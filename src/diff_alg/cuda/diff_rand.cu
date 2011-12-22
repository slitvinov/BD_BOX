/**************************************************************************
 *                                                                        *
 *   Copyright (C) 2010 - 2011: University of Warsaw                      *
 *   Authors:            Pawel Zielinski <pzielins@icm.edu.pl>            *
 *                       Maciej Dlugosz  <mdlugosz@icm.edu.pl>            *
 *                                                                        *
 *   This program is free software; you can redistribute it and/or modify *
 *   it under the terms of the GNU General Public License as published by *
 *   the Free Software Foundation; either version 3 of the License, or    *
 *   (at your option) any later version.                                  *
 *                                                                        *
 *   This program is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
 *   GNU General Public License for more details.                         *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program; if not, see http://www.gnu.org/licenses     *
 *   or write to the Free Software Foundation,Inc., 51 Franklin Street,   *
 *   Fifth Floor, Boston, MA 02110-1301  USA                              *
 *                                                                        *
 **************************************************************************/

#include "cuda.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

extern "C"
{
#include "../diff_tensor.h"
#include "../../cuda.h"
#include "../../err.h"
#include "../../trans.h"
#include "../../rand_move.h"
#include "../../myblas.h"
}

#include "diff_tensor.cuh"

void mstart();
void mstop();

DOUBLE* d_eps = NULL;
DOUBLE* h_eps = NULL;
DOUBLE* d_eff = NULL;
DOUBLE* d_C   = NULL;

template<int cuda_block>
__device__ void eps_reduction( float e_curr, float* out )
{
    extern __shared__ float sh[];
    int tid = threadIdx.x;

    __syncthreads( );
    sh[tid] = e_curr;
    __syncthreads( );
    if ( cuda_block >= 512 )
    {
        if ( tid < 256 ) sh[tid] += sh[tid + 256];
        __syncthreads( );
    }
    if ( cuda_block >= 256 )
    {
        if ( tid < 128 ) sh[tid] += sh[tid + 128];
        __syncthreads( );
    }
    if ( cuda_block >= 128 )
    {
        if ( tid < 64 ) sh[tid] += sh[tid + 64];
        __syncthreads( );
    }
    volatile float *redE = sh;
    if ( tid < 32 )
    {
        if ( cuda_block >= 64 ) redE[tid] += redE[tid + 32];
        if ( cuda_block >= 32 ) redE[tid] += redE[tid + 16];
        if ( cuda_block >= 16 ) redE[tid] += redE[tid + 8];
        if ( cuda_block >= 8 ) redE[tid] += redE[tid + 4];
        if ( cuda_block >= 4 ) redE[tid] += redE[tid + 2];
        if ( cuda_block >= 2 ) redE[tid] += redE[tid + 1];
    }
    if ( tid == 0 )
    {
        out[blockIdx.x] = redE[0];
    }
}

template<int cuda_block >
__global__ void comp_geyer_beg( DOUBLE* ptr_d, DOUBLE* ptr_x, DOUBLE* ptr_eff, DOUBLE* ptr_c, DOUBLE* ptr_eps )
{
    float eps_curr = 0;
    float Feff = 0;
    float C = 0;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int LD = (gridDim.x * blockDim.x/3)*4;

    for ( int k = 0; k < gridDim.x*blockDim.x; ++k )
    {
        float dk = ptr_d[LD*k+k];
        float dik = 0;
        if(id<k)
            dik = ptr_d[LD*id + k];
        else
            dik = ptr_d[LD*k + id];
        if(k!=id)
        {
            Feff += ptr_x[k] * dik;
            C += dik*dik/dk;
            eps_curr+=dik;
        }
    }
    ptr_eff[id] = Feff;
    ptr_c[id] = C;
    eps_reduction<cuda_block>( eps_curr/2, ptr_eps );
}

template<int cuda_block >
__global__ void comp_geyer_end( DOUBLE* ptr_d, DOUBLE* ptr_x, DOUBLE* ptr_eff, DOUBLE* ptr_c, DOUBLE beta, DOUBLE* out)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    float Feff = ptr_eff[id];
    float C = ptr_c[id];
    int LD = (gridDim.x * blockDim.x/DIMS0)*DIMS1;
    float dii = ptr_d[LD*id+id];
    C *= beta*beta/dii;
    C = rsqrtf( 1 + C );
    out[id] += C*(beta*Feff + dii*ptr_x[id]);
}

extern "C"
void compute_r( DOUBLE* _D )
{
    INT i;
    DOUBLE a, eps;

    gauss_rand( X, size * DIMS0 );
    switch ( hydro )
    {
        case DIFF_NONE:
        {
            for ( i = 0; i < size * DIMS0; ++i )
            {
                X[ i ] = X[ i ] * SQRTD( diag_D[i] * 2 * curr_dt );
            }
            cudaMemcpy( d_X, X, sizeof(DOUBLE)*size*DIMS0, cudaMemcpyHostToDevice ); CCERR
            CUBLAS_AXPYD( size*DIMS0, 1.0f, d_X, 1, d_Rs0, 1 ); cublasCheckError(__FILE__,__LINE__);
            break;
        }
        case DIFF_CHOLESKY:
        {
            /* computation on CPU */

            a = SQRTD( 2 * curr_dt );
            for ( i = 0; i < size * DIMS0; ++i )
            {
                X[i] *= a;
            }
            cudaMemcpy( d_X, X, sizeof(DOUBLE)*size*DIMS0, cudaMemcpyHostToDevice ); CCERR
            if( curr_iter != save_iter )
            {
                if ( my_dpotrf( _D, size*DIMS0, size*DIMS1 ) )
                    UNERR("Cholesky decomposition error")
            }
            my_dtrmv( _D, d_X, d_Rs0, size*DIMS0, size*DIMS1 );
            break;
        }
        case DIFF_GEYER:
        {
            for ( i = 0; i < size * DIMS0; ++i )
            {
                X[ i ] *= SQRTD( 2 * curr_dt * inv_diag_D[i] );
            }
            cudaMemcpy( d_X, X, sizeof(DOUBLE)*size*DIMS0, cudaMemcpyHostToDevice ); CCERR

            int to_get = DIMS0*blocks;
            if ( geyer_on_the_fly )
            {
                comp_geyer_beg_in_fly( d_coord, d_X, d_eff, d_C, d_eps );
                to_get = blocks;
            }
            else
            {
                switch(cuda_block)
                {
                    case(1024):comp_geyer_beg<1024><<< DIMS0*blocks, cuda_block, 4*cuda_block >>>( d_D, d_X, d_eff, d_C, d_eps );break;
                    case( 512):comp_geyer_beg< 512><<< DIMS0*blocks, cuda_block, 4*cuda_block >>>( d_D, d_X, d_eff, d_C, d_eps );break;
                    case( 256):comp_geyer_beg< 256><<< DIMS0*blocks, cuda_block, 4*cuda_block >>>( d_D, d_X, d_eff, d_C, d_eps );break;
                    case( 128):comp_geyer_beg< 128><<< DIMS0*blocks, cuda_block, 4*cuda_block >>>( d_D, d_X, d_eff, d_C, d_eps );break;
                    case(  64):comp_geyer_beg<  64><<< DIMS0*blocks, cuda_block, 4*cuda_block >>>( d_D, d_X, d_eff, d_C, d_eps );break;
                    case(  32):comp_geyer_beg<  32><<< DIMS0*blocks, cuda_block, 4*cuda_block >>>( d_D, d_X, d_eff, d_C, d_eps );break;
                }
            }
            cudaThreadSynchronize(); CCERR
            cudaMemcpy( h_eps, d_eps, to_get*sizeof(float),cudaMemcpyDeviceToHost); CCERR
            eps=0;
            for( int i = 0; i < to_get; ++i )
            {
                eps += h_eps[i];
            }
            //eps /= (size * DIMS0 * size * DIMS0 - size * DIMS0) / 2;
            eps /= size * DIMS0;
            eps /= size * DIMS0 - 1;
            eps *= 2;
            if ( eps > 1.0 )
            {
                eps = 1.0;
                warning( "epsilon truncated to 1.0", __FILE__, __LINE__ );
            }
            a = (3 * size - 1) * eps * eps-(3 * size - 2) * eps;
            DOUBLE beta = 0.5;
            if ( abs( a ) < 0.0000001f )/*we assume no correlation*/
            {
                warning( "Epsilon == 0.0", __FILE__, __LINE__ );
            }
            else
            {
                beta = (1 - SQRTD( 1 - a )) / a;
            }
            if ( hydro == DIFF_GEYER && algorithm < DIFF_ALG_IGT_CONST )
            {
                comp_geyer_end_in_fly( d_coord, d_X, d_eff, d_C, beta, d_Rs0 );
            }
            else
            {
                switch(cuda_block)
                {
                    case(1024):comp_geyer_end<1024><<< DIMS0*blocks, cuda_block >>>( d_D, d_X, d_eff, d_C, beta, d_Rs0 );break;
                    case( 512):comp_geyer_end< 512><<< DIMS0*blocks, cuda_block >>>( d_D, d_X, d_eff, d_C, beta, d_Rs0 );break;
                    case( 256):comp_geyer_end< 256><<< DIMS0*blocks, cuda_block >>>( d_D, d_X, d_eff, d_C, beta, d_Rs0 );break;
                    case( 128):comp_geyer_end< 128><<< DIMS0*blocks, cuda_block >>>( d_D, d_X, d_eff, d_C, beta, d_Rs0 );break;
                    case(  64):comp_geyer_end<  64><<< DIMS0*blocks, cuda_block >>>( d_D, d_X, d_eff, d_C, beta, d_Rs0 );break;
                    case(  32):comp_geyer_end<  32><<< DIMS0*blocks, cuda_block >>>( d_D, d_X, d_eff, d_C, beta, d_Rs0 );break;
                }
            }
            cudaThreadSynchronize(); CCERR
            break;
        }
        case DIFF_CHEBYSHEV:
        {
            UNERR( "Unimplemented CHEBYSHEV method")
        }
    }
    save_iter = curr_iter;
}

extern "C"
void init_cuda_rand()
{
    if ( hydro == DIFF_GEYER )
    {
        cudaMalloc( (void**) &d_eps, sizeof(DOUBLE)*DIMS0*size); CCERR
        h_eps = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS0*size ); CHMEM(h_eps);
        cudaMalloc( (void**) &d_eff, sizeof(DOUBLE)*DIMS0*size); CCERR
        cudaMalloc( (void**) &d_C, sizeof(DOUBLE)*DIMS0*size); CCERR
    }
}

extern "C"
void free_cuda_rand()
{
    cudaFree( d_eps ); CCERR
    free( h_eps );
    cudaFree( d_eff ); CCERR
    cudaFree( d_C ); CCERR
}
