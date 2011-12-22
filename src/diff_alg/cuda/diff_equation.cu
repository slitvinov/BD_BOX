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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cublas.h>

extern "C"
{
#include "../diff_tensor.h"
#include "../../cuda.h"
#include "../../input.h"
#include "../../main_loop.h"
#include "../../math_help.h"
#include "../../myblas.h"
#include "../../rand_move.h"
#include "../../trans.h"

#include "../../potentials/calc_func.h"
#include "../../potentials/LJ.h"
#include "../../potentials/electro.h"
void mstart();
void mstop();
}

#include "diff_tensor.cuh"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

__constant__ float dev_box0_e;
__constant__ float dev_box1_e;
__constant__ float dev_box2_e;
__constant__ float dev_invL0_e;
__constant__ float dev_invL1_e;
__constant__ float dev_invL2_e;
__constant__ float dev_grad[9];

DOUBLE* d_FS = NULL;

__global__ void vec_mult( void* ptr_a, void* ptr_b, void* ptr_c )
{
    float* a = (float*)ptr_a;
    float* b = (float*)ptr_b;
    float* c = (float*)ptr_c;
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/

    c[id] = a[id] * b[id];
}

void compute_Rf( DOUBLE** _Rf, DOUBLE* _D )
{
    DOUBLE pre = curr_dt / (kB * T);
    INT LD = size*DIMS1;
    
    DOUBLE dzero = 0;
    DOUBLE done = 1;
    int ione = 1;
    if( curr_iter != save_iter )
    {
        if( geyer_on_the_fly )
        {
            comp_rinfly( d_coord, d_forces, _Rf[0] );
            cudaThreadSynchronize(); CCERR
        }
        else if( hydro )
        {
            compute_D( _D ); CCERR
            if (d_forces)
            {
                CUBLAS_SYMVD( 'L', size*DIMS0, done, _D, LD, d_forces, ione, dzero, _Rf[0], ione ); CCERR
            }
            else
            {
                cudaMemset( _Rf[0], 0, sizeof(float)*size*DIMS0 );
            }
        }
        else
        {
            if ( d_forces )
            {
                vec_mult<<<blocks*DIMS0, cuda_block>>>( _D, d_forces, _Rf[0] );
                cudaThreadSynchronize(); CCERR
            }
            else
            {
                cudaMemset( _Rf[0], 0, sizeof(float)*size*DIMS0 ); CCERR
            }
        }
        CUBLAS_SSCALD( size * DIMS0, pre, _Rf[0], 1); CCERR
        cudaMemcpy( d_sRs0, _Rf[0], size*DIMS0*sizeof(float), cudaMemcpyDeviceToDevice ); CCERR
    }
    else
    {
        cudaMemcpy( _Rf[0], d_sRs0, size*DIMS0*sizeof(float), cudaMemcpyDeviceToDevice ); CCERR
    }
}

template<bool pbc, bool grad>
__global__ void apply_move( void* ptr_c, void* ptr_v )
{
    float4* coord = (float4*)ptr_c;
    float3* move = (float3*)ptr_v;
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/

    float4 curr_c = coord[id];
    float3 curr_m = move[id];

    if ( grad )
    {
        float tmp1 = dev_grad[0]*curr_c.x + dev_grad[1]*curr_c.y + dev_grad[2]*curr_c.z;
        float tmp2 = dev_grad[3]*curr_c.x + dev_grad[4]*curr_c.y + dev_grad[5]*curr_c.z;
        float tmp3 = dev_grad[6]*curr_c.x + dev_grad[7]*curr_c.y + dev_grad[8]*curr_c.z;
        curr_c.x += tmp1;
        curr_c.y += tmp2;
        curr_c.z += tmp3;
    }

    curr_c.x += curr_m.x;
    curr_c.y += curr_m.y;
    curr_c.z += curr_m.z;

    if ( pbc )
    {
        curr_c.x -= dev_box0_e * floor( curr_c.x * dev_invL0_e + 0.5f );
        curr_c.y -= dev_box1_e * floor( curr_c.y * dev_invL1_e + 0.5f );
        curr_c.z -= dev_box2_e * floor( curr_c.z * dev_invL2_e + 0.5f );
    }
    coord[id] = curr_c;
}

template<bool pbc, bool first >
__global__ void apply_move_grad( void* ptr_c, void* ptr_v, float3* FS )
{
    float4* coord = (float4*)ptr_c;
    float3* move = (float3*)ptr_v;
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/

    float4 curr_c = coord[id];
    float3 curr_m = move[id];

    float3 tmp;
    if ( first )
    {
        tmp.x = dev_grad[0]*curr_c.x + dev_grad[1]*curr_c.y + dev_grad[2]*curr_c.z;
        tmp.y = dev_grad[3]*curr_c.x + dev_grad[4]*curr_c.y + dev_grad[5]*curr_c.z;
        tmp.z = dev_grad[6]*curr_c.x + dev_grad[7]*curr_c.y + dev_grad[8]*curr_c.z;
        curr_c.x += tmp.x;
        curr_c.y += tmp.y;
        curr_c.z += tmp.z;
    }

    curr_c.x += curr_m.x;
    curr_c.y += curr_m.y;
    curr_c.z += curr_m.z;
    if ( first )
    {
        tmp.x += dev_grad[0]*curr_c.x + dev_grad[1]*curr_c.y + dev_grad[2]*curr_c.z;
        tmp.y += dev_grad[3]*curr_c.x + dev_grad[4]*curr_c.y + dev_grad[5]*curr_c.z;
        tmp.z += dev_grad[6]*curr_c.x + dev_grad[7]*curr_c.y + dev_grad[8]*curr_c.z;
        FS[id] = tmp;
    }
    else
    {
        tmp = FS[id];
        curr_c.x += tmp.x * 0.5f;
        curr_c.y += tmp.y * 0.5f;
        curr_c.z += tmp.z * 0.5f;
    }
    if ( pbc )
    {
        curr_c.x -= dev_box0_e * floor( curr_c.x * dev_invL0_e + 0.5f );
        curr_c.y -= dev_box1_e * floor( curr_c.y * dev_invL1_e + 0.5f );
        curr_c.z -= dev_box2_e * floor( curr_c.z * dev_invL2_e + 0.5f );
    }
    coord[id] = curr_c;
}

void compute_hydro_ermak( DOUBLE* coord_out )
{
    compute_Rf( &d_Rs0, d_D );

    if ( (!geyer_on_the_fly && algorithm > DIFF_ALG_ERMAK_NEWTON) )
    {
        cudaMemcpy( d_aD, d_D, size_D, cudaMemcpyDeviceToDevice ); CCERR
    }
    if ( algorithm > DIFF_ALG_ERMAK_NEWTON )
    {
        cudaMemcpy( d_aRs0, d_Rs0, sizeof(float)*size*DIMS0, cudaMemcpyDeviceToDevice ); CCERR
    }
    compute_r( d_D );
    if ( !vel_grad_tensor )
    {
        if(bc==BC_BOX)
            apply_move<  true, false ><<<blocks, cuda_block>>>( coord_out, d_Rs0 );
        else
            apply_move< false, false ><<<blocks, cuda_block>>>( coord_out, d_Rs0 );
    }
    else
    {
        if ( algorithm > DIFF_ALG_ERMAK_NEWTON ) //for IGT
        {
            if(bc==BC_BOX)
                apply_move_grad< true, true ><<<blocks, cuda_block>>>( coord_out, d_Rs0, (float3*)d_FS );
            else
                apply_move_grad< false, true ><<<blocks, cuda_block>>>( coord_out, d_Rs0, (float3*)d_FS );
        }
        else
        {
            if(bc==BC_BOX)
                apply_move<  true, true ><<<blocks, cuda_block>>>( coord_out, d_Rs0 );
            else
                apply_move< false, true ><<<blocks, cuda_block>>>( coord_out, d_Rs0 );
        }
    }
    cudaThreadSynchronize(); CCERR
}

INT compute_hydro_IGT_second_phase( DOUBLE* _E )
{
    cudaMemcpy( d_Rs0, d_aRs0, sizeof(float)*size*DIMS0, cudaMemcpyDeviceToDevice );
    compute_r( d_D );
    if ( !vel_grad_tensor )
    {
        if(bc==BC_BOX)
            apply_move<  true, false ><<<blocks, cuda_block>>>( d_coord, d_Rs0 );
        else
            apply_move< false, false ><<<blocks, cuda_block>>>( d_coord, d_Rs0 );
    }
    else
    {
        if(bc==BC_BOX)
            apply_move_grad< true, false ><<<blocks, cuda_block>>>( d_coord, d_Rs0, (float3*)d_FS );
        else
            apply_move_grad< false, false ><<<blocks, cuda_block>>>( d_coord, d_Rs0, (float3*)d_FS );
    }
    cudaThreadSynchronize(); CCERR
    return compute_forces( _E, curr_time + curr_dt ) << 1;
}

__global__ void comp_mean_D( float* D, float* Ds )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    const int LD = gridDim.x * blockDim.x * DIMS1;

    if ( idx <= idy )//lower traingle
    {
        D[ DIMS0 * idx * LD + DIMS0 * idy ] =
            (D [ DIMS0 * idx * LD + DIMS0 * idy ] +
             Ds[ DIMS0 * idx * LD + DIMS0 * idy ])/2;
    }
}

INT compute_hydro_IGT( DOUBLE* _E )
{
    float* tmp_coord;
    cudaMemcpy( d_acoord, d_coord, sizeof(float)*size_coord, cudaMemcpyDeviceToDevice ); CCERR
    compute_hydro_ermak( d_acoord );
    
    tmp_coord = d_coord;
    d_coord = d_acoord;

    if ( compute_forces( NULL, curr_time + curr_dt ) )
    {
        d_coord = tmp_coord;
        return 1;
    }
    curr_iter++;

    compute_Rf( &d_Rs0, d_D );
    d_coord = tmp_coord;

    cublasSaxpy( size*DIMS0, 1.0f, d_aRs0, 1, d_Rs0, 1 ); cublasCheckError(__FILE__,__LINE__);
    cublasSscal( size*DIMS0, 0.5f,  d_Rs0, 1 );           cublasCheckError(__FILE__,__LINE__);

    if ( hydro )
    {
        const int xtdim = 1;
        const int ytdim = 128;
        dim3 dim_block( xtdim, ytdim );
        dim3 dim_grid( size / xtdim, size / ytdim );
        comp_mean_D<<< dim_grid, dim_block >>>( d_D, d_aD );
        cudaThreadSynchronize(); CCERR
    }
    cudaMemcpy( d_aRs0, d_Rs0, sizeof(float)*size*DIMS0, cudaMemcpyDeviceToDevice );

    return compute_hydro_IGT_second_phase( _E );
}

void init_cuda_equation()
{
    if ( bc == BC_BOX )
    {
        float bx[3];
        float iL[3];
        bx[0] = box[0], bx[1] = box[1], bx[2] = box[2];
        iL[0] = inv_L[0], iL[1] = inv_L[1], iL[2] = inv_L[2];
        cudaMemcpyToSymbol( "dev_box0_e", bx, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
        cudaMemcpyToSymbol( "dev_box1_e", bx + 1, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
        cudaMemcpyToSymbol( "dev_box2_e", bx + 2, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
        cudaMemcpyToSymbol( "dev_invL0_e", iL, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
        cudaMemcpyToSymbol( "dev_invL1_e", iL + 1, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
        cudaMemcpyToSymbol( "dev_invL2_e", iL + 2, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    }
    if ( vel_grad_tensor )
    {
        float tmp[9];
        for ( int i = 0; i < 9; ++i )
        {
            tmp[i] = vel_grad_tensor[i];
        }
        cudaMemcpyToSymbol ( "dev_grad", tmp, sizeof(tmp), 0, cudaMemcpyHostToDevice ); CCERR
        if ( algorithm > DIFF_ALG_ERMAK_NEWTON ) //for IGT
        {
            cudaMalloc( (void**)&d_FS, sizeof(float)*size*DIMS0 ); CCERR
        }
    }
}

void free_cuda_equation()
{
    if ( d_FS ) { cudaFree( d_FS ); CCERR }
}
