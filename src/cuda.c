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

#include "err.h"

#include "diff_alg/diff_tensor.h"
#include "potentials/electro.h"
#include "potentials/LJ.h"

MAKE_STR_IN(PTR_INT,cuda_devices,NULL,"Cuda devices used")
MAKE_STR_IN(INT,cuda_block,256,"The dimension of the thread block")

#if USE_CUDA || USE_MAGMA
#include <cuda_runtime_api.h>
#include <cublas.h>

INT sm_red;
DOUBLE* d_coord = NULL;
DOUBLE* d_scoord = NULL;
DOUBLE* d_D = NULL;
INT    size_D = 0;
INT    size_coord = 0;
DOUBLE* d_X = NULL;
DOUBLE* h_X = NULL;
DOUBLE* d_Rs0 = NULL;
DOUBLE* h_Rs0 = NULL;
DOUBLE* d_sRs0 = NULL;
DOUBLE* d_sforces = NULL;
DOUBLE* h_sforces = NULL;
DOUBLE* d_forces = NULL;
DOUBLE* h_forces = NULL;

DOUBLE* d_acoord = NULL;
DOUBLE* d_aRs0 = NULL;
DOUBLE* d_aD = NULL;

INT blocks;

void cudaCheckError( const char* file, int line )
{
    cudaError_t error = cudaGetLastError();
    if( error != 0 )
    {
        warning( cudaGetErrorString(error), file, line );
        unexpected_error( "End execution due to cuda error", file, line );
    }
}

void cublasCheckError( const char* file, int line )
{
    INT error = cublasGetError();
    if( error != 0 )
    {
        warning( "CUBLAS error", file, line );
        UNERR("End execution due to cuda error");
    }
}

void init_cuda()
{
    INT id = 0;
    struct cudaDeviceProp deviceProp;
    if( cuda_devices == NULL )
    {
        g_devices = 1;
        cudaSetDevice( 0 ); CCERR
        cublasInit(); CCERR
    }
    else
    {
        if( variable_cuda_devices.size == 1 )
        {
            id = cuda_devices[0];
            cudaSetDevice( cuda_devices[0] ); CCERR
            cublasInit(); CCERR
        }
        else
        {
            UNERR("Unimplemented");
        }
    }
    cudaGetDeviceProperties( &deviceProp, id ); CCERR
    sm_red = deviceProp.multiProcessorCount;
    
    size_D = size * size * DIMS1 * DIMS0 * sizeof( float );
    size_coord = size * DIMS1 * sizeof( float );
    if( hydro && !geyer_on_the_fly )
    {
        cudaMalloc( (void**)&d_D, size_D ); CCERR
    }
    else
    {
        cudaMalloc( (void**)&d_D, size * DIMS0 * sizeof( float ) ); CCERR
    }

    cudaMalloc( (void**)&d_coord, size_coord ); CCERR
    cudaMalloc( (void**)&d_scoord, size_coord ); CCERR

    h_X = (DOUBLE*) malloc( size_coord ); CHMEM(h_X);
    h_Rs0 = (DOUBLE*) malloc( size_coord ); CHMEM(h_Rs0);

    cudaMalloc( (void**)&d_X, size_coord ); CCERR
    cudaMalloc( (void**)&d_Rs0, size_coord ); CCERR
    cudaMalloc( (void**)&d_sRs0, size_coord ); CCERR

    blocks = size / cuda_block;
    if ( blocks * cuda_block != size )
    {
        printf( "%d (%d != %d)\n", blocks, blocks*cuda_block, size );
        UNERR("Wrong size, size should be multiplication of cuda_block");
    }
    if ( algorithm > DIFF_ALG_ERMAK_NEWTON )
    {
        cudaMalloc( (void**)&d_acoord, size_coord ); CCERR
        cudaMalloc( (void**)&d_aRs0, size_coord ); CCERR
        cudaMalloc( (void**)&d_aD, size_D ); CCERR
    }
    cudaMemcpy( d_coord, coord, size_coord, cudaMemcpyHostToDevice ); CCERR
}

void free_cuda()
{
    cudaFree( d_coord ); CCERR
    if ( hydro )
    {
        cudaFree( d_D ); CCERR
        if ( d_acoord ) { cudaFree( d_acoord ); CCERR }
        if ( d_aRs0 )   { cudaFree( d_aRs0   ); CCERR }
        if ( d_aD )     { cudaFree( d_aD     ); CCERR }
    }
    free( h_X );
    cudaFree( d_X ); CCERR
    free( h_Rs0 );
    cudaFree( d_Rs0 ); CCERR

    cublasShutdown(); CCERR
    cudaThreadExit(); CCERR
}

void init_cuda_after_restart()
{
    init_cuda();
}

void show_array( float* array, int count, const char* name )
{
    if ( name ) printf( "%s ", name );
    if ( array )
    {
        int csf = count * sizeof(float);
        int i;
        float* td = (float*) malloc( csf ); CHMEM(td);
        cudaMemcpy( td, array, csf, cudaMemcpyDeviceToHost );
        for ( i = 0; i < count; ++i )
        {
            printf( "%g ", td[i] );
        }
        free(td);
    }
    else
    {
        printf("NULL");
    }
    printf("\n");
}

#endif
