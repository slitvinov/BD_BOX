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
#include <stdlib.h>
#include <stdio.h>

extern "C"
{
#include "../newton.h"
#include "../diff_tensor.h"
#include "../../err.h"
#include "../../input.h"
#include "../../math_help.h"
#include "../../trans.h"
#include "../../cuda.h"
}

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

DOUBLE* d_velocity = NULL;
DOUBLE* h_velocity = NULL;
STR d_rev = NULL;
STR h_rev = NULL;
DOUBLE* d_times = NULL;
DOUBLE* h_times = NULL;
DOUBLE s_time;

__constant__ float dev_box0;
__constant__ float dev_box1;
__constant__ float dev_box2;
__constant__ float dev_invL0;
__constant__ float dev_invL1;
__constant__ float dev_invL2;

template<bool pbc>
__device__ float single_time( float4 ci, float4 vi, float4 cj, float4 vj )
{
    float3 r;
    r.x = ci.x - cj.x;
    r.y = ci.y - cj.y;
    r.z = ci.z - cj.z;

    if ( pbc )
    {
        r.x -= dev_box0 * floorf( r.x * dev_invL0 + 0.5f );
        r.y -= dev_box1 * floorf( r.y * dev_invL1 + 0.5f );
        r.z -= dev_box2 * floorf( r.z * dev_invL2 + 0.5f );
    }
    float3 v;
    v.x = vi.x - vj.x;
    v.y = vi.y - vj.y;
    v.z = vi.z - vj.z;
    float rv = r.x * v.x + r.y * v.y + r.z * v.z;
    if ( rv < 0 )
    {
        float l = vi.w + vj.w + NEWTON_OFFSET; /* velocity.w contain radii_LJ */
        float vv = v.x*v.x + v.y*v.y + v.z*v.z;
        float rr = r.x*r.x + r.y*r.y + r.z*r.z;
        float delta = rv * rv - vv * ( rr - l * l);
        if ( delta > 0 )
        {
            float inv_vv = 1.0f / vv;
            delta = sqrtf( delta );
            float t2 = -(rv + delta) * inv_vv;
            return (t2 < 0.0f) ? 0.0f : t2 ;
        }
    }
    return FLT_MAX;
}

template<bool pbc>
__global__ void comp_timeD( void* ptr_coord, void* ptr_velocity, void* ptr_times )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/
    extern __shared__ float sh[];
    float4* shCoord = (float4*) sh;
    float4* shVel = ((float4*) sh)+ blockDim.x;
    float4* coord = (float4*)ptr_coord;
    float4* vel = (float4*)ptr_velocity;
    
    float4 currCoord = coord[id];
    float4 currVel = vel[id];
    
    float out_t = FLT_MAX;
    float out_pos = -1;
    for ( int i = 0; i < blockIdx.x; ++i )
    {
        shCoord[ threadIdx.x ] = coord[ i * blockDim.x + threadIdx.x ];
        shVel[ threadIdx.x ] = vel[ i * blockDim.x + threadIdx.x ];
        __syncthreads( );
        for ( int j = 0; j < blockDim.x; ++j )
        {
            float t = single_time<pbc> (currCoord, currVel, shCoord[j], shVel[j] );
            if( out_t > t )
            {
                out_pos = i * blockDim.x + j;
                out_t = t;
            }
        }
        __syncthreads( );
    }

    shCoord[ threadIdx.x ] = currCoord;
    shVel[threadIdx.x] = currVel;
    __syncthreads( );
    for ( int j = 0; j < threadIdx.x; ++j )
    {
        float t = single_time<pbc> (currCoord, currVel, shCoord[j], shVel[j] );
        if( out_t > t )
        {
            out_pos = blockIdx.x * blockDim.x + j;
            out_t = t;
        }
    }
    float* outTab = (float*)ptr_times;
    outTab[id] = out_t;
    outTab[id + blockDim.x*gridDim.x] = out_pos;//blockDim.x*gridDim.x == size
}

void comp_time( )
{
    int i;
    if( bc == BC_BOX )
        comp_timeD<true> <<< blocks, cuda_block, cuda_block * 4* 4 *2 >>>( d_coord, d_velocity, d_times );
    else
        comp_timeD<false> <<< blocks, cuda_block, cuda_block * 4* 4 *2 >>>( d_coord, d_velocity, d_times );
    cudaThreadSynchronize(); CCERR
    cudaMemcpy( h_times, d_times, sizeof(DOUBLE) * size * 2, cudaMemcpyDeviceToHost ); CCERR
    s_time = DOUBLE_MAX;
    for( i = 0; i < size; ++i )
    {
        if( s_time > h_times[i] ) s_time = h_times[i];
    }
}

template<bool pbc>
__global__ void apply_velocities( void* ptr_c, void* ptr_v, DOUBLE t )
{
    float4* coord = (float4*)ptr_c;
    float4* velocity = (float4*)ptr_v;
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/

    float4 curr_c = coord[id];
    float4 curr_v = velocity[id];

    curr_c.x += curr_v.x * t;
    curr_c.y += curr_v.y * t;
    curr_c.z += curr_v.z * t;
    if ( pbc )
    {
        curr_c.x -= dev_box0 * floor( curr_c.x * dev_invL0 + 0.5f );
        curr_c.y -= dev_box1 * floor( curr_c.y * dev_invL1 + 0.5f );
        curr_c.z -= dev_box2 * floor( curr_c.z * dev_invL2 + 0.5f );
    }
    coord[id] = curr_c;
}

template<int cuda_block>
__global__ void apply_velocities( void* ptr_c, void* ptr_v, DOUBLE t, DOUBLE radius, void* ptr_rev  )
{
    float4* coord = (float4*)ptr_c;
    float4* velocity = (float4*)ptr_v;
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/
    int tid = threadIdx.x;

    float4 curr_c = coord[id];
    float4 curr_v = velocity[id];

    curr_c.x += curr_v.x * t;
    curr_c.y += curr_v.y * t;
    curr_c.z += curr_v.z * t;

    coord[id] = curr_c;
    
    char rev = (curr_c.x*curr_c.x + curr_c.y*curr_c.y + curr_c.z*curr_c.z) > radius;
    extern __shared__ char rev_tab[];
    rev_tab[threadIdx.x] = rev;
    __syncthreads();
    if( cuda_block >= 512 ) { if( tid < 256 ) rev_tab[tid] |= rev_tab[tid+256]; __syncthreads(); }
    if( cuda_block >= 256 ) { if( tid < 128 ) rev_tab[tid] |= rev_tab[tid+128]; __syncthreads(); }
    if( cuda_block >= 128 ) { if( tid <  64 ) rev_tab[tid] |= rev_tab[tid+ 64]; __syncthreads(); }
    if( tid < 32 )
    {
        if( cuda_block >= 64 ) rev_tab[tid] |= rev_tab[tid+32];
        if( cuda_block >= 32 ) rev_tab[tid] |= rev_tab[tid+16];
        if( cuda_block >= 16 ) rev_tab[tid] |= rev_tab[tid+ 8];
        if( cuda_block >=  8 ) rev_tab[tid] |= rev_tab[tid+ 4];
        if( cuda_block >=  4 ) rev_tab[tid] |= rev_tab[tid+ 2];
        if( cuda_block >=  2 ) rev_tab[tid] |= rev_tab[tid+ 1];
    }
    if( tid == 0 ) ((char*)ptr_rev)[blockIdx.x] = rev_tab[0];
}

INT newton_moves( INT* pcount )
{
    DOUBLE next_move_time;
    DOUBLE elapsed = 0.0;
    INT ret = 0,p2,p1;
    while ( elapsed < dt - 2 * NEWTON_EPSILON )
    {
		ret = 0;
        comp_time( );
        next_move_time = s_time;
        if ( next_move_time < 0 || next_move_time > dt - elapsed )/* we can move*/
        {
            next_move_time = dt - elapsed;
        }
        next_move_time = (next_move_time > NEWTON_EPSILON) ? next_move_time - NEWTON_EPSILON : 0;
        elapsed += next_move_time;
        if ( next_move_time > 0.0 )
        {
            /*apply velocities*/
            if( bc == BC_SPHERE )
            {
                INT i;
                DOUBLE sr = sphere_radius * sphere_radius;
                switch (cuda_block)
                {
                    case 1024: apply_velocities<1024> <<<blocks, cuda_block, cuda_block>>>( d_coord, d_velocity, next_move_time, sr, d_rev ); break;
                    case  512: apply_velocities< 512> <<<blocks, cuda_block, cuda_block>>>( d_coord, d_velocity, next_move_time, sr, d_rev ); break;
                    case  256: apply_velocities< 256> <<<blocks, cuda_block, cuda_block>>>( d_coord, d_velocity, next_move_time, sr, d_rev ); break;
                    case  128: apply_velocities< 128> <<<blocks, cuda_block, cuda_block>>>( d_coord, d_velocity, next_move_time, sr, d_rev ); break;
                    case   64: apply_velocities<  64> <<<blocks, cuda_block, cuda_block>>>( d_coord, d_velocity, next_move_time, sr, d_rev ); break;
                    case   32: apply_velocities<  32> <<<blocks, cuda_block, cuda_block>>>( d_coord, d_velocity, next_move_time, sr, d_rev ); break;
                }
                cudaMemcpy( h_rev, d_rev, blocks,cudaMemcpyDeviceToHost ); CCERR
                cudaThreadSynchronize(); CCERR
                for( i = 0; i < blocks; ++i )
                {
                    ret |= h_rev[i];
                }
            }
            else
            {
                if( bc == BC_BOX )
                    apply_velocities<true><<<blocks, cuda_block>>>( d_coord, d_velocity, next_move_time );
                else
                    apply_velocities<false><<<blocks, cuda_block>>>( d_coord, d_velocity, next_move_time );
                cudaThreadSynchronize(); CCERR
            }
        }
        if ( ret )/*from pbc==sphere*/
        {
            break;
        }
        
        /*revert speeds */
        cudaMemcpy( h_velocity, d_velocity, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
        if( !d_bforces )
        {
            cudaMemcpy( coord, d_coord, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
        }//else copied for bonded interactions
        for( p1 = 0; p1 < size; ++p1 )
        {
            if( h_times[p1] == s_time )
            {
                p2 = (INT)h_times[size+p1];
                DOUBLE r_norm, r_sq;
                DOUBLE e[3], v[3];
                DOUBLE ev;
                DOUBLE m1 = masses[p1];
                DOUBLE m2 = masses[p2];
                get_v2_norms( p2, p1, e, &r_norm, &r_sq );
                dist_vec( h_velocity + DIMS1*p2, h_velocity + DIMS1*p1, v );
                ev = e[0] * v[0] + e[1] * v[1] + e[2] * v[2];
                ev *= (1 + e_collision) / (m1 + m2);
                h_velocity[DIMS1 * p1 ] -= m2 * ev * e[0];
                h_velocity[DIMS1 * p1 + 1] -= m2 * ev * e[1];
                h_velocity[DIMS1 * p1 + 2] -= m2 * ev * e[2];

                h_velocity[DIMS1 * p2 ] += m1 * ev * e[0];
                h_velocity[DIMS1 * p2 + 1] += m1 * ev * e[1];
                h_velocity[DIMS1 * p2 + 2] += m1 * ev * e[2];
            }
        }
        cudaMemcpy( d_velocity, h_velocity,
                    sizeof(DOUBLE)*size*DIMS1, cudaMemcpyHostToDevice ); CCERR
        ++(*pcount);
        if ( *pcount > move_attempts ) UNERR( "Ermak Newton error")
    }/* while */
    return ret;
}

template<bool pbc>
__global__ void comp_velocities( void* ptr_c, void* ptr_sc, void* ptr_v, DOUBLE inv_t )
{
    float4* coord = (float4*)ptr_c;
    float4* scoord = (float4*)ptr_sc;
    float4* velocity = (float4*)ptr_v;
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/

    float4 curr_c = coord[id];
    float4 curr_sc = scoord[id];

    float3 r_ij;
    r_ij.x = curr_c.x - curr_sc.x;
    r_ij.y = curr_c.y - curr_sc.y;
    r_ij.z = curr_c.z - curr_sc.z;

    if ( pbc )
    {
        r_ij.x -= dev_box0 * floor( r_ij.x * dev_invL0 + 0.5f );
        r_ij.y -= dev_box1 * floor( r_ij.y * dev_invL1 + 0.5f );
        r_ij.z -= dev_box2 * floor( r_ij.z * dev_invL2 + 0.5f );
    }
    r_ij.x *= inv_t;
    r_ij.y *= inv_t;
    r_ij.z *= inv_t;

    velocity[id].x = r_ij.x;
    velocity[id].y = r_ij.y;
    velocity[id].z = r_ij.z;
}

void newton_init_vel( )
{
    if( bc == BC_BOX )
        comp_velocities<true ><<<blocks,cuda_block>>>( d_coord, d_scoord, d_velocity, 1.0f / curr_dt);
    else
        comp_velocities<false><<<blocks,cuda_block>>>( d_coord, d_scoord, d_velocity, 1.0f / curr_dt);
    cudaThreadSynchronize( ); CCERR
}

void init_newton( )
{
    INT i;
    
    cudaMalloc( (void**) &d_velocity, sizeof (DOUBLE) * DIMS1 * size ); CCERR
    h_velocity = (DOUBLE*) malloc( sizeof(DOUBLE) * DIMS1 * size); CHMEM(h_velocity);
    cudaMalloc( (void**) &d_rev, blocks ); CCERR
    h_rev = (STR) malloc( blocks ); CHMEM(h_rev);
    cudaMalloc( (void**) &d_times, sizeof (DOUBLE)*size*2 ); CCERR
    h_times = (DOUBLE*) malloc( sizeof (DOUBLE)*size*2 ); CHMEM(h_times);

    for( i = 0; i < size; ++i )
    {
        h_velocity[DIMS1*i + 3] = LJ[i*2];
    }
    cudaMemcpy( d_velocity, h_velocity, sizeof(float)*size*DIMS1, cudaMemcpyHostToDevice ); CCERR

    float bx[3];
    float iL[3];
    bx[0] = box[0], bx[1] = box[1], bx[2] = box[2];
    iL[0] = inv_L[0], iL[1] = inv_L[1], iL[2] = inv_L[2];
    cudaMemcpyToSymbol( "dev_box0", bx, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_box1", bx + 1, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_box2", bx + 2, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_invL0", iL, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_invL1", iL + 1, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_invL2", iL + 2, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
}

void free_newton( )
{
    cudaFree( d_velocity ); CCERR
    free( h_velocity );
    cudaFree( d_rev ); CCERR
    free( h_rev );
    cudaFree( d_times ); CCERR
    free( h_times );
}
