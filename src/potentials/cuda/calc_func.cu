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

#include <cuda_runtime_api.h>
#include <cublas.h>

#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <memory.h>

using namespace std;
extern "C"
{
#include "../calc_func.h"

#include "../../cuda.h"
#include "../../err.h"
#include "../../input.h"
#include "../../math_help.h"
#include "../../trans.h"

#include "../angle.h"
#include "../angle_cos.h"
#include "../bond.h"
#include "../dihe.h"
#include "../dihe_angle.h"
#include "../electro.h"
#include "../electro_ext.h"
#include "../LJ.h"

#include "../bucket.h"
}

#include "../../netdecl.cuh"
#include "sorter.cuh"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

__constant__ int dev_magic_tab[1<<BUCKET_DIM][1<<BUCKET_DIM];
__constant__ float dev_kappa;
__constant__ float dev_gamma_eps;
__constant__ float dev_cutoffc;
__constant__ float dev_cutoffLJ;
__constant__ float dev_alpha_lj;
__constant__ float dev_lj_6_term;
__constant__ float dev_box0_f;
__constant__ float dev_box1_f;
__constant__ float dev_box2_f;
__constant__ float dev_invL0_f;
__constant__ float dev_invL1_f;
__constant__ float dev_invL2_f;
__constant__ float dev_bond_lj_scale;
__constant__ float dev_bond_c_scale;
__constant__ int dev_size;
__constant__ float dev_offset[3];
__constant__ float dev_blen[3];
__constant__ int dev_bonds;

texture< float4, 1, cudaReadModeElementType> tr_coords;
texture< float, 1, cudaReadModeElementType> tr_q;
texture< float2, 1, cudaReadModeElementType> tr_lj;

INT max_gs = 0;
INT* graph = NULL;
INT* d_gconns = NULL;
INT* kkeys = NULL;
DOUBLE* d_out = NULL;
DOUBLE* d_bforces = NULL;
DOUBLE* d_help = NULL;
INT* d_indxs = NULL;
DOUBLE* d_ljcoord;
#if USE_CUDPP
StaticSorter* sorter = NULL;
#endif

extern "C"
{
    DOUBLE* E;

    MAKE_STR_IN( DOUBLE, bond_lj_scale, 1.0, "scaling factor for Lennard–Jones interactions between bonded partners 1–2" )
    MAKE_STR_IN( DOUBLE, bond_c_scale, 1.0, "scaling factor for electrostatic interactions between bonded partners 1–2" )
    MAKE_STR_IN( YESNO, check_overlap, 1, "checking overlap" )

    char* d_overlap = NULL;
    DOUBLE* d_lj = NULL;
    DOUBLE* d_Q = NULL;
}

void mstart( );

double mstop( );

template<bool pbc>
__device__ float4 get_dist( float4 ci, float4 cj )
{
    float4 r_ij;
    r_ij.x = ci.x - cj.x;
    r_ij.y = ci.y - cj.y;
    r_ij.z = ci.z - cj.z;

    if ( pbc )
    {
        r_ij.x -= dev_box0_f * floor( r_ij.x * dev_invL0_f + 0.5f );
        r_ij.y -= dev_box1_f * floor( r_ij.y * dev_invL1_f + 0.5f );
        r_ij.z -= dev_box2_f * floor( r_ij.z * dev_invL2_f + 0.5f );
    }
    return r_ij;
}

template<bool pbc>
__device__ char isoverlap( float4 ci, float4 cj )
{
    float4 r_ij = get_dist<pbc>( ci, cj );

    float r_sq = r_ij.x*r_ij.x + r_ij.y*r_ij.y + r_ij.z*r_ij.z;

    float sigma = (ci.w + cj.w)*(ci.w + cj.w);

    return (char) (r_sq < sigma);
}

template<bool pbc>
__device__ int isoverlap( float4 ci, float x, float y, float z, float w )
{
    float4 cj = {x,y,z,w};
    float4 r_ij = get_dist<pbc>( ci, cj );

    float r_sq = r_ij.x*r_ij.x + r_ij.y*r_ij.y + r_ij.z*r_ij.z;

    float sigma = (ci.w + cj.w)*(ci.w + cj.w);

    return (r_sq < sigma);
}

template<bool pbc>
__global__ void check_beads_overlap( void* ptr, void* ptr_lj, char* overlap )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float4 shPosition[];
    float4* tab = (float4*) ptr;
    float4 curr = tab[id];
    float2* ljtab = (float2*) ptr_lj;
    curr.w = ljtab[id].x;
    char rep = 0;
    for ( int i = 0; i < blockIdx.x; ++i )
    {
        shPosition[ threadIdx.x ] = tab[ i * blockDim.x + threadIdx.x ];
        float2 lw = ljtab[i * blockDim.x + threadIdx.x];
        shPosition[ threadIdx.x ].w = lw.x;
        __syncthreads( );
        for ( int j = 0; j < blockDim.x; ++j )
        {
            rep = rep | isoverlap<pbc > (curr, shPosition[j]);
        }
        __syncthreads( );
    }
    shPosition[ threadIdx.x ] = curr;
    __syncthreads( );
    for ( int j = 0; j < threadIdx.x; ++j )
    {
        rep = rep | isoverlap<pbc > (curr, shPosition[j]);
    }
    if( rep ) overlap[ 0 ] = 1;
}

template<bool pbc, int uncount>
__global__ void check_beads_overlap_tex( char* overlap )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float4 shPosition[];
    float4 curr = tex1Dfetch( tr_coords, id );
    int rep = 0;
    for ( int i = 0; i < blockIdx.x; ++i )
    {
        shPosition[ threadIdx.x ] = tex1Dfetch( tr_coords, i * blockDim.x + threadIdx.x );
        __syncthreads( );
#pragma unroll
        for ( int j = 0; j < uncount; ++j )
        {
            rep = rep | isoverlap<pbc > (curr, shPosition[j]);
        }
        __syncthreads( );
    }
    shPosition[ threadIdx.x ] = curr;
    __syncthreads( );
    for ( int j = 0; j < threadIdx.x; ++j )
    {
        rep = rep | isoverlap<pbc > (curr, shPosition[j]);
    }
    if( rep ) overlap[ 0 ] = 1;
}

__global__ void make_coords( float4* coords, float4* src_coords, float2* src_lj )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    coords[id] = src_coords[id];
    coords[id].w = src_lj[id].x;
}

template<int cuda_block>
__device__ void energy_reduction( float* out, float e_curr )
{
    extern __shared__ float sh[];
    int tid = threadIdx.x;
    volatile float* redE = sh;

    __syncthreads( );
    redE[tid] = e_curr;
    __syncthreads( );
    if ( cuda_block >= 1024 )
    {
        if ( tid < 512 ) redE[tid] += redE[tid + 512];
        __syncthreads( );
    }
    if ( cuda_block >= 512 )
    {
        if ( tid < 256 ) redE[tid] += redE[tid + 256];
        __syncthreads( );
    }
    if ( cuda_block >= 256 )
    {
        if ( tid < 128 ) redE[tid] += redE[tid + 128];
        __syncthreads( );
    }
    if ( cuda_block >= 128 )
    {
        if ( tid < 64 ) redE[tid] += redE[tid + 64];
        __syncthreads( );
    }
    if ( tid < 32 )
    {
        if ( cuda_block >= 64 ) redE[tid] += redE[tid + 32];
        if ( cuda_block >= 32 ) redE[tid] += redE[tid + 16];
        if ( cuda_block >= 16 ) redE[tid] += redE[tid + 8];
        if ( cuda_block >=  8 ) redE[tid] += redE[tid + 4];
        if ( cuda_block >=  4 ) redE[tid] += redE[tid + 2];
        if ( cuda_block >=  2 ) redE[tid] += redE[tid + 1];
    }
    if ( tid == 0 )
    {
        *out = redE[0];
    }
}

template<bool checkOverlap>
__device__ float4 electro( float4 r_ij, float qj, float r_sq, float li, float lj )
{
    float r_norm = sqrtf( r_sq );
    float E = qj * expf( -dev_kappa * r_norm ) / r_norm;
    float f = E * (1 + dev_kappa * r_norm) / r_sq;

    r_ij.x *= f;
    r_ij.y *= f;
    r_ij.z *= f;
    if ( checkOverlap ) r_ij.w = (li + lj) > r_norm ? nanf( "" ) : E;
    else r_ij.w = E;
    return r_ij;
}

template<bool pbc,bool checkOverlap>
__device__ float4 electro_s( float4 ci, float4 cj, float li, float lj )
{
    float4 r_ij = get_dist<pbc>( ci, cj );

    float r_sq = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    if ( r_sq <= dev_cutoffc )
    {
        return electro<checkOverlap>( r_ij, cj.w, r_sq, li, lj );
    }
    float4 ret = {0.0f, 0.0f, 0.0f, 0.0f};
    return ret;
}

template<int cg>
__device__ bool gconnected( int id, int jd, int* gconns )
{
    bool conn = 0;
    if ( cg >= 9 )
    {
		int* gcurr = gconns + id * dev_bonds;
        for ( int i = 0; i < dev_bonds; ++i )
            conn = conn || (jd == gcurr[i]);
    }
    else
    {
		int* gcurr = gconns + id * cg;
        if ( cg >= 1 ) conn = conn || (jd == gcurr[0]);
        if ( cg >= 2 ) conn = conn || (jd == gcurr[1]);
        if ( cg >= 3 ) conn = conn || (jd == gcurr[2]);
        if ( cg >= 4 ) conn = conn || (jd == gcurr[3]);
        if ( cg >= 5 ) conn = conn || (jd == gcurr[4]);
        if ( cg >= 6 ) conn = conn || (jd == gcurr[5]);
        if ( cg >= 7 ) conn = conn || (jd == gcurr[6]);
        if ( cg >= 8 ) conn = conn || (jd == gcurr[7]);
    }
    return conn;
}

template<bool pbc, int cg,bool checkOverlap>
__device__ float4 electro_g( float4 ci, float4 cj, float li, float lj, int id, int jd, int* gconns )
{
    float4 r_ij = get_dist<pbc>( ci, cj );

    float r_sq = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    if ( r_sq <= dev_cutoffc )
    {
        r_ij = electro<checkOverlap>( r_ij, cj.w, r_sq, li, lj );
        bool conn = gconnected<cg>( id, jd, gconns );
        if ( conn )
        {
            r_ij.x *= dev_bond_c_scale;
            r_ij.y *= dev_bond_c_scale;
            r_ij.z *= dev_bond_c_scale;
            r_ij.w *= dev_bond_c_scale;
            if ( dev_bond_c_scale != 0 ) return r_ij;
        }
        else
        {
            return r_ij;
        }
    }
    float4 ret = {0.0f, 0.0f, 0.0f, 0.0f};
    return ret;
}

template<int cuda_block, bool pbc, int cg,bool checkOverlap>
__global__ void comp_electro( float4* tab, float* qtab, float2* ljtab, void* forces, int* gconns )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/
    extern __shared__ float sh[];
    float4* shPosition = (float4*)sh;
    float* shLJ = sh + 4*blockDim.x;
    float3* ftab = (float3*) forces;
    
    float4 curr = tab[id];
    curr.w = qtab[id];
    float currRadiiLJ = ljtab[id].x;
    
    float3 f_curr = {0.0f, 0.0f, 0.0f};
    float e_curr = 0.0f;

    for ( int i = 0; i < gridDim.x; ++i )
    {
        float2 lw = ljtab[i * blockDim.x + threadIdx.x];
        shLJ[ threadIdx.x ] = lw.x;
        shPosition[ threadIdx.x ] = tab[ i * blockDim.x + threadIdx.x ];
        shPosition[ threadIdx.x ].w = qtab[i * blockDim.x + threadIdx.x];
        __syncthreads( );
        for ( int j = 0; j < blockDim.x; ++j )
        {
            float4 f;
            int jd = i * blockDim.x + j;
            if ( cg == 0 )
                f = electro_s<pbc,checkOverlap> (curr, shPosition[j],currRadiiLJ,shLJ[j]);
            else
                f = electro_g<pbc, cg, checkOverlap> (curr, shPosition[j],currRadiiLJ,shLJ[j], id, jd, gconns );
            f_curr.x += (id!=jd) ? f.x : 0.0f;
            f_curr.y += (id!=jd) ? f.y : 0.0f;
            f_curr.z += (id!=jd) ? f.z : 0.0f;
            e_curr +=   (id!=jd) ? f.w : 0.0f;
        }
        __syncthreads( );
    }
    f_curr.x *= curr.w * dev_gamma_eps; /* curr.w = Qi*/
    f_curr.y *= curr.w * dev_gamma_eps;
    f_curr.z *= curr.w * dev_gamma_eps;
    e_curr *= curr.w * dev_gamma_eps;
    ftab[id] = f_curr;
    energy_reduction<cuda_block>( (float*)forces + 3 * blockDim.x * gridDim.x + blockIdx.x, e_curr );
}

typedef struct
{
    float4 a;
    float4 b;
} float8;

template<bool chOverlap>
__device__ float4 LJ_s( float4 r_ij, float r_sq , float2 li, float2 lj, float flag )
{
    float r_norm = sqrtf( r_sq );
    float epsilon_ij = sqrtf( li.y * lj.y );
    float sig_r_6 = (li.x + lj.x) / r_norm;

    sig_r_6 *= sig_r_6;
    sig_r_6 *= sig_r_6 * sig_r_6;

    float El = epsilon_ij * (sig_r_6 * (sig_r_6 - dev_lj_6_term));
    float fl = epsilon_ij * 6 * (sig_r_6 * (2 * sig_r_6 - dev_lj_6_term)) / r_sq;

    r_ij.x = r_ij.x * fl;
    r_ij.y = r_ij.y * fl;
    r_ij.z = r_ij.z * fl;
    if ( chOverlap ) r_ij.w = flag + El;
    else  r_ij.w = El;

    return r_ij;
}

template<bool pbc,bool dyn_cutoff,bool chOverlap>
__device__ float8 electroLJ( float4 ci, float4 cj, float2 li, float2 lj )
{
    float8 r_ij;
    r_ij.a = get_dist<pbc>( ci, cj );
    r_ij.b.x = r_ij.b.y = r_ij.b.z = r_ij.b.w = 0;
    
    float r_sq = r_ij.a.x * r_ij.a.x + r_ij.a.y * r_ij.a.y + r_ij.a.z * r_ij.a.z;
    
    float flag = ((li.x + lj.x)*(li.x + lj.x) > r_sq && chOverlap) ? nanf( "" ) : 0;
    if ( r_sq <= (dyn_cutoff ? (1.259920941444f*(li.x + lj.x)*(li.x + lj.x)) : dev_cutoffLJ) )
    {
        r_ij.b = LJ_s<chOverlap>( r_ij.a, r_sq , li, lj, flag );
    }
    if ( r_sq <= dev_cutoffc )
    {
        r_ij.a = electro<chOverlap>( r_ij.a, cj.w, r_sq, li.x, lj.x );
    }
    else
    {
        r_ij.a.x = r_ij.a.y = r_ij.a.z = r_ij.a.w = 0;
    }
    if (chOverlap) r_ij.a.w += flag;
    return r_ij;
}

template<bool pbc,bool dyn_cutoff,int cg,bool chOverlap>
__device__ float8 electroLJ_g( float4 ci, float4 cj, float2 li, float2 lj, int id, int jd, int* gconns )
{
    float8 r_ij;
    r_ij.a = get_dist<pbc>( ci, cj );
    r_ij.b.x = r_ij.b.y = r_ij.b.z = r_ij.b.w = 0;

    float r_sq = r_ij.a.x * r_ij.a.x + r_ij.a.y * r_ij.a.y + r_ij.a.z * r_ij.a.z;

    float flag = ((li.x + lj.x)*(li.x + lj.x) > r_sq && chOverlap) ? nanf( "" ) : 0;
    int conn = -1;
    if ( r_sq <= (dyn_cutoff ? (1.259920941444f*(li.x + lj.x)*(li.x + lj.x)) : dev_cutoffLJ) )
    {
        r_ij.b = LJ_s<chOverlap>( r_ij.a, r_sq , li, lj, flag );
        conn = gconnected<cg>( id, jd, gconns );
        if ( conn )
        {
            r_ij.b.x *= dev_bond_lj_scale;
            r_ij.b.y *= dev_bond_lj_scale;
            r_ij.b.z *= dev_bond_lj_scale;
            r_ij.b.w *= dev_bond_lj_scale;
            if ( dev_bond_lj_scale == 0 )
                r_ij.b.x = r_ij.b.y = r_ij.b.z = r_ij.b.w = 0;
        }
    }
    if ( r_sq <= dev_cutoffc )
    {
        r_ij.a = electro<chOverlap>( r_ij.a, cj.w, r_sq, li.x, lj.x );
        if (conn==-1) conn = gconnected<cg>( id, jd, gconns );
        if ( conn )
        {
            r_ij.a.x *= dev_bond_c_scale;
            r_ij.a.y *= dev_bond_c_scale;
            r_ij.a.z *= dev_bond_c_scale;
            r_ij.a.w *= dev_bond_c_scale;
            if ( dev_bond_c_scale == 0 )
                r_ij.a.x = r_ij.a.y = r_ij.a.z = r_ij.a.w = 0;
        }
    }
    else
    {
        r_ij.a.x = r_ij.a.y = r_ij.a.z = r_ij.a.w = 0;
    }
    if (chOverlap) r_ij.a.w += flag;
    return r_ij;
}

template<int cuda_block, bool pbc, bool dyn_cutoff, int cg, bool chOverlap>
__global__ void comp_electroLJ( float4* tab, float* qtab, float2* ljtab, void* forces, int* gconns )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/
    extern __shared__ float sh[];
    float4* shPosition = (float4*)sh;
    float2* shLJ = (float2*)(sh + 4*blockDim.x);
    float3* ftab = (float3*) forces;

    float4 curr = tab[id];
    curr.w = qtab[id];
    float2 currLJ = ljtab[id];

    float8 f_curr = {0.0f, 0.0f, 0.0f, 0.0f,0.0f, 0.0f, 0.0f, 0.0f};
    
    for ( int i = 0; i < gridDim.x; ++i )
    {
        shLJ[ threadIdx.x ] = ljtab[ i * blockDim.x + threadIdx.x];
        shPosition[ threadIdx.x ] = tab[ i * blockDim.x + threadIdx.x ];
        shPosition[ threadIdx.x ].w = qtab[ i * blockDim.x + threadIdx.x];
        __syncthreads( );
        for ( int j = 0; j < blockDim.x; ++j )
        {
            float8 f;
            int jd = i * blockDim.x + j;
            if ( cg == 0 ) f = electroLJ<pbc,dyn_cutoff,chOverlap> (curr, shPosition[j],currLJ,shLJ[j]);
            else f = electroLJ_g<pbc,dyn_cutoff,cg,chOverlap> (curr, shPosition[j],currLJ,shLJ[j],id,jd,gconns);
            f_curr.a.x += (id!=jd) ? f.a.x : 0.0f;
            f_curr.a.y += (id!=jd) ? f.a.y : 0.0f;
            f_curr.a.z += (id!=jd) ? f.a.z : 0.0f;
            f_curr.a.w += (id!=jd) ? f.a.w : 0.0f;

            f_curr.b.x += (id!=jd) ? f.b.x : 0.0f;
            f_curr.b.y += (id!=jd) ? f.b.y : 0.0f;
            f_curr.b.z += (id!=jd) ? f.b.z : 0.0f;
            f_curr.b.w += (id!=jd) ? f.b.w : 0.0f;
        }
        __syncthreads( );
    }
    f_curr.a.x *= curr.w * dev_gamma_eps; /* curr.w = Qi*/
    f_curr.a.y *= curr.w * dev_gamma_eps;
    f_curr.a.z *= curr.w * dev_gamma_eps;
    f_curr.a.w *= curr.w * dev_gamma_eps;

    f_curr.b.x *= dev_alpha_lj;
    f_curr.b.y *= dev_alpha_lj;
    f_curr.b.z *= dev_alpha_lj;
    f_curr.b.w *= dev_alpha_lj;

    float e_cu = f_curr.a.w;
    ftab[id].x = f_curr.a.x + f_curr.b.x;
    ftab[id].y = f_curr.a.y + f_curr.b.y;
    ftab[id].z = f_curr.a.z + f_curr.b.z;
    
    energy_reduction<cuda_block>( (float*)forces + 3 * blockDim.x * gridDim.x + blockIdx.x, e_cu );
    __syncthreads();
    e_cu = f_curr.b.w;
    energy_reduction<cuda_block>( (float*)forces + 3 * blockDim.x * gridDim.x + blockIdx.x + gridDim.x, e_cu );
}

template<bool pbc, bool dyn_cutoff,bool chOverlap>
__device__ float4 sLJ( float4 ci, float4 cj, float eps1, float eps2 )
{
    float4 r_ij = get_dist<pbc>( ci, cj );

    float r_sq = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    float flag = ((ci.w + cj.w)*(ci.w + cj.w) > r_sq && chOverlap) ? nanf( "" ) : 0;
    
    if ( r_sq <= (dyn_cutoff ? (1.259920941444f * (ci.w + cj.w)*(ci.w + cj.w)) : dev_cutoffLJ ) )
    {
        float2 li = {ci.w, eps1};
        float2 lj = {cj.w, eps2};
        return LJ_s<chOverlap>( r_ij, r_sq , li, lj, flag );
    }
    float4 ret = {0.0f, 0.0f, 0.0f, flag};
    return ret;
}

template<bool pbc, bool dyn_cutoff, int cg, bool chOverlap>
__device__ float4 sLJ_g( float4 ci, float4 cj, float eps1, float eps2, int id, int jd, int* gconns )
{
    float4 r_ij = get_dist<pbc>( ci, cj );

    float r_sq = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    float flag = ((ci.w + cj.w)*(ci.w + cj.w) > r_sq && chOverlap) ? nanf( "" ) : 0;

    if ( r_sq <= (dyn_cutoff ? (1.259920941444f * (ci.w + cj.w)*(ci.w + cj.w)) : dev_cutoffLJ ) )
    {
        float2 li = {ci.w, eps1};
        float2 lj = {cj.w, eps2};
        r_ij = LJ_s<chOverlap>( r_ij, r_sq , li, lj, flag );
        bool conn = gconnected<cg>( id, jd, gconns );
        if ( conn )
        {
            r_ij.x *= dev_bond_lj_scale;
            r_ij.y *= dev_bond_lj_scale;
            r_ij.z *= dev_bond_lj_scale;
            r_ij.w *= dev_bond_lj_scale;
            if ( dev_bond_lj_scale != 0 ) return r_ij;
        }
        else
        {
            return r_ij;
        }
    }
    float4 ret = {0.0f, 0.0f, 0.0f, flag};
    return ret;
}

template<int cuda_block, bool pbc, bool dyn_cutoff, int cg, bool chOverlap>
__global__ void comp_LJ( void* ptr, void* ptr_lj, void* forces, int* gconns )
{
    int id = blockIdx.x * blockDim.x + threadIdx.x; /*id*/
    extern __shared__ float sh[];
    float4* shPosition = (float4*) sh;
    float* shEps = sh + 4 * blockDim.x;
    float4* tab = (float4*) ptr;
    float3* ftab = (float3*) forces;
    float4 curr = tab[id];
    float3 f_curr = {0.0f, 0.0f, 0.0f};
    float2* vLJ = (float2*) ptr_lj;
    float eps_curr = 0;
    float e_curr = 0;
    float2 currLJ = vLJ[id];
    curr.w = currLJ.x;
    eps_curr = currLJ.y;
    for ( int i = 0; i < gridDim.x; ++i )
    {
        float4 posi = tab[ i * blockDim.x + threadIdx.x ];
        float2 LJi = vLJ[ i * blockDim.x + threadIdx.x ];
        posi.w = LJi.x;
        shPosition[ threadIdx.x ] = posi;
        shEps[ threadIdx.x ] = LJi.y;
        __syncthreads( );
        for ( int j = 0; j < blockDim.x; ++j )
        {
            float4 f;
            int jd = i * blockDim.x + j;
            if ( cg == 0 ) f = sLJ<pbc, dyn_cutoff,chOverlap> ( curr, shPosition[j], shEps[j], eps_curr );
            else f = sLJ_g<pbc, dyn_cutoff,cg,chOverlap> ( curr, shPosition[j], shEps[j], eps_curr, id, jd, gconns );
            f_curr.x += (id!=jd) ? f.x : 0.0f;
            f_curr.y += (id!=jd) ? f.y : 0.0f;
            f_curr.z += (id!=jd) ? f.z : 0.0f;
            e_curr += (id!=jd) ? f.w : 0.0f;
        }
        __syncthreads( );
    }
    f_curr.x *= dev_alpha_lj;
    f_curr.y *= dev_alpha_lj;
    f_curr.z *= dev_alpha_lj;
    e_curr   *= dev_alpha_lj;

    ftab[id] = f_curr;

    energy_reduction<cuda_block>( (float*)forces + 3 * blockDim.x * gridDim.x + blockIdx.x, e_curr );
}

INT check_E( DOUBLE* _E, int energy_type, int iters, DOUBLE* es = h_forces )
{
    INT rep = 0;
    INT i;
    if ( _E )
    {
        _E[ energy_type ] = 0;
        for ( i = 0; i < iters; ++i )
        {
            _E[ energy_type ] += es[ i ];
        }
        _E[ energy_type ] *= 0.5;
        rep = isnan( _E[ energy_type ] );
    }
    else
    {
        for ( i = 0; i < iters; ++i )
        {
            rep |= isnan( es[ i ] );
        }
    }
    return rep;
}

template<bool pbc, int cg, bool checkOverlap>
void run_elec_bs()
{
    switch(cuda_block)
    {
        case(1024): comp_electro<1024,pbc,cg,checkOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces, d_gconns); break;
        case( 512): comp_electro< 512,pbc,cg,checkOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces, d_gconns); break;
        case( 256): comp_electro< 256,pbc,cg,checkOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces, d_gconns); break;
        case( 128): comp_electro< 128,pbc,cg,checkOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces, d_gconns); break;
        case(  64): comp_electro<  64,pbc,cg,checkOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces, d_gconns); break;
        case(  32): comp_electro<  32,pbc,cg,checkOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces, d_gconns); break;
        default: UNERR( "Unimplemented");
    }
}

template<bool pbc,bool checkOverlap>
void run_elec_g()
{
	if ( max_gs > 8 )
	{
		run_elec_bs<pbc,9,checkOverlap>();
	}
	else
	{
		switch(max_gs)
		{
			case(8): run_elec_bs<pbc,8,checkOverlap>(); break;
			case(7): run_elec_bs<pbc,7,checkOverlap>(); break;
			case(6): run_elec_bs<pbc,6,checkOverlap>(); break;
			case(5): run_elec_bs<pbc,5,checkOverlap>(); break;
			case(4): run_elec_bs<pbc,4,checkOverlap>(); break;
			case(3): run_elec_bs<pbc,3,checkOverlap>(); break;
			case(2): run_elec_bs<pbc,2,checkOverlap>(); break;
			case(1): run_elec_bs<pbc,1,checkOverlap>(); break;
			case(0): run_elec_bs<pbc,0,checkOverlap>(); break;
			default: UNERR("Unimplemented");
		}
	}
}

template<bool checkOverlap>
INT calc_elec( DOUBLE* _E )
{
    if ( bc == BC_BOX )
    {
        run_elec_g<true,checkOverlap>();
    }
    else
    {
        run_elec_g<false,checkOverlap>();
    }
    cudaThreadSynchronize( ); CCERR

    cudaMemcpy( h_forces, d_forces + DIMS0*size, blocks * sizeof (float), cudaMemcpyDeviceToHost ); CCERR

    return check_E( _E, ENERGY_COULOMB, blocks );
}

template <bool pbc, bool dyn_cutoff,int cg,bool chOverlap>
void run_LJ_g()
{
    switch(cuda_block)
    {
        case(1024): comp_LJ <1024, pbc, dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>(d_coord, d_lj, d_forces,d_gconns);break;
        case( 512): comp_LJ < 512, pbc, dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>(d_coord, d_lj, d_forces,d_gconns);break;
        case( 256): comp_LJ < 256, pbc, dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>(d_coord, d_lj, d_forces,d_gconns);break;
        case( 128): comp_LJ < 128, pbc, dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>(d_coord, d_lj, d_forces,d_gconns);break;
        case(  64): comp_LJ <  64, pbc, dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>(d_coord, d_lj, d_forces,d_gconns);break;
        case(  32): comp_LJ <  32, pbc, dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 5, 0 >>>(d_coord, d_lj, d_forces,d_gconns);break;
        default: UNERR("Unimplemented");
    }
}

template <bool pbc, bool dyn_cutoff,bool chOverlap>
void run_LJ()
{
	if ( max_gs > 8 )
	{
		run_LJ_g<pbc,dyn_cutoff,9,chOverlap>();
	}
	else
	{
		switch(max_gs)
		{
			case(8): run_LJ_g<pbc,dyn_cutoff,8,chOverlap>(); break;
			case(7): run_LJ_g<pbc,dyn_cutoff,7,chOverlap>(); break;
			case(6): run_LJ_g<pbc,dyn_cutoff,6,chOverlap>(); break;
			case(5): run_LJ_g<pbc,dyn_cutoff,5,chOverlap>(); break;
			case(4): run_LJ_g<pbc,dyn_cutoff,4,chOverlap>(); break;
			case(3): run_LJ_g<pbc,dyn_cutoff,3,chOverlap>(); break;
			case(2): run_LJ_g<pbc,dyn_cutoff,2,chOverlap>(); break;
			case(1): run_LJ_g<pbc,dyn_cutoff,1,chOverlap>(); break;
			case(0): run_LJ_g<pbc,dyn_cutoff,0,chOverlap>(); break;
			default: UNERR("Unimplemented");
		}
	}
}

template<bool chOverlap>
INT calc_LJ( DOUBLE* _E )
{
    if ( bc == BC_BOX )
    {
        if ( cutoff_lj <= 0.0f )
        {
            run_LJ<true,true,chOverlap>();
        }
        else
        {
            run_LJ<true,false,chOverlap>();
        }
    }
    else
    {
        if ( cutoff_lj <= 0.0f )
        {
            run_LJ<false,true,chOverlap>();
        }
        else
        {
            run_LJ<false,false,chOverlap>();
        }
    }
    cudaThreadSynchronize( ); CCERR

    cudaMemcpy( h_forces, d_forces + DIMS0*size, blocks * sizeof (float), cudaMemcpyDeviceToHost ); CCERR

    return check_E( _E, ENERGY_LJ, blocks );
}

template <bool pbc, bool dyn_cutoff,int cg,bool chOverlap>
void run_electro_LJ_g()
{
    switch(cuda_block)
    {
        case(1024): comp_electroLJ<512,pbc,dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 6, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces,d_gconns); break;
        case( 512): comp_electroLJ<512,pbc,dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 6, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces,d_gconns); break;
        case( 256): comp_electroLJ<256,pbc,dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 6, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces,d_gconns); break;
        case( 128): comp_electroLJ<128,pbc,dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 6, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces,d_gconns); break;
        case(  64): comp_electroLJ< 64,pbc,dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 6, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces,d_gconns); break;
        case(  32): comp_electroLJ< 32,pbc,dyn_cutoff,cg,chOverlap> <<<blocks, cuda_block, cuda_block * 4 * 6, 0 >>>((float4*)d_coord, d_Q, (float2*)d_lj, d_forces,d_gconns); break;
        default: UNERR("Unimplemented");
    }
}

template <bool pbc, bool dyn_cutoff, bool chOverlap>
void run_electro_LJ()
{
	if ( max_gs > 8 )
	{
		run_electro_LJ_g<pbc,dyn_cutoff,9,chOverlap>();
	}
	else
	{
		switch(max_gs)
		{
			case(8): run_electro_LJ_g<pbc,dyn_cutoff,8,chOverlap>(); break;
			case(7): run_electro_LJ_g<pbc,dyn_cutoff,7,chOverlap>(); break;
			case(6): run_electro_LJ_g<pbc,dyn_cutoff,6,chOverlap>(); break;
			case(5): run_electro_LJ_g<pbc,dyn_cutoff,5,chOverlap>(); break;
			case(4): run_electro_LJ_g<pbc,dyn_cutoff,4,chOverlap>(); break;
			case(3): run_electro_LJ_g<pbc,dyn_cutoff,3,chOverlap>(); break;
			case(2): run_electro_LJ_g<pbc,dyn_cutoff,2,chOverlap>(); break;
			case(1): run_electro_LJ_g<pbc,dyn_cutoff,1,chOverlap>(); break;
			case(0): run_electro_LJ_g<pbc,dyn_cutoff,0,chOverlap>(); break;
			default: UNERR("Unimplemented");
		}
	}
}

template<bool chOverlap>
INT calc_electro_LJ( DOUBLE* _E )
{
    if ( bc == BC_BOX )
    {
        if ( cutoff_lj <= 0.0f )
        {
            run_electro_LJ<true,true,chOverlap>();
        }
        else
        {
            run_electro_LJ<true,false,chOverlap>();
        }
    }
    else
    {
        if ( cutoff_lj <= 0.0f )
        {
            run_electro_LJ<false,true,chOverlap>();
        }
        else
        {
            run_electro_LJ<false,false,chOverlap>();
        }
    }
    cudaThreadSynchronize( ); CCERR

    cudaMemcpy( h_forces, d_forces + DIMS0*size, 2*blocks*sizeof(float), cudaMemcpyDeviceToHost ); CCERR

    INT rep = check_E( _E, ENERGY_COULOMB, blocks );
    rep |= check_E( _E, ENERGY_LJ, blocks, h_forces + blocks );
    return rep;
}

INT calc_overlap()
{
    INT rep;
    rep = 0;
    cudaMemset( d_overlap, 0, sizeof (char) ); CCERR

    if ( bc == BC_BOX )
    {
        check_beads_overlap<true> <<<blocks, cuda_block, cuda_block * 4 * 4 >>>(d_coord, d_lj, d_overlap);
    }
    else
    {
        make_coords<<<blocks,cuda_block>>>( (float4*)d_ljcoord, (float4*)d_coord, (float2*)d_lj );
        cudaThreadSynchronize( ); CCERR
        switch( cuda_block )
        {
            case(512): check_beads_overlap_tex<false,512><<<blocks, cuda_block, cuda_block * 4 * 4 >>>( d_overlap ); break;
            case(256): check_beads_overlap_tex<false,256><<<blocks, cuda_block, cuda_block * 4 * 4 >>>( d_overlap ); break;
            case(128): check_beads_overlap_tex<false,128><<<blocks, cuda_block, cuda_block * 4 * 4 >>>( d_overlap ); break;
            case( 64): check_beads_overlap_tex<false, 64><<<blocks, cuda_block, cuda_block * 4 * 4 >>>( d_overlap ); break;
            case( 32): check_beads_overlap_tex<false, 32><<<blocks, cuda_block, cuda_block * 4 * 4 >>>( d_overlap ); break;
            default: UNERR("Wrong cuda_block value");
        }
    }
    cudaThreadSynchronize( ); CCERR
    char overlap;
    cudaMemcpy( &overlap, d_overlap, sizeof (char), cudaMemcpyDeviceToHost ); CCERR
    rep = rep || (int)(overlap);
    return rep;
}

template<bool chOverlap>
INT brute_force( DOUBLE* _E )
{
    INT rep = 0;
    if ( elec || alpha_lj != 0.0 ) /* if elec or lj turn on*/
    {
        if ( elec && alpha_lj == 0.0 )/* only electro static */
        {
            rep = calc_elec<chOverlap>( _E );
        }
        else if ( !elec && alpha_lj != 0.0 )/* only LJ */
        {
            rep = calc_LJ<chOverlap>( _E );
        }
        else /* both */
        {
            rep = calc_electro_LJ<chOverlap>( _E );
        }
    }
    else
    {
        if ( chOverlap ) rep = calc_overlap();
    }
    return rep;
}

__device__ void red_single( int tid, int curr_size, float* sh )
{
    sh[tid             ] = fmax( sh[tid             ], sh[tid+curr_size             ] );
    sh[tid+  blockDim.x] = fmax( sh[tid+  blockDim.x], sh[tid+curr_size+  blockDim.x] );
    sh[tid+2*blockDim.x] = fmax( sh[tid+2*blockDim.x], sh[tid+curr_size+2*blockDim.x] );
    sh[tid+3*blockDim.x] = fmin( sh[tid+3*blockDim.x], sh[tid+curr_size+3*blockDim.x] );
    sh[tid+4*blockDim.x] = fmin( sh[tid+4*blockDim.x], sh[tid+curr_size+4*blockDim.x] );
    sh[tid+5*blockDim.x] = fmin( sh[tid+5*blockDim.x], sh[tid+curr_size+5*blockDim.x] );
}

template<int cuda_block,int curr_size>
__device__ void red6( int tid )
{
    extern __shared__ float sh[];
    if ( cuda_block >= 2*curr_size )
    {
        if ( tid < curr_size )
        {
            red_single( tid, curr_size, sh );
        }
    }
}

template<int cuda_block>
__global__ void comp_h_l( void* coord, float* fout )
{
    int tid = threadIdx.x;
    int id = blockDim.x*blockIdx.x + tid;
    float4* cds = (float4*)coord;
    extern __shared__ float sh[];
    float hx, hy, hz, lx, ly, lz;
    hx = hy = hz = -FLT_MAX;
    lx = ly = lz = FLT_MAX;
    while ( id < dev_size )
    {
        float x, y, z;
        x = cds[id].x;
        y = cds[id].y;
        z = cds[id].z;
        hx = fmax(x,hx); hy = fmax(y,hy); hz = fmax(z,hz);
        lx = fmin(x,lx); ly = fmin(y,ly); lz = fmin(z,lz);
        id += blockDim.x*gridDim.x;
    }
    sh[tid             ] = hx;
    sh[tid+  blockDim.x] = hy;
    sh[tid+2*blockDim.x] = hz;
    sh[tid+3*blockDim.x] = lx;
    sh[tid+4*blockDim.x] = ly;
    sh[tid+5*blockDim.x] = lz;

    __syncthreads();
    red6<cuda_block,512>( tid ); __syncthreads();
    red6<cuda_block,256>( tid ); __syncthreads();
    red6<cuda_block,128>( tid ); __syncthreads();
    red6<cuda_block, 64>( tid ); __syncthreads();

    if( tid < 32 )
    {
        if( cuda_block >= 64 ) red_single( tid, 32, sh );
        if( cuda_block >= 32 ) red_single( tid, 16, sh );
        if( cuda_block >= 16 ) red_single( tid,  8, sh );
        if( cuda_block >=  8 ) red_single( tid,  4, sh );
        if( cuda_block >=  4 ) red_single( tid,  2, sh );
        if( cuda_block >=  2 ) red_single( tid,  1, sh );
    }
    if ( tid == 0 )
    {
#pragma unroll 6
        for ( int i = 0; i < 6; ++i )
            fout[6*blockIdx.x+i]=sh[i*blockDim.x];
    }
}

template<int last>
__device__ INT correct_indx( int indx )
{
    return ( indx < 0 ) ? (last + indx) : (indx)%last;
}

__device__
int make_key( int x, int y, int z )
{
    int cx = correct_indx<BUCKET_SPATIAL_M + 1>( x );
    int cy = correct_indx<BUCKET_SPATIAL_M + 1>( y );
    int cz = correct_indx<BUCKET_SPATIAL_M + 1>( z );
    int key = ( cx << BUCKET_SPATIAL_X) |
              ( cy << BUCKET_SPATIAL_Y) |
              ( cz << BUCKET_SPATIAL_Z);
    return key;
}

template<int i>
__device__ void add_key( int* kkeys, int key, int id, int kx, int ky, int kz )
{
    key = make_key( kx+(i>>2), ky+((i>>1)&1), kz+(i&1) );
    kkeys[ 8*id + i ] = key;
    kkeys[ 8*id + i + dev_size*8 ] = (id*8) + i;
}

__global__ void compute_keys_value( float4* coords, int* kkeys )
{
    int id = blockDim.x*blockIdx.x + threadIdx.x;

    if ( id < dev_size )
    {
        int kx, ky, kz, key = 0;
        kx = (INT)floor( (coords[id].x + dev_offset[0]) / dev_blen[0] );
        ky = (INT)floor( (coords[id].y + dev_offset[1]) / dev_blen[1] );
        kz = (INT)floor( (coords[id].z + dev_offset[2]) / dev_blen[2] );
        add_key<0>( kkeys, key, id, kx, ky, kz);
        add_key<1>( kkeys, key, id, kx, ky, kz);
        add_key<2>( kkeys, key, id, kx, ky, kz);
        add_key<3>( kkeys, key, id, kx, ky, kz);
        add_key<4>( kkeys, key, id, kx, ky, kz);
        add_key<5>( kkeys, key, id, kx, ky, kz);
        add_key<6>( kkeys, key, id, kx, ky, kz);
        add_key<7>( kkeys, key, id, kx, ky, kz);
    }
}

#define MASK (~15)
#define COUNT 16

template<int op>
__device__ int check_s_overlap( float4 curr, volatile float4* cop, int true_id, volatile int* iop)
{
    int res = isoverlap<false>( curr, cop[ op ].x, cop[ op ].y
                                    , cop[ op ].z, cop[ op ].w);
    res = (true_id != iop[ op ]) * res;
    return res;
}

__global__ void checking_overlap( float4* coords, float2* lj, int* kkeys, char* out )
{
    int tid = threadIdx.x;
    int id = blockDim.x*blockIdx.x + tid;
    int mid = tid&(COUNT-1);
    extern __shared__ float sh[];
    volatile float4* sh_coord = (float4*)sh  + (tid&MASK);
    volatile int* sh_ids = (int*)sh + 4*blockDim.x  + (tid&MASK);
    int last_warp_key = kkeys[ (id & MASK) + COUNT - 1 ];//mask of last in wrap
    int curr_id = kkeys[id+dev_size*8]>>3;
    int rep = 0;
    float4 curr = coords[curr_id];
    curr.w = lj[curr_id].x;
    sh_coord[mid].x = curr.x;
    sh_coord[mid].y = curr.y;
    sh_coord[mid].z = curr.z;
    sh_coord[mid].w = curr.w;
    sh_ids[mid] = curr_id;

    while(1)
    {
        rep = rep || check_s_overlap< 0>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 1>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 2>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 3>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 4>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 5>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 6>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 7>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 8>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap< 9>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap<10>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap<11>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap<12>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap<13>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap<14>( curr, sh_coord, curr_id, sh_ids );
        rep = rep || check_s_overlap<15>( curr, sh_coord, curr_id, sh_ids );

        id += COUNT;
        if ( id < dev_size*8 )
        {
            int ckey = kkeys[ id & MASK ];
            int  cid = kkeys[ id+dev_size*8 ]>>3;
            sh_coord[mid].x = coords[cid].x;
            sh_coord[mid].y = coords[cid].y;
            sh_coord[mid].z = coords[cid].z;
            sh_coord[mid].w = lj[cid].x;
            sh_ids[mid] = cid;

            if ( ckey != last_warp_key ) break;
        }
        else
        {
            break;
        }
    }
    if ( rep ) *out = rep;
}

template<bool checkOverlap>
__device__ float4 compute_electro( float4 curr, float c_lj_r , int c_key,
                                   float4* cop, float* lj_r,   int* kop,
                                   int true_id, int* iop,
                                   int op, int p1 )
{
    float4 res = {0.0f,0.0f,0.0f,0.0f};
    int p2 = iop[ op ] & 7;
    int pred_m = dev_magic_tab[p1][p2];
    int pred = (true_id != (iop[ op ]>>3)) && (c_key == kop[op]) && pred_m;
    if(pred)
    {
        float4 r_op = { cop[ op ].x, cop[ op ].y, cop[ op ].z, cop[ op ].w };
        res = electro_s<false,checkOverlap>( curr, r_op, c_lj_r, lj_r[op] );
    }
    return res;
}

__device__ float4 add4(float4 a, float4 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
    return a;
}

__device__ float8 add8(float8 a, float8 b)
{
    a.a = add4( a.a, b.a );
    a.b = add4( a.b, b.b );
    return a;
}

template<bool checkOverlap>
__global__ void compute_electro_spatial_indxs( float* forces, int* kkeys, int* indxs )
{
    int tid = threadIdx.x;
    int mid = tid & (COUNT-1);
    int id = blockDim.x*blockIdx.x + tid;
    extern __shared__ float sh[];
    float4* sh_coord = ((float4*)sh) + (tid&MASK);//begin of wrap memory
    float* sh_lj_r = sh + 4*blockDim.x  + (tid&MASK);
    int* sh_kids = (int*)sh + 5*blockDim.x  + (tid&MASK);
    int* sh_keys = (int*)sh + 6*blockDim.x  + (tid&MASK);
    int curr_key = kkeys[id];
    int curr_kid = kkeys[id+dev_size*8];
    float curr_lj_r = (tex1Dfetch( tr_lj, curr_kid>>3 )).x;
    float4 res = {0.0f,0.0f,0.0f,0.0f};
    float4 curr = tex1Dfetch( tr_coords, curr_kid>>3 );
    curr.w = tex1Dfetch( tr_q, curr_kid>>3 );
    int pos = id/COUNT;
    int beg = indxs[ 2*pos   ];
    int end = indxs[ 2*pos +1];

    for ( int id = beg; id < end; id += COUNT )
    {
        int cid = kkeys[ id+mid+dev_size*8 ];
        float4 sc = tex1Dfetch( tr_coords, cid>>3 );
        float2 slj = tex1Dfetch( tr_lj, cid>>3 );
        sh_coord[mid].x = sc.x;
        sh_coord[mid].y = sc.y;
        sh_coord[mid].z = sc.z;
        sh_coord[mid].w = tex1Dfetch( tr_q, cid>>3 );
        sh_lj_r[mid] = slj.x;
        sh_kids[mid] = cid;
        sh_keys[mid] = kkeys[ id+mid ];

        for( int i = 0; i < COUNT; ++i )
        {
            res = add4(res,compute_electro<checkOverlap>( curr, curr_lj_r, curr_key, sh_coord,
                sh_lj_r, sh_keys, curr_kid>>3, sh_kids, i, curr_kid&7 ) );
        }
    }

    int out_id =  curr_kid &7;
    int curr_id = curr_kid>>3;
    forces[ DIMS0*curr_id + DIMS1*dev_size*out_id   ] = res.x * curr.w * dev_gamma_eps;
    forces[ DIMS0*curr_id + DIMS1*dev_size*out_id +1] = res.y * curr.w * dev_gamma_eps;
    forces[ DIMS0*curr_id + DIMS1*dev_size*out_id +2] = res.z * curr.w * dev_gamma_eps;

    forces[ curr_id + DIMS0*dev_size + DIMS1*dev_size*out_id ] = res.w * curr.w * dev_gamma_eps;
}

template<bool dyn_cutoff,bool chOverlap>
__device__ float4 compute_LJ( float4 curr, float curr_eps, int c_key,
                              float4* cop, float* lj_eps, int* kop,
                              int ckid, int* iop, int op )
{
    float4 res = {0.0f,0.0f,0.0f,0.0f};
    int p1 = ckid &7;
    int p2 = iop[ op ]&7;
    int pred_m = dev_magic_tab[p1][p2];
    int pred = ((ckid>>3) != (iop[ op ]>>3)) && (c_key == kop[op]) && pred_m;
    if(pred)
    {
        float4 r_op = { cop[ op ].x, cop[ op ].y, cop[ op ].z, cop[ op ].w };
        res = sLJ<false,dyn_cutoff,chOverlap>( curr, r_op, curr_eps, lj_eps[op] );
    }
    return res;
}

template<bool dyn_cutoff, bool chOverlap>
__global__ void compute_LJ_spatial_indxs(float* forces, int* kkeys, int* indxs )
{
    int tid = threadIdx.x;
    int mid = tid & (COUNT-1);
    int id = blockDim.x*blockIdx.x + tid;
    extern __shared__ float sh[];
    float4* sh_coord = ((float4*)sh) + (tid&MASK);//begin of wrap memory
    float* sh_lj_eps = sh + 4*blockDim.x  + (tid&MASK);
    int* sh_kids = (int*)sh + 5*blockDim.x  + (tid&MASK);
    int* sh_keys = (int*)sh + 6*blockDim.x  + (tid&MASK);
    int curr_key = kkeys[id];
    int curr_kid = kkeys[id+dev_size*8];
    float4 res = {0.0f,0.0f,0.0f,0.0f};
    float4 curr = tex1Dfetch( tr_coords, curr_kid>>3 );
    float2 curr_slj = tex1Dfetch( tr_lj, curr_kid>>3 );
    float curr_eps = curr_slj.y;
    curr.w = curr_slj.x;
    int pos = id/COUNT;
    int beg = indxs[ 2*pos   ];
    int end = indxs[ 2*pos +1];
    for ( int id = beg; id < end; id += COUNT )
    {
        int  cid = kkeys[ id+mid+dev_size*8 ];
        float4 sc = tex1Dfetch( tr_coords, cid>>3 );
        float2 slj = tex1Dfetch( tr_lj, cid>>3 );
        sh_coord[mid].x = sc.x;
        sh_coord[mid].y = sc.y;
        sh_coord[mid].z = sc.z;
        sh_coord[mid].w = slj.x;
        sh_lj_eps[mid] = slj.y;
        sh_kids[mid] = cid;
        sh_keys[mid] = kkeys[ id+mid ];

        for( int i = 0; i < COUNT; ++i )
        {
            res = add4( res, compute_LJ<dyn_cutoff,chOverlap>( curr, curr_eps,
                                                   curr_key, sh_coord, sh_lj_eps,
                                                   sh_keys, curr_kid, sh_kids, i) );
        }
    }
    id = blockDim.x*blockIdx.x + tid;
    int out_id = curr_kid &7;
    int curr_id = curr_kid>>3;
    forces[ DIMS0*curr_id + DIMS1*dev_size*out_id   ] = res.x * dev_alpha_lj;
    forces[ DIMS0*curr_id + DIMS1*dev_size*out_id +1] = res.y * dev_alpha_lj;
    forces[ DIMS0*curr_id + DIMS1*dev_size*out_id +2] = res.z * dev_alpha_lj;

    forces[ curr_id + DIMS0*dev_size + DIMS1*dev_size*out_id ] = res.w * dev_alpha_lj;
}

template<bool dyn_cutoff,bool chOverlap>
__device__ float8 compute_electroLJ( float4 curr, float2 c_lj , int c_key, 
                                     int ckid, float4* cop, float2* lj,
                                     int* kop, int* iop, int op )
{
    float8 res = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
    float2 r_lj = { lj[op].x, lj[op].y };
    int p1 = ckid &7;
    int p2 = iop[op] &7;
    int pred_m = dev_magic_tab[p1][p2];
    int pred = ((ckid>>3) != (iop[ op ]>>3)) && (c_key == kop[op]) && pred_m;
    if(pred)
    {
        float4 r_op = { cop[ op ].x, cop[ op ].y, cop[ op ].z, cop[ op ].w };
        res = electroLJ<false,dyn_cutoff,chOverlap>( curr, r_op, c_lj, r_lj );
    }
    return res;
}

template<bool dyn_cutoff,bool chOverlap>
__global__ void compute_electroLJ_spatial_indxs( float* forces, int* kkeys,
                                                 int* indxs )
{
    int tid = threadIdx.x;
    int mid = tid & (COUNT-1);
    int id = blockDim.x*blockIdx.x + tid;
    extern __shared__ float sh[];
    float4* sh_coord = ((float4*)sh) + (tid&MASK);//begin of wrap memory
    float2* sh_lj = (float2*)(sh + 4*blockDim.x) + (tid&MASK);
    int* sh_ids = (int*)sh + 6*blockDim.x  + (tid&MASK);
    int* sh_keys = (int*)sh + 7*blockDim.x  + (tid&MASK);
    int curr_key = kkeys[id];
    int curr_kid = kkeys[id+dev_size*8];
    float2 curr_lj = tex1Dfetch( tr_lj, curr_kid>>3 );
    float4 curr = tex1Dfetch( tr_coords, curr_kid>>3 );
    curr.w = tex1Dfetch( tr_q, curr_kid>>3 );
    int pos = id/COUNT;
    int beg = indxs[ 2*pos   ];
    int end = indxs[ 2*pos +1];
    float8 res = {0.0f, 0.0f, 0.0f, 0.0f,0.0f, 0.0f, 0.0f, 0.0f};
    for ( int id = beg; id < end; id += COUNT )
    {
        int  cid = kkeys[ id+mid+dev_size*8 ];
        float4 sc = tex1Dfetch( tr_coords, cid>>3 );
        float2 slj = tex1Dfetch( tr_lj, cid>>3 );
        sh_coord[mid].x = sc.x;
        sh_coord[mid].y = sc.y;
        sh_coord[mid].z = sc.z;
        sh_coord[mid].w = tex1Dfetch( tr_q, cid>>3 );
        sh_lj[mid] = slj;
        sh_ids[mid] = cid;
        sh_keys[mid] = kkeys[ id+mid ];

        for( int i = 0; i < COUNT; ++i )
        {
            res = add8( res,compute_electroLJ<dyn_cutoff,chOverlap>( curr, curr_lj,
                        curr_key, curr_kid, sh_coord, sh_lj, sh_keys, sh_ids, i ) );
        }
    }

    int out_id = curr_kid & 7;
    int curr_id = curr_kid>>3;

    id = DIMS0*curr_id + (DIMS1+1)*dev_size*out_id;
    forces[ id    ] = res.a.x * curr.w * dev_gamma_eps + res.b.x * dev_alpha_lj;
    forces[ id +1 ] = res.a.y * curr.w * dev_gamma_eps + res.b.y * dev_alpha_lj;
    forces[ id +2 ] = res.a.z * curr.w * dev_gamma_eps + res.b.z * dev_alpha_lj;

    id = curr_id + DIMS0*dev_size + (DIMS1+1)*dev_size*out_id;
    forces[ id ] = res.a.w * curr.w * dev_gamma_eps;

    id = curr_id + DIMS1*dev_size + (DIMS1+1)*dev_size*out_id;
    forces[ id ] = res.b.w * dev_alpha_lj;
}

__global__ void get_indx_spatial( int* kkeys, int* indxs )
{
    int tid = threadIdx.x;
    if( (blockDim.x*blockIdx.x + tid)*COUNT < 8*dev_size )
    {
        int id_end = (blockDim.x*blockIdx.x + tid) * COUNT;/*first in warp*/
        int last_warp_key = kkeys[ id_end + COUNT - 1 ];/*last in wrap*/
        while(1)
        {
            id_end += COUNT;
            if ( id_end < dev_size*8 )
            {
                int ckey = kkeys[ id_end ];
                if ( ckey != last_warp_key ) break;
            }
            else
            {
                break;
            }
        }
        int id_beg = (blockDim.x*blockIdx.x + tid) * COUNT;
        int first_warp_key = kkeys[ id_beg ];//first in wrap
        id_beg += COUNT-1; //move to last in warp
        while(1)
        {
            id_beg -= COUNT;
            if ( id_beg >= 0 )
            {
                int ckey = kkeys[ id_beg ];
                if ( ckey != first_warp_key ) break;
            }
            else
            {
                break;
            }
        }
        int id = blockDim.x*blockIdx.x + tid;
        indxs[ id*2   ] = (id_beg & MASK)+COUNT;// >=
        indxs[ id*2 +1] = id_end;// <
    }
}

template<int cuda_block>
__global__ void reduce_forces( float* forces, float* help )
{
    int tid = threadIdx.x;
    int id = blockDim.x*blockIdx.x + tid;
    float3 f = {0.0f, 0.0f, 0.0f};
#pragma unroll
    for ( int i = 0; i < 8; ++i )
    {
        f.x += forces[ DIMS0*id + DIMS1*dev_size*i   ];
        f.y += forces[ DIMS0*id + DIMS1*dev_size*i +1];
        f.z += forces[ DIMS0*id + DIMS1*dev_size*i +2];
    }
    forces[ DIMS0*id   ] = f.x;
    forces[ DIMS0*id +1] = f.y;
    forces[ DIMS0*id +2] = f.z;

    float e_curr = 0.0f;
#pragma unroll
    for ( int i = 0; i < 8; ++i )
    {
        e_curr += forces[ id + DIMS0*dev_size + DIMS1*dev_size*i ];
    }
    energy_reduction<cuda_block>( help + blockIdx.x, e_curr );
}

template<int cuda_block>
__global__ void reduce_forcesB( float* forces, float* help )
{
    int tid = threadIdx.x;
    int id = blockDim.x*blockIdx.x + tid;
    float3 f = {0.0f, 0.0f, 0.0f};
#pragma unroll
    for ( int i = 0; i < 8; ++i )
    {
        f.x += forces[ DIMS0*id + (DIMS1+1)*dev_size*i   ];
        f.y += forces[ DIMS0*id + (DIMS1+1)*dev_size*i +1];
        f.z += forces[ DIMS0*id + (DIMS1+1)*dev_size*i +2];
    }
    forces[ DIMS0*id   ] = f.x;
    forces[ DIMS0*id +1] = f.y;
    forces[ DIMS0*id +2] = f.z;

    float e_curr = 0.0f;
#pragma unroll
    for ( int i = 0; i < 8; ++i )
    {
        e_curr += forces[ id + DIMS0*dev_size + (DIMS1+1)*dev_size*i ];
    }
    energy_reduction<cuda_block>( help + blockIdx.x, e_curr );

    e_curr = 0.0f;
#pragma unroll
    for ( int i = 0; i < 8; ++i )
    {
        e_curr += forces[ id + DIMS1*dev_size + (DIMS1+1)*dev_size*i ];
    }
    energy_reduction<cuda_block>( help + blockIdx.x + gridDim.x, e_curr );
}

template<bool chOverlap>
INT spatial_force( DOUBLE* _E )
{
    INT rep = 0;
    float* hl = new float[sm_red*6];
    INT i;
    switch(cuda_block)
    {
        case(1024): comp_h_l<1024><<<sm_red, cuda_block,cuda_block*4*6>>>( d_coord, d_out ); break;
        case( 512): comp_h_l< 512><<<sm_red, cuda_block,cuda_block*4*6>>>( d_coord, d_out ); break;
        case( 256): comp_h_l< 256><<<sm_red, cuda_block,cuda_block*4*6>>>( d_coord, d_out ); break;
        case( 128): comp_h_l< 128><<<sm_red, cuda_block,cuda_block*4*6>>>( d_coord, d_out ); break;
        case(  64): comp_h_l<  64><<<sm_red, cuda_block,cuda_block*4*6>>>( d_coord, d_out ); break;
        case(  32): comp_h_l<  32><<<sm_red, cuda_block,cuda_block*4*6>>>( d_coord, d_out ); break;
        default: UNERR("Unimplemented");
    }
    cudaThreadSynchronize(); CCERR
    cudaMemcpy( hl, d_out, sm_red*sizeof(float)*6, cudaMemcpyDeviceToHost ); CCERR
    for ( i = 1; i < sm_red; ++i )
    {
        hl[0]=MAX(hl[6*i  ],hl[0]);
        hl[1]=MAX(hl[6*i+1],hl[1]);
        hl[2]=MAX(hl[6*i+2],hl[2]);
        hl[3]=MIN(hl[6*i+3],hl[3]);
        hl[4]=MIN(hl[6*i+4],hl[4]);
        hl[5]=MIN(hl[6*i+5],hl[5]);
    }
    hl[0] = MAXD( (hl[0]-hl[3])/BUCKET_SPATIAL_M, max_force_cutoff );
    hl[1] = MAXD( (hl[1]-hl[4])/BUCKET_SPATIAL_M, max_force_cutoff );
    hl[2] = MAXD( (hl[2]-hl[5])/BUCKET_SPATIAL_M, max_force_cutoff );
    hl[3]=-hl[3];
    hl[4]=-hl[4];
    hl[5]=-hl[5];
    cudaMemcpyToSymbol( "dev_offset", hl+3, sizeof(float)*3 ); CCERR
    cudaMemcpyToSymbol( "dev_blen", hl, sizeof(float)*3 ); CCERR
    compute_keys_value<<<blocks,cuda_block>>>( (float4*) d_coord, kkeys ); cudaThreadSynchronize(); CCERR
#if USE_CUDPP
    sorter->sort();
#endif
    delete[] hl;
    if ( alpha_lj != 0.0 || elec )
    {
        int indx_blocks = (8*size/COUNT + cuda_block - 1)/cuda_block;
        get_indx_spatial<<<indx_blocks,cuda_block>>>( kkeys, d_indxs );
        cudaThreadSynchronize(); CCERR
        if ( elec && alpha_lj == 0.0 )
        {
            compute_electro_spatial_indxs<chOverlap><<<blocks*8,cuda_block,cuda_block*4*7>>>( d_forces, kkeys, d_indxs );
        }
        else if( !elec )
        {
            if ( cutoff_lj <= 0.0f )
            {
                compute_LJ_spatial_indxs<true,chOverlap><<<blocks*8,cuda_block,cuda_block*4*7>>>
                    ( d_forces, kkeys, d_indxs );
            }
            else
            {
                compute_LJ_spatial_indxs<false,chOverlap><<<blocks*8,cuda_block,cuda_block*4*7>>>
                    ( d_forces, kkeys, d_indxs );
            }
        }
        else
        {
            if ( cutoff_lj <= 0.0f )
            {
                compute_electroLJ_spatial_indxs<true,chOverlap><<<blocks*8,cuda_block,cuda_block*4*8>>>
                    ( d_forces, kkeys, d_indxs );
            }
            else
            {
                compute_electroLJ_spatial_indxs<false,chOverlap><<<blocks*8,cuda_block,cuda_block*4*8>>>
                    ( d_forces, kkeys, d_indxs );
            }
        }
        cudaThreadSynchronize(); CCERR
        if ( alpha_lj != 0.0 && elec )
        {
            switch( cuda_block )
            {
                case(1024):reduce_forcesB<1024><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case( 512):reduce_forcesB< 512><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case( 256):reduce_forcesB< 256><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case( 128):reduce_forcesB< 128><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case(  64):reduce_forcesB<  64><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case(  32):reduce_forcesB<  32><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                default: UNERR("Unimplemented");
            }
            cudaThreadSynchronize(); CCERR
            cudaMemcpy( h_forces, d_help, 2*blocks * sizeof (float), cudaMemcpyDeviceToHost ); CCERR
            rep = rep || check_E( _E, ENERGY_COULOMB, blocks );
            rep = rep || check_E( _E,      ENERGY_LJ, blocks, h_forces + blocks );
        }
        else
        {
            switch( cuda_block )
            {
                case(1024):reduce_forces<1024><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case( 512):reduce_forces< 512><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case( 256):reduce_forces< 256><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case( 128):reduce_forces< 128><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case(  64):reduce_forces<  64><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                case(  32):reduce_forces<  32><<<blocks,cuda_block,cuda_block*4>>>( d_forces, d_help );break;
                default: UNERR("Unimplemented");
            }
            cudaThreadSynchronize(); CCERR
            cudaMemcpy( h_forces, d_help, blocks * sizeof (float), cudaMemcpyDeviceToHost ); CCERR
            rep = rep || check_E( _E, elec ? ENERGY_COULOMB: ENERGY_LJ, blocks );
        }
    }
    else
    {
        //checking overlap
        cudaMemset( d_overlap, 0, sizeof(char) ); CCERR

        checking_overlap<<<blocks*8,cuda_block,cuda_block*4*5>>>( (float4*) d_coord, (float2*)d_lj, kkeys, d_overlap);
        cudaThreadSynchronize(); CCERR

        char overlap;
        cudaMemcpy( &overlap, d_overlap, sizeof (char), cudaMemcpyDeviceToHost ); CCERR
        rep = rep || (int)(overlap);
        return rep;
    }
    return rep;
}

extern "C"
INT compute_forces( DOUBLE* _E, DOUBLE curr_time )
{
    INT rep = 0;

    if( d_bforces || bc == BC_SPHERE )
    {
        cudaMemcpy( coord, d_coord, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
        init_iter_forces( _E );
    }
    else
    {
        if ( _E ) memset( _E, 0, sizeof (DOUBLE) * (ENERGY_LAST + 1) );
    }

    if ( nb_list == BUCKET_NONE )
    {
        if ( check_overlap ) rep = rep || brute_force<true>( _E );
        else rep = rep || brute_force<false>( _E );
    }
    else if ( nb_list == BUCKET_SPATIAL && bc != BC_BOX )
    {
#if !USE_CUDPP
        UNERR("Compile with CUDPP");
#endif
        if ( check_overlap ) rep = rep || spatial_force<true>( _E );
        else rep = rep || spatial_force<false>( _E );
    }
    else
    {
        UNERR("Unimplemented");
    }
    if( d_bforces )
    {
        if (!rep)
        {
            rep |= bonded_inter_parallel( _E, curr_time );
            if( !rep )
            {
                INT j;
                cudaMemcpy( d_bforces, F[0], sizeof(float)*size*DIMS0, cudaMemcpyHostToDevice ); CCERR
                CUBLAS_AXPYD( size*DIMS0, 1.0f, d_bforces, 1, d_forces, 1 ); cublasCheckError(__FILE__,__LINE__);
                if( _E )
                {
                    for ( j = 0; j < ENERGY_LAST; ++j )
                    {
                        _E[ j ] += ES[0][j];
                    }
                }
            }
        }
    }
    return rep;
}

extern "C"
void init_forces( )
{
    int fbonded = bond_size || angle_cos_size || angle_size ||
                  dihe_size || dihe_angle_size || E_ext;
    E = (DOUBLE*) malloc( sizeof (DOUBLE) * (ENERGY_LAST + 1) ); CHMEM(E)

    cudaMalloc( (void**) &d_overlap, sizeof (char) ); CCERR
    
    if ( elec || alpha_lj != 0 || fbonded )
    {
        int add = ( nb_list == BUCKET_SPATIAL ) * 7;
        int addelj = add && elec && (alpha_lj != 0);
        h_forces = (float*) malloc( sizeof (float) * size * DIMS1 );
        cudaMalloc( (void**) &d_forces, ( 1 + add ) * size * (DIMS1 + addelj) * sizeof (float) ); CCERR
        cudaMalloc( (void**) &d_sforces, size * DIMS1 * sizeof (float) ); CCERR
        if ( elec )
        {
            float src_kappa = kappa_c;
            float src_gamma_eps = gamma_c / epsilon_c;
            float src_cutoff = cutoff_c*cutoff_c;
            float src_bondQs = bond_c_scale;
            cudaMemcpyToSymbol( "dev_kappa", &src_kappa, sizeof(float) ); CCERR
            cudaMemcpyToSymbol( "dev_gamma_eps", &src_gamma_eps, sizeof(float) ); CCERR
            cudaMemcpyToSymbol( "dev_cutoffc", &src_cutoff, sizeof(float) ); CCERR
            cudaMemcpyToSymbol( "dev_bond_c_scale", &src_bondQs, sizeof(float) ); CCERR
            float* h_Q = (float*) malloc( sizeof(float)*size ); CHMEM(h_Q)
            for( int i = 0; i < size; ++ i )
            {
                h_Q[ i ] = Q[ i ];
            }
            cudaMalloc( (void**)&d_Q,sizeof(float)*size); CCERR
            cudaMemcpy(d_Q,h_Q,sizeof(float)*size,cudaMemcpyHostToDevice); CCERR
            free(h_Q);
        }
        if ( alpha_lj != 0.0 )
        {
            float src_cutoff = cutoff_lj*cutoff_lj;
            float src_alpha_lj = alpha_lj;
            float src_bondLJs = bond_lj_scale;
            cudaMemcpyToSymbol( "dev_cutoffLJ", &src_cutoff, sizeof(float) ); CCERR
            cudaMemcpyToSymbol( "dev_alpha_lj", &src_alpha_lj, sizeof(float) ); CCERR
            float src_6term = (float)lj_6_term;
            cudaMemcpyToSymbol( "dev_lj_6_term", &src_6term, sizeof(float) ); CCERR
            cudaMemcpyToSymbol( "dev_bond_lj_scale", &src_bondLJs, sizeof(float) ); CCERR
        }
    }
    /* for all we check overlap*/
    float* lj_data = (float*) malloc( size*sizeof(float)*2 ); CHMEM(lj_data)
    for ( int i = 0; i < 2 * size; ++i )
    {
        lj_data[i] = LJ[i];
    }
    cudaMalloc( (void**) &d_lj, sizeof(float)*size*2 ); CCERR
    cudaMemcpy( d_lj, lj_data, sizeof(float)*size*2, cudaMemcpyHostToDevice ); CCERR
    free( lj_data );
    
    if ( bc == BC_BOX )
    {
        float bx[3];
        float iL[3];
        bx[0] = box[0], bx[1] = box[1], bx[2] = box[2];
        iL[0] = inv_L[0], iL[1] = inv_L[1], iL[2] = inv_L[2];
        cudaMemcpyToSymbol( "dev_box0_f", bx, sizeof (float) ); CCERR
        cudaMemcpyToSymbol( "dev_box1_f", bx + 1, sizeof (float) ); CCERR
        cudaMemcpyToSymbol( "dev_box2_f", bx + 2, sizeof (float) ); CCERR
        cudaMemcpyToSymbol( "dev_invL0_f", iL, sizeof (float) ); CCERR
        cudaMemcpyToSymbol( "dev_invL1_f", iL + 1, sizeof (float) ); CCERR
        cudaMemcpyToSymbol( "dev_invL2_f", iL + 2, sizeof (float) ); CCERR
    }
    max_gs = 0;
    
    if( bond_size > 0 && (( elec && bond_c_scale != 1 ) || ( alpha_lj != 0.0 && bond_lj_scale != 1 )) )
    {
        graph = get_connection_graph( &max_gs );
	printf("max connections: %d\n", max_gs);
	cudaMemcpyToSymbol( "dev_bonds", &max_gs, sizeof (int) ); CCERR
        cudaMalloc( (void**) &d_gconns, sizeof(int)*max_gs*size ); CCERR
        cudaMemcpy( d_gconns, graph, sizeof(int)*max_gs*size, cudaMemcpyHostToDevice ); CCERR
        if( nb_list != BUCKET_NONE )
        {
            nb_list = BUCKET_NONE;
            warning( "Connections detected, brute force method only", __FILE__, __LINE__ );
        }
    }

    if( nb_list == BUCKET_SPATIAL )
    {
#if USE_CUDPP
        if ( bc != BC_BOX )
        {
            cudaMalloc( (void**) &kkeys, sizeof(int)*2*8*size ); CCERR
            cudaMalloc( (void**) &d_help, sizeof(float)*2*blocks ); CCERR
            cudaThreadSynchronize(); CCERR
            sorter = new StaticSorter( kkeys, kkeys+8*size, 8*size );
            int magic_tab[1<<BUCKET_DIM][1<<BUCKET_DIM];
            memset( magic_tab, 0, sizeof( magic_tab ) );
            for ( int i = 0; i < magic_pos_size; ++i )
            {
                magic_tab[magic_pos[i][0]][magic_pos[i][1]] = 1;
                magic_tab[magic_pos[i][1]][magic_pos[i][0]] = 1;
            }
            cudaMemcpyToSymbol( "dev_magic_tab", magic_tab, sizeof(magic_tab) ); CCERR
            if ( dev_size*8 % COUNT )
            {
                UNERR("Wrong data size");
            }
            cudaMalloc( (void**) &d_indxs, sizeof(int)*2*8*size/COUNT ); CCERR
            
            //bind textures
            tr_coords.addressMode[0] = cudaAddressModeClamp;
            tr_coords.addressMode[1] = cudaAddressModeClamp;
            tr_coords.addressMode[2] = cudaAddressModeClamp;
            tr_coords.filterMode = cudaFilterModePoint;
            tr_coords.normalized = false;
            cudaChannelFormatDesc channel1 = cudaCreateChannelDesc<float4>(); CCERR
            cudaBindTexture( 0, tr_coords, d_coord, channel1, size*DIMS1*sizeof(float) ); CCERR

            if ( d_Q )
            {
                tr_q.addressMode[0] = cudaAddressModeClamp;
                tr_q.addressMode[1] = cudaAddressModeClamp;
                tr_q.addressMode[2] = cudaAddressModeClamp;
                tr_q.filterMode = cudaFilterModePoint;
                tr_q.normalized = false;
                cudaChannelFormatDesc channel2 = cudaCreateChannelDesc<float>(); CCERR
                cudaBindTexture( 0, tr_q, d_Q, channel2, size*sizeof(float) ); CCERR
            }

            if ( d_lj )
            {
                tr_lj.addressMode[0] = cudaAddressModeClamp;
                tr_lj.addressMode[1] = cudaAddressModeClamp;
                tr_lj.addressMode[2] = cudaAddressModeClamp;
                tr_lj.filterMode = cudaFilterModePoint;
                tr_lj.normalized = false;
                cudaChannelFormatDesc channel3 = cudaCreateChannelDesc<float2>(); CCERR
                cudaBindTexture( 0, tr_lj, d_lj, channel3, size*sizeof(float2) ); CCERR
            }
        }
        else
#endif
        {
            UNERR("Unimplemented method");
        }
        cudaMalloc( (void**) &d_out, sm_red*sizeof(float)*6 ); CCERR
    }
    else if ( nb_list != BUCKET_NONE )
    {
        UNERR("Only nb_list=brute and nb_list=bucket_list are supported.");
    }
    else
    {
        cudaMalloc( (void**) &d_ljcoord, size * DIMS1 * sizeof (float) ); CCERR
        tr_coords.addressMode[0] = cudaAddressModeClamp;
        tr_coords.addressMode[1] = cudaAddressModeClamp;
        tr_coords.addressMode[2] = cudaAddressModeClamp;
        tr_coords.filterMode = cudaFilterModePoint;
        tr_coords.normalized = false;
        cudaChannelFormatDesc channel1 = cudaCreateChannelDesc<float4>(); CCERR
        cudaBindTexture( 0, tr_coords, d_ljcoord, channel1, size*DIMS1*sizeof(float) ); CCERR
    }
    init_cutoff();
    if( fbonded )
    {
        init_Fs();
        cudaMalloc( (void**) &d_bforces, sizeof(float)*size*DIMS0 ); CCERR
    }
    int h_size = size;
    cudaMemcpyToSymbol( "dev_size", &h_size, sizeof(int) ); CCERR
}

extern "C"
void free_forces( )
{
    free( E );
    if ( h_forces ) free( h_forces );
    if ( d_forces ) { cudaFree( d_forces ); CCERR }
    if ( d_ljcoord ) { cudaFree( d_ljcoord ); CCERR }
    free( graph );
    if( d_gconns ) { cudaFree( d_gconns ); CCERR }
    if( kkeys ) { cudaFree( kkeys ); CCERR }
    cudaFree( d_overlap ); CCERR
    if( d_bforces )
    {
        free_Fs();
        cudaFree( d_bforces ); CCERR
    }
    if ( d_indxs ) { cudaFree( d_indxs ); CCERR }
#if USE_CUDPP
    if ( sorter ) delete sorter;
#endif
}
