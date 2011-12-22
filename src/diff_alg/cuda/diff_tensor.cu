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

extern "C"
{
#include "../diff_tensor.h"
#include "../../cuda.h"
#include "../../err.h"
#include "../../trans.h"
}

void mstart( );
void mstop( );

__constant__ float dev_preDii;
__constant__ float dev_preDij;
__constant__ float dev_ewald_alpha;
__constant__ float dev_ewald_alpha2;
__constant__ float dev_ewald_alpha3;
__constant__ float dev_M_PI_ALPHA_SQ;
__constant__ float dev_max_sq;
__constant__ float dev_box0_t;
__constant__ float dev_box1_t;
__constant__ float dev_box2_t;
__constant__ float dev_invL0_t;
__constant__ float dev_invL1_t;
__constant__ float dev_invL2_t;


__device__ static void tensor_calc_pos_diag( float4 a, float* D, int LD )
{
    float3 tmp = {0.0f, 0.0f, 0.0f};
    tmp.x = dev_preDii / a.w;
    *((float3*) D) = tmp;
    tmp.y = tmp.x;
    tmp.x = 0;
    *((float3*) (D + LD)) = tmp;
    tmp.z = tmp.y;
    tmp.y = 0;
    *((float3*) (D + 2 * LD)) = tmp;
}

__device__ float3 add3(float3 a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

__device__ float3 sub3(float3 a, float3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

__device__ static void tensor_calc_pos( float4 a, float4 b, float* D, int LD )
{
    float4 r_ij;
    r_ij.x = a.x - b.x;
    r_ij.y = a.y - b.y;
    r_ij.z = a.z - b.z;

    float r_sq = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    float r_norm = sqrtf( r_sq );
    float3 d1,d2,d3;

    float sigma = (a.w + b.w) / 2;
    float pre = (r_norm < a.w + b.w) ? (dev_preDii / sigma) : (dev_preDij / r_norm);
    float sigm_r = (a.w * a.w + b.w * b.w) / r_sq;

    float divr = (r_norm < a.w + b.w) ? (sigma * r_norm * 32) : (r_sq);
    float mult = (r_norm < a.w + b.w) ? (3) : (1 - sigm_r);
    float part1 = (r_norm < a.w + b.w) ? (1 - (9 * r_norm / (32 * sigma))) : (1 + (sigm_r / 3));

    d3.z = d2.y = d1.x = part1;

    d1.x += mult * r_ij.x * r_ij.x / divr;
    d1.x *= pre;

    d2.y += mult * r_ij.y * r_ij.y / divr;
    d2.y *= pre;

    d3.z += mult * r_ij.z * r_ij.z / divr;
    d3.z *= pre;

    d2.x = d1.y = pre * mult * r_ij.x * r_ij.y / divr;
    d3.x = d1.z = pre * mult * r_ij.x * r_ij.z / divr;
    d2.z = d3.y = pre * mult * r_ij.z * r_ij.y / divr;

    *((float3*) D) = d1;
    *((float3*) (D + LD)) = d2;
    *((float3*) (D + 2 * LD)) = d3;
}

__device__ void get_ev2_norms( float x1, float y1, float z1, float3& e, float& inv_r, float& r_norm, float& r_sq )
{
    inv_r = 0.0;
    e.x = x1;
    e.y = y1;
    e.z = z1;

    r_sq = (x1)*(x1) + (y1)*(y1) + (z1)*(z1);
    r_norm = sqrtf( r_sq );
    inv_r = 1 / r_norm;

    e.x *= inv_r;
    e.y *= inv_r;
    e.z *= inv_r;
}

__device__ void comp_OQ( int f, int g, int h, float4 v, float3* O1, float3* O2, float3* Q1, float3* Q2 )
{
    float r_norm, r_sq, inv_r;
    float3 e;

    get_ev2_norms( v.x + h, v.y + f, v.z + g, e, inv_r, r_norm, r_sq );
    if ( r_norm != 0.0f )
    {
        float inv_r3 = 1 / (r_norm * r_sq);
        float sc = erfc( r_norm * dev_ewald_alpha ) * inv_r;
        float sc2 = 2 * dev_ewald_alpha * M1_SQRTPIF * expf( -dev_ewald_alpha2 * r_sq );
        float sc_ek = e.x * (sc + sc2);

        (*O1).x += sc + sc_ek * e.x;
        sc_ek = e.y * (sc + sc2);
        (*O2).y += sc + sc_ek * e.y;
        (*O1).y += e.x * sc_ek;
        sc_ek = e.z * (sc + sc2);
        (*O2).x += sc + sc_ek * e.z;
        (*O1).z += e.x * sc_ek;
        (*O2).z += e.y * sc_ek;

        sc = inv_r3 * (erfc( r_norm * dev_ewald_alpha ) + 2 * dev_ewald_alpha * r_norm * M1_SQRTPIF * expf( -dev_ewald_alpha2 * r_sq ));
        sc2 = -4 * dev_ewald_alpha3 * M1_SQRTPIF * expf( -dev_ewald_alpha2 * r_sq );
        sc_ek = e.x * (-3 * sc + sc2);
        (*Q1).x += sc + sc_ek * e.x;
        sc_ek = e.y * (-3 * sc + sc2);
        (*Q2).y += sc + sc_ek * e.y;
        (*Q1).y += e.x * sc_ek;
        sc_ek = e.z * (-3 * sc + sc2);
        (*Q2).x += sc + sc_ek * e.z;
        (*Q1).z += e.x * sc_ek;
        (*Q2).z += e.y * sc_ek;
    }

    get_ev2_norms( h, f, g, e, inv_r, r_norm, r_sq );
    if ( r_norm != 0.0f )
    {
        float sc = 2 * expf( dev_M_PI_ALPHA_SQ * r_sq ) / (M_PIF * r_sq);
        float sc2 = -1 + dev_M_PI_ALPHA_SQ * r_sq;
        sc *= cos( 2 * M_PIF * (v.x * h + v.y * f + v.z * g) );
        sc2 *= sc;
        float sc_ek = sc2 * e.x;
        (*O1).x += sc + sc_ek * e.x;
        sc_ek = e.y * sc2;
        (*O2).y += sc + sc_ek * e.y;
        (*O1).y += e.x * sc_ek;
        sc_ek = e.z * sc2;
        (*O2).x += sc + sc_ek * e.z;
        (*O1).z += e.x * sc_ek;
        (*O2).z += e.y * sc_ek;

        sc = 4 * M_PIF * expf( dev_M_PI_ALPHA_SQ * r_sq );
        sc *= cosf( 2 * M_PIF * (v.x * h + v.y * f + v.z * g) );
        sc_ek = e.x * sc;
        (*Q1).x += sc_ek * e.x;
        sc_ek = e.y * sc;
        (*Q2).y += sc_ek * e.y;
        (*Q1).y += e.x * sc_ek;
        sc_ek = e.z * sc;
        (*Q2).x += sc_ek * e.z;
        (*Q1).z += e.x * sc_ek;
        (*Q2).z += e.y * sc_ek;
    }
}

template<int ewald_max>
__device__ void tensor_calc_pos_ewald( float4 a, float4 b, float* D, int LD )
{
    float4 v;
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    v.z = a.z - b.z;

    v.x *= dev_invL0_t;
    v.y *= dev_invL1_t;
    v.z *= dev_invL2_t;

    float3 O1, O2, Q1, Q2 = {0.0f, 0.0f, 0.0f};
    O1 = O2 = Q1 = Q2;
    int f, g, h;
    comp_OQ( ewald_max, 0, 0, v, &O1, &O2, &Q1, &Q2 );
    if ( ewald_max )comp_OQ( -ewald_max, 0, 0, v, &O1, &O2, &Q1, &Q2 );
    for ( h = -ewald_max + 1; h < ewald_max; ++h )
    {
        int next_f = sqrtf( dev_max_sq - h * h );
        for ( f = -next_f; f <= next_f; ++f )
        {
            int next_g = sqrtf( dev_max_sq - h * h - f * f );
            for ( g = -next_g; g <= next_g; ++g )
            {
                comp_OQ( h, f, g, v, &O1, &O2, &Q1, &Q2 );
            }
        }
    }
    const float preO = 0.75f * dev_invL0_t;
    const float sigma_sq = 0.5f * (a.w * a.w + b.w * b.w);
    const float preQ = 0.5f * sigma_sq * dev_invL0_t * dev_invL0_t * dev_invL0_t;
    float3 d1, d2, d3;
    d1.x = dev_preDii * (preO * O1.x + preQ * Q1.x);
    d1.y = dev_preDii * (preO * O1.y + preQ * Q1.y);
    d1.z = dev_preDii * (preO * O1.z + preQ * Q1.z);

    d2.x = d1.y;
    d2.y = dev_preDii * (preO * O2.y + preQ * Q2.y);
    d2.z = dev_preDii * (preO * O2.z + preQ * Q2.z);

    d3.x = d1.z;
    d3.y = d2.z;
    d3.z = dev_preDii * (preO * O2.x + preQ * Q2.x);

    *((float3*) D) = d1;
    *((float3*) (D + LD)) = d2;
    *((float3*) (D + 2 * LD)) = d3;
}

__device__ void tensor_calc_pos_correction( float4 a, float4 b, float* D, int LD )
{
    float4 r_ij;
    r_ij.x = a.x - b.x;
    r_ij.y = a.y - b.y;
    r_ij.z = a.z - b.z;

    float r_sq = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    if( r_sq < (a.w + b.w)*(a.w + b.w) )
    {
        float3 bd1,bd2,bd3;
        {
            float r_norm = sqrtf( r_sq );

            float sigma = (a.w + b.w) / 2;
            float pre = dev_preDii / sigma;

            float divr = sigma * r_norm * 32;

            bd3.z = bd2.y = bd1.x = 1 - (9 * r_norm / (32 * sigma));

            bd1.x += 3 * r_ij.x * r_ij.x / divr;
            bd1.x *= pre;

            bd2.y += 3 * r_ij.y * r_ij.y / divr;
            bd2.y *= pre;

            bd3.z += 3 * r_ij.z * r_ij.z / divr;
            bd3.z *= pre;

            bd2.x = bd1.y = pre * 3 * r_ij.x * r_ij.y / divr;
            bd3.x = bd1.z = pre * 3 * r_ij.x * r_ij.z / divr;
            bd2.z = bd3.y = pre * 3 * r_ij.z * r_ij.y / divr;
        }
        {
            float r_norm = sqrtf( r_sq );
            float3 d1,d2,d3;

            float pre = dev_preDij / r_norm;
            float sigm_r = (a.w * a.w + b.w * b.w) / r_sq;

            float divr = r_sq;
            float mult = 1 - sigm_r;
            float part1 = 1 + (sigm_r / 3);

            d3.z = d2.y = d1.x = part1;

            d1.x += mult * r_ij.x * r_ij.x / divr;
            d1.x *= pre;

            d2.y += mult * r_ij.y * r_ij.y / divr;
            d2.y *= pre;

            d3.z += mult * r_ij.z * r_ij.z / divr;
            d3.z *= pre;

            d2.x = d1.y = pre * mult * r_ij.x * r_ij.y / divr;
            d3.x = d1.z = pre * mult * r_ij.x * r_ij.z / divr;
            d2.z = d3.y = pre * mult * r_ij.z * r_ij.y / divr;

            bd1 = sub3( bd1, d1 );
            bd2 = sub3( bd2, d2 );
            bd3 = sub3( bd3, d3 );
        }

        *((float3*) D) = add3( bd1, *((float3*) D) );
        *((float3*) (D + LD)) = add3( bd2, *((float3*) (D + LD)) );
        *((float3*) (D + 2 * LD)) = add3( bd3, *((float3*) (D + 2 * LD)) ) ;
    }
}

template<int ewald_real>
__device__ void tensor_calc_pos_ewald_diag( float4 a, float* D, int LD )
{
    float3 tab;
    tensor_calc_pos_ewald<ewald_real > (a, a, (float*) &tab, 0);
    float sc3 = dev_ewald_alpha * a.w * dev_invL0_t;
    sc3 = sc3 * sc3 * sc3 * M1_SQRTPIF / 3;
    sc3 = 0.75f * dev_invL0_t * 1.5f * dev_ewald_alpha * a.w * M1_SQRTPIF * dev_invL0_t +
        0.5f * a.w * a.w * dev_invL0_t * dev_invL0_t * dev_invL0_t * sc3;
    sc3 *= dev_preDii;
    float3 tmp = {0.0f, 0.0f, 0.0f};
    tmp.x = dev_preDii / a.w - sc3 + tab.z; // z is only correct due to 0 in LD
    *((float3*) D) = tmp;
    tmp.y = tmp.x;
    tmp.x = 0;
    *((float3*) (D + LD)) = tmp;
    tmp.z = tmp.y;
    tmp.y = 0;
    *((float3*) (D + 2 * LD)) = tmp;
}

__device__ void comp_Dt( void* ptr_coord, void* D )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    float* Dtab = (float*) D;
    float4* ctab = (float4*) ptr_coord;
    const int LD = gridDim.x * blockDim.x * 4;

    if ( idx <= idy )
    {
        float4 curr_x = ctab[ idx ];
        float4 curr_y = ctab[ idy ];
        //__syncthreads( );
        if ( idx != idy )
            tensor_calc_pos( curr_y, curr_x, Dtab + 3 * idx * LD + 3 * idy, LD );
        else
            tensor_calc_pos_diag( curr_x, Dtab + idx * LD * 3 + idx * 3, LD );
    }
}

template<int ewald_real>
__device__ void comp_Dt_ewald_sq( void* ptr_coord, void* D )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    float* Dtab = (float*) D;
    float4* ctab = (float4*) ptr_coord;
    const int LD = gridDim.x * blockDim.x * 4;

    if ( idx <= idy )
    {
        float4 curr_x = ctab[ idx ];
        float4 curr_y = ctab[ idy ];
        //__syncthreads( );
        if ( idx != idy )
            tensor_calc_pos_ewald<ewald_real > (curr_y, curr_x, Dtab + 3 * idx * LD + 3 * idy, LD);
        else
            tensor_calc_pos_ewald_diag<ewald_real > (curr_x, Dtab + idx * LD * 3 + idx * 3, LD);
    }
}

__device__ void comp_Dt_ewald_correction( void* ptr_coord, void* D )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    float* Dtab = (float*) D;
    float4* ctab = (float4*) ptr_coord;
    const int LD = gridDim.x * blockDim.x * 4;

    if ( idx < idy )
    {
        float4 curr_x = ctab[ idx ];
        float4 curr_y = ctab[ idy ];
        tensor_calc_pos_correction(curr_y, curr_x, Dtab + 3 * idx * LD + 3 * idy, LD);
    }
}

template<bool ewald, int ewald_real >
__global__ void choose_comp_D( void* coord_ptr, void* _D )
{
    if ( !ewald )
    {
        comp_Dt( coord_ptr, _D );
    }
    else
    {
        comp_Dt_ewald_sq< ewald_real >( coord_ptr, _D );
    }
}

__global__ void choose_comp_D_correction( void* coord_ptr, void* _D )
{
    comp_Dt_ewald_correction( coord_ptr, _D );
}

extern "C"
void compute_D( DOUBLE* _D )
{
    if ( !hydro ) return; /* End if no hydro */

    /*double mintime = 111111;
    int mx,my;
    for ( int xtdim = 1; xtdim <= 256; xtdim*=2)
    for ( int ytdim = 1; ytdim <= 256; ytdim*=2)
    if(xtdim*ytdim<=128)
    {
    printf("%d %d\n",xtdim,ytdim);
    mstart();*/
    if ( !ewald )
    {
        const int xtdim = 1;
        const int ytdim = (128 < size ) ? 128 : size;
        dim3 dim_block( xtdim, ytdim );
        dim3 dim_grid( size / xtdim, size / ytdim );
        if ( bc == BC_BOX ) //4 float4, 4 sizeof(float)
        {
            choose_comp_D< false, 0 > <<< dim_grid, dim_block >>>(d_coord, _D);
        }
        else
        {
            choose_comp_D< false, 0 > <<< dim_grid, dim_block >>>(d_coord, _D);
        }
    }
    else
    {
        int xtdim = 2;
        int ytdim = 64;
        if ( xtdim*ytdim > size )
        {
            xtdim = 1;
            ytdim = size;
        }
        dim3 dim_block( xtdim, ytdim );
        dim3 dim_grid( size / xtdim, size / ytdim );
        if( ewald_method == EWALD_METHOD_SMITH )
        {
            if ( ewald_real == 0 && ewald_recip == 0 )
                choose_comp_D<true, 0 > <<< dim_grid, dim_block >>>(d_coord, _D);
            else if ( ewald_real == 1 && ewald_recip == 1 )
                choose_comp_D<true, 1 > <<< dim_grid, dim_block >>>(d_coord, _D);
            else if ( ewald_real == 2 && ewald_recip == 2 )
                choose_comp_D<true, 2 > <<< dim_grid, dim_block >>>(d_coord, _D);
            else if ( ewald_real == 3 && ewald_recip == 3 )
                choose_comp_D<true, 3 > <<< dim_grid, dim_block >>>(d_coord, _D);
            else if ( ewald_real == 4 && ewald_recip == 4 )
                choose_comp_D<true, 4 > <<< dim_grid, dim_block >>>(d_coord, _D);
            else if ( ewald_real == 5 && ewald_recip == 5 )
                choose_comp_D<true, 5 > <<< dim_grid, dim_block >>>(d_coord, _D);
            else
                UNERR( "Wrong ewald values")
        }
        else
        {
            UNERR( "Unimplemented ewald method beenakker")
        }
        choose_comp_D_correction <<< dim_grid, dim_block >>>( d_coord, _D );//
    }
    cudaThreadSynchronize( ); CCERR
    /*mstop();
    extern double g_start;
    if( g_start < mintime )
    {
        mintime = g_start;
        my = ytdim;
        mx = xtdim;
    }
    }
    printf("Winner %g, %d %d\n",mintime, mx,my);*/
}

extern "C"
void init_tensor()
{
    float fDii = preDii, fDij = preDij, ewa = ewald_alpha;
    cudaMemcpyToSymbol( "dev_preDii", &fDii, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_preDij", &fDij, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_ewald_alpha", &ewa, sizeof ( float), 0, cudaMemcpyHostToDevice ); CCERR

    float ewa2 = ewa * ewa;
    cudaMemcpyToSymbol( "dev_ewald_alpha2", &ewa2, sizeof ( float), 0, cudaMemcpyHostToDevice ); CCERR
    float ewa3 = ewa * ewa2;
    cudaMemcpyToSymbol( "dev_ewald_alpha3", &ewa3, sizeof ( float), 0, cudaMemcpyHostToDevice ); CCERR
    float M_PI_ALPHA_SQ = -M_PIF * M_PIF / ewa2;
    cudaMemcpyToSymbol( "dev_M_PI_ALPHA_SQ", &M_PI_ALPHA_SQ, sizeof ( float), 0, cudaMemcpyHostToDevice ); CCERR
    float max_sq = (float)ewald_real*ewald_real;
    cudaMemcpyToSymbol( "dev_max_sq", &max_sq, sizeof ( float), 0, cudaMemcpyHostToDevice ); CCERR

    float bx[3];
    float iL[3];
    bx[0] = box[0], bx[1] = box[1], bx[2] = box[2];
    iL[0] = inv_L[0], iL[1] = inv_L[1], iL[2] = inv_L[2];
    cudaMemcpyToSymbol( "dev_box0_t", bx, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_box1_t", bx + 1, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_box2_t", bx + 2, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_invL0_t", iL, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_invL1_t", iL + 1, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    cudaMemcpyToSymbol( "dev_invL2_t", iL + 2, sizeof (float), 0, cudaMemcpyHostToDevice ); CCERR
    if( hydro && !( hydro == DIFF_GEYER && algorithm < DIFF_ALG_IGT_CONST ) )
    {
        cudaMemset( d_D, 0, size * DIMS1 * size * DIMS0 * sizeof (float) ); CCERR
    }
}

void cm_dh_tensor( DOUBLE* _D )
{
    INT i, j;
    DOUBLE* h_D = (DOUBLE*) malloc( size_D ); CHMEM(h_D);
    cudaMemcpy( h_D, d_D, size_D, cudaMemcpyDeviceToHost ); CCERR
    memset( _D, 0, DIMS0 * size * DIMS0 * size * sizeof (_D[0]) );
    for ( i = 0; i < DIMS0 * size; ++i )
    {
        for ( j = 0; j < DIMS0 * size; ++j )
        {
            if ( i <= j )
            {
                _D[ j * DIMS0 * size + i ] = h_D[ i * DIMS1 * size + j ];
            }
        }
    }
    free( h_D );
}

void compute_single_D( DOUBLE* _D )
{
    compute_D( NULL );
    cm_dh_tensor( _D );
}

template<bool ewald, int ewald_real>
__global__ void rinfly( float4* coords, float3* forces, float3* r )
{
    int tid = threadIdx.x;
    int id = blockIdx.x*blockDim.x + tid;
    extern __shared__ float sh[];
    float4* sh_coords = (float4*) sh;
    float3* sh_forces = (float3*)( sh_coords + blockDim.x );

    float4 curr_coords = coords[id];
    float3 out_r = {0.0f, 0.0f, 0.0f};

    for ( int i = 0; i < gridDim.x; ++i )
    {
        sh_coords[tid] = coords[ i*blockDim.x + tid ];
        sh_forces[tid] = forces[ i*blockDim.x + tid ];
        __syncthreads();
        for ( int k = 0; k < blockDim.x; ++k )
        {
            int jd = i * blockDim.x + k;
            if ( jd == id )//branch divergent on diagonal but only n-times
            {
                float3 tmp;
                if ( ewald ) tensor_calc_pos_ewald_diag<ewald_real>( curr_coords, (float*) &tmp, 0 );
                else tensor_calc_pos_diag( curr_coords, (float*) &tmp, 0 );
                out_r.x += tmp.z * sh_forces[k].x;
                out_r.y += tmp.z * sh_forces[k].y;
                out_r.z += tmp.z * sh_forces[k].z;
            }
            else
            {
                float t[9];
                if ( ewald ) tensor_calc_pos_ewald<ewald_real>( curr_coords, sh_coords[k], t, 3 );
                else tensor_calc_pos( curr_coords, sh_coords[k], t, 3 );
                float3 f = sh_forces[k];
                out_r.x += f.x * t[0] + f.y * t[1] + f.z * t[2];
                out_r.y += f.x * t[3] + f.y * t[4] + f.z * t[5];
                out_r.z += f.x * t[6] + f.y * t[7] + f.z * t[8];
            }
        }
        __syncthreads();
    }
    r[id] = out_r;
}

void comp_rinfly( float* coord, float* forces, float* r )
{
    if( !ewald )
    {
        rinfly<false,0><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
    }
    else
    {
        if ( ewald_real == 0 && ewald_recip == 0 )
            rinfly<true,0><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
        else if ( ewald_real == 1 && ewald_recip == 1 )
            rinfly<true,1><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
        else if ( ewald_real == 2 && ewald_recip == 2 )
            rinfly<true,2><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
        else if ( ewald_real == 3 && ewald_recip == 3 )
            rinfly<true,3><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
        else if ( ewald_real == 4 && ewald_recip == 4 )
            rinfly<true,4><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
        else if ( ewald_real == 5 && ewald_recip == 5 )
            rinfly<true,5><<<blocks, cuda_block, cuda_block*4*7>>>
                    ( (float4*)coord, (float3*)forces, (float3*)r );
        else
            UNERR( "Wrong ewald values")
    }
}

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

__device__ float3 add( float3 a, float3 b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }

template<int cuda_block, bool ewald, int ewald_real>
__global__ void comp_geyer_beg_in_fly_kernel( float4* coords, float3* ptr_x,
                                              float3* ptr_eff, float3* ptr_c,
                                              float* ptr_eps )
{
    float eps_curr = 0;
    float3 Feff = {0.0f,0.0f,0.0f};
    float3 C = {0.0f,0.0f,0.0f};
    int tid = threadIdx.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float sh[];
    float4* sh_coords = (float4*)sh;
    float3* sh_X = (float3*)( sh + 4*blockDim.x );
    float4 curr_coords = coords[ id ];
    for ( int k = 0; k < gridDim.x; ++k )
    {
        sh_coords[tid] = coords[k*blockDim.x+tid];
        sh_X[tid] = ptr_x[k*blockDim.x+tid];
        __syncthreads();
        for ( int j = 0; j < blockDim.x; ++j )
        {
            int jd = k*blockDim.x + j;
            if( jd != id )
            {
                float4 k_coord = sh_coords[j];
                float3 tmp;
                if ( ewald ) tensor_calc_pos_ewald_diag<ewald_real>( k_coord, (float*) &tmp, 0 );
                else tensor_calc_pos_diag( k_coord, (float*) &tmp, 0 );
                float dk = tmp.z;
                float t[9]; 
                if ( ewald ) tensor_calc_pos_ewald<ewald_real>( curr_coords, k_coord, t, 3 );
                else tensor_calc_pos( curr_coords, k_coord, t, 3 );
                float3 X = sh_X[j];
                Feff.x += X.x*t[0] + X.y*t[3] + X.z*t[6];
                Feff.y += X.x*t[1] + X.y*t[4] + X.z*t[7];
                Feff.z += X.x*t[2] + X.y*t[5] + X.z*t[8];
                //Feff += ptr_x[k] * dik;
                C.x += (t[0]*t[0] + t[3]*t[3] + t[6]*t[6])/dk;
                C.y += (t[1]*t[1] + t[4]*t[4] + t[7]*t[7])/dk;
                C.z += (t[2]*t[2] + t[5]*t[5] + t[8]*t[8])/dk;
                //C += dik*dik/dk;
                eps_curr += t[0]+t[1]+t[2]+t[3]+t[4]+t[5]+t[6]+t[7]+t[8];
            }
        }
        __syncthreads();
    }
    ptr_eff[id] = Feff;
    ptr_c[id] = C;
    eps_reduction<cuda_block>( eps_curr/2, ptr_eps );
}

template<int cuda_block, bool ewald, int ewald_real>
__global__ void comp_geyer_end_in_fly( float4* coords, float* ptr_x, float* ptr_eff, float* ptr_c, float beta, float* out)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    float Feff = ptr_eff[id];
    float C = ptr_c[id];

    float dii;
    float3 tmp;
    if ( ewald ) tensor_calc_pos_ewald_diag<ewald_real>( coords[id/3], (float*) &tmp, 0 );
    else tensor_calc_pos_diag( coords[id/3], (float*) &tmp, 0 );
    dii = tmp.z;
    C *= beta*beta/dii;
    C = rsqrtf( 1 + C );
    out[id] += C*(beta*Feff + dii*ptr_x[id]);
}

template<int cuda_block>
void comp_geyer_beg_in_fly_ewald( float* coord, float* x, float* eff, float* c, float* eps )
{
    if ( ewald_real == ewald_recip )
    {
        switch( ewald_real )
        {
            case 0: comp_geyer_beg_in_fly_kernel<cuda_block,true,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );return;
            case 1: comp_geyer_beg_in_fly_kernel<cuda_block,true,1><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );return;
            case 2: comp_geyer_beg_in_fly_kernel<cuda_block,true,2><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );return;
            case 3: comp_geyer_beg_in_fly_kernel<cuda_block,true,3><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );return;
            case 4: comp_geyer_beg_in_fly_kernel<cuda_block,true,4><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );return;
            case 5: comp_geyer_beg_in_fly_kernel<cuda_block,true,5><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );return;
        }
    }
    UNERR("Unimplemented")
}

void comp_geyer_beg_in_fly( float* coord, float* x, float* eff, float* c, float* eps )
{
    if( !ewald )
    {
        switch(cuda_block)
        {
            case(1024):comp_geyer_beg_in_fly_kernel<1024,false,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );break;
            case( 512):comp_geyer_beg_in_fly_kernel< 512,false,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );break;
            case( 256):comp_geyer_beg_in_fly_kernel< 256,false,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );break;
            case( 128):comp_geyer_beg_in_fly_kernel< 128,false,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );break;
            case(  64):comp_geyer_beg_in_fly_kernel<  64,false,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );break;
            case(  32):comp_geyer_beg_in_fly_kernel<  32,false,0><<< blocks, cuda_block, 7*4*cuda_block >>>( (float4*)coord, (float3*)x, (float3*)eff, (float3*)c, eps );break;
            default: UNERR("Unimplemented")
        }
    }
    else
    {
        switch(cuda_block)
        {
            case(1024): comp_geyer_beg_in_fly_ewald<1024>( coord, x, eff, c, eps );break;
            case( 512): comp_geyer_beg_in_fly_ewald< 512>( coord, x, eff, c, eps );break;
            case( 256): comp_geyer_beg_in_fly_ewald< 256>( coord, x, eff, c, eps );break;
            case( 128): comp_geyer_beg_in_fly_ewald< 128>( coord, x, eff, c, eps );break;
            case(  64): comp_geyer_beg_in_fly_ewald<  64>( coord, x, eff, c, eps );break;
            case(  32): comp_geyer_beg_in_fly_ewald<  32>( coord, x, eff, c, eps );break;
            default: UNERR("Unimplemented")
        }
    }
}

template<int cuda_block>
void comp_geyer_end_in_fly_ewald( float* coord, float* ptr_x, float* ptr_eff,
                                       float* ptr_c, float beta, float* out)
{
    if ( ewald_real == ewald_recip )
    {
        switch( ewald_real )
        {
            case 0: comp_geyer_end_in_fly<cuda_block,true,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case 1: comp_geyer_end_in_fly<cuda_block,true,1><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case 2: comp_geyer_end_in_fly<cuda_block,true,2><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case 3: comp_geyer_end_in_fly<cuda_block,true,3><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case 4: comp_geyer_end_in_fly<cuda_block,true,4><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case 5: comp_geyer_end_in_fly<cuda_block,true,5><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
        }
    }
    UNERR("Unimplemented ewald")
}

void comp_geyer_end_in_fly( float* coord, float* ptr_x, float* ptr_eff,
                                       float* ptr_c, float beta, float* out)
{
    if( !ewald )
    {
        switch(cuda_block)
        {
            case(1024):comp_geyer_end_in_fly<1024,false,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case( 512):comp_geyer_end_in_fly< 512,false,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case( 256):comp_geyer_end_in_fly< 256,false,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case( 128):comp_geyer_end_in_fly< 128,false,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case(  64):comp_geyer_end_in_fly<  64,false,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case(  32):comp_geyer_end_in_fly<  32,false,0><<< DIMS0*blocks, cuda_block >>>( (float4*)coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
        }
    }
    else
    {
        switch(cuda_block)
        {
            case(1024):comp_geyer_end_in_fly_ewald<1024>( coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case( 512):comp_geyer_end_in_fly_ewald< 512>( coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case( 256):comp_geyer_end_in_fly_ewald< 256>( coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case( 128):comp_geyer_end_in_fly_ewald< 128>( coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case(  64):comp_geyer_end_in_fly_ewald<  64>( coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
            case(  32):comp_geyer_end_in_fly_ewald<  32>( coord, ptr_x, ptr_eff, ptr_c, beta, out );return;
        }
    }
    UNERR("Unimplemented in fly")
}
