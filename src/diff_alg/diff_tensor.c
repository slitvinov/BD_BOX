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

#include "diff_tensor.h"

#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../err.h"
#include "../input.h"
#include "../main_loop.h"
#include "../math_help.h"
#include "../myblas.h"
#include "../rand_move.h"
#include "../trans.h"

void diff_tensor_calc_pos( INT a, INT b, DOUBLE* tD )
{
    DOUBLE sigma_i = coord[ DIMS1 * b + 3 ];
    DOUBLE sigma_j = coord[ DIMS1 * a + 3 ];
    DOUBLE r_norm;
    DOUBLE r_sq;
    DOUBLE e[3];
    int i, j;

    get_v2_norms_no( b, a, e, &r_norm, &r_sq );
    if ( r_norm < sigma_i + sigma_j )
    {
        DOUBLE sigma = (sigma_i + sigma_j) / 2;
        DOUBLE pre = preDii / sigma;
        DOUBLE part1 = (1 - (9 * r_norm / (32 * sigma)));
        tD[0] = part1;
        tD[4] = part1;
        tD[8] = part1;
        for ( i = 0; i < 3; ++i )
        {
            tD[ i + 3 * i ] += 3 * e[i] * e[i] * r_norm / (32 * sigma);
            tD[ i + 3 * i ] *= pre;
            for ( j = 0; j < i; ++j )
            {
                tD[ j + 3 * i ] = pre * 3 * e[i] * e[j] * r_norm / (32 * sigma);
                tD[ i + 3 * j ] = tD[ j + 3 * i ];
            }
        }
    }
    else
    {
        DOUBLE pre = preDij / r_norm;
        DOUBLE sigm_r = (sigma_i * sigma_i + sigma_j * sigma_j) / r_sq;
        tD[0] = 1 + sigm_r / 3;
        tD[4] = 1 + sigm_r / 3;
        tD[8] = 1 + sigm_r / 3;
        for ( i = 0; i < 3; ++i )
        {
            tD[ i + 3 * i ] += e[i] * e[i]*(1 - sigm_r);
            tD[ i + 3 * i ] *= pre;
            for ( j = 0; j < i; ++j )
            {
                tD[ j + 3 * i ] = pre * e[i] * e[j]*(1 - sigm_r);
                tD[ i + 3 * j ] = tD[ j + 3 * i ];
            }
        }
    }
}

DOUBLE preO = 0;
INT ewald_rr = 0;
INT rr_sq = 0;
DOUBLE alpha2 = 0;
DOUBLE alpha3 = 0;
DOUBLE inv_alpha2 = 0;
DOUBLE M_PI_ALPHA_SQ = 0;

void diff_tensor_calc_ewald_smith( INT i, INT j, DOUBLE* _tD )
{
    DOUBLE sigma_i = coord[ DIMS1 * j + 3 ];
    DOUBLE sigma_j = coord[ DIMS1 * i + 3 ];
    DOUBLE v[3];
    const DOUBLE sigma_sq = 0.5 * (sigma_i * sigma_i + sigma_j * sigma_j);
    INT k, l;
    DOUBLE Op[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    DOUBLE Qp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    INT f, g, h;
    const DOUBLE preQ = 0.5 * sigma_sq * inv_L[0] * inv_L[0] * inv_L[0];
    dist_vec( coord + j*DIMS1, coord + i*DIMS1, v );

    if ( (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) < (sigma_i + sigma_j)*(sigma_i + sigma_j) && i != j )
    {
        DOUBLE r_norm, r_sq, e[3];
        get_ev2_norms( v[0], v[1], v[2], e, &r_norm, &r_sq );
        {/*add overlapping Dij*/
            DOUBLE sigma = (sigma_i + sigma_j) / 2;
            DOUBLE pre = preDii / sigma;
            DOUBLE part1 = (1 - (9 * r_norm / (32 * sigma)));
            for ( i = 0; i < 3; ++i )
            {
                _tD[ i + 3 * i ] += pre * (part1 + 3 * e[i] * e[i] * r_norm / (32 * sigma));
                for ( j = 0; j < i; ++j )
                {
                    DOUBLE vval = pre * 3 * e[i] * e[j] * r_norm / (32 * sigma);
                    _tD[ j + 3 * i ] += vval;
                    _tD[ i + 3 * j ] += vval;
                }
            }
        }
        {/*remove none overlapping Dij*/
            DOUBLE pre = preDij / r_norm;
            DOUBLE sigm_r = (sigma_i * sigma_i + sigma_j * sigma_j) / r_sq;
            DOUBLE part1 = 1 + sigm_r / 3;
            for ( i = 0; i < 3; ++i )
            {
                _tD[ i + 3 * i ] -= pre*(part1 + e[i] * e[i]*(1 - sigm_r));
                for ( j = 0; j < i; ++j )
                {
                    DOUBLE vval = pre * e[i] * e[j] * (1 - sigm_r);
                    _tD[ j + 3 * i ] -= vval;
                    _tD[ i + 3 * j ] -= vval;
                }
            }
        }
    }

    v[0] *= inv_L[0];
    v[1] *= inv_L[1];
    v[2] *= inv_L[2];

    for ( h = -ewald_rr; h <= ewald_rr; ++h )
    {
        int next_f = (int)SQRTD( rr_sq - h*h );
        for ( f = -next_f; f <= next_f; ++f )
        {
            int next_g = (int)SQRTD( rr_sq - h*h - f*f );
            for ( g = -next_g; g <= next_g; ++g )
            {
                
                    DOUBLE r_norm, r_sq, e[3];
                    DOUBLE inv_r;
                    DOUBLE inv_r3;
                    DOUBLE sc;
                    DOUBLE sc2;

                    get_ev2_norms( v[0] + h, v[1] + f, v[2] + g, e, &r_norm, &r_sq );
                    if ( r_norm != 0.0f )
                    {
                        inv_r = 1 / r_norm;
                        inv_r3 = 1 / (r_norm * r_sq);
                        sc = ERFCD( r_norm * ewald_alpha ) * inv_r;
                        sc2 = 2 * ewald_alpha * M1_SQRTPI * EXPD( -alpha2 * r_sq );
                        for ( k = 0; k < 3; ++k )
                        {
                            DOUBLE sc_ek = e[k] * (sc + sc2);
                            Op[ k * 4 ] += sc + sc_ek * e[k];
                            for ( l = 0; l < k; ++l )
                            {
                                Op[ k + l * 3 ] += e[l] * sc_ek;
                            }
                        }
                        sc = inv_r3 * (ERFCD( r_norm * ewald_alpha ) + 2 * ewald_alpha * r_norm * M1_SQRTPI * EXPD( -alpha2 * r_sq ));
                        sc2 = -4 * alpha3 * M1_SQRTPI * EXPD( -alpha2 * r_sq );
                        for ( k = 0; k < 3; ++k )
                        {
                            DOUBLE sc_ek = (-3 * sc + sc2) * e[k];
                            Qp[ k * 4 ] += sc + sc_ek * e[k];
                            for ( l = 0; l < k; ++l )
                            {
                                Qp[ k + l * 3 ] += e[l] * sc_ek;
                            }
                        }
                    }
                /*recip*/
                    get_ev2_norms( h, f, g, e, &r_norm, &r_sq );
                    if ( r_norm != 0.0f )
                    {
                        sc = 2 * EXPD( M_PI_ALPHA_SQ * r_sq ) / (M_PIF * r_sq);
                        sc2 = -1 + M_PI_ALPHA_SQ * r_sq;
                        sc *= cos( 2 * M_PIF * (v[0] * h + v[1] * f + v[2] * g) );
                        sc2 *= sc;
                        for ( k = 0; k < 3; ++k )
                        {
                            DOUBLE sc_ek = sc2 * e[k];
                            Op[ k * 4 ] += sc + sc_ek * e[k];
                            for ( l = 0; l < k; ++l )
                            {
                                Op[ k + l * 3 ] += e[l] * sc_ek;
                            }
                        }
                        sc = 4 * M_PIF * EXPD( M_PI_ALPHA_SQ * r_sq );
                        sc *= cos( 2 * M_PIF * (v[0] * h + v[1] * f + v[2] * g) );
                        for ( k = 0; k < 3; ++k )
                        {
                            DOUBLE sc_ek = sc * e[k];
                            Qp[ k * 4 ] += sc_ek * e[k];
                            for ( l = 0; l < k; ++l )
                            {
                                Qp[ k + l * 3 ] += e[l] * sc_ek;
                            }
                        }
                    }
            }
        }
    }
    for ( k = 0; k < 3; ++k )
    {
        _tD[ k * 4 ] += preDii * (preO * Op[k * 4] + preQ * Qp[k * 4]);
        for ( l = 0; l < k; ++l )
        {
            _tD[ k + l * 3 ] += preDii * (preO * Op[ k + l * 3 ] + preQ * Qp[ k + l * 3 ]);
            _tD[ l + k * 3 ] += preDii * (preO * Op[ k + l * 3 ] + preQ * Qp[ k + l * 3 ]);
        }
    }
}

void diff_tensor_calc_ewald_smith_diag( INT i, DOUBLE* tD )
{
    INT k, l;
    DOUBLE aw = coord[i * DIMS1 + 3];
    DOUBLE sc3;

    diff_tensor_calc_ewald_smith( i, i, tD );

    sc3 = ewald_alpha * aw * inv_L[0];
    sc3 = sc3 * sc3 * sc3 * M1_SQRTPI / 3;
    sc3 = 0.75f * inv_L[0] * 1.5f * ewald_alpha * aw * M1_SQRTPI * inv_L[0] +
        0.5f * aw * aw * inv_L[0] * inv_L[0] * inv_L[0] * sc3;
    sc3 *= preDii;
    for ( k = 0; k < 3; ++k )
    {
        tD[ k * 4 ] = preDii / aw - sc3 + tD[ k * 4 ];
        for ( l = 0; l < k; ++l )
        {
            tD[ l + 3 * k ] = 0;
            tD[ k + 3 * l ] = 0;
        }
    }
}

void diff_tensor_calc_ewald_beenakker( INT i, INT j, DOUBLE* _tD )
{
    DOUBLE sigma_i = coord[ DIMS1 * j + 3 ];
    DOUBLE sigma_j = coord[ DIMS1 * i + 3 ];
    DOUBLE v[3];
    DOUBLE r[3];
    const DOUBLE sigma_sq = 0.5 * (sigma_i * sigma_i + sigma_j * sigma_j);
    INT k, l;
    DOUBLE Op[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    DOUBLE Qp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    INT f, g, h;
    const DOUBLE preQ = inv_L[0] * inv_L[0] * inv_L[0];
    dist_vec( coord + i*DIMS1, coord + j*DIMS1, r );

    if ( (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) < (sigma_i + sigma_j)*(sigma_i + sigma_j) && i != j )
    {
        DOUBLE r_norm, r_sq, e[3];
        get_ev2_norms( r[0], r[1], r[2], e, &r_norm, &r_sq );
        {/*add overlapping Dij*/
            DOUBLE sigma = (sigma_i + sigma_j) / 2;
            DOUBLE pre = preDii / sigma;
            DOUBLE part1 = (1 - (9 * r_norm / (32 * sigma)));
            for ( i = 0; i < 3; ++i )
            {
                _tD[ i + 3 * i ] += pre * (part1 + 3 * e[i] * e[i] * r_norm / (32 * sigma));
                for ( j = 0; j < i; ++j )
                {
                    DOUBLE vval = pre * 3 * e[i] * e[j] * r_norm / (32 * sigma);
                    _tD[ j + 3 * i ] += vval;
                    _tD[ i + 3 * j ] += vval;
                }
            }
        }
        {/*remove none overlapping Dij*/
            DOUBLE pre = preDij / r_norm;
            DOUBLE sigm_r = (sigma_i * sigma_i + sigma_j * sigma_j) / r_sq;
            DOUBLE part1 = 1 + sigm_r / 3;
            for ( i = 0; i < 3; ++i )
            {
                _tD[ i + 3 * i ] -= pre*(part1 + e[i] * e[i]*(1 - sigm_r));
                for ( j = 0; j < i; ++j )
                {
                    DOUBLE vval = pre * e[i] * e[j] * (1 - sigm_r);
                    _tD[ j + 3 * i ] -= vval;
                    _tD[ i + 3 * j ] -= vval;
                }
            }
        }
    }

    v[0] = r[0] * inv_L[0];
    v[1] = r[1] * inv_L[1];
    v[2] = r[2] * inv_L[2];

    for ( h = -ewald_rr; h <= ewald_rr; ++h )
    {
        int next_f = (int)SQRTD( rr_sq - h*h );
        for ( f = -next_f; f <= next_f; ++f )
        {
            int next_g = (int)SQRTD( rr_sq - h*h - f*f );
            for ( g = -next_g; g <= next_g; ++g )
            {

                    DOUBLE r_norm, r_sq, e[3];
                    DOUBLE inv_r, inv_r2, inv_r3;
                    DOUBLE sc, sc2;

                    get_ev2_norms( r[0] + box[0]*h, r[1] + box[0]*f, r[2] + box[0]*g, e, &r_norm, &r_sq );
                    if ( r_norm != 0.0f )
                    {
                        inv_r = 1 / r_norm;
                        inv_r2 = 1 / r_sq;
                        inv_r3 = inv_r*inv_r2;
                        sc = ewald_alpha*sigma_sq*inv_r2;
                        sc += 14*alpha3*sigma_sq;
                        sc += -4.5f*ewald_alpha;
                        sc += -20*alpha2*alpha3*sigma_sq*r_sq;
                        sc += 3*alpha3*r_sq;
                        sc += 4*alpha3*alpha2*alpha2*sigma_sq*r_sq*r_sq;
                        sc *= EXPD(-alpha2*r_sq);
                        sc /= SQRTPI;
                        sc += ERFCD( r_norm * ewald_alpha ) * ( 0.75f*inv_r + 0.5f*sigma_sq*inv_r3);

                        sc2 = -3*ewald_alpha*sigma_sq*inv_r2;
                        sc2 += -2*alpha3*sigma_sq;
                        sc2 += 1.5f*ewald_alpha;
                        sc2 += 16*alpha2*alpha3*sigma_sq*r_sq;
                        sc2 += -3*alpha3*r_sq;
                        sc2 += -4*alpha3*alpha2*alpha2*sigma_sq*r_sq*r_sq;
                        sc2 *= EXPD(-alpha2*r_sq);
                        sc2 /= SQRTPI;
                        sc2 += ERFCD( r_norm * ewald_alpha ) * ( 0.75f*inv_r - 1.5f*sigma_sq*inv_r3);

                        for ( k = 0; k < 3; ++k )
                        {
                            DOUBLE sc_ek = e[k] * sc2;
                            Op[ k * 4 ] += sc + sc_ek * e[k];
                            for ( l = 0; l < k; ++l )
                            {
                                Op[ k + l * 3 ] += e[l] * sc_ek;
                            }
                        }
                    }
                /*recip*/
                    get_ev2_norms( h, f, g, e, &r_norm, &r_sq );
                    if ( r_norm != 0.0f )
                    {
                        r_norm *= 2*M_PIF;
                        r_sq *= 4*M_PIF*M_PIF;
                        sc = 1 - sigma_sq*r_sq/3;
                        sc *= 1 + (inv_alpha2*r_sq/4)*(1+inv_alpha2*r_sq/2);
                        sc *= 6* M_PIF *EXPD(-inv_alpha2*r_sq/4)/r_sq;
                        sc *= cos( 2 * M_PIF * (v[0] * h + v[1] * f + v[2] * g) );
                        for ( k = 0; k < 3; ++k )
                        {
                            DOUBLE sc_ek = sc * e[k];
                            Qp[ k * 4 ] += sc - sc_ek * e[k];
                            for ( l = 0; l < k; ++l )
                            {
                                Qp[ k + l * 3 ] -= e[l] * sc_ek;
                            }
                        }
                    }
            }
        }
    }
    for ( k = 0; k < 3; ++k )
    {
        _tD[ k * 4 ] += preDii * ( Op[k * 4] + preQ * Qp[k * 4] );
        for ( l = 0; l < k; ++l )
        {
            _tD[ k + l * 3 ] += preDii * ( Op[ k + l * 3 ] + preQ * Qp[ k + l * 3 ]);
            _tD[ l + k * 3 ] += preDii * ( Op[ k + l * 3 ] + preQ * Qp[ k + l * 3 ]);
        }
    }
}


void diff_tensor_calc_ewald_beenakker_diag( INT i, DOUBLE* tD )
{
    INT k, l;
    DOUBLE aw = coord[i * DIMS1 + 3];
    DOUBLE sc3;

    diff_tensor_calc_ewald_beenakker( i, i, tD );

    sc3 = -6*ewald_alpha/SQRTPI + 40*alpha3*aw*aw/(3*SQRTPI);
    sc3 *= preDii;
    for ( k = 0; k < 3; ++k )
    {
        tD[ k * 4 ] = preDii / aw + sc3 + tD[ k * 4 ];
        for ( l = 0; l < k; ++l )
        {
            tD[ l + 3 * k ] = 0;
            tD[ k + 3 * l ] = 0;
        }
    }
}


void insert9( INT i, INT j, DOUBLE* _D, DOUBLE* tD, INT LD )
{
    INT k, l;
    for ( k = 0; k < 3; ++k )
    {
        for ( l = 0; l < 3; ++l )
        {
            _D[ k + l * LD + (j + i * LD) * DIMS0 ] = tD[ k + l * 3 ];
        }
    }
}

void compute_Dij( INT i, INT j, DOUBLE* tD )
{
    if( i == j )
    {
        if ( !ewald )
        {
            tD[1]=tD[2]=tD[3]=tD[5]=tD[6]=tD[7]=0.0;
            tD[0]=tD[4]=tD[8]=diag_D[i * DIMS0];
        }
        else if( ewald_method == EWALD_METHOD_SMITH )
        {
            diff_tensor_calc_ewald_smith_diag( i, tD );
        }
        else if( ewald_method == EWALD_METHOD_BEENAKKER )
        {
            diff_tensor_calc_ewald_beenakker_diag( i, tD );
        }
    }
    else
    {
        if( !ewald )
        {
            diff_tensor_calc_pos( i, j, tD );
        }
        else if( ewald_method == EWALD_METHOD_SMITH )
        {
            diff_tensor_calc_ewald_smith( i, j, tD );
        }
        else if( ewald_method == EWALD_METHOD_BEENAKKER )
        {
            diff_tensor_calc_ewald_beenakker( i, j, tD );
        }
    }
}

void compute_D( DOUBLE* _D )
{
    INT i;
    INT LD = size*DIMS0;
    if ( !hydro ) return; /* End if no hydro */

    if ( !ewald )
    {
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
        for ( i = 0; i < size; ++i )
        {
            INT j;
            for ( j = g_id; j < i; j += g_numprocs )
            {
                DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                diff_tensor_calc_pos( i, j, tD );
                insert9( i, j, _D, tD, LD );
            }
            _D[ i * DIMS0 * LD + i * DIMS0 ] = diag_D[i * DIMS0];
            _D[ (i * DIMS0 + 1) * LD + i * DIMS0 ] = 0.0;
            _D[ (i * DIMS0 + 1) * LD + i * DIMS0 + 1 ] = diag_D[i * DIMS0+1];
            _D[ (i * DIMS0 + 2) * LD + i * DIMS0 ] = 0.0;
            _D[ (i * DIMS0 + 2) * LD + i * DIMS0 + 1 ] = 0.0;
            _D[ (i * DIMS0 + 2) * LD + i * DIMS0 + 2 ] = diag_D[i * DIMS0+2];
        }
    }
    else if( ewald_method == EWALD_METHOD_SMITH )
    {
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
        for ( i = 0; i < size; ++i )
        {
            INT j;
            for ( j = g_id; j < i; j += g_numprocs )
            {
                DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                diff_tensor_calc_ewald_smith( i, j, tD );
                insert9( i, j, _D, tD, LD );
            }
            /* diag elements */
            {
                DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                diff_tensor_calc_ewald_smith_diag( i, tD );
                insert9( i, i, _D, tD, LD );
            }
        }
    }
    else if( ewald_method == EWALD_METHOD_BEENAKKER )
    {
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
        for ( i = 0; i < size; ++i )
        {
            INT j;
            for ( j = g_id; j < i; j += g_numprocs )
            {
                DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                diff_tensor_calc_ewald_beenakker( i, j, tD );
                insert9( i, j, _D, tD, LD );
            }
            /* diag elements */
            {
                DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                diff_tensor_calc_ewald_beenakker_diag( i, tD );
                insert9( i, i, _D, tD, LD );
            }
        }
    }
}

void init_tensor()
{
    preO = 0.75 * inv_L[0];
    ewald_rr = ewald_real;
    rr_sq = ewald_real * ewald_real;
    if( ewald_method == EWALD_METHOD_BEENAKKER )
    {
        ewald_alpha = SQRTPI * inv_L[0];
    }
    alpha2 = ewald_alpha*ewald_alpha;
    alpha3 = ewald_alpha*alpha2;
    inv_alpha2 = 1.0f / alpha2;
    M_PI_ALPHA_SQ = -M_PIF * M_PIF / alpha2;
}

void compute_single_D( DOUBLE* _D )
{
    compute_D( _D );
}
