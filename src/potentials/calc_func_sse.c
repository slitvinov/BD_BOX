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

#if USE_SSE

#include "calc_func.h"

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>

#include "../input.h"
#include "../math_help.h"
#include "../trans.h"

#include "angle.h"
#include "angle_cos.h"
#include "bond.h"
#include "dihe.h"
#include "dihe_angle.h"
#include "electro.h"
#include "electro_ext.h"
#include "LJ.h"

#include "sse_mathfun.h"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

DOUBLE* E;
DOUBLE** ES;
DOUBLE** F; /*Forces*/
INT* repeat;

MAKE_STR_IN( DOUBLE, bond_lj_scale, 1.0, "scaling factor for Lennard–Jones interactions between bonded partners (1–2)" )
MAKE_STR_IN( DOUBLE, bond_c_scale, 1.0, "scaling factor for electrostatic interactions between bonded partners (1–2)" )
MAKE_STR_IN( YESNO, check_overlap, 1, "whether to detect overlaps after each BD step" )

DOUBLE* h_coord; /*host coord + LJ or Q*/
DOUBLE** h_f; /* forces */

inline
__m128d exp_pd( __m128d in )
{
    __m128d ret;
    __m128 s1, s2;
    /* mixed arithmetic method */
    s1 = _mm_cvtpd_ps( in );
    s2 = exp_ps( s1 );
    ret = _mm_cvtps_pd( s2 );
    return ret;
}

INT compute_forces( DOUBLE* _E, DOUBLE curr_time )
{
    int i;
#ifdef _OPENMP
    INT tid = 0;
#endif
    INT rep = 0;

    if ( !iscorrect( ) ) return 1;
    memset( repeat, 0, sizeof (INT) * g_threads );
    for ( i = 0; i < g_threads; ++i )
    {
        memset( F[i], 0, sizeof (DOUBLE) * DIMS0 * size );
        memset( h_f[i], 0 , sizeof (DOUBLE) * (size * DIMS0 + 16) );
    }
    if ( _E )
    {
        memset( _E, 0, sizeof (DOUBLE) * (ENERGY_LAST + 1) );
        for ( i = 0; i < g_threads; ++i )
        {
            memset( ES[i], 0, sizeof (DOUBLE) * (ENERGY_LAST + 1) );
        }
    }

    if ( elec || alpha_lj != 0.0 ) /* if elec or lj turn on*/
    {
        if ( elec && alpha_lj == 0.0 )
        {
            int al_size = size & (~(SSECOUNT - 1));
            M128* tab;
            M128 m_kappac = MMSET1( kappa_c );
            M128 m_mkappac = MMSET1( -kappa_c );
            M128 m_ones = MMSET1( 1 );
            M128 gE = MMSETZERO( );
            DOUBLE* hc;
            DOUBLE sE[SSECOUNT];
            INT k;
            memset( F[0], 0, sizeof (DOUBLE) * DIMS0 * size );
            for ( i = al_size; i < size; ++i )
            {
                INT j;
                for ( j = i + 1 + g_id; j < size && !rep; j += g_numprocs )
                {
                    DOUBLE r_norm;
                    DOUBLE r_sq;
                    DOUBLE e[3];
                    get_v2_norms( i, j, e, &r_norm, &r_sq );
                    rep = rep || electro_force( r_norm, r_sq, e, F[0] + i*DIMS0, F[0] + j*DIMS0, Q[i], Q[j], _E ? (_E + ENERGY_COULOMB) : _E );
                }
            }
            for ( i = 0; i < al_size; ++i )
            {
                INT j;
                INT end = (i + SSECOUNT) - (i & (SSECOUNT - 1));
                if ( end > size ) end = size;
                for ( j = i + g_id + 1; j < end && !rep; j += g_numprocs )
                {
                    DOUBLE r_norm;
                    DOUBLE r_sq;
                    DOUBLE e[3];
                    get_v2_norms( i, j, e, &r_norm, &r_sq );
                    rep = rep || electro_force( r_norm, r_sq, e, F[0] + i*DIMS0, F[0] + j*DIMS0, Q[i], Q[j], _E ? (_E + ENERGY_COULOMB) : _E );
                }
            }

            hc = (DOUBLE*) (16 * (((size_t) h_coord + 15) / 16)); /*make alignment*/
            for ( i = 0; i < al_size / SSECOUNT; ++i )
            {
                INT j = 0;
                for ( j = 0; j < SSECOUNT; ++j )
                {
                    hc[SSECOUNT * i * DIMS1 + j] = coord[ (SSECOUNT * i + j) * DIMS1 ];
                    hc[SSECOUNT * (i * DIMS1 + 1) + j] = coord[ (SSECOUNT * i + j) * DIMS1 + 1];
                    hc[SSECOUNT * (i * DIMS1 + 2) + j] = coord[ (SSECOUNT * i + j) * DIMS1 + 2];
                    hc[SSECOUNT * (i * DIMS1 + 3) + j] = Q[SSECOUNT * i + j];
                }
            }
            tab = (M128*) hc;
            al_size /= SSECOUNT;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic ) private(rep)
#endif
            for ( i = 0; i < al_size; ++i )
            {
                INT j;
                M128 x, y, z, q;
                M128 Fxi, Fyi, Fzi;
                M128* out_F = (M128*) (16 * (((size_t) h_f[0] + 15) / 16));
                M128 s10, s11;
                M128 mE;
#ifdef _OPENMP
                out_F = (M128*) (16 * (((size_t) h_f[omp_get_thread_num( )] + 15) / 16));
#endif
                mE  = MMSETZERO( );
                Fxi = MMSETZERO( );
                Fyi = MMSETZERO( );
                Fzi = MMSETZERO( );

                x = tab[ i * DIMS1 ];
                y = tab[ i * DIMS1 + 1];
                z = tab[ i * DIMS1 + 2];
                q = tab[ i * DIMS1 + 3];

                for ( j = 0; j < i; ++j )
                {
                    M128 xj, yj, zj, qj;
                    M128 rx, ry, rz, Fxj, Fyj, Fzj;
                    M128 s8, s9;
                    INT k;
                    Fxj = MMSETZERO( );
                    Fyj = MMSETZERO( );
                    Fzj = MMSETZERO( );
                    xj = tab[ j * DIMS1    ];
                    yj = tab[ j * DIMS1 + 1];
                    zj = tab[ j * DIMS1 + 2];
                    qj = tab[ j * DIMS1 + 3];

                    for ( k = 0; k < SSECOUNT; ++k )
                    {
                        M128 rxx, ryy, rzz;
                        M128 s1, s2, s3, s4, s5, s6, s7, r3, s12, invr, res;
                        M128 s14, s16;
                        M128 kappac_rnorm;
                        M128 mkappac_rnorm;

                        rx = MMSUB( x, xj );
                        ry = MMSUB( y, yj );
                        rz = MMSUB( z, zj );

                        rxx = MMMUL( rx, rx );
                        ryy = MMMUL( ry, ry );
                        rzz = MMMUL( rz, rz );

                        s1 = MMADD( rxx, ryy );
                        s2 = MMADD( s1, rzz );
                        s3 = MMSQRT( s2 );
                        r3 = MMMUL( s3, s2 );
                        invr = MMDIV( m_ones, r3 );

                        kappac_rnorm = MMMUL( m_kappac, s3 );
                        mkappac_rnorm = MMMUL( m_mkappac, s3 );
                        s4 = MMADD( kappac_rnorm, m_ones );
                        s5 = MEXP( mkappac_rnorm );
                        s6 = MMMUL( s5, qj );
                        s7 = MMMUL( s6, q );

                        mE = MMADD( MMDIV( s7, s3 ), mE );
                        res = MMMUL( s4, s7 );

                        res = MMMUL( res, invr );
                        s12 = MMMUL( res, rx );
                        Fxi = MMADD( Fxi, s12 );
                        Fxj = MMSUB( Fxj, s12 );

                        s14 = MMMUL( res, ry );
                        Fyi = MMADD( Fyi, s14 );
                        Fyj = MMSUB( Fyj, s14 );

                        s16 = MMMUL( res, rz );
                        Fzi = MMADD( Fzi, s16 );
                        Fzj = MMSUB( Fzj, s16 );

                        /* Rotate to get all pairs */
                        xj = MMSHUFFLE( xj, xj, SSESHUFFLE );
                        yj = MMSHUFFLE( yj, yj, SSESHUFFLE );
                        zj = MMSHUFFLE( zj, zj, SSESHUFFLE );
                        qj = MMSHUFFLE( qj, qj, SSESHUFFLE );
                        Fxj = MMSHUFFLE( Fxj, Fxj, SSESHUFFLE );
                        Fyj = MMSHUFFLE( Fyj, Fyj, SSESHUFFLE );
                        Fzj = MMSHUFFLE( Fzj, Fzj, SSESHUFFLE );
                    }
                    s8 = out_F[3 * j];
                    s9 = MMADD( s8, Fxj );
                    out_F[3 * j] = s9;

                    s8 = out_F[3 * j + 1];
                    s9 = MMADD( s8, Fyj );
                    out_F[3 * j + 1] = s9;

                    s8 = out_F[3 * j + 2];
                    s9 = MMADD( s8, Fzj );
                    out_F[3 * j + 2] = s9;
                }
                s10 = out_F[3 * i];
                s11 = MMADD( s10, Fxi );
                out_F[3 * i] = s11;

                s10 = out_F[3 * i + 1];
                s11 = MMADD( s10, Fyi );
                out_F[3 * i + 1] = s11;

                s10 = out_F[3 * i + 2];
                s11 = MMADD( s10, Fzi );
                out_F[3 * i + 2] = s11;
#ifdef _OPENMP
#pragma omp critical
                {
#endif
                    gE = MMADD( gE, mE );
#ifdef _OPENMP
                }
#endif
            }
            if ( _E )
            {
                MMSTORE( sE, gE );
                for ( i = 1; i < SSECOUNT; ++i )
                {
                    sE[0] += sE[i];
                }
                _E[ ENERGY_COULOMB ] += (gamma_c / epsilon_c) * sE[0];
            }
            for ( k = 0; k < g_threads; ++k )
            {
                DOUBLE* hF = (DOUBLE*) (16 * (((size_t) h_f[k] + 15) / 16));
                for ( i = 0; i < al_size; ++i )
                {
                    INT j = 0;
                    for ( j = 0; j < SSECOUNT; ++j )
                    {
                        F[0][ (SSECOUNT * i + j) * DIMS0 ] += (gamma_c / epsilon_c) * hF[SSECOUNT * i * DIMS0 + j];
                        F[0][ (SSECOUNT * i + j) * DIMS0 + 1] += (gamma_c / epsilon_c) * hF[SSECOUNT * (i * DIMS0 + 1) + j];
                        F[0][ (SSECOUNT * i + j) * DIMS0 + 2] += (gamma_c / epsilon_c) * hF[SSECOUNT * (i * DIMS0 + 2) + j];
                    }
                }
            }
        }
        else
        {
            UNERR("Can't use LJ and electro in SSE mode");
        }
    }
    /* always check overlap */
    {
        int al_size = size & (~(SSECOUNT - 1));
        M128* tab;
        DOUBLE* hc;
        rep = 0;
        if ( al_size & 1 ) al_size--;
        if ( al_size < 0 ) al_size = 0;
        for ( i = al_size; i < size; ++i )
        {
            INT j;
            for ( j = i + 1 + g_id; j < size && !rep; j += g_numprocs )
            {
                DOUBLE r_sq;
                get_v2_sq( i, j, &r_sq );
                rep = rep || LJ_force_check( r_sq, LJ + 2 * i, LJ + 2 * j );
            }
        }
        for ( i = 0; i < al_size; ++i )
        {
            INT j;
            INT end = (i + SSECOUNT) - (i & (SSECOUNT - 1));
            for ( j = i + g_id + 1; j < end && !rep; j += g_numprocs )
            {
                DOUBLE r_sq;
                get_v2_sq( i, j, &r_sq );
                rep = rep || LJ_force_check( r_sq, LJ + 2 * i, LJ + 2 * j );
            }
        }
        hc = (DOUBLE*) (16 * (((size_t) h_coord + 15) / 16)); /*make alignment*/
        for ( i = 0; i < al_size / SSECOUNT; ++i )
        {
            INT j = 0;
            for ( j = 0; j < SSECOUNT; ++j )
            {
                hc[SSECOUNT * i * DIMS1 + j] = coord[ (SSECOUNT * i + j) * DIMS1 ];
                hc[SSECOUNT * (i * DIMS1 + 1) + j] = coord[ (SSECOUNT * i + j) * DIMS1 + 1];
                hc[SSECOUNT * (i * DIMS1 + 2) + j] = coord[ (SSECOUNT * i + j) * DIMS1 + 2];
                hc[SSECOUNT * (i * DIMS1 + 3) + j] = LJ[(SSECOUNT * i + j)*2 + LJ_SIGMA_OFF];
            }
        }
        tab = (M128*) hc;
        al_size /= SSECOUNT;
#ifdef _OPENMP
#pragma omp parallel reduction(||:rep)
#endif
{
    M128 res_sv = MMSETZERO();
#ifdef _MSC_VER
    __declspec(align(16)) DOUBLE sv[SSECOUNT];
#else
    __attribute__((aligned(16))) DOUBLE sv[SSECOUNT];
#endif
    int h;
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
    for ( i = 0; i < al_size; i+=2 )
    {
        INT j,k;
        M128 X[2],Y[2],Z[2],W[2];

        X[0] = tab[ i * DIMS1];
        Y[0] = tab[ i * DIMS1 + 1];
        Z[0] = tab[ i * DIMS1 + 2];
        W[0] = tab[ i * DIMS1 + 3];
        X[1] = tab[ i * DIMS1 + 4];
        Y[1] = tab[ i * DIMS1 + 5];
        Z[1] = tab[ i * DIMS1 + 6];
        W[1] = tab[ i * DIMS1 + 7];
        for ( j = 0; j < i; j+=2 )
        {
            M128 XJ[2], YJ[2], ZJ[2], WJ[2];
            XJ[0] = tab[ j * DIMS1];
            YJ[0] = tab[ j * DIMS1 + 1];
            ZJ[0] = tab[ j * DIMS1 + 2];
            WJ[0] = tab[ j * DIMS1 + 3];
            XJ[1] = tab[ j * DIMS1 + 4];
            YJ[1] = tab[ j * DIMS1 + 5];
            ZJ[1] = tab[ j * DIMS1 + 6];
            WJ[1] = tab[ j * DIMS1 + 7];
            for ( k = 0; k < SSECOUNT; ++k )
            {
                int l,m;
                M128 rx, ry, rz, rw;
                M128 rxx, ryy, rzz;
                M128 s1, s2, res;
                
                for ( l = 0; l < 2; ++l )
                for ( m = 0; m < 2; ++m )
                {
                    rx = MMSUB( X[l], XJ[m] );
                    ry = MMSUB( Y[l], YJ[m] );
                    rz = MMSUB( Z[l], ZJ[m] );

                    rw = MMADD( W[l], WJ[m] );
                    rw = MMMUL( rw, rw );

                    rxx = MMMUL( rx, rx );
                    ryy = MMMUL( ry, ry );
                    rzz = MMMUL( rz, rz );

                    s1 = MMADD( rxx, ryy );
                    s2 = MMADD( s1, rzz );

                    res = MMCMPLT( s2, rw );

                    res_sv = MMADD(res_sv,res);
                }

                XJ[0] = MMSHUFFLE( XJ[0], XJ[0], SSESHUFFLE );
                YJ[0] = MMSHUFFLE( YJ[0], YJ[0], SSESHUFFLE );
                ZJ[0] = MMSHUFFLE( ZJ[0], ZJ[0], SSESHUFFLE );
                WJ[0] = MMSHUFFLE( WJ[0], WJ[0], SSESHUFFLE );
                
                XJ[1] = MMSHUFFLE( XJ[1], XJ[1], SSESHUFFLE );
                YJ[1] = MMSHUFFLE( YJ[1], YJ[1], SSESHUFFLE );
                ZJ[1] = MMSHUFFLE( ZJ[1], ZJ[1], SSESHUFFLE );
                WJ[1] = MMSHUFFLE( WJ[1], WJ[1], SSESHUFFLE );
            }
        }
        for ( k = 0; k < SSECOUNT; ++k )
        {
            M128 rx, ry, rz, rw;
            M128 rxx, ryy, rzz;
            M128 s1, s2, res;

            rx = MMSUB( X[0], X[1] );
            ry = MMSUB( Y[0], Y[1] );
            rz = MMSUB( Z[0], Z[1] );

            rw = MMADD( W[0], W[1] );
            rw = MMMUL( rw, rw );

            rxx = MMMUL( rx, rx );
            ryy = MMMUL( ry, ry );
            rzz = MMMUL( rz, rz );

            s1 = MMADD( rxx, ryy );
            s2 = MMADD( s1, rzz );

            res = MMCMPLT( s2, rw );

            res_sv = MMADD(res_sv,res);
            X[1] = MMSHUFFLE( X[1], X[1], SSESHUFFLE );
            Y[1] = MMSHUFFLE( Y[1], Y[1], SSESHUFFLE );
            Z[1] = MMSHUFFLE( Z[1], Z[1], SSESHUFFLE );
            W[1] = MMSHUFFLE( W[1], W[1], SSESHUFFLE );
        }
    }
    MMSTORE( sv, res_sv );
    for ( h = 0; h < SSECOUNT; ++h )
    {
        rep |= (sv[h] != 0);
    }
}
        repeat[0] = (rep != 0);
    }
    if ( _E )
    {
        for ( i = 0; i < g_threads; ++i )
        {
            INT j;
            for ( j = 0; j < ENERGY_LAST; ++j )
            {
                _E[ j ] += ES[i][j];
            }
        }
    }
    for ( i = 1; i < g_threads; ++i )
    {
        INT j;
        for ( j = 0; j < size * DIMS0; ++j )
        {
            F[0][j] += F[i][j];
        }
        repeat[0] = repeat[0] || repeat[i];
    }
#if USE_MPI
    F[0][size * DIMS0] = repeat[0];
    if ( _E )
    {
        for ( i = 1; i <= ENERGY_LAST; ++i )
        {
            F[0][ size * DIMS0 + i ] = _E[ i - 1 ];
        }
        MPI_Reduce( F[0], F[1], size * DIMS0 + ENERGY_LAST + 1, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    }
    else
    {
        MPI_Reduce( F[0], F[1], size * DIMS0 + 1, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    }
    if ( g_id == 0 )
    {
        /* copy 'recv buffor' to 'send send' */
        memcpy( F[0], F[1], sizeof (DOUBLE)*(size * DIMS0 + ENERGY_LAST + 1) );
    }
    MPI_Bcast( F[0], size * DIMS0 + 1, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
    if ( _E && g_id == 0 ) /*Only master saves results*/
    {
        for ( i = 1; i <= ENERGY_LAST; ++i )
        {
            _E[ i - 1 ] = F[0][ size * DIMS0 + i ];
        }
    }
    repeat[0] = F[0][size * DIMS0] != 0.0;
#endif
    return repeat[0];
}

void init_forces( )
{
    INT i;
    repeat = (INT*) malloc( sizeof (INT) * g_threads ); CHMEM(repeat);
    E = (DOUBLE*) malloc( sizeof (DOUBLE) * (ENERGY_LAST + 1) ); CHMEM(E);
#if USE_MPI
    F = (DOUBLE**) malloc( sizeof (DOUBLE*) * (g_threads + 1) ); CHMEM(F);
#else
    F = (DOUBLE**) malloc( sizeof (DOUBLE*) * g_threads ); CHMEM(F);
#endif
    for ( i = 0; i < g_threads; ++i )
    {
        F[i] = (DOUBLE*) malloc( sizeof (DOUBLE) * (DIMS1 * size + ENERGY_LAST + 2) ); /*We will be using it for sending in MPI*/
        CHMEM(F[i]);
    }
#if USE_MPI
    if ( g_id == 0 )
    {
        F[ g_threads ] = (DOUBLE*) malloc( sizeof (DOUBLE) * (DIMS1 * size + ENERGY_LAST + 2) );
        CHMEM(F[ g_threads ]);
    }
    else
    {
        F[ g_threads ] = NULL;
    }
#endif
    ES = (DOUBLE**) malloc( sizeof (DOUBLE*) * g_threads ); CHMEM(ES);
    for ( i = 0; i < g_threads; ++i )
    {
        ES[i] = (DOUBLE*) malloc( sizeof (DOUBLE) * (ENERGY_LAST + 1) );
        CHMEM(ES[i]);
    }
    h_coord = (DOUBLE*) malloc( sizeof (DOUBLE)*(size * DIMS1 + 16) );
    CHMEM(h_coord);
    h_f = (DOUBLE**) malloc( sizeof (DOUBLE*) * g_threads ); CHMEM(h_f);
    for ( i = 0; i < g_threads; ++i )
    {
        h_f[ i ] = (DOUBLE*) malloc( sizeof (DOUBLE) * (size * DIMS0 + 16) );
        memset( h_f[i], 0 , sizeof (DOUBLE) * (size * DIMS0 + 16) );
    }
}

void free_forces( )
{
    INT i;
    for ( i = 0; i < g_threads; ++i )
    {
        free( F[i] );
        free( ES[i] );
        free( h_f[i] );
    }
    free( h_f );
    free( F );
    free( ES );
    free( E );
    free( repeat );
    free( h_coord );
}

#else
void dummy_sse(){}
#endif
