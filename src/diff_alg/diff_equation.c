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

#include "../input.h"
#include "../main_loop.h"
#include "../math_help.h"
#include "../myblas.h"
#include "../rand_move.h"
#include "../trans.h"

#include "../potentials/calc_func.h"
#include "../potentials/LJ.h"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

void compute_Rf_geyer_in_fly( DOUBLE** _Rf )
{
    INT i;
    DOUBLE pre = curr_dt / (kB * T);
    for ( i = 0; i < g_threads; ++i )
    {
        memset( _Rf[i], 0, sizeof (DOUBLE) * DIMS0 * size );
    }

#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
    for ( i = 0; i < size; ++i )
    {
        INT k, j, l;
        INT tid = 0;
        DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#ifdef _OPENMP
        tid = omp_get_thread_num( );
#endif
        for ( k = g_id; k < i; k += g_numprocs )
        {
            compute_Dij( i , k, tD );
            for ( j = 0; j < 3; ++j )
            {
                for ( l = 0; l < 3; ++l )
                {
                    _Rf[tid][ i * DIMS0 + j ] += F[0][ k * DIMS0 + l ] * tD[ j * DIMS0 + l ];/* symetric */
                    _Rf[tid][ k * DIMS0 + l ] += F[0][ i * DIMS0 + j ] * tD[ j * DIMS0 + l ];/* tD[j*3+l]==tD[l*3+j] */
                }
            }
        }
        if ( g_id == 0 )
        {
            compute_Dij( i , i, tD );
            _Rf[tid][ i * DIMS0    ] += F[0][ i * DIMS0    ] * tD[0];
            _Rf[tid][ i * DIMS0 + 1] += F[0][ i * DIMS0 + 1] * tD[4];
            _Rf[tid][ i * DIMS0 + 2] += F[0][ i * DIMS0 + 2] * tD[8];
        }
    }
    /* Reduce forces */
    for ( i = 1; i < g_threads; ++i )
    {
        INT j;
        for ( j = 0; j < size * DIMS0; ++j )
        {
            _Rf[0][j] += _Rf[i][j];
        }
    }

    for ( i = 0; i < size * DIMS0; ++i )
    {
        _Rf[0][ i ] *= pre;
    }
#if USE_MPI
    MPI_Reduce( _Rf[0], _Rf[1], size*DIMS0, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    if ( g_id == 0 )
    {
        memcpy( _Rf[0], _Rf[1], sizeof (DOUBLE)*(size * DIMS0) );
    }
    MPI_Bcast( _Rf[0], size*DIMS0, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
#endif
}

void compute_Rf( DOUBLE** _Rf, DOUBLE* _D )
{
    INT i;
    DOUBLE pre = curr_dt / (kB * T);
    INT LD = size*DIMS0;
#if USE_MPI
    if ( hydro == DIFF_CHOLESKY )
    {
        compute_Rf_cholesky_mpi( _Rf, _D );
        return;
    }
#endif
    if( geyer_on_the_fly )
    {
        compute_Rf_geyer_in_fly( _Rf );
        return;
    }
    compute_D( _D );

    if ( !hydro )
    {
        for ( i = 0; i < size; ++i )
        {
            _Rf[0][ i * DIMS0 ] = pre * diag_D[ i * DIMS0 ] * F[0][ i * DIMS0 ];
            _Rf[0][ i * DIMS0 + 1 ] = pre * diag_D[ i * DIMS0 + 1 ] * F[0][ i * DIMS0 + 1 ];
            _Rf[0][ i * DIMS0 + 2 ] = pre * diag_D[ i * DIMS0 + 2 ] * F[0][ i * DIMS0 + 2 ];
        }
        return;
    }
    for ( i = 0; i < g_threads; ++i )
    {
        memset( _Rf[i], 0, sizeof (DOUBLE) * DIMS0 * size );
    }
#if USE_LAPACK && !USE_MPI
    {
        DOUBLE dzero = 0;
        DOUBLE done = 1;
        int ione = 1;
        const char* uplo = "U";
        SYMVD( uplo, &LD, &done, _D, &LD, F[0], &ione, &dzero, _Rf[0], &ione );
    }
#else

#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
    for ( i = 0; i < size; ++i )
    {
        INT k, j, l;
        INT tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num( );
#endif
        for ( j = 0; j < 3; ++j )
        {
            for ( k = g_id; k < i; k += g_numprocs )
            {
                for ( l = 0; l < 3; ++l )
                {
                    _Rf[tid][ i * DIMS0 + j ] += F[0][ k * DIMS0 + l ] * _D[ (i * DIMS0 + j) * LD + k * DIMS0 + l ];
                    _Rf[tid][ k * DIMS0 + l ] += F[0][ i * DIMS0 + j ] * _D[ (i * DIMS0 + j) * LD + k * DIMS0 + l ];
                }
            }
        }
        if ( g_id == 0 )
        {
            _Rf[tid][ i * DIMS0 ] += F[0][ i * DIMS0 ] * _D[ i * DIMS0 * LD + i * DIMS0 ];
            _Rf[tid][ i * DIMS0 + 1] += F[0][ i * DIMS0 + 1] * _D[ i * DIMS0 * LD + i * DIMS0 ];
            _Rf[tid][ i * DIMS0 + 2] += F[0][ i * DIMS0 + 2] * _D[ i * DIMS0 * LD + i * DIMS0 ];
        }
    }
    /* Reduce forces */
    for ( i = 1; i < g_threads; ++i )
    {
        INT j;
        for ( j = 0; j < size * DIMS0; ++j )
        {
            _Rf[0][j] += _Rf[i][j];
        }
    }
#endif

    for ( i = 0; i < size * DIMS0; ++i )
    {
        _Rf[0][ i ] *= pre;
    }

#if USE_MPI
    MPI_Reduce( _Rf[0], _Rf[1], size*DIMS0, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    if ( g_id == 0 )
    {
        memcpy( _Rf[0], _Rf[1], sizeof (DOUBLE)*(size * DIMS0) );
    }
    MPI_Bcast( _Rf[0], size*DIMS0, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
#endif
}

void compute_hydro_ermak( DOUBLE* coord_out )
{
    INT i;
    compute_Rf( Rf, Ds );

    if ( (hydro == DIFF_CHOLESKY && algorithm > DIFF_ALG_ERMAK_NEWTON) )
    {
#ifdef _OPENMP
        if ( size < 50 )
#endif
        {
            memcpy( D, Ds, sizeof (DOUBLE) * DIMS0 * DIMS0 * size * size );
        }
#ifdef _OPENMP
        else
        {
            INT LD = size*DIMS0;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
            for ( i = 0; i < LD; ++i )
            {
                INT j;
                for ( j = 0; j <= i; ++j )
                {
                    D[ i * LD + j ] = Ds[ i * LD + j ];
                }
            }
        }
#endif
    }

    memcpy( Rf0, Rf[0], sizeof (DOUBLE) * DIMS0 * size );

    compute_r( Ds );

    memcpy( Rf[0], Rf0, sizeof (DOUBLE) * DIMS0 * size );

    if ( !vel_grad_tensor )
    {
        for ( i = 0; i < size; ++i )
        {
            coord_out[ i * DIMS1 ] = Rf[0][ i * DIMS0 ] + Rs[0][ i * DIMS0 ] + coord[ i * DIMS1 ];
            coord_out[ i * DIMS1 + 1] = Rf[0][ i * DIMS0 + 1] + Rs[0][ i * DIMS0 + 1] + coord[ i * DIMS1 + 1];
            coord_out[ i * DIMS1 + 2] = Rf[0][ i * DIMS0 + 2] + Rs[0][ i * DIMS0 + 2] + coord[ i * DIMS1 + 2];
            coord_out[ i * DIMS1 + 3] = coord[ i * DIMS1 + 3];
            trans_to_box( coord_out + DIMS1 * i );
        }
    }
    else
    {
        DOUBLE* G = vel_grad_tensor;
        for ( i = 0; i < size; ++i )
        {
            DOUBLE tmp0 = curr_dt * (G[0] * coord[ i * DIMS1 ] + G[1] * coord[ i * DIMS1 + 1] + G[2] * coord[ i * DIMS1 + 2]);
            DOUBLE tmp1 = curr_dt * (G[3] * coord[ i * DIMS1 ] + G[4] * coord[ i * DIMS1 + 1] + G[5] * coord[ i * DIMS1 + 2]);
            DOUBLE tmp2 = curr_dt * (G[6] * coord[ i * DIMS1 ] + G[7] * coord[ i * DIMS1 + 1] + G[8] * coord[ i * DIMS1 + 2]);

            coord_out[ i * DIMS1 ] = Rf[0][ i * DIMS0 ] + Rs[0][ i * DIMS0 ] + coord[ i * DIMS1 ] + tmp0;
            coord_out[ i * DIMS1 + 1] = Rf[0][ i * DIMS0 + 1] + Rs[0][ i * DIMS0 + 1] + coord[ i * DIMS1 + 1] + tmp1;
            coord_out[ i * DIMS1 + 2] = Rf[0][ i * DIMS0 + 2] + Rs[0][ i * DIMS0 + 2] + coord[ i * DIMS1 + 2] + tmp2;
            coord_out[ i * DIMS1 + 3] = coord[ i * DIMS1 + 3];

            FS[ i * DIMS1 ] = tmp0 + curr_dt * (G[0] * coord_out[ i * DIMS1 ] + G[1] * coord_out[ i * DIMS1 + 1] + G[2] * coord_out[ i * DIMS1 + 2]);
            FS[ i * DIMS1 + 1] = tmp1 + curr_dt * (G[3] * coord_out[ i * DIMS1 ] + G[4] * coord_out[ i * DIMS1 + 1] + G[5] * coord_out[ i * DIMS1 + 2]);
            FS[ i * DIMS1 + 2] = tmp2 + curr_dt * (G[6] * coord_out[ i * DIMS1 ] + G[7] * coord_out[ i * DIMS1 + 1] + G[8] * coord_out[ i * DIMS1 + 2]);

            trans_to_box( coord_out + DIMS1 * i );
        }
    }
}

INT compute_hydro_IGT_second_phase( DOUBLE* _E )
{
    INT i;

    compute_r( D ); /* in Rf0 random move */

    if ( !vel_grad_tensor )
    {
        for ( i = 0; i < size; ++i )
        {
            coord[ i * DIMS1 ] += Rf0[ i * DIMS0 ] + Rs[0][ i * DIMS0 ];
            coord[ i * DIMS1 + 1] += Rf0[ i * DIMS0 + 1] + Rs[0][ i * DIMS0 + 1];
            coord[ i * DIMS1 + 2] += Rf0[ i * DIMS0 + 2] + Rs[0][ i * DIMS0 + 2];
            trans_to_box( coord + DIMS1 * i );
        }
    }
    else
    {
        for ( i = 0; i < size; ++i )
        {
            coord[ i * DIMS1 ] += Rf0[ i * DIMS0 ] + Rs[0][ i * DIMS0 ] + 0.5 * FS[ i * DIMS1 ];
            coord[ i * DIMS1 + 1] += Rf0[ i * DIMS0 + 1] + Rs[0][ i * DIMS0 + 1] + 0.5 * FS[ i * DIMS1 + 1];
            coord[ i * DIMS1 + 2] += Rf0[ i * DIMS0 + 2] + Rs[0][ i * DIMS0 + 2] + 0.5 * FS[ i * DIMS1 + 2];
            trans_to_box( coord + DIMS1 * i );
        }
    }
    return compute_forces( _E, curr_time + curr_dt ) << 1;
}

INT compute_hydro_IGT( DOUBLE* _E )
{
    INT i;
    DOUBLE* tmp_coord;

    compute_hydro_ermak( Rf0 ); /* in Rf determistic move, Rf0 coord for predictor step */

    tmp_coord = coord;
    coord = Rf0;

    if ( compute_forces( NULL, curr_time + curr_dt ) )
    {
        coord = tmp_coord;
        return 1;
    }

    compute_Rf( Rs, (hydro == DIFF_CHOLESKY) ? Ds : D );
    coord = tmp_coord;

    for ( i = 0; i < size * DIMS0; ++i )
    {
        Rf0[i] = 0.5 * (Rf[0][i] + Rs[0][i]); /*determistic move*/
    }

    while ( hydro )
    {
        INT LD = size*DIMS0;
#if USE_MPI
        if ( hydro == DIFF_CHOLESKY )
        {
            sum_tensors_cholesky_mpi( );
            break;
        }
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
        for ( i = 0; i < size; ++i )
        {
            INT j, k, l;
            for ( j = 0; j < 3; ++j )
            {
                for ( k = g_id; k < i; k += g_numprocs )
                {
                    for ( l = 0; l < 3; ++l )
                    {
                        D[ (i * DIMS0 + j) * LD + k * DIMS0 + l ] =
                            (D [ (i * DIMS0 + j) * LD + k * DIMS0 + l ] +
                            Ds[ (i * DIMS0 + j) * LD + k * DIMS0 + l ]) / 2;
                    }
                }
            }
        }
        break; /*only once*/
    }
    curr_iter++;
    return compute_hydro_IGT_second_phase( _E );
}
