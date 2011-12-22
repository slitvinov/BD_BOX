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

#include <math.h>
#include <memory.h>

#include "../err.h"
#include "../math_help.h"
#include "../myblas.h"
#include "../rand_move.h"

#if USE_MPI
#include "mpi.h"
#include "diff_alg/cholesky_mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

void compute_r( DOUBLE* _D )
{
    INT i;
    DOUBLE a;

    gauss_rand( X, size * DIMS0 );
    switch ( hydro )
    {
        case DIFF_NONE:
        {
            /*We ignore _D*/
            for ( i = 0; i < size * DIMS0; ++i )
            {
                Rs[0][ i ] = X[ i ] * SQRTD( diag_D[i] * 2 * curr_dt );
            }
            break;
        }
        case DIFF_CHOLESKY:
        {
            a = SQRTD( 2 * curr_dt );
            for ( i = 0; i < size * DIMS0; ++i )
            {
                X[i] *= a;
            }
            memset( Rs[0], 0, sizeof (DOUBLE) * DIMS1 * size );
#if USE_MPI
            compute_r_cholesky_mpi( _D );
            break;
#endif
            if ( !curr_iter || save_iter != curr_iter )
            {
                if ( my_dpotrf( _D, size * DIMS0, size * DIMS0 ) )
                    UNERR( "Cholesky decomposition error")
                save_iter = curr_iter;
                memcpy( S[ save_iter & 1 ], _D, sizeof (DOUBLE) * DIMS0 * size * DIMS0 * size );
            }
            my_dtrmv( S[ save_iter & 1 ], X, Rs[0], size*DIMS0, size * DIMS0 );
            break;
        }
        case DIFF_GEYER:
        {
            DOUBLE beta = 0.0;
            INT j;
            for ( i = 0; i < g_threads; ++i )
            {
                memset( Rs[i], 0, sizeof (DOUBLE) * DIMS1 * size );
                memset( Rf[i], 0, sizeof (DOUBLE) * DIMS1 * size );
            }
            for ( i = 0; i < size * DIMS0; ++i )
            {
                X[ i ] *= SQRTD( 2 * curr_dt * inv_diag_D[i] );
            }
            memset( epsilon, 0, sizeof (DOUBLE) * g_threads );
            if( geyer_on_the_fly )
            {
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
                for ( i = 0; i < size; ++i )
                {
                    int k, j, l;
                    INT tid = 0;
                    DOUBLE* C;
                    DOUBLE* Feff;
                    DOUBLE eps = 0.0;
                    DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#ifdef _OPENMP
                    tid = omp_get_thread_num( );
#endif
                    C = Rf[tid];
                    Feff = Rs[tid];
                    for ( k = g_id; k < i; k += g_numprocs )
                    {
                        compute_Dij( i , k, tD );
                        for ( j = 0; j < 3; ++j )
                        {
                            for ( l = 0; l < 3; ++l )
                            {
                                C[i * DIMS0 + j] += tD[ j * DIMS0 + l ] * tD[ j * DIMS0 + l ] * inv_diag_D[k * DIMS0 + l]; /* Rf - C[i] */
                                C[k * DIMS0 + l] += tD[ j * DIMS0 + l ] * tD[ j * DIMS0 + l ] * inv_diag_D[i * DIMS0 + j];
                                Feff[i * DIMS0 + j] += X[k * DIMS0 + l] * tD[ j * DIMS0 + l ];
                                Feff[k * DIMS0 + l] += X[i * DIMS0 + j] * tD[ j * DIMS0 + l ];
                                eps += tD[ j * DIMS0 + l ];
                            }
                        }
                    }
                    epsilon[tid] += eps;
                }
            }
            else
            {
                int LD = size*DIMS0;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
                for ( i = 0; i < size*DIMS0; ++i )
                {
                    int k;
                    INT tid = 0;
                    DOUBLE* C;
                    DOUBLE* Feff;
                    DOUBLE eps = 0.0;
#ifdef _OPENMP
                    tid = omp_get_thread_num( );
#endif
                    C = Rf[tid];
                    Feff = Rs[tid];
                    for ( k = g_id; k < i; k += g_numprocs )
                    {
                        C[i] += _D[ i * LD + k ] * _D[ i * LD + k ] * inv_diag_D[k];
                        C[k] += _D[ i * LD + k ] * _D[ i * LD + k ] * inv_diag_D[i];
                        Feff[i] += X[k] * _D[ i * LD + k ];
                        Feff[k] += X[i] * _D[ i * LD + k ];
                        eps += _D[ i * LD + k ];
                    }
                    epsilon[tid] += eps;
                }
            }
            /* sum after parallel computations*/
            for ( j = 1; j < g_threads; ++j )
            {
                for ( i = 0; i < size * DIMS0; ++i )
                {
                    Rf[0][i] += Rf[j][i];
                }
            }
            for ( i = 1; i < g_threads; ++i )
            {
                for ( j = 0; j < size * DIMS0; ++j )
                {
                    Rs[0][j] += Rs[i][j];
                }
            }
            for ( i = 1; i < g_threads; ++i )
            {
                epsilon[0] += epsilon[i];
            }
#if USE_MPI
            Rf[0][size * DIMS0] = epsilon[0];
            MPI_Reduce( Rf[0], Rf[1], size * DIMS0 + 1, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce( Rs[0], Rs[1], size * DIMS0, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
            if ( g_id == 0 )
            {
                /* copy 'recv buffor' to 'send buffor' */
                memcpy( Rf[0], Rf[1], sizeof (DOUBLE)*(size * DIMS0 + 1) );
                memcpy( Rs[0], Rs[1], sizeof (DOUBLE)*(size * DIMS0) );
            }
            MPI_Bcast( Rf[0], size * DIMS0 + 1, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
            MPI_Bcast( Rs[0], size * DIMS0, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
            epsilon[0] = Rf[0][size * DIMS0];
#endif
            /* compute epsilon */
            //epsilon[0] /= (size * DIMS0 * size * DIMS0 - size * DIMS0) / 2;
            epsilon[0] /= size * DIMS0;
            epsilon[0] /= (size * DIMS0 - 1);
            epsilon[0] *= 2;
            if ( epsilon[0] > 1.0 )
            {
                epsilon[0] = 1.0;
                warning( "epsilon truncated to 1.0", __FILE__, __LINE__ );
            }
            a = (3 * size - 1) * epsilon[0] * epsilon[0]-(3 * size - 2) * epsilon[0];
            if ( ABSD( a ) < 0.0000001f )/*we assume no correlation*/
            {
                beta = 0.5;
                warning( "Epsilon == 0.0", __FILE__, __LINE__ );
            }
            else
            {
                beta = (1 - SQRTD( 1 - a )) / a;
            }

            for ( i = 0; i < size * DIMS0; ++i )
            {
                Rf[0][i] *= beta * beta * inv_diag_D[i];
                Rf[0][i] = 1 / (1 + Rf[0][i]);
                Rf[0][i] = SQRTD( Rf[0][i] );
            }
            for ( i = 0; i < size * DIMS0; ++i )
            {
                Rs[0][i] = Rf[0][i]*(beta * Rs[0][i] + diag_D[ i ] * X[i]);
            }
            break;
        }
        case DIFF_CHEBYSHEV:
        {
            UNERR("Unimplemented CHEBYSHEV method")
        }
    }
}
