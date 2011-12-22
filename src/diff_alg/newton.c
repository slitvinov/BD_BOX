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

#ifndef USE_CUDA

#include "newton.h"

#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>

#include "diff_tensor.h"

#include "../err.h"
#include "../input.h"
#include "../math_help.h"
#include "../trans.h"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

DOUBLE* velocity = NULL;
INT size_pairs;
INT** pairs = NULL;
DOUBLE s_time;

void comp_time( )
{
    INT i;
#ifdef _OPENMP
#pragma omp parallel
    {
#endif
        INT tid = 0;
        INT* curr_pairs;
        INT j;
        DOUBLE curr_stime = DOUBLE_MAX;
        INT curr_spairs = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num( );
#endif
        curr_pairs = pairs[tid];
#ifdef _OPENMP
#pragma omp for schedule( dynamic ) nowait
#endif
        for ( i = 0; i < size; ++i )
        {
            for ( j = i + 1 + g_id; j < size; j += g_numprocs )
            {
                DOUBLE r[3];
                DOUBLE v[3];
                DOUBLE rv;
                DOUBLE l;
                trans_dist_vec( coord + DIMS1*j, coord + DIMS1*i, r );
                l = LJ[ 2*j ] + LJ[ 2*i ] + NEWTON_OFFSET;
                dist_vec( velocity + DIMS0*j, velocity + DIMS0*i, v );
                rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
                if ( rv < 0 )
                {
                    DOUBLE vv = sq3v( v );
                    DOUBLE delta = rv * rv - vv * (sq3v( r ) - l * l);
                    if ( delta > 0 )
                    {
                        DOUBLE t2;
                        DOUBLE inv_vv = 1 / vv;
                        delta = SQRTD( delta );
                        t2 = -(rv + delta) * inv_vv;

                        if ( t2 < 0 ) t2 = 0;
                        if ( t2 <= curr_stime )
                        {
                            curr_spairs = (t2 < curr_stime) ? 0 : curr_spairs;
                            if ( size_pairs + 1 < 2 * size )
                            {
                                curr_pairs[ curr_spairs ] = i;
                                curr_pairs[ curr_spairs + 1 ] = j;
                                curr_spairs += 2;
                            }
                            curr_stime = t2;
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
#pragma omp master
        {
#endif
            size_pairs = curr_spairs;
            s_time = curr_stime;
#ifdef _OPENMP
        }
#pragma omp barrier
#pragma omp critical
#endif
        if ( tid > 0 && curr_stime <= s_time )
        {
            if ( curr_stime < s_time ) size_pairs = 0;
            for ( i = 0; i < curr_spairs / 2; ++i )
            {
                if ( size_pairs + 1 < 2 * size )
                {
                    pairs[0][ size_pairs ] = pairs[tid][2 * i ];
                    pairs[0][ size_pairs + 1] = pairs[tid][2 * i + 1];
                    size_pairs += 2;
                }
                else
                {
                    break;
                }
            }
            s_time = curr_stime;
        }
#ifdef _OPENMP
    }
#endif
#if USE_MPI
    {
        double tab[2];
        double* recv_tab = NULL;
        tab[0] = curr_stime;
        tab[1] = curr_spairs;
        if ( g_id == 0 )
        {
            recv_tab = (double*) malloc( sizeof (double) * g_numprocs * 2 );
        }
        MPI_Gather( tab, 2, MPI_DOUBLE, recv_tab, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        if ( g_id == 0 )
        {
            for ( i = 1; i < g_numprocs; ++i )
            {
                if ( recv_tab[2 * i] <= tab[0] )
                {
                    if ( recv_tab[2 * i] < tab[0] ) tab[1] = 0;
                    tab[0] = recv_tab[2 * i];
                    tab[1] += recv_tab[2 * i + 1];
                }
            }
        }
        MPI_Bcast( tab, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        s_time = tab[0];
        size_pairs = tab[1] > 2 * size ? 2 * size : tab[1];
        if ( g_id == 0 )
        {
            INT* recv_pairs = (INT*) malloc( sizeof (INT) * tab[1] );
            MPI_Status status;
            INT j;
            if ( curr_stime != tab[0] )
            {
                curr_spairs = 0;
            }
            for ( i = 1; i < g_numprocs; ++i )
            {
                if ( recv_tab[2 * i] == tab[0] )
                {
                    MPI_Recv( recv_pairs, recv_tab[2 * i + 1], MPI_INT, i, 1, MPI_COMM_WORLD, &status );
                    for ( j = 0; j < recv_tab[2 * i + 1]; ++j, ++curr_spairs )
                    {
                        if ( curr_spairs < 2 * size )
                        {
                            pairs[0][ curr_spairs ] = recv_pairs[ curr_spairs ];
                        }
                    }
                }
            }
            if ( size_pairs != curr_spairs ) UNERR("Assertion fault")
            free( recv_pairs );
        }
        else
        {
            if ( tab[ 0 ] == curr_stime ) /*send to root*/
            {
                MPI_Send( pairs[0], curr_spairs, MPI_INT, 0, 1, MPI_COMM_WORLD );
            }
        }
        MPI_Bcast( pairs[0], size_pairs, MPI_INT, 0, MPI_COMM_WORLD );
        free( recv_tab );
    }
#endif
}

INT newton_moves( INT* pcount )
{
    DOUBLE next_move_time;
    DOUBLE elapsed = 0.0;
    INT ret, i, j;
	ret = 0;
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
            for ( j = 0; j < size; ++j )
            {
                coord[ DIMS1 * j ] += velocity[j * DIMS0 ] * next_move_time;
                coord[ DIMS1 * j + 1] += velocity[j * DIMS0 + 1] * next_move_time;
                coord[ DIMS1 * j + 2] += velocity[j * DIMS0 + 2] * next_move_time;
                trans_to_box( coord + DIMS1 * j );
            }
            ret = !iscorrect( );
        }
        if ( ret )
        {
            break;
        }
        for ( i = 0; i < size_pairs / 2; ++i )
        {
            int p1 = pairs[0][ i * 2 ];
            int p2 = pairs[0][ i * 2 + 1];
            DOUBLE r_norm, r_sq;
            DOUBLE e[3], v[3];
            DOUBLE ev;
            DOUBLE m1 = masses[p1];
            DOUBLE m2 = masses[p2];
            get_v2_norms( p2, p1, e, &r_norm, &r_sq );
            dist_vec( velocity + DIMS0*p2, velocity + DIMS0*p1, v );
            ev = e[0] * v[0] + e[1] * v[1] + e[2] * v[2];
            ev *= (1 + e_collision) / (m1 + m2);
            velocity[DIMS0 * p1 ] -= m2 * ev * e[0];
            velocity[DIMS0 * p1 + 1] -= m2 * ev * e[1];
            velocity[DIMS0 * p1 + 2] -= m2 * ev * e[2];

            velocity[DIMS0 * p2 ] += m1 * ev * e[0];
            velocity[DIMS0 * p2 + 1] += m1 * ev * e[1];
            velocity[DIMS0 * p2 + 2] += m1 * ev * e[2];
        }
        ++(*pcount);
        if ( *pcount > move_attempts ) UNERR("Ermak Newton error")
    }/* while */
    return ret;
}

void newton_init_vel( )
{
    INT i;
    DOUBLE inv_t = 1.0 / curr_dt;
    for ( i = 0; i < size; ++i )
    {
        trans_dist_vec( save_coord + DIMS1*i, coord + DIMS1*i, velocity + i * DIMS0 );
        velocity[i * DIMS0 ] *= inv_t;
        velocity[i * DIMS0 + 1] *= inv_t;
        velocity[i * DIMS0 + 2] *= inv_t;
    }
}

void init_newton( )
{
    INT i;
    velocity = (DOUBLE*) malloc( sizeof (DOUBLE) * DIMS1 * size );
    pairs = (INT**) malloc( sizeof (INT*) * g_threads );
    for ( i = 0; i < g_threads; ++i )
    {
        pairs[i] = (INT*) malloc( sizeof (INT) * 2 * size );
    }
}

void free_newton( )
{
    INT i;
    free( velocity );
    for ( i = 0; i < g_threads; ++i )
    {
        free( pairs[i] );
    }
    free( pairs );
}

#endif
