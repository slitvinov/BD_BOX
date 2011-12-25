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

#include "calc_func.h"

#include <memory.h>
#include <stdlib.h>

#include "../err.h"
#include "../trans.h"

#include "angle.h"
#include "angle_cos.h"
#include "bond.h"
#include "bucket.h"
#include "dihe.h"
#include "dihe_angle.h"
#include "electro.h"
#include "electro_ext.h"
#include "LJ.h"
#include "sboundary.h"
#include "pwell.h"

#if _OPENMP
#include <omp.h>
#endif

DOUBLE** ES;
DOUBLE** F; /*Forces*/
INT* repeat;

INT init_iter_forces( DOUBLE* _E )
{
    INT i;
    if ( !iscorrect( ) ) return 1;
    memset( repeat, 0, sizeof (INT) * g_threads );
    for ( i = 0; i < g_threads; ++i )
    {
        memset( F[i], 0, sizeof (DOUBLE) * DIMS0 * size );
    }
    if ( _E )
    {
        memset( _E, 0, sizeof (DOUBLE) * (ENERGY_LAST + 1) );
        for ( i = 0; i < g_threads; ++i )
        {
            memset( ES[i], 0, sizeof (DOUBLE) * (ENERGY_LAST + 1) );
        }
    }
	return 0;
}

INT collect_results( DOUBLE* _E )
{
    INT i;
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
    return repeat[0];
}

INT bonded_inter( DOUBLE* _E, INT tid, DOUBLE curr_time )
{
    INT i, rep=0;
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
    for ( i = g_id; i < angle_cos_size; i += g_numprocs )
    {
        rep |= angle_cos_calc( F[tid], i, _E ? ES[tid] + ENERGY_ANGLE_COS : _E );
    }
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
    for ( i = g_id; i < angle_size; i += g_numprocs )
    {
        rep |= angle_calc( F[tid], i, _E ? ES[tid] + ENERGY_ANGLE : _E );
    }
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
    for ( i = g_id; i < bond_size; i += g_numprocs )
    {
        rep |= bond_calc( F[tid], i, _E ? ES[tid] + ENERGY_BOND : _E );
    }
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
    for ( i = g_id; i < dihe_size; i += g_numprocs )
    {
        rep |= dihe_calc( F[tid], i, _E ? ES[tid] + ENERGY_DIHE : _E );
    }
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
    for ( i = g_id; i < dihe_angle_size; i += g_numprocs )
    {
        rep |= dihe_calc_angle( F[tid], i, _E ? ES[tid] + ENERGY_DIHE_ANGLE : _E );
    }
    if ( E_ext )
    {
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
        for ( i = g_id; i < size; i += g_numprocs )
        {
            if ( Q[i] != 0.0 )
            {
                rep |= electro_ext_force( Q[i], F[tid] + DIMS0*i, curr_time );
            }
        }
    }
    if( sboundary && bc == BC_SPHERE )
    {
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
        for ( i = g_id; i < size; i += g_numprocs )
        {
            rep |= sboundary_ext( coord + i * DIMS1, F[tid] + DIMS0*i );
        }
    }
    if( pwell && bc == BC_PWELL )
    {
#ifdef _OPENMP
#pragma omp for schedule( static ) nowait
#endif
        for ( i = g_id; i < size; i += g_numprocs )
        {
            rep |= pwell_ext( coord + i * DIMS1, F[tid] + DIMS0*i );
        }
    }
    return rep;
}


INT bonded_inter_parallel( DOUBLE* _E, DOUBLE curr_time )
{
#if _OPENMP
    int i;
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        repeat[tid] = bonded_inter( _E, tid, curr_time );
    }
    if ( _E )
    {
        for ( i = 1; i < g_threads; ++i )
        {
            INT j;
            for ( j = 0; j < ENERGY_LAST; ++j )
            {
                ES[0][j] += ES[i][j];
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
    return repeat[0];
#else
    return bonded_inter( _E, 0, curr_time );
#endif
}

void init_cutoff()
{
    if ( nb_list )
    {
        INT i;
        max_force_cutoff = 0;
        if ( elec ) max_force_cutoff = cutoff_c + 0.01f;
        if ( alpha_lj != 0 )
        {
            if ( cutoff_lj == -1 )
            {
                for ( i = 0; i < size; ++i )
                {
                    if ( 2 * 1.1225f * LJ[2 * i] > max_force_cutoff )
                    {
                        max_force_cutoff = 2 * 1.1225f * LJ[2 * i] + 0.01f;
                    }
                }
            }
            else
            {
                if ( cutoff_lj > max_force_cutoff )
                {
                    max_force_cutoff = cutoff_lj + 0.01f;
                }
            }
        }
        for ( i = 0; i < size; ++i )
        {
            if ( 2 * LJ[2 * i] > max_force_cutoff )
            {
                max_force_cutoff = 2 * LJ[2 * i] + 0.01f;/*because it's radius*/
            }
        }
    }
}

void init_Fs()
{
    INT i;
    repeat = (INT*) malloc( sizeof(INT)*g_threads ); CHMEM(repeat);
    E = (DOUBLE*) malloc( sizeof(DOUBLE)*(ENERGY_LAST + 1) ); CHMEM(E);
#if USE_MPI
    F = (DOUBLE**) malloc( sizeof (DOUBLE*) * (g_threads + 1) ); CHMEM(F);
#else
    F = (DOUBLE**) malloc( sizeof(DOUBLE*)*g_threads ); CHMEM(F);
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
}

void free_Fs()
{
    INT i;
    for ( i = 0; i < g_threads; ++i )
    {
        free( F[i] );
        free( ES[i] );
    }
    free( F );
    free( ES );
    free( repeat );
}
