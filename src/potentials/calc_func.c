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

#include "calc_func.h"

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "../err.h"
#include "../input.h"
#include "../math_help.h"
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

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

DOUBLE* E;

MAKE_STR_IN( DOUBLE, bond_lj_scale, 1.0, "scaling factor for Lennard-Jones interactions between bonded partners 1-2" )
MAKE_STR_IN( DOUBLE, bond_c_scale, 1.0, "scaling factor for electrostatic interactions between bonded partners 1-2" )
MAKE_STR_IN( YESNO, check_overlap, 1, "whether to detect overlaps after each step" )

/*!
 * Function computes nonbonded interactions beetwen two beads
 * \param i first bead id
 * \param j second bead id
 * \param Fi pointer to output force array
 * \param lEc pointer to output variable containnig summary energy of electrostatic potential
 * \param lElj pointer to output variable containnig summary energy of LJ potential
 * \return one if beads overlap, zero otherwise
 */
INT force_single( INT i, INT j, DOUBLE* Fi, DOUBLE* lEc, DOUBLE* lElj );

/*!
 * Function computes interaction beetwen beads in system.
 * Use spatial subdivision and bucket counting for it.
 * \param _E Output energy, can be null
 * \param tid Id of computing thread
 * \param br Additional help structure
 * \return one if beads overlap, zero otherwise
 */
INT bucket_eval( DOUBLE* _E, INT tid, bucket_f* br );

/*!
 * Function computes interaction beetwen beads in system.
 * Use spatial subdivision and sortin for it.
 * \param _E Output energy, can be null
 * \param tid Id of computing thread
 * \param bs Additional help structure
 * \return one if beads overlap, zero otherwise
 */
INT spatial_eval( DOUBLE* _E, INT tid, bucket_s* bs );

/*!
 * Function computes interaction beetwen beads in system.
 * Use naive method, checks all pairs
 * \param _E Output energy, can be null
 * \param tid Id of computing thread
 * \return one if beads overlap, zero otherwise
 */
INT brute_force_val( DOUBLE* _E, INT tid );

INT compute_forces( DOUBLE* _E, DOUBLE curr_time )
{
    INT tid = 0;
    INT rep = 0;
    bucket_f* br = NULL;
    bucket_s* bs = NULL;
    /** Init phase */
    init_iter_forces( _E );
    br = (nb_list == BUCKET_BUCKET) ? (init_bucket_f( size, (DOUBLE4*) coord, NULL, max_force_cutoff )) : NULL;
    bs = (nb_list == BUCKET_SPATIAL) ? (init_keys( size, (DOUBLE4*) coord, NULL, max_force_cutoff )) : NULL;
    /** Computation phase */
#ifdef _OPENMP
#pragma omp parallel private( tid, rep )
    {
#endif
#ifdef _OPENMP
        rep = 0;
        tid = omp_get_thread_num( );
#endif
        /**-- Bonded interactions*/
        rep |= bonded_inter( _E, tid, curr_time );
        /**-- Non-Bonded interactions*/
        if ( br )
        {
            rep |= bucket_eval( _E, tid, br );
        }
        else if ( bs )
        {
            rep |= spatial_eval( _E, tid, bs );
        }
        else
        {
            rep |= brute_force_val( _E, tid );
        }
        repeat[tid] = rep;
#ifdef _OPENMP
    }
#endif

    /** OpenMP - Integrating phase, summarizing forces and energies  */
    rep |= collect_results( _E );

    if ( br ) free_bucket_f( br, size, (DOUBLE4*) coord );
    if ( bs ) free_keys( bs );

    /** MPI - Integrating phase, send results */
#if USE_MPI
    /**-- Reduce forces and energies over nodes*/
    F[0][size * DIMS0] = repeat[0];
    if ( _E )
    {
        int i;
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
    /**-- Broadcast from root node to all*/
    MPI_Bcast( F[0], size * DIMS0 + 1, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
    if ( _E && g_id == 0 ) /*Only master saves results*/
    {
        int i;
        for ( i = 1; i <= ENERGY_LAST; ++i )
        {
            _E[ i - 1 ] = F[0][ size * DIMS0 + i ];
        }
    }
    rep = repeat[0] = F[0][size * DIMS0] != 0.0;
#endif
    return rep;
}

INT brute_force_val( DOUBLE* _E, INT tid )
{
    INT i, rep = 0;
    if ( elec || alpha_lj != 0.0 ) /* if elec or lj turn on*/
    {
        DOUBLE* lEc;
        DOUBLE* lElj;
        DOUBLE* Fi;
        lEc = _E ? ES[ tid ] + ENERGY_COULOMB : _E;
        lElj = _E ? ES[ tid ] + ENERGY_LJ : _E;
        Fi = F[ tid ];
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
        for ( i = 0; i < size; ++i )
        {
            INT j;
            for ( j = i + 1 + g_id; j < size && !rep; j += g_numprocs )
            {
                rep = rep || force_single( i, j, Fi, lEc, lElj );
            }
        }
    }
    else if( check_overlap )
    {
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
        for ( i = 0; i < size; ++i )
        {
            INT j;
            for ( j = i + 1 + g_id; j < size && !rep; j += g_numprocs )
            {
                DOUBLE r_sq;
                get_v2_sq( i, j, &r_sq );
                rep = rep || LJ_force_check( r_sq, LJ + 2 * i, LJ + 2 * j );
            }
        }
    }
    return rep;
}

INT spatial_eval( DOUBLE* _E, INT tid, bucket_s* bs )
{
    INT rep, k;
    rep = 0;
    if ( elec || alpha_lj != 0.0 ) /* if elec or lj turn on*/
    {
        DOUBLE* lEc;
        DOUBLE* lElj;
        DOUBLE* Fi;
        lEc = _E ? ES[ tid ] + ENERGY_COULOMB : _E;
        lElj = _E ? ES[ tid ] + ENERGY_LJ : _E;
        Fi = F[ tid ];
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
        for ( k = g_id; k < bs->size_ids; k += g_numprocs )
        {
            INT b = bs->ids[2 * k];
            INT e = bs->ids[2 * k + 1];
            INT i;
            for ( i = b; i < e && !rep; ++i )
            {
                INT j;
                INT n = bs->keys[i].value;
                if (n>=size)
                {
                    if ( n-size >= bs->size_phan ) UNERR("Assertion fault");
                    n = bs->phantoms[n-size];
                }
                for ( j = i + 1; j < e && !rep; ++j )
                {
                    INT m = bs->keys[j].value;
                    if (m>=size)
                    {
                        if ( m-size >= bs->size_phan ) UNERR("Assertion fault");
                        m = bs->phantoms[m-size];
                    }
                    if ( n != m && can_interact( bs, bs->keys + i, bs->keys + j ) )
                    {
                        rep = rep || force_single( n, m, Fi, lEc, lElj );
                    }
                }
            }
        }
    }
    else if( check_overlap )
    {
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
        for ( k = g_id; k < bs->size_ids; k += g_numprocs )
        {
            INT b = bs->ids[2 * k];
            INT e = bs->ids[2 * k + 1];
            INT i;
            for ( i = b; i < e && !rep; ++i )
            {
                INT j;
                INT n = bs->keys[i].value;
                for ( j = i + 1; j < e && !rep; ++j )
                {
                    DOUBLE r_sq;
                    INT m = bs->keys[j].value;
                    if ( n != m )
                    {
                        get_v2_sq( n, m, &r_sq );
                        rep = rep || LJ_force_check( r_sq, LJ + 2 * n, LJ + 2 * m );
                    }
                }
            }
        }
    }
    return rep;
}

INT bucket_eval( DOUBLE* _E, INT tid, bucket_f* br )
{
    INT i, rep;
    rep = 0;
    if ( elec || alpha_lj != 0.0 ) /* if elec or lj turn on*/
    {
        DOUBLE* lEc;
        DOUBLE* lElj;
        DOUBLE* Fi;
        lEc = _E ? ES[ tid ] + ENERGY_COULOMB : _E;
        lElj = _E ? ES[ tid ] + ENERGY_LJ : _E;
        Fi = F[ tid ];
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
        for ( i = g_id; i < size; i += g_numprocs )
        {
            int sizes[64];
            int* ids[64];
            int k;
            get_neighbours( br, (DOUBLE4*) coord + i, sizes, ids );
            for ( k = 0; k < sizeof (sizes) / sizeof (sizes[0]) && sizes[k]; ++k )
            {
                int l;
                for ( l = 0; l < sizes[k]; ++l )
                {
                    int j = ids[k][l];
                    if ( i < j )
                    {
                        rep = rep || force_single( i, j, Fi, lEc, lElj );
                    }
                }
            }
        }
    }
    else if( check_overlap )
    {
#ifdef _OPENMP
#pragma omp for schedule( dynamic )
#endif
        for ( i = g_id; i < size; i += g_numprocs )
        {
            int sizes[64];
            int* ids[64];
            int k;
            get_neighbours( br, (DOUBLE4*) coord + i, sizes, ids );
            for ( k = 0; k < sizeof (sizes) / sizeof (sizes[0]) && sizes[k]; ++k )
            {
                int l;
                for ( l = 0; l < sizes[k]; ++l )
                {
                    int j = ids[k][l];
                    DOUBLE r_sq;
                    if ( i < j )
                    {
                        get_v2_sq( i, j, &r_sq );
                        rep = rep || LJ_force_check( r_sq, LJ + 2 * i, LJ + 2 * j );
                    }
                }
            }
        }
    }
    return rep;
}

INT force_single( INT i, INT j, DOUBLE* Fi, DOUBLE* lEc, DOUBLE* lElj )
{
    DOUBLE r_norm;
    DOUBLE r_sq;
    DOUBLE e[3];
    int rep = 0;
    get_v2_norms( i, j, e, &r_norm, &r_sq );

    if ( elec && r_norm <= cutoff_c )
    {
        if ( !are_conn( i, j ) )
        {
            electro_force( r_norm, r_sq, e, Fi + i*DIMS0, Fi + j*DIMS0, Q[i], Q[j], lEc );
        }
        else if ( bond_c_scale != 0.0 )
        {
            electro_force_scale( r_norm, r_sq, e, Fi + i*DIMS0, Fi + j*DIMS0, Q[i], Q[j], bond_c_scale, lEc );
        }
    }
    if ( alpha_lj != 0.0 && r_norm <= (cutoff_lj > 0.0 ? cutoff_lj : 1.122462f * (LJ[2 * i + LJ_SIGMA_OFF] + LJ[2 * j + LJ_SIGMA_OFF])) )
    {
        if ( !are_conn( i, j ) )
        {
            rep = LJ_force( r_norm, e, Fi + i*DIMS0, Fi + j*DIMS0, LJ + 2 * i, LJ + 2 * j, lElj );
        }
        else if ( bond_lj_scale != 0.0 )
        {
            rep = LJ_force_scale( r_norm, e, Fi + i*DIMS0, Fi + j*DIMS0, LJ + 2 * i, LJ + 2 * j, bond_lj_scale, lElj );
        }
        else
        {
            rep = LJ_force_check( r_sq, LJ + 2 * i, LJ + 2 * j );
        }
    }
    else
    {
        rep = LJ_force_check( r_sq, LJ + 2 * i, LJ + 2 * j );
    }
    return check_overlap && rep;
}

void init_forces( )
{
    init_cutoff();
    init_Fs( );
}

void free_forces( )
{
    free_Fs( );
    free( E );
}

#endif
