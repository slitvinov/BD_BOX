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

#include <stdlib.h>
#include <memory.h>

#include "../err.h"
#include "../input.h"
#include "../math_help.h"
#include "../myblas.h"
#include "../potentials/calc_func.h"
#include "newton.h"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

DOUBLE* D = NULL;
DOUBLE* diag_D = NULL;
DOUBLE* inv_diag_D = NULL;

DOUBLE* S[2];
INT save_iter = -1; /*time stamp for S*/
INT curr_iter = 0;

DOUBLE preDij;
DOUBLE preDii;
DOUBLE* X = NULL;
DOUBLE** Rf = NULL;
DOUBLE** Rs = NULL;
DOUBLE* Ds = NULL;
DOUBLE* epsilon = NULL;
DOUBLE* FS = NULL;
DOUBLE* Rf0 = NULL;

MAKE_STR_IN( DOUBLE, visc, 0.01, "viscovity" )
MAKE_STR_IN( DOUBLE, vfactor, 14.4, "factor for viscovity units conversion" )
MAKE_STR_IN( DOUBLE, curr_dt, 0, "inner value" )
MAKE_STR_IN( NO_CHOLS_GEYER, hydro, 0, "method used to evaluate hydrodynamic interactions" )
MAKE_STR_IN( ERMAK_IGTCONST_IGTVAR, algorithm, 2, "integration algorithm" )
MAKE_STR_IN( INT, move_attempts, 10000, "number of attempts to draw a new random vector in case of an overlap" )
MAKE_STR_IN( PTR_DOUBLE9, vel_grad_tensor, NULL, "the velocity gradient tensor matrix (row,column)" )
MAKE_STR_IN( DOUBLE, e_collision, 1, "restitution parameter used in the elastic collision method" )

INT ewald = 0;
MAKE_STR_IN( INT, ewald_real, 10, "magnitude of the real lattice vectors" )
MAKE_STR_IN( INT, ewald_recip, 10, "magnitude of the reciprocal lattice vectors" )
MAKE_STR_IN( DOUBLE, ewald_alpha, 1.77245385090551602729, "Ewald summation convergence control" )
MAKE_STR_IN( EWALD_METHOD, ewald_method, 0, "method used to built the Ewald-summed tensor" )

MAKE_STR_IN( YESNO, geyer_on_the_fly, 1, "whether to use 3x3 submatrices of D (yes) or the whole D during evaluation of HI (no) with TEA-HI" )

void init_diff_tensor( )
{
    int i = 0;
    DOUBLE pre;
    if ( hydro && algorithm > DIFF_ALG_ERMAK_NEWTON )
    {
        D = (DOUBLE*) malloc( sizeof(DOUBLE)*size*size*DIMS0*DIMS0 ); CHMEM(D);
        memset( D, 0, sizeof(DOUBLE)*DIMS0*DIMS0*size*size );
    }
    diag_D = (DOUBLE*) malloc( sizeof(DOUBLE)*size*DIMS0 ); CHMEM(diag_D);
    inv_diag_D = (DOUBLE*) malloc( sizeof(DOUBLE)*size*DIMS0 ); CHMEM(diag_D);
    pre = kB_1_6PI * T / (visc * vfactor);
    for ( i = 0; i < size; ++i )
    {
        diag_D[DIMS0 * i + 2] = diag_D[DIMS0 * i + 1] = diag_D[DIMS0 * i] = pre / coord[ DIMS1 * i + 3 ];
        inv_diag_D[DIMS0 * i + 2] = inv_diag_D[DIMS0 * i + 1] = inv_diag_D[DIMS0 * i] = 1 / diag_D[DIMS0 * i];
    }
    preDij = kB_1_8PI * T / (visc * vfactor);
    preDii = pre;
    X = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS0*size ); CHMEM(X);

#if USE_MPI
    Rf = (DOUBLE**) malloc( sizeof (DOUBLE*) * (g_threads + 1) ); CHMEM(Rf);
#else
    Rf = (DOUBLE**) malloc( sizeof (DOUBLE*) * g_threads ); CHMEM(Rf);
#endif

    for ( i = 0; i < g_threads; ++i )
    {
        Rf[i] = (DOUBLE*) malloc( sizeof (DOUBLE) * (DIMS1 * size + 1) ); CHMEM(Rf[i]);
    }

#if USE_MPI
    if ( g_id == 0 )
    {
        Rf[ g_threads ] = (DOUBLE*) malloc( sizeof (DOUBLE) * (DIMS1 * size + 1) );
        CHMEM(Rf[i]);
    }
    else
    {
        Rf[ g_threads ] = NULL;
    }
#endif

#if USE_MPI
    Rs = (DOUBLE**) malloc( sizeof (DOUBLE*) * (g_threads + 1) ); CHMEM(Rs);
#else
    Rs = (DOUBLE**) malloc( sizeof (DOUBLE*) * g_threads ); CHMEM(Rs);
#endif

    for ( i = 0; i < g_threads; ++i )
    {
        Rs[i] = (DOUBLE*) malloc( sizeof (DOUBLE) * (DIMS1 * size + 1) ); CHMEM(Rs[i]);
    }

#if USE_MPI
    if ( g_id == 0 )
    {
        Rs[g_threads] = (DOUBLE*) malloc( sizeof (DOUBLE) * (DIMS1 * size + 1) );
        CHMEM(Rs[i]);
    }
    else
    {
        Rs[g_threads] = NULL;
    }
#endif

    Rf0 = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS1*size ); CHMEM(Rf0);
    if ( vel_grad_tensor )
    {
        FS = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS1*size ); CHMEM(FS);
    }
    if ( hydro && !geyer_on_the_fly )
    {
        Ds = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS0*DIMS0*size*size ); CHMEM(Ds);
        memset( Ds, 0, sizeof (DOUBLE)*DIMS0*DIMS0*size*size );
    }
    if ( hydro == DIFF_CHOLESKY )
    {
        S[0] = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS0*DIMS0*size*size ); CHMEM(S[0]);
        S[1] = (DOUBLE*) malloc( sizeof (DOUBLE)*DIMS0*DIMS0*size*size ); CHMEM(S[1]);
        init_cholesky( );
#if USE_MPI
        init_cholesky_mpi( );
#endif
    }
    epsilon = (DOUBLE*) malloc( sizeof (DOUBLE) * g_threads ); CHMEM(epsilon);
    if ( algorithm == DIFF_ALG_ERMAK_NEWTON )
    {
        init_newton();
    }
    init_tensor();
}

void free_diff_tensor( )
{
    INT i;
    free( epsilon );
    if ( hydro == DIFF_CHOLESKY )
    {
        free_cholesky( );
        free( S[0] );
        free( S[1] );
#if USE_MPI
        free_cholesky_mpi( );
#endif
    }
    if ( Rf )
    {
        for ( i = 0; i < g_threads; ++i )
        {
            free( Rf[i] );
        }
        free( Rf );
    }
    free( Rf0 );
    if ( FS ) free( FS );
    if ( Rs )
    {
        for ( i = 0; i < g_threads; ++i )
        {
            free( Rs[i] );
        }
        free( Rs );
    }
    if ( Ds ) free( Ds );
    free( X );
    free( diag_D );
    free( inv_diag_D );
    free( D );
    if ( algorithm == DIFF_ALG_ERMAK_NEWTON )
    {
        free_newton();
    }
}

void save_curr_coord()
{
    memcpy( save_coord, coord, sizeof (DOUBLE) * DIMS1 * size );
}

void save_curr_forces()
{
    memcpy( save_forces, F[0], sizeof (DOUBLE) * DIMS0 * size );
}

void load_coord()
{
    memcpy( coord, save_coord, sizeof (DOUBLE) * DIMS1 * size );
}

void load_forces()
{
    memcpy( F[0], save_forces, sizeof (DOUBLE) * DIMS0 * size );
}
