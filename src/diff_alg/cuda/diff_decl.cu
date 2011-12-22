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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

extern "C"
{
#include "../diff_tensor.h"
#include "../newton.h"
#include "../../cuda.h"
#include "../../err.h"
#include "../../input.h"
#include "../../math_help.h"
#include "../../myblas.h"
#include "../../potentials/calc_func.h"
}

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif


DOUBLE* diag_D = NULL;
DOUBLE* inv_diag_D = NULL;

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

MAKE_STR_IN( DOUBLE, visc, 0.01f, "viscovity[P]" )
MAKE_STR_IN( DOUBLE, vfactor, 14.4f, "factor for viscovity units conversion" )
MAKE_STR_IN( DOUBLE, curr_dt, 0, "inner value" )
MAKE_STR_IN( NO_CHOLS_GEYER, hydro, 0, "method used to evaluate hydrodynamic interactions" )
MAKE_STR_IN( ERMAK_IGTCONST_IGTVAR, algorithm, 2, "integration algorithm" )
MAKE_STR_IN( INT, move_attempts, 10000, "number of attempts to draw a new random vector in case of overlap" )
MAKE_STR_IN( PTR_DOUBLE9, vel_grad_tensor, NULL, "the velocity gradient tensor matrix (row,column)" )
MAKE_STR_IN( DOUBLE, e_collision, 1.0, "coefficient applied in collision for the velocity" )

INT ewald = 0;
MAKE_STR_IN( INT, ewald_real, 10, "an integer which deffines the range of the real-space sum" )
MAKE_STR_IN( INT, ewald_recip, 10, "an integer deffining the summation range in the reciprocal-space and its number of vectors" )
MAKE_STR_IN( DOUBLE, ewald_alpha, 1.77245385090551602729f, "alfa" )
MAKE_STR_IN( EWALD_METHOD, ewald_method, 0, "Ewald method" )

MAKE_STR_IN( YESNO, geyer_on_the_fly, 1, "Geyer on the fly version on/off" )

extern "C"
void init_diff_tensor( )
{
    int i = 0;
    DOUBLE pre;

    diag_D = (DOUBLE*) malloc( sizeof(DOUBLE)*size*DIMS0 ); CHMEM( diag_D);
    inv_diag_D = (DOUBLE*) malloc( sizeof(DOUBLE)*size*DIMS0 );CHMEM(diag_D);
    pre = kB_1_6PI * T / (visc * vfactor);
    for ( i = 0; i < size; ++i )
    {
        diag_D[DIMS0 * i + 2] = diag_D[DIMS0 * i + 1] = diag_D[DIMS0 * i] = pre / coord[ DIMS1 * i + 3 ];
        inv_diag_D[DIMS0 * i + 2] = inv_diag_D[DIMS0 * i + 1] = inv_diag_D[DIMS0 * i] = 1 / diag_D[DIMS0 * i];
    }

    if( !hydro )
    {
        cudaMemcpy( d_D, diag_D, sizeof(DOUBLE)*size*DIMS0, cudaMemcpyHostToDevice ); CCERR
    }

    preDij = kB_1_8PI * T / (visc * vfactor);
    preDii = pre;
    X = (DOUBLE*) malloc( sizeof (DOUBLE)*DIMS0*size );
    CHMEM(X);

    Rf = (DOUBLE**) malloc( sizeof(DOUBLE*)*g_threads ); CHMEM(Rf);

    for ( i = 0; i < g_threads; ++i )
    {
        Rf[i] = (DOUBLE*) malloc( sizeof(DOUBLE)*(DIMS1*size + 1) ); CHMEM(Rf[i]);
    }

    Rs = (DOUBLE**) malloc( sizeof(DOUBLE*)*g_threads );CHMEM(Rs);

    for ( i = 0; i < g_threads; ++i )
    {
        Rs[i] = (DOUBLE*) malloc( sizeof(DOUBLE)*(DIMS1*size + 1) ); CHMEM(Rs[i]);
    }

    init_cholesky( );

    init_tensor();

    epsilon = (DOUBLE*) malloc( sizeof(DOUBLE)*g_threads ); CHMEM(epsilon);
    if ( algorithm == DIFF_ALG_ERMAK_NEWTON )
    {
        init_newton();
    }
    init_cuda_equation();
    init_cuda_rand();
}

void free_diff_tensor( )
{
    INT i;
    free( epsilon );
    if ( hydro == DIFF_CHOLESKY )
    {
        free_cholesky( );
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
    free( X );
    free( diag_D );
    free( inv_diag_D );
    if ( algorithm == DIFF_ALG_ERMAK_NEWTON )
    {
        free_newton();
    }
    free_cuda_rand();
}

extern "C"
void save_curr_coord()
{
    cudaMemcpy( d_scoord, d_coord, sizeof(DOUBLE) * DIMS1 * size, cudaMemcpyDeviceToDevice );
    CCERR
}

extern "C"
void save_curr_forces()
{
    if( d_sforces )
    {
        cudaMemcpy( d_sforces, d_forces, sizeof(DOUBLE) * DIMS1 * size, cudaMemcpyDeviceToDevice );
        CCERR
    }
}

extern "C"
void load_coord()
{
    cudaMemcpy( d_coord, d_scoord, sizeof(DOUBLE) * DIMS1 * size, cudaMemcpyDeviceToDevice );
    CCERR
}

extern "C"
void load_forces()
{
    if( d_sforces )
    {
        cudaMemcpy( d_forces, d_sforces, sizeof(DOUBLE) * DIMS1 * size, cudaMemcpyDeviceToDevice );
        CCERR
    }
}
