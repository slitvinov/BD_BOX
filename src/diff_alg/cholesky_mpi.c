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

#include "cholesky_mpi.h"

#include <stdio.h>
#include <memory.h>
#include <float.h>

#include "diff_tensor.h"
#include "../err.h"
#include "../math_help.h"

MAKE_STR_IN(INT,MPI_nprow,0,"Number of rows in the processors grid")
MAKE_STR_IN(INT,MPI_npcol,0,"Number of cols in the processors grid")
MAKE_STR_IN(INT,MPI_block,0,"Size of a block in the processors grid")

/*Cblas and ScaLAPACK functions*/
extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, const char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);

extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);

extern void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
            double *b, int *ib, int *jb, int *descb);
extern void pdpotrf_( const char* UPLO, int *n, double *A, int *ia, int *ja, int *desca, int *info);
extern int  indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int  indxl2g_( int *indxloc, int *nb, int *iproc, int *isrcproc, int *nprocs);

#if USE_SCALAPACK && USE_MPI
static int iam, nprocs, ictxt;
static int myrow, mycol, np, nq;
static int izero = 0;
static int ione = 1;
static int descA[9];
#endif

extern INT curr_iter;/* from diff_tensor.c*/
extern INT save_iter;
extern DOUBLE* S[2];
extern DOUBLE* X;
extern DOUBLE** Rs;
extern DOUBLE** F;
extern DOUBLE preDii;
extern DOUBLE preDij;
extern DOUBLE T;
extern DOUBLE curr_dt;

void init_cholesky_mpi()
{
#if USE_SCALAPACK && USE_MPI
    int n, itemp, info;
    if( MPI_nprow * MPI_npcol > g_numprocs )
    {
        UNERR("Not enough processes available to make a process grid")
    }
    Cblacs_pinfo( &iam, &nprocs ); /*get info about id of current proces, number of process*/
    Cblacs_get( -1, 0, &ictxt );/*get current contex*/
    Cblacs_gridinit( &ictxt, "Row", MPI_nprow, MPI_npcol ); /*init current contex*/
    Cblacs_gridinfo( ictxt, &MPI_nprow, &MPI_npcol, &myrow, &mycol );/*get info about contex*/
    n = size * DIMS0;
    np    = numroc_( &n, &MPI_block, &myrow, &izero, &MPI_nprow );
    nq    = numroc_( &n, &MPI_block, &mycol, &izero, &MPI_npcol );
    itemp = np > 1 ? np : 1;/* LD */
    descinit_( descA, &n, &n, &MPI_block, &MPI_block, &izero, &izero, &ictxt, &itemp, &info );
    if( info < 0 ) UNERR( "Error in descinit")
#else
    UNERR("Unimplemented init_cholesky_mpi, use --with-scalapack and --enable-mpi")
#endif
}

void free_cholesky_mpi()
{
#if USE_SCALAPACK && USE_MPI
    Cblacs_gridexit( ictxt );
#else
    UNERR( "Unimplemented free_cholesky_mpi")
#endif
}

void compute_r_cholesky_mpi( DOUBLE* _D )
{
#if USE_SCALAPACK && USE_MPI
    int n = size*DIMS0;
    int info = 0;
    int k, i, j;
    DOUBLE* cD = NULL;
    if ( !curr_iter || save_iter != curr_iter )
    {
        pdpotrf_( "Uplo", &n, _D, &ione, &ione, descA, &info );
        if ( info ) UNERR("Cholesky decomposition error")
        save_iter = curr_iter;
        memcpy( S[ save_iter & 1 ], _D, sizeof (DOUBLE) * np * nq );
    }
    cD = S[ save_iter & 1 ];
    k = 0;
    for (j = 1; j <= nq; j++)
    {
        for (i = 1; i <= np; i++)
        {
            int ic = indxl2g_( &i, &MPI_block, &myrow, &izero, &MPI_nprow ) - 1;
            int jc = indxl2g_( &j, &MPI_block, &mycol, &izero, &MPI_npcol ) - 1;
            if( ic <= jc )
            {
                Rs[0][jc] += cD[k]*X[ic];
            }
            k++;
        }
    }
    MPI_Reduce( Rs[0], Rs[1], size*DIMS0, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    if ( g_id == 0 )
    {
        /* copy 'recv buffor' to 'send buffor' */
        memcpy( Rs[0], Rs[1], sizeof (DOUBLE)*(size * DIMS0) );
    }
    MPI_Bcast( Rs[0], size*DIMS0, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
#else
	_D = NULL;
#endif
}

void compute_Rf_cholesky_mpi( DOUBLE** _Rf, DOUBLE* _D )
{
#if USE_SCALAPACK && USE_MPI
    int k, i, j;
    DOUBLE pre = curr_dt / (kB * T);
    k=0;
    for (j = 1; j <= nq; j++)
    {
        for (i = 1; i <= np; i++)
        {
            int ic = indxl2g_( &i, &MPI_block, &myrow, &izero, &MPI_nprow ) - 1;
            int jc = indxl2g_( &j, &MPI_block, &mycol, &izero, &MPI_npcol ) - 1;
            if( ic <= jc )
            {
                DOUBLE tD[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                if( !ewald )
                {
                    if( (ic/3) == (jc/3) )
                    {
                        _D[k] = ( (ic%3) == (jc%3) ) ? diag_D[ic/3] : 0.0;
                    }
                    else
                    {
                        diff_tensor_calc_pos( ic/3, jc/3, tD );
                        _D[k] = tD[ (ic%3) + 3*(jc%3) ];
                    }
                }
                else if( ewald_method == EWALD_METHOD_SMITH )
                {
                    DOUBLE val = 0.0;
                    if( (ic/3) == (jc/3) )
                    {
                        if( (ic%3) == (jc%3) )
                        {
                            diff_tensor_calc_ewald_smith_diag( ic/3, tD );
                            val = tD[ (ic%3)*4 ];
                        }
                    }
                    else
                    {
                        diff_tensor_calc_ewald_smith( ic/3, jc/3, tD );
                        val = tD[ (ic%3) + 3*(jc%3) ];
                    }
                    _D[k] = val;
                }
            }
            k++;
        }
    }
    memset( _Rf[0], 0, sizeof (DOUBLE) * DIMS0 * size );
    k = 0;
    for (j = 1; j <= nq; j++)
    {
        for (i = 1; i <= np; i++)
        {
            int ic = indxl2g_( &i, &MPI_block, &myrow, &izero, &MPI_nprow ) - 1;
            int jc = indxl2g_( &j, &MPI_block, &mycol, &izero, &MPI_npcol ) - 1;
            if( ic <= jc )
            {
                if( (ic/3) != (jc/3) )
                {
                    _Rf[0][ic] += _D[k]*F[0][jc];
                }
                _Rf[0][jc] += _D[k]*F[0][ic];
            }
            k++;
        }
    }
    for ( i = 0; i < size * DIMS0; ++i )
    {
        _Rf[0][ i ] *= pre;
    }
    MPI_Reduce( _Rf[0], _Rf[1], size*DIMS0, MPI_VALUE_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    if ( g_id == 0 )
    {
        memcpy( _Rf[0], _Rf[1], sizeof (DOUBLE)*(size * DIMS0) );
    }
    MPI_Bcast( _Rf[0], size*DIMS0, MPI_VALUE_TYPE, 0, MPI_COMM_WORLD );
#else
	_Rf = NULL;
	_D = NULL;
#endif
}

void sum_tensors_cholesky_mpi()
{
#if USE_SCALAPACK && USE_MPI
    int k, i, j;
    k = 0;
    for( i = 0; i < np; ++i )
    {
        for( j = 0; j < nq; ++j )
        {
            D[k] = (D[k]+Ds[k])/2;
            k++;
        }
    }
#endif
}
