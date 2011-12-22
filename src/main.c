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

#include "cuda.h"
#include "data.h"
#include "err.h"
#include "input.h"
#include "main_loop.h"
#include "myblas.h"
#include "output.h"
#include "rand_move.h"
#include "restart.h"

#include "diff_alg/diff_tensor.h"
#include "potentials/bond.h"
#include "potentials/calc_func.h"

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

int main( int argc, char* argv[] )
{
    /*Init phase */
#if USE_MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &g_numprocs );/* get number of process */
    MPI_Comm_rank( MPI_COMM_WORLD, &g_id );/* get id == iam*/
    IFP printf( "MPI version, numprocs %d\n", g_numprocs );
#else
#ifdef _OPENMP
    g_threads = omp_get_max_threads();
    printf( "OpenMP version, threads %d\n", g_threads );
#endif
#endif
#ifdef __GNUC__
    IFP printf( "gcc %d.%d\n", __GNUC__, __GNUC_MINOR__ );
#endif
#ifdef __INTEL_COMPILER
    IFP printf( "icc %d\n", __INTEL_COMPILER );
#endif
#ifdef __CYGWIN__
    IFP printf( "CYGWIN\n" );
#endif
#ifdef _MSC_VER
    IFP printf( "cl %d\n", _MSC_VER );
#endif
#if USE_CUDA || USE_MAGMA
    IFP printf( "CUDA version, devices %d\n", g_devices );
#endif
#if USE_SSE
    IFP printf( "SSE used\n");
#endif
    IFP printf( "%s precision, bytes %lu\n",
		(sizeof(DOUBLE)==4) ? "Single" :
                ((sizeof(DOUBLE)==8) ? "Double" : "Unknow" ), sizeof(DOUBLE) );
    /* Options */
    if ( parse_args( argc, argv ) )
    {
        if( !restart )
        {
            read_prm_file();
            apply_options(0);
            init_rand( rand_seed );
            read_str_file();
            IFP init_save();
        }
        else
        {
            restart_read( restart );
            alloc_connections();
            if( prm_filename )
            {
                read_prm_file();
            }
            apply_options(1);/* values from prm or command line have greater priority */
            IFP init_save_after_restart();
        }
#if USE_CUDA || USE_MAGMA
        if( !restart )
        {
            init_cuda();
        }
        else
        {
            init_cuda_after_restart();
        }
#endif
        init_diff_tensor();
        init_forces();

        main_loop();

        /*Clean up phase*/
        free_forces();
        IFP free_save();
        free_diff_tensor();
        free_data();
        free_input();
        free_restart();
        free_logfile();

#if USE_CUDA || USE_MAGMA
        free_cuda();
#endif
    }
#if USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
