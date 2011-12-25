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

#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

extern "C"
{
#include "../diff_tensor.h"
#include "../../input.h"
#include "../../main_loop.h"
#include "../../math_help.h"
#include "../../myblas.h"
#include "../../rand_move.h"
#include "../../trans.h"

#include "../../potentials/calc_func.h"
#include "../../potentials/LJ.h"
#include "../../err.h"
#include "../../output.h"
#include "../newton.h"
#include "../../cuda.h"
#include "../../restart.h"
}

#if USE_MPI
#include "mpi.h"
#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

/*double g_start = -1;
DOUBLE g_start_c = -1;*/

void mstart( )
{
    /*g_start_c = 1.*clock()/CLOCKS_PER_SEC;
    g_start = omp_get_wtime();*/
}

double mstop( )
{
    /*printf( "%g\n", 1.*clock()/CLOCKS_PER_SEC - g_start_c );
    g_start = omp_get_wtime() - g_start;
    printf( "%g\n", g_start );
    return g_start;*/
    return 0.0;
}

DOUBLE comp_step_IGT_const( DOUBLE* _E )
{
    INT count = 0;
    curr_dt = dt;
    curr_iter++;
    while ( 1 )
    {
        INT res = compute_hydro_IGT( _E );
        if ( res )
        {
            if ( res == 2 )/* if first bas been accepted, revoke only second */
            {
                load_coord();
                while ( compute_hydro_IGT_second_phase( _E ) )
                {
                    load_coord();
                    ++count;
                    if ( count > move_attempts )
                        UNERR("Wrong conformation, can't move")
                }
                break;
            }
            else
            {
                /* revoke predictor step, only forces */
                load_forces();
                curr_iter++;
            }
        }
        else
        {
            break;
        }
        ++count;
        if ( count > move_attempts ) UNERR("Wrong conformation, can't move")
    }
    if ( !d_bforces && ( ( xyz_file && (step % save_xyz_freq ) == 0 ) ||
                         ( dcd_file && (step % save_dcd_freq ) == 0 ) ||
                         ( enr_file && (step % save_enr_freq ) == 0 ) ||
                         (             (step % save_rst_freq ) == 0 ) ) )
    {
        cudaMemcpy( coord, d_coord, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
    }
    return curr_dt;
}

DOUBLE comp_step_IGT_var( DOUBLE* _E )
{
    INT count = 0;
    curr_dt = dt;
    curr_iter++;
    while ( 1 )
    {
        INT res = compute_hydro_IGT( _E );
        if ( res )
        {
            if ( count < 16 )
            {
                curr_dt /= 2;
            }
            load_coord();
            load_forces();
            if ( res == 2 )/* revoke both steps */
            {
                save_iter--;
                curr_iter--;
            }
        }
        else
        {
            break;
        }
        ++count;
        if ( count > move_attempts ) UNERR( "Wrong conformation, can't move")
    }
    if ( !d_bforces && ( ( xyz_file && (step % save_xyz_freq ) == 0 ) ||
                         ( dcd_file && (step % save_dcd_freq ) == 0 ) ||
                         ( enr_file && (step % save_enr_freq ) == 0 ) ||
                         (             (step % save_rst_freq ) == 0 ) ) )
    {
        cudaMemcpy( coord, d_coord, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
    }
    return curr_dt;
}

DOUBLE comp_step_rev_IGT( DOUBLE* _E )
{
    UNERR( "Unimplemented method: comp_step_rev_IGT")
    return -1;
}

DOUBLE comp_step_ermak( DOUBLE* _E )
{
    INT count = 0;
    curr_dt = dt;
    curr_iter++;
    while ( 1 )
    {
        compute_hydro_ermak( d_coord );
        if ( compute_forces( _E, curr_time + curr_dt ) )
        {
            if ( algorithm == DIFF_ALG_ERMAK_VAR && count < 16 )
            {
                curr_dt /= 2;
                --save_iter;
            }
            load_coord();
            load_forces();
        }
        else
        {
            break;
        }
        ++count;
        if ( count > move_attempts ) UNERR( "Wrong conformation, can't move")
    }
    if ( !d_bforces && ( ( xyz_file && (step % save_xyz_freq ) == 0 ) ||
                         ( dcd_file && (step % save_dcd_freq ) == 0 ) ||
                         ( enr_file && (step % save_enr_freq ) == 0 ) ||
                         (             (step % save_rst_freq ) == 0 ) ) )
    {
        cudaMemcpy( coord, d_coord, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
    }
    return curr_dt;
}

DOUBLE comp_step_ermak_newton( DOUBLE* _E )
{
    INT count = 0;
    INT ret = 1;

    curr_iter++;
    while(ret)
    {
        curr_dt = dt;
        compute_hydro_ermak( d_coord );
        ret = compute_forces( _E, curr_time + curr_dt );
        if ( ret )
        {
            ret = 0;/* clear ret */
            newton_init_vel();
            load_coord();
            ret = newton_moves( &count );
            ret = ret || compute_forces( _E, curr_time + curr_dt );
            if( ret )
            {
                load_coord();
                load_forces();
                if( bc != BC_SPHERE ) warning("Revert move in none sphere case",__FILE__,__LINE__);
                if( bc != BC_PWELL ) warning("Revert move in none pwell case",__FILE__,__LINE__);
            }
        }
        else
        {
            break;
        }
    }
    if ( !d_bforces && ( ( xyz_file && (step % save_xyz_freq ) == 0 ) ||
                         ( dcd_file && (step % save_dcd_freq ) == 0 ) ||
                         ( enr_file && (step % save_enr_freq ) == 0 ) ||
                         (             (step % save_rst_freq ) == 0 ) ) )
    {
        cudaMemcpy( coord, d_coord, sizeof(float)*size*DIMS1, cudaMemcpyDeviceToHost ); CCERR
    }
    return dt;
}

extern "C"
void compute_hydro( DOUBLE* _E )
{
    switch ( algorithm )
    {
        case DIFF_ALG_ERMAK_CONST:
        case DIFF_ALG_ERMAK_VAR:
            curr_time += comp_step_ermak( _E );
            break;
        case DIFF_ALG_ERMAK_NEWTON:
            curr_time += comp_step_ermak_newton( _E );
            break;
        case DIFF_ALG_IGT_CONST:
            curr_time += comp_step_IGT_const( _E );
            break;
        case DIFF_ALG_IGT_VAR:
            curr_time += comp_step_IGT_var( _E );
            break;
        case DIFF_ALG_IGT_VAR_REV:
            UNERR( "Unimplemented algorithm")
            curr_time += comp_step_rev_IGT( _E );
            break;
        default: UNERR( "Unimplemented algorithm")
    }
}
