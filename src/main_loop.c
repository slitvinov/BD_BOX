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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "main_loop.h"

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "err.h"
#include "input.h"
#include "math_help.h"
#include "output.h"
#include "trans.h"

#include "diff_alg/diff_tensor.h"
#include "potentials/calc_func.h"

#define STARS 25

MAKE_STR_IN(INT,step,0,"inner variable")
MAKE_STR_IN(DOUBLE,curr_time,0,"inner variable")

void single_point_print();

/**
 * Serving main loop.
 * Computations and saving of intermediate results.
 */
void main_loop( )
{
    int delta_step = bdsteps / STARS;
    int pr = delta_step;
    int printed=0;
    int i;

    for ( i = 0; i < DIMS1 * size; i += DIMS1 )
    {
        trans_to_box( coord + i );
    }
    if( single_point )
    {
        single_point_print();
        return;
    }

    IFP printf( "%s", STARS ? "[         Start        ]\n" : "" );
    
    IFP
    {
        if( !step ) /* only first time */
        {
            save_pqr( );
        }
        else
        {
            for ( i = 0; i < step; i += delta_step ? delta_step : 1 )
            {
                if( printed + 1 < STARS ) printf( "*" );
                pr += delta_step;
                printed++;
            }
        }
    }
    IFNP enr_file = enr_filename != NULL ? stderr : NULL;
    if ( compute_forces( ( enr_file && !step ) ? E : NULL, 0.0 ) )
    {
        UNERR("Wrong starting position");
    }
    if( !step )
    {
        IFP save_output();
        step = 1;
    }
    for ( ; step < bdsteps; ++step )
    {
        int comp_E = enr_file && ((step % save_enr_freq) == 0);
        IFP
        if ( step >= pr )
        {
            if( printed + 1 < STARS ) printf( "*" );
            fflush( stdout );
            pr += delta_step;
            printed++;
        }

        save_curr_coord();
        save_curr_forces();
        compute_hydro( comp_E ? E : NULL );

        IFP save_output();
    }
    IFP
    {
        while( printed + 1 < STARS )
        {
            printf( "*" );
            printed++;
        }
        printf( "\n" );
    }
}

void single_point_print()
{
    extern DOUBLE* Ds;
    INT i, j;
    FILE* f;
    if ( Ds == NULL )
    {
        Ds = (DOUBLE*) malloc(size*size*DIMS0*DIMS0*sizeof(DOUBLE)); CHMEM(Ds);
    }
    compute_single_D( Ds );
    f = fopen( single_point, "w" );
    if ( hydro )
    {
        for ( i = 0; i < DIMS0*size; ++i )
        {
            for ( j = 0; j < DIMS0*size; ++j )
            {
                if ( i <= j )
                {
                    fprintf( f, "%" FORMAT_DOUBLEG " ", Ds[ j * DIMS0*size + i ] );
                }
                else
                {
                    fprintf( f, "%" FORMAT_DOUBLEG " ", Ds[ i * DIMS0*size + j ] );
                }
            }
            fprintf( f, "\n" );
        }
    }
    else
    {
        for ( i = 0; i < DIMS0*size; ++i )
        {
            for ( j = 0; j < DIMS0*size; ++j )
            {
                if ( i == j )
                {
                    fprintf( f, "%" FORMAT_DOUBLEG " ", diag_D[i] );
                }
                else
                {
                    fprintf( f, "0 " );
                }
            }
            fprintf( f, "\n" );
        }
    }
    fclose(f);
}
