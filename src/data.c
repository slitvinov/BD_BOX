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

#include "data.h"

#include <malloc.h>
#include <string.h>

#include "err.h"
#include "input.h"

#include "potentials/angle.h"
#include "potentials/angle_cos.h"
#include "potentials/bond.h"
#include "potentials/dihe.h"
#include "potentials/dihe_angle.h"

/*Global data*/
MAKE_STR_IN(INT,size,0,"inner variable")

MAKE_STR_IN(PTR_DOUBLE,coord,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,Q,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,LJ,NULL,"inner variable")
MAKE_STR_IN(PTR_STR,names,NULL,"inner variable")
MAKE_STR_IN(PTR_INT,ids,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,masses,NULL,"inner variable")

DOUBLE box[4];
INT g_devices = 1;
INT g_threads = 1;
INT g_id = 0;
INT g_numprocs = 1;
DOUBLE* save_coord;
DOUBLE* save_forces;

void free_data()
{
    int i;
    for( i = 0; i < size; ++i )
    {
        free( names[i] );
    }
    free(ids);
    free(names);
    free(LJ);
    free(Q);
    free(coord);
    free(masses);
    
    free( save_coord );
    free( save_forces );

    free_angle();
    free_angle_cos();
    free_bond();
    free_dihe();
    free_dihe_angle();
}

static INT capacity = 0;
void add_sub(CSTR namei, INT idi, DOUBLE xi, DOUBLE yi, DOUBLE zi, DOUBLE sigmai, DOUBLE Qi, DOUBLE sigma_LJi, DOUBLE epsilon_LJi, DOUBLE mass )
{
    if( capacity == size )
    {
        capacity = capacity ? 2*capacity : 128;
        coord = (DOUBLE*) realloc( coord, sizeof(DOUBLE) * 4 * capacity ); CHMEM(coord);
        Q = (DOUBLE*) realloc( Q, sizeof(DOUBLE) * capacity ); CHMEM(Q);
        LJ = (DOUBLE*) realloc( LJ, sizeof(DOUBLE) * 2 * capacity ); CHMEM(LJ);
        names = (STR*) realloc( names, sizeof(STR*) * capacity );  CHMEM(names);
        ids = (INT*) realloc( ids, sizeof(INT) * capacity );  CHMEM(ids);
        masses = (DOUBLE*) realloc( masses, sizeof(DOUBLE) * capacity ); CHMEM(masses);
    }
    coord[ 4 * size ] = xi;
    coord[ 4 * size + 1 ] = yi;
    coord[ 4 * size + 2 ] = zi;
    coord[ 4 * size + 3 ] = sigmai;
    Q[size] = Qi;
    LJ[2*size + 0 ] = sigma_LJi/2;
    LJ[2*size + 1 ] = epsilon_LJi;
    ids[ size ] = idi;
    masses[ size ] = mass;
    names[ size ] = (STR) malloc( (strlen(namei) + 1) * sizeof(CHAR) );
    strcpy( names[size], namei );
    size++;
}

void end_sub()
{
    if( size )
    {
        int i;
        coord = (DOUBLE*) realloc( coord, sizeof(DOUBLE) * 4 * size );
        MAKE_STR_SIZE(coord,4 * size)
        Q = (DOUBLE*) realloc( Q, sizeof(DOUBLE) * size );
        MAKE_STR_SIZE(Q,size)
        LJ = (DOUBLE*) realloc( LJ, sizeof(DOUBLE) * 2 * size );
        MAKE_STR_SIZE(LJ,2*size)
        names = (STR*) realloc( names, sizeof(STR*) * size );
        MAKE_STR_SIZE(names,size)
        ids = (INT*) realloc( ids, sizeof(INT) * size );
        MAKE_STR_SIZE(ids,size)
        masses = (DOUBLE*) realloc( masses, sizeof(DOUBLE) * capacity );
        MAKE_STR_SIZE(masses,size)
        save_coord = (DOUBLE*) malloc( sizeof (DOUBLE) * DIMS1 * size ); CHMEM(save_coord);
        save_forces = (DOUBLE*) malloc( sizeof (DOUBLE) * DIMS1 * size ); CHMEM(save_forces);

        /* check LJ sigma */
        for ( i = 0; i < size; ++i )
        {
            if( LJ[ 2*i ] == 0.0 )
            {
                warning( "LJ radii is zero", __FILE__, __LINE__ );
                break;
            }
        }
    }
    else
    {
        UNERR("No beads");
    }
}

INT resolve_pos(INT id)
{
    int i;
    for( i = 0; i < size; ++i )
    {
        if( ids[i] == id )
        {
            return i;
        }
    }
    UNERR("resolve_pos unknown id");
    return -1;
}
