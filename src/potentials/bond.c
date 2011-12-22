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

#include "bond.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "calc_func.h"
#include "electro.h"
#include "LJ.h"

#include "../err.h"
#include "../math_help.h"

MAKE_STR_IN(INT,bond_size,0,"inner variable")
MAKE_STR_IN(PTR_INT,bond_pairs,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,bond_parms,NULL,"inner variable")
UCHAR* bond_conns;

INT bond_calc(DOUBLE* F, int pos, DOUBLE* _E )
{
    DOUBLE r_norm;
    DOUBLE r_sq;
    DOUBLE e[3];
    DOUBLE r_0 = bond_parms[ 3*pos ];
    double r_max = bond_parms[ 3*pos + 1 ];
    DOUBLE H = bond_parms[ 3*pos + 2 ];
    DOUBLE Ff = 0;
    DOUBLE dist_rmax;
    DOUBLE* F1 = F + DIMS0 * bond_pairs[2*pos];
    DOUBLE* F2 = F + DIMS0 * bond_pairs[2*pos + 1];

    get_v2_norms( bond_pairs[2*pos], bond_pairs[2*pos+1], e, &r_norm, &r_sq );
    if ( r_norm < r_max )
    {
        if(_E) _E[0] += -0.5 * H * (
            r_max*r_max * log((r_max*r_max-r_sq)/(r_max*r_max-r_0*r_0))
          + r_0 * r_max * log((r_max+r_norm)*(r_max-r_0)/((r_max-r_norm)*(r_max+r_0)) ) );
        dist_rmax = 1 - r_sq/(r_max*r_max);
        Ff = H * ( r_0 - r_norm ) / dist_rmax;

        F1[0] -= Ff * e[0];
        F1[1] -= Ff * e[1];
        F1[2] -= Ff * e[2];

        F2[0] += Ff * e[0];
        F2[1] += Ff * e[1];
        F2[2] += Ff * e[2];
    }
    else
    {
        CHAR buff[1024];
        sprintf( buff, "Overstretched bond between subunits %d and %d", bond_pairs[2*pos], bond_pairs[2*pos+1] );
        UNERR( buff );
    }

    return isnan(Ff) || isinf(Ff) || 
           isnan(e[0]) || isinf(e[0]) ||
           isnan(e[1]) || isinf(e[1]) ||
           isnan(e[2]) || isinf(e[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

int are_conn(INT p1,INT p2)
{
    int a2 = p2/sizeof(UCHAR);
    return bond_conns && bond_conns[ p1 * (size+sizeof(UCHAR))/sizeof(UCHAR) + a2 ] &
           ( 1 << (p2 % sizeof(UCHAR)) );
}

static INT capacity = 0;
void add_bond( INT id1, INT id2, DOUBLE r_0, DOUBLE r_max, DOUBLE H )
{
    if( capacity == bond_size )
    {
        capacity = capacity ? 2*capacity : 128;
        bond_pairs = (INT*) realloc( bond_pairs, sizeof(INT) * 2 * capacity ); 
        CHMEM(bond_pairs);
        bond_parms = (DOUBLE*) realloc( bond_parms, sizeof(DOUBLE) * 3 * capacity );
        CHMEM(bond_parms);
    }
    bond_pairs[2*bond_size] = id1;
    bond_pairs[2*bond_size + 1] = id2;
    bond_parms[3*bond_size] = r_0;
    bond_parms[3*bond_size + 1] = r_max;
    bond_parms[3*bond_size + 2] = H;
    bond_size++;
}

void alloc_connections()
{
    int i;
    if( bond_size && (( elec && bond_c_scale != 1 ) || ( alpha_lj != 0.0 && bond_lj_scale != 1 )) )
    {
        bond_conns = (UCHAR*) malloc( sizeof(UCHAR) * size * ( size + sizeof(UCHAR) ) / sizeof(UCHAR) );
        memset( bond_conns, 0, sizeof(UCHAR) * size * ( size + sizeof(UCHAR) ) / sizeof(UCHAR) );
        for( i = 0; i < bond_size; ++i )
        {
            int p1 = bond_pairs[ 2 * i ];
            int p2 = bond_pairs[ 2 * i + 1 ];
            int a1 = p1/sizeof(UCHAR);
            int a2 = p2/sizeof(UCHAR);
            bond_conns[ p1 * (size+sizeof(UCHAR))/sizeof(UCHAR) + a2 ] |= 1 << (p2 % sizeof(UCHAR));
            bond_conns[ p2 * (size+sizeof(UCHAR))/sizeof(UCHAR) + a1 ] |= 1 << (p1 % sizeof(UCHAR));
        }
    }
}

INT* get_connection_graph( INT* max_group )
{
    INT* counts;
    INT i;
    INT max_g = 0;
    INT* ret;

    *max_group = 0;
    if( bond_size == 0 ) return NULL;
    counts = (INT*) malloc( sizeof(INT) * size ); CHMEM(counts);
    memset( counts, 0, sizeof(INT)*size );

    for( i = 0; i < bond_size; ++i ) /* count connection in graph*/
    {
        counts[ bond_pairs[2*i  ] ]++;
        counts[ bond_pairs[2*i+1] ]++;
    }
    for( i = 0; i < size; ++i )
    {
        if( max_g < counts[i] ) max_g = counts[i];
    }
    ret = (INT*) malloc( sizeof(INT) * size * max_g ); CHMEM(ret);
    memset( ret, -1, sizeof(INT)*size );/* set unused edge to -1 */
    memset( counts, 0, sizeof(INT)*size );
    for( i = 0; i < bond_size; ++i )
    {
        ret[ max_g*bond_pairs[2*i  ] + counts[ bond_pairs[2*i  ] ] ] = bond_pairs[2*i+1];
        ret[ max_g*bond_pairs[2*i+1] + counts[ bond_pairs[2*i+1] ] ] = bond_pairs[2*i  ];
        counts[ bond_pairs[2*i  ] ]++;
        counts[ bond_pairs[2*i+1] ]++;
    }
    *max_group =  max_g;
    return ret;
}

void end_bond()
{
    if( bond_size )
    {
        bond_parms = (DOUBLE*) realloc( bond_parms, sizeof(DOUBLE) * 3 * bond_size ); CHMEM(bond_parms);
        MAKE_STR_SIZE(bond_parms,3*bond_size)
        bond_pairs = (INT*) realloc( bond_pairs, sizeof(INT) * 2 * bond_size ); CHMEM(bond_pairs);
        MAKE_STR_SIZE(bond_pairs,2*bond_size)
    }
    alloc_connections();
}

void free_bond()
{
    if( bond_conns )
    {
        free( bond_conns );
        bond_conns = NULL;
    }
    if( bond_pairs )
    {
        free( bond_pairs );
        bond_pairs = NULL;
    }
    if( bond_parms )
    {
        free( bond_parms );
        bond_parms = NULL;
    }
}
