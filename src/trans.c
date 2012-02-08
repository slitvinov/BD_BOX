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

#include "trans.h"

#include <math.h>
#include <stdio.h>

#include "err.h"
#include "input.h"
#include "math_help.h"

DOUBLE inv_L[4];
MAKE_STR_IN(NONE_CUBIC_SPHERE,bc,1,"periodic boundary conditions, currently none, a rectangular box or a reflective sphere")
MAKE_STR_IN(DOUBLE,sphere_radius,0,"radius of the reflective sphere enclosing the studied system")
#define EPSILON 0.000001f

void center_coords()
{
    INT i, j;
    DOUBLE m[3][2];
    DOUBLE g[3] = {0.0, 0.0, 0.0};
    DOUBLE begin[3] = {0.0, 0.0, 0.0};
    if( size == 0 ) UNERR( "Empty str file, no particles");
    if( !bc ) return;
    if( bc == BC_PWELL ) return;
    for( i = 0; i < 3; ++i )
    {
        m[i][0] = coord[i];
        m[i][1] = coord[i];
        g[i] += coord[i];
        for( j = 1; j < size; ++j )
        {
            if( coord[ j*DIMS1 + i ] < m[i][0] )
            {
                m[i][0] = coord[ j*DIMS1 + i];
            }
            else if( coord[ j*DIMS1 + i ] > m[i][1] )
            {
                m[i][1] = coord[ j*DIMS1 + i];
            }
            g[i] += coord[ j*DIMS1 + i ];
        }
    }
    g[0] /= size;
    g[1] /= size;
    g[2] /= size;
    if( bc == BC_BOX || bc == BC_SPHERE )
    {
        begin[0] = begin[1] = begin[2] = 0;
    }
    else
    {
        UNERR( "Wrong BC value");
    }
    for( i = 0; i < 3; ++i )
    {
        if( bc == BC_BOX && (m[i][1]-m[i][0]) > box[i] )
        {
            UNERR("Particles don't fit in the box");
        }
        m[i][0] = (m[i][1]+m[i][0])/2;
        m[i][0] = begin[i] - g[i];
    }
    for( j = 0; j < size; ++j )
    {
        for( i = 0; i < 3; ++i )
        {
            coord[ j * DIMS1 + i ] += m[i][0];
            if( bc == BC_BOX )
            {
                if( ABSD( coord[ j * DIMS1 + i ] ) > box[i]/2 + EPSILON  )
                {
                    UNERR("Wrong correction");
                }
            }
        }
    }
    if( !iscorrect() ) UNERR("Wrong correction");
}

void trans_to_box( DOUBLE* v )
{
    if( bc != BC_BOX ) return;
    v[0] -= box[0] * FLOORD( v[0] * inv_L[0] + 0.5f );
    v[1] -= box[1] * FLOORD( v[1] * inv_L[1] + 0.5f );
    v[2] -= box[2] * FLOORD( v[2] * inv_L[2] + 0.5f );
}

DOUBLE ABSD(DOUBLE);

void trans_dist_vec( DOUBLE* v1, DOUBLE* v2, DOUBLE* v )
{
    v[0] = v2[0] - v1[0];
    v[1] = v2[1] - v1[1];
    v[2] = v2[2] - v1[2];
    if( bc != BC_BOX ) return;
    /* function driven */
    v[0] -= box[0] * FLOORD( v[0] * inv_L[0] + 0.5f );
    v[1] -= box[1] * FLOORD( v[1] * inv_L[1] + 0.5f );
    v[2] -= box[2] * FLOORD( v[2] * inv_L[2] + 0.5f );
}

void dist_vec( DOUBLE* v1, DOUBLE* v2, DOUBLE* out )
{
    out[0] = v2[0] - v1[0];
    out[1] = v2[1] - v1[1];
    out[2] = v2[2] - v1[2];
}

INT iscorrect()
{
    if( bc == BC_SPHERE )
    {
        int i;
        DOUBLE sr = sphere_radius*sphere_radius;
        for( i = 0; i < size; ++i )
        {
            DOUBLE distance = sq3v( coord + i * DIMS1 ); /* from (0,0,0) */
            if( distance > sr )
            {
                return 0; /* not correct if outside sphere */
            }
        }
    }

    return 1;
}
