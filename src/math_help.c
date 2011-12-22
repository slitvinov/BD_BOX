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

#include "math_help.h"

#include <math.h>

#include "trans.h"

DOUBLE sq3v(DOUBLE* v)
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

DOUBLE norm3v(DOUBLE* v)
{
    return SQRTD( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

void get_v2_norms( int p1, int p2, DOUBLE* e1, DOUBLE* r_norm, DOUBLE* r_sq )
{
    DOUBLE* v1 = coord + p1*DIMS1;
    DOUBLE* v2 = coord + p2*DIMS1;
    DOUBLE inv = 0.0;
    
    trans_dist_vec( v1, v2, e1 );

    *r_sq = sq3v(e1);
    *r_norm = SQRTD( *r_sq );
    inv = 1 / *r_norm;

    e1[0] *= inv;
    e1[1] *= inv;
    e1[2] *= inv;
}

void get_v2_sq( int p1, int p2, DOUBLE* r_sq )
{
    DOUBLE* v1 = coord + p1*DIMS1;
    DOUBLE* v2 = coord + p2*DIMS1;
    DOUBLE e1[3];
    trans_dist_vec( v1, v2, e1 );

    *r_sq = sq3v(e1);
}

void get_v2_norms_no( int p1, int p2, DOUBLE* e1, DOUBLE* r_norm, DOUBLE* r_sq )
{
    DOUBLE* v1 = coord + p1*DIMS1;
    DOUBLE* v2 = coord + p2*DIMS1;
    DOUBLE inv = 0.0;

    e1[0] = v2[0] - v1[0];
    e1[1] = v2[1] - v1[1];
    e1[2] = v2[2] - v1[2];

    *r_sq = sq3v(e1);
    *r_norm = SQRTD( *r_sq );
    inv = 1 / *r_norm;

    e1[0] *= inv;
    e1[1] *= inv;
    e1[2] *= inv;
}

void get_ev2_norms( DOUBLE x1, DOUBLE y1, DOUBLE z1,
                    DOUBLE* e1, DOUBLE* r_norm, DOUBLE* r_sq )
{
    DOUBLE inv = 0.0;
    e1[0] = x1;
    e1[1] = y1;
    e1[2] = z1;
    
    *r_sq = (x1)*(x1) + (y1)*(y1) + (z1)*(z1);
    *r_norm = SQRTD( *r_sq );
    inv = 1 / *r_norm;

    e1[0] *= inv;
    e1[1] *= inv;
    e1[2] *= inv;
}

void get_v3_norms( int p1, int p2, int p3, DOUBLE* e1, DOUBLE* e2 )
{
    DOUBLE* v1 = coord + p1*DIMS1;
    DOUBLE* v2 = coord + p2*DIMS1;
    DOUBLE* v3 = coord + p3*DIMS1;
    DOUBLE nrm1;
    DOUBLE nrm2;

    trans_dist_vec( v1, v2, e1 );

    trans_dist_vec( v3, v2, e2 );

    nrm1 = 1.0f / norm3v(e1);
    nrm2 = 1.0f / norm3v(e2);
    e1[0] *= nrm1;
    e1[1] *= nrm1;
    e1[2] *= nrm1;

    e2[0] *= nrm2;
    e2[1] *= nrm2;
    e2[2] *= nrm2;
}

DOUBLE cos_vs( DOUBLE* e1, DOUBLE* e2 )
{
    DOUBLE cos_v = (e1[0]*e2[0]+e1[1]*e2[1]+e1[2]*e2[2]);
    if( cos_v <= -1.0 ) cos_v = -1.0;
    if( cos_v >=  1.0 ) cos_v =  1.0;
    return cos_v;
}

DOUBLE ABSD(DOUBLE val)
{
    return val < 0.0 ? -val : val;
}

DOUBLE MAXD(DOUBLE v1, DOUBLE v2)
{
    return v1 < v2 ? v2 : v1;
}

DOUBLE MIND(DOUBLE v1, DOUBLE v2)
{
    return v1 > v2 ? v2 : v1;
}

void cross_product(DOUBLE*a,DOUBLE*b,DOUBLE*o)
{
    o[0] = a[1]*b[2]-a[2]*b[1];
    o[1] = a[2]*b[0]-a[0]*b[2];
    o[2] = a[0]*b[1]-a[1]*b[0];
}

DOUBLE SINGD(DOUBLE a)
{
    if( a == 0.0 ) return 0.0;
    if( a < 0.0 ) return -1.0;
    return 1.0;
}
