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

#include "angle.h"

#include <stdlib.h>
#include <math.h>

#include "../err.h"
#include "../math_help.h"

#define EPSILON 0.000001

MAKE_STR_IN(INT,angle_size,0,"inner variable")
MAKE_STR_IN(PTR_INT,angle_pairs,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,angle_phis,NULL,"inner variable")

INT angle_calc( DOUBLE* F, int pos, DOUBLE* _E  )
{
    DOUBLE phi_0 = angle_phis[ 2*pos ];
    DOUBLE k_phi = angle_phis[ 2*pos + 1 ];
    DOUBLE* F1 = F + DIMS0 * angle_pairs[ 3 * pos ];
    DOUBLE* F2 = F + DIMS0 * angle_pairs[ 3 * pos + 1 ];
    DOUBLE* F3 = F + DIMS0 * angle_pairs[ 3 * pos + 2 ];
    DOUBLE Fi[3];
    DOUBLE Fk[3];
    DOUBLE e1[3];
    DOUBLE e2[3];
    DOUBLE cos_phi;
    DOUBLE phi;
    get_v3_norms( angle_pairs[ 3 * pos ], angle_pairs[ 3 * pos + 1 ], angle_pairs[ 3 * pos + 2 ], e1, e2 );
    cos_phi = cos_vs( e1, e2 );
    if(cos_phi>=1-EPSILON || cos_phi<=-1+EPSILON) return 0;
    phi  = ACOSD( cos_phi );
    k_phi *= ( phi - phi_0 );
    if(_E) _E[0] += 0.5f * k_phi * ( phi - phi_0 );
    k_phi /= SQRTD( 1 - cos_phi*cos_phi );

    Fi[0] = -k_phi * (e2[0] - e1[0] * cos_phi);
    Fi[1] = -k_phi * (e2[1] - e1[1] * cos_phi);
    Fi[2] = -k_phi * (e2[2] - e1[2] * cos_phi);

    Fk[0] = -k_phi * (e1[0] - e2[0] * cos_phi);
    Fk[1] = -k_phi * (e1[1] - e2[1] * cos_phi);
    Fk[2] = -k_phi * (e1[2] - e2[2] * cos_phi);

    F1[0] += Fi[0];
    F1[1] += Fi[1];
    F1[2] += Fi[2];

    F2[0] += -Fi[0]-Fk[0];
    F2[1] += -Fi[1]-Fk[1];
    F2[2] += -Fi[2]-Fk[2];

    F3[0] += Fk[0];
    F3[1] += Fk[1];
    F3[2] += Fk[2];

    return isnan(Fk[0]) || isinf(Fk[1]) ||
           isnan(Fk[1]) || isinf(Fk[1]) ||
           isnan(Fk[2]) || isinf(Fk[2]) ||
           isnan(Fi[0]) || isinf(Fi[0]) ||
           isnan(Fi[1]) || isinf(Fi[1]) ||
           isnan(Fi[2]) || isinf(Fi[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

static INT capacity = 0;
void add_angle( INT id1, INT id2, INT id3, DOUBLE phi_0, DOUBLE k_phi )
{
    if( capacity == angle_size )
    {
        capacity = capacity ? 2*capacity : 128;
        angle_pairs = (INT*) realloc( angle_pairs, sizeof(INT) * 3 * capacity );  CHMEM(angle_pairs);
        angle_phis = (DOUBLE*) realloc( angle_phis, sizeof(DOUBLE) * 2 * capacity );  CHMEM(angle_phis);
    }
    angle_pairs[3*angle_size] = id1;
    angle_pairs[3*angle_size + 1] = id2;
    angle_pairs[3*angle_size + 2] = id3;
    angle_phis[2*angle_size] = phi_0 * M_PIF / 180;
    angle_phis[2*angle_size + 1] = k_phi;
    angle_size++;
}

void end_angle()
{
    if( angle_size )
    {
        angle_phis = (DOUBLE*) realloc( angle_phis, sizeof(DOUBLE) * 2 * angle_size );
        MAKE_STR_SIZE(angle_phis,2 * angle_size)
        angle_pairs = (INT*) realloc( angle_pairs, sizeof(INT) * 3 * angle_size );
        MAKE_STR_SIZE(angle_pairs,3 * angle_size)
    }
}
void free_angle()
{
    if( angle_phis )
    {
        free( angle_phis );
        angle_phis = NULL;
    }
    if( angle_pairs )
    {
        free( angle_pairs );
        angle_pairs = NULL;
    }
}
