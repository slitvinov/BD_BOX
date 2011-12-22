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

#include "dihe.h"

#include <stdlib.h>
#include <math.h>

#include "../err.h"
#include "../math_help.h"
#include "../trans.h"

#define EPSILON 0.000001

MAKE_STR_IN(INT,dihe_size,0,"inner variable")
MAKE_STR_IN(PTR_INT,dihe_pairs,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,dihe_parms,NULL,"inner variable")

DOUBLE get_delta_cos( int m, DOUBLE cos_theta);
DOUBLE get_m_cos( int m, DOUBLE cos_theta);

INT dihe_calc( DOUBLE* F, int pos, DOUBLE* _E )
{
    DOUBLE* F1 = F + DIMS0 * dihe_pairs[ 4 * pos ];
    DOUBLE* F2 = F + DIMS0 * dihe_pairs[ 4 * pos + 1 ];
    DOUBLE* F3 = F + DIMS0 * dihe_pairs[ 4 * pos + 2 ];
    DOUBLE* F4 = F + DIMS0 * dihe_pairs[ 4 * pos + 3 ];
    DOUBLE* ri = coord + 4 * dihe_pairs[ 4 * pos ];
    DOUBLE* rj = coord + 4 * dihe_pairs[ 4 * pos + 1];
    DOUBLE* rk = coord + 4 * dihe_pairs[ 4 * pos + 2];
    DOUBLE* rl = coord + 4 * dihe_pairs[ 4 * pos + 3];
    DOUBLE rim[3];
    DOUBLE inv_nrim;
    DOUBLE rln[3];
    DOUBLE inv_nrln;
    
    DOUBLE sijkj;
    DOUBLE sklkj;
    DOUBLE nkjsq;
    DOUBLE rkj[3];
    DOUBLE rij[3];
    DOUBLE rkl[3];

    DOUBLE cos_theta;
    DOUBLE delta_cos;

    DOUBLE Fi[3];
    DOUBLE Fj[3];
    DOUBLE Fk[3];
    DOUBLE Fl[3];

    DOUBLE minus_k_cos_delta = dihe_parms[ 2*pos ];
    DOUBLE tmp;

    trans_dist_vec( rk, rj, rkj );
    nkjsq = sq3v(rkj);
    
    trans_dist_vec( ri, rj, rij );

    trans_dist_vec( rk, rl, rkl );

    sijkj = (rij[0] * rkj[0] + rij[1] * rkj[1] + rij[2] * rkj[2] ) / nkjsq;
    sklkj = (rkl[0] * rkj[0] + rkl[1] * rkj[1] + rkl[2] * rkj[2] ) / nkjsq;
    
    /*(33)*/
    rim[0] = -rij[0] + sijkj * rkj[0];
    rim[1] = -rij[1] + sijkj * rkj[1];
    rim[2] = -rij[2] + sijkj * rkj[2];
    tmp = norm3v(rim);
    if( ABSD(tmp) < EPSILON ) return 0;
    inv_nrim = 1.0f / tmp;
    /*(34)*/
    rln[0] = rkl[0] - sklkj * rkj[0];
    rln[1] = rkl[1] - sklkj * rkj[1];
    rln[2] = rkl[2] - sklkj * rkj[2];
    tmp = norm3v(rln);
    if( ABSD(tmp) < EPSILON ) return 0;
    inv_nrln = 1.0f / tmp;
    /*(35)*/
    cos_theta = (rim[0]*rln[0]+rim[1]*rln[1]+rim[2]*rln[2])*inv_nrln*inv_nrim;
    if(cos_theta < -1) cos_theta = -1;
    if(cos_theta > 1) cos_theta = 1;
    /*V*/
    if(_E) _E[0] += ABSD(minus_k_cos_delta) * ( 1 + ABSD(-minus_k_cos_delta)*get_m_cos( (int)dihe_parms[ 2*pos + 1], cos_theta) );
    delta_cos = get_delta_cos( (int)dihe_parms[ 2*pos + 1], cos_theta);

    /*(51)*/
    tmp = minus_k_cos_delta * delta_cos * inv_nrim;
    Fi[0] = tmp * ( rln[0]*inv_nrln - rim[0]*inv_nrim*cos_theta );
    Fi[1] = tmp * ( rln[1]*inv_nrln - rim[1]*inv_nrim*cos_theta );
    Fi[2] = tmp * ( rln[2]*inv_nrln - rim[2]*inv_nrim*cos_theta );
    /*(52)*/
    tmp = minus_k_cos_delta * delta_cos * inv_nrln;
    Fl[0] = tmp * ( rim[0]*inv_nrim - rln[0]*inv_nrln*cos_theta );
    Fl[1] = tmp * ( rim[1]*inv_nrim - rln[1]*inv_nrln*cos_theta );
    Fl[2] = tmp * ( rim[2]*inv_nrim - rln[2]*inv_nrln*cos_theta );
    /*(53)*/
    Fj[0] = (sijkj - 1)*Fi[0] - sklkj * Fl[0];
    Fj[1] = (sijkj - 1)*Fi[1] - sklkj * Fl[1];
    Fj[2] = (sijkj - 1)*Fi[2] - sklkj * Fl[2];
    /*(54)*/
    Fk[0] = -Fi[0]-Fj[0]-Fl[0];
    Fk[1] = -Fi[1]-Fj[1]-Fl[1];
    Fk[2] = -Fi[2]-Fj[2]-Fl[2];

    F1[0] += Fi[0];
    F1[1] += Fi[1];
    F1[2] += Fi[2];

    F2[0] += Fj[0];
    F2[1] += Fj[1];
    F2[2] += Fj[2];

    F3[0] += Fk[0];
    F3[1] += Fk[1];
    F3[2] += Fk[2];

    F4[0] += Fl[0];
    F4[1] += Fl[1];
    F4[2] += Fl[2];

    return isnan(Fk[0]) || isinf(Fk[0]) ||
           isnan(Fk[1]) || isinf(Fk[1]) ||
           isnan(Fk[2]) || isinf(Fk[2]) ||
           isnan(Fi[0]) || isinf(Fi[0]) ||
           isnan(Fi[1]) || isinf(Fi[1]) ||
           isnan(Fi[2]) || isinf(Fi[2]) ||
           isnan(Fj[0]) || isinf(Fj[0]) ||
           isnan(Fj[1]) || isinf(Fj[1]) ||
           isnan(Fj[2]) || isinf(Fj[2]) ||
           isnan(Fl[0]) || isinf(Fl[0]) ||
           isnan(Fl[1]) || isinf(Fl[1]) ||
           isnan(Fl[2]) || isinf(Fl[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

DOUBLE get_delta_cos( int m, DOUBLE cos_theta)
{
    switch(m)
    {
        /* equations (44) - (50) from manual */
        case 1: return 1;
        case 2: return  4 * cos_theta;
        case 3: return 12 * cos_theta*cos_theta - 3;
        case 4: return cos_theta*(32 * cos_theta*cos_theta - 16);
        case 5: return cos_theta*cos_theta*(80*cos_theta*cos_theta-60)+5;
        case 6: return cos_theta*(192*cos_theta*cos_theta*(cos_theta*cos_theta-1) + 36);
        default: return 0;
    }
}

DOUBLE get_m_cos( int m, DOUBLE cos_theta)
{
    switch(m)
    {
        case 1: return cos_theta;
        case 2: return 2 * cos_theta*cos_theta - 1;
        case 3: return 4 * cos_theta*cos_theta*cos_theta - 3*cos_theta;
        case 4: return 8 * cos_theta*cos_theta*(cos_theta*cos_theta-1)+1;
        case 5: return cos_theta*(cos_theta*cos_theta*(16*cos_theta*cos_theta-20) + 5);
        case 6: return cos_theta*cos_theta*(cos_theta*cos_theta*(32*cos_theta*cos_theta-48) + 18)-1;
        default: return 1;
    }
}

static INT capacity = 0;
void add_dihe( INT id1, INT id2, INT id3, INT id4, DOUBLE k_theta, DOUBLE m, DOUBLE delta )
{
    if( capacity == dihe_size )
    {
        capacity = capacity ? 2*capacity : 128;
        dihe_pairs = (INT*) realloc( dihe_pairs, sizeof(INT) * 4 * capacity );
        CHMEM(dihe_pairs);
        dihe_parms = (DOUBLE*) realloc( dihe_parms, sizeof(DOUBLE) * 2 * capacity );
        CHMEM(dihe_parms);
    }
    dihe_pairs[4*dihe_size] = id1;
    dihe_pairs[4*dihe_size + 1] = id2;
    dihe_pairs[4*dihe_size + 2] = id3;
    dihe_pairs[4*dihe_size + 3] = id4;
    dihe_parms[2*dihe_size] = -k_theta * delta;
    dihe_parms[2*dihe_size + 1] = m;
    dihe_size++;
}

void end_dihe()
{
    if( dihe_size )
    {
        dihe_parms = (DOUBLE*) realloc( dihe_parms, sizeof(DOUBLE) * 2 * dihe_size );
        MAKE_STR_SIZE(dihe_parms,2*dihe_size)
        dihe_pairs = (INT*) realloc( dihe_pairs, sizeof(INT) * 4 * dihe_size );
        MAKE_STR_SIZE(dihe_pairs,4*dihe_size)
    }
}

void free_dihe()
{
    if( dihe_parms )
    {
        free( dihe_parms );
        dihe_parms = NULL;
    }
    if( dihe_pairs )
    {
        free( dihe_pairs );
        dihe_pairs = NULL;
    }
}
