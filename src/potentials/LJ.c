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

#include "LJ.h"

#include <stdio.h>
#include <math.h>

#include "../err.h"
#include "../input.h"
#include "../math_help.h"

MAKE_STR_IN(DOUBLE,alpha_lj,4.0,"factor for L-J interactions scaling")
MAKE_STR_IN(DOUBLE,cutoff_lj,0.0,"cutoff for L-J interactions")
MAKE_STR_IN(YESNO,lj_6_term,1,"switch the attractive part of the L-J potential on/off")

INT LJ_force(DOUBLE r_norm, DOUBLE* e, DOUBLE* F1, DOUBLE* F2, DOUBLE* lj1, DOUBLE* lj2, DOUBLE* _E )
{
    DOUBLE sigma_ij = lj1[ LJ_SIGMA_OFF ] + lj2[ LJ_SIGMA_OFF ];
    DOUBLE epsilon_ij = SQRTD(lj1[ LJ_EPS_OFF ] * lj2[ LJ_EPS_OFF ]);
    DOUBLE sig_r_6 = sigma_ij / r_norm;
    DOUBLE F;

    sig_r_6 *= sig_r_6;
    sig_r_6 *= sig_r_6 * sig_r_6; /* much faster than pow( x, 6 ) */

    if(_E) _E[0] += alpha_lj * epsilon_ij * ( sig_r_6 * ( sig_r_6 - lj_6_term ) );
    F = alpha_lj * epsilon_ij * 6 * ( sig_r_6 * ( 2 * sig_r_6 - lj_6_term ) ) / r_norm;

    F1[0] -= e[0] * F;
    F1[1] -= e[1] * F;
    F1[2] -= e[2] * F;

    F2[0] += e[0] * F;
    F2[1] += e[1] * F;
    F2[2] += e[2] * F;

    return (sigma_ij > r_norm) ||
           isnan(F)    || isinf(F)    ||
           isnan(e[0]) || isinf(e[0]) ||
           isnan(e[1]) || isinf(e[1]) ||
           isnan(e[2]) || isinf(e[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

INT LJ_force_scale(DOUBLE r_norm, DOUBLE* e, DOUBLE* F1, DOUBLE* F2, DOUBLE* lj1, DOUBLE* lj2, DOUBLE scale, DOUBLE* _E )
{
    DOUBLE sigma_ij = lj1[ LJ_SIGMA_OFF ] + lj2[ LJ_SIGMA_OFF ];
    DOUBLE epsilon_ij = SQRTD(lj1[ LJ_EPS_OFF ] * lj2[ LJ_EPS_OFF ]);
    DOUBLE sig_r_6 = sigma_ij / r_norm;
    DOUBLE F = scale;

    sig_r_6 *= sig_r_6;
    sig_r_6 *= sig_r_6 * sig_r_6; /* much faster than pow( x, 6 ) */

    if(_E) _E[0] += scale*alpha_lj * epsilon_ij * ( sig_r_6 * ( sig_r_6 - lj_6_term ) );
    F *= alpha_lj * epsilon_ij * 6 * ( sig_r_6 * ( 2 * sig_r_6 - lj_6_term ) ) / r_norm;

    F1[0] -= e[0] * F;
    F1[1] -= e[1] * F;
    F1[2] -= e[2] * F;

    F2[0] += e[0] * F;
    F2[1] += e[1] * F;
    F2[2] += e[2] * F;

    return (sigma_ij > r_norm)  ||
           isnan(F)    || isinf(F)    ||
           isnan(e[0]) || isinf(e[0]) ||
           isnan(e[1]) || isinf(e[1]) ||
           isnan(e[2]) || isinf(e[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

INT LJ_force_check(DOUBLE r_sq, DOUBLE* lj1, DOUBLE* lj2 )
{
    DOUBLE sigma_ij = lj1[ LJ_SIGMA_OFF ] + lj2[ LJ_SIGMA_OFF ];
    return (sigma_ij*sigma_ij > r_sq);
}
