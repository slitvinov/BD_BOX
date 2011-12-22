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

#include "electro.h"

#include <math.h>

#include "../math_help.h"

MAKE_STR_IN(DOUBLE,epsilon_c,78.54f,"solvent dielectric constant")
MAKE_STR_IN(DOUBLE,kappa_c,0.1f,"inverse of the screening length")
MAKE_STR_IN(DOUBLE,gamma_c,331.842f,"scaling factor for electrostatic interactions")
MAKE_STR_IN(YESNO,elec,1,"switch electrostatic interactions on/off")
MAKE_STR_IN(DOUBLE,cutoff_c,0,"cutoff for electrostatic interactions")

INT electro_force( DOUBLE r_norm, DOUBLE r_sq, DOUBLE* e,
                    DOUBLE* F1, DOUBLE* F2, DOUBLE Q1, DOUBLE Q2, DOUBLE* _E )
{
    DOUBLE F = (gamma_c / epsilon_c) * Q1 * Q2 * ( 1 + kappa_c * r_norm );
    F *= EXPD( - kappa_c * r_norm ) / r_sq;
    if(_E) _E[0] += (gamma_c / epsilon_c) * Q1 * Q2 * EXPD( - kappa_c * r_norm ) / r_norm;
    
    F1[0] -= e[0] * F;
    F1[1] -= e[1] * F;
    F1[2] -= e[2] * F;

    F2[0] += e[0] * F;
    F2[1] += e[1] * F;
    F2[2] += e[2] * F;
    
    return isnan(F)    || isinf(F)    ||
           isnan(e[0]) || isinf(e[0]) ||
           isnan(e[1]) || isinf(e[1]) ||
           isnan(e[2]) || isinf(e[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

INT electro_force_scale( DOUBLE r_norm, DOUBLE r_sq, DOUBLE* e,
                    DOUBLE* F1, DOUBLE* F2, DOUBLE Q1, DOUBLE Q2, DOUBLE scale, DOUBLE* _E )
{
    DOUBLE F = (gamma_c / epsilon_c) * Q1 * Q2 * ( 1 + kappa_c * r_norm );
    F *= EXPD( - kappa_c * r_norm ) / r_sq;
    F *= scale;
    if(_E) _E[0] += scale*(gamma_c / epsilon_c) * Q1 * Q2 * EXPD( - kappa_c * r_norm ) / r_norm;
    F1[0] -= e[0] * F;
    F1[1] -= e[1] * F;
    F1[2] -= e[2] * F;

    F2[0] += e[0] * F;
    F2[1] += e[1] * F;
    F2[2] += e[2] * F;

    return isnan(F)    || isinf(F)    ||
           isnan(e[0]) || isinf(e[0]) ||
           isnan(e[1]) || isinf(e[1]) ||
           isnan(e[2]) || isinf(e[2]) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}
