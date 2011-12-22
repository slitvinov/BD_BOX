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

#include "electro_ext.h"

#include <math.h>

#include "../math_help.h"

MAKE_STR_IN(YESNO,E_ext,0,"yes/no - switch the external electric field on/off")
MAKE_STR_IN(DOUBLE,E_factor,1.0,"electrostatic field, units conversion factor")
MAKE_STR_IN(DOUBLE,E_magn,1,"the magnitude of the external electric field")
MAKE_STR_IN(DC_AC_RF,E_type,0,"DC/AC/RF - a choice between AC, DC or RF electric fields")
MAKE_STR_IN(DOUBLE,E_freq,0.0,"frequency of the external electric field, applicable in case of AC or RF")
MAKE_STR_IN(CHAR,E_dir1,'x',"x/y/z - the direction in which the external electric field is applied")
MAKE_STR_IN(CHAR,E_dir2,'y',"x/y/z - the second direction in which the external electric field (only for RF) is applied")

INT electro_ext_force( DOUBLE Q, DOUBLE* F, DOUBLE t )
{
    int pos = E_dir1 - 'x';
    int pos2 = 0;
    switch(E_type)
    {
        case(0): F[pos] += Q * E_factor * E_magn; break;
        case(1): F[pos] += Q * E_factor * E_magn * COSD( E_freq * t); break;
        case(2): pos2 = E_dir2 - 'x';
                 F[pos] += Q * E_factor * E_magn * COSD( E_freq * t);
                 F[pos2] += Q * E_factor * E_magn * SIND( E_freq * t); break;
    }
    return isnan(F[pos]) || isinf(F[pos])
           || (E_type==2 && (isnan(F[pos2]) || isinf(F[pos2])) );
}
