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

#include "pwell.h"

#include "../trans.h"
#include "../math_help.h"

#define EPSILON 0.000001

MAKE_STR_IN(YESNO,  pwell,   0, "yes/no - switch the potential well on/off")
MAKE_STR_IN(DOUBLE, pwell_A, 1, "the magnitude of the potential well")

INT pwell_ext( DOUBLE* coord, DOUBLE* F )
{
  DOUBLE e[3];
  e[0] = 0.0;
  e[1] = pwell_A * coord[1];
  e[2] = 0.0;
  
  F[0] -= e[0];
  F[1] -= e[1];
  F[2] -= e[2];
  return 0;
}
