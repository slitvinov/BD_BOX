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

#ifndef _ELECTRO_EXT_H
#define	_ELECTRO_EXT_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Switch the external electric field on/off. */
extern YESNO E_ext;
MAKE_STR_DEC(YESNO,E_ext)

/*! Electrostatic filed, units conversion factor, for example,
 *  if input units are V/m, this value should be 0.230710 * 10^-8*/
extern DOUBLE E_factor;
MAKE_STR_DEC(DOUBLE,E_factor)

/*! A choice between AC, DC or RF electric fields. */
extern DOUBLE E_magn;
MAKE_STR_DEC(DOUBLE,E_magn)

/*! A choice between AC, DC or RF electric fields */
extern DC_AC_RF E_type;
MAKE_STR_DEC(DC_AC_RF,E_type)

/*! Frequency of the external electric field, applicable in case of AC or RF. */
extern DOUBLE E_freq;
MAKE_STR_DEC(DOUBLE,E_freq)

/*! The direction in which the external electric field is applied (DC/AC/RF)*/
extern CHAR E_dir1;
MAKE_STR_DEC(CHAR,E_dir1)

/*! The direction in which the external electric field is applied (RF)*/
extern CHAR E_dir2;
MAKE_STR_DEC(CHAR,E_dir2)

/*! Computes Forces acting on charged subunits due to the external electric field.
    \param Q charge of bead
    \param F [out] computed force
    \param t current time */
INT electro_ext_force( DOUBLE Q, DOUBLE* F, DOUBLE t );

#ifdef	__cplusplus
}
#endif

#endif	/* ELECTRO_EXT_H */
