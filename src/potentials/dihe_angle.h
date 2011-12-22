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

#ifndef DIHE_ANGLE_H
#define	DIHE_ANGLE_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Number of entries in dihe_angle_parms and dihe_angle_pairs array. */
extern INT dihe_angle_size;
MAKE_STR_DEC(INT,dihe_angle_size)

/*! Array of interacting beads. */
extern PTR_INT dihe_angle_pairs;
MAKE_STR_DEC(PTR_INT,dihe_angle_pairs)

/*! Array of interacting beads parameters. */
extern PTR_DOUBLE dihe_angle_parms;
MAKE_STR_DEC(PTR_DOUBLE,dihe_angle_parms)

/*! Compute forces on beads.
 *  \param F pointer to output forces array
 *  \param pos position of interacting pair in dihe_phis
    \param E [out] energy of dihedral deformation */
INT dihe_calc_angle( DOUBLE* F, int pos, DOUBLE* E );

/*! Add angle desciption to dihe_parms array.
    \param id1 id of first bead
    \param id2 id of second bead
    \param id3 id of third bead
    \param id4 id of third bead
    \param k_theta force constant of dihedral angle deformation
 *  \param theta value of equilibrium dihedral angle value */
void add_dihe_angle( INT id1, INT id2, INT id3, INT id4, DOUBLE k_theta, DOUBLE theta );

/*! Trim size of angle_phis. */
void end_dihe_angle();

/*! Free resources. */
void free_dihe_angle();

#ifdef	__cplusplus
}
#endif

#endif	/* DIHE_ANGLE_H */
