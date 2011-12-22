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

#ifndef _ANGLE_H
#define	_ANGLE_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Number of entries in angle_phis and angle_pairs array.*/
extern INT angle_size;
MAKE_STR_DEC_INNER(INT,angle_size)

/*! Array of interacting beads. */
extern PTR_INT angle_pairs;
MAKE_STR_DEC_INNER(PTR_INT,angle_pairs)

/*! Array of interacting beads parameters. */
extern PTR_DOUBLE angle_phis;
MAKE_STR_DEC_INNER(PTR_DOUBLE,angle_phis)

/*! Compute forces on beads.
 *  \param F pointer to output forces array
 *  \param pos position of interacting pair in angle_phis
    \param E [out] energy of angle deformation */
INT angle_calc( DOUBLE* F, int pos, DOUBLE* E  );

/*! Add angle desciption to angle_phis array.
    \param id1 id of first bead
    \param id2 id of second(middle) bead
    \param id3 id of third bead
    \param phi_0 value of angle
    \param k_phi force constant of angle deformation */
void add_angle( INT id1, INT id2, INT id3, DOUBLE phi_0, DOUBLE k_phi );

/*! Trim size of angle_phis. */
void end_angle();

/*! Free resources. */
void free_angle();

#ifdef	__cplusplus
}
#endif

#endif	/* ANGLE_H */
