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

#ifndef LJ_H
#define	LJ_H

#include "../data.h"

/*! Offset of sigma in LJ array. */
#define LJ_SIGMA_OFF 0

/*! Offset of epsilon in LJ array. */
#define LJ_EPS_OFF 1

#ifdef	__cplusplus
extern "C" {
#endif

/*! Scaling factor for excluded-volume interactions. */
extern DOUBLE alpha_lj;
MAKE_STR_DEC(DOUBLE,alpha_lj)

/*! cutoff for excluded-volume interactions. */
extern DOUBLE cutoff_lj;
MAKE_STR_DEC(DOUBLE,cutoff_lj)

/*! Switch the attractive part of the L–J potential on/off. */
extern YESNO lj_6_term;
MAKE_STR_DEC(YESNO,lj_6_term)

/*! Compute excluded-volume interaction between pair of beads.
    \param r_norm distance between beads
    \param r_sq square of distance
    \param e unit vector of distance vector
    \param F1 force 3D vector of first bead
    \param F2 force 3D vector of second bead
    \param lj1 sigma and epsilon of first bead
    \param lj2 sigma and epsilon of second bead
    \param E [out] value of energy
    \return true if beads overlap or if arithmetic error occured */
INT LJ_force(DOUBLE r_norm, DOUBLE* e, DOUBLE* F1, DOUBLE* F2, DOUBLE* lj1, DOUBLE* lj2, DOUBLE* E );

/*! Compute scaled excluded-volume interaction between pair of beads.
    \param r_norm distance between beads
    \param r_sq square of distance
    \param e unit vector of distance vector
    \param F1 force 3D vector of first bead
    \param F2 force 3D vector of second bead
    \param lj1 sigma and epsilon of first bead
    \param lj2 sigma and epsilon of second bead
    \param scale scale factor for connected beads
    \param E [out] value of energy
    \return true if beads overlap or if arithmetic error occured */
INT LJ_force_scale(DOUBLE r_norm, DOUBLE* e, DOUBLE* F1, DOUBLE* F2, DOUBLE* lj1, DOUBLE* lj2, DOUBLE scale, DOUBLE* E );

/*! Check overlap.
 *  \param r_norm distance between beads
 *  \param lj1 sigma and epsilon of first bead
 *  \param lj2 sigma and epsilon of second bead
    \return true if beads overlap */
INT LJ_force_check(DOUBLE r_norm, DOUBLE* lj1, DOUBLE* lj2 );

#ifdef	__cplusplus
}
#endif

#endif	/* LJ_H */
