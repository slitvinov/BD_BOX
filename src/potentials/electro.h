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

#ifndef ELECTRO_H
#define	ELECTRO_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Solvent dielectric constant. */
extern DOUBLE epsilon_c;
MAKE_STR_DEC(DOUBLE,epsilon_c)

/*! Inverse of screening length. */
extern DOUBLE kappa_c;
MAKE_STR_DEC(DOUBLE,kappa_c)
    
/*! Scaling factor electrostatic interactions. */
extern DOUBLE gamma_c;
MAKE_STR_DEC(DOUBLE,gamma_c)

/* Switch electrostatic interactions on/off. */
extern YESNO elec;
MAKE_STR_DEC(YESNO,elec)

/*! Cutoff for electrostatic interactions. */
MAKE_STR_DEC(DOUBLE,cutoff_c)

/*! Computes eletrostatic interaction.
    \param r_norm distance beetwen beads
    \param r_sq square of r_norm
    \param e
    \param F1 [out] force vector for first bead
    \param F2 [out] force vector for second bead
    \param Q1 charge of first bead
    \param Q2 charge of second bead
    \param _E [out] energy value */
INT electro_force( DOUBLE r_norm, DOUBLE r_sq, DOUBLE* e, DOUBLE* F1, DOUBLE* F2, DOUBLE Q1, DOUBLE Q2, DOUBLE* _E );

/*! Computes eletrostatic interaction with scaling.
    \param r_norm distance beetwen beads
    \param r_sq square of r_norm
    \param e
    \param F1 [out] force vector for first bead
    \param F2 [out] force vector for second bead
    \param Q1 charge of first bead
    \param Q2 charge of second bead
    \param _E [out] energy value
    \param scale scaling factor for connected beads */
INT electro_force_scale( DOUBLE r_norm, DOUBLE r_sq, DOUBLE* e, DOUBLE* F1, DOUBLE* F2, DOUBLE Q1, DOUBLE Q2, DOUBLE scale, DOUBLE* _E );

#ifdef	__cplusplus
}
#endif

#endif	/* _ELECTRO_H */
