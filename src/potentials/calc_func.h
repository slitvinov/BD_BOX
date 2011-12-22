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

#ifndef _CALC_FUNC_H
#define	_CALC_FUNC_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Array of energies. */
extern DOUBLE* E;
extern DOUBLE** ES;

/*! Array of forces. */
extern DOUBLE** F;
extern INT* repeat;

/*! Scaling factor for Lennard-Jones interactions between bonded partners 1-2. */
extern DOUBLE bond_lj_scale;
MAKE_STR_DEC(DOUBLE,bond_lj_scale)

/*! Scaling factor for electrostatic interactions between bonded partners 1-2. */
extern DOUBLE bond_c_scale;
MAKE_STR_DEC(DOUBLE,bond_c_scale)

/*! Switch on/off checking overlap. */
extern YESNO check_overlap;
MAKE_STR_DEC( YESNO, check_overlap )

/*!
 * Compute forces and energies for beads.
 * Can use CUDA, OpenMP, MPI technology, depends on compilation.
 * \param E optional array of energies
 * \param curr_time current time
 * \return nonzero value for wrong conformation
 */
INT compute_forces( DOUBLE* E, DOUBLE curr_time );
void init_forces();
void free_forces();

/*!
 * Compute bonded interactions.
 * \param _E energy array
 * \param tid id of thread
 * \param curr_time current time
 * \return one if overlap, zero otherwise
 */
INT bonded_inter( DOUBLE* _E, INT tid, DOUBLE curr_time );

/*!
 * Compute bonded interactions in all threads.
 * \param _E energy array
 * \param tid id of thread
 * \param curr_time current time
 * \return one if overlap, zero otherwise
 */
INT bonded_inter_parallel( DOUBLE* _E, DOUBLE curr_time );

/*!
 * Initialization of interactions on CPU.
 * Cleaning energies and forces.
 * \param _E energy array
 * \return one if exists bead outside sphere, zero otherwise
 */
INT init_iter_forces( DOUBLE* _E );

/*!
 * Sumarize energies, forces and repeats in OpenMP version.
 * \param _E energy array
 * \return one if step should be repeated
 */
INT collect_results( DOUBLE* _E );

/*!
 * Compute optimal cutt-off. Used in subdivion algorithms.
 */
void init_cutoff();

/*! Init forces and energies. */
void init_Fs();

/*! Free resources. */
void free_Fs();

#ifdef	__cplusplus
}
#endif

#endif	/* CALC_FUNC_H */
