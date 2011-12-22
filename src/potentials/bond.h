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

#ifndef BOND_H
#define	BOND_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Number of entries in bond_parms and bond_pairs array.*/
extern INT bond_size;
MAKE_STR_DEC(INT,bond_size)

/*! Array of interacting beads. */
extern PTR_INT bond_pairs;
MAKE_STR_DEC(PTR_INT,bond_pairs)

/*! Array of interacting beads parameters. */
extern PTR_DOUBLE bond_parms;
MAKE_STR_DEC(PTR_DOUBLE,bond_parms)

/*! Compute forces on beads.
 *  \param F pointer to output forces array
 *  \param pos position of interacting pair in bond_parms
    \param E interaction energy */
INT bond_calc(DOUBLE* F, int pos, DOUBLE* E );

/*! Check if beads are conneted.
    \return true if connected*/
int are_conn(INT,INT);

/*! Alloc table for connections. */
void alloc_connections();

/*! Create connection graph of beads.
    \param max_group */
INT* get_connection_graph( INT* max_group );

/*! Add bond desciption to bond_parms array.
    \param id1 id of first bead
    \param id2 id of second bead
    \param r_0 value of equilibrium bond length
    \param r_max maximum bond length
    \param H force constant */
void add_bond( INT id1, INT id2, DOUBLE r_0, DOUBLE r_max, DOUBLE H );

/*! Trim size of angle_phis. */
void end_bond();

/*! Free resources. */
void free_bond();

#ifdef	__cplusplus
}
#endif

#endif	/* BOND_H */
