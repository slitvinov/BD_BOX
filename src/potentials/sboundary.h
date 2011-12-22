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

#ifndef SBOUNDARY_H
#define	SBOUNDARY_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Switch on/off whether the bounding sphere force should be used. */
extern YESNO sboundary;
MAKE_STR_DEC(YESNO,sboundary)

/*! The amplitude of the bounding sphere force. */
extern DOUBLE sboundary_A;
MAKE_STR_DEC(DOUBLE,sboundary_A)

/*! The bounding sphere force, n value. */
extern INT sboundary_n;
MAKE_STR_DEC(INT,sboundary_n)

/*! The bounding sphere force, cutoff. */
extern DOUBLE sboundary_cutoff;
MAKE_STR_DEC(DOUBLE, sboundary_cutoff)

/*! Compute sphere boundary forces values.
 *  \param coord coordinate array of bead
 *  \param F [out] forces value */
INT sboundary_ext( DOUBLE* coord, DOUBLE* F );

#ifdef	__cplusplus
}
#endif

#endif	/* SBOUNDARY_H */
