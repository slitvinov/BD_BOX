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

#ifndef PWELL_H
#define	PWELL_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Switch on/off whether the bounding sphere force should be used. */
extern YESNO pwell;
MAKE_STR_DEC(YESNO,pwell)

/*! The amplitude of the bounding sphere force. */
extern DOUBLE pwell_A;
MAKE_STR_DEC(DOUBLE,pwell_A)

/*! Compute sphere boundary forces values.
 *  \param coord coordinate array of bead
 *  \param F [out] forces value */
INT pwell_ext( DOUBLE* coord, DOUBLE* F );

#ifdef	__cplusplus
}
#endif

#endif	/* PWELL_H */
