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

#ifndef TRANS_H
#define	TRANS_H

#include "data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Boundary conditions, currently none, a rectangular box or a reflective sphere. */
extern NONE_CUBIC_SPHERE bc;
MAKE_STR_DEC(NONE_CUBIC_SPHERE,bc)

/*! Radius of the reflective sphere enclosing the studied system. */
extern DOUBLE sphere_radius;
MAKE_STR_DEC(DOUBLE,sphere_radius)

/*! None periodic boundary condition.
    \sa bs */
#define BC_NONE 0

/*! Cubic periodic boundary condition.
    \sa bs
    \sa box */
#define BC_BOX 1

/*! Reflective sphere.
    \sa bs
    \sa sphere_radius
    \sa iscorrect */
#define BC_SPHERE 2

/*! Potential well.
    \sa bs
    \sa sphere_radius
    \sa iscorrect */
#define BC_PWELL 3


/*! String value in prm or argmunent line option indicate 'none' bs. */
#define BC_NONE_STR "none"

/*! String value in prm or argmunent line option indicate 'periodic' bs. */
#define BC_BOX_STR "periodic"

/*! String value in prm or argmunent line option indicate 'sphere' bs. */
#define BC_SPHERE_STR "sphere"

/*! String value in prm or argmunent line option indicate 'sphere' bs. */
#define BC_PWELL_STR "pwell"

/*! Center coordinations in box. */
void center_coords();

/*! Vector of inverted box sizes. */
extern DOUBLE inv_L[4];

/*! Function keeps beads in the primary box.
 *  \param v vector with 3 coordinates to transform in box, for cubic bc */
void trans_to_box( DOUBLE* v );

/*! Compute vector connecting nearest images.
 *  \param v1 [in] position vector of first bead
 *  \param v2 [in] position vector of second bead
 *  \param out [out] transformed vector connecting beads */
void trans_dist_vec( DOUBLE* v1, DOUBLE* v2, DOUBLE* out );

/*! Computes distance between points.
 *  \param v1 [in] position vector of first bead
 *  \param v2 [in] position vector of second bead
 *  \param out [out] distance vector*/
void dist_vec( DOUBLE* v1, DOUBLE* v2, DOUBLE* out );

/*! Check for sphere BC if all beads are in sphere of sphere_radius. */
INT iscorrect();

#ifdef	__cplusplus
}
#endif

#endif	/* TRANS_H */
