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

#ifndef MATH_HELP_H
#define	MATH_HELP_H

#include "data.h"

#ifndef M_PI
#define M_PI 3.14159265358
#endif

#ifdef	__cplusplus
extern "C" {
#endif

#ifndef HAVE_ISNAN
int isnan(DOUBLE a);
#endif

#ifndef HAVE_ISINF
int isinf(DOUBLE a);
#endif

/*! Minimum of two numbers. */
#define MIN(a,b) (((a)<(b))?(a):(b))

/*! Maximum of two numbers. */
#define MAX(a,b) (((a)>(b))?(a):(b))

/*! Square of vector. */
DOUBLE sq3v(DOUBLE* v);

/*! The Euclidean norm of the vector.
 *  \return Norm of 3D vector. */
DOUBLE norm3v(DOUBLE* v);

/*! Computes distance between beads.
 *  \param p1 id of first bead
 *  \param p2 id of seconf bead
 *  \param r_sq square of distance */
void get_v2_sq( int p1, int p2, DOUBLE* r_sq );

/*! Computes distance between beads.
 *  \param p1 id of first bead
 *  \param p2 id of seconf bead
 *  \param e1 unit vector
 *  \param r_norm [out] distance between beads
 *  \param r_sq [out] square of distance */
void get_v2_norms( int p1, int p2, DOUBLE* e1, DOUBLE* r_norm, DOUBLE* r_sq );

/*! Computes normalized distance vectors between beads.
 *  \param p1 id of first bead
 *  \param p2 id of second bead
 *  \param p3 id of third bead
 *  \param e1 unit vector P2P1
 *  \param e2 unit vector P2P3 */
void get_v3_norms( int p1, int p2, int p3, DOUBLE* e1, DOUBLE* e2 );

/*! Normalize vector and computes it length.
 *  \param x1 x coordinate of 3D vector
 *  \param y1 y coordinate of 3D vector
 *  \param z1 z coordinate of 3D vector
 *  \param e1 [out] unit vector of (x1,y1,z1) vector
 *  \param r_norm [out] distance between beads
 *  \param r_sq [out] square of distance */
void get_ev2_norms( DOUBLE x1, DOUBLE y1, DOUBLE z1, DOUBLE* e1, DOUBLE* r_norm, DOUBLE* r_sq );

/*! Normalize vector and computes it length.
 *  \param p1 id of first bead
 *  \param p2 id of seconf bead
 *  \param e1 [out] unit vector of (x1,y1,z1) vector
 *  \param r_norm [out] distance between beads
 *  \param r_sq [out] square of distance */
void get_v2_norms_no( int p1, int p2, DOUBLE* e1, DOUBLE* r_norm, DOUBLE* r_sq );

/*! Calculate the cosine of the angle between the vectors.
    \param e1 unit vector
    \param e2 unit vector
    \return the cosine value */
DOUBLE cos_vs( DOUBLE* e1, DOUBLE* e2 );

/*! Absolute value function.
    \return the absolute value of parameter*/
DOUBLE ABSD(DOUBLE);

/*! Maximum function
 *  \return max of two numbers. */
DOUBLE MAXD(DOUBLE,DOUBLE);

/*! Minimum function
 *  \return min of two numbers. */
DOUBLE MIND(DOUBLE,DOUBLE);

/*! Cross product of two vectors. */
void cross_product(DOUBLE*,DOUBLE*,DOUBLE*);

/*! Sign function.
 *  \return Sign of number. */
DOUBLE SINGD(DOUBLE);

#ifdef	__cplusplus
}
#endif

#endif	/* MATH_HELP_H */
