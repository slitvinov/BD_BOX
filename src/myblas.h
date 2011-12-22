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

#ifndef _CHOLESKY_H
#define	_CHOLESKY_H

#include "data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Wrapper for dpotrf from diffrent implementation and libraries.
    \param A matrix in column-wise order
    \param size dimension of matrix
    \param LD leading deminsion of matrix */
INT my_dpotrf( DOUBLE* A, INT size, INT LD );

/*! Wrapper for dtrmv from diffrent implementation and libraries.
    \param A matrix in column-wise order
    \param v vector to multiplicate
    \param o output vector A*x
    \param n dimension of matrix
    \param LD leading deminsion of matrix */
void my_dtrmv( DOUBLE* A, DOUBLE* v, DOUBLE* o, INT n, INT LD );

/*! Init cholesky decomposition. Alloc additional vector. */
void init_cholesky();

/*! Free resource used in cholesky decomposition. */
void free_cholesky();

/*! Stub, only declaration. */
void dpotrf_( const char* uplo, int* n, double* a, int* lda, int* info );

/*! Stub, only declaration. */
void dtrmv_(const char *uplo, const char *transa, char *diag, int *n, double *a, int *lda, double *b, int *incx);

/*! Stub, only declaration. */
void dsymv_(const char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

#ifdef	__cplusplus
}
#endif

#endif	/* CHOLESKY_H */
