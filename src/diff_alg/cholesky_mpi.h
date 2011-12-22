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

#ifndef CHOLESKY_MPI_H
#define CHOLESKY_MPI_H

#include "../data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Specifies the number of rows in the process grid. */
extern INT MPI_nprow;
MAKE_STR_DEC(INT,MPI_nprow)

/*! Specifies the number of columns in the process grid. */
extern INT MPI_npcol;
MAKE_STR_DEC(INT,MPI_npcol)

/*! Size of block. */
extern INT MPI_block;
MAKE_STR_DEC(INT,MPI_block)

/*! Initialization of cholesky version using mpi. */
void init_cholesky_mpi();

/*! Free resources. */
void free_cholesky_mpi();

/*! Compute cholesky decomposition of diffusion tensor matrix. */
void compute_r_cholesky_mpi( DOUBLE* _D );

/*! Compute deterministic bead displacement. */
void compute_Rf_cholesky_mpi( DOUBLE** _Rf, DOUBLE* _D );

/*! Sum diffusion tensor matrix in IG-T algorithm. */
void sum_tensors_cholesky_mpi();

#ifdef	__cplusplus
}
#endif

#endif /* CHOLESKY_MPI_H */
