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

#ifndef CUDA_H
#define	CUDA_H

#include "data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Array of devices in OS to use. */
extern PTR_INT cuda_devices;
MAKE_STR_DEC(PTR_INT,cuda_devices)

/*! Size of block of threads in kernel. */
extern INT cuda_block;
MAKE_STR_DEC(INT,cuda_block)

#if USE_CUDA || USE_MAGMA

/*! Number of blocks in kernel. */
extern INT blocks;

/*! Number of streaming processors in device. */
extern INT sm_red;

/*! Array of coordinates on device. */
extern float* d_coord;

/*! Copy of coordinates on device. */
extern float* d_scoord;

/*! Copy of forces on device. */
extern float* d_sforces;

/*! Copy of forces on host. */
extern float* h_sforces;

/*! Array of forces on device. */
extern float* d_forces;

/*! Array of forces on host. */
extern float* h_forces;

/*! Diffusion tensor matrix on device.
 *  May contain only diagonal elements in non-hi case. */
extern float* d_D;

/*! Diffusion tensor matrix on host. */
extern float* h_D;

/*! Vector of random numbers on device. */
extern float* d_X;

/*! Vector of random numbers on host. */
extern float* h_X;

/*! Vector of beads displacements on device. */
extern float* d_Rs0;

/*! Vector of beads displacements on host. */
extern float* h_Rs0;

/*! Copy of beads displacements on device. */
extern float* d_sRs0;

/*! Sizeof diffusion tensor matrix. \sa d_D */
extern int    size_D;

/*! Size of coordinates array. \sa d_coord */
extern int    size_coord;

/*! Array of bonded forces values on device. */
extern float* d_bforces;

/*! Duplicate of coordinates for IGT algorithm on device. */
extern float* d_acoord;

/*! Duplicate of beads displacments for IGT algorithm on device. */
extern float* d_aRs0;

/*! Duplicate of beads displacments for IGT algorithm on device. */
extern float* d_aD;

/*! Initialize data for random number generation on device. */
void init_cuda_rand();

/*! Free resources on device. */
void free_cuda_rand();

/*! Initialize data for position computation on device. */
void init_cuda_equation();

/*! Free resources on device. */
void free_cuda_equation();

/*! Print on stdout float array from CUDA device.
    \param array array to print
    \param number of elements to print
    \param name optional name of array */
void show_array( float* array, int count, const char*name );
#endif

/*! Macro for checking errors from device. */
#define CCERR cudaCheckError( __FILE__, __LINE__ );

/*! Help function to check if cuda error ocurred */
void cudaCheckError( const char* file, int line );

/*! Help function to check if cudablas error ocurred */
void cublasCheckError( const char* file, int line );

/*! Initializes CUDA device and allocs memory after restart. */
void init_cuda_after_restart();

/*! Initializes CUDA device and allocs memory. */
void init_cuda();

/*! Free memory and other resources used in CUDA. */
void free_cuda();

#ifdef	__cplusplus
}
#endif

#endif	/* CUDA_H */
