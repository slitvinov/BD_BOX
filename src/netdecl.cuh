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

#ifndef NETDECL_CUH
#define	NETDECL_CUH

#ifdef __NETBEANS__

/* Declarations for Netbeans ONLY!!! */

#define __CUDACC__ 1
#define __cplusplus 1
#include <vector_functions.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <math_functions.h>
#include <channel_descriptor.h>
#include <cuda_texture_types.h>
#include <texture_fetch_functions.h>

uint3 threadIdx;
uint3 blockIdx;
dim3 blockDim;
dim3 gridDim;
int warpSize;

#endif  /* __NETBEANS__ */

#endif	/* NETDECL_CUH */
