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

#ifndef DIFF_TENSOR_CUH
#define DIFF_TENSOR_CUH

/*! Computes deterministic bead displacement in fly */
void comp_rinfly( float* coords, float* forces, float* r );

/*! Computes random displacement using Geyer-Winter approach in fly. First phase. */
void comp_geyer_beg_in_fly( float* coord, float* ptr_x, float* ptr_eff,
                                       float* ptr_c, float* ptr_eps );

/*! Computes random displacement using Geyer-Winter approach in fly. Second phase. */
void comp_geyer_end_in_fly( float* coord, float* ptr_x, float* ptr_eff,
                                       float* ptr_c, float beta, float* out);

#endif /* DIFF_TENSOR_CUH */
