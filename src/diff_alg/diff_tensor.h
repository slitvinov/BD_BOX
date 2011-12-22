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

#ifndef DIFF_TENSOR_H
#define	DIFF_TENSOR_H

#include "../data.h"

#define DIFF_ALG_ERMAK_CONST 0
#define DIFF_ALG_ERMAK_VAR 1
#define DIFF_ALG_ERMAK_NEWTON 2
#define DIFF_ALG_IGT_CONST 3
#define DIFF_ALG_IGT_VAR 4
#define DIFF_ALG_IGT_VAR_REV 5

#define DIFF_ALG_ERMAK_CONST_STR "ermak_const"
#define DIFF_ALG_ERMAK_VAR_STR "ermak_var"
#define DIFF_ALG_ERMAK_NEWTON_STR "ermak_newton"
#define DIFF_ALG_IGT_CONST_STR "igt_const"
#define DIFF_ALG_IGT_VAR_STR "igt_var"
#define DIFF_ALG_IGT_VAR_REV_STR "igt_var_rev"

#define DIFF_NONE 0
#define DIFF_CHOLESKY 1
#define DIFF_GEYER 2
#define DIFF_CHEBYSHEV 3

#define DIFF_NONE_STR "none"
#define DIFF_CHOLESKY_STR "cholesky"
#define DIFF_GEYER_STR "geyer"
#define DIFF_CHEBYSHEV_STR "chebyshev"

#define EWALD_METHOD_SMITH 0
#define EWALD_METHOD_BEENAKKER 1

#define EWALD_METHOD_SMITH_STR "smith"
#define EWALD_METHOD_BEENAKKER_STR "beenakker"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Diffusion tensor matrix. */
    extern DOUBLE* D;

/*! Second diffusion tensor matrix for IG-T algorithm. */
    extern DOUBLE* Ds;

/*! Diagonal elements of diffusion tensor matrix. */
    extern DOUBLE* diag_D;

/*! Inverted diagonal elements of diffusion tensor matrix. */
    extern DOUBLE* inv_diag_D;

/* Viscovity[P] */
extern DOUBLE visc;
MAKE_STR_DEC(DOUBLE,visc)

/*! Factor for viscosity units conversion. */
extern DOUBLE vfactor;
MAKE_STR_DEC(DOUBLE,vfactor)

/* Method used to evaluate hydrodynamic interactions. */
extern NO_CHOLS_GEYER hydro;
MAKE_STR_DEC(NO_CHOLS_GEYER,hydro)

/*! Algorithm used, IG-T with a constant/variable
    time step or the algorithm of Ermak, with a constant/variable time step or
    the algorithm of Ermak with a Newtonian correction for overlapping beads */
extern ERMAK_IGTCONST_IGTVAR algorithm;
MAKE_STR_DEC(ERMAK_IGTCONST_IGTVAR,algorithm)

/*! Number of attempts to draw a new random vector in case of overlap*/
extern INT move_attempts;
MAKE_STR_DEC(INT,move_attempts)

/*! The velocity gradient tensor matrix (row,column) */
extern PTR_DOUBLE9 vel_grad_tensor;
MAKE_STR_DEC(PTR_DOUBLE9,vel_grad_tensor)

/*! The coefficient of restitution for the ermak newton algorithm, falls between 0 and 1*/
extern DOUBLE e_collision;
MAKE_STR_DEC(DOUBLE,e_collision)

/*! Current dt. */
extern DOUBLE curr_dt;
MAKE_STR_DEC(DOUBLE,curr_dt)

/*! Ewald method for PBC. */
extern EWALD_METHOD ewald_method;
MAKE_STR_DEC(EWALD_METHOD,ewald_method)

/*! Switch on/off to use geyer on the fly, without diffusion tensor matrix. */
extern YESNO geyer_on_the_fly;
MAKE_STR_DEC(YESNO,geyer_on_the_fly)

/*! Flag indicates ewald summation. */
extern INT      ewald;

/*! Precomputed constant coefficient of diagonal element in diffusion tensor matrix. */
extern DOUBLE   preDii;

/*! Precomputed constant coefficient of off diagonal element in diffusion tensor matrix. */
extern DOUBLE   preDij;

/*! Array of random numbers. */
extern DOUBLE*  X;

/*! Arrays for deterministic displacements. */
extern DOUBLE** Rf;

/*! Arrays for random displacements. */
extern DOUBLE** Rs;

/*! Help array used in Geyer-Winter approach. */
extern DOUBLE*  epsilon;

/*! Help array used in homogeneous flow in IG-T algorithm. */
extern DOUBLE*  FS;

/*! Help array. */
extern DOUBLE*  Rf0;

/*! Saved decomposed diffusion tensor matrix. */
extern DOUBLE* S[2];

/*! Time stamp for S. */
extern INT save_iter;

/*! Current iteration. */
extern INT curr_iter;

/*! An integer which deffines the range of the real-space sum and controls its maximum
    number of vectors (i.e. image cells)*/
extern INT ewald_real;
MAKE_STR_DEC(INT,ewald_real)

/*! An integer deffining the summation range in the reciprocal-space and its number of vectors.*/
extern INT ewald_recip;
MAKE_STR_DEC(INT,ewald_recip)

/*! Parameter that controls the convergence of Ewald sums */
extern DOUBLE ewald_alpha;
MAKE_STR_DEC(DOUBLE,ewald_alpha)
    
/*! Initialization of diffusion tensor matrix. */
void init_diff_tensor();

/*! Initialization of diffusion tensor matrix using Ewald summation. */
void init_tensor();

/*! Free resources. */
void free_diff_tensor();

/*! Computes single cell of diffusion matrix tensor without periodic condition.
    \param i id of first bead
    \param j id of second bead
    \param tD nine values of Dij  */
void diff_tensor_calc_pos( INT i, INT j, DOUBLE* tD );

/*! Computes single cell of diffusion matrix tensor
    using Ewald summation for periodic condition proposed by Smith.
    \param i id of first bead
    \param j id of second bead
    \param _tD nine values of Dij  */
void diff_tensor_calc_ewald_smith( INT i, INT j, DOUBLE* _tD );

/*! Computes single diagonal cell of diffusion matrix tensor
    using Ewald summation for periodic condition proposed by Smith.
    \param i id of bead
    \param tD nine values of Dij  */
void diff_tensor_calc_ewald_smith_diag( INT i, DOUBLE* tD );

/*! Computes single cell of diffusion matrix tensor.
    using Ewald summation for periodic condition proposed by Beenakker.
    \param i id of first bead
    \param j id of second bead
    \param _tD nine values of Dij  */
void diff_tensor_calc_ewald_beenakker( INT i, INT j, DOUBLE* _tD );

/*! Computes single diagonal cell of diffusion matrix tensor
    using Ewald summation for periodic condition proposed by Beenakker.
    \param i id of bead
    \param tD nine values of Dij  */
void diff_tensor_calc_ewald_beenakker_diag( INT i, DOUBLE* tD );

/*! Dispath computation algorithm.
    \sa algorithm
    \param E nullable, array of energies */
void compute_hydro( DOUBLE* E );

/*! Generate Brownian dynamics trajectories using the original Ermak-McCammon algorithm.
    \param coord_out [out] new trajectory */
void compute_hydro_ermak( DOUBLE* coord_out );

/*! Computes correctorstep of IG-T algorithm. */
INT compute_hydro_IGT_second_phase( DOUBLE* _E );

/*! Computes diffusion tensor matrix.
 *  \param _D [out] matrix in column-wise for diffusion tensor*/
void compute_D( DOUBLE* _D );

/*! Computes single cell of diffusion matrix tensor using appropriate algorithm.
 *  \sa diff_tensor_calc_pos
 *  \sa diff_tensor_calc_ewald_smith
 *  \sa diff_tensor_calc_ewald_beenakker
 *  \param i id of first bead
 *  \param j id of second bead
    \param tD nine values of Dij */
void compute_Dij( INT i, INT j, DOUBLE* tD );

/*! Computes predictor step of IG-T algorithm. */
INT compute_hydro_IGT( DOUBLE* _E );

/*! Computes random displacement.
    \param _D diffusion matrix tensor used to compute corellation*/
void compute_r( DOUBLE* _D );

/*! Single diffusion matrix tensor computation.
 *  \sa compute_D*/
void compute_single_D( DOUBLE* _D );

/*! Save current coordinates. Used in step reverting after beads overlap. */
void save_curr_coord();

/*! Save current forces. Used in step reverting after beads overlap. */
void save_curr_forces();

/*! Load saved coordinates. Used in step reverting after beads overlap. */
void load_coord();

/*! Load saved forces. Used in step reverting after beads overlap. */
void load_forces();

#ifdef	__cplusplus
}
#endif

#endif	/* DIFF_TENSOR_H */
