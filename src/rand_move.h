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

#ifndef RAND_MOVE_H
#define RAND_MOVE_H

#include "data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Initialization value for random generator. */
extern INT rand_seed;
MAKE_STR_DEC(INT,rand_seed)

/*! Dummy variable, used for polymorphic state of diffrent random generators. */
extern RAND_STATE rand_state;
MAKE_STR_DEC_INNER(RAND_STATE,rand_state)

/*! Generate variable from normal distribution of mean equal 0 and variance equal 1.
 *  \param tabl output array
 *  \param tabl_size number of random number to generate */
void gauss_rand(DOUBLE *tabl, INT tabl_size);

/*! Generate random number.
    \return single value from uniform distribution in [0,1) inteval. */
double rand_d48( );

/*! Initialization of random generator.
 *  \param seed initialization parameter. */
void init_rand( unsigned int seed );

/*! Free all resources used by random generator. */
void free_rand();

/*! Save state to argument.
    \param state [out] variable containing state of random generator.
    \param size [out] sizeof state in bytes */
void get_state( char** state, INT* size );

/*! Set state of random generator.
    \param state new state of random generator
    \param size size of state */
void set_state( char* state, INT size );

#ifdef	__cplusplus
}
#endif

#endif	/* RAND_MOVE_H */
