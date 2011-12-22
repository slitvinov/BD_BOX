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

#ifndef INPUT_H
#define INPUT_H

#include "data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/**********************
 *    Module data     *
 *********************/

/*! Input *.prm filename. In this file the following syntax is used <br> <code>keyword</code> <i>value</i>. */
extern STR prm_filename;
MAKE_STR_DEC(STR,prm_filename)

/*! Input *.str filename. Contain beads and interactions description. */
extern STR str_filename;
MAKE_STR_DEC(STR,str_filename)

/*! Filename of restart file. */
extern STR restart;
MAKE_STR_DEC(STR,restart)

/*! Undocumented feature. */
extern STR single_point;
MAKE_STR_DEC(STR,single_point)

/*! Delta t, time, - ps */
extern DOUBLE dt;
MAKE_STR_DEC(DOUBLE,dt)

/*! Temperature */
extern DOUBLE T;
MAKE_STR_DEC(DOUBLE,T)

/*! Number of trajectory step. */
extern INT bdsteps;

/*! Computional (periodic box ) - x size. */
extern DOUBLE xbox;
MAKE_STR_DEC(DOUBLE,xbox)

/*! Computional (periodic box ) - y size. */
extern DOUBLE ybox;
MAKE_STR_DEC(DOUBLE,ybox)

/*! Computional (periodic box ) - z size. */
extern DOUBLE zbox;
MAKE_STR_DEC(DOUBLE,zbox)

/*! Global table of input wrappers of variables*/
extern struct variable* input_prm[];

/*! Size of input_prm.
    \sa input_prm */
extern int size_input_prm;/*size of input_prm*/

/*! Command line arguments for variables.
    i-th position contain string value for i-th position in input_prm.
    \sa input_prm */
extern STR* input_strs_args;

/*! PRM input file arguments for variables.
    i-th position contain string value for i-th position in input_prm.
    \sa input_prm */
extern STR* input_strs_file;

/********************
 * Module functions *
*********************/

/*! Parsing program options from arguments */
INT parse_args( int argc, char* argv[] );

/*! Read prm option file */
void read_prm_file();

/*! Apply options.
    \param res_mode seted value indicates restart mode */
void apply_options( int res_mode );

/*! Read beads structure file.
    \sa parse_sub
    \sa parse_bond
    \sa parse_angle
    \sa parse_dihe*/
void read_str_file();

/*! Parse bead description from global variable initialized by strtok*/
void parse_sub();

/*! Parse bond description from global variable initialized by strtok*/
void parse_bond();

/*! Parse angle description from global variable initialized by strtok*/
void parse_angle();

/*! Parse dihe-angle description from global variable initialized by strtok*/
void parse_dihe();

/*! Free resources */
void free_input();

/*! Help function to parse ints, from global variable initialized by strtok */
void parse_ints( INT* tab ,INT count );

/*! Help function to parse double, from global variable initialized by strtok */
void parse_doubles( DOUBLE* tab ,INT count );

/* Type conversion */

/*! Print single input value, depends on it's type*/
void print_type( struct variable* in );

/*! Check if value is of type
    \param type INT, DOUBLE, YESNO application defined types
    \param value string representing value of type
    \return true if can convert string value to type */
int check_type( CSTR type, STR value );

/*! Apply value from desciption to wrapper structure containing value.
    \param in wrapper containing destination variable
    \param value string value to be set to variable */
void apply_value( struct variable* in, STR value );

/*! Prints help on stdout */
void print_help();

/*! Prints usage on stdout */
void print_usage();

/*! Prints version on stdout */
void print_version();

#ifdef	__cplusplus
}
#endif

#endif	/* INPUT_H */
