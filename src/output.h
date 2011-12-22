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

#ifndef _OUTPUT_H
#define	_OUTPUT_H

#include <stdio.h>

#include "data.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Save frequency of xyz trajectory file. */
extern INT save_xyz_freq;
MAKE_STR_DEC(INT,save_xyz_freq)

/*! Save frequency of dcd trajectory file. */
extern INT save_dcd_freq;
MAKE_STR_DEC(INT,save_dcd_freq)

/*! Save frequency of restart file, containing current system state. */
extern INT save_rst_freq;
MAKE_STR_DEC(INT,save_rst_freq)

/*! Save frequency of energy file. */
extern INT save_enr_freq;
MAKE_STR_DEC(INT,save_enr_freq)
        
/*! Filename of xyz trajectory file. */
extern STR xyz_filename;
MAKE_STR_DEC(STR,xyz_filename)

/*! Filename of dcd trajectory file. */
extern STR dcd_filename;
MAKE_STR_DEC(STR,dcd_filename)

/*! Filename of restart file, containing current system state. */
extern STR rst_filename;
MAKE_STR_DEC(STR,rst_filename)

/*! Filename of file containing energy of beads. */
extern STR enr_filename;
MAKE_STR_DEC(STR,enr_filename)

/*! Filename of pqr file. */
extern STR pqr_filename;
MAKE_STR_DEC(STR,pqr_filename)

/*! Index of xyz file size in file_sizes.
    \sa file_sizes */
#define XYZ_FILE 0

/*! Index of dcd file size in file_sizes.
    \sa file_sizes */
#define DCD_FILE 1

/*! Index of enr file size in file_sizes.
    \sa file_sizes */
#define ENR_FILE 2

/*! Suffix for pqr-file containing hydrodynamic radiis. */
#define PQR_SUFFIX "_hi"

/*! Dot */
#define DOT '.'

/*! Size of files xyz, dcd and enr. */
extern INT file_sizes[3];

/*! Save pqr file of system. */
void save_pqr();

/*! Save cordinates and energies to output files. */
void save_output();

/*! Init save files after restart.
    Trim file with file_sizes. */
void init_save_after_restart();

/*! Initialize files.*/
void init_save();

/*! Close and flush files. */
void free_save();

/*! xyz file descriptor. */
extern FILE* xyz_file;

/*! dcd file descriptor. */
extern FILE* dcd_file;

/*! Energy file descriptor. */
extern FILE* enr_file;

/*! pqr file descriptor. */
extern FILE* pqr_file;

/*! Help function to add suffix on end of copy of filename.
    \param filename to copy
    \param suffix string to add to filename, bu before extension*/
char* add_suffix( const char* filename, const char* suffix );

#ifdef	__cplusplus
}
#endif

#endif	/* OUTPUT_H */
