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

#ifndef ERR_H
#define ERR_H

#include "data.h"

/*! Default name of error log file. */
#define ERR_FILENAME "err.log"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Name of plain text log file (output). */
extern STR out_filename;
MAKE_STR_DEC(STR,out_filename)

/*! Print error to log file and end program with error message.
    \param mssg Message to print
    \param file name of source file
    \param line number of line in source file*/
void unexpected_error( CSTR mssg, CSTR file, INT line );

#define UNERR(mssg) unexpected_error( mssg, __FILE__, __LINE__ );

/*! Print warning to log file.
    \param mssg Message to print
    \param file name of source file
    \param line number of line in source file */
void warning( CSTR mssg, CSTR file, INT line );

/*! Print log-message to log file.
    \param mssg Message to print
    \param file name of source file
    \param line number of line in source file*/
void logmssg( CSTR mssg, CSTR file, INT line );

/*! Check if memory was allocated, if not end program.
    \param ptr pointer to check
    \param file name of source file
    \param line number of line in source file*/
void checkMem( void* ptr, CSTR file, INT line );

#define CHMEM(ptr) checkMem( ptr, __FILE__, __LINE__ );

/*! Free resources, flush and close log file. */
void free_logfile();

#ifdef	__cplusplus
}
#endif

#endif	/* ERR_H */
