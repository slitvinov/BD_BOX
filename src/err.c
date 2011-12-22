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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "err.h"

#include <stdio.h>
#include <stdlib.h>

#if USE_MPI
#include "mpi.h"
#endif

#include "output.h"

MAKE_STR_IN(STR,out_filename,ERR_FILENAME,"plain text log file (output)")

FILE* getLogFile();

void checkMem( void* ptr, CSTR file, INT line)
{
    if( ptr == NULL ) unexpected_error("Cann't alloc memory",file,line);
}

void unexpected_error( CSTR mssg, CSTR file, INT line )
{
    static INT lock = 0;
    FILE* log_file;
#if _OPENMP
#pragma omp critical
{
#endif
    if( lock ) /* because getLogFile can throw unexpected_error */
    {
        fprintf( stderr, "Double fault, terminating\n" );
        free_logfile();
#if USE_MPI
        MPI_Finalize();
#endif
        exit( 2 );
    }
    lock = 1;
    log_file = getLogFile();
    if( log_file )
    {
        fprintf( log_file, "ERROR: [%s:%d]\n%s\n", file, line, mssg );
        fprintf( log_file, "***********************************\n" );
        fprintf( stderr, "Unexpected error, check %s for more information\n", out_filename );
        free_logfile();
#if USE_MPI
        MPI_Finalize();
#endif
        exit( 1 );/*terminate*/
    }
    else
    {
        fprintf( stderr, "Can't open log file.\n" );
        free_logfile();
#if USE_MPI
        MPI_Finalize();
#endif
        exit( 3 );    
    }
#if _OPENMP
}
#endif
}

void warning( CSTR mssg, CSTR file, INT line )
{
    FILE* log_file = getLogFile();
    fprintf( log_file, "WARNING: [%s:%d]\n%s\n", file, line, mssg );
    fprintf( log_file, "***********************************\n" );
    fflush( log_file );
}

void logmssg( CSTR mssg, CSTR file, INT line )
{
    FILE* log_file = getLogFile();
    fprintf( log_file, "LOG: [%s:%d]\n%s\n", file, line, mssg );
    fprintf( log_file, "***********************************\n" );
    fflush( log_file );
}

static FILE* log_file = NULL;
FILE* getLogFile()
{
    if( log_file == NULL )
    {
        char* filename = out_filename;
        IFNP
        {
            char buff[64];
            sprintf( buff, "_%d", g_id );
            filename = add_suffix( out_filename, buff );    
        }
        log_file = fopen( filename, "w" );/* overwrite in all cases, even after restart */
        if( log_file == NULL )
        {
            log_file = stderr;
            warning( "Can't open log file", __FILE__, __LINE__ );
        }
        IFNP out_filename = filename;
    }
    return log_file;
}

void free_logfile()
{
    if( log_file )
    {
        fflush(log_file);
        fclose(log_file);
    }
    IFNP free( out_filename );
}
