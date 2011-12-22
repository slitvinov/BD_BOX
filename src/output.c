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

#include "output.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "err.h"
#include "input.h"
#include "main_loop.h"
#include "restart.h"

#include "potentials/LJ.h"
#include "potentials/calc_func.h"

MAKE_STR_IN(STR,xyz_filename,NULL,"xyz–formatted trajectory file (output)")
MAKE_STR_IN(STR,dcd_filename,NULL,"dcd–formatted trajectory file (output)")
MAKE_STR_IN(STR,rst_filename,NULL,"name of the restart file")
MAKE_STR_IN(STR,enr_filename,NULL,"plain-text file with energies (output)")
MAKE_STR_IN(STR,pqr_filename,NULL,"pqr-formatted template file (output)")

MAKE_STR_IN(INT,save_xyz_freq,1,"frequency (number of steps) for writting to the xyz trajectory")
MAKE_STR_IN(INT,save_dcd_freq,1,"frequency (number of steps) for writting to the dcd trajectory")
MAKE_STR_IN(INT,save_rst_freq,1,"frequency (number of steps) for writting to the restart file")
MAKE_STR_IN(INT,save_enr_freq,1,"frequency (number of steps) for writting to the energy file")

FILE* xyz_file;
FILE* dcd_file;
FILE* enr_file;
FILE* pqr_file;
FILE* pqr_file_hi;

INT file_sizes[3];

static float* dcd_tab;

void save_xyz()
{
    int i;
    fprintf( xyz_file, "%d\n", size );
    fprintf( xyz_file, "%s time [ps] %" FORMAT_DOUBLEF "\n",xyz_filename, curr_time );
    for( i = 0; i < size; ++i )
    {
        fprintf( xyz_file, "%s %" FORMAT_DOUBLEF " %" FORMAT_DOUBLEF " %" FORMAT_DOUBLEF "\n", names[i],
                 coord[ i * DIMS1 ], coord[ i * DIMS1 + 1], coord[ i * DIMS1 + 2 ]);
    }
    fflush( xyz_file );
}

void save_dcd()
{
    int i;
    unsigned blocksize = size*sizeof(FLOAT);
    for( i = 0; i < size* DIMS1; ++i )
    {
        int pos = (i % DIMS1)*size + i / DIMS1;
        if( pos < 3 * size )
        {
            dcd_tab[ pos ] = (FLOAT)coord[i];
        }
    }
    for( i = 0; i < 3; ++i )
    {
        fwrite( &blocksize, sizeof(unsigned), 1, dcd_file );
        fwrite( dcd_tab + size * i, sizeof(FLOAT), size, dcd_file );
        fwrite( &blocksize, sizeof(unsigned), 1, dcd_file );
    }
    fflush( dcd_file );
}

void init_dcd( int numberOfFrames )
{
    unsigned int blockSize;
    char HDR[5];
    int ICNTRL[9];
    double DELTA;
    int ICNTRL2[10];
    char TITLE[80];
    int NTITL;
    int i;

    strcpy( HDR, "CORD" );
    ICNTRL[0] = numberOfFrames;
    ICNTRL[1] = 0;
    ICNTRL[2] = 1;
    for ( i = 3; i <= 8; i++ )
        ICNTRL[i] = 0;
    /*ICNTRL[9] = 1;*/
    DELTA = 1.0;
    for ( i = 0; i <= 9; i++ )
        ICNTRL2[i] = 0;

    NTITL = 1;
    memset( TITLE, 0, sizeof(char) * 80 );
    strcpy( TITLE, "bd_box DCD file" );

    blockSize = 84;

    /*doing header of dcd file*/
    fwrite( &blockSize, sizeof (unsigned int), 1, dcd_file );
    fwrite( HDR, sizeof (char), 4, dcd_file );
    fwrite( ICNTRL, sizeof (int), 9, dcd_file );
    fwrite( & DELTA, sizeof (float), 1, dcd_file );
    fwrite( ICNTRL2, sizeof (int), 10, dcd_file );
    fwrite( & blockSize, sizeof (int), 1, dcd_file );
    
    fwrite( & blockSize, sizeof (int), 1, dcd_file );
    fwrite( & NTITL, sizeof (int), 1, dcd_file );
    fwrite( TITLE, sizeof (char), 80, dcd_file );
    fwrite( & blockSize, sizeof (int), 1, dcd_file );
    blockSize = 4;
    fwrite( & blockSize, sizeof (int), 1, dcd_file );
    fwrite( & size, sizeof (int), 1, dcd_file );
    fwrite( & blockSize, sizeof (int), 1, dcd_file );

    fflush( dcd_file );
}

void save_energy()
{
    DOUBLE sum = E[ENERGY_BOND]+E[ENERGY_ANGLE]+E[ENERGY_ANGLE_COS]+E[ENERGY_DIHE_ANGLE]+E[ENERGY_DIHE]+E[ENERGY_LJ]+E[ENERGY_COULOMB];
    fprintf(enr_file, "%15" FORMAT_DOUBLEF "%15" FORMAT_DOUBLEF " %15" FORMAT_DOUBLEG " %15" FORMAT_DOUBLEG
        " %15" FORMAT_DOUBLEG " %15" FORMAT_DOUBLEG " %15" FORMAT_DOUBLEG " %15" FORMAT_DOUBLEG " %15" FORMAT_DOUBLEG "\n",
        curr_time,sum,E[ENERGY_BOND],E[ENERGY_ANGLE],E[ENERGY_ANGLE_COS],E[ENERGY_DIHE_ANGLE],E[ENERGY_DIHE],E[ENERGY_LJ],E[ENERGY_COULOMB]);
    fflush( enr_file );
}

void init_enr()
{
    fprintf(enr_file, "#%14s%15s %15s %15s %15s %15s %15s %15s %15s\n","TIME","ALL","BOND","ANG","ANG_COS","DIHE","DIHE_COS","LJ","COULOMB");
    fflush( enr_file );
}

void save_pqr()
{
    if( pqr_file )
    {
        int i;
        for( i = 0; i < size; ++i )
        {
            fprintf( pqr_file,
                "%-6s%5d %4s sub  %4d    %8.3" FORMAT_DOUBLEF "%8.3" FORMAT_DOUBLEF
                "%8.3" FORMAT_DOUBLEF "%6.2" FORMAT_DOUBLEF "%6.2" FORMAT_DOUBLEF "\n",
                "ATOM  ",ids[i],names[i],i+1,
                coord[DIMS1 * i], coord[DIMS1 * i + 1], coord[DIMS1 * i + 2],
                Q[i], LJ[ 2 * i  + LJ_SIGMA_OFF ]);
        }
        fclose(pqr_file);
        pqr_file = NULL;
    }
    if( pqr_file_hi )
    {
        int i;
        for( i = 0; i < size; ++i )
        {
            fprintf( pqr_file_hi,
                "%-6s%5d %4s sub  %4d    %8.3" FORMAT_DOUBLEF "%8.3" FORMAT_DOUBLEF
                "%8.3" FORMAT_DOUBLEF "%6.2" FORMAT_DOUBLEF "%6.2" FORMAT_DOUBLEF "\n",
                "ATOM  ",ids[i],names[i],i+1,
                coord[DIMS1 * i], coord[DIMS1 * i + 1], coord[DIMS1 * i + 2],
                Q[i], coord[DIMS1 * i + 3]);
        }
        fclose(pqr_file_hi);
        pqr_file_hi = NULL;
    }
}

void init_out( STR filename, FILE** out_file, const char* open_type )
{
    if( filename )
    {
        *out_file = fopen( filename, open_type );
        if( *out_file == NULL )
        {
            UNERR("Cann't open out-file");
        }
    }
    else
    {
        *out_file = NULL;
    }
}

void init_out_ar( STR filename, FILE** out_file, const char* open_type, int prev_size )
{
    if( !prev_size )/* in previous session wasn't defined*/
    {
        init_out( filename, out_file, open_type );
    }
    else
    {
        if( filename )
        {
            FILE* f = fopen( filename, "rb" );
            long int size_f;
            if( !f ) UNERR("Can't read file");
            /* get file size*/
            if( fseek( f, 0, SEEK_END ) ) UNERR("Can't read file");
            size_f = ftell( f );
            if( fseek( f, 0, SEEK_SET ) ) UNERR("Can't read file");
            if( size_f >= prev_size )
            {
                char* buffor = (char*)malloc( prev_size );
                if( fread( buffor, 1, prev_size, f ) != prev_size ) UNERR("Can't read file");
                fclose(f);
                *out_file = fopen( filename, open_type );
                if( !*out_file ) UNERR("Can't open file");
                if( fwrite( buffor, 1, prev_size, *out_file ) != prev_size ) UNERR("Can't write to file");
                fflush( *out_file );
                free(buffor);
            }
            else
            {
                fclose( f );
                warning( filename, __FILE__, __LINE__ );
                warning( "File too small, no output", __FILE__, __LINE__ );
            }
        }
        else
        {
            warning( "Filename was override by null", __FILE__, __LINE__ );
        }
    }
}

void save_output()
{
    if( xyz_file && (step % save_xyz_freq) == 0 ) save_xyz();
    if( dcd_file && (step % save_dcd_freq) == 0 ) save_dcd();
    if( enr_file && (step % save_enr_freq) == 0  ) save_energy();

    if( (step % save_rst_freq) == 0  ) restart_write( rst_filename );/* must be last */
}

void init_save()
{
    init_out(xyz_filename, &xyz_file, "w" );
    init_out(dcd_filename, &dcd_file, "wb" );
    if( dcd_file )
    {
        dcd_tab = (float*) malloc( sizeof(float) * size * DIMS0 ); CHMEM(dcd_tab);
        init_dcd( bdsteps / save_dcd_freq );
    }
    init_out(enr_filename, &enr_file, "w" );
    if( enr_file ) init_enr();
    init_out( pqr_filename, &pqr_file, "w" );/*pqr is closed after save*/
    if( pqr_filename )
    {
        char* buff = add_suffix( pqr_filename, PQR_SUFFIX);
        init_out( buff, &pqr_file_hi, "w" );
        free(buff);
    }
}

void init_save_after_restart()
{
    init_out_ar( xyz_filename, &xyz_file, "w", file_sizes[XYZ_FILE] );
    init_out_ar( dcd_filename, &dcd_file, "wb", file_sizes[DCD_FILE] );
    if( dcd_file )
    {
        int frames;
        dcd_tab = (float*) malloc( sizeof(float) * size * DIMS0 );CHMEM(dcd_tab)
        if( fseek(dcd_file,8,SEEK_SET) )
        {
            warning( "Can't seek in dcd file", __FILE__, __LINE__ );
        }
        else
        {
            frames = bdsteps / save_dcd_freq;/* because it may be new value */
            fwrite( &frames, sizeof(int), 1, dcd_file);
            if( fseek(dcd_file,0,SEEK_END) )  warning( "Can't seek in dcd file to end", __FILE__, __LINE__ );
        }
    }
    init_out_ar( enr_filename, &enr_file, "w", file_sizes[ENR_FILE] );
}

void free_save()
{
    if( dcd_file ) free( dcd_tab );
    if( dcd_file ) fclose( dcd_file );
    if( xyz_file ) fclose( xyz_file );
    if( enr_file ) fclose( enr_file );
}

char* add_suffix( const char* filename, const char* suffix )
{
    char* find;
    char* buff = (STR) malloc( sizeof(CHAR) * (strlen(filename) + strlen(suffix) + 2 ) );
    CHMEM(buff);
    memcpy( buff, filename, sizeof(CHAR) * (strlen(filename) + 1) );
    find = strrchr( buff, DOT );
    if( !find )
    {
        strcat( buff, suffix );
    }
    else
    {
        *find = '\0';
        strcat( buff, suffix );
        find[ strlen(suffix) ] = '\0';
        strcat( buff, filename + (find-buff) );
    }
    return buff;
}
