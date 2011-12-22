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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "restart.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "err.h"
#include "input.h"
#include "main_loop.h"
#include "output.h"
#include "rand_move.h"

#include "potentials/angle.h"
#include "potentials/angle_cos.h"
#include "potentials/bond.h"
#include "potentials/dihe.h"
#include "potentials/dihe_angle.h"

#if CHAR_BIT != 8
#warning "Restart file might be incorrect"
#endif

#define HEADER_SIZE 8
#define FILES_SIZE 3
#define REVERSED 0x7F

#if USE_LAPACK
CHAR use_lapack = 1;
#else
CHAR use_lapack = 0;
#endif

#if USE_GSL
CHAR use_gsl = 1;
#else
CHAR use_gsl = 0;
#endif

#if HAVE_DRAND48
CHAR have_drand48 = 1;
#else
CHAR have_drand48 = 0;
#endif

CHAR size_double = sizeof(DOUBLE);
CHAR size_int = sizeof(INT);
CHAR endian = 0;/*little*/

/*! Check byte order.
 *  \return 1 if bigendian */
CHAR is_bigendian();

/*! Struct for inner variables, that should be saved. */
struct variable* input_rst[] = {
MAKE_STR_PTR(size),
MAKE_STR_PTR(curr_time),
MAKE_STR_PTR(step),
MAKE_STR_PTR(coord),
MAKE_STR_PTR(masses),
MAKE_STR_PTR(names),
MAKE_STR_PTR(ids),
MAKE_STR_PTR(Q),
MAKE_STR_PTR(LJ),
MAKE_STR_PTR(angle_cos_size),
MAKE_STR_PTR(angle_cos_pairs),
MAKE_STR_PTR(angle_cos_phis),
MAKE_STR_PTR(angle_size),
MAKE_STR_PTR(angle_pairs),
MAKE_STR_PTR(angle_phis),
MAKE_STR_PTR(bond_size),
MAKE_STR_PTR(bond_pairs),
MAKE_STR_PTR(bond_parms),
MAKE_STR_PTR(dihe_size),
MAKE_STR_PTR(dihe_pairs),
MAKE_STR_PTR(dihe_parms),
MAKE_STR_PTR(dihe_angle_size),
MAKE_STR_PTR(dihe_angle_pairs),
MAKE_STR_PTR(dihe_angle_parms),
MAKE_STR_PTR(rand_state)
};

INT size_input_rst = sizeof(input_rst)/sizeof(struct variable*);

/*! Read variable from file.
    \param var wrapper of variable in system
    \param f restart file*/
void read_value( struct variable* var, FILE* f );

/*! Read string from file.
    \param str [out] pointer of new string read from file
    \param f restart file */
void read_str( STR* str, FILE* f );

void restart_read( CSTR filename )
{
    CHAR header[HEADER_SIZE];
    INT i, j;
    STR var_name = NULL;
    FILE* file = fopen( filename, "rb");
    INT step1, step2;
    if( file == NULL ) UNERR( "Can't open restart file to read");

    endian = is_bigendian();
    /* Read header */
    if( fread( header, sizeof(char), HEADER_SIZE, file ) != HEADER_SIZE )
        UNERR("Wrong restart file header" );
    if( header[0] != use_lapack )
        warning( header[0] ? "LAPACK was used" : "Using LAPACK", __FILE__, __LINE__ );
    if( header[1] != use_gsl )
        warning( header[1] ? "Gsl was used" : "Using Gsl", __FILE__, __LINE__ );
    else if( !header[1] && header[2] != have_drand48 )
        warning( header[2] ? "drand48 was used" : "using drand48", __FILE__, __LINE__ );
    if( header[3] != size_double ) UNERR("double size differ");
    if( header[4] != size_int ) UNERR("int size differ");
    if( header[5] != endian ) UNERR("endian differ");
    if( header[6] != REVERSED ) UNERR("REVERSED value wrong");
    if( header[7] != REVERSED ) UNERR("REVERSED value wrong");

    if( fread( &step1, sizeof(INT), 1, file ) != 1 )
        UNERR("Wrong restart file header");

    if( fread( file_sizes, sizeof(INT), FILES_SIZE, file ) != FILES_SIZE )
        UNERR("Wrong restart file header");

    for( i = 0; i < size_input_prm + size_input_rst; ++i ) /* For all variables */
    {
        INT find = 0;
        read_str( &var_name, file );
        for( j = 0; j < size_input_prm; ++j )
        {
            if( strcmp( input_prm[j]->name, var_name ) == 0 )
            {
                find = 1;
                read_value( input_prm[j], file );
                break;
            }
        }
        if( !find )
        {
            for( j = 0; j < size_input_prm; ++j ) /* Check inner variables */
            {
                if( strcmp( input_rst[j]->name, var_name ) == 0 )
                {
                    find = 1;
                    read_value( input_rst[j], file );
                    break;
                }
            }
        }
        if( !find ) UNERR("Unknown tag in restart file");

        free( var_name ); var_name = NULL;
    }
    save_coord = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS1*size ); CHMEM(save_coord)
    save_forces = (DOUBLE*) malloc( sizeof(DOUBLE)*DIMS1*size ); CHMEM(save_forces)
    if( var_name ) free(var_name);
    step2 = -1;
    if( fread( &step2, sizeof(INT), 1, file ) != 1 || step1 != step2 || step2 != step )
    {
        printf( "%d %d %d\n", step, step1, step2 );
        UNERR("Restart file corrupted");
    }
    fclose( file );
    step++;/* snapshot was made at the end of loop */
}

/*! Reads variable from binary file. */
void read_value( struct variable* var, FILE* f)
{
    if( strcmp( var->type, TO_STR(DOUBLE) ) == 0 )
    {
        if( fread( var->ptr, sizeof(DOUBLE), 1, f ) != 1 ) UNERR("Can't read restart file");
    }
    else if( strcmp( var->type, TO_STR(INT) ) == 0 ||
             strcmp( var->type, TO_STR(YESNO) ) == 0 ||
             strcmp( var->type, TO_STR(CHAR) ) == 0 ||
             strcmp( var->type, TO_STR(NO_CHOLS_GEYER) ) == 0 ||
             strcmp( var->type, TO_STR(ERMAK_IGTCONST_IGTVAR) ) == 0 ||
             strcmp( var->type, TO_STR(EWALD_METHOD) ) == 0 ||
             strcmp( var->type, TO_STR(NONE_CUBIC_SPHERE) ) == 0 ||
             strcmp( var->type, TO_STR(BUCKET_NONE_SPATIAL) ) == 0 ||
             strcmp( var->type, TO_STR(DC_AC_RF) ) == 0 )
    {
        if( fread( var->ptr, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
    }
    else if( strcmp( var->type, TO_STR(STR) ) == 0 )
    {
        STR* val = (STR*)var->ptr;
        read_str( val, f );
        if( *val ) var->size = -1;/*to free flag*/
    }
    else if( strcmp( var->type, TO_STR(PTR_DOUBLE ) ) == 0 ||
             strcmp( var->type, TO_STR(PTR_DOUBLE9) ) == 0)
    {
        if( fread( &var->size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
        if( var->size > 0 )
        {
            *((PTR_DOUBLE*)var->ptr) = (DOUBLE*) malloc( sizeof(DOUBLE) * var->size ); checkMem(*((PTR_DOUBLE*)var->ptr),__FILE__,__LINE__);
            if( fread( *((PTR_DOUBLE*)var->ptr), sizeof(DOUBLE), var->size, f ) != var->size ) UNERR("Can't read to file");
        }
        else
        {
            *((PTR_DOUBLE*)var->ptr) = NULL;
        }
    }
    else if( strcmp( var->type, TO_STR(PTR_INT) ) == 0 )
    {
        if( fread( &var->size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
        if( var->size > 0 )
        {
            *((PTR_INT*)var->ptr) = (INT*) malloc( sizeof(INT) * var->size ); checkMem(*((PTR_INT*)var->ptr),__FILE__,__LINE__);
            if( fread( *((PTR_INT*)var->ptr), sizeof(INT), var->size, f ) != var->size ) UNERR("Can't read restart file");
        }
        else
        {
            *((PTR_INT*)var->ptr) = NULL;
        }
    }
    else if( strcmp( var->type, TO_STR(PTR_STR) ) == 0 )
    {
        PTR_STR strs = NULL;
        int i;
        if( fread( &var->size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
        if( var->size )
        {
            strs = (CHAR**) malloc( sizeof(PTR_STR) * var->size ); CHMEM(strs);
            for( i = 0; i < var->size; ++i )
            {
                read_str( strs + i, f );
            }
        }
        *((PTR_STR*)var->ptr) = strs;
    }
    else if( strcmp( var->type, TO_STR(RAND_STATE) ) == 0 )
    {
        INT size;
        char* buff;
        if( fread(&size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
        if( size )
        {
            buff = (char*) malloc( sizeof(char)*size ); CHMEM(buff);
            if( fread( buff, size, 1, f) != 1 ) UNERR("Can't read restart file");
            set_state( buff, size );
            free( buff );
        }
    }
    else
    {
        UNERR( "Unknown type to read");
    }
}

void read_str( STR* str, FILE* f)
{
    INT size;
    if( fread( &size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
    if( size > 0 )
    {
        *str = (STR) malloc( sizeof(CHAR) * (size + 2) ); CHMEM(str);
        if( fread( *str, sizeof(CHAR), size, f ) != size ) UNERR("Can't read restart file");
        (*str)[size] = '\0';
    }
    else
    {
        *str = NULL;
    }
}

static INT pos = 0;
static STR* filenames = NULL;
void check_filenames(CSTR);

/*! Write variable to file.
    \param var wrapper of variable in system
    \param f restart file*/
void write_value( struct variable* var, FILE* f );

/*! Write string to file.
    \param var string to write
    \param f restart file */
void write_str( CSTR var, FILE* f );

void restart_write( CSTR filename )
{
    CHAR header[HEADER_SIZE];
    INT sizes[FILES_SIZE];
    FILE* file = NULL;
    INT i;
    if( filename == NULL ) return;
    
    check_filenames( filename );
    file = fopen( filenames[pos], "wb"); if( !file ) UNERR("Can't open restart file to write");
    header[0] = use_lapack;
    header[1] = use_gsl;
    header[2] = have_drand48;
    header[3] = size_double;
    header[4] = size_int;
    header[5] = is_bigendian();
    header[6] = REVERSED;
    header[7] = REVERSED;

    if( fwrite( header, sizeof(CHAR), HEADER_SIZE, file ) != HEADER_SIZE ) UNERR("Can't write to file");

    if( fwrite( &step, sizeof(INT), 1, file ) != 1 ) UNERR("Can't write to file");

    memset( sizes, 0, sizeof(INT)*FILES_SIZE );

    sizes[XYZ_FILE] = xyz_file ? ftell(xyz_file) : 0;
    sizes[DCD_FILE] = dcd_file ? ftell(dcd_file) : 0;
    sizes[ENR_FILE] = enr_file ? ftell(enr_file) : 0;
    if( fwrite(sizes, sizeof(INT), FILES_SIZE, file ) != FILES_SIZE ) UNERR("Can't write to file");

    for( i = 0; i < size_input_prm; ++i )
    {
        write_value( input_prm[i], file );
    }

    for( i = 0; i < size_input_rst; ++i )
    {
        write_value( input_rst[i], file );
    }

    if( fwrite( &step, sizeof(INT), 1, file ) != 1 ) UNERR("Can't write to file");

    fflush( file );
    fclose( file );
}

void check_filenames( CSTR filename )
{
    if( filenames == NULL )
    {
        int i;
        int size_filename = (int)strlen(filename);
        filenames = (STR*) malloc( sizeof(STR) * 2 ); CHMEM(filenames);
        for( i = 0; i < 2; ++i )
        {
            filenames[i] = (STR) malloc( sizeof(CHAR)*(size_filename + 2) ); CHMEM(filenames[0]);
            strcpy( filenames[i], filename ); CHMEM(filenames[i]);
            filenames[i][ size_filename ] = (char)(i+'1');
            filenames[i][ size_filename + 1 ] = '\0';
        }
    }
    pos = pos ^ 1;
}

CHAR is_bigendian()
{
    union
    {
#ifdef u_int32_t
        u_int32_t a;
#else
        int a;
#endif
        CHAR b[4];
    } uint_char = {0x0100000F};
    const char* msg = (sizeof(uint_char.a) != 4 ||
        	      sizeof(uint_char.b) != 4) ?
		      "Wrong integer value": NULL;
    if ( msg ) UNERR(msg);
    return (CHAR)( uint_char.b[0] == 0x01 );
}

void write_str( CSTR str, FILE* f)
{
    if( str )
    {
        INT size_str = (INT)strlen( str );
        if( fwrite( &size_str, sizeof(INT), 1, f ) != 1 ) UNERR("Can't write to file");
        if( fwrite( str, sizeof(CHAR), size_str, f ) != size_str ) UNERR("Can't write to file");
    }
    else
    {
        INT size_str = 0; /* Zero size indicates NULL pointer */
        if( fwrite( &size_str, sizeof(INT), 1, f ) != 1 ) UNERR("Can't write to file");
    }
}

void write_value( struct variable* var, FILE* f)
{
    write_str( var->name, f );
    if( strcmp( var->type, TO_STR(DOUBLE) ) == 0 )
    {
        if( fwrite( var->ptr, sizeof(DOUBLE), 1, f ) != 1 ) UNERR("Can't write to file");
    }
    else if( strcmp( var->type, TO_STR(INT) ) == 0 ||
             strcmp( var->type, TO_STR(YESNO) ) == 0 ||
             strcmp( var->type, TO_STR(CHAR) ) == 0 ||
             strcmp( var->type, TO_STR(NO_CHOLS_GEYER) ) == 0 ||
             strcmp( var->type, TO_STR(ERMAK_IGTCONST_IGTVAR) ) == 0 ||
             strcmp( var->type, TO_STR(EWALD_METHOD) ) == 0 ||
             strcmp( var->type, TO_STR(NONE_CUBIC_SPHERE) ) == 0 ||
             strcmp( var->type, TO_STR(BUCKET_NONE_SPATIAL) ) == 0 ||
             strcmp( var->type, TO_STR(DC_AC_RF) ) == 0 )
    {
        if( fwrite( var->ptr, sizeof(INT), 1, f ) != 1 ) UNERR("Can't write to file");
    }
    else if( strcmp( var->type, TO_STR(STR) ) == 0 )
    {
        STR* val = (STR*)var->ptr;
        write_str( *val, f );
    }
    else if( strcmp( var->type, TO_STR(PTR_DOUBLE ) ) == 0 ||
             strcmp( var->type, TO_STR(PTR_DOUBLE9) ) == 0)
    {
        if( fwrite( &var->size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't write to file");
        if( var->size > 0 )
        {
            if( fwrite( *((PTR_DOUBLE*)var->ptr), sizeof(DOUBLE), var->size, f ) != var->size ) UNERR("Can't write to file");
        }
    }
    else if( strcmp( var->type, TO_STR(PTR_INT) ) == 0 )
    {
        if( fwrite( &var->size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't write to file");
        if( var->size > 0 )
        {
            if( fwrite( *((PTR_INT*)var->ptr), sizeof(INT), var->size, f ) != var->size ) UNERR("Can't write to file");
        }
    }
    else if( strcmp( var->type, TO_STR(PTR_STR) ) == 0 )
    {
        PTR_STR strs = *((PTR_STR*)var->ptr);
        int i;
        if( fwrite( &var->size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't write to file");
        for( i = 0; i < var->size; ++i )
        {
            write_str( strs[i], f );
        }
    }
    else if( strcmp( var->type, TO_STR(RAND_STATE) ) == 0 )
    {
        INT size;
        char* buff;
        get_state( &buff, &size );
        if( fwrite( &size, sizeof(INT), 1, f ) != 1 ) UNERR("Can't read restart file");
        if( fwrite( buff,  size,        1, f ) != 1 ) UNERR("Can't read restart file");
    }
    else
    {
        UNERR("Unknown type to write");
    }
}

void free_restart()
{
    if( filenames )
    {
        free( filenames[0] );
        free( filenames[1] );
        free( filenames );
        filenames = NULL;
    }
}
