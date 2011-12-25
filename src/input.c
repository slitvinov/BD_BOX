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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "input.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cuda.h"
#include "err.h"
#include "output.h"
#include "rand_move.h"
#include "trans.h"

#include "diff_alg/cholesky_mpi.h"
#include "diff_alg/diff_tensor.h"
#include "potentials/angle.h"
#include "potentials/angle_cos.h"
#include "potentials/bond.h"
#include "potentials/bucket.h"
#include "potentials/calc_func.h"
#include "potentials/dihe.h"
#include "potentials/dihe_angle.h"
#include "potentials/electro.h"
#include "potentials/electro_ext.h"
#include "potentials/LJ.h"
#include "potentials/sboundary.h"
#include "potentials/pwell.h"

#define STR_HELP "help"
#define STR_HELPs "h"
#define STR_USAGE "usage"
#define STR_USAGEs "u"
#define STR_VER "version"
#define STR_VERs "v"
#define CHECK_OPT2( src, dst1, dst2 ) strncmp( src, "--" dst1, strlen("--" dst1) ) == 0 ||\
strncmp( src, "-" dst2, strlen("-" dst2) ) == 0

#ifdef _MSC_VER
#define STRDUP _strdup
#else
#define STRDUP strdup
#endif

/**************************************
 * Input format data declaration      *
***************************************/
STR prm_filename = NULL;

MAKE_STR_IN(DOUBLE,dt,0.0,"timestep")
MAKE_STR_IN(DOUBLE,T,298.15f,"temperature")
MAKE_STR_IN(INT,bdsteps,0,"number of simulation's steps")
MAKE_STR_IN(DOUBLE,xbox,0.0,"computational box, x-size")
MAKE_STR_IN(DOUBLE,ybox,0.0,"computational box, y-size")
MAKE_STR_IN(DOUBLE,zbox,0.0,"computational box, z-size")
MAKE_STR_IN(STR,str_filename,NULL,"structure file (input)")
MAKE_STR_IN(STR,restart,NULL,"apply restart from filename")
MAKE_STR_IN(STR,single_point,NULL,"Undocumented")

struct variable* input_prm[] = {
MAKE_STR_PTR(out_filename),
MAKE_STR_PTR(dt),
MAKE_STR_PTR(visc),
MAKE_STR_PTR(vfactor),
MAKE_STR_PTR(T),
MAKE_STR_PTR(epsilon_c),
MAKE_STR_PTR(kappa_c),
MAKE_STR_PTR(cutoff_c),
MAKE_STR_PTR(bdsteps),
MAKE_STR_PTR(gamma_c),
MAKE_STR_PTR(alpha_lj),
MAKE_STR_PTR(cutoff_lj),
MAKE_STR_PTR(lj_6_term),
MAKE_STR_PTR(save_xyz_freq),
MAKE_STR_PTR(save_dcd_freq),
MAKE_STR_PTR(save_rst_freq),
MAKE_STR_PTR(save_enr_freq),
MAKE_STR_PTR(xbox),
MAKE_STR_PTR(ybox),
MAKE_STR_PTR(zbox),
MAKE_STR_PTR(elec),
MAKE_STR_PTR(hydro),
MAKE_STR_PTR(algorithm),
MAKE_STR_PTR(str_filename),
MAKE_STR_PTR(xyz_filename),
MAKE_STR_PTR(dcd_filename),
MAKE_STR_PTR(enr_filename),
MAKE_STR_PTR(rst_filename),
MAKE_STR_PTR(pqr_filename),
MAKE_STR_PTR(rand_seed),
MAKE_STR_PTR(E_ext),
MAKE_STR_PTR(E_magn),
MAKE_STR_PTR(E_type),
MAKE_STR_PTR(E_freq),
MAKE_STR_PTR(E_dir1),
MAKE_STR_PTR(E_dir2),
MAKE_STR_PTR(E_factor),
MAKE_STR_PTR(bc),
MAKE_STR_PTR(sphere_radius),
MAKE_STR_PTR(bond_lj_scale),
MAKE_STR_PTR(bond_c_scale),
MAKE_STR_PTR(restart),
MAKE_STR_PTR(move_attempts),
MAKE_STR_PTR(vel_grad_tensor),
MAKE_STR_PTR(ewald_real),
MAKE_STR_PTR(ewald_recip),
MAKE_STR_PTR(ewald_alpha),
MAKE_STR_PTR(ewald_method),
MAKE_STR_PTR(MPI_nprow),
MAKE_STR_PTR(MPI_npcol),
MAKE_STR_PTR(MPI_block),
MAKE_STR_PTR(cuda_devices),
MAKE_STR_PTR(e_collision),
MAKE_STR_PTR(single_point),
MAKE_STR_PTR(nb_list),
MAKE_STR_PTR(sboundary),
MAKE_STR_PTR(sboundary_A),
MAKE_STR_PTR(sboundary_n),
MAKE_STR_PTR(sboundary_cutoff),
MAKE_STR_PTR(pwell),
MAKE_STR_PTR(pwell_A),
MAKE_STR_PTR(check_overlap),
MAKE_STR_PTR(geyer_on_the_fly),
MAKE_STR_PTR(cuda_block)};

INT size_input_prm = sizeof(input_prm)/sizeof(struct variable*);
STR* input_strs_args = NULL;/* on 'i' position str value for input_prm[i]->name from prm-file*/
STR* input_strs_file = NULL;/* on 'i' position str value for input_prm[i]->name from line arguments*/

/********************
 * Module functions *
*********************/
/*parsing*/
INT parse_args( int argc, char* argv[] )
{
    int i;
    if( argc == 1 )
    {
        IFP print_help();
        return 0;
    }
    else
    {
        int i;
        int files = 0;/*number of input files*/

        /*prepare place for arguments*/
        input_strs_args = (STR*) malloc( sizeof(STR) * size_input_prm );
        CHMEM(input_strs_args);
        memset( input_strs_args, 0, sizeof(STR) * size_input_prm );

        /*parse arguments*/
        for( i = 1; i < argc; ++i )
        {
            STR cmd = argv[i];
            if( CHECK_OPT2( cmd, STR_HELP, STR_HELPs) )
            {
                IFP print_help();
                return 0;
            }
            else if( CHECK_OPT2( cmd, STR_USAGE, STR_USAGEs) )
            {
                IFP print_usage();
                return 0;
            }
            else if( CHECK_OPT2( cmd, STR_VER, STR_VERs) )
            {
                IFP print_version();
                return 0;
            }
            else
            {
                int find_prm_args = 0;
                int j;
                size_t len_cmd = strlen(cmd);
                if( len_cmd > 3 && cmd[0] =='-'  && cmd[1] =='-' )
                {
                    for( j = 0; j < size_input_prm; ++j )
                    {
                        size_t name_len = strlen(input_prm[j]->name);
                        if( name_len + 3 < len_cmd  &&
                            strncmp( input_prm[j]->name, cmd + 2, name_len ) == 0 &&
                            cmd[ name_len + 2 ] == '=' )
                        {
                            STR tmp_str;
                            char buff[64];
                            find_prm_args = 1;
                            if( input_strs_args[j] != NULL )
                            {
                                sprintf( buff, "Arg option duplicated %s\n", input_prm[j]->name );
                                warning( buff, __FILE__, __LINE__ );
                            }
                            tmp_str = STRDUP( cmd + name_len + 3 );
                            input_strs_args[j] = STRDUP( cmd + name_len + 3 );
                            if( !check_type( input_prm[j]->type, tmp_str ) )
                            {
                                sprintf( buff, "Arg option has wrong value type %s\n", input_prm[j]->name );
                                UNERR(buff);
                            }
                            free( tmp_str );
                            break;
                        }
                    }
                }

                if( find_prm_args )
                {
                    continue;
                }
                else if( strlen( cmd ) && cmd[0] == '-' )
                {
                    IFP printf( "Unknown option: %s\n", cmd );
                    IFP print_usage();
                    return 0;
                }
                else
                {
                    if( files )
                    {
                        IFP printf("Too many files\n");
                        IFP print_usage();
                        return 0;
                    }
                    else
                    {
                        files++;
                        prm_filename = cmd;
                    }
                }
            }
        }
    }
    for( i = 0; i < size_input_prm; ++i )
    {
        char* value = input_strs_args[i];
        if( value )
        {
            apply_value( input_prm[i], value );
        }
    }
    return 1;
}

void read_prm_file()
{
    int len = 1024;
    int readed = -1;
    STR buff;
    int file_line = 0;
    FILE* prm_file;
    if( prm_filename == NULL ) UNERR( "prm_filename not defined");

    buff = (STR) malloc( sizeof(char) * len ); CHMEM(buff);
    input_strs_file = (STR*) malloc( sizeof(STR) * size_input_prm );
    CHMEM(input_strs_file);
    memset( input_strs_file, 0, sizeof(STR) * size_input_prm );
    
    prm_file = fopen( prm_filename, "r" );
    if( prm_file == NULL ) UNERR("Can't open prm-file");

    while( fgets( buff, len, prm_file ) )
    {
        if( readed ) /*Non empty line*/
        {
            int i;
            char* pos = strchr( buff, COMMENT_SIGN );
            char* keyword;
            char* value;
            char* end;
            file_line++;
            if( pos ) *pos = '\0'; /* 'remove' comment */
            if( !strlen(buff) ) continue;

            keyword = strtok( buff, " \t\n" );
            if( !keyword ) continue; /* ok we have only spaces and tabs */

            value = strtok( NULL, "\n" );
            if( !value || !value[0] ) UNERR("Wrong prm file format, no value");
            end = value + strlen(value)-1;
            while( *end == '\t' || *end == ' ' )
                *end-- = 0;

            for( i = 0; i < size_input_prm; ++i )
            {
                if( strcmp( input_prm[i]->name, keyword ) == 0 )
                {
                    if( input_strs_file[i] )
                    {
                        warning( "Arg in file duplicated", __FILE__, __LINE__ );
                    }

                    input_strs_file[i] = STRDUP( value );
                    
                    if( !check_type( input_prm[i]->type, value ) )
                    {
                        sprintf( buff, "File option has wrong value type %s\n", input_prm[i]->name );
                        UNERR(buff);
                    }

                    break;
                }
            }
            if( i >= size_input_prm )
            {
                char buff2[32];
                sscanf( keyword, "%31s", buff2 );
                sprintf( buff, "Wrong prm file format, unknown token, prm file line %d, token %s", file_line, buff2 );
                UNERR(buff);
            }
        }
    }
    free( buff );
    fclose( prm_file );
}

int simple_int_cmp(const void *a, const void *b)
{
    const int* aa = (const int*)a;
    const int* bb = (const int*)b;
    return *aa - *bb;
}

/*Apply options, check data*/
void apply_options( int res_mode )
{
    int i = 0;
    int j;
    int* cid;
    for( i = 0; i < size_input_prm; ++i )
    {
        char* value = input_strs_file ? input_strs_file[i] : NULL;
        value = input_strs_args[i] ? input_strs_args[i] : value;
        if( value )
        {
            /* from command line are setted in parse_args at the end,
             * but in restart mode their are overwrited in read of rst_file */
            if( res_mode || ( input_strs_file && value == input_strs_file[i] ) )
            {
                apply_value( input_prm[i], value );
            }
        }
        else if( !res_mode ) /* in case of restart we don't show warning */
        {
            char buff[128];
            sprintf( buff, "Parameter not set %s", input_prm[i]->name );
            warning( buff, __FILE__, __LINE__ );
        }
    }
    
    if( bc == BC_BOX )
    {
        if( xbox <= 0.0 || ybox <= 0.0 || zbox <= 0.0)
        {
            UNERR("Wrong box size");
        }
        inv_L[0] = 1 / xbox;
        inv_L[1] = 1 / ybox;
        inv_L[2] = 1 / zbox;
        inv_L[3] = 0.0;

        ewald = (hydro != DIFF_NONE);/* Run Ewald summation */
        if( ewald )
        {
            if( xbox!=ybox || ybox!=zbox ) UNERR("In Ewald box must be cube");

            for( i = 0; i < size; ++i )
            {
                if( LJ[ i*2 ] < coord[ i *DIMS1 + 3 ] )
                {
                    UNERR("LJ radii too small, should be at least 2 * bead radii");
                }
            }
            if( ewald_real != ewald_recip )
            {
                UNERR("ewald_real != ewald_recip");
            }
            if( ewald_real < 0 || ewald_real > 5 )
            {
                UNERR("ewald_real wrong value, only 0..5");
            }
        }

        box[0] = xbox;
        box[1] = ybox;
        box[2] = zbox;
        box[3] = 0.0;
        
        if( cutoff_c > xbox/2 ) warning("cutoff_c greater than half of xbox", __FILE__, __LINE__ );
        if( cutoff_c > ybox/2 ) warning("cutoff_c greater than half of ybox", __FILE__, __LINE__ );
        if( cutoff_c > zbox/2 ) warning("cutoff_c greater than half of zbox", __FILE__, __LINE__ );

        if( cutoff_lj > xbox/2 ) warning("cutoff_lj greater than half of xbox", __FILE__, __LINE__ );
        if( cutoff_lj > ybox/2 ) warning("cutoff_lj greater than half of ybox", __FILE__, __LINE__ );
        if( cutoff_lj > zbox/2 ) warning("cutoff_lj greater than half of zbox", __FILE__, __LINE__ );
        if( cutoff_lj < 0.0 && cutoff_lj != -1 ) UNERR("Wrong cutoff_lj value");
        if( cutoff_c < 0.0 ) UNERR("Wrong cutoff_c value");
#if USE_SSE
        UNERR("Currently, bc periodic is not implemented for SSE");
#endif
    }
    cid = (int*) malloc( sizeof(int)*size );
    if ( cid )
    {
        for( i = 0; i < size; ++i )
        {
            cid[i] = ids[i];
        }
        qsort( cid, size, sizeof(int), simple_int_cmp );
        for( i = 0; i < size-1; ++i )
        {
            if ( cid[i] == cid[i+1] ) UNERR( "Wrong atom numeration");
        }
        free( cid );
    }
    else
    {
        for( i = 0; i < size; ++i )
        {
            for( j = i + 1; j < size; ++j )
            {
                if( ids[i] == ids[j] ) UNERR("Wrong atom numeration");
            }
        }
    }
    for( i = 0; i < size; ++i )
    {
        if( coord[DIMS1*i +3] <= 0.0 ) UNERR( "Wrong atom radii");
    }

    if( vfactor == 0.0 ) UNERR("Wrong value for vfactor");
    if( visc == 0.0 ) UNERR( "Wrong value for visc");
    if ( alpha_lj < 0.0 ) UNERR( "Wrong value for alpha_lj");
    if( vel_grad_tensor && bc == BC_SPHERE )
    {
        warning( "Simulations in a flowing solvent turned off because of sphere bc", __FILE__, __LINE__ );
        free(vel_grad_tensor);
        vel_grad_tensor = NULL;
        MAKE_STR_SIZE(vel_grad_tensor,0);
    }
    if(lj_6_term) lj_6_term=1;/* 'y' - 'n' = 11 */
    if(E_ext)
    {
        if(E_dir1!='x'&&E_dir1!='y'&&E_dir1!='z')
            UNERR("Wrong value for E_dir1");
        if(E_type==2) if(E_dir2!='x'&&E_dir2!='y'&&E_dir2!='z')
            UNERR("Wrong value for E_dir2");
    }
#if USE_SSE
    if( ( 4 - 2*(sizeof(DOUBLE) == sizeof(double)) ) != SSECOUNT )
    {
        UNERR("SSE wrong data.h header");
    }
    if ( nb_list != BUCKET_NONE )
    {
        UNERR("nb_list not implemented in SSE mode, use brute");
    }
#endif
    if( elec && kappa_c < 0.0) UNERR("Wrong kappa_c value");
    if( algorithm == DIFF_ALG_ERMAK_NEWTON )
    {
        int i;
        for( i = 0; i < size; ++i )
        {
            if( masses[i] <= 0.0 ) UNERR("Wrong mass value");
            if( LJ[2*i] <= 0.0 ) UNERR("Wrong LJ value");
        }
    }
    if(hydro != DIFF_GEYER || algorithm >= DIFF_ALG_IGT_CONST )
    {
        geyer_on_the_fly = 0; /* force turn off */
        warning("Geyer on the fly turned off", __FILE__, __LINE__ );
    }
    if( sboundary_n < 0 ) UNERR("Wrong value for sboundary_n > 0");
}

void read_str_file()
{
    int len = 1024;
    int readed = -1;
    STR buff;
    int file_line = 0;
    FILE* str_file;
    if( str_filename == NULL ) UNERR("str_filename not defined");

    buff = (STR) malloc( sizeof(char) * len ); CHMEM(buff);

    str_file = fopen( str_filename, "r" );
    if( str_file == NULL ) UNERR("Can't open str-file");
    while( fgets( buff, len, str_file ) )
    {
        if( readed ) /*Non empty line*/
        {
            char* pos = strchr( buff, COMMENT_SIGN );
            char* keyword;
            file_line++;
            if( pos ) *pos = '\0'; /* 'remove' comment */
            if( !strlen(buff) ) continue;

            keyword = strtok( buff, " \t\n" );
            if( !keyword ) continue; /* ok we have only spaces and tabs */

            if( strcmp( keyword , SUB_TAG ) == 0 )
            {
                parse_sub();
            }
            else if( strcmp( keyword, BOND_TAG ) == 0 )
            {
                parse_bond();
            }
            else if( strcmp( keyword, ANGLE_TAG ) == 0 )
            {
                parse_angle();
            }
            else if( strcmp( keyword, DIHE_TAG ) == 0 )
            {
                parse_dihe();
            }
            else
            {
                warning( keyword, __FILE__, __LINE__ );
                UNERR( "Unknow tag in str file" );
            }
        }
    }
    end_sub();
    end_angle();
    end_angle_cos();
    end_dihe();
    end_dihe_angle();
    end_bond();

    free( buff );
    fclose( str_file );

    center_coords();/* Center coordinates */
}

void parse_sub()
{
    STR namei;
    INT idi;
    DOUBLE temp[8];
    namei = strtok( NULL, " \t\n" );
    if( namei == NULL ) UNERR("Wrong line in " SUB_TAG " no name");

    parse_ints( &idi, 1 );
    parse_doubles( temp, 8 );
    if( strtok( NULL, " \t\n" ) )
        UNERR("Wrong line in " SUB_TAG " too many columns");
    add_sub(namei,idi,temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],temp[7]);
}

void parse_bond()
{
    INT ids[2];
    DOUBLE temp[3];
    parse_ints( ids, 2 );
    parse_doubles( temp, 3 );
    if( strtok( NULL, " \t\n" ) )
        UNERR("Wrong line in " BOND_TAG " too many columns");
    add_bond( resolve_pos(ids[0]), resolve_pos(ids[1]),
              temp[0], temp[1], temp[2] );
}

void parse_angle()
{
    INT ids[3];
    DOUBLE temp[2];
    STR typea;
    typea = strtok( NULL, " \t\n" );
    if( typea == NULL ) UNERR("Wrong line in " ANGLE_TAG " no name");

    parse_ints( ids, 3 );
    parse_doubles( temp, 2 );
    if( strtok( NULL, " \t\n" ) ) UNERR("Wrong line in " ANGLE_TAG " too many columns");
    if( strcmp( typea, "cosine") == 0 )
    {
        add_angle_cos( resolve_pos(ids[0]), resolve_pos(ids[1]),
               resolve_pos(ids[2]), temp[0], temp[1] );
    }
    else
    {
        add_angle( resolve_pos(ids[0]), resolve_pos(ids[1]),
               resolve_pos(ids[2]), temp[0], temp[1] );
    }
}

void parse_dihe()
{
    INT ids[4];
    DOUBLE temp[3];
    STR typea;
    typea = strtok( NULL, " \t\n" );
    if( typea == NULL ) UNERR("Wrong line in " DIHE_TAG " no name");
    
    parse_ints( ids, 4 );
    parse_doubles( temp, 2 );
    if( strtok( NULL, " \t\n" ) ) UNERR("Wrong line in " DIHE_TAG " too many columns");
    if( strcmp( typea, "cosine" ) == 0 )
    {
        parse_doubles( temp+2, 1 );
        add_dihe( resolve_pos(ids[0]), resolve_pos(ids[1]),
              resolve_pos(ids[2]), resolve_pos(ids[3]),
              temp[0], temp[1], temp[2] );
    }
    else
    {
        add_dihe_angle( resolve_pos(ids[0]), resolve_pos(ids[1]),
              resolve_pos(ids[2]), resolve_pos(ids[3]),
              temp[1], temp[0] );
    }
}

void free_input()
{
    int i;
    if( input_strs_args )
    {
        for( i = 0; i < size_input_prm; ++i )
        {
            if( input_strs_args[i] )
            {
                free( input_strs_args[i] );
            }
        }
        free( input_strs_args );
    }
    if( input_strs_file )
    {
        for( i = 0; i < size_input_prm; ++i )
        {
            if( input_strs_file[i] )
            {
                free( input_strs_file[i] );
            }
        }
        free( input_strs_file );
    }
    for( i = 0; i < size_input_prm; ++i )
    {
        if( input_prm[i]->ptr )
        {
            if( strcmp( input_prm[i]->type, TO_STR(PTR_DOUBLE9) ) == 0 )
            {
                free(*(PTR_DOUBLE9*)input_prm[i]->ptr);
            }
            else if( strcmp( input_prm[i]->type, TO_STR(PTR_INT) ) == 0 )
            {
                free(*(PTR_INT*)input_prm[i]->ptr);
            }
            else if( strcmp( input_prm[i]->type, TO_STR(PTR_DOUBLE) ) == 0 )
            {
                free(*(PTR_DOUBLE*)input_prm[i]->ptr);
            }
            else if( strcmp( input_prm[i]->type, TO_STR(PTR_STR) ) == 0 )
            {
                free(*(PTR_STR*)input_prm[i]->ptr);
            }
            else if( strcmp( input_prm[i]->type, TO_STR(STR) ) == 0 && input_prm[i]->size == -1 )
            {
                free(*(STR*)input_prm[i]->ptr);
            }
        }
    }
    free_rand();
}

void print_type( struct variable* in )
{
    if( strcmp( in->type, TO_STR(DOUBLE) ) == 0 )
    {
        DOUBLE* val = (DOUBLE*)in->ptr;
        printf( "%" FORMAT_DOUBLEG "", *val );
    }
    else if( strcmp( in->type, TO_STR(PTR_INT) ) == 0 )
    {
        PTR_INT* val = (PTR_INT*)in->ptr;
        INT i;
        for ( i = 0; i < in->size; ++i )
        {
            printf( "%d", (*val)[i] );
        }
    }
	else if( strcmp( in->type, TO_STR(PTR_DOUBLE9) ) == 0 )
    {
        printf( "\"0 0 0 0 0 0 0 0 0\"" );
    }
    else if( strcmp( in->type, TO_STR(INT) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        printf( "%d", *val );
    }
    else if( strcmp( in->type, TO_STR(YESNO) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        if(*val)
        {
            printf( "yes" );
        }
        else
        {
            printf( "no" );
        }
    }
    else if( strcmp( in->type, TO_STR(BUCKET_NONE_SPATIAL) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        switch(*val)
        {
            case BUCKET_NONE:       printf(BUCKET_NONE_STR); break;
            case BUCKET_BUCKET:     printf(BUCKET_BUCKET_STR); break;
            case BUCKET_SPATIAL:    printf(BUCKET_SPATIAL_STR); break;
            default: UNERR("Unknow BUCKET_NONE_SPATIAL");
        }
    }
    else if( strcmp( in->type, TO_STR(CHAR) ) == 0 )
    {
        CHAR* val = (CHAR*)in->ptr;
        printf("%c",*val);
    }
    else if( strcmp( in->type, TO_STR(STR) ) == 0 )
    {
        STR* val = (STR*)in->ptr;
        printf("%s",*val);
    }
    else if( strcmp( in->type, TO_STR(NO_CHOLS_GEYER) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        switch(*val)
        {
            case DIFF_NONE:      printf(DIFF_NONE_STR); break;
            case DIFF_CHOLESKY:  printf(DIFF_CHOLESKY_STR); break;
            case DIFF_GEYER:     printf(DIFF_GEYER_STR); break;
            case DIFF_CHEBYSHEV: printf(DIFF_CHEBYSHEV_STR); break;
            default: UNERR("Unknow NO_CHOLS_GEYER");
        }
    }
    else if( strcmp( in->type, TO_STR(ERMAK_IGTCONST_IGTVAR) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        switch(*val)
        {
            case DIFF_ALG_ERMAK_CONST:    printf( DIFF_ALG_ERMAK_CONST_STR ); break;
            case DIFF_ALG_ERMAK_VAR:      printf( DIFF_ALG_ERMAK_VAR_STR ); break;
            case DIFF_ALG_ERMAK_NEWTON:   printf( DIFF_ALG_ERMAK_NEWTON_STR ); break;
            case DIFF_ALG_IGT_CONST:      printf( DIFF_ALG_IGT_CONST_STR   ); break;
            case DIFF_ALG_IGT_VAR:        printf( DIFF_ALG_IGT_VAR_STR     ); break;
            case DIFF_ALG_IGT_VAR_REV:    printf( DIFF_ALG_IGT_VAR_REV_STR ); break;
            default: UNERR("Unknow ERMAK_IGTCONST_IGTVAR");
        }
    }
    else if( strcmp( in->type, TO_STR(NONE_CUBIC_SPHERE) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        switch(*val)
        {
            case BC_NONE:   printf(BC_NONE_STR); break;
            case BC_BOX:    printf(BC_BOX_STR); break;
            case BC_SPHERE: printf(BC_SPHERE_STR); break;
            case BC_PWELL: printf(BC_PWELL_STR); break;
            default: UNERR("Unknow NONE_CUBIC_SPHERE");
        }
    }
    else if( strcmp( in->type, TO_STR(EWALD_METHOD) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        switch(*val)
        {
            case EWALD_METHOD_SMITH:        printf(EWALD_METHOD_SMITH_STR); break;
            case EWALD_METHOD_BEENAKKER:    printf(EWALD_METHOD_BEENAKKER_STR); break;
            default: UNERR("Unknow EWALD_METHOD");
        }
    }
    else if( strcmp( in->type, TO_STR(DC_AC_RF) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        switch(*val)
        {
            case 0: printf("DC"); break;
            case 1: printf("AC"); break;
            case 2: printf("RF"); break;
            default: UNERR("Unknow DC_AC_RF");
        }
    }
    else
    {
        UNERR("Unknown data type");
    }
}

int check_type( CSTR type, STR value )
{
    if( !value || value[0] == '\0' )
    {
        return 0;
    }
    if( strcmp( type, TO_STR(DOUBLE) ) == 0 )
    {
        DOUBLE val;
        return sscanf( value, "%" FORMAT_DOUBLESG "", &val );
    }
    else if( strcmp( type, TO_STR(INT) ) == 0 )
    {
        INT val;
        return sscanf( value, "%d", &val );
    }
    else if( strcmp( type, TO_STR(PTR_INT) ) == 0 )
    {
        CSTR cstr_len = strtok( value, " \t\n" );
        INT ilen;        
        if( sscanf( cstr_len, "%d", &ilen ) != 1 ) return 0;
        if( ilen < 0 ) return 0;
        while(ilen--)
        {
            INT tmp;
            CSTR str = strtok( NULL, " \t\n" );
            if( !str && sscanf( str, "%d", &tmp ) != 1 ) return 0;
        }
        if( strtok( NULL, "\n" ) ) return 0;
        return 1;
    }
    else if( strcmp( type, TO_STR(YESNO) ) == 0 )
    {
        return strcmp( value, "yes" )==0 || strcmp( value, "no" )==0;
    }
    else if( strcmp( type, TO_STR(NO_CHOLS_GEYER) ) == 0 )
    {
        int tmp = strcmp( value, DIFF_NONE_STR ) == 0 ||
                  strcmp( value, DIFF_CHOLESKY_STR ) == 0 ||
                  strcmp( value, DIFF_GEYER_STR ) == 0 ||
                  strcmp( value, DIFF_CHEBYSHEV_STR) == 0;
        return tmp;
    }
    else if( strcmp( type, TO_STR(ERMAK_IGTCONST_IGTVAR) ) == 0 )
    {
        int tmp = strcmp( value, DIFF_ALG_ERMAK_CONST_STR ) == 0 ||
                  strcmp( value, DIFF_ALG_ERMAK_VAR_STR ) == 0 ||
                  strcmp( value, DIFF_ALG_ERMAK_NEWTON_STR ) == 0 ||
                  strcmp( value, DIFF_ALG_IGT_CONST_STR   ) == 0 ||
                  strcmp( value, DIFF_ALG_IGT_VAR_STR     ) == 0 ||
                  strcmp( value, DIFF_ALG_IGT_VAR_REV_STR ) == 0;
        return tmp;
    }
    else if( strcmp( type, TO_STR(NONE_CUBIC_SPHERE) ) == 0 )
    {
        int tmp = strcmp( value, BC_NONE_STR ) == 0 ||
                  strcmp( value, BC_BOX_STR ) == 0 ||
                  strcmp( value, BC_SPHERE_STR ) == 0 ||
	          strcmp( value, BC_PWELL_STR ) == 0;
        return tmp;
    }
    else if( strcmp( type, TO_STR(BUCKET_NONE_SPATIAL) ) == 0 )
    {
        int tmp = strcmp( value, BUCKET_NONE_STR ) == 0 ||
                  strcmp( value, BUCKET_BUCKET_STR ) == 0 ||
                  strcmp( value, BUCKET_SPATIAL_STR ) == 0;
        return tmp;
    }
    else if( strcmp( type, TO_STR(EWALD_METHOD) ) == 0 )
    {
        int tmp = strcmp( value, EWALD_METHOD_SMITH_STR ) == 0 ||
                  strcmp( value, EWALD_METHOD_BEENAKKER_STR ) == 0;
        return tmp;
    }
    else if( strcmp( type, TO_STR(DC_AC_RF) ) == 0 )
    {
        int tmp = strcmp( value, "DC" ) == 0 ||
                  strcmp( value, "AC" ) == 0 ||
                  strcmp( value, "RF" ) == 0;
        return tmp;
    }
    else if( strcmp( type, TO_STR(CHAR) ) == 0 )
    {
        return value[0] && !value[1];/*Only one*/
    }
    else if( strcmp( type, TO_STR(STR) ) == 0 )
    {
        return 1;
    }
    else if( strcmp( type, TO_STR(PTR_DOUBLE9) ) == 0 )
    {
        DOUBLE val[9];
        return sscanf( value, "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG ""
            , val,val+1,val+2,val+3,val+4,val+5,val+6,val+7,val+8 ) == 9;
    }
    else
    {
        UNERR("Unknown data type");
        return 0;
    }
}

void apply_value( struct variable* in, STR value )
{
    if( strcmp( in->type, TO_STR(DOUBLE) ) == 0 )
    {
        DOUBLE* val = (DOUBLE*)in->ptr;
        *val = (FLOAT)atof(value);
    }
    else if( strcmp( in->type, TO_STR(INT) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        *val = atoi( value );
    }
    else if( strcmp( in->type, TO_STR(PTR_INT) ) == 0 )
    {
        PTR_INT* val = (PTR_INT*)in->ptr;
        STR cstr_len = strtok( value, " \t\n" );
        INT ilen, i;
        sscanf( cstr_len, "%d", &ilen );
        *val = (PTR_INT) malloc( sizeof(INT) * ilen );
        in->size = ilen;
        for ( i = 0; i < ilen; ++i )
        {
            CSTR cstr = strtok( NULL, " \t\n" );
            if( cstr )
            {
                sscanf( cstr, "%d", (*val) + i );
            }
            else
            {
                UNERR("Wrong token");
            }
        }
    }
    else if( strcmp( in->type, TO_STR(YESNO) ) == 0 )
    {
        INT* val = (INT*)in->ptr;
        *val = tolower(value[0]) - 'n';
    }
    else if( strcmp( in->type, TO_STR(CHAR) ) == 0 )
    {
        CHAR* val = (CHAR*)in->ptr;
        *val = value[0];
    }
    else if( strcmp( in->type, TO_STR(STR) ) == 0 )
    {
        STR* val = (STR*)in->ptr;
        if( in->size == -1 && (*val) ) free( *val );
        *val = value;
        in->size = 0;
    }
    else if( strcmp( in->type, TO_STR(DC_AC_RF) ) == 0 )
    {
        DC_AC_RF* val = (DC_AC_RF*)in->ptr;
        if( strcmp( value, "DC" ) == 0  )
        {
            *val = 0;
        }
        else if( strcmp( value, "AC" ) == 0  )
        {
            *val = 1;
        }
        else if( strcmp( value, "RF" ) == 0  )
        {
            *val = 2;
        }
        else
        {
            UNERR("Unknow value for hydro");
        }
    }
    else if( strcmp( in->type, TO_STR(EWALD_METHOD) ) == 0 )
    {
        EWALD_METHOD* val = (EWALD_METHOD*)in->ptr;
        if( strcmp( value, EWALD_METHOD_BEENAKKER_STR ) == 0  )
        {
            *val = EWALD_METHOD_BEENAKKER;
        }
        else if( strcmp( value, EWALD_METHOD_SMITH_STR ) == 0  )
        {
            *val = EWALD_METHOD_SMITH;
        }
        else
        {
            UNERR("Unknow value for ewald_method");
        }
    }
    else if( strcmp( in->type, TO_STR(NO_CHOLS_GEYER) ) == 0 )
    {
        NO_CHOLS_GEYER* val = (NO_CHOLS_GEYER*)in->ptr;
        if( strcmp( value, DIFF_NONE_STR ) == 0  )
        {
            *val = DIFF_NONE;
        }
        else if( strcmp( value, DIFF_CHOLESKY_STR ) == 0  )
        {
            *val = DIFF_CHOLESKY;
        }
        else if( strcmp( value, DIFF_GEYER_STR ) == 0  )
        {
            *val = DIFF_GEYER;
        }
        else if( strcmp( value, DIFF_CHEBYSHEV_STR ) == 0  )
        {
            *val = DIFF_CHEBYSHEV;
        }
        else
        {
            UNERR("Unknow value for hydro");
        }
    }
    else if( strcmp( in->type, TO_STR(ERMAK_IGTCONST_IGTVAR) ) == 0 )
    {
        ERMAK_IGTCONST_IGTVAR* val = (ERMAK_IGTCONST_IGTVAR*)in->ptr;
        if( strcmp( value, DIFF_ALG_ERMAK_CONST_STR ) == 0  )
        {
            *val = DIFF_ALG_ERMAK_CONST;
        }
        else if( strcmp( value, DIFF_ALG_ERMAK_VAR_STR ) == 0  )
        {
            *val = DIFF_ALG_ERMAK_VAR;
        }
        else if( strcmp( value, DIFF_ALG_ERMAK_NEWTON_STR ) == 0  )
        {
            *val = DIFF_ALG_ERMAK_NEWTON;
        }
        else if( strcmp( value, DIFF_ALG_IGT_CONST_STR ) == 0  )
        {
            *val = DIFF_ALG_IGT_CONST;
        }
        else if( strcmp( value, DIFF_ALG_IGT_VAR_STR ) == 0  )
        {
            *val = DIFF_ALG_IGT_VAR;
        }
        else if( strcmp( value, DIFF_ALG_IGT_VAR_REV_STR ) == 0  )
        {
            *val = DIFF_ALG_IGT_VAR_REV;
        }
        else
        {
            UNERR("Unknow value for hydro");
        }
    }
    else if( strcmp( in->type, TO_STR(NONE_CUBIC_SPHERE) ) == 0 )
    {
        NONE_CUBIC_SPHERE* val = (NONE_CUBIC_SPHERE*)in->ptr;
        if( strcmp( value, BC_NONE_STR ) == 0  )
        {
            *val = BC_NONE;
        }
        else if( strcmp( value, BC_BOX_STR ) == 0  )
        {
            *val = BC_BOX;
        }
        else if( strcmp( value, BC_SPHERE_STR ) == 0  )
        {
            *val = BC_SPHERE;
        }
        else if( strcmp( value, BC_PWELL_STR ) == 0  )
        {
            *val = BC_PWELL;
        }
        else
        {
            UNERR("Unknow value for bc");
        }
    }
    else if( strcmp( in->type, TO_STR(BUCKET_NONE_SPATIAL) ) == 0 )
    {
        BUCKET_NONE_SPATIAL* val = (BUCKET_NONE_SPATIAL*)in->ptr;
        if( strcmp( value, BUCKET_NONE_STR ) == 0  )
        {
            *val = BUCKET_NONE;
        }
        else if( strcmp( value, BUCKET_BUCKET_STR ) == 0  )
        {
            *val = BUCKET_BUCKET;
        }
        else if( strcmp( value, BUCKET_SPATIAL_STR ) == 0  )
        {
            *val = BUCKET_SPATIAL;
        }
        else
        {
            UNERR("Unknow value for bc");
        }
    }
    else if( strcmp( in->type, TO_STR(PTR_DOUBLE9) ) == 0 )
    {
        PTR_DOUBLE9* val = (PTR_DOUBLE9*)in->ptr;
        in->size = 9;
        *val = (DOUBLE*) malloc( sizeof( DOUBLE ) * in->size );
        sscanf( value, "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG "" \
            "%" FORMAT_DOUBLESG ""
            , *val,*val+1,*val+2,*val+3,*val+4,*val+5,*val+6,*val+7,*val+8 );
    }
    else
    {
        UNERR("Unknown data type");
    }
}

void parse_ints( INT* tab ,INT count )
{
    INT i;
    STR str_temp;
    for( i = 0; i < count; ++i )
    {
        str_temp = strtok( NULL, " \t" );
        if( str_temp == NULL || !check_type( TO_STR(INT), str_temp )  ) UNERR("Wrong line in sub");
        tab[ i ] = atoi( str_temp );
    }
}

void parse_doubles( DOUBLE* tab ,INT count )
{
    INT i;
    STR str_temp;
    for( i = 0; i < count; ++i )
    {
        str_temp = strtok( NULL, " \t" );
        if( str_temp == NULL || !check_type( TO_STR(DOUBLE), str_temp )  )  UNERR("Wrong line in sub");
        tab[ i ] = (DOUBLE)atof( str_temp );
    }
}

void print_help()
{
    int i;
    printf("bd_box help\n");
    for( i = 0; i < size_input_prm; ++i )
    {
        printf( "  --%s=", input_prm[i]->name );
        print_type(input_prm[i]);
        printf( ", %s\n", input_prm[i]->doc );
    }
    printf( "Report bugs to: " PACKAGE_BUGREPORT "\n" );
    printf( PACKAGE_NAME " home page: <" PACKAGE_URL ">\n" );
}

void print_version()
{
#ifdef PACKAGE
    printf(PACKAGE);
#endif
#ifdef PACKAGE_VERSION
    printf(" " PACKAGE_VERSION);
#endif
    printf("\n");
}

void print_usage()
{
    int i;
    printf("bd_box help\n");
    for( i = 1; i <= size_input_prm; ++i )
    {
        printf( "[ --%s=", input_prm[i-1]->name );
        print_type(input_prm[i-1]);
        printf(" ] ");
        if( (i % 5) == 0 ) printf("\n");
    }
    printf("\n");
}
