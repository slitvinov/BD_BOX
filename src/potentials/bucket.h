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

#ifndef BUCKET_H
#define	BUCKET_H

#define BUCKET_N 16
#define BUCKET_CAP_MLT 1
#define BUCKET_CAP_ADD 8

#define BUCKET_SPATIAL_BITS 10
#define BUCKET_SPATIAL_X  0
#define BUCKET_SPATIAL_Y 10
#define BUCKET_SPATIAL_Z 20
#define BUCKET_SPATIAL_M 1000

#define BUCKET_DIM 3

#include "../data.h"

#define BUCKET_NONE 0
#define BUCKET_BUCKET 1
#define BUCKET_SPATIAL 2

#define BUCKET_NONE_STR "brute"
#define BUCKET_BUCKET_STR "bucket_space"
#define BUCKET_SPATIAL_STR "bucket_list"

#ifdef	__cplusplus
extern "C" {
#endif

/*! Algorith to generate interaction list in spatial subdivision. */
extern BUCKET_NONE_SPATIAL nb_list;
MAKE_STR_DEC(BUCKET_NONE_SPATIAL,nb_list)
extern DOUBLE max_force_cutoff;

typedef union
{
    struct
    {
        DOUBLE x,y,z,w;
    } c;
    DOUBLE t[DIMS1];
} DOUBLE4;

typedef struct
{
    int x,y,z,size,capacity,*ids;
} reg_bucket;

typedef struct
{
    int size;
    reg_bucket* regs;
} list_reg;

typedef struct
{
    DOUBLE blen  [DIMS0];
    DOUBLE offset[DIMS0];
    INT mlocal[DIMS0];
    list_reg lreg[BUCKET_N][BUCKET_N][BUCKET_N];
} bucket_f;

typedef struct
{
    INT key, value;
} key_value;

typedef struct
{
    DOUBLE blen[3];
    INT size_keys;
    key_value* keys;
    INT* ids;
    INT size_ids;
    CHAR magic_tab[1<<BUCKET_DIM][1<<BUCKET_DIM];
    INT* coords;
    INT* phantoms;
    INT size_phan;
    INT outers;
} bucket_s;

extern INT magic_pos[][2];
extern INT magic_pos_size;


INT is_phantom( INT i );
INT can_gen_phantom(DOUBLE x,DOUBLE y,DOUBLE z);

bucket_f* init_bucket_f( int size, DOUBLE4* coords, DOUBLE* radii, DOUBLE blen );
void get_bid_f( bucket_f* br, DOUBLE4* bead, int* p, int* local );
void add_coord_f( bucket_f* br, int id, DOUBLE4* bead );
void free_bucket_f( bucket_f* br, int size, DOUBLE4* coord );
void print_statistic_f( bucket_f* br );
void get_neighbours( bucket_f* br, DOUBLE4* bead, int sizes[64], int* ids[64] );

bucket_s* init_keys( int size, DOUBLE4* coords, DOUBLE* radii, DOUBLE blen );
INT make_mag_pos( const bucket_s* bs, INT key, INT id );
INT can_interact( const bucket_s* bs, const key_value* i, const key_value* j);
void free_keys( bucket_s* keys );

#ifdef	__cplusplus
}
#endif

#endif/* BUCKET_H */
