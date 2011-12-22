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

#include "bucket.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "../input.h"
#include "../err.h"
#include "../math_help.h"
#include "../trans.h"
#include "LJ.h"
#include "electro.h"
#include "calc_func.h"

MAKE_STR_IN(BUCKET_NONE_SPATIAL,nb_list,0,"generation of nonbonded interactions lists algorithm")
DOUBLE max_force_cutoff;

INT in_phan_reg( DOUBLE x, INT dim )
{
    switch(dim)
    {
        case 0: return x <= -xbox/2 + max_force_cutoff + 0.01f;
        case 1: return x <= -ybox/2 + max_force_cutoff + 0.01f;
        case 2: return x <= -zbox/2 + max_force_cutoff + 0.01f;
        default: UNERR("Wrong dim"); return 0;
    }
}

INT in_phan_reg_up( DOUBLE x, INT dim )
{
    switch(dim)
    {
        case 0: return x >= xbox/2 - max_force_cutoff - 0.01f;
        case 1: return x >= ybox/2 - max_force_cutoff - 0.01f;
        case 2: return x >= zbox/2 - max_force_cutoff - 0.01f;
        default: UNERR("Wrong dim"); return 0;
    }
}

INT can_gen_phantom(DOUBLE x,DOUBLE y,DOUBLE z)
{
    if ( bc == BC_BOX && ( in_phan_reg(x,0) || in_phan_reg(y,1) || in_phan_reg(z,2) ))
        return 1;
    return 0;
}

INT is_phantom( INT i )
{
    return i >= size;
}

void get_h_l( INT size, DOUBLE4* coords, DOUBLE* radii, DOUBLE blen, DOUBLE h[3], DOUBLE l[3] )
{
    INT i, j;
    for(j=0;j<3;++j) h[j] = -FLT_MAX;
    for(j=0;j<3;++j) l[j] = FLT_MAX;
    for( i = 0; i < size; ++i )
    {
        for(j=0;j<3;++j)
        {
            if( coords[i].t[j] > h[j] )
            {
                h[j]=coords[i].t[j];
            }
            if( coords[i].t[j] < l[j] )
            {
                l[j]=coords[i].t[j];
            }
        }
        if( radii && 2*radii[2*i] > blen ) blen = 2*radii[2*i] + 0.01f;
    }
}

INT correct_indx( int indx, int last )
{
    if( indx < 0 ) return last + indx;
    return (indx)%last;
}

int cmp_keys( const void* a, const void* b )
{
    const key_value* aa = (const key_value*)a;
    const key_value* bb = (const key_value*)b;
    if( aa->key - bb->key == 0 )
        return aa->value - bb->value;
    else
        return aa->key - bb->key;
}

INT can_interact( const bucket_s* bs, const key_value* i, const key_value* j)
{
    INT m1 = make_mag_pos( bs, i->key, i->value );
    INT m2 = make_mag_pos( bs, j->key, j->value );
    if ( bc == BC_BOX )
    {
        if ( is_phantom( i->value ) &&
             is_phantom( j->value ) ) return 0;
    }
    if ( i->value >= bs->outers && j->value >= bs->phantoms[i->value - size] ) return 0;
    if ( j->value >= bs->outers && i->value >= bs->phantoms[j->value - size] ) return 0;
    if ( m1 == -1 || m2 == -1 )
        UNERR( "Wrong cell");
    return bs->magic_tab[m1][m2];
}

INT make_mag_pos( const bucket_s* bs, INT key, INT id )
{
    INT x,y,z, kx,ky,kz;
    INT mask;
    /* primary cell */
    x = bs->coords[DIMS0*id  ];
    y = bs->coords[DIMS0*id+1];
    z = bs->coords[DIMS0*id+2];
    /* current cell */
    mask = (1 << BUCKET_SPATIAL_BITS)-1;
    kx = (key & (mask << BUCKET_SPATIAL_X))>>BUCKET_SPATIAL_X;
    ky = (key & (mask << BUCKET_SPATIAL_Y))>>BUCKET_SPATIAL_Y;
    kz = (key & (mask << BUCKET_SPATIAL_Z))>>BUCKET_SPATIAL_Z;
    if( ((kx-x)!=0 && (kx-x)!=1) || ((ky-y)!=0 && (ky-y)!=1) || ((kz-z)!=0 && (kz-z)!=1) ) return -1;
    return ((kx-x)<<2) | ((ky-y)<<1) | (kz-z);
}

INT magic_pos[][2] = { {0,0}, {0,1}, {0,2}, {0,3}, {0,4}, {0,5}, {0,6},
                       {0,7}, {1,2}, {1,4}, {1,6}, {2,4}, {2,5}, {3,4} };
INT magic_pos_size = sizeof(magic_pos)/sizeof(magic_pos[0]);

INT get_coord( INT key, INT p )
{
    INT mask = (1 << BUCKET_SPATIAL_BITS)-1;
    if ( p == 0 ) return (key & (mask << BUCKET_SPATIAL_X))>>BUCKET_SPATIAL_X;
    if ( p == 1 ) return (key & (mask << BUCKET_SPATIAL_Y))>>BUCKET_SPATIAL_Y;
    if ( p == 2 ) return (key & (mask << BUCKET_SPATIAL_Z))>>BUCKET_SPATIAL_Z;
    UNERR("Unknown dimension");
    return -1;
}

INT make_keyb(INT x,INT y,INT z)
{
    INT cx = correct_indx( x, BUCKET_SPATIAL_M + 1 );
    INT cy = correct_indx( y, BUCKET_SPATIAL_M + 1 );
    INT cz = correct_indx( z, BUCKET_SPATIAL_M + 1 );
    INT key = ( cx << BUCKET_SPATIAL_X) |
              ( cy << BUCKET_SPATIAL_Y) |
              ( cz << BUCKET_SPATIAL_Z);
    return key;
}

bucket_s* init_keys( int size, DOUBLE4* coords, DOUBLE* radii, DOUBLE blen )
{
    DOUBLE h[3], l[3];
    INT i;
    INT prev;
    INT keys_size = 0;
    INT phan_count = 0;
    bucket_s* ret = (bucket_s*) malloc( sizeof(bucket_s) ); CHMEM(ret);

    get_h_l( size, coords, radii, blen, h, l );
    if ( bc == BC_BOX )
    {
        h[0] = MAXD( max_force_cutoff + xbox/2, h[0] );
        h[1] = MAXD( max_force_cutoff + ybox/2, h[1] );
        h[2] = MAXD( max_force_cutoff + zbox/2, h[2] );
        l[0] = MIND( xbox/2 - 3*max_force_cutoff, l[0]);
        l[1] = MIND( ybox/2 - 3*max_force_cutoff, l[1]);
        l[2] = MIND( zbox/2 - 3*max_force_cutoff, l[2]);
    }
    ret->blen[0] = MAXD( (h[0]-l[0])/BUCKET_SPATIAL_M, blen);
    ret->blen[1] = MAXD( (h[1]-l[1])/BUCKET_SPATIAL_M, blen);
    ret->blen[2] = MAXD( (h[2]-l[2])/BUCKET_SPATIAL_M, blen);
    if ( bc == BC_BOX )
    {
        if ( xbox / ret->blen[0] <= 3.01 ||
             ybox / ret->blen[1] <= 3.01 ||
             zbox / ret->blen[2] <= 3.01 ||
             l[0] < -xbox/2 || l[1] < -ybox/2 || l[2] < -zbox/2 )
        {
            free(ret);
            return NULL;
        }
    }
    memset( ret->magic_tab, 0, sizeof( ret->magic_tab ) );
    for ( i = 0; i < sizeof(magic_pos)/sizeof(magic_pos[0]); ++i )
    {
        ret->magic_tab[magic_pos[i][0]][magic_pos[i][1]] = 1;
        ret->magic_tab[magic_pos[i][1]][magic_pos[i][0]] = 1;
    }
    ret->phantoms = NULL;
    ret->outers = size;
    if ( bc == BC_BOX )
    {
        INT ctab[] = {0,4,10,19};
        for ( i = 0; i < size; ++i )
        {
            DOUBLE x,y,z;
            INT c,d;
            x = coords[i].c.x;
            y = coords[i].c.y;
            z = coords[i].c.z;
            c =    in_phan_reg(x,0) +    in_phan_reg(y,1) +    in_phan_reg(z,2);
            d = in_phan_reg_up(x,0) + in_phan_reg_up(y,1) + in_phan_reg_up(z,2);
            phan_count += ctab[c]*( 1 + d );
        }
        ret->phantoms = (INT*)malloc( sizeof(INT) * phan_count );
        ret->size_phan = phan_count;
    }
    ret->keys = (key_value*) malloc( sizeof(key_value)*(size*8 + phan_count+1) ); CHMEM(ret->keys);
    ret->size_keys = size*8 + phan_count;
    ret->coords = (INT*) malloc( sizeof(INT)*(size+phan_count+1)*DIMS0 ); CHMEM(ret->coords);
    for ( i = 0; i < size; ++i )
    {
        INT x,y,z,k;

        x = (INT)FLOORD((coords[i].c.x-l[0])/ret->blen[0]);
        y = (INT)FLOORD((coords[i].c.y-l[1])/ret->blen[1]);
        z = (INT)FLOORD((coords[i].c.z-l[2])/ret->blen[2]);
        ret->coords[i*DIMS0  ]=x;
        ret->coords[i*DIMS0+1]=y;
        ret->coords[i*DIMS0+2]=z;
        
        for ( k = 0; k < 8; ++k )
        {
            ret->keys[ keys_size ].key = make_keyb( x+(k&1), y+((k&2)>>1), z+((k&4)>>2) );
            ret->keys[ keys_size ].value = i;
            keys_size++;
        }
    }
    if ( bc == BC_BOX )
    {
        INT id = size;
        for ( i = 0; i < size; ++i )
        {
            INT k, n, m;
            const DOUBLE xx = coords[i].c.x;
            const DOUBLE yy = coords[i].c.y;
            const DOUBLE zz = coords[i].c.z;
            for ( k = 0; k <= in_phan_reg(xx,0); ++k )
                for ( n = 0; n <= in_phan_reg(yy,1); ++n )
                    for ( m = 0; m <= in_phan_reg(zz,2); ++m )
                    {
                        if( k || n || m )
                        {
                            INT ix,iy,iz;
                            INT a,b,c;
                            DOUBLE dx,dy,dz;
                            INT added = 0;
                            dx = xx + k * xbox;
                            dy = yy + n * ybox;
                            dz = zz + m * zbox;

                            ix = (INT)FLOORD((dx-l[0])/ret->blen[0]);
                            iy = (INT)FLOORD((dy-l[1])/ret->blen[1]);
                            iz = (INT)FLOORD((dz-l[2])/ret->blen[2]);
                            ret->coords[id*DIMS0]=ix;
                            ret->coords[id*DIMS0+1]=iy;
                            ret->coords[id*DIMS0+2]=iz;

                            for ( a = 0; a <= 1-k; ++a )
                                for ( b = 0; b <= 1-n; ++b )
                                    for ( c = 0; c <= 1-m; ++c )
                                    {
                                        ret->keys[ keys_size ].key = make_keyb( ix+a, iy+b, iz+c );
                                        ret->keys[ keys_size ].value = id;
                                        keys_size++;
                                        added = 1;
                                    }
                            ret->phantoms[id-size] = i;
                            id++;
                        }
                    }
        }
        ret->outers = id;
        for ( i = 0; i < size; ++i )
        {
            INT k, n, m;
            const DOUBLE xx = coords[i].c.x;
            const DOUBLE yy = coords[i].c.y;
            const DOUBLE zz = coords[i].c.z;
            for ( k = 0; k <= in_phan_reg(xx,0); ++k )
                for ( n = 0; n <= in_phan_reg(yy,1); ++n )
                    for ( m = 0; m <= in_phan_reg(zz,2); ++m )
                    {
                        if( k || n || m )
                        {
                            INT ix,iy,iz;
                            INT d,e,f;
                            DOUBLE dx,dy,dz;
                            dx = xx + k * xbox;
                            dy = yy + n * ybox;
                            dz = zz + m * zbox;

                            ix = (INT)FLOORD((dx-l[0])/ret->blen[0]);
                            iy = (INT)FLOORD((dy-l[1])/ret->blen[1]);
                            iz = (INT)FLOORD((dz-l[2])/ret->blen[2]);

                            for ( d = 0; d <= in_phan_reg_up(xx,0); ++d )
                                for ( e = 0; e <= in_phan_reg_up(yy,0); ++e )
                                    for ( f = 0; f <= in_phan_reg_up(zz,0); ++f )
                                    {
                                        if( d || e || f )
                                        {
                                            INT added = 0;
                                            INT cx, cy, cz;
                                            cx = ret->coords[id*DIMS0  ] = (INT)FLOORD((dx-d*xbox-l[0])/ret->blen[0]);
                                            cy = ret->coords[id*DIMS0+1] = (INT)FLOORD((dy-e*ybox-l[1])/ret->blen[1]);
                                            cz = ret->coords[id*DIMS0+2] = (INT)FLOORD((dz-f*zbox-l[2])/ret->blen[2]);
                                            if ( cx>=-1 && cy >=-1 && cz>=-1)
                                            {
												INT a, b, c;
                                                cx = d ? 0 : ix;
                                                cy = e ? 0 : iy;
                                                cz = f ? 0 : iz;

                                                for ( a = 0; a <= (1-k)*(1-d); ++a )
                                                    for ( b = 0; b <= (1-n)*(1-e); ++b )
                                                        for ( c = 0; c <= (1-m)*(1-f); ++c )
                                                        {
                                                            INT X = cx+a;
                                                            INT Y = cy+b;
                                                            INT Z = cz+c;
                                                            ret->keys[ keys_size ].key = make_keyb(X,Y,Z);
                                                            ret->keys[ keys_size ].value = id;
                                                            keys_size++;
                                                            added = 1;
                                                        }
                                                if ( added )
                                                {
                                                    ret->phantoms[id-size] = i;
                                                    id++;
                                                }
                                            }
                                        }
                                    }
                        }
                    }
        }
    }
    ret->size_keys = keys_size;
    qsort( ret->keys, ret->size_keys, sizeof(key_value), cmp_keys );
    prev = -1;
    ret->size_ids = 0;
    ret->ids = (INT*) malloc( sizeof(INT) * ret->size_keys*2 );
    for( i = 0; i < ret->size_keys; ++i )
    {
        if( prev != ret->keys[i].key )
        {
            ret->ids[2*ret->size_ids] = i;
            if ( ret->size_ids ) ret->ids[2*(ret->size_ids-1)+1] = i;
            ret->size_ids++;
        }
        prev = ret->keys[i].key;
    }
    if ( ret->size_ids ) ret->ids[2*(ret->size_ids-1)+1] = ret->size_keys;
    return ret;
}

void free_keys( bucket_s* br )
{
    free( br->ids );
    free( br->keys );
    free( br->coords );
    free( br->phantoms );
    free( br );
}

void get_neighbours( bucket_f* br, DOUBLE4* bead, int sizes[64], int* ids[64] )
{
    int count = 0;
    int i,j,k,m;
    int bp[3];
    int p[3];
    int blocal[3];
    int local[3];
    list_reg* plr;
    int li, lj, lk, hi, hj, hk;
    int to_corr = br->mlocal[0] >= 0;

    get_bid_f( br, bead, bp, blocal );
    li = lj = lk = -1;
    hi = hj = hk = 1;
    if ( to_corr )
    {
        if ( br->mlocal[0] > 2 && blocal[0] == 0 ) li = -2;
        if ( br->mlocal[1] > 2 && blocal[1] == 0 ) lj = -2;
        if ( br->mlocal[2] > 2 && blocal[2] == 0 ) lk = -2;
        
        if ( br->mlocal[0] > 2 && blocal[0] == br->mlocal[0]-1 ) hi = 2;
        if ( br->mlocal[1] > 2 && blocal[1] == br->mlocal[1]-1 ) hj = 2;
        if ( br->mlocal[2] > 2 && blocal[2] == br->mlocal[2]-1 ) hk = 2;
        if ( br->mlocal[0] == 0 ) li = hi = 0;
        if ( br->mlocal[1] == 0 ) lj = hj = 0;
        if ( br->mlocal[2] == 0 ) lk = hk = 0;

        if ( br->mlocal[0] == 1 ) li = 0;
        if ( br->mlocal[1] == 1 ) lj = 0;
        if ( br->mlocal[2] == 1 ) lk = 0;
    }
    
    for( i = li; i <= hi; ++i)
    {
        for( j = lj; j <= hj; ++j)
        {
            for( k = lk; k <= hk; ++k)
            {
                if( to_corr )
                {
                    local[0] = correct_indx(blocal[0]+i,br->mlocal[0]+1);
                    local[1] = correct_indx(blocal[1]+j,br->mlocal[1]+1);
                    local[2] = correct_indx(blocal[2]+k,br->mlocal[2]+1);
                    p[0] = correct_indx( local[0], MIN(BUCKET_N,br->mlocal[0]+1) );
                    p[1] = correct_indx( local[1], MIN(BUCKET_N,br->mlocal[1]+1) );
                    p[2] = correct_indx( local[2], MIN(BUCKET_N,br->mlocal[2]+1) );
                }
                else
                {
                    p[0] = correct_indx( bp[0]+i, BUCKET_N );
                    p[1] = correct_indx( bp[1]+j, BUCKET_N );
                    p[2] = correct_indx( bp[2]+k, BUCKET_N );
                    local[0] = blocal[0]+i;
                    local[1] = blocal[1]+j;
                    local[2] = blocal[2]+k;
                }            
                if( local[0] >= 0 && local[1] >= 0 && local[2] >= 0 )
                {
                    plr = &br->lreg[p[0]][p[1]][p[2]];
                    sizes[count]=0;
                    ids[count]=NULL;
                    for( m = 0; m < plr->size; ++m )
                    {
                        if( plr->regs[m].x == local[0] &&
                            plr->regs[m].y == local[1] &&
                            plr->regs[m].z == local[2])
                        {
                            sizes[count] = plr->regs[m].size;
                            ids[count] = plr->regs[m].ids;
                            count ++;
                            break;
                        }
                    }
                }
            }
        }
    }
    if( count < 64 ) sizes[count] = 0;
}

bucket_f* init_bucket_f( int size, DOUBLE4* coords, DOUBLE* radii, DOUBLE blen )
{
    int i,j;
    DOUBLE h[3], l[3];
    int count = 0;
    bucket_f* br = (bucket_f*) malloc(sizeof(bucket_f)); CHMEM(br);
    memset( br, 0, sizeof(bucket_f) );

    get_h_l( size, coords, radii, blen, h, l );

    for(j=0;j<3;++j)
    {
        count += (INT)((h[j]-l[j])/blen);
    }
    if( count <= 6 )/*no sense*/
    {
        free( br );
        return NULL;
    }
    for(j=0;j<3;++j) br->offset[j] = -l[j];
    br->blen[0] = br->blen[1] = br->blen[2] = blen;
    br->mlocal[0] = br->mlocal[1] = br->mlocal[2] = -1;
    if ( bc == BC_BOX )
    {
        for ( i = 0; i < size; ++i )
        {
            for ( j = 0; j < DIMS0; ++j )
            {
                INT tmp = (INT)((coords[i].t[j]+br->offset[j])/br->blen[j]);
                if ( tmp > br->mlocal[j] ) br->mlocal[j] = tmp;
            }
        }
    }
    for ( i = 0; i < size; ++i )
    {
        add_coord_f( br, i, coords+i);
    }
    return br;
}

void get_bid_f( bucket_f* br, DOUBLE4* bead, int* p, int* local )
{
    int j;
    for(j=0;j<3;++j)
    {
        local[j] = (INT)((bead->t[j]+br->offset[j])/br->blen[j]);
        p[j] = local[j] % BUCKET_N;
    }
}

void add_coord_f( bucket_f* br, int id, DOUBLE4* pbead )
{
    INT p[3], j, local[3];
    list_reg* plr;
    get_bid_f( br, pbead, p, local);
    plr = &br->lreg[p[0]][p[1]][p[2]];
    if( plr->size )
    {
        for( j = 0; j < plr->size; ++j )
        {
            if( plr->regs[j].x == local[0] && plr->regs[j].y == local[1] && plr->regs[j].z == local[2] )
                break;
        }
        if( j < plr->size )
        {
            if( plr->regs[j].size + 1 < plr->regs[j].capacity )
            {
                plr->regs[j].ids[plr->regs[j].size++]=id;
            }
            else
            {
                plr->regs[j].capacity = plr->regs[j].capacity*BUCKET_CAP_MLT + BUCKET_CAP_ADD;
                plr->regs[j].ids = (int*)realloc( plr->regs[j].ids, sizeof(int)*plr->regs[j].capacity );
                plr->regs[j].ids[plr->regs[j].size++]=id;
            }
            j = 0;
        }
        else
        {
            j = 1;
        }
    }
    else
    {
        j = 1;
    }
    if( j )
    {
        plr->size++;
        plr->regs = (reg_bucket*)realloc( plr->regs, sizeof(reg_bucket)*plr->size );
        plr->regs[plr->size-1].x = local[0];
        plr->regs[plr->size-1].y = local[1];
        plr->regs[plr->size-1].z = local[2];
        plr->regs[plr->size-1].size=1;
        plr->regs[plr->size-1].capacity=2;
        plr->regs[plr->size-1].ids= (int*) malloc(sizeof(int)*2);
        plr->regs[plr->size-1].ids[0] = id;
    }
}

void free_bucket_f( bucket_f* br, int size, DOUBLE4* coord )
{
    int i;
    if( size > 0 && coord )
    {
        int p[3], local[3];
        for( i = 0; i < size; ++i )
        {
            list_reg* plr;
            get_bid_f(br,coord+i,p,local);
            plr = &br->lreg[p[0]][p[1]][p[2]];
            if( plr->size )
            {
                int j;
                for( j = 0; j < plr->size; ++j )
                {
                    free(plr->regs[j].ids);
                }
                free( plr->regs );
                plr->size = 0;
            }
        }
    }
    else
    {
        int j, k, l;
        for( i = 0; i < BUCKET_N; ++i )
        {
            for( j = 0; j < BUCKET_N; ++j )
            {
                for( k = 0; k < BUCKET_N; ++k )
                {
                    list_reg* plr;
                    plr = &br->lreg[i][j][k];
                    if( plr->size )
                    {
                        for( l = 0; l < plr->size; ++l )
                        {
                            free(plr->regs[l].ids);
                        }
                        free( plr->regs );
                        plr->size = 0;
                    }
                }
            }
        }
    }
    free(br);
}

void print_statistic_f( bucket_f* br )
{
    int i, j, k, l;
    int count = 0;
    int max_reg = 0;
    int max_ids = 0;
    for( i = 0; i < BUCKET_N; ++i )
    {
        for( j = 0; j < BUCKET_N; ++j )
        {
            for( k = 0; k < BUCKET_N; ++k )
            {
                list_reg* plr;
                plr = &br->lreg[i][j][k];
                if( plr->size )
                {
                    for( l = 0; l < plr->size; ++l )
                    {
                        count += plr->regs[l].size;
                        if( plr->regs[l].size > max_ids )
                            max_ids = plr->regs[l].size;
                    }
                    if( plr->size > max_reg )
                        max_reg = plr->size;
                }
            }
        }
    }
    printf("C: %d, mreg: %d, mi: %d\n",count,max_reg,max_ids);
}
