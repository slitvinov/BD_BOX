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

#include "rand_move.h"
#include "err.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#if USE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

MAKE_STR_IN(INT,rand_seed, 12345, "random generator seed")
MAKE_STR_IN(RAND_STATE,rand_state,0,"inner variable")

#if USE_GSL
    gsl_rng* r_gsl;
#else
#if HAVE_DRAND48_R
    struct drand48_data __seed;
#else
    unsigned int __seed;
#endif
#endif

void init_rand( unsigned int seed )
{
#if USE_GSL
    const gsl_rng_type* T;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    
    r_gsl = gsl_rng_alloc(T);
    gsl_rng_set( r_gsl, seed );
#else

#if HAVE_DRAND48_R
    srand48_r( seed, &__seed );
#else

#if HAVE_RAND_R
    __seed = seed;
#else
    __seed = seed;
    srand(seed);
#endif

#endif

#endif
}

void free_rand()
{
#if USE_GSL
    gsl_rng_free(r_gsl);
#endif
}

double rand_d48()
{
#if USE_GSL
    return gsl_rng_uniform(r_gsl);
#else
#if HAVE_DRAND48_R
    double val;
    drand48_r( &__seed, &val );
    return val;
#else
#if HAVE_RAND_R
    long long int x = rand_r( &__seed );
    long long int mx = RAND_MAX;
    double result = 0.0;
    x <<= 16;
    x ^= rand_r( &__seed );
    x <<= 16;
    x ^= rand_r( &__seed ); /*XOR if RAND_MAX > 1<<16 - 1*/
    mx <<= 16;
    mx |= RAND_MAX;
    mx <<= 16;
    mx |= RAND_MAX;
    result = ((double) x / (double) mx);
    return result;
#else
    long long int x = rand( );
    long long int mx = RAND_MAX;
    double result = 0.0;
    x <<= 16;
    x ^= rand( );
    x <<= 16;
    x ^= rand( ); /*XOR if RAND_MAX > 1<<16 - 1*/
    mx <<= 16;
    mx |= RAND_MAX;
    mx <<= 16;
    mx |= RAND_MAX;
    result = ((double) x / (double) mx);
    return result;
#endif
#endif
#endif
}

void gauss_rand( DOUBLE *tabl, INT tabl_size )
{
    int i;
#if USE_GSL
    for( i = 0; i < tabl_size; ++i )
    {
        tabl[i] = gsl_ran_gaussian(r_gsl,1.0);
    }
#else
    int iset = 0;
    DOUBLE gset;
    DOUBLE fac, rsq, v1, v2;
    v1 = v2 = gset = 0.0;
    for (i = 0; i < tabl_size; i++)
    {
        if (iset == 0)
        {
            rsq = 0.0;
            while ((rsq >= 1) || (rsq == 0.0))
            {
                v1 = 2.0f * (DOUBLE)rand_d48() - 1.0f;
                v2 = 2.0f * (DOUBLE)rand_d48() - 1.0f;
                rsq = v1 * v1 + v2*v2;
            }
            fac = SQRTD(-2.0f * LOGD(rsq) / rsq);
            gset = v1*fac;
            iset = 1;
            tabl[i] = v2*fac;
        }
        else
        {
            iset = 0;
            tabl[i] = gset;
        }
    }
#endif
}

void get_state( char** state, INT* size )
{
#if USE_GSL
    *state = gsl_rng_state( r_gsl );
    *size = gsl_rng_size( r_gsl );
#else
    *state = (char*)&__seed;
    *size = sizeof( __seed );
#endif
}

void set_state( char* state, INT size )
{
#if USE_GSL
    const gsl_rng_type* T;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r_gsl = gsl_rng_alloc(T);
    if( size != gsl_rng_size( r_gsl ) )
        UNERR( "random generators state sizes differ");
    memcpy( gsl_rng_state( r_gsl ), state, size );
#else
    if( sizeof( __seed ) != size )
        UNERR( "random generators state sizes differ");
    memcpy( &__seed, state, size );
#endif
}
