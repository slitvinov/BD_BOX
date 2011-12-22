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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "myblas.h"
#include "err.h"

#include <memory.h>
#include <math.h>
#include <stdlib.h>

#if USE_CUDA || USE_MAGMA
#include <cuda_runtime_api.h>
#include <cublas.h>
#include "cuda.h"
#endif

#if !USE_LAPACK
static DOUBLE* p;
#endif

INT my_dpotrf( DOUBLE* A, INT n, INT LD )
{
#if USE_CUDA
#if USE_MAGMA
    int INFO[1];
    INFO[0] = -1234567;
    MAGMA_POTRFD( 'L', n, A, LD, INFO );
    return INFO[0];
#else
    UNERR("CUDA version without MAGMA - can't use cholesky decomposition");
    return 1;
#endif
#else
#if USE_LAPACK
    int INFO[1];
    POTRFD( "U", &n, A, &LD, INFO );
    return INFO[0];
#else
    int i, j, k;
    DOUBLE sum;
    for (i = 0; i < n; ++i)
    {
        for (j = i; j < n; ++j)
        {
            for (sum = A[ i + j * LD ], k = i - 1; k >= 0; k--)
                sum -= A[ i + k * LD ] * A[ j + k * n ];
            if (i == j)
            {
                if (sum <= 0.0)
                {
                    return -1;
                }
                p[i] = SQRTD(sum);
            }
            else
            {
                A[ j + n * i ] = sum / p[i];
            }
        }
    }
    for (i = 0; i < n; ++i)
    {
        A[ i * LD + i ] = p[i];
    }
    return 0;
#endif
#endif
}

void my_dtrmv( DOUBLE* A, DOUBLE* v, DOUBLE* o, INT n, INT LD )
{
#if USE_CUDA
    CUBLAS_TRMVD( 'L', 'N', 'N', n, A, LD, v, 1 ); cublasCheckError(__FILE__,__LINE__);
    CUBLAS_AXPYD( n, 1.0f, v, 1, o, 1 ); cublasCheckError(__FILE__,__LINE__);
#else
#if USE_LAPACK
    int ione = 1;
    TRMVD( "U", "T", "N", &n, A, &LD, v, &ione );
    memcpy( o, v, n*sizeof(DOUBLE) );
#else
    DOUBLE sum;
    INT i, j;
    for( i = 0; i < n; ++i)
    {
        sum = 0.0;
        for( j = 0; j <= i; ++j )
        {
            sum += A[ LD*j + i ] * v[ j ];
        }
        o[ i ] = sum;
    }
#endif
#endif
}

void init_cholesky()
{
#if !USE_LAPACK
    p = (DOUBLE*) malloc( sizeof(DOUBLE) * DIMS0 * size );
#endif
}

void free_cholesky()
{
#if !USE_LAPACK
    free(p);
#endif
}
