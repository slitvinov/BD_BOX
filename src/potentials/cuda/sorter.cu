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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef USE_CUDPP

#include <cudpp.h>
extern "C"
{
#include "../../cuda.h"
#include "../../err.h"
}
#include "sorter.cuh"

#include <iostream>
using namespace std;

#define CERR cudaCheckError(__FILE__,__LINE__);

#ifdef WIN32
typedef unsigned uint;
#endif

void cudppCheck( CUDPPResult res, const char* file, int line );

/*! Intenal implementation of sorting wrapper */
class StaticSorter_impl
{
private:

    /*! */
    uint* d_keys_;
    uint* d_vals_;
    int size_;
    CUDPPHandle plan_;
    CUDPPConfiguration config_;
public:
    StaticSorter_impl( void* d_keys, void* d_vals, int size ): size_(size)
    {
        d_keys_ = static_cast<uint*>(d_keys);
        d_vals_ = static_cast<uint*>(d_vals);
        memset( (void*) &config_, 0, sizeof(CUDPPConfiguration) );
        config_.algorithm = CUDPP_SORT_RADIX;
        config_.op = CUDPP_MIN;
        config_.datatype = CUDPP_UINT;

        cudppCheck( cudppPlan( &plan_, config_, size, 1, 0 ), __FILE__, __LINE__ );
    }

    void sort( )
    {
        cudppSort( plan_, d_keys_, d_vals_, 32, size_ );
        cudaThreadSynchronize( );
        CERR
    }

    ~StaticSorter_impl( )
    {
        cudppCheck( cudppDestroyPlan( plan_ ), __FILE__, __LINE__ );
    }
};

StaticSorter::StaticSorter( void* d_keys, void* d_vals, int size )
{
    o_ = new StaticSorter_impl( d_keys, d_vals, size );
}

StaticSorter::~StaticSorter( )
{
    delete o_;
}

void StaticSorter::sort( )
{
    o_->sort( );
}

/*! Check error of cudpp and raise unexpected_error.
    \sa unexpected_error
    \param res result from cudpp routine
    \param file source of error
    \param line number of line that caused error */
void cudppCheck( CUDPPResult res, const char* file, int line )
{
    if ( res != CUDPP_SUCCESS )
    {
        const char* reason = "Unknown";
        switch ( res )
        {
            case ( CUDPP_ERROR_INVALID_HANDLE ): reason = "Invalid handle";
                break;
            case ( CUDPP_ERROR_ILLEGAL_CONFIGURATION ): reason = "Illegal configuration";
                break;
        }
        UNERR(reason);
    }
}

#endif
