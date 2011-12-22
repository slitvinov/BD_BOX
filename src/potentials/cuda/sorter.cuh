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

#ifndef SORTER_CUH
#define SORTER_CUH

class StaticSorter_impl;

/*! Sorte class for device arrays */
class StaticSorter
{
private:
    
    /*! Internal object */
    StaticSorter_impl* o_;
public:

    /*! Initialization of object. Bind array to sorter.
        \param d_keys keys to sort of type unsigned int
        \param d_vals data binded to keys of 'int' size
        \param size number of elements */
    StaticSorter( void* d_keys, void* d_vals, int size );
    
    /*! Free resources*/
    ~StaticSorter();

    /*! Sort arrays by keys. */
    void sort();
};

#endif
