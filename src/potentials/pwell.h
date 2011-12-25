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

#include "pwell.h"

#include "../trans.h"
#include "../math_help.h"

#define EPSILON 0.000001

MAKE_STR_IN(YESNO,  pwell,   0, "yes/no - switch the bounding sphere field on/off")
MAKE_STR_IN(DOUBLE, pwell_A, 1, "the magnitude of the bounding sphere force")
MAKE_STR_IN(INT,    pwell_n, 2, "the power of the radial distance dependence of the bounding sphere force")
MAKE_STR_IN(DOUBLE, pwell_cutoff, 0, "the bounding sphere force is applied outside this cutoff radius" )

INT pwell_ext( DOUBLE* coord, DOUBLE* F )
{
    DOUBLE r = norm3v( coord );
    if ( r >= pwell_cutoff )
    {
        DOUBLE R = sphere_radius - r;
        if( R >= EPSILON )
        {
            DOUBLE e[3];
            DOUBLE inr = 1.0f / r;
            DOUBLE inR = 1.0f / R;
            INT i;
            for ( i = 1; i <= pwell_n; i*=2 )
            {
                if ( pwell_n & i )
                    inr *= inR;
                inR *= inR;
            }
            e[0] = pwell_A * coord[0] * inr;
            e[1] = pwell_A * coord[1] * inr;
            e[2] = pwell_A * coord[2] * inr;

            F[0] -= e[0];
            F[1] -= e[1];
            F[2] -= e[2];
            
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else
    {
        return 0;
    }
}
