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

#include "dihe_angle.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../err.h"
#include "../math_help.h"
#include "../trans.h"

#define EPSILON 0.000001

MAKE_STR_IN(INT,dihe_angle_size,0,"inner variable")
MAKE_STR_IN(PTR_INT,dihe_angle_pairs,NULL,"inner variable")
MAKE_STR_IN(PTR_DOUBLE,dihe_angle_parms,NULL,"inner variable")

void measureAndCalcGradDihedral(DOUBLE* points1, DOUBLE* points2, DOUBLE* points3,
    DOUBLE* points4,DOUBLE* outDihedral,DOUBLE* outGradDihedral);

INT dihe_calc_angle( DOUBLE* F, int pos, DOUBLE* _E )
{
    DOUBLE* F1 = F + DIMS0 * dihe_angle_pairs[ 4 * pos ];
    DOUBLE* F2 = F + DIMS0 * dihe_angle_pairs[ 4 * pos + 1 ];
    DOUBLE* F3 = F + DIMS0 * dihe_angle_pairs[ 4 * pos + 2 ];
    DOUBLE* F4 = F + DIMS0 * dihe_angle_pairs[ 4 * pos + 3 ];
    DOUBLE* ri = coord + DIMS1 * dihe_angle_pairs[ 4 * pos ];
    DOUBLE* rj = coord + DIMS1 * dihe_angle_pairs[ 4 * pos + 1];
    DOUBLE* rk = coord + DIMS1 * dihe_angle_pairs[ 4 * pos + 2];
    DOUBLE* rl = coord + DIMS1 * dihe_angle_pairs[ 4 * pos + 3];
    DOUBLE oD[12];

    DOUBLE k_theta = dihe_angle_parms[ 2*pos    ];
    DOUBLE theta_0 = dihe_angle_parms[ 2*pos + 1];
    DOUBLE theta = 0;
    DOUBLE tmp;
    
    measureAndCalcGradDihedral(ri,rj,rk,rl,&theta,oD);
    tmp = (theta - theta_0);
    if( tmp < -M_PIF) tmp += 2*M_PIF;
    if( tmp > M_PIF) tmp -= 2*M_PIF;
    k_theta *= -tmp;
    /*V*/
    if(_E) _E[0] -= 0.5f * k_theta * tmp;

    F1[0] += k_theta * oD[0];
    F1[1] += k_theta * oD[1];
    F1[2] += k_theta * oD[2];

    F2[0] += k_theta * oD[3];
    F2[1] += k_theta * oD[4];
    F2[2] += k_theta * oD[5];

    F3[0] += k_theta * oD[6];
    F3[1] += k_theta * oD[7];
    F3[2] += k_theta * oD[8];

    F4[0] += k_theta * oD[9];
    F4[1] += k_theta * oD[10];
    F4[2] += k_theta * oD[11];

    return isnan(oD[0]) || isinf(oD[0]) ||
           isnan(oD[1]) || isinf(oD[1]) ||
           isnan(oD[2]) || isinf(oD[2]) ||
           isnan(oD[3]) || isinf(oD[3]) ||
           isnan(oD[4]) || isinf(oD[4]) ||
           isnan(oD[5]) || isinf(oD[5]) ||
           isnan(oD[6]) || isinf(oD[6]) ||
           isnan(oD[7]) || isinf(oD[7]) ||
           isnan(oD[8]) || isinf(oD[8]) ||
           isnan(oD[9]) || isinf(oD[9]) ||
           isnan(oD[10]) || isinf(oD[10]) ||
           isnan(oD[11]) || isinf(oD[11]) ||
           isnan(k_theta) || isinf(k_theta) ||
           ( _E ? (isnan( _E[0] ) || isinf( _E[0] )) : 0);
}

/* From RedMD */
void measureAndCalcGradDihedral(DOUBLE* points1, DOUBLE* points2, DOUBLE* points3,
    DOUBLE* points4,DOUBLE* outDihedral,DOUBLE* outGradDihedral)
{
    DOUBLE K1=1;

    DOUBLE r12[3];
    DOUBLE r23[3];
    DOUBLE r34[3];

    DOUBLE A[3],B[3],C[3];
    DOUBLE rAinv;
    DOUBLE rBinv;
    DOUBLE rCinv;

    DOUBLE cos_phi;
    DOUBLE sin_phi;

    DOUBLE f1[3],f2[3],f3[3];
    DOUBLE phi;
    r12[0] = points2[0]-points1[0];
    r12[1] = points2[1]-points1[1];
    r12[2] = points2[2]-points1[2];

    r23[0] = points3[0]-points2[0];
    r23[1] = points3[1]-points2[1];
    r23[2] = points3[2]-points2[2];

    r34[0] = points4[0]-points3[0];
    r34[1] = points4[1]-points3[1];
    r34[2] = points4[2]-points3[2];

    cross_product(r12,r23,A);
    cross_product(r23,r34,B);
    cross_product(r23,A,C);
    rAinv = 1.0f / norm3v(A);
    rBinv = 1.0f / norm3v(B);
    rCinv = 1.0f / norm3v(C);

    cos_phi = (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])*(rAinv*rBinv);
    sin_phi = (C[0]*B[0]+C[1]*B[1]+C[2]*B[2])*(rCinv*rBinv);
    if( cos_phi < -1 ) cos_phi = -1;
    if( cos_phi >  1 ) cos_phi =  1;
    if( sin_phi < -1 ) sin_phi = -1;
    if( sin_phi >  1 ) sin_phi =  1;
    phi= -ATAN2D(sin_phi,cos_phi);
    (*outDihedral)=phi;

    B[0] *= rBinv;
    B[1] *= rBinv;
    B[2] *= rBinv;

    if ( ABSD(sin_phi) > 0.1f )
    {
        DOUBLE dcosdA[3];
        DOUBLE dcosdB[3];
        
        A[0] *= rAinv;
        A[1] *= rAinv;
        A[2] *= rAinv;

        dcosdA[0] = rAinv*(cos_phi*A[0]-B[0]);
        dcosdA[1] = rAinv*(cos_phi*A[1]-B[1]);
        dcosdA[2] = rAinv*(cos_phi*A[2]-B[2]);

        dcosdB[0] = rBinv*(cos_phi*B[0]-A[0]);
        dcosdB[1] = rBinv*(cos_phi*B[1]-A[1]);
        dcosdB[2] = rBinv*(cos_phi*B[2]-A[2]);

        K1 = K1/sin_phi;
        f1[0] = K1*(r23[1]*dcosdA[2] - r23[2]*dcosdA[1]);
        f1[1] = K1*(r23[2]*dcosdA[0] - r23[0]*dcosdA[2]);
        f1[2] = K1*(r23[0]*dcosdA[1] - r23[1]*dcosdA[0]);

        f3[0] = K1*(r23[2]*dcosdB[1] - r23[1]*dcosdB[2]);
        f3[1] = K1*(r23[0]*dcosdB[2] - r23[2]*dcosdB[0]);
        f3[2] = K1*(r23[1]*dcosdB[0] - r23[0]*dcosdB[1]);

        f2[0] = K1*(r12[2]*dcosdA[1] - r12[1]*dcosdA[2]
                + r34[1]*dcosdB[2] - r34[2]*dcosdB[1]);
        f2[1] = K1*(r12[0]*dcosdA[2] - r12[2]*dcosdA[0]
                + r34[2]*dcosdB[0] - r34[0]*dcosdB[2]);
        f2[2] = K1*(r12[1]*dcosdA[0] - r12[0]*dcosdA[1]
                + r34[0]*dcosdB[1] - r34[1]*dcosdB[0]);
  }
  else
  {
    /*  This angle is closer to 0 or 180 than it is to
      90, so use the cos version to avoid 1/sin terms

      Normalize C
      rC = 1.0/rC;*/

    DOUBLE dsindC[3];
    DOUBLE dsindB[3];
    
    C[0] *= rCinv;
    C[1] *= rCinv;
    C[2] *= rCinv;

    dsindC[0] = rCinv*(sin_phi*C[0]-B[0]);
    dsindC[1] = rCinv*(sin_phi*C[1]-B[1]);
    dsindC[2] = rCinv*(sin_phi*C[2]-B[2]);

    dsindB[0] = rBinv*(sin_phi*B[0]-C[0]);
    dsindB[1] = rBinv*(sin_phi*B[1]-C[1]);
    dsindB[2] = rBinv*(sin_phi*B[2]-C[2]);

    K1 = -K1/cos_phi;

    f1[0] = K1*((r23[1]*r23[1] + r23[2]*r23[2])*dsindC[0]
                - r23[0]*r23[1]*dsindC[1]
                - r23[0]*r23[2]*dsindC[2]);
    f1[1] = K1*((r23[2]*r23[2] + r23[0]*r23[0])*dsindC[1]
                - r23[1]*r23[2]*dsindC[2]
                - r23[1]*r23[0]*dsindC[0]);
    f1[2] = K1*((r23[0]*r23[0] + r23[1]*r23[1])*dsindC[2]
                - r23[2]*r23[0]*dsindC[0]
                - r23[2]*r23[1]*dsindC[1]);

    cross_product(dsindB,r23,f3);

    f3[0]*=K1;
    f3[1]*=K1;
    f3[2]*=K1;

    f2[0] = K1*(-(r23[1]*r12[1] + r23[2]*r12[2])*dsindC[0]
                +(2.0f*r23[0]*r12[1] - r12[0]*r23[1])*dsindC[1]
                +(2.0f*r23[0]*r12[2] - r12[0]*r23[2])*dsindC[2]
                +dsindB[2]*r34[1] - dsindB[1]*r34[2]);
    f2[1] = K1*(-(r23[2]*r12[2] + r23[0]*r12[0])*dsindC[1]
                +(2.0f*r23[1]*r12[2] - r12[1]*r23[2])*dsindC[2]
                +(2.0f*r23[1]*r12[0] - r12[1]*r23[0])*dsindC[0]
                +dsindB[0]*r34[2] - dsindB[2]*r34[0]);
    f2[2] = K1*(-(r23[0]*r12[0] + r23[1]*r12[1])*dsindC[2]
                +(2.0f*r23[2]*r12[0] - r12[2]*r23[0])*dsindC[0]
                +(2.0f*r23[2]*r12[1] - r12[2]*r23[1])*dsindC[1]
                +dsindB[1]*r34[0] - dsindB[0]*r34[1]);
  }

  outGradDihedral[0]=f1[0];
  outGradDihedral[1]=f1[1];
  outGradDihedral[2]=f1[2];

  outGradDihedral[3] = f2[0] - f1[0];
  outGradDihedral[4] = f2[1] - f1[1];
  outGradDihedral[5] = f2[2] - f1[2];

  outGradDihedral[6] = f3[0] - f2[0];
  outGradDihedral[7] = f3[1] - f2[1];
  outGradDihedral[8] = f3[2] - f2[2];

  outGradDihedral[9] = -f3[0];
  outGradDihedral[10] = -f3[1];
  outGradDihedral[11] = -f3[2];
}

static INT capacity = 0;
void add_dihe_angle( INT id1, INT id2, INT id3, INT id4, DOUBLE k_theta, DOUBLE theta )
{
    if( capacity == dihe_angle_size )
    {
        capacity = capacity ? 2*capacity : 128;
        dihe_angle_pairs = (INT*) realloc( dihe_angle_pairs, sizeof(INT)*4*capacity );
        CHMEM(dihe_angle_pairs);
        dihe_angle_parms = (DOUBLE*) realloc( dihe_angle_parms, sizeof(DOUBLE)*2*capacity );
        CHMEM(dihe_angle_parms);
    }
    dihe_angle_pairs[4*dihe_angle_size] = id1;
    dihe_angle_pairs[4*dihe_angle_size + 1] = id2;
    dihe_angle_pairs[4*dihe_angle_size + 2] = id3;
    dihe_angle_pairs[4*dihe_angle_size + 3] = id4;
    dihe_angle_parms[2*dihe_angle_size] = k_theta;
    dihe_angle_parms[2*dihe_angle_size + 1] = theta*M_PIF/180;
    dihe_angle_size++;
}

void end_dihe_angle()
{
    if( dihe_angle_size )
    {
        dihe_angle_parms = (DOUBLE*) realloc( dihe_angle_parms, sizeof(DOUBLE) * 2 * dihe_angle_size );
        MAKE_STR_SIZE(dihe_angle_parms,2*dihe_angle_size)
        dihe_angle_pairs = (INT*) realloc( dihe_angle_pairs, sizeof(INT) * 4 * dihe_angle_size );
        MAKE_STR_SIZE(dihe_angle_pairs,4*dihe_angle_size)
    }
}

void free_dihe_angle()
{
    if( dihe_angle_parms )
    {
        free( dihe_angle_parms );
        dihe_angle_parms = NULL;
    }
    if( dihe_angle_pairs )
    {
        free( dihe_angle_pairs );
        dihe_angle_pairs = NULL;
    }
}
