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

/**
 * \mainpage
 *
 * \section introduction Introduction
 *
 * BD_BOX is a scalable Brownian dynamics package for UNIX/LINUX platforms.
 * BD BOX is written in C and uses modern computer architectures and technologies.
 *
 * \section overview Overview Presentation
 *
 * A brief set of description is available here:
 * <a href="http://bionano.icm.edu.pl">Bionano</a>.
 *
 * \section homepage Homepage
 * Homepage for BD_BOX: http://bionano.icm.edu.pl
 *
 *
 * \section getting-started Getting Started with BD_BOX
 * 
 * The BD BOX package is distributed as source code only.
 * A UNIX/LINUX make tool is needed to build a working binary from source code (see below).
 * The CPU version of BD BOX can be run either in a serial mode or in parallel,
 * either on shared-memory machines using OpenMP and MPI libraries or on 
 * architectures with distributed memory using the MPI library.
 * The GPU version of BD BOX requires the CUDA Toolkit.
 *
 * \section getting-help Getting Help and Reporting Problems
 *
 * To get help using BD_BOX, please use the manual.
 *
 * \section opSys Operating System Support
 *
 * This release (1.0) has been thoroughly tested on the following OSes.
 * - Windows XP (32-bit) (CUDA 3.2)
 * - Windows 7 (64-bit) (CUDA 3.2)
 * - Ubuntu Linux 10.04 (32-bit and 64-bit) (CUDA 3.2)
 * - Ubuntu Linux  9.10 (32-bit and 64-bit) (CUDA 3.2)
 *
 * We expect BD_BOX to build and run correctly on other flavors of Linux, but these
 * are not actively tested by the developers at this time.
 *
 * Note: BD_BOX is not compatible with previous CUDA releases.
 *
 * \section citing Citing BD_BOX
 *
 *
 * \section credits Credits
 * \subsection managers  BD_BOX Project manager
 * - Maciej Dlugosz,    University of Warsaw, Warsaw
 * \subsection developers BD_BOX Developers
 * - Pawel Zielinski,   University of Warsaw, Warsaw
 *
 * \subsection acknowledgments Acknowledgments
 *
 * Thanks to .
 *
 * BD_BOX Developers from University of Warsaw thank their funding agencies:
 * -
 *
 * \section license-overview BD_BOX Copyright and Software License
 *   This program is free software; you can redistribute it and/or modify <BR>
 *   it under the terms of the GNU General Public License as published by <BR>
 *   the Free Software Foundation; either version 3 of the License, or <BR>
 *   (at your option) any later version. <BR>
 *   <BR>
 *   This program is distributed in the hope that it will be useful,  <BR>
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of  <BR>
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  <BR>
 *   GNU General Public License for more details.  <BR>
 *
 *   You should have received a copy of the GNU General Public License  <BR>
 *   along with this program; if not, see http://www.gnu.org/licenses  <BR>
 *   or write to the Free Software Foundation,Inc., 51 Franklin Street,  <BR>
 *   Fifth Floor, Boston, MA 02110-1301  USA  <BR>
 */

#ifndef DATA_H
#define DATA_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if USE_MPI
#include "mpi.h"
#endif

#if !HAVE_ERFC

/*! Error function, double precision version. Used mainly for windows platform. */
double erfc(double);

/*! Error function, single precision version. Used mainly for windows platform. */
float  erfcf(float);
#endif

#if USE_FLOAT
typedef float DOUBLE;
#define FORMAT_DOUBLEG "g"
#define FORMAT_DOUBLESG "g"
#define FORMAT_DOUBLEF "f"
#define SQRTD sqrtf
#define LOGD logf
#define EXPD expf
#define ERFCD erfcf
#define FLOORD floorf
#define COSD cosf
#define ACOSD acosf
#define SIND sinf
#define ATAN2D atan2f
#define POTRFD spotrf_
#define TRMVD strmv_
#define SYMVD ssymv_

#define SSECOUNT 4
#define M128 __m128
#define SSESHUFFLE 57
#define MMSET1 _mm_set_ps1
#define MMSETZERO _mm_setzero_ps
#define MMSQRT _mm_sqrt_ps
#define MMRSQRT _mm_rsqrt_ps
#define MMSHUFFLE _mm_shuffle_ps
#define MMSUB _mm_sub_ps
#define MMADD _mm_add_ps
#define MMMUL _mm_mul_ps
#define MMDIV _mm_div_ps
#define MMCMPLT _mm_cmplt_ps
#define MMSTORE _mm_store_ps
#define MEXP exp_ps

#define DOUBLE_MAX FLT_MAX

#define MAGMA_POTRFD magma_spotrf_gpu
#define CUBLAS_TRMVD cublasStrmv
#define CUBLAS_SYMVD cublasSsymv
#define CUBLAS_SSCALD cublasSscal
#define CUBLAS_AXPYD cublasSaxpy

#define NEWTON_EPSILON 0.00001f
#define NEWTON_OFFSET  0.01f

#if USE_MPI
#define MPI_VALUE_TYPE MPI_FLOAT
#endif

#define kB 0.001986248532617f
#define kB_1_8PI 0.00007903031804378703f
#define kB_1_6PI 0.00010537375739171604f
#define SQRTPI 1.77245385090551602729f
#define M1_SQRTPI 0.5641895835f

#else
typedef double DOUBLE;
#define FORMAT_DOUBLEG "g"
#define FORMAT_DOUBLESG "lg"
#define FORMAT_DOUBLEF "f"
#define SQRTD sqrt
#define LOGD log
#define EXPD exp
#define ERFCD erfc
#define FLOORD floor
#define COSD cos
#define ACOSD acos
#define SIND sin
#define ATAN2D atan2
#define POTRFD dpotrf_
#define TRMVD dtrmv_
#define SYMVD dsymv_

#define SSECOUNT 2
#define M128 __m128d
#define SSESHUFFLE 1
#define MMSET1 _mm_set1_pd
#define MMSETZERO _mm_setzero_pd
#define MMSQRT _mm_sqrt_pd
#define MMRSQRT _mm_rsqrt_pd
#define MMSHUFFLE _mm_shuffle_pd
#define MMSUB _mm_sub_pd
#define MMADD _mm_add_pd
#define MMMUL _mm_mul_pd
#define MMDIV _mm_div_pd
#define MMCMPLT _mm_cmplt_pd
#define MMSTORE _mm_store_pd
#define MEXP exp_pd

#define DOUBLE_MAX DBL_MAX

#define MAGMA_POTRFD magma_dpotrf_gpu
#define CUBLAS_TRMVD cublasDtrmv
#define CUBLAS_SYMVD cublasDsymv
#define CUBLAS_SSCALD cublasDscal
#define CUBLAS_AXPYD cublasDaxpy

#define NEWTON_EPSILON 0.000001
#define NEWTON_OFFSET  0.0001

#if USE_MPI
#define MPI_VALUE_TYPE MPI_DOUBLE
#endif

#define kB 0.001986248532617
#define kB_1_8PI 0.00007903031804378703
#define kB_1_6PI 0.00010537375739171604
#define SQRTPI 1.77245385090551602729
#define M1_SQRTPI 0.5641895835

#endif

typedef float FLOAT;
typedef int INT;
typedef const char* CSTR;
typedef char CHAR;
typedef unsigned char UCHAR;
typedef CHAR* STR;
typedef INT YESNO;
typedef INT DC_AC_RF;
typedef INT NO_CHOLS_GEYER;
typedef INT NONE_CUBIC_SPHERE;
typedef INT ERMAK_IGTCONST_IGTVAR;
typedef INT EWALD_METHOD;
typedef INT RAND_STATE;
typedef INT BUCKET_NONE_SPATIAL;
typedef DOUBLE* PTR_DOUBLE;
typedef DOUBLE* PTR_DOUBLE9;
typedef INT* PTR_INT;
typedef STR* PTR_STR;
typedef void VOID;

#define DIMS0 3
#define DIMS1 4
#define M_PIF 3.14159265358979323846f
#define M1_SQRTPIF 0.5641895835f

#define SUB_TAG "sub"
#define BOND_TAG "bond"
#define ANGLE_TAG "angle"
#define DIHE_TAG "dihe"

/*! Structure for single variable, wrapper. */
struct variable
{
    /*! Pointer to variable*/
    void* ptr;
    /*! Name in program of variable */
    CSTR name;
    /*! Type of variable */
    CSTR type;
    /*! Documentation of variable, used in print_help */
    CSTR doc;
    /*! Size of variable, used in array variables */
    INT size;
};

/*! Types of energies in system. */
enum ENERGY_TYPE
{
    ENERGY_NONE,
    ENERGY_BOND,
    ENERGY_ANGLE,
    ENERGY_ANGLE_COS,
    ENERGY_DIHE,
    ENERGY_DIHE_ANGLE,
    ENERGY_LJ,
    ENERGY_COULOMB,
    ENERGY_LAST
};

/*! Macro for declaring variable, that can be parsed from string and saved in restart file. */
#define MAKE_STR_DEC(TYPE,NAME) extern TYPE NAME; extern struct variable variable_##NAME;

/*! Macro for declaring variable, saved in restart file. */
#define MAKE_STR_DEC_INNER(TYPE,NAME) MAKE_STR_DEC(TYPE,NAME)
#define COMMENT_SIGN '#'
#define MAKE_STR_IN(TYPE,NAME,DEFVAL,DOC) TYPE NAME=DEFVAL; struct variable variable_##NAME = { & NAME, #NAME, #TYPE, DOC, 0 };
#define MAKE_STR_SIZE(NAME,SIZE) variable_##NAME.size = SIZE;

/*! Make pointer to wrapper of variable. */
#define MAKE_STR_PTR(NAME) &variable_##NAME
#define TO_STR( a ) #a

/*! for MPI - IF Primary process*/
#define IFP if(g_id==0)

/*! for MPI - IF NOT Primary process*/
#define IFNP if(g_id!=0)

#ifdef	__cplusplus
extern "C" {
#endif

/*! Number of beads */
extern INT size;
MAKE_STR_DEC_INNER(INT,size)

/*! Coordinates array - x,y,z,radii */
extern PTR_DOUBLE coord;
MAKE_STR_DEC_INNER(PTR_DOUBLE,coord)

/*! Optional masses of beads used in newton for collision. */
extern PTR_DOUBLE masses;
MAKE_STR_DEC_INNER(PTR_DOUBLE,masses)

/*! Charges of beads */
extern PTR_DOUBLE Q;
MAKE_STR_DEC_INNER(PTR_DOUBLE,Q)

/*! Lennard Jones variables */
extern PTR_DOUBLE LJ;
MAKE_STR_DEC_INNER(PTR_DOUBLE,LJ)

/*! Names of beads */
extern PTR_STR names;
MAKE_STR_DEC_INNER(PTR_STR,names)

/*! Ids of beads */
extern PTR_INT ids;
MAKE_STR_DEC_INNER(PTR_INT,ids)

/*! Box size. */
extern DOUBLE box[4];

/*! Number of threads in single process, for OpenMP. */
extern INT g_threads;

/*! Id of process, for MPI. */
extern INT g_id;

/*! Number of process in MPI. */
extern INT g_numprocs;

/*! Number of cuda devices. */
extern INT g_devices;

/*! Copy of coordinates. */
extern DOUBLE* save_coord;

/*! Copy of forces. */
extern DOUBLE* save_forces;

/*! Convert id to array offset. */
INT resolve_pos(INT id);

/*! Add beads to coordinates, Q, LJ, names and ids arrays. */
void add_sub(CSTR namei, INT idi, DOUBLE xi, DOUBLE yi, DOUBLE zi, DOUBLE sigmai, DOUBLE Qi, DOUBLE sigma_LJi, DOUBLE epsilon_LJi, DOUBLE mass );

/*! Trim arrays to number of beads. */
void end_sub();

/*! Data and reqources free. */
void free_data( );

#ifdef	__cplusplus
}
#endif

#endif	/* DATA_H */
