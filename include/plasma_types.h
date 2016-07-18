/**
 *
 * @file plasma_types.h
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver.
 *
 * @version 3.0.0
 * @author  Jakub Kurzak
 * @author  Maksims Abalenkovs
 * @date    2016-07-18
 *
 **/
#ifndef ICL_PLASMA_TYPES_H
#define ICL_PLASMA_TYPES_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
#ifdef PLASMA_WITH_MKL
#define lapack_complex_float PLASMA_Complex32_t
#define lapack_complex_double PLASMA_Complex64_t
#endif

/******************************************************************************/
#ifndef CBLAS_SADDR
#define CBLAS_SADDR(var) &(var)
#endif

/******************************************************************************/
enum {
    PlasmaByte          = 0,
    PlasmaInteger       = 1,
    PlasmaRealFloat     = 2,
    PlasmaRealDouble    = 3,
    PlasmaComplexFloat  = 4,
    PlasmaComplexDouble = 5
};

/***************************************************************************//**
 *
 *  PLASMA constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
enum {
    PlasmaNoTrans    = 111,
    PlasmaTrans      = 112,
    PlasmaConjTrans  = 113,

    PlasmaUpper      = 121,
    PlasmaLower      = 122,
    PlasmaFull       = 123,  // formerly PlasmaUpperLower

    PlasmaNonUnit    = 131,
    PlasmaUnit       = 132,

    PlasmaLeft       = 141,
    PlasmaRight      = 142,

    PlasmaInfNorm    = 175
};

enum {
    PLASMA_SUCCESS              =    0,
    PLASMA_ERR_NOT_INITIALIZED  = -101,
    PLASMA_ERR_REINITIALIZED    = -102,
    PLASMA_ERR_NOT_SUPPORTED    = -103,
    PLASMA_ERR_ILLEGAL_VALUE    = -104,
    PLASMA_ERR_NOT_FOUND        = -105,
    PLASMA_ERR_OUT_OF_RESOURCES = -106,
    PLASMA_ERR_INTERNAL_LIMIT   = -107,
    PLASMA_ERR_UNALLOCATED      = -108,
    PLASMA_ERR_FILESYSTEM       = -109,
    PLASMA_ERR_UNEXPECTED       = -110,
    PLASMA_ERR_SEQUENCE_FLUSHED = -111
};

enum {
    PLASMA_INPLACE,
    PLASMA_OUTOFPLACE
};

enum {
    PLASMA_TILE_SIZE,
    PLASMA_TRANSLATION_MODE
};

/******************************************************************************/
typedef int PLASMA_enum;
typedef int PLASMA_bool;

typedef float  _Complex PLASMA_Complex32_t;
typedef double _Complex PLASMA_Complex64_t;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_TYPES_H
