/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef ICL_CORE_BLAS_H
#define ICL_CORE_BLAS_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

static const char *lapack_constants[] = {
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",

    "", "", "", "", "", "", "", "", "", "",
    "",
    "NoTrans",                              ///< 111: PlasmaNoTrans
    "Trans",                                ///< 112: PlasmaTrans
    "ConjTrans",                            ///< 113: PlasmaConjTrans

    "", "", "", "", "", "", "",
    "Upper",                                ///< 121: PlasmaUpper
    "Lower",                                ///< 122: PlasmaLower
    "General",                              ///< 123: PlasmaFull

    "", "", "", "", "", "", "",
    "NonUnit",                              ///< 131: PlasmaNonUnit
    "Unit",                                 ///< 132: PlasmaUnit

    "", "", "", "", "", "", "", "",
    "Left",                                 ///< 141: PlasmaLeft
    "Right",                                ///< 142: PlasmaRight
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "",
    "Forward",                             ///< 391: PlasmaForward
    "Backward",                            ///< 392: PlasmaBackward
    "", "", "", "", "", "", "", "",
    "Columnwise",                          ///< 401: PlasmaColumnwise
    "Rowwise"                              ///< 402: PlasmaRowwise
};

static inline char lapack_const(int plasma_const) {
    return lapack_constants[plasma_const][0];
}

#define coreblas_error(msg) \
        coreblas_error_func_line_file(__func__, __LINE__, __FILE__, msg)

static inline void coreblas_error_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "COREBLAS ERROR at %d of %s() in %s: %s\n",
            line, func, file, msg);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#include "core_blas_s.h"
#include "core_blas_d.h"
#include "core_blas_c.h"
#include "core_blas_z.h"

#endif // ICL_CORE_BLAS_H
