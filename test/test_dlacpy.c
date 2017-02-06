/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlacpy.c, normal z -> d, Mon Feb  6 14:07:08 2017
 *
 **/

#include "core_blas.h"
#include "core_lapack.h"
#include "flops.h"
#include "plasma.h"
#include "test.h"

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***************************************************************************//**
 *
 * @brief Tests DLACPY
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_dlacpy(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info
            print_usage(PARAM_UPLO);
            print_usage(PARAM_DIM);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "UpLo",
                     InfoSpacing, "m",
                     InfoSpacing, "n",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadB",
                     InfoSpacing, "nb");
        }
        return;
    }
    // Return column values
    snprintf(info, InfoLen,
             "%*c %*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_UPLO].c,
             InfoSpacing, param[PARAM_DIM].dim.m,
             InfoSpacing, param[PARAM_DIM].dim.n,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters
    //================================================================
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);

    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;

    int lda = imax(1, m + param[PARAM_PADA].i);
    int ldb = imax(1, m + param[PARAM_PADB].i);

    int    test = param[PARAM_TEST].c == 'y';
    double eps  = LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    double *A =
        (double*)malloc((size_t)lda*n*sizeof(double));
    assert(A != NULL);

    double *B =
        (double*)malloc((size_t)ldb*n*sizeof(double));
    assert(B != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    lapack_int mtrxLayout = LAPACK_COL_MAJOR;

    // LAPACKE_[zcds]lantr_work has a bug (returns 0)
    // in MKL <= 11.3.3 (at least). Fixed in LAPACK 3.6.1.
    // For now, zero out the opposite triangle and use lange.
    // @sa test_dtrmm

    // Enforce zeroes in general rectangle or upper or lower triangle of B
    switch (uplo) {
        case PlasmaLower:
            LAPACKE_dlaset_work(
                mtrxLayout, 'U', m-1, n-1, 0.0, 0.0, &B[m], ldb);
            break;
        case PlasmaUpper:
            LAPACKE_dlaset_work(
                mtrxLayout, 'L', m-1, n-1, 0.0, 0.0, &B[1], ldb);
            break;
        default:
            LAPACKE_dlaset_work(
                mtrxLayout, 'G', m,   n,   0.0, 0.0,  B,    ldb);
    }

    double *Bref = NULL;
    if (test) {
        Bref = (double*)malloc(
            (size_t)ldb*n*sizeof(double));
        assert(Bref != NULL);

        memcpy(Bref, B, (size_t)ldb*n*sizeof(double));
    }

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();

    retval = plasma_dlacpy(uplo, m, n, A, lda, B, ldb);

    plasma_time_t stop = omp_get_wtime();

    param[PARAM_GFLOPS].d = 0.0;

    if (retval != PlasmaSuccess) {
        plasma_error("plasma_dlacpy() failed");
        param[PARAM_TIME].d    = 0.0;
        param[PARAM_ERROR].d   = 1.0;
        param[PARAM_SUCCESS].i = false;
        return;
    }
    else {
        param[PARAM_TIME].d = stop-start;
    }

    //================================================================
    // Test results by comparing to result of core_dlacpy function
    //================================================================
    if (test) {
        // Calculate relative error |B_ref - B|_F / |B_ref|_F < 3*eps
        // Using 3*eps covers complex arithmetic

        retval = LAPACKE_dlacpy_work(
                    mtrxLayout, lapack_const(uplo), m, n, A, lda, Bref, ldb);

        if (retval != PlasmaSuccess) {
            coreblas_error("LAPACKE_dlacpy_work() failed");
            param[PARAM_ERROR].d   = 1.0;
            param[PARAM_SUCCESS].i = false;
            return;
        }

        double work[1];

        // Calculate Frobenius norm of reference result B_ref
        double BnormRef  = LAPACKE_dlange_work(
                               mtrxLayout, 'F', m, n, Bref, ldb, work);

        // Calculate difference B_ref-B
        double zmone = -1.0;
        cblas_daxpy((size_t)ldb*n, (zmone), B, 1, Bref, 1);

        // Calculate Frobenius norm of B_ref-B
        double BnormDiff = LAPACKE_dlange_work(
                               mtrxLayout, 'F', m, n, Bref, ldb, work);

        // Calculate relative error |B_ref-B|_F / |B_ref|_F
        double error = BnormDiff/BnormRef;

        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A);
    free(B);

    if (test) {
        free(Bref);
    }
}
