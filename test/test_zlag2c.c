/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions mixed zc -> ds
 *
 **/

#include "test.h"
#include "flops.h"
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>

/***************************************************************************//**
 *
 * @brief Tests ZLAG2C
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_zlag2c(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s",
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
             "%*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters
    //================================================================
    int m = param[PARAM_M].i;
    int n = param[PARAM_N].i;

    int lda  = imax(1, m + param[PARAM_PADA].i);
    int ldas = imax(1, m + param[PARAM_PADB].i);

    int    test = param[PARAM_TEST].c == 'y';
    double eps  = LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    plasma_complex64_t *A =
        (plasma_complex64_t*)malloc((size_t)lda*n*sizeof(plasma_complex64_t));
    assert(A != NULL);

    plasma_complex32_t *As =
        (plasma_complex32_t*)malloc((size_t)ldas*n*sizeof(plasma_complex32_t));
    assert(As != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_zlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    plasma_complex32_t *AsRef = NULL;
    if (test) {
        AsRef = (plasma_complex32_t*)malloc(
            (size_t)ldas*n*sizeof(plasma_complex32_t));
        assert(AsRef != NULL);
    }

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();

    PLASMA_zlag2c(m, n, A, lda, As, ldas);

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d   = time;
    param[PARAM_GFLOPS].d = NAN;

    //================================================================
    // Test results by comparing to result of core_zlag2c function
    //================================================================
    if (test) {
        // Calculate relative error |As_ref - As|_F / |As_ref|_F < 3*eps
        // Using 3*eps covers complex arithmetic

        core_zlag2c(m, n, A, lda, AsRef, ldas);

        float work[1];

        // Calculate Frobenius norm of reference result As_ref
        double AsNormRef = LAPACKE_clange_work(
                               LAPACK_COL_MAJOR, 'F', ldas, n, AsRef, ldas, work);

        // Calculate difference As_ref-As
        plasma_complex32_t cmone = -1.0;
        cblas_caxpy((size_t)ldas*n, CBLAS_SADDR(cmone), As, 1, AsRef, 1);

        // Calculate Frobenius norm of As_ref-As
        double AsNormDiff = LAPACKE_clange_work(
                               LAPACK_COL_MAJOR, 'F', ldas, n, AsRef, ldas, work);

        // Calculate relative error |As_ref-As|_F / |As_ref|_F
        double error = AsNormDiff/AsNormRef;

        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A);
    free(As);

    if (test)
        free(AsRef);
}
