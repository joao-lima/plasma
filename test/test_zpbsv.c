/**
 *
 * @file test_zpbtrf.c
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define COMPLEX

#define A(i,j)  A[i + j*lda]

/***************************************************************************//**
 *
 * @brief Tests ZPBSV.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_zpbsv(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            //  pbtrf params
            print_usage(PARAM_UPLO);
            print_usage(PARAM_N);
            print_usage(PARAM_KL);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
            //  gbtrs params for check
            print_usage(PARAM_NRHS);
            print_usage(PARAM_PADB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                "%*s %*s %*s %*s %*s %*s %*s ",
                InfoSpacing, "UpLo",
                InfoSpacing, "N",
                InfoSpacing, "KD",
                InfoSpacing, "PadA",
                InfoSpacing, "NB",
                InfoSpacing, "NRHS",
                InfoSpacing, "PadB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
        "%*c %*d %*d %*d %*d %*d %*d",
        InfoSpacing, param[PARAM_UPLO].c,
        InfoSpacing, param[PARAM_N].i,
        InfoSpacing, param[PARAM_KL].i,
        InfoSpacing, param[PARAM_PADA].i,
        InfoSpacing, param[PARAM_NB].i,
        InfoSpacing, param[PARAM_NRHS].i,
        InfoSpacing, param[PARAM_PADB].i);

    //================================================================
    // Set parameters.
    //================================================================
    PLASMA_enum uplo = PLASMA_uplo_const(param[PARAM_UPLO].c);
    int pada = param[PARAM_PADA].i;
    int n    = param[PARAM_N].i;
    int lda  = imax(1, n + pada);

    int kd   = (uplo == PlasmaUpper ? param[PARAM_KU].i : param[PARAM_KL].i);

    int test = param[PARAM_TEST].c == 'y';
    double tol = param[PARAM_TOL].d * LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    PLASMA_Set(PLASMA_TILE_SIZE, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================

    // band matrix A in full storage (also used for solution check)
    PLASMA_Complex64_t *A =
        (PLASMA_Complex64_t*)malloc((size_t)lda*n*sizeof(PLASMA_Complex64_t));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_zlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);
    // make it SPD
    int i,j;
    for (i = 0; i < n; ++i) {
        A(i,i) = creal(A(i,i)) + n;
        for (j = 0; j < i; ++j) {
            A(j,i) = conj(A(i,j));
        }
    }
    // zero out elements outside the band
    PLASMA_Complex64_t zzero = 0.0;
    for (i = 0; i < n; i++) {
        for (j = i+kd+1; j < n; j++) A(i, j) = zzero;
    }
    for (j = 0; j < n; j++) {
        for (i = j+kd+1; i < n; i++) A(i, j) = zzero;
    }

    // band matrix A in LAPACK storage
    int ldab = imax(1, kd+1+pada);
    PLASMA_Complex64_t *AB = NULL;
    AB = (PLASMA_Complex64_t*)malloc((size_t)ldab*n*sizeof(PLASMA_Complex64_t));
    assert(AB != NULL);

    // convert into LAPACK's symmetric band storage
    for (j = 0; j < n; j++) {
        for (i = 0; i < ldab; i++) AB[i + j*ldab] = zzero;
        if (uplo == PlasmaUpper) {
            for (i = imax(0, j-kd); i <= j; i++) AB[i-j+kd + j*ldab] = A(i, j);
        } else {
            for (i = j; i <= imin(n-1, j+kd); i++) AB[i-j + j*ldab] = A(i, j);
        }
    }

    // set up solution vector X with right-hand-side B
    int nrhs = param[PARAM_NRHS].i;
    int ldx = imax(1, n + param[PARAM_PADB].i);
    PLASMA_Complex64_t *X = NULL;
    X = (PLASMA_Complex64_t*)malloc((size_t)ldx*nrhs*sizeof(PLASMA_Complex64_t));
    assert(X != NULL);

    retval = LAPACKE_zlarnv(1, seed, (size_t)ldx*nrhs, X);
    assert(retval == 0);
    
    // copy B to X
    int ldb = ldx;
    PLASMA_Complex64_t *B = NULL;
    if (test) {
        B = (PLASMA_Complex64_t*)malloc((size_t)ldb*nrhs*sizeof(PLASMA_Complex64_t));
        assert(B != NULL);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'F', n, nrhs, X, ldx, B, ldb);
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    int iinfo;

    plasma_time_t start = omp_get_wtime();
    iinfo = PLASMA_zpbsv(uplo, n, kd, nrhs, AB, ldab, X, ldx);
    if (iinfo != 0) printf( " zpbsv failed with info=%d\n",iinfo );
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = 0.0 / time / 1e9;

    //================================================================
    // Test results by computing residual norm.
    //================================================================
    if (test) {
        PLASMA_Complex64_t zone =   1.0;
        PLASMA_Complex64_t zmone = -1.0;

        // compute residual vector
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, nrhs, n,
                    CBLAS_SADDR(zmone), A, lda,
                                        X, ldx,
                    CBLAS_SADDR(zone),  B, ldb);

        // compute various norms
        double *work = NULL;
        work = (double*)malloc((size_t)n*sizeof(double));
        assert(work != NULL);

        double Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'F', n, n,    A, lda, work);
        double Xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', n, nrhs, X, ldb, work);
        double Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', n, nrhs, B, ldb, work);
        double residual = Rnorm/(n*Anorm*Xnorm);

        param[PARAM_ERROR].d = residual;
        param[PARAM_SUCCESS].i = residual < tol;

        // free arrays
        free(work);
        free(B);
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(X);
    free(AB);
    free(A);
}
