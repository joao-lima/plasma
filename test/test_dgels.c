/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgels.c, normal z -> d, Mon Feb  6 14:07:04 2017
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

#define REAL

/***************************************************************************//**
 *
 * @brief Tests DGELS.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_dgels(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_DIM);
            print_usage(PARAM_NRHS);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
            print_usage(PARAM_IB);
            print_usage(PARAM_HMODE);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "NRHS",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadB",
                     InfoSpacing, "NB",
                     InfoSpacing, "IB",
                     InfoSpacing, "Hous. mode");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d %*d %*d %*d %*c",
             InfoSpacing, param[PARAM_DIM].dim.m,
             InfoSpacing, param[PARAM_DIM].dim.n,
             InfoSpacing, param[PARAM_NRHS].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_IB].i,
             InfoSpacing, param[PARAM_HMODE].c);

    //================================================================
    // Set parameters.
    //================================================================
    int m    = param[PARAM_DIM].dim.m;
    int n    = param[PARAM_DIM].dim.n;
    int nrhs = param[PARAM_NRHS].i;

    int lda = imax(1, m + param[PARAM_PADA].i);
    int ldb = imax(1, imax(m, n) + param[PARAM_PADB].i);

    int test = param[PARAM_TEST].c == 'y';
    double tol = param[PARAM_TOL].d * LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    if (param[PARAM_HMODE].c == 't') {
        plasma_set(PlasmaHouseholderMode, PlasmaTreeHouseholder);
    }
    else {
        plasma_set(PlasmaHouseholderMode, PlasmaFlatHouseholder);
    }

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    double *A =
        (double*)malloc((size_t)lda*n*sizeof(double));
    assert(A != NULL);

    double *B =
        (double*)malloc((size_t)ldb*nrhs*
                                    sizeof(double));
    assert(B != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    retval = LAPACKE_dlarnv(1, seed, (size_t)ldb*nrhs, B);
    assert(retval == 0);

    // store the original arrays if residual is to be evaluated
    double *Aref = NULL;
    double *Bref = NULL;
    if (test) {
        Aref = (double*)malloc(
            (size_t)lda*n*sizeof(double));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(double));

        Bref = (double*)malloc(
            (size_t)ldb*nrhs*sizeof(double));
        assert(Bref != NULL);

        memcpy(Bref, B, (size_t)ldb*nrhs*sizeof(double));
    }

    //================================================================
    // Prepare the descriptor for matrix T.
    //================================================================
    plasma_desc_t T;

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_dgels(PlasmaNoTrans, m, n, nrhs,
                 A, lda,
                 &T,
                 B, ldb);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;

    double flop;
    if (m >= n) {
        // cost of QR-based factorization and solve
        flop = flops_dgeqrf(m, n) + flops_dgeqrs(m, n, nrhs);
    }
    else {
        // cost of LQ-based factorization, triangular solve, and Q^T application
        flop = flops_dgelqf(m, n) +
               flops_dtrsm(PlasmaLeft, m, nrhs) +
               flops_dormlq(PlasmaLeft, n, nrhs, m);
    }
    param[PARAM_GFLOPS].d = flop / time / 1e9;

    //================================================================
    // Test results by solving a linear system.
    //================================================================
    if (test) {
        // |A|_F
        double work[1];
        double Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n,
                                           Aref, lda, work);

        // |B|_F
        double Bnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, nrhs,
                                           Bref, ldb, work);

        // |X|_F, solution X is now stored in the n-by-nrhs part of B
        double Xnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', n, nrhs,
                                           B, ldb, work);

        // compute residual and store it in B = A*X - B
        double zone  =  1.0;
        double zmone = -1.0;
        double zzero =  0.0;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, nrhs, n,
                    (zone), Aref, lda, B, ldb,
                    (zmone), Bref, ldb);

        // Compute B = A^T * (A*X - B)
        cblas_dgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, nrhs, m,
                    (zone), Aref, lda, Bref, ldb,
                    (zzero), B, ldb);

        // |RES|_F
        double Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', n, nrhs,
                                           B, ldb, work);

        // normalize the result
        double result = Rnorm / ((Anorm*Xnorm+Bnorm)*n);

        param[PARAM_ERROR].d = result;
        param[PARAM_SUCCESS].i = result < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    plasma_desc_destroy(&T);
    free(A);
    free(B);
    if (test) {
        free(Aref);
        free(Bref);
    }
}
