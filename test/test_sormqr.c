/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zunmqr.c, normal z -> s, Mon Feb  6 14:07:17 2017
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
#include <math.h>

#include <omp.h>

#define REAL

/***************************************************************************//**
 *
 * @brief Tests SORMQR.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_sormqr(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_SIDE);
            print_usage(PARAM_TRANS);
            print_usage(PARAM_DIM);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
            print_usage(PARAM_IB);
            print_usage(PARAM_HMODE);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "side",
                     InfoSpacing, "trans",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "K",
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
             "%*c %*c %*d %*d %*d %*d %*d %*d %*d %*c",
             InfoSpacing, param[PARAM_SIDE].c,
             InfoSpacing, param[PARAM_TRANS].c,
             InfoSpacing, param[PARAM_DIM].dim.m,
             InfoSpacing, param[PARAM_DIM].dim.n,
             InfoSpacing, param[PARAM_DIM].dim.k,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_IB].i,
             InfoSpacing, param[PARAM_HMODE].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t trans = plasma_trans_const(param[PARAM_TRANS].c);
    plasma_enum_t side  = plasma_side_const(param[PARAM_SIDE].c);

    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;

    // Number of Householder reflectors to use.
    int k = param[PARAM_DIM].dim.k;

    // Dimensions of matrix A differ for different combinations of
    // side and trans.
    int am, an;
    if (side == PlasmaLeft) {
        am = m;
        if (trans == PlasmaNoTrans) {
            an = m;
        }
        else {
            an = k;
        }
    }
    else {
        am = n;
        if (trans == PlasmaNoTrans) {
            an = k;
        }
        else {
            an = n;
        }
    }
    int lda = imax(1, am + param[PARAM_PADA].i);

    // Dimensions of matrix B.
    int bm = m;
    int bn = n;
    int ldb = imax(1, bm  + param[PARAM_PADB].i);

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

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
    // Allocate and initialize array A for construction of matrix Q as
    // A = Q*R.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*an*sizeof(float));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*an, A);
    assert(retval == 0);

    //================================================================
    // Prepare factorization of matrix A.
    //================================================================
    plasma_desc_t T;
    plasma_sgeqrf(am, an, A, lda, &T);

    //================================================================
    // Prepare m-by-n matrix B.
    //================================================================
    float *B =
        (float*)malloc((size_t)ldb*bn*
                                            sizeof(float));
    assert(B != NULL);

    retval = LAPACKE_slarnv(1, seed, (size_t)ldb*bn, B);
    assert(retval == 0);

    float *Bref = NULL;
    if (test) {
        // Store the original array if residual is to be evaluated.
        Bref = (float*)malloc((size_t)ldb*bn*
                                           sizeof(float));
        assert(Bref != NULL);

        memcpy(Bref, B, (size_t)ldb*bn*sizeof(float));
    }

    //================================================================
    // Prepare explicit matrix Q.
    //================================================================
    // Number of Householder reflectors to be used depends on
    // side and trans combination.
    int qk = an;

    int qm, qn, ldq;
    float *Q = NULL;
    if (test) {
        qm  = am;
        qn  = an;
        ldq = qm;
        Q = (float *)malloc((size_t)ldq*qn*
                                         sizeof(float));
        // Build explicit Q.
        plasma_sorgqr(qm, qn, qk, A, lda, T, Q, ldq);
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_sormqr(side, trans,
                  bm, bn, qk,
                  A, lda, T,
                  B, ldb);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_sormqr(side, bm, bn, qk) /
                            time / 1e9;

    //================================================================
    // Test results by comparing implicit and explicit actions of Q.
    //================================================================
    if (test) {
        // Set dimensions of the resulting matrix C = op(Q)*B or C = B*op(Q).
        int cm, cn;
        if (side == PlasmaLeft) {
            cn = bn;
            if (trans == PlasmaNoTrans) {
                cm = qm;
            }
            else {
                cm = qn;
            }
        }
        else {
            cm = bm;
            if (trans == PlasmaNoTrans) {
                cn = qn;
            }
            else {
                cn = qm;
            }
        }

        // Apply explicit Q and compute the difference. For example, for
        // PlasmaLeft and PlasmaNoTrans, B <- implicit(Q)*B - Q*Bref.
        if (side == PlasmaLeft) {
            plasma_sgemm(trans, PlasmaNoTrans,
                         cm, cn, m,
                         -1.0, Q, ldq,
                               Bref, ldb,
                          1.0, B, ldb);
        }
        else {
            plasma_sgemm(PlasmaNoTrans, trans,
                         cm, cn, n,
                         -1.0, Bref, ldb,
                               Q, ldq,
                          1.0, B, ldb);
        }

        // |B|_1
        float work[1];
        float normB = LAPACKE_slange_work(LAPACK_COL_MAJOR, '1', bm, bn,
                                           Bref, ldb, work);

        // Compute error in the difference.
        // |implicit(Q)*B - Q*Bref|_1
        float error = LAPACKE_slange_work(LAPACK_COL_MAJOR, '1', cm, cn,
                                           B, ldb, work);

        // Normalize the result.
        error /= (cm * normB);

        // Store the results.
        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = (error < tol);
    }


    //================================================================
    // Free arrays.
    //================================================================
    plasma_desc_destroy(&T);
    free(A);
    free(B);
    if (test) {
        free(Bref);
        free(Q);
    }
}
