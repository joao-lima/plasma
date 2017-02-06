/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgbsv.c, normal z -> c, Mon Feb  6 14:07:03 2017
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define A(m, n) (plasma_complex32_t*)plasma_tile_addr(A, m, n)

#define COMPLEX

/***************************************************************************//**
 *
 * @brief Tests CGBSV.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_cgbsv(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_DIM);
            print_usage(PARAM_KL);
            print_usage(PARAM_KU);
            print_usage(PARAM_NRHS);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
            print_usage(PARAM_IB);
            print_usage(PARAM_NTPF);
            print_usage(PARAM_ZEROCOL);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "N",
                     InfoSpacing, "KL",
                     InfoSpacing, "KU",
                     InfoSpacing, "NRHS",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB",
                     InfoSpacing, "IB",
                     InfoSpacing, "NTPF",
                     InfoSpacing, "ZeroCol");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d %*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_DIM].dim.n,
             InfoSpacing, param[PARAM_KL].i,
             InfoSpacing, param[PARAM_KU].i,
             InfoSpacing, param[PARAM_KU].i,
             InfoSpacing, param[PARAM_NRHS].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_IB].i,
             InfoSpacing, param[PARAM_NTPF].i,
             InfoSpacing, param[PARAM_ZEROCOL].i);

    //================================================================
    // Set parameters.
    //================================================================
    int n    = param[PARAM_DIM].dim.n;
    int kl   = param[PARAM_KL].i;
    int ku   = param[PARAM_KU].i;
    int lda  = imax(1, n+param[PARAM_PADA].i);
    int nrhs = param[PARAM_NRHS].i;
    int ldb  = imax(1, n + param[PARAM_PADB].i);
    int ldx  = ldb;

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    plasma_set(PlasmaNumPanelThreads, param[PARAM_NTPF].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    plasma_complex32_t *A = (plasma_complex32_t*)malloc(
        (size_t)lda*n*sizeof(plasma_complex32_t));
    assert(A != NULL);
    plasma_complex32_t *X = (plasma_complex32_t*)malloc(
        (size_t)ldx*nrhs*sizeof(plasma_complex32_t));
    assert(X != NULL);
    int *ipiv = (int*)malloc((size_t)n*sizeof(int));
    assert(ipiv != NULL);

    // set up right-hand-sides X
    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_clarnv(1, seed, (size_t)ldb*nrhs, X);
    assert(retval == 0);
    // copy X to B for test
    plasma_complex32_t *B = NULL;
    if (test) {
        B = (plasma_complex32_t*)malloc(
            (size_t)ldb*nrhs*sizeof(plasma_complex32_t));
        assert(B != NULL);
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'F', n, nrhs, X, ldx, B, ldb);
    }

    // set up matrix A
    retval = LAPACKE_clarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);
    // zero out elements outside the band
    for (int i = 0; i < n; i++) {
        for (int j = i+ku+1; j < n; j++) A[i + j*lda] = 0.0;
    }
    for (int j = 0; j < n; j++) {
        for (int i = j+kl+1; i < n; i++) A[i + j*lda] = 0.0;
    }

    int zerocol = param[PARAM_ZEROCOL].i;
    if (zerocol >= 0 && zerocol < n)
        memset(&A[zerocol*lda], 0, n*sizeof(plasma_complex32_t));

    // save A for test
    plasma_complex32_t *Aref = NULL;
    if (test) {
        Aref = (plasma_complex32_t*)malloc(
            (size_t)lda*n*sizeof(plasma_complex32_t));
        assert(Aref != NULL);
        memcpy(Aref, A, (size_t)lda*n*sizeof(plasma_complex32_t));
    }

    int nb = param[PARAM_NB].i;
    // band matrix A in skewed LAPACK storage
    int kut  = (ku+kl+nb-1)/nb; // # of tiles in upper band (not including diagonal)
    int klt  = (kl+nb-1)/nb;    // # of tiles in lower band (not including diagonal)
    int ldab = (kut+klt+1)*nb;  // since we use cgetrf on panel, we pivot back within panel.
                                // this could fill the last tile of the panel,
                                // and we need extra NB space on the bottom
    plasma_complex32_t *AB = NULL;
    AB = (plasma_complex32_t*)malloc((size_t)ldab*n*sizeof(plasma_complex32_t));
    assert(AB != NULL);
    // convert into LAPACK's skewed storage
    for (int j = 0; j < n; j++) {
        int i_kl = imax(0,   j-ku);
        int i_ku = imin(n-1, j+kl);
        for (int i = 0; i < ldab; i++)
            AB[i + j*ldab] = 0.0;
        for (int i = i_kl; i <= i_ku; i++)
            AB[kl + i-(j-ku) + j*ldab] = A[i + j*lda];
    }
    plasma_complex32_t *ABref = NULL;
    if (test) {
        ABref = (plasma_complex32_t*)malloc(
            (size_t)ldab*n*sizeof(plasma_complex32_t));
        assert(ABref != NULL);

        memcpy(ABref, AB, (size_t)ldab*n*sizeof(plasma_complex32_t));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    int plainfo = plasma_cgbsv(n, kl, ku, nrhs, AB, ldab, ipiv, X, ldx);

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = 0.0;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        if (plainfo == 0) {
            // compute residual vector
            plasma_complex32_t zone  =  1.0;
            plasma_complex32_t zmone = -1.0;
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, nrhs, n,
                        CBLAS_SADDR(zmone), Aref, lda,
                                            X, ldx,
                        CBLAS_SADDR(zone),  B, ldb);

            // compute various norms
            float *work = NULL;
            work = (float*)malloc((size_t)n*sizeof(float));
            assert(work != NULL);

            float Anorm = LAPACKE_clange_work(
                LAPACK_COL_MAJOR, 'F', n, n,    A, lda, work);
            float Xnorm = LAPACKE_clange_work(
                LAPACK_COL_MAJOR, 'I', n, nrhs, X, ldb, work);
            float Rnorm = LAPACKE_clange_work(
                LAPACK_COL_MAJOR, 'I', n, nrhs, B, ldb, work);
            float residual = Rnorm/(n*Anorm*Xnorm);

            param[PARAM_ERROR].d = residual;
            param[PARAM_SUCCESS].i = residual < tol;

            // free workspaces
            free(work);
        }
        else {
            int lapinfo = LAPACKE_cgbsv(
                              LAPACK_COL_MAJOR,
                              n, kl, ku, nrhs, ABref, ldab, ipiv, X, ldx);
            if (plainfo == lapinfo) {
                param[PARAM_ERROR].d = 0.0;
                param[PARAM_SUCCESS].i = 1;
            }
            else {
                param[PARAM_ERROR].d = INFINITY;
                param[PARAM_SUCCESS].i = 0;
            }
        }
        free(B);
        free(Aref);
        free(ABref);
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(AB);
    free(ipiv);
    free(X);
}
