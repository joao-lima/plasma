/**
 *
 * @file
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
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define A(m, n) (plasma_complex64_t*)plasma_tile_addr(A, m, n)

#define COMPLEX

/***************************************************************************//**
 *
 * @brief Tests ZGBTRF.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_zgbtrf(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_KL);
            print_usage(PARAM_KU);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
            print_usage(PARAM_IB);
            print_usage(PARAM_NTPF);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "KL",
                     InfoSpacing, "KU",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB",
                     InfoSpacing, "IB",
                     InfoSpacing, "NTPF");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_KL].i,
             InfoSpacing, param[PARAM_KU].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_IB].i,
             InfoSpacing, param[PARAM_NTPF].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_complex64_t zzero =  0.0;
    plasma_complex64_t zone  =  1.0;
    plasma_complex64_t zmone = -1.0;

    int m  = param[PARAM_M].i;
    int n  = param[PARAM_N].i;
    int kl = param[PARAM_KL].i;
    int ku = param[PARAM_KU].i;
    int lda = imax(1, m+param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    double tol = param[PARAM_TOL].d * LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    plasma_set(PlasmaNumPanelThreads, param[PARAM_NTPF].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    plasma_complex64_t *A =
        (plasma_complex64_t*)malloc((size_t)lda*n*sizeof(plasma_complex64_t));
    assert(A != NULL);

    int *ipiv = (int*)malloc((size_t)m*sizeof(int));
    assert(ipiv != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_zlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);
    // zero out elements outside the band
    for (int i = 0; i < m; i++) {
        for (int j = i+ku+1; j < n; j++) A[i + j*lda] = zzero;
    }
    for (int j = 0; j < n; j++) {
        for (int i = j+kl+1; i < m; i++) A[i + j*lda] = zzero;
    }

    // save A for test
    plasma_complex64_t *Aref = NULL;
    if (test) {
        Aref = (plasma_complex64_t*)malloc(
            (size_t)lda*n*sizeof(plasma_complex64_t));
        assert(Aref != NULL);
        memcpy(Aref, A, (size_t)lda*n*sizeof(plasma_complex64_t));
    }

    int nb = param[PARAM_NB].i;
    // band matrix A in skewed LAPACK storage
    int kut  = (ku+kl+nb-1)/nb; // # of tiles in upper band (not including diagonal)
    int klt  = (kl+nb-1)/nb;    // # of tiles in lower band (not including diagonal)
    int ldab = (kut+klt+1)*nb;  // since we use zgetrf on panel, we pivot back within panel.
                                // this could fill the last tile of the panel,
                                // and we need extra NB space on the bottom
    plasma_complex64_t *AB = NULL;
    AB = (plasma_complex64_t*)malloc((size_t)ldab*n*sizeof(plasma_complex64_t));
    assert(AB != NULL);
    // convert into LAPACK's skewed storage
    for (int j = 0; j < n; j++) {
        int i_kl = imax(0,   j-ku);
        int i_ku = imin(m-1, j+kl);
        for (int i = 0; i < ldab; i++) AB[i + j*ldab] = zzero;
        for (int i = i_kl; i <= i_ku; i++) AB[kl + i-(j-ku) + j*ldab] = A[i + j*lda];
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_zgbtrf(m, n, kl, ku, AB, ldab, ipiv);

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = 0.0;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        if (m == n) {
            // compute the residual norm ||A-bx||
            int nrhs = param[PARAM_NRHS].i;
            int ldb = imax(1, n + param[PARAM_PADB].i);

            // set up right-hand-side B
            plasma_complex64_t *B = NULL;
            B = (plasma_complex64_t*)malloc((size_t)ldb*nrhs*sizeof(plasma_complex64_t));
            assert(B != NULL);

            retval = LAPACKE_zlarnv(1, seed, (size_t)ldb*nrhs, B);
            assert(retval == 0);

            // copy B to X
            int ldx = ldb;
            plasma_complex64_t *X = NULL;
            X = (plasma_complex64_t*)malloc((size_t)ldx*nrhs*sizeof(plasma_complex64_t));
            assert(X != NULL);
            LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'F', n, nrhs, B, ldb, X, ldx);

            // solve for X
            int iinfo = plasma_zgbtrs(PlasmaNoTrans, n, kl, ku, nrhs, AB, ldab, ipiv, X, ldb);
            if (iinfo != 0) printf( " zpbtrs failed with info = %d\n", iinfo );

            // compute residual vector
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, nrhs, n,
                        CBLAS_SADDR(zmone), Aref, lda,
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

            // free workspaces
            free(work);
            free(X);
            free(B);
        }
        else {
            // compute the factorization error norm ||A-LU||
            plasma_complex64_t *LU = NULL, *work = NULL;
            double Anorm, Enorm = 0.0, temp;
            LU = (plasma_complex64_t*)malloc((size_t)n*lda*sizeof(plasma_complex64_t));
            work = (plasma_complex64_t*)malloc((size_t)m*sizeof(plasma_complex64_t));
            Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, '1', m, n, A, lda, &temp);
            for (int j = 1; j <= n; j++) {
                // compute L*U(:,j)
                int kd = kl + ku + 1;
                int ju = imin(kl+ku, j-1);
                int jl = imin(kl, m-j);
                int lenj = imin(m, j) - j + ju + 1;
                if (lenj > 0) {
                    int iw;
                    plasma_complex64_t alpha;
                    // reverse the piovot applied back within the panel
                    int jnb = imin(nb*(1+(j-1)/nb), imin(m,n));
                    for (int i = jnb; i > j; i--) {
                        iw = kd - (j-i);
                        alpha = AB[iw-1 + (j-1)*ldab];
                        int ip = ipiv[i-1];
                        if (i != ip) {
                            ip = kd - (j-ip);
                            AB[iw-1 + (j-1)*ldab] = AB[ip-1 + (j-1)*ldab];
                            AB[ip-1 + (j-1)*ldab] = alpha;
                        }
                    }
                    // compute L*U(:,j)
                    // copy U(:,j) into work, i.e., multiply with diagonals of L
                    cblas_zcopy(lenj, &AB[kd-ju-1 + (j-1)*ldab], 1, work, 1);
                    for (int i = lenj; i <= ju+jl; i++) {
                       work[i] = zzero;
                    }
                    // sum up U(i,j)*L(:,i)
                    for (int i = imin(m-1, j); i >= j-ju; i--) {
                        int il = imin(kl, m-i);
                        if (il > 0) {
                            iw = i - j + ju + 1;
                            alpha = work[iw-1];
                            cblas_zaxpy(il, CBLAS_SADDR(alpha), &AB[kd + (i-1)*ldab], 1, &work[iw], 1);
                            // revert the i-th pivot
                            int ip = ipiv[i-1];
                            if (i != ip) {
                                ip = ip - j + ju + 1;
                                work[iw-1] = work[ip-1];
                                work[ip-1] = alpha;
                            }
                        }
                    }
                    // subtract A(:,j), and compute 1-norm
                    cblas_zaxpy(ju+jl+1, CBLAS_SADDR(zmone), &A[(j-ju-1)+(j-1)*lda], 1, work, 1);
                    double Enormj = cblas_dzasum(ju+jl+1, work, 1);
                    if (Enormj > Enorm ) Enorm = Enormj;
                }
            }
            param[PARAM_ERROR].d = Enorm / (n*Anorm);
            param[PARAM_SUCCESS].i = (Enorm / (n*Anorm)) < tol;
            free(LU); free(work);
        }
        // free arrays
        free(Aref);
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(AB);
    free(ipiv);
}
