/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlacpy.c, normal z -> d, Sat Feb  4 20:58:54 2017
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_lacpy
 *
 *  Copies all or part of a two-dimensional matrix A to another matrix B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaGeneral: entire A,
 *          - PlasmaUpper:   upper triangle,
 *          - PlasmaLower:   lower triangle.
 *
 * @param[in] m
 *          The number of rows of the matrices A and B.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices A and B.
 *          n >= 0.
 *
 * @param[in] A
 *          The m-by-n matrix to copy.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          The m-by-n copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          ldb >= max(1,m).
 *
 ******************************************************************************/
void core_dlacpy(plasma_enum_t uplo,
                 int m, int n,
                 const double *A, int lda,
                       double *B, int ldb)
{
    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,
                        lapack_const(uplo),
                        m, n,
                        A, lda,
                        B, ldb);
}

/******************************************************************************/
void core_omp_dlacpy(plasma_enum_t uplo,
                     int m, int n,
                     const double *A, int lda,
                           double *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
#if defined(USE_OMPEXT)
omp_set_task_name("dlacpy");
#endif
    #pragma omp task depend(in:A[0:lda*n]) \
                     depend(out:B[0:ldb*n])
    {
        if (sequence->status == PlasmaSuccess)
            core_dlacpy(uplo,
                        m, n,
                        A, lda,
                        B, ldb);
    }
}
