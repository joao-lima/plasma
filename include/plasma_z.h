/**
 *
 * @file
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester, Univ. of California Berkeley and
 *  Univ. of Colorado Denver.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef ICL_PLASMA_Z_H
#define ICL_PLASMA_Z_H

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Standard interface.
 **/
int PLASMA_zgelqf(int m, int n,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT);

int PLASMA_zgelqs(int m, int n, int nrhs,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT,
                  plasma_complex64_t *B, int ldb);

int PLASMA_zgels(PLASMA_enum trans, int m, int n, int nrhs,
                 plasma_complex64_t *A, int lda,
                 plasma_desc_t *descT,
                 plasma_complex64_t *B, int ldb);

int PLASMA_zgemm(PLASMA_enum transA, PLASMA_enum transB,
                 int m, int n, int k,
                 plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                                           plasma_complex64_t *B, int ldb,
                 plasma_complex64_t beta,  plasma_complex64_t *C, int ldc);

int PLASMA_zgeqrf(int m, int n,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT);

int PLASMA_zgeqrs(int m, int n, int nrhs,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT,
                  plasma_complex64_t *B, int ldb);

int PLASMA_zhemm(PLASMA_enum side, PLASMA_enum uplo, int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                                           plasma_complex64_t *B, int ldb,
                 plasma_complex64_t beta,  plasma_complex64_t *C, int ldc);

int PLASMA_zher2k(PLASMA_enum uplo, PLASMA_enum trans,
                  int n, int k,
                  plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                                            plasma_complex64_t *B, int ldb,
                               double beta, plasma_complex64_t *C, int ldc);

int PLASMA_zherk(PLASMA_enum uplo, PLASMA_enum trans,
                 int n, int k,
                 double alpha, plasma_complex64_t *A, int lda,
                 double beta,  plasma_complex64_t *C, int ldc);

int PLASMA_zpbsv(PLASMA_enum uplo, int n, int kd, int nrhs,
                 plasma_complex64_t *AB, int ldab,
                 plasma_complex64_t *B, int ldb);

int PLASMA_zpbtrs(PLASMA_enum uplo, int n, int kd, int nrhs,
                  plasma_complex64_t *AB, int ldab,
                  plasma_complex64_t *B, int ldb);

int PLASMA_zpbtrf(PLASMA_enum uplo,
                  int n, int kd,
                  plasma_complex64_t *AB, int ldab);

int PLASMA_zposv(PLASMA_enum uplo, int n, int nrhs,
                 plasma_complex64_t *A, int lda,
                 plasma_complex64_t *B, int ldb);

int PLASMA_zpotrf(PLASMA_enum uplo,
                  int n,
                  plasma_complex64_t *A, int lda);

int PLASMA_zpotrs(PLASMA_enum uplo,
                  int n, int nrhs,
                  plasma_complex64_t *A, int lda,
                  plasma_complex64_t *B, int ldb);

int PLASMA_zsymm(PLASMA_enum side, PLASMA_enum uplo, int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                                           plasma_complex64_t *B, int ldb,
                 plasma_complex64_t beta,  plasma_complex64_t *C, int ldc);

int PLASMA_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                  int n, int k,
                  plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                                            plasma_complex64_t *B, int ldb,
                  plasma_complex64_t beta,  plasma_complex64_t *C, int ldc);

int PLASMA_zsyrk(PLASMA_enum uplo, PLASMA_enum trans, int n, int k,
                 plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                 plasma_complex64_t beta,  plasma_complex64_t *C, int ldc);

int PLASMA_ztradd(PLASMA_enum uplo, PLASMA_enum transA, int m, int n,
                  plasma_complex64_t  alpha,
                  plasma_complex64_t *A, int lda,
                  plasma_complex64_t  beta,
                  plasma_complex64_t *B, int ldb);

int PLASMA_ztrmm(PLASMA_enum side, PLASMA_enum uplo,
                 PLASMA_enum transA, PLASMA_enum diag,
                 int m, int n, plasma_complex64_t alpha,
                 plasma_complex64_t *A, int lda,
                 plasma_complex64_t *B, int ldb);

int PLASMA_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                 PLASMA_enum transA, PLASMA_enum diag,
                 int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *A, int lda,
                                           plasma_complex64_t *B, int ldb);

int PLASMA_zunglq(int m, int n, int k,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT,
                  plasma_complex64_t *Q, int ldq);

int PLASMA_zungqr(int m, int n, int k,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT,
                  plasma_complex64_t *Q, int ldq);

int PLASMA_zunmlq(PLASMA_enum side, PLASMA_enum trans, int m, int n, int k,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT,
                  plasma_complex64_t *C, int ldc);

int PLASMA_zunmqr(PLASMA_enum side, PLASMA_enum trans, int m, int n, int k,
                  plasma_complex64_t *A, int lda,
                  plasma_desc_t *descT,
                  plasma_complex64_t *C, int ldc);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/

void PLASMA_zccrb2cm_Async(plasma_desc_t *A, plasma_complex64_t *Af77, int lda,
                           plasma_sequence_t *sequence,
                           plasma_request_t *request);

void PLASMA_zccrb2cm_band_Async(PLASMA_enum uplo,
                                plasma_desc_t *A,
                                plasma_complex64_t *Af77, int lda,
                                plasma_sequence_t *sequence,
                                plasma_request_t *request);

void PLASMA_zcm2ccrb_Async(plasma_complex64_t *Af77, int lda,
                           plasma_desc_t *A,
                           plasma_sequence_t *sequence,
                           plasma_request_t *request);

void PLASMA_zcm2ccrb_band_Async(PLASMA_enum uplo,
                                plasma_complex64_t *Af77, int lda,
                                plasma_desc_t *A,
                                plasma_sequence_t *sequence,
                                plasma_request_t *request);

void plasma_omp_zgelqf(plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgelqs(plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_desc_t *descB, plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgels(PLASMA_enum trans,
                      plasma_desc_t *descA, plasma_desc_t *descT,
                      plasma_desc_t *descB, plasma_workspace_t *work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgemm(PLASMA_enum transA, PLASMA_enum transB,
                      plasma_complex64_t alpha, plasma_desc_t *A,
                                                plasma_desc_t *B,
                      plasma_complex64_t beta,  plasma_desc_t *C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgeqrf(plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgeqrs(plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_desc_t *descB, plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zhemm(PLASMA_enum side, PLASMA_enum uplo,
                      plasma_complex64_t alpha, plasma_desc_t *A,
                                                plasma_desc_t *B,
                      plasma_complex64_t beta,  plasma_desc_t *C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zher2k(PLASMA_enum uplo, PLASMA_enum trans,
                       plasma_complex64_t alpha, plasma_desc_t *A,
                                                 plasma_desc_t *B,
                       double beta,              plasma_desc_t *C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zherk(PLASMA_enum uplo, PLASMA_enum trans,
                      double alpha, plasma_desc_t *A,
                      double beta,  plasma_desc_t *C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpbsv(PLASMA_enum uplo,
                      plasma_desc_t *AB,
                      plasma_desc_t *B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpbtrf(PLASMA_enum uplo, plasma_desc_t *AB,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpbtrs(PLASMA_enum uplo, plasma_desc_t *AB, plasma_desc_t *B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zposv(PLASMA_enum uplo, plasma_desc_t *A, plasma_desc_t *B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpotrf(PLASMA_enum uplo, plasma_desc_t *A,
                       plasma_sequence_t *sequence,
                       plasma_request_t *request);

void plasma_omp_zpotrs(PLASMA_enum uplo, plasma_desc_t *A, plasma_desc_t *B,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zsymm(PLASMA_enum side, PLASMA_enum uplo,
                      plasma_complex64_t alpha, plasma_desc_t *A,
                                                plasma_desc_t *B,
                      plasma_complex64_t beta,  plasma_desc_t *C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                       plasma_complex64_t alpha, plasma_desc_t *A,
                                                  plasma_desc_t *B,
                       plasma_complex64_t beta,  plasma_desc_t *C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zsyrk(PLASMA_enum uplo, PLASMA_enum trans,
                      plasma_complex64_t alpha, plasma_desc_t *A,
                      plasma_complex64_t beta,  plasma_desc_t *C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ztradd(PLASMA_enum uplo, PLASMA_enum transA,
                       plasma_complex64_t alpha, plasma_desc_t *A,
                       plasma_complex64_t beta,  plasma_desc_t *B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_ztrmm(PLASMA_enum side, PLASMA_enum uplo,
                      PLASMA_enum transA, PLASMA_enum diag,
                      plasma_complex64_t alpha, plasma_desc_t *A,
                                                plasma_desc_t *B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                      PLASMA_enum transA, PLASMA_enum diag,
                      plasma_complex64_t alpha, plasma_desc_t *A,
                                                plasma_desc_t *B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zunglq(plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_desc_t *descQ, plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zungqr(plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_desc_t *descQ, plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zunmlq(PLASMA_enum side, PLASMA_enum trans,
                       plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_desc_t *descC, plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zunmqr(PLASMA_enum side, PLASMA_enum trans,
                       plasma_desc_t *descA, plasma_desc_t *descT,
                       plasma_desc_t *descC, plasma_workspace_t *work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_Z_H
