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
int plasma_zgeadd(plasma_enum_t transa,
                  int m, int n,
                  plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                  plasma_complex64_t beta,  plasma_complex64_t *pB, int ldb);

int plasma_zgelqf(int m, int n,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T);

int plasma_zgelqs(int m, int n, int nrhs,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex64_t *pB, int ldb);

int plasma_zgels(plasma_enum_t trans,
                 int m, int n, int nrhs,
                 plasma_complex64_t *pA, int lda,
                 plasma_desc_t T,
                 plasma_complex64_t *pB, int ldb);

int plasma_zgemm(plasma_enum_t transa, plasma_enum_t transb,
                 int m, int n, int k,
                 plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                           plasma_complex64_t *pB, int ldb,
                 plasma_complex64_t beta,  plasma_complex64_t *pC, int ldc);

int plasma_zgeqrf(int m, int n,
                  plasma_complex64_t *f77A, int lda,
                  plasma_desc_t T);

int plasma_zgeqrs(int m, int n, int nrhs,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex64_t *pB, int ldb);

int plasma_zhemm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                           plasma_complex64_t *pB, int ldb,
                 plasma_complex64_t beta,  plasma_complex64_t *pC, int ldc);

int plasma_zher2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                            plasma_complex64_t *pB, int ldb,
                  double beta,              plasma_complex64_t *pC, int ldc);

int plasma_zherk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 double alpha, plasma_complex64_t *pA, int lda,
                 double beta,  plasma_complex64_t *pC, int ldc);

int plasma_zlacpy(plasma_enum_t uplo,
                  int m, int n,
                  plasma_complex64_t *pA, int lda,
                  plasma_complex64_t *pB, int ldb);

int plasma_zpbsv(plasma_enum_t uplo,
                 int n, int kd, int nrhs,
                 plasma_complex64_t *pAB, int ldab,
                 plasma_complex64_t *pB,  int ldb);

int plasma_zpbtrs(plasma_enum_t uplo,
                  int n, int kd, int nrhs,
                  plasma_complex64_t *pAB, int ldab,
                  plasma_complex64_t *pB,  int ldb);

int plasma_zpbtrf(plasma_enum_t uplo,
                  int n, int kd,
                  plasma_complex64_t *pAB, int ldab);

int plasma_zposv(plasma_enum_t uplo,
                 int n, int nrhs,
                 plasma_complex64_t *pA, int lda,
                 plasma_complex64_t *pB, int ldb);

int plasma_zpotrf(plasma_enum_t uplo,
                  int n,
                  plasma_complex64_t *pA, int lda);

int plasma_zpotrs(plasma_enum_t uplo,
                  int n, int nrhs,
                  plasma_complex64_t *pA, int lda,
                  plasma_complex64_t *pB, int ldb);

int plasma_zsymm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                           plasma_complex64_t *pB, int ldb,
                 plasma_complex64_t beta,  plasma_complex64_t *pC, int ldc);

int plasma_zsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                            plasma_complex64_t *pB, int ldb,
                  plasma_complex64_t beta,  plasma_complex64_t *pC, int ldc);

int plasma_zsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                 plasma_complex64_t beta,  plasma_complex64_t *pC, int ldc);

int plasma_ztradd(plasma_enum_t uplo, plasma_enum_t transa,
                  int m, int n,
                  plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                  plasma_complex64_t beta,  plasma_complex64_t *pB, int ldb);

int plasma_ztrmm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                           plasma_complex64_t *pB, int ldb);

int plasma_ztrsm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 plasma_complex64_t alpha, plasma_complex64_t *pA, int lda,
                                           plasma_complex64_t *pB, int ldb);

int plasma_zunglq(int m, int n, int k,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex64_t *Qf77, int ldq);

int plasma_zungqr(int m, int n, int k,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex64_t *Qf77, int ldq);

int plasma_zunmlq(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex64_t *pC, int ldc);

int plasma_zunmqr(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  plasma_complex64_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex64_t *pC, int ldc);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/
void plasma_omp_zgeadd(plasma_enum_t transa,
                       plasma_complex64_t alpha, plasma_desc_t A,
                       plasma_complex64_t beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_zdesc2ge(plasma_desc_t A,
                         plasma_complex64_t *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_zdesc2pb(plasma_desc_t A,
                         plasma_complex64_t *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_zge2desc(plasma_complex64_t *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_zge2desc(plasma_complex64_t *pA, int lda,
                           plasma_desc_t A,
                           plasma_sequence_t *sequence,
                           plasma_request_t *request);

void plasma_omp_zgelqf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgelqs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgels(plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T,
                      plasma_desc_t B, plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgemm(plasma_enum_t transa, plasma_enum_t transb,
                      plasma_complex64_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_complex64_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgeqrf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zgeqrs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zhemm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_complex64_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_complex64_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zher2k(plasma_enum_t uplo, plasma_enum_t trans,
                       plasma_complex64_t alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       double beta,              plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zherk(plasma_enum_t uplo, plasma_enum_t trans,
                      double alpha, plasma_desc_t A,
                      double beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zlacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpb2desc(plasma_complex64_t *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_zpbsv(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpbtrf(plasma_enum_t uplo, plasma_desc_t AB,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpbtrs(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zposv(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpotrf(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zpotrs(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zsymm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_complex64_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_complex64_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                       plasma_complex64_t alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       plasma_complex64_t beta,  plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                      plasma_complex64_t alpha, plasma_desc_t A,
                      plasma_complex64_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ztradd(plasma_enum_t uplo, plasma_enum_t transa,
                       plasma_complex64_t alpha, plasma_desc_t A,
                       plasma_complex64_t beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_ztrmm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      plasma_complex64_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ztrsm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      plasma_complex64_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zunglq(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zungqr(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zunmlq(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_zunmqr(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_Z_H
