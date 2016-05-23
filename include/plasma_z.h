/**
 *
 * @file plasma_z.h
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley, Univ. of Colorado Denver and
 *  Univ. of Manchester.
 *
 * @version 3.0.0
 * @author Jakub Kurzak
 * @date 2016-01-01
 * @precisions normal z -> s d c
 *
 **/
#ifndef ICL_PLASMA_Z_H
#define ICL_PLASMA_Z_H

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Standard interface.
 **/
int PLASMA_zgemm(
    PLASMA_enum transA, PLASMA_enum transB,
    int m, int n, int k,
    PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int lda,
                              PLASMA_Complex64_t *B, int ldb,
    PLASMA_Complex64_t beta,  PLASMA_Complex64_t *C, int ldc);

int PLASMA_zsyrk(
    PLASMA_enum uplo, PLASMA_enum trans, 
    int n, int k,
    PLASMA_Complex64_t alpha, 
    PLASMA_Complex64_t *A, int lda,
    PLASMA_Complex64_t beta,  
    PLASMA_Complex64_t *C, int ldc);

int PLASMA_zherk(
    PLASMA_enum uplo, PLASMA_enum trans, 
    int n, int k, 
    double alpha, 
    PLASMA_Complex64_t *A, int lda, 
    double beta, 
    PLASMA_Complex64_t *C, int ldc);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/
void PLASMA_zgemm_Tile_Async(
    PLASMA_enum transA, PLASMA_enum transB,
    PLASMA_Complex64_t alpha, PLASMA_desc *descA,
                              PLASMA_desc *descB,
    PLASMA_Complex64_t beta,  PLASMA_desc *descC,
    PLASMA_sequence *sequence, PLASMA_request *request);

void PLASMA_zsyrk_Tile_Async(
     PLASMA_enum uplo, PLASMA_enum trans,
     PLASMA_Complex64_t alpha, PLASMA_desc *A,
     PLASMA_Complex64_t beta,  PLASMA_desc *C,
     PLASMA_sequence *sequence, PLASMA_request *request);

void PLASMA_zherk_Tile_Async(
    PLASMA_enum uplo, PLASMA_enum trans, 
    double alpha, PLASMA_desc *A, 
    double beta, PLASMA_desc *C, 
    PLASMA_sequence *sequence, PLASMA_request *request);

/***************************************************************************//**
 *  Layout translation async.
 **/
void PLASMA_zcm2ccrb_Async(
    PLASMA_Complex64_t *Af77, int lda,
    PLASMA_desc *A,
    PLASMA_sequence *sequence, PLASMA_request *request);

void PLASMA_zccrb2cm_Async(
    PLASMA_desc *A,
    PLASMA_Complex64_t *Af77, int lda,
    PLASMA_sequence *sequence, PLASMA_request *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_Z_H
