/**
 *
 * @file zgemm.c
 *
 *  PLASMA computational routine.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver.
 *
 * @version 3.0.0
 * @author Jakub Kurzak
 * @date 2016-01-01
 * @precisions normal z -> s d c
 *
 **/

#include "../control/async.h"
#include "../control/context.h"
#include "../control/descriptor.h"
#include "../control/internal.h"
#include "../include/plasma_z.h"
#include "../include/plasmatypes.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  Performs one of the matrix-matrix operations
 *
 *          \f[ C = \alpha [op( A )\times op( B )] + \beta C, \f]
 *
 *  where op( X ) is one of:
 *          - op( X ) = X  or
 *          - op( X ) = X' or
 *          - op( X ) = conjg( X' ),
 *
 *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          - PlasmaNoTrans:   B is not transposed,
 *          - PlasmaTrans:     B is transposed,
 *          - PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrix op( A ) and of the matrix C.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix op( B ) and of the matrix C.
 *          n >= 0.
 *
 * @param[in] k
 *          The number of columns of the matrix op( A ) and the number of rows
 *          of the matrix op( B ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix, where ka is k when transA = PlasmaNoTrans,
 *          and is m otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[in] B
 *          An ldb-by-kb matrix, where kb is n when transB = PlasmaNoTrans,
 *          and is k otherwise.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix. On exit, the array is overwritten by the m-by-n
 *          matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 *******************************************************************************
 *
 * @retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgemm_Tile
 * @sa PLASMA_zgemm_Tile_Async
 * @sa PLASMA_cgemm
 * @sa PLASMA_dgemm
 * @sa PLASMA_sgemm
 *
 ******************************************************************************/
int PLASMA_zgemm(PLASMA_enum transA, PLASMA_enum transB,
                 int m, int n, int k,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int lda,
                                           PLASMA_Complex64_t *B, int ldb,
                 PLASMA_Complex64_t beta,  PLASMA_Complex64_t *C, int ldc)
{
    int Am, An;
    int Bm, Bn;
    int nb;
    int retval;
    int status;

    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descC;

    PLASMA_Complex64_t zzero = (PLASMA_Complex64_t)0.0;
    PLASMA_Complex64_t zone = (PLASMA_Complex64_t)1.0;

    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    // Check input arguments.
    if ((transA != PlasmaNoTrans) &&
        (transA != PlasmaTrans) &&
        (transA != PlasmaConjTrans)) {
        plasma_error("illegal value of transA");
        return -1;
    }
    if ((transB != PlasmaNoTrans) &&
        (transB != PlasmaTrans) &&
        (transB != PlasmaConjTrans)) {
        plasma_error("illegal value of transB");
        return -2;
    }
    if (m < 0) {
        plasma_error("illegal value of m");
        return -3;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -4;
    }
    if (k < 0) {
        plasma_error("illegal value of k");
        return -5;
    }
    if (A == NULL) {
        plasma_error("NULL A");
        return -7;
    }

    if (transA == PlasmaNoTrans) {
        Am = m;
        An = k;
    }
    else {
        Am = k;
        An = m;
    }
    if (transB == PlasmaNoTrans) {
        Bm = k;
        Bn = n;
    }
    else {
        Bm = n;
        Bn = k;
    }

    if (lda < imax(1, Am)) {
        plasma_error("illegal value of lda");
        return -8;
    }
    if (B == NULL) {
        plasma_error("NULL B");
        return -9;
    }
    if (ldb < imax(1, Bm)) {
        plasma_error("illegal value of ldb");
        return -10;
    }
    if (C == NULL) {
        plasma_error("NULL C");
        return -12;
    }
    if (ldc < imax(1, m)) {
        plasma_error("illegal value of ldc");
        return -13;
    }

    // quick return
    if (m == 0 || n == 0 || ((alpha == zzero || k == 0) && beta == zone))
        return PLASMA_SUCCESS;

    // Tune.
    // if (plasma_tune(PLASMA_FUNC_ZGEMM, m, n, 0) != PLASMA_SUCCESS) {
    //     plasma_error("plasma_tune() failed");
    //     return status;
    // }
    nb = plasma->nb;

    // Initialize tile matrix descriptors.
    descA = plasma_desc_init(PlasmaComplexDouble, nb, nb,
                             nb*nb, Am, An, 0, 0, Am, An);

    descB = plasma_desc_init(PlasmaComplexDouble, nb, nb,
                             nb*nb, Bm, Bn, 0, 0, Bm, Bn);

    descC = plasma_desc_init(PlasmaComplexDouble, nb, nb,
                             nb*nb, Cm, An, 0, 0, Am, An);

    // Allocate matrices in tile layout.
    retval = plasma_desc_mat_alloc(&descA);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("plasma_desc_mat_alloc() failed");
        return retval;
    }

    retval = plasma_desc_mat_alloc(&descB);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("plasma_desc_mat_alloc() failed");
        plasma_desc_mat_free(&descA);
        return retval;
    }

    retval = plasma_desc_mat_alloc(&descC);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("plasma_desc_mat_alloc() failed");
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
        return retval;
    }

    // Create sequence.
    PLASMA_sequence *sequence = NULL;
    retval = plasma_sequence_create(&sequence);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("plasma_sequence_create() failed");
        return retval;
    }
    // Initialize request.
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    // Translate to tile layout.
    retavl = PLASMA_zcm2ccrb_Async(A, lda, &descA, sequence, &request);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcm2ccrb_Async() failed");
        return retval;
    }
    retval = PLASMA_zcm2ccrb_Async(B, ldb, &descB, sequence, &request);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcm2ccrb_Async() failed");
        return retval;
    }
    retval = PLASMA_zcm2ccrb_Async(C, ldc, &descC, sequence, &request);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcm2ccrb_Async() failed");
        return retval;
    }

    // Call the tile async function.
    retval = PLASMA_zgemm_Tile_Async(transA, transB,
                                     alpha, &descA,
                                            &descB,
                                      beta, &descC,
                                     sequence, &request);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgemm_Tile_Async() failed");
        return retval;
    }

    // Translate back to LAPACK layout.
    retval = PLASMA_zccrb2cm_Async(C, ldc, &descC, sequence, &request);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zccrb2cm_Async() failed");
        return retval;
    }

    // Free matrices in tile layout.
    plasma_desc_mat_free(&descA);
    plasma_desc_mat_free(&descB);
    plasma_desc_mat_free(&descC);

    // Return status.
    status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  Performs matrix multiplication.
 *  Tile equivalent of PLASMA_zgemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of matrix A.
 *
 * @param[in] B
 *          Descriptor of matrix B.
 *
 * @param[in,out] C
 *          Descriptor of matrix C.
 *
 *******************************************************************************
 *
 * @retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgemm
 * @sa PLASMA_zgemm_Tile_Async
 * @sa PLASMA_cgemm_Tile
 * @sa PLASMA_dgemm_Tile
 * @sa PLASMA_sgemm_Tile
 *
 ******************************************************************************/
int PLASMA_zgemm_Tile(PLASMA_enum transA, PLASMA_enum transB,
                      PLASMA_Complex64_t alpha, PLASMA_desc *A,
                                                PLASMA_desc *B,
                      PLASMA_Complex64_t beta,  PLASMA_desc *C)
{
    int retval;
    int status;

    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    // Create sequence.
    PLASMA_sequence *sequence = NULL;
    retval = plasma_sequence_create(&sequence);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("plasma_sequence_create() failed");
        return retval;
    }
    // Initialize request.
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    // Call the tile async function.
    retval = PLASMA_zgemm_Tile_Async(transA, transB,
                                     alpha, A,
                                            B,
                                      beta, C,
                                     sequence, &request);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgemm_Tile_Async() failed");
        return retval;
    }

    // Return status.
    status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  Performs matrix multiplication.
 *  Non-blocking equivalent of PLASMA_zgemm_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgemm
 * @sa PLASMA_zgemm_Tile
 * @sa PLASMA_cgemm_Tile_Async
 * @sa PLASMA_dgemm_Tile_Async
 * @sa PLASMA_sgemm_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB,
                            PLASMA_Complex64_t alpha, PLASMA_desc *A,
                                                      PLASMA_desc *B,
                            PLASMA_Complex64_t beta,  PLASMA_desc *C,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    // Check input arguments.
    if ((transA != PlasmaNoTrans) &&
        (transA != PlasmaTrans) &&
        (transA != PlasmaConjTrans)) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("illegal value of transA");
        return -1;
    }
    if ((transB != PlasmaNoTrans) &&
        (transB != PlasmaTrans) &&
        (transB != PlasmaConjTrans)) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("illegal value of transB");
        return -2;
    }
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("invalid A");
        return -4;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("invalid B");
        return -5;
    }
    if (plasma_desc_check(C) != PLASMA_SUCCESS) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("invalid C");
        return -7;
    }
    if (sequence == NULL) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("NULL sequence");
        return -8;
    }
    if (request == NULL) {
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        plasma_error("NULL request");
        return -9;
    }

    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    if (transA == PlasmaNoTrans) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
        Anb = A->nb;
        Ai  = A->i;
        Aj  = A->j;
    }
    else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
        Anb = A->mb;
        Ai  = A->j;
        Aj  = A->i;
    }
    if (transB == PlasmaNoTrans) {
        Bm  = B->m;
        Bn  = B->n;
        Bmb = B->mb;
        Bnb = B->nb;
        Bi  = B->i;
        Bj  = B->j;
    }
    else {
        Bm  = B->n;
        Bn  = B->m;
        Bmb = B->nb;
        Bnb = B->mb;
        Bi  = B->j;
        Bj  = B->i;
    }

    if (Amb != C->mb || Anb != Bmb || Bnb != C->nb) {
        plasma_error("tile size mismatch");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (Am != C->m || An != Bm || Bn != C->n) {
        plasma_error("matrix size mismatch");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (Ai%Amb != C->i%C->mb ||
        Bj%Bnb != C->j%C->nb || Aj%Anb != Bi%Bmb) {
        plasma_error("start indexes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    // Check sequence status.
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request,
                                   PLASMA_ERR_SEQUENCE_FLUSHED);

    // quick return
    if (C->m == 0 || C->n == 0 || An == 0 ||
        (alpha == (PLASMA_Complex64_t)0.0 && beta == (PLASMA_Complex64_t)1.0))
        return PLASMA_SUCCESS;

    // Call the parallel function.
    plasma_pzgemm(transA, transB,
                  alpha, A,
                         B,
                   beta, C,
                  sequence, request);

    return PLASMA_SUCCESS;
}
