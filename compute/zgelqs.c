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

#include "plasma_types.h"
#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_z.h"

/***************************************************************************//**
 *
 * @ingroup plasma_gelqs
 *
 *  Computes a minimum-norm solution min | A*X - B | using the
 *  LQ factorization A = L*Q computed by PLASMA_zgelqf.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix A. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= m >= 0.
 *
 * @param[in] nrhs
 *          The number of columns of B. nrhs >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned
 *          by PLASMA_zgelqf.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= m.
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by PLASMA_zgelqf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= n.
 *
 *******************************************************************************
 *
 * @retval PLASMA_SUCCESS successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa plasma_omp_zgelqs
 * @sa PLASMA_cgelqs
 * @sa PLASMA_dgelqs
 * @sa PLASMA_sgelqs
 * @sa PLASMA_zgelqf
 *
 ******************************************************************************/
int PLASMA_zgelqs(int m, int n, int nrhs,
                  PLASMA_Complex64_t *A, int lda,
                  PLASMA_desc *descT,
                  PLASMA_Complex64_t *B, int ldb)
{
    int ib, nb;
    int retval;
    int status;

    PLASMA_desc descA, descB;

    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    // Check input arguments
    if (m < 0) {
        plasma_error("illegal value of m");
        return -1;
    }
    if (n < 0 || m > n) {
        plasma_error("illegal value of n");
        return -2;
    }
    if (nrhs < 0) {
        plasma_error("illegal value of nrhs");
        return -3;
    }
    if (lda < imax(1, m)) {
        plasma_error("illegal value of lda");
        return -5;
    }
    if (ldb < imax(1, imax(1, n))) {
        plasma_error("illegal value of ldb");
        return -8;
    }
    // Quick return
    if (m == 0 || n == 0 || nrhs == 0)
        return PLASMA_SUCCESS;

    // Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE
    //status = plasma_tune(PLASMA_FUNC_ZGELS, M, N, NRHS);
    //if (status != PLASMA_SUCCESS) {
    //    plasma_error("plasma_tune() failed");
    //    return status;
    //}
    ib = plasma->ib;
    nb = plasma->nb;

    // Initialize tile matrix descriptors.
    descA = plasma_desc_init(PlasmaComplexDouble, nb, nb,
                             nb*nb, lda, n, 0, 0, m, n);

    descB = plasma_desc_init(PlasmaComplexDouble, nb, nb,
                             nb*nb, ldb, nrhs, 0, 0, n, nrhs);

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

    // Allocate workspace.
    PLASMA_workspace work;
    size_t lwork = ib*nb;  // unmlq: work
    retval = plasma_workspace_alloc(&work, lwork, PlasmaComplexDouble);
    if (retval != PLASMA_SUCCESS) {
        plasma_error("plasma_workspace_alloc() failed");
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

    #pragma omp parallel
    #pragma omp master
    {
        // the Async functions are submitted here.  If an error occurs
        // (at submission time or at run time) the sequence->status
        // will be marked with an error.  After an error, the next
        // Async will not _insert_ more tasks into the runtime.  The
        // sequence->status can be checked after each call to _Async
        // or at the end of the parallel region.

        // Translate to tile layout.
        PLASMA_zcm2ccrb_Async(A, lda, &descA, sequence, &request);
        PLASMA_zcm2ccrb_Async(B, ldb, &descB, sequence, &request);

        // Call the tile async function.
        plasma_omp_zgelqs(&descA, descT, &descB, &work, sequence, &request);

        // Translate back to LAPACK layout.
        // It is not needed to translate the descriptor back
        // for out-of-place storage.
        //PLASMA_zccrb2cm_Async(&descA, A, lda, sequence, &request);
        PLASMA_zccrb2cm_Async(&descB, B, ldb, sequence, &request);
    }
    // implicit synchronization

    plasma_workspace_free(&work);

    // Free matrices in tile layout.
    plasma_desc_mat_free(&descA);
    plasma_desc_mat_free(&descB);

    // Return status.
    status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup plasma_gelqs
 *
 *  Computes a minimum-norm solution using previously computed LQ factorization.
 *  Non-blocking tile version of PLASMA_zgelqs().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of matrix A.
 *          A is stored in the tile layout.
 *
 * @param[in] T
 *          Descriptor of matrix T.
 *          Auxiliary factorization data, computed by PLASMA_zgelqf.
 *
 * @param[in,out] B
 *          Descriptor of matrix B.
 *          On entry, right-hand side matrix B in the tile layout.
 *          On exit, solution matrix X in the tile layout.
 *
 * @param[in] work
 *          Workspace for the auxiliary arrays needed by some coreblas kernels.
 *          For multiplication by Q contains preallocated space for WORK
 *          arrays. Allocated by the plasma_workspace_alloc function.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 * @retval void
 *          Errors are returned by setting sequence->status and
 *          request->status to error values.  The sequence->status and
 *          request->status should never be set to PLASMA_SUCCESS (the
 *          initial values) since another async call may be setting a
 *          failure value at the same time.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgelqs
 * @sa plasma_omp_cgelqs
 * @sa plasma_omp_dgelqs
 * @sa plasma_omp_sgelqs
 * @sa plasma_omp_zgelqf
 *
 ******************************************************************************/
void plasma_omp_zgelqs(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B,
                       PLASMA_workspace *work,
                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }

    // Check input arguments.
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("invalid descriptor A");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("invalid descriptor T");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("invalid descriptor B");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }
    if (A->nb != A->mb || B->nb != B->mb) {
        plasma_error("only square tiles supported");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }
    if (sequence == NULL) {
        plasma_error("NULL sequence");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }
    if (request == NULL) {
        plasma_error("NULL request");
        plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        return;
    }

    // Check sequence status.
    if (sequence->status != PLASMA_SUCCESS) {
        plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);
        return;
    }

    // Quick return
    // (   m == 0 ||    n == 0 || nrhs == 0)
    if (A->m == 0 || A->n == 0 || B->n == 0)
        return;

    // zero the trailing block of the right-hand side matrix
    // (B has less rows than X)
    plasma_pzlaset(PlasmaFull, 0., 0.,
                   plasma_desc_submatrix(*B, A->m, 0,
                                         A->n - A->m, B->n),
                   sequence, request);

    // Solve L * Y = B
    PLASMA_Complex64_t zone  =  1.0;
    plasma_pztrsm(
        PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit,
        zone, plasma_desc_submatrix(*A, 0, 0, A->m, A->m),
        plasma_desc_submatrix(*B, 0, 0, A->m, B->n),
        sequence, request);

    // Find X = Q^H * Y
    // Plasma_ConjTrans will be converted to PlasmaTrans by the
    // automatic datatype conversion, which is what we want here.
    // Note that PlasmaConjTrans is protected from this conversion.
    plasma_pzunmlq(PlasmaLeft, Plasma_ConjTrans,
                   *A, *B, *T,
                   work, sequence, request);
}
