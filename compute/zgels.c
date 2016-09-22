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

#include "plasma.h"
#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

/***************************************************************************//**
 *
 * @ingroup plasma_gels
 *
 *  Solves overdetermined or underdetermined linear systems
 *  involving an m-by-n matrix A using a QR or LQ factorization of A.  It
 *  is assumed that A has full rank.  The following options are provided:
 *
 *  # trans = PlasmaNoTrans and m >= n: find the least squares solution of an
 *    overdetermined system, i.e., solve the least squares problem:
 *    minimize || B - A*X ||.
 *
 *  # trans = PlasmaNoTrans and m < n: find the minimum norm solution of an
 *    underdetermined system A * X = B.
 *
 *  Several right-hand side vectors B and solution vectors X can be handled in a
 *  single call; they are stored as the columns of the m-by-nrhs right-hand side
 *  matrix B and the n-by-nrhs solution matrix X.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          - PlasmaNoTrans:  the linear system involves A
 *                            (the only supported option for now).
 *
 * @param[in] m
 *          The number of rows of the matrix A. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0.
 *
 * @param[in] nrhs
 *          The number of right hand sides, i.e., the number of columns of the
 *          matrices B and X.  nrhs >= 0.
 *
 * @param[in,out] A
 *          On entry, the m-by-n matrix A.
 *          On exit,
 *          if m >= n, A is overwritten by details of its QR factorization as
 *                     returned by PLASMA_zgeqrf;
 *          if m < n, A is overwritten by details of its LQ factorization as
 *                      returned by PLASMA_zgelqf.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[out] descT
 *          On exit, auxiliary factorization data.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs matrix B of right-hand side vectors, stored
 *          columnwise;
 *          On exit, if return value = 0, B is overwritten by the solution
 *          vectors, stored columnwise:
 *          if m >= n, rows 1 to N of B contain the least squares solution
 *          vectors; the residual sum of squares for the solution in each column
 *          is given by the sum of squares of the modulus of elements n+1 to m
 *          in that column;
 *          if m < n, rows 1 to n of B contain the minimum norm solution
 *          vectors;
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m,n).
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa plasma_omp_zgels
 * @sa PLASMA_cgels
 * @sa PLASMA_dgels
 * @sa PLASMA_sgels
 * @sa PLASMA_zgeqrf
 * @sa PLASMA_zgeqrs
 *
 ******************************************************************************/
int PLASMA_zgels(plasma_enum_t trans, int m, int n, int nrhs,
                 plasma_complex64_t *A, int lda,
                 plasma_desc_t *descT,
                 plasma_complex64_t *B, int ldb)
{
    int ib, nb;
    int retval;
    int status;

    plasma_desc_t descA;
    plasma_desc_t descB;

    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }

    // Check input arguments
    if (trans != PlasmaNoTrans) {
        plasma_error("only PlasmaNoTrans supported");
        return PlasmaErrorNotSupported;
    }
    if (m < 0) {
        plasma_error("illegal value of m");
        return -2;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -3;
    }
    if (nrhs < 0) {
        plasma_error("illegal value of nrhs");
        return -4;
    }
    if (lda < imax(1, m)) {
        plasma_error("illegal value of lda");
        return -6;
    }
    if (ldb < imax(1, imax(m, n))) {
        plasma_error("illegal value of ldb");
        return -9;
    }
    // Quick return
    if (imin(m, imin(n, nrhs)) == 0) {
        for (int i = 0; i < imax(m, n); i++)
            for (int j = 0; j < nrhs; j++)
                B[j*ldb+i] = 0.0;
        return PlasmaSuccess;
    }

    // Tune NB & IB depending on M, N & NRHS; Set NBNB
    //status = plasma_tune(PLASMA_FUNC_ZGELS, M, N, NRHS);
    //if (status != PlasmaSuccess) {
    //    plasma_error("plasma_tune() failed");
    //    return status;
    //}
    ib = plasma->ib;
    nb = plasma->nb;

    // Create tile matrices.
    retval = plasma_desc_general_create(PlasmaComplexDouble, nb, nb,
                                        lda, n, 0, 0, m, n, &descA);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        return retval;
    }
    retval = plasma_desc_general_create(PlasmaComplexDouble, nb, nb,
                                        ldb, nrhs, 0, 0, imax(m,n), nrhs,
                                        &descB);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        plasma_desc_destroy(&descA);
        return retval;
    }

    // Allocate workspace.
    plasma_workspace_t work;
    size_t lwork = nb + ib*nb;  // geqrt/gelqt: tau + work
    retval = plasma_workspace_alloc(&work, lwork, PlasmaComplexDouble);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_workspace_alloc() failed");
        return retval;
    }

    // Create sequence.
    plasma_sequence_t *sequence = NULL;
    retval = plasma_sequence_create(&sequence);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_sequence_create() failed");
        return retval;
    }

    // Initialize request.
    plasma_request_t request = PlasmaRequestInitializer;

    // asynchronous block
    #pragma omp parallel
    #pragma omp master
    {
        // Translate to tile layout.
        PLASMA_zcm2ccrb_Async(A, lda, &descA, sequence, &request);
        PLASMA_zcm2ccrb_Async(B, ldb, &descB, sequence, &request);

        // Call the tile async function.
        plasma_omp_zgels(PlasmaNoTrans, &descA, descT, &descB,
                         &work, sequence, &request);

        // Translate back to LAPACK layout.
        PLASMA_zccrb2cm_Async(&descA, A, lda, sequence, &request);
        PLASMA_zccrb2cm_Async(&descB, B, ldb, sequence, &request);
    }
    // implicit synchronization

    plasma_workspace_free(&work);

    // Free matrices in tile layout.
    plasma_desc_destroy(&descA);
    plasma_desc_destroy(&descB);

    // Return status.
    status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup plasma_gels
 *
 *  Solves overdetermined or underdetermined linear
 *  system of equations using the tile QR or the tile LQ factorization.
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          - PlasmaNoTrans:  the linear system involves A
 *                            (the only supported option for now).
 *
 * @param[in,out] A
 *          Descriptor of matrix A stored in the tile layout.
 *          On exit,
 *          if m >= n, A is overwritten by details of its QR factorization
 *                     as returned by PLASMA_zgeqrf;
 *          if m < n,  A is overwritten by details of its LQ factorization
 *                     as returned by PLASMA_zgelqf.
 *
 * @param[out] T
 *          Descriptor of matrix T.
 *          Auxiliary factorization data, computed by
 *          PLASMA_zgeqrf or PLASMA_zgelqf.
 *
 * @param[in,out] B
 *          Descriptor of matrix B.
 *          On entry, right-hand side matrix B in the tile layout.
 *          On exit, solution matrix X in the tile layout.
 *
 * @param[in] work
 *          Workspace for the auxiliary arrays needed by some coreblas kernels.
 *          For QR/LQ factorizations used in GELS, it contains preallocated
 *          space for TAU and WORK arrays.
 *          Allocated by the plasma_workspace_alloc function.
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
 *          request->status should never be set to PlasmaSuccess (the
 *          initial values) since another async call may be setting a
 *          failure value at the same time.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgels
 * @sa plasma_omp_cgels
 * @sa plasma_omp_dgels
 * @sa plasma_omp_sgels
 *
 ******************************************************************************/
void plasma_omp_zgels(plasma_enum_t trans, plasma_desc_t *A,
                      plasma_desc_t *T, plasma_desc_t *B,
                      plasma_workspace_t *work,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // Check input arguments
    if (trans != PlasmaNoTrans) {
        plasma_error("only PlasmaNoTrans supported");
        plasma_request_fail(sequence, request, PlasmaErrorNotSupported);
        return;
    }
    if (plasma_desc_check(A) != PlasmaSuccess) {
        plasma_error("invalid descriptor A");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (A->mb != plasma->nb || A->nb != plasma->nb) {
        plasma_error("wrong tile dimensions of A");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (plasma_desc_check(T) != PlasmaSuccess) {
        plasma_error("invalid descriptor T");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (T->mb != plasma->ib || T->nb != plasma->nb) {
        plasma_error("wrong tile dimensions of T");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (plasma_desc_check(B) != PlasmaSuccess) {
        plasma_error("invalid descriptor B");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (B->mb != plasma->nb || B->nb != plasma->nb) {
        plasma_error("wrong tile dimensions of B");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (sequence == NULL) {
        plasma_error("NULL sequence");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (request == NULL) {
        plasma_error("NULL request");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // Check sequence status.
    if (sequence->status != PlasmaSuccess) {
        plasma_request_fail(sequence, request, PlasmaErrorSequence);
        return;
    }

    // Quick return
    // (imin(m, imin(n, nrhs)) == 0)
    if (imin(A->m, imin(A->n, B->n)) == 0) {
        // zero matrix B
        plasma_pzlaset(PlasmaGeneral, 0., 0., *B, sequence, request);
        return;
    }

    if (A->m >= A->n) {
        // solution based on QR factorization of A
        plasma_pzgeqrf(*A, *T, work, sequence, request);

        // Plasma_ConjTrans will be converted to PlasmaTrans by the
        // automatic datatype conversion, which is what we want here.
        // Note that PlasmaConjTrans is protected from this conversion.
        plasma_pzunmqr(PlasmaLeft, Plasma_ConjTrans,
                       *A, *B, *T,
                       work, sequence, request);

        plasma_pztrsm(PlasmaLeft, PlasmaUpper,
                      PlasmaNoTrans, PlasmaNonUnit,
                      1.0,
                      plasma_desc_view(*A, 0, 0, A->n, A->n),
                      plasma_desc_view(*B, 0, 0, A->n, B->n),
                      sequence, request);
    }
    else {
        // solution based on LQ factorization of A
        plasma_pzgelqf(*A, *T, work, sequence, request);

        // zero the trailing block of the right-hand side matrix
        // (B has less rows than X)
        plasma_pzlaset(PlasmaGeneral, 0., 0.,
                       plasma_desc_view(*B, A->m, 0,
                                        A->n - A->m, B->n),
                       sequence, request);

        // Solve L * Y = B
        plasma_complex64_t zone  =  1.0;
        plasma_pztrsm(
            PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit,
            zone, plasma_desc_view(*A, 0, 0, A->m, A->m),
                  plasma_desc_view(*B, 0, 0, A->m, B->n),
            sequence, request);

        // Find X = Q^H * Y
        // Plasma_ConjTrans will be converted to PlasmaTrans by the
        // automatic datatype conversion, which is what we want here.
        // Note that PlasmaConjTrans is protected from this conversion.
        plasma_pzunmlq(PlasmaLeft, Plasma_ConjTrans,
                       *A, *B, *T,
                       work, sequence, request);
    }
}
