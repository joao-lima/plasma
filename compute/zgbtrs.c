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
 * @ingroup plasma_gbtrs
 *
 *  Solves a system of linear equations A * X = B with a Hermitian positive definite
 *  matrix A using the Cholesky factorization of A (i.e., A = L*L^T or A = U^T*U)
 *  computed by plasma_zpbtrf.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of A is stored;
 *          - PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in] kd
 *          The number of suuperdiagonals within the band of A if uplo is upper,
 *          or the number of suuperdiagonals if uplo is lower. kd >= 0.
 *
 * @param[in] nrhs
 *          The number of right hand sides, i.e., the number of
 *          columns of the matrix B. nrhs >= 0.
 *
 * @param[in,out] AB
 *          The triangular factor U or L from the Cholesky
 *          factorization A = U^H*U or A = L*L^H, computed by
 *          plasma_zpotrf.
 *          Remark: If out-of-place layout translation is used, the
 *          matrix A can be considered as input, however if inplace
 *          layout translation is enabled, the content of A will be
 *          reordered for computation and restored before exiting the
 *          function.
 *
 * @param[in] ldab
 *          The leading dimension of the array AB.
 *
 * @param[in,out] B
 *          On entry, the n-by-nrhs right hand side matrix B.
 *          On exit, if return value = 0, the n-by-nrhs solution matrix X.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,n).
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval  < 0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa plasma_omp_zgbtrs
 * @sa plasma_cgbtrs
 * @sa plasma_dgbtrs
 * @sa plasma_sgbtrs
 * @sa plasma_zpbtrf
 *
 ******************************************************************************/
int plasma_zgbtrs(int m, int n, int kl, int ku, int nrhs,
                  plasma_complex64_t *pAB, int ldab,
                  int *IPIV,
                  plasma_complex64_t *pB,  int ldb)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }

    // Check input arguments.
    if (m < 0) {
        plasma_error("illegal value of m");
        return -1;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -2;
    }
    if (kl < 0) {
        plasma_error("illegal value of kd");
        return -3;
    }
    if (ku < 0) {
        plasma_error("illegal value of ku");
        return -4;
    }
    if (nrhs < 0) {
        plasma_error("illegal value of nrhs");
        return -5;
    }
    if (ldab < imax(1, 1+kl+ku)) {
        plasma_error("illegal value of ldab");
        return -7;
    }
    if (ldb < imax(1, n)) {
        plasma_error("illegal value of ldb");
        return -10;
    }

    // quick return
    if (imax(n, nrhs) == 0)
        return PlasmaSuccess;

    // Set tiling parameters.
    int nb = plasma->nb;

    // Initialize tile matrix descriptors.
    plasma_desc_t AB;
    plasma_desc_t B;
    int tku = (ku+kl+nb-1)/nb; // number of tiles in upper band (not including diagonal)
    int tkl = (kl+nb-1)/nb;    // number of tiles in lower band (not including diagonal)
    int lm  = (tku+tkl+1)*nb;  // since we use zgetrf on panel, we pivot back within panel.
                               // this could fill the last tile of the panel,
                               // and we need extra NB space on the bottom
    int retval;
    retval = plasma_desc_general_band_create(PlasmaComplexDouble, PlasmaGeneral,
                                             nb, nb, lm, n, 0, 0, m, n, kl, ku,
                                             &AB);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_desc_general_band_create() failed");
        return retval;
    }
    retval = plasma_desc_general_create(PlasmaComplexDouble, nb, nb,
                                        m, nrhs, 0, 0, m, nrhs, &B);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        plasma_desc_destroy(&AB);
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
        plasma_omp_zpb2desc(pAB, ldab, AB, sequence, &request);
        plasma_omp_zge2desc(pB, ldb, B, sequence, &request);

        // Call the tile async function.
        plasma_omp_zgbtrs(AB, IPIV, B, sequence, &request);

        // Translate back to LAPACK layout.
        plasma_omp_zdesc2ge(B, pB, ldb, sequence, &request);
    }
    // implicit synchronization

    // Free matrix A in tile layout.
    plasma_desc_destroy(&AB);
    plasma_desc_destroy(&B);

    // Return status.
    int status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup plasma_gbtrs
 *
 *  Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Non-blocking tile version of plasma_zgbtrs().
 *  May return before the computation is finished.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of A is stored;
 *          - PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] AB
 *          The triangular factor U or L from the Cholesky factorization
 *          A = U^H*U or A = L*L^H, computed by plasma_zpotrf.
 *
 * @param[in,out] B
 *          On entry, the n-by-nrhs right hand side matrix B.
 *          On exit, if return value = 0, the n-by-nrhs solution matrix X.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).  Check
 *          the sequence->status for errors.
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
 * @sa plasma_zgbtrs
 * @sa plasma_omp_zgbtrs
 * @sa plasma_omp_cgbtrs
 * @sa plasma_omp_dgbtrs
 * @sa plasma_omp_sgbtrs
 * @sa plasma_omp_zgbtrf
 *
 ******************************************************************************/
void plasma_omp_zgbtrs(plasma_desc_t AB, int *IPIV, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // Check input arguments.
    if (plasma_desc_check(AB) != PlasmaSuccess) {
        plasma_error("invalid A");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (plasma_desc_check(B) != PlasmaSuccess) {
        plasma_error("invalid B");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (sequence == NULL) {
        plasma_fatal_error("NULL sequence");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (request == NULL) {
        plasma_fatal_error("NULL request");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // quick return
    if (AB.n == 0 || B.n == 0)
        return;

    // Call the parallel functions.
    plasma_pztbsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans,
                  PlasmaUnit,
                  1.0, AB,
                       B,
                  IPIV,
                  sequence, request);
    plasma_pztbsm(PlasmaLeft, PlasmaUpper, PlasmaNoTrans,
                  PlasmaNonUnit,
                  1.0, AB,
                       B,
                  IPIV,
                  sequence, request);
}