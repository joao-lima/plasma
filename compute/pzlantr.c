/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee,  US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m, n) (plasma_complex64_t*)plasma_tile_addr(A, m, n)

/***************************************************************************//**
 *  Parallel tile calculation of max, one, infinity or Frobenius matrix norm
 *  for a triangular matrix.
 ******************************************************************************/
void plasma_pzlantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A, double *work, double *value,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Check sequence status.
    if (sequence->status != PlasmaSuccess) {
        plasma_request_fail(sequence, request, PlasmaErrorSequence);
        return;
    }

    switch (norm) {
    double stub;
    double *workspace;
    double *scale;
    double *sumsq;
    //================
    // PlasmaMaxNorm
    //================
    case PlasmaMaxNorm:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            if (uplo == PlasmaLower) {
                for (int n = 0; n < imin(m, A.nt); n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zlange(PlasmaMaxNorm,
                                    mvam, nvan,
                                    A(m, n), ldam, 
                                    &stub, &work[A.mt*n+m],
                                    sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zlange(PlasmaMaxNorm,
                                    mvam, nvan,
                                    A(m, n), ldam,
                                    &stub, &work[A.mt*n+m],
                                    sequence, request);
                }
            }
            if (m < A.nt) {
                int nvam = plasma_tile_nview(A, m);
                core_omp_zlantr(PlasmaMaxNorm, uplo, diag,
                                mvam, nvam,
                                A(m, m), ldam, 
                                &stub, &work[A.mt*m+m],
                                sequence, request);
            }
        }
        #pragma omp taskwait
        core_omp_dlantr(PlasmaMaxNorm, uplo, PlasmaNonUnit,
                        A.mt, A.nt,
                        work, A.mt,
                        &stub, value,
                        sequence, request);
        break;
    //================
    // PlasmaOneNorm
    //================
    case PlasmaOneNorm:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            if (uplo == PlasmaLower) {
                for (int n = 0; n < imin(m, A.nt); n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zlange_aux(PlasmaOneNorm,
                                        mvam, nvan,
                                        A(m, n), ldam, 
                                        &work[A.n*m+n*A.nb],
                                        sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zlange_aux(PlasmaOneNorm,
                                        mvam, nvan,
                                        A(m, n), ldam, 
                                        &work[A.n*m+n*A.nb],
                                        sequence, request);
                }
            }
            if (m < A.nt) {
                int nvam = plasma_tile_nview(A, m);
                core_omp_zlantr_aux(PlasmaOneNorm, uplo, diag,
                                    mvam, nvam,
                                    A(m, m), ldam, 
                                    &work[A.n*m+m*A.nb],
                                    sequence, request);
            }

        }
        #pragma omp taskwait
        workspace = work + A.mt*A.n;
        core_omp_dlange(PlasmaInfNorm,
                        A.n, A.mt,
                        work, A.n,
                        workspace, value,
                        sequence, request);
        break;
    //================
    // PlasmaInfNorm
    //================
    case PlasmaInfNorm:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);

            if (uplo == PlasmaLower) {
                for (int n = 0; n < imin(m, A.nt); n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zlange_aux(PlasmaInfNorm,
                                        mvam, nvan,
                                        A(m, n), ldam, 
                                        &work[A.m*n+m*A.mb],
                                        sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zlange_aux(PlasmaInfNorm,
                                        mvam, nvan,
                                        A(m, n), ldam, 
                                        &work[A.m*n+m*A.mb],
                                        sequence, request);
                }
            }
            if (m < A.nt) {
                int nvam = plasma_tile_nview(A, m);
                core_omp_zlantr_aux(PlasmaInfNorm, uplo, diag,
                                    mvam, nvam,
                                    A(m, m), ldam, 
                                    &work[A.m*m+m*A.nb],
                                    sequence, request);
            }
        }
        #pragma omp taskwait
        workspace = work + A.nt*A.m;
        core_omp_dlange(PlasmaInfNorm,
                        A.m, A.nt,
                        work, A.m,
                        workspace, value,
                        sequence, request);
        break;
    //======================
    // PlasmaFrobeniusNorm
    //======================
    case PlasmaFrobeniusNorm:
        scale = work;
        sumsq = work + A.mt*A.nt;
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            if (uplo == PlasmaLower) {
                for (int n = 0; n < imin(m, A.nt); n++) {
                    int nvan = plasma_tile_nview(A, n);                    
                    core_omp_zgessq(mvam, nvan,
                                    A(m, n), ldam,
                                    &scale[A.mt*n+m], &sumsq[A.mt*n+m],
                                    sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_zgessq(mvam, nvan,
                                    A(m, n), ldam,
                                    &scale[A.mt*n+m], &sumsq[A.mt*n+m],
                                    sequence, request);
                }
            }
            if (m < A.nt) {
                int nvam = plasma_tile_nview(A, m);
                core_omp_ztrssq(uplo, diag,
                                mvam, nvam,
                                A(m, m), ldam, 
                                &scale[A.mt*m+m], &sumsq[A.mt*m+m],
                                sequence, request);
            }
        }
        #pragma omp taskwait
        core_omp_dgessq_aux(A.mt*A.nt,
                            scale, sumsq,
                            value,
                            sequence, request);
        break;
    }
}
