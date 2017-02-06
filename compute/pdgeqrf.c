/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzgeqrf.c, normal z -> d, Mon Feb  6 14:06:15 2017
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m, n) (double*)plasma_tile_addr(A, m, n)
#define T(m, n) (double*)plasma_tile_addr(T, m, n)

/***************************************************************************//**
 *  Parallel tile QR factorization - dynamic scheduling
 * @see plasma_omp_dgeqrf
 **/
void plasma_pdgeqrf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    // Set inner blocking from the T tile row-dimension.
    int ib = T.mb;

    for (int k = 0; k < imin(A.mt, A.nt); k++) {
        int mvak = plasma_tile_mview(A, k);
        int nvak = plasma_tile_nview(A, k);
        int ldak = plasma_tile_mmain(A, k);
        core_omp_dgeqrt(
            mvak, nvak, ib,
            A(k, k), ldak,
            T(k, k), T.mb,
            work,
            sequence, request);

        for (int n = k+1; n < A.nt; n++) {
            int nvan = plasma_tile_nview(A, n);
            core_omp_dormqr(
                PlasmaLeft, PlasmaTrans,
                mvak, nvan, imin(mvak, nvak), ib,
                A(k, k), ldak,
                T(k, k), T.mb,
                A(k, n), ldak,
                work,
                sequence, request);
        }
        for (int m = k+1; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            core_omp_dtsqrt(
                mvam, nvak, ib,
                A(k, k), ldak,
                A(m, k), ldam,
                T(m, k), T.mb,
                work,
                sequence, request);

            for (int n = k+1; n < A.nt; n++) {
                int nvan = plasma_tile_nview(A, n);
                core_omp_dtsmqr(
                    PlasmaLeft, PlasmaTrans,
                    A.mb, nvan, mvam, nvan, nvak, ib,
                    A(k, n), ldak,
                    A(m, n), ldam,
                    A(m, k), ldam,
                    T(m, k), T.mb,
                    work,
                    sequence, request);
            }
        }
    }
}
