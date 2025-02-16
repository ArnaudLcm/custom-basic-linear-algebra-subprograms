/**
 *
 * @file perf_dgemm.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Binary to assess the performance of a GEMM implementation
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-30
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include "algonum.h"

int main( int argc, char **argv )
{
    dplrnt_tiled_fct_t tested_tiled_dplrnt = dplrnt_tiled;
    option_t *options = &global_options;
    int i;

    algonum_init( argc, argv, options, ALGO_GEMM );

#if defined(ENABLE_STARPU)
    if ( options->fct->starpu ) {
        tested_tiled_dplrnt = dplrnt_tiled_starpu;
    }
#endif

    for( i=0; i<options->iter; i++ ) {
        if ( options->fct->tiled ) {
            testone_dgemm_tiled( tested_tiled_dplrnt,
                                 options->fct->fctptr,
                                 options->transA, options->transB,
                                 options->M, options->N, options->K,
                                 options->b, options->check );
        }
        else {
            testone_dgemm( options->fct->fctptr,
                           options->transA, options->transB,
                           options->M, options->N, options->K, options->check );
        }
    }

    algonum_exit( options, ALGO_GEMM );

    return EXIT_SUCCESS;
}
