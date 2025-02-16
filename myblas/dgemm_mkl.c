/**
 *
 * @file dgemm_mkl.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Register the MKL function for references
 *
 * This version of the DGEMM should not be modified. It is provided to
 * you in order to have a reference point for the sequential
 * algorithm.
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-30
 *
 */
#include "myblas.h"

int dgemm_mkl( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
               CBLAS_TRANSPOSE transB, const int M, const int N,
               const int K, const double alpha, const double *A,
               const int lda, const double *B, const int ldb,
               const double beta, double *C, const int ldc )
{
    cblas_dgemm( layout, transA, transB, M, N, K,
                 alpha, A, lda, B, ldb, beta, C, ldc );

    return ALGONUM_SUCCESS;
}

/* To make sure we use the right prototype */
static dgemm_fct_t valid_dgemm_mkl __attribute__ ((unused)) = dgemm_mkl;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgemm_mkl;

/**
 * @brief Registration function
 */
void dgemm_mkl_init( void ) __attribute__( ( constructor ) );
void
dgemm_mkl_init( void )
{
    fct_dgemm_mkl.mpi    = 0;
    fct_dgemm_mkl.tiled  = 0;
    fct_dgemm_mkl.starpu = 0;
    fct_dgemm_mkl.name   = "mkl";
    fct_dgemm_mkl.helper = "MKL implementation of the dgemm";
    fct_dgemm_mkl.fctptr = dgemm_mkl;
    fct_dgemm_mkl.next   = NULL;

    register_fct( &fct_dgemm_mkl, ALGO_GEMM );
}
