/**
 *
 * @file dgemm_gpu.cu
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Template for the initial sequential GEMM function
 *
 * This file is a template that you can copy/paste as many times as
 * you like to create new versions of the GEMM implementation.
 *
 * To do that, replace all occurence of the TEMPLATE keyword by the shortname
 * you like, and save the file under the name dgemm_TEMPLATE.c
 * Then, add the file to the CMakeLists.txt gemm list, compile and enjoy.
 *
 * @version 0.2.0
 * @author YOURSELF
 * @date 2021-09-30
 *
 */
#include "myblas.h"

int dgemm_gpu( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                    CBLAS_TRANSPOSE transB, const int M, const int N,
                    const int K, const double alpha, const double *A,
                    const int lda, const double *B, const int ldb,
                    const double beta, double *C, const int ldc )
{
    /* Here is where you put your own code */

    return ALGONUM_NOT_IMPLEMENTED; /* To return in all cases that are not implemented */
    return ALGONUM_SUCCESS; /* To return in all implemented cases */
}

/* To make sure we use the right prototype */
static dgemm_fct_t valid_dgemm_gpu __attribute__ ((unused)) = dgemm_gpu;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgemm_gpu;

/**
 * @brief Registration function
 */
void dgemm_gpu_init( void ) __attribute__( ( constructor ) );
void
dgemm_gpu_init( void )
{
    fct_dgemm_TEMPLATE.tiled  = 0;
    fct_dgemm_TEMPLATE.starpu = 0;
    fct_dgemm_TEMPLATE.name   = "gpu";
    fct_dgemm_TEMPLATE.helper = "CUDA version of DGEMM";
    fct_dgemm_TEMPLATE.fctptr = dgemm_gpu;
    fct_dgemm_TEMPLATE.next   = NULL;

    register_fct( &fct_dgemm_TEMPLATE, ALGO_GEMM );
}
