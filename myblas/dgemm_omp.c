/**
 *
 * @file dgemm_omp.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Template to develop the OpenMP version
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-30
 *
 */
#include "myblas.h"
#include <omp.h>
#include "algonum.h"

static int dgemm_seq_block_size = -1;

int dgemm_seq_omp(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
              CBLAS_TRANSPOSE transB, const int M, const int N,
              const int K, const double alpha, const double *A,
              const int lda, const double *B, const int ldb,
              const double beta, double *C, const int ldc)
{
    int m, n, k;
#pragma omp parallel for private(m, n, k) collapse(2)
for (m = 0; m < M; m++)
{   
    for (n = 0; n < N; n++)
    {
        double acc = 0.0;

        for (k = 0; k < K; k++)
        {
            double a_val, b_val;

            if (transA == CblasNoTrans)
                a_val = A[lda * k + m];
            else
                a_val = A[lda * m + k];

            if (transB == CblasNoTrans)
                b_val = B[ldb * n + k];
            else
                b_val = B[ldb * k + n];

            acc += a_val * b_val;
        }

        C[ldc * n + m] = beta * C[ldc * n + m] + alpha * acc;
    }
}

    return ALGONUM_SUCCESS;
}










    

int dgemm_block_omp(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
              CBLAS_TRANSPOSE transB, const int M, const int N,
              const int K, const double alpha, const double *A,
              const int lda, const double *B, const int ldb,
              const double beta, double *C, const int ldc)
{
    int block_size = dgemm_seq_block_size;

    if (block_size >= N)
    {
        return dgemm_base(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    }

#pragma omp parallel for collapse(2)
    for (int i = 0; i < M; i += block_size)
    {
        for (int j = 0; j < N; j += block_size)
        {
            for (int k = 0; k < K; k += block_size)
            {
                // Each thread works on a   different block
                dgemm_base(layout, transA, transB, block_size, block_size, block_size,
                           alpha, &A[i + k * lda], lda, &B[k + j * ldb], ldb,
                           (k == 0) ? beta : 1, &C[i + j * ldc], ldc);
            }
        }
    }
}

int dgemm_omp(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
              CBLAS_TRANSPOSE transB, const int M, const int N,
              const int K, const double alpha, const double *A,
              const int lda, const double *B, const int ldb,
              const double beta, double *C, const int ldc)
{
    return dgemm_seq_omp(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

/* To make sure we use the right prototype */
static dgemm_fct_t valid_dgemm_omp __attribute__ ((unused)) = dgemm_omp;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgemm_omp;

/**
 * @brief Registration function
 */
void dgemm_omp_init( void ) __attribute__( ( constructor ) );
void
dgemm_omp_init( void )
{
    fct_dgemm_omp.mpi    = 0;
    fct_dgemm_omp.tiled  = 0;
    fct_dgemm_omp.starpu = 0;
    fct_dgemm_omp.name   = "omp";
    fct_dgemm_omp.helper = "Scalar OpenMP implementation of the dgemm";
    fct_dgemm_omp.fctptr = dgemm_omp;
    fct_dgemm_omp.next   = NULL;

    register_fct( &fct_dgemm_omp, ALGO_GEMM );
    dgemm_seq_block_size = myblas_getenv_value_int("BLOCKSIZE", 32);

}
