/**
 *
 * @file dgemm_seq.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Sequential version of the Matrix-Matrix multiply operation
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-30
 *
 */
#include "myblas.h"
#include "assert.h"
#include <assert.h>
#include <immintrin.h>

/* Size in bytes of an SIMD register */
#define REG_BYTES (sizeof(__m256d))

/* Number of elements in an SIMD register */
#define REG_NB_ELEMENTS (REG_BYTES / sizeof(double))

// Exemple of ways to add additionnal parameters to your kernel
// See the registration function to change its value
static int dgemm_seq_block_size = -1;

int dgemm_base(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
               CBLAS_TRANSPOSE transB, const int M, const int N,
               const int K, const double alpha, const double *A,
               const int lda, const double *B, const int ldb,
               const double beta, double *C, const int ldc)
{

    int m, n, k;
    if (transA == CblasNoTrans)
    {
        if (transB == CblasNoTrans)
        {
            for (n = 0; n < N; n++)
            {
                if (beta != 1.)
                {
                    for (m = 0; m < M; m++)
                    {
                        C[ldc * n + m] = beta * C[ldc * n + m];
                    }
                }
                for (k = 0; k < K; k++)
                {
                    for (m = 0; m < M; m++)
                    {
                        C[ldc * n + m] += alpha * A[lda * k + m] * B[ldb * n + k];
                    }
                }
            }
        }
        else
        {
            for (m = 0; m < M; m++)
            {
                for (n = 0; n < N; n++)
                {
                    if (beta != 1.)
                    {
                        C[ldc * n + m] = beta * C[ldc * n + m];
                    }
                    for (k = 0; k < K; k++)
                    {
                        C[ldc * n + m] += alpha * A[lda * k + m] * B[ldb * k + n];
                    }
                }
            }
        }
    }
    else
    {
        if (transB == CblasNoTrans)
        {
            for (m = 0; m < M; m++)
            {
                for (n = 0; n < N; n++)
                {
                    if (beta != 1.)
                    {
                        C[ldc * n + m] = beta * C[ldc * n + m];
                    }
                    for (k = 0; k < K; k++)
                    {
                        C[ldc * n + m] += alpha * A[lda * m + k] * B[ldb * n + k];
                    }
                }
            }
        }
        else
        {
            for (m = 0; m < M; m++)
            {
                for (n = 0; n < N; n++)
                {
                    if (beta != 1.)
                    {
                        C[ldc * n + m] = beta * C[ldc * n + m];
                    }
                    for (k = 0; k < K; k++)
                    {
                        C[ldc * n + m] += alpha * A[lda * m + k] * B[ldb * k + n];
                    }
                }
            }
        }
    }

    return ALGONUM_SUCCESS;
}



#ifdef __AVX2__

int dgemm_vectorized(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                     CBLAS_TRANSPOSE transB, const int M, const int N,
                     const int K, const double alpha, const double *A,
                     const int lda, const double *B, const int ldb,
                     const double beta, double *C, const int ldc)
{

    // assert((M % REG_NB_ELEMENTS) == 0);
    // assert((N % REG_NB_ELEMENTS) == 0);
    // assert((K % REG_NB_ELEMENTS) == 0);

    if (((M % REG_NB_ELEMENTS) != 0) || ((N % REG_NB_ELEMENTS) != 0) || ((K % REG_NB_ELEMENTS) != 0))
    {
        printf("Wrong matrix size\n");
        return ALGONUM_FAIL;
    }

    int m, n, k;

    __m256d reg_alpha;
    reg_alpha = _mm256_set1_pd(alpha);

    __m256d reg_beta;
    reg_beta = _mm256_set1_pd(beta);

    if (transA == CblasNoTrans)
    {
        if (transB == CblasNoTrans)
        {

            for (n = 0; n < N; n++)
            {

                if (beta != 1.)
                {
                    for (m = 0; m < M; m += REG_NB_ELEMENTS)
                    {
                        __m256d reg_C;
                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_C = _mm256_mul_pd(reg_beta, reg_C);

                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);

                        // C[ ldc * n + m ] = beta * C[ ldc * n + m ];
                    }
                }

                for (k = 0; k < K; k++)
                {

                    __m256d reg_B;
                    reg_B = _mm256_set1_pd(B[ldb * n + k]);
                    for (m = 0; m < M; m += REG_NB_ELEMENTS)
                    {
                        __m256d reg_C;
                        __m256d reg_A;

                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_A = _mm256_loadu_pd(&A[lda * k + m]);

                        reg_A = _mm256_mul_pd(reg_alpha, reg_A);
                        reg_C = _mm256_fmadd_pd(reg_A, reg_B, reg_C);
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);

                        // C[ ldc * n + m ] += alpha * A[ lda * k + m ] * B[ ldb * n + k ];
                    }
                }
            }
        }
        else
        {
            for (m = 0; m < M; m += REG_NB_ELEMENTS)
            {
                for (n = 0; n < N; n++)
                {

                    if (beta != 1.)
                    {
                        __m256d reg_C;
                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_C = _mm256_mul_pd(reg_beta, reg_C);
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);
                        // C[ ldc * n + m ] = beta * C[ ldc * n + m ];
                    }
                    for (k = 0; k < K; k++)
                    {

                        __m256d reg_C;
                        __m256d reg_A;
                        __m256d reg_B;

                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_A = _mm256_loadu_pd(&A[lda * k + m]);
                        reg_B = _mm256_set1_pd(B[ldb * k + n]);

                        reg_A = _mm256_mul_pd(reg_alpha, reg_A);
                        reg_C = _mm256_fmadd_pd(reg_A, reg_B, reg_C);
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);
                        // C[ ldc * n + m ] += alpha * A[ lda * k + m ] * B[ ldb * k + n ];
                    }
                }
            }
        }
    }
    else
    {
        if (transB == CblasNoTrans)
        {
            for (m = 0; m < M; m += REG_NB_ELEMENTS)
            {

                for (n = 0; n < N; n++)
                {

                    if (beta != 1.)
                    {
                        __m256d reg_C;
                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_C = _mm256_mul_pd(reg_beta, reg_C);
                        // if ((ldc * n_col + m )>=(N*M))
                        // {
                        //     printf("n : %d, m : %d, ldc : %d, reg_nb_elt : %ld ",n_col,m,ldc,REG_NB_ELEMENTS);
                        //     printf("Fail 1\n");
                        //     return ALGONUM_FAIL;
                        // }
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);
                        // C[ ldc * n + m ] = beta * C[ ldc * n + m ];
                    }
                    for (k = 0; k < K; k++)
                    {
                        int k_col = k / REG_NB_ELEMENTS;
                        __m256d reg_C;
                        __m256d reg_A;
                        __m256d reg_B;

                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_A = _mm256_set1_pd(A[lda * m + k]);
                        reg_B = _mm256_loadu_pd(&B[ldb * n + k]);

                        reg_A = _mm256_mul_pd(reg_alpha, reg_A);
                        reg_C = _mm256_fmadd_pd(reg_A, reg_B, reg_C);
                        // if ((ldc * n_col + m)>=(N*M))
                        // {
                        //     printf("n : %d, m : %d, ldc : %d, reg_nb_elt : %ld ",n_col,m,ldc,REG_NB_ELEMENTS);
                        //     printf("Fail 2\n");
                        //     return ALGONUM_FAIL;
                        // }
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);
                        // C[ ldc * n + m ] += alpha * A[ lda * m + k ] * B[ ldb * n + k ];
                    }
                }
            }
        }
        else
        {
            for (m = 0; m < M; m += REG_NB_ELEMENTS)
            {

                for (n = 0; n < N; n++)
                {

                    if (beta != 1.)
                    {
                        // printf("n : %d, m : %d, ldc : %d, reg_nb_elt : %ld\n",n,m,ldc,REG_NB_ELEMENTS);
                        // printf("indice : %d\n", ldc * n + m);
                        __m256d reg_C;
                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_C = _mm256_mul_pd(reg_beta, reg_C);
                        // if ((ldc * n_col + m)>=(N*ldc))
                        // {
                        //     printf("n : %d, m : %d, ldc : %d, reg_nb_elt : %ld ",n_col,m,ldc,REG_NB_ELEMENTS);
                        //     printf("Fail 3\n");
                        //     return ALGONUM_FAIL;
                        // }
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);
                        // C[ ldc * n + m ] = beta * C[ ldc * n + m ];
                    }
                    for (k = 0; k < K; k++)
                    {

                        __m256d reg_C;
                        __m256d reg_A;
                        __m256d reg_B;

                        reg_C = _mm256_loadu_pd(&C[ldc * n + m]);
                        reg_A = _mm256_set1_pd(A[lda * m + k]);
                        reg_B = _mm256_loadu_pd(&B[ldb * k + n]);

                        reg_A = _mm256_mul_pd(reg_alpha, reg_A);
                        reg_C = _mm256_fmadd_pd(reg_A, reg_B, reg_C);
                        // if ((ldc * n_col + m)>=(N*M))
                        // {
                        //     printf("n : %d, m : %d, ldc : %d, reg_nb_elt : %ld ",n_col,m,ldc,REG_NB_ELEMENTS);
                        //     printf("Fail 4\n");
                        //     return ALGONUM_FAIL;
                        // }
                        _mm256_storeu_pd(&C[ldc * n + m], reg_C);
                        // C[ ldc * n + m ] += alpha * A[ lda * m + k ] * B[ ldb * k + n ];
                    }
                }
            }
        }
    }

    return ALGONUM_SUCCESS;
}

#endif



void display_matrix(const double *A, int M, int N, int lda)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f ", A[j * lda + i]);
        }
        printf("\n");
    }
}

int dgemm_bloc(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
               CBLAS_TRANSPOSE transB, const int M, const int N,
               const int K, const double alpha, const double *A,
               const int lda, const double *B, const int ldb,
               const double beta, double *C, const int ldc)
{
    int block_size = dgemm_seq_block_size;

    if(block_size>=N){
        return dgemm_base(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    }
   

    for (size_t i = 0; i < M; i += block_size)
    {
        for (int j=0; j < N; j += block_size)
            dgemm_base(layout, transA, transB, block_size, block_size, block_size, alpha, &A[i], lda, &B[j*ldb], ldb, beta, &C[i + j * ldc ], ldc);

            
        for (size_t k = block_size; k < N; k += block_size)
            
            {
                

            for (int j  = 0; j < N; j += block_size)

                {

                    dgemm_base(layout, transA, transB, block_size, block_size, block_size, alpha, &A[i + k * lda], lda, &B[k + j * ldb], ldb, 1, &C[i + j * ldc], ldc);
                }
            }
    }

}

int dgemm_bloc_vecto(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
               CBLAS_TRANSPOSE transB, const int M, const int N,
               const int K, const double alpha, const double *A,
               const int lda, const double *B, const int ldb,
               const double beta, double *C, const int ldc)
{
    int block_size = dgemm_seq_block_size;

    if(block_size>=N){
        return dgemm_vectorized(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    for (size_t i = 0; i < M; i += block_size)
    {
        for (int j=0; j < N; j += block_size)
            dgemm_vectorized(layout, transA, transB, block_size, block_size, block_size, alpha, &A[i], lda, &B[j*ldb], ldb, beta, &C[i + j * ldc ], ldc);

            
    for (size_t k = block_size; k < N; k += block_size)
        
        {
            

        for (int j  = 0; j < N; j += block_size)

            {

                dgemm_vectorized(layout, transA, transB, block_size, block_size, block_size, alpha, &A[i + k * lda], lda, &B[k + j * ldb], ldb, 1, &C[i + j * ldc], ldc);
            }
        }
    }

    // printf("\n \nC=");
    // display_matrix(C, M, N);
}


/* To make sure we use the right prototype */
static dgemm_fct_t valid_dgemm_seq __attribute__((unused)) = dgemm_seq;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgemm_seq;

/**
 * @brief Registration function
 */

int dgemm_seq(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
              CBLAS_TRANSPOSE transB, const int M, const int N,
              const int K, const double alpha, const double *A,
              const int lda, const double *B, const int ldb,
              const double beta, double *C, const int ldc)
{
    return dgemm_bloc(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void dgemm_seq_init(void) __attribute__((constructor));
void dgemm_seq_init(void)
{
    fct_dgemm_seq.mpi = 0;
    fct_dgemm_seq.tiled = 0;
    fct_dgemm_seq.starpu = 0;
    fct_dgemm_seq.name = "seq";
    fct_dgemm_seq.helper = "Sequential version of DGEMM";
    fct_dgemm_seq.fctptr = dgemm_seq;
    fct_dgemm_seq.next = NULL;

    register_fct(&fct_dgemm_seq, ALGO_GEMM);

    /* Read the value of dgemm_block_size */
    dgemm_seq_block_size = myblas_getenv_value_int("BLOCKSIZE", 32);
}
