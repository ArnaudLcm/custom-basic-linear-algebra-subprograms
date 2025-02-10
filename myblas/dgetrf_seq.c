/**
 *
 * @file dgetrf_seq.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Basic sequential implementation of the dgetrf.
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-21
 *
 */
#include "myblas.h"
#include <assert.h>
#define Up 0    
#define Low 1
#include <immintrin.h>

static int dgemm_seq_block_size = -1;

int dgetrf_base(CBLAS_LAYOUT layout, int M, int N, double *A, int lda)
{

    int m, n, k;
    int K = (M > N) ? N : M;
    for (k = 0; k < K; k++)
    {
        for (m = k + 1; m < M; m++)
        {
            A[lda * k + m] = A[lda * k + m] / A[lda * k + k];
            for (n = k + 1; n < N; n++)
            {
                A[lda * n + m] = A[lda * n + m] - A[lda * k + m] * A[lda * n + k];
            }
        }
    }

        // display_matrix(A, M, N, lda);

    return ALGONUM_SUCCESS; /* Success */
}




int dgetrf_bloc(CBLAS_LAYOUT layout, int M, int N, double *A, int lda)
{
    int mm, nn, k, kk;
    int K = ( M > N ) ? N : M;
    int block_size = dgemm_seq_block_size;
    printf("Blocksize : %d\n", block_size);
    int nbr_rows = (M + block_size -1) / block_size;
    int nbr_columns = (N + block_size -1) / block_size;
    if(block_size > M)
        return dgetrf_base(layout, M, N, A, lda);
    for (k = 0; k < nbr_rows; k++)
    {
        kk = k == (nbr_rows-1) ? N - k * block_size : block_size;
        dgetrf_base(layout, kk, kk, A + (lda * (k * kk) + k * kk), lda);
        for (int i = k + 1; i < nbr_rows; i++)
        {
            mm = i == (nbr_rows-1) ? M - i * block_size : block_size;
            cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, mm, block_size, 1,  &A[k * block_size + (k * block_size) * lda], lda,  &A[i * block_size + (k * block_size) * lda], lda);
        }

        for(int j=k+1; j<nbr_columns; j++){
            nn = j == (nbr_rows-1) ? N - j * block_size : block_size;
            cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, block_size, nn, 1,  &A[k * block_size + (k * block_size) * lda],lda ,  &A[k * block_size + (j * block_size) * lda], lda);
        }
    
        for (int j = k + 1; j < nbr_columns; j++)
        { 
            for (int i = k + 1; i < nbr_rows; i++)
            {
                mm = i == (nbr_rows-1) ? M - i * block_size : block_size;
                nn = j == (nbr_rows-1) ? N - j * block_size : block_size;
                dgemm_base(layout, CblasNoTrans, CblasNoTrans, mm, nn, block_size, -1, &A[i * block_size + (k * block_size) * lda], lda, &A[k * block_size + (j * block_size) * lda], lda, 1, &A[i * block_size + (j * block_size) * lda], lda);
            }
        }
    }
    // display_matrix(A, M, N, lda);

    return ALGONUM_SUCCESS;
}

/* Size in bytes of an SIMD register */
#define REG_BYTES (sizeof(__m256d))

/* Number of elements in an SIMD register */
#define REG_NB_ELEMENTS (REG_BYTES / sizeof(double))

#ifdef __AVX2__

int dgetrf_vectorized(CBLAS_LAYOUT layout, int M, int N, double *A, int lda)
{

    int m, n, k;
    int K = (M > N) ? N : M;
    __m256d diag_k, k_m, n_m, n_k;

    for (k = 0; k < K; k++)
    {

        diag_k = _mm256_set1_pd(A[lda * k + k]);
        for (m = k + 1; m + REG_NB_ELEMENTS <= M; m += REG_NB_ELEMENTS)
        {
            // A[ lda * k + m ] = A[ lda * k + m ] / A[ lda * k + k ];

            k_m = _mm256_loadu_pd(&A[lda * k + m]);
            k_m = _mm256_div_pd(k_m, diag_k);
            _mm256_storeu_pd(&A[lda * k + m], k_m);

            for (n = k + 1; n < N; n++)
            {
                // A[ lda * n + m ] = A[ lda * n + m ] - A[ lda * k + m ] * A[ lda * n + k ];
                n_k = _mm256_set1_pd(-A[lda * n + k]);

                n_m = _mm256_loadu_pd(&A[lda * n + m]);
                n_m = _mm256_fmadd_pd(n_k, k_m, n_m);
                _mm256_storeu_pd(&A[lda * n + m], n_m);
            }
        }
        for (; m < M; m++)
        {
            A[lda * k + m] = A[lda * k + m] / A[lda * k + k];

            for (n = k + 1; n < N; n++)
            {

                A[lda * n + m] = A[lda * n + m] - A[lda * k + m] * A[lda * n + k];
            }
        }
    }

    return ALGONUM_SUCCESS; /* Success */
}

#endif

int dgetrf_seq(CBLAS_LAYOUT layout, int M, int N, double *A, int lda)
{
    return dgetrf_bloc(layout, M, N, A, lda);
}

/* To make sure we use the right prototype */
static dgetrf_fct_t valid_dgetrf_seq __attribute__((unused)) = dgetrf_seq;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgetrf_seq;

/**
 * @brief Registration function
 */
void dgetrf_seq_init(void) __attribute__((constructor));
void dgetrf_seq_init(void)
{
    fct_dgetrf_seq.mpi = 0;
    fct_dgetrf_seq.tiled = 0;
    fct_dgetrf_seq.starpu = 0;
    fct_dgetrf_seq.name = "seq";
    fct_dgetrf_seq.helper = "Basic sequential implementation of the dgetrf";
    fct_dgetrf_seq.fctptr = dgetrf_seq;
    fct_dgetrf_seq.next = NULL;

    register_fct(&fct_dgetrf_seq, ALGO_GETRF);

    dgemm_seq_block_size = myblas_getenv_value_int("BLOCKSIZE", 10);

}
