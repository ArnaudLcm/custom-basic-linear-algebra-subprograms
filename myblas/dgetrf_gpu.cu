/**
 *
 * @file dgetrf_gpu.cu
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Template for the initial sequential GETRF function
 *
 * This file is a template that you can copy/paste as many times as
 * you like to create new versions of the GETRF implementation.
 *
 * To do that, replace all occurence of the TEMPLATE keyword by the shortname
 * you like, and save the file under the name dgetrf_TEMPLATE.c
 * Then, add the file to the CMakeLists.txt GETRF list, compile and enjoy.
 *
 * @version 0.2.0
 * @author YOURSELF
 * @date 2021-09-30
 *
 */
#include "myblas.h"
#include <cuda_runtime.h>

#define REG_NB_ELEMENTS 4

__global__ void dgetrf_gpu_kernel(int M, int N, double *A, int lda) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (tid < M) {
        int k;
        __m256d reg_minus1 = _mm256_set1_pd(-1.0);

        for (k = 0; k < min(M, N); k++) {
            if (tid == k) {
                // Division
                A[lda * k + tid] /= A[lda * k + k];
            }
            __syncthreads();

            if (tid > k && tid < M) {
                // Calculs vectorisés
                __m256d reg_A = _mm256_loadu_pd(&A[lda * k + tid]);
                __m256d reg_KK = _mm256_set1_pd(A[lda * k + k]);
                reg_A = _mm256_div_pd(reg_A, reg_KK);
                _mm256_storeu_pd(&A[lda * k + tid], reg_A);
            }
            __syncthreads();

            if (tid > k && tid < M) {
                // Calculs vectorisés
                __m256d reg_AK = _mm256_loadu_pd(&A[lda * k + tid]);
                __m256d reg_AKN = _mm256_set1_pd(A[lda * tid + k]);
                __m256d reg_AN = _mm256_loadu_pd(&A[lda * tid + k]);
                reg_AK = _mm256_fmadd_pd(reg_minus1, reg_AK, _mm256_mul_pd(reg_AK, reg_AKN));
                reg_AN = _mm256_sub_pd(reg_AN, _mm256_mul_pd(reg_AK, reg_AN));
                _mm256_storeu_pd(&A[lda * tid + k], reg_AN);
            }
            __syncthreads();
        }
    }
}

int dgetrf_gpu( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                    CBLAS_TRANSPOSE transB, const int M, const int N,
                    const int K, const double alpha, const double *A,
                    const int lda, const double *B, const int ldb,
                    const double beta, double *C, const int ldc )
{
    /* Here is where you put your own code */

    double *d_A;

    // Allouer de la mémoire sur le GPU
    cudaMalloc((void**)&d_A, M * lda * sizeof(double));

    // Copier les données depuis le CPU vers le GPU
    cudaMemcpy(d_A, A, M * lda * sizeof(double), cudaMemcpyHostToDevice);

    // Définir la grille et la taille des blocs pour les threads du GPU
    dim3 gridSize((M + REG_NB_ELEMENTS - 1) / REG_NB_ELEMENTS, 1, 1);
    dim3 blockSize(REG_NB_ELEMENTS, 1, 1);

    // Appeler le noyau GPU
    dgetrf_gpu_kernel<<<gridSize, blockSize>>>(M, N, d_A, lda);

    // Attendre la fin des calculs du GPU
    cudaDeviceSynchronize();

    // Copier les résultats depuis le GPU vers le CPU
    cudaMemcpy(A, d_A, M * lda * sizeof(double), cudaMemcpyDeviceToHost);

    // Libérer la mémoire allouée sur le GPU
    cudaFree(d_A);

    return ALGONUM_SUCCESS;
}

/* To make sure we use the right prototype */
static dgetrf_fct_t valid_dgetrf_gpu __attribute__ ((unused)) = dgetrf_gpu;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgetrf_gpu;

/**
 * @brief Registration function
 */
void dgetrf_gpu_init( void ) __attribute__( ( constructor ) );
void
dgetrf_gpu_init( void )
{
    fct_dgetrf_TEMPLATE.tiled  = 0;
    fct_dgetrf_TEMPLATE.starpu = 0;
    fct_dgetrf_TEMPLATE.name   = "gpu";
    fct_dgetrf_TEMPLATE.helper = "CUDA version of dgetrf";
    fct_dgetrf_TEMPLATE.fctptr = dgetrf_gpu;
    fct_dgetrf_TEMPLATE.next   = NULL;

    register_fct( &fct_dgetrf_TEMPLATE, ALGO_GETRF );
}
