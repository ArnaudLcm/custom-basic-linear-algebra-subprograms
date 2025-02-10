/**
 *
 * @file dgetrf_omp.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief OpenMP implementation of the dgetrf.
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-21
 *
 */
#include "myblas.h"
#include <assert.h>
#include <omp.h>
#include "algonum.h"

static int dgetrf_seq_block_size = 2;



int
dgetrf_omp( CBLAS_LAYOUT layout, int M, int N, double *A, int lda )
{
    int m, n, k;
    int K = ( M > N ) ? N : M;
    int block_size = dgetrf_seq_block_size;
    printf("Blocksize : %d\n", block_size);
    assert(lda/block_size == 0);
    int nbr_rows = (M + block_size -1) / block_size;
    int nbr_columns = (N + block_size -1) / block_size;
    for (k = 0; k < nbr_rows; k++)
    {
        dgetrf_base(layout, block_size, block_size, &A[k * block_size + (k * block_size) * lda] , lda);
        #pragma omp parallel for
        for (int i = k + 1; i < nbr_rows; i++)
        {
            cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, block_size, block_size, 1, &A[k * block_size + (k * block_size) * lda] , lda, &A[i * block_size + (k * block_size) * lda], lda);
        }

        #pragma omp parallel for
        for(int j=k+1; j<nbr_columns; j++){
            cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, block_size, block_size, 1,&A[k * block_size + (k * block_size) * lda] ,block_size , &A[k * block_size + (j * block_size) * lda], lda);
        }

        #pragma omp parallel for collapse(2)
        for (int j = k + 1; j < nbr_columns; j++)
        {
            for (int i = k + 1; i < nbr_rows; i++)
            {
                dgemm_base(layout, CblasNoTrans, CblasNoTrans, block_size, block_size, block_size, -1, &A[i * block_size + (k * block_size) * lda], lda, &A[k * block_size + (j * block_size) * lda], lda, 1, &A[i * block_size + (j * block_size) * lda], lda);
            }
        }
    }
    return ALGONUM_SUCCESS;
}








/* To make sure we use the right prototype */
static dgetrf_fct_t valid_dgetrf_omp __attribute__ ((unused)) = dgetrf_omp;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgetrf_omp;

/**
 * @brief Registration function
 */
void dgetrf_omp_init( void ) __attribute__( ( constructor ) );
void
dgetrf_omp_init( void )
{
    fct_dgetrf_omp.mpi    = 0;
    fct_dgetrf_omp.tiled  = 0;
    fct_dgetrf_omp.starpu = 0;
    fct_dgetrf_omp.name   = "omp";
    fct_dgetrf_omp.helper = "OpenMP implementation of the dgetrf";
    fct_dgetrf_omp.fctptr = dgetrf_omp;
    fct_dgetrf_omp.next   = NULL;

    register_fct( &fct_dgetrf_omp, ALGO_GETRF );
    dgetrf_seq_block_size = myblas_getenv_value_int("BLOCKSIZE", 32);

}
