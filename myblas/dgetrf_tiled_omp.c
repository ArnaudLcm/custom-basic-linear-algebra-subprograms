/**
 *
 * @file dgetrf_tiled_omp.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief OpenMP tile-based implementation of the dgetrf.
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-21
 *
 */
#include "myblas.h"

int
dgetrf_tiled_omp( CBLAS_LAYOUT layout,
                  int M, int N, int b, double **Atiled)
{

    
    int MT = my_iceil( M, b );
    int NT = my_iceil( N, b );
    int KT = my_imin( MT, NT );
    int m, n, k;
    

    for (k = 0; k < KT; k++)
    {
        dgetrf_base(layout, b, b, Atiled[k*MT+k], b);
        #pragma omp parallel for
        for (int m = k + 1; m < MT; m++)
        {
            cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, b, b, 1,Atiled[k*MT+k], b, Atiled[k*MT+m], b);
        }

        #pragma omp parallel for
        for(int n=k+1; n<NT; n++){
            cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, b, b, 1, Atiled[k*MT+k],b , Atiled[n*MT+k], b);
        }

        #pragma omp parallel for collapse(2)
        for (int n = k + 1; n < NT; n++)
        {
            for (int m = k + 1; m < MT; m++)
            {
                dgemm_base(layout, CblasNoTrans, CblasNoTrans, b, b, b, -1, Atiled[k*MT+m], b, Atiled[m*MT+k], b, 1, Atiled[m*MT+n], b);
            }
        }
    }
    

    return ALGONUM_SUCCESS;
}


/* To make sure we use the right prototype */
static dgetrf_tiled_fct_t valid_dgetrf_tiled_omp __attribute__ ((unused)) = dgetrf_tiled_omp;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgetrf_tiled_omp;

/**
 * @brief Registration function
 */
void dgetrf_tiled_omp_init( void ) __attribute__( ( constructor ) );
void
dgetrf_tiled_omp_init( void )
{
    fct_dgetrf_tiled_omp.mpi    = 0;
    fct_dgetrf_tiled_omp.tiled  = 1;
    fct_dgetrf_tiled_omp.starpu = 0;
    fct_dgetrf_tiled_omp.name   = "omp-t";
    fct_dgetrf_tiled_omp.helper = "OpenMP tile-based implementation of the dgetrf";
    fct_dgetrf_tiled_omp.fctptr = dgetrf_tiled_omp;
    fct_dgetrf_tiled_omp.next   = NULL;

    register_fct( &fct_dgetrf_tiled_omp, ALGO_GETRF );
}
