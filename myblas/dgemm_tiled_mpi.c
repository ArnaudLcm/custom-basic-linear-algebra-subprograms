/**
 *
 * @file dgemm_tiled_mpi.c
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

int dgemm_tiled_mpi( CBLAS_LAYOUT layout,
                     CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                     int M, int N, int K, int b,
                     double alpha, const double **A,
                                   const double **B,
                     double beta,        double **C )
{
#if !defined(ENABLE_MPI)
    return ALGONUM_NOT_IMPLEMENTED;
#else
    MPI_Status    status;
    const double *Aptr, *Bptr;
    double        wsA[b*b], wsB[b*b];
    double        lbeta;

    /* Let's compute the total number of tiles with a *ceil* */
    int MT = my_iceil( M, b );
    int NT = my_iceil( N, b );
    int KT = my_iceil( K, b );

    int m, n, k, mm, nn, kk;
    int ownerA, ownerB, ownerC;

    if ( transA == CblasNoTrans ) {
        if ( transB == CblasNoTrans ) {
            for( m=0; m<MT; m++ ) {
                mm = m == (MT-1) ? M - m * b : b;

                for( n=0; n<NT; n++ ) {
                    nn = n == (NT-1) ? N - n * b : b;
                    ownerC = get_rank_of( m, n );

                    if ( global_options.mpirank == ownerC ) {
                        for( k=0; k<KT; k++ ) {
                            kk = k == (KT-1) ? K - k * b : b;
                            lbeta = (k == 0) ? beta : 1.;

                            ownerA = get_rank_of( m, k );
                            ownerB = get_rank_of( k, n );   

                            if ( global_options.mpirank == ownerA ) {
                                Aptr = A[MT * k + m];
                            }
                            else {
                                MPI_Recv( wsA, b*kk, MPI_DOUBLE,
                                          ownerA, MT * k + m, MPI_COMM_WORLD, &status );
                                Aptr = wsA;
                            }

                            if ( global_options.mpirank == ownerB ) {
                                Bptr = B[KT * n + k];
                            }
                            else {
                                MPI_Recv( wsB, b*nn, MPI_DOUBLE,
                                          ownerB, KT * n + k, MPI_COMM_WORLD, &status );
                                Bptr = wsB;
                            }

                            dgemm_seq( CblasColMajor, transA, transB,
                                       mm, nn, kk,
                                       alpha, Aptr, b, Bptr, b,
                                       lbeta, C[ MT * n + m ], b );
                        }
                    }
                    else {
                        for( k=0; k<KT; k++ ) {
                            kk = k == (KT-1) ? K - k * b : b;
                            lbeta = (k == 0) ? beta : 1.;

                            ownerA = get_rank_of( m, k );
                            ownerB = get_rank_of( k, n );

                            if ( global_options.mpirank == ownerA ) {
                                MPI_Send( A[MT * k + m], b * kk, MPI_DOUBLE,
                                          ownerC, MT * k + m, MPI_COMM_WORLD );
                            }

                            if ( global_options.mpirank == ownerB ) {
                                MPI_Send( B[KT * n + k], b * nn, MPI_DOUBLE,
                                          ownerC, KT * n + k, MPI_COMM_WORLD );
                            }
                        }
                    }
                }
            }
        }
        else {
            return ALGONUM_NOT_IMPLEMENTED;
        }
    }
    else {
        return ALGONUM_NOT_IMPLEMENTED;
    }

    return ALGONUM_SUCCESS;
#endif
}


void matrix_sum(void *in, void *inout, int *len, MPI_Datatype *dptr) {
    int i, j;
    int ld = sqrt(len[0]);
    double *matrix_in = (double *)in;
    double *matrix_inout = (double *)inout;
    for (i = 0; i < ld; i++) {
        for (j = 0; j < ld; j++) {
            matrix_inout[j * ld+ i] += matrix_in[j * ld + i];
        }
    }
}

void matrix_somme(double *in, double *inout, int len) {
    int i, j;
    for (i = 0; i < len; i++) {
        for (j = 0; j < len; j++) {
            inout[j * len+ i] += in[j * len + i];
        }
    }
}

void matrix_multiply(double *inout, int len, double beta) {
    int i, j;
    for (i = 0; i < len; i++) {
        for (j = 0; j < len; j++) {
            inout[j * len+ i] *= beta;
        }
    }
}





int dgemm_dist_compute_centralize(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
               CBLAS_TRANSPOSE transB, const int M, const int N,
               const int K, const double alpha, const double *A,
               const int lda, const double *B, const int ldb,
               const double beta, double *C, const int ldc)
{
    
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
   
    int b = M/10;
    int mm, nn, kk;
    double* C_local ;
    double* C_global;
    MPI_Status    status;

    
    double ** Atiled = lapack2tile(M, N, b, A, lda);
    double ** Btiled = lapack2tile(M, N, b, B, ldb);
    double ** Ctiled = lapack2tile(M, N, b, C, ldc);


    int MT = my_iceil( M, b );
    int NT = my_iceil( N, b );
    int KT = my_iceil( K, b );

    MPI_Op mat_add_op;
    MPI_Op_create(&matrix_sum, 1, &mat_add_op);
    if(b>=N){
            return dgemm_base(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
        }


    for (int m = 0; m < MT; m += 1)
    {       
        mm = m == (MT-1) ? M - m * b : b;
        for (int n = 0; n < NT; n += 1)
            {
            nn = n == (NT-1) ? N - n * b : b;
            if(rank == n%size){
                matrix_multiply(Ctiled[n*NT+m], b, beta);
            }
            C_global  = calloc(b*b, sizeof(double));
            C_local  = calloc(b*b, sizeof(double));
            double* tile_recv = calloc(b*b, sizeof(double)) ;
            for (int k = 0; k < KT; k += 1)    
                {
                    if(k%size == rank){
                        if(k%size == n%size){
                            tile_recv = Btiled[n*NT+k];
                        }
                        else{
                            MPI_Recv( tile_recv, b*b, MPI_DOUBLE,
                                          n%size, KT * n + k, comm, &status );
                        }
                        kk = k == (KT-1) ? K - k * b : b;
                        dgemm_base(layout, transA, transB, mm, nn, kk, alpha, Atiled[m + k * NT], b, tile_recv, b, 1, C_local, b);
                    }
                    else{
                        if(n%size == rank){
                            MPI_Send( Btiled[n*NT + k], b * b, MPI_DOUBLE,
                                          k%size, KT*n+k, comm );

                        }
                    }
                }
            MPI_Reduce(C_local, C_global, b*b, MPI_DOUBLE, mat_add_op, n%size, comm );
            
            
            free(C_local);
            if(rank == n%size)
                matrix_somme(C_global, Ctiled[n*NT+m], b); 
            
            free(C_global);
            }
    }
    
    tile2lapack( M,  N, b ,Ctiled , C,  ldc );

    tileFree(M, N, b, Atiled);
    tileFree(M, N, b, Btiled);
    tileFree(M, N, b, Ctiled);


    MPI_Bcast(C, M*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      
    return ALGONUM_SUCCESS;

}

/* To make sure we use the right prototype */
static dgemm_tiled_fct_t valid_dgemm_tiled_mpi __attribute__ ((unused)) = dgemm_tiled_mpi;

/* Declare the variable that will store the information about this version */
fct_list_t fct_dgemm_tiled_mpi;

/**
 * @brief Registration function
 */
void dgemm_tiled_mpi_init( void ) __attribute__( ( constructor ) );
void
dgemm_tiled_mpi_init( void )
{
    fct_dgemm_tiled_mpi.mpi    = 1;
    fct_dgemm_tiled_mpi.tiled  = 0;
    fct_dgemm_tiled_mpi.starpu = 0;
    fct_dgemm_tiled_mpi.name   = "mpi";
    fct_dgemm_tiled_mpi.helper = "MPI tiled implementation of the dgemm";
    fct_dgemm_tiled_mpi.fctptr = dgemm_dist_compute_centralize;
    fct_dgemm_tiled_mpi.next   = NULL;

    register_fct( &fct_dgemm_tiled_mpi, ALGO_GEMM );
}
