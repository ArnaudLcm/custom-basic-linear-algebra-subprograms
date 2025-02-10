#include "myblas.h"

void affiche(int M, int N, double *A, int lda, FILE *flux)
{
    int m, n, k;
    k = (M > N) ? N : M;

    fprintf(flux, "\n");

    for (n = 0; n < N; n++)
    {
        for (m = 0; m < M; m++)
        {
            fprintf(flux, "%lf ", A[n * lda + m]);
        }
        fprintf(flux, "\n");
    }
    fprintf(flux, "\n");
}