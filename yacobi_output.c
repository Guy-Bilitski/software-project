#include "spkmeans.h"
#include <stdlib.h>
#ifndef YACOBI_OUTPUT_IS_DEFINED
#define YACOBI_OUTPUT_IS_DEFINED
typedef struct YacobiOutput
{
    Matrix *A;
    Matrix *V;
} YacobiOutput;
#endif


/* Yacobi Output API */
YacobiOutput *create_empty_yacobi_output() {
    YacobiOutput *yacobi_output = (YacobiOutput *)malloc(sizeof(YacobiOutput));
    return yacobi_output;
}

/* Getters */
Matrix *yacobi_output_get_A(YacobiOutput *yacobi_output) {
    return yacobi_output->A;
}

Matrix *yacobi_output_get_V(YacobiOutput *yacobi_output) {
    return yacobi_output->V;
}

/* Setters */
void set_yacobi_output_values(YacobiOutput *yacobi_output, Matrix *A, Matrix *V) {
    yacobi_output->A = A;
    yacobi_output->V = V;
}

/* cleanup */
void free_yacobi_output(YacobiOutput *yacobi_output) {
    free(yacobi_output->A);
    free(yacobi_output->V);
    free(yacobi_output);
}

/* debugging */
void print_jacobi_output(YacobiOutput *J) {
    print_matrix_diag(J->A);
    print_matrix(J->V);
}