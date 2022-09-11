#include "spkmeans.h"
#include <stdlib.h>
#ifndef JACOBI_OUTPUT_IS_DEFINED
#define JACOBI_OUTPUT_IS_DEFINED
typedef struct JacobiOutput
{
    Matrix *A;
    Matrix *V;
} JacobiOutput;
#endif


/* Jacobi Output API */
JacobiOutput *create_empty_jacobi_output() {
    JacobiOutput *jacobi_output = (JacobiOutput *)malloc(sizeof(JacobiOutput));
    return jacobi_output;
}

/* Getters */
Matrix *jacobi_output_get_A(JacobiOutput *jacobi_output) {
    return jacobi_output->A;
}

Matrix *jacobi_output_get_V(JacobiOutput *jacobi_output) {
    return jacobi_output->V;
}

/* Setters */
void set_jacobi_output_values(JacobiOutput *jacobi_output, Matrix *A, Matrix *V) {
    jacobi_output->A = A;
    jacobi_output->V = V;
}

/* cleanup */
void free_jacobi_output(JacobiOutput *jacobi_output) {
    free_matrix(jacobi_output->A);
    free_matrix(jacobi_output->V);
    free(jacobi_output);
}

/* debugging */
void print_jacobi_output(JacobiOutput *J) {
    print_matrix_diag(J->A);
    print_matrix(J->V);
}