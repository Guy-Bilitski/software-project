#include "spkmeans.h"
#include <stdlib.h>
#ifndef YACOBI_OUTPUT_IS_DEFINED
#define YACOBI_OUTPUT_IS_DEFINED
typedef struct YacobiOutput
{
    Matrix *A;
    Matrix *V;
    int k;
} YacobiOutput;
#endif


YacobiOutput *create_empty_yacobi_output() {
    YacobiOutput *yacobi_output = (YacobiOutput *)malloc(sizeof(YacobiOutput));
    return yacobi_output;
}

void set_yacobi_output_values(YacobiOutput *yacobi_output, Matrix *A, Matrix *V, int k) {
    yacobi_output->A = A;
    yacobi_output->V = V;
    yacobi_output->k = k;
}

void free_yacobi_output(YacobiOutput *yacobi_output) {
    free(yacobi_output->A);
    free(yacobi_output->V);
    free(yacobi_output);
}