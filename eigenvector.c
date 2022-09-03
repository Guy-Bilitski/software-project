#include "spkmeans.h"
#include <stdlib.h>
#ifndef EIGENVECTOR_IS_DEFINED
#define EIGENVECTOR_IS_DEFINED
typedef struct Eigenvector {
    Point *point;
    double eigen_value;
} Eigenvector;
#endif

Eigenvector *create_empty_eigen_vector() {
    Eigenvector *eigen_vector = (Eigenvector *)malloc(sizeof(Eigenvector));
    return eigen_vector;
}


Eigenvector *create_eigen_vector(Point *point, double eigen_value) {
    Eigenvector *eigen_vector = create_empty_eigen_vector();
    eigen_vector->point = point;
    eigen_vector->eigen_value = eigen_value;
    return eigen_vector;
}

Point *eigen_vector_get_point(Eigenvector *eigen_vector) {
    return eigen_vector->point;
}

double eigen_vector_get_eigen_value(Eigenvector *eigen_vector) {
    return eigen_vector->eigen_value;
}

int compare_eigenvectors(const void *p1, const void *p2)
{
const Eigenvector *v1 = p1, *v2 = p2;
double diff = (v1->eigen_value) - (v2->eigen_value); /* if diff > 0 then v1 should be first, return -1 */
if (diff < 0) return 1;
if (diff > 0) return -1;
return 0;
}

void sort_eigenvectors_array(Eigenvector *array, size_t n) {
    qsort(array, n, sizeof(Eigenvector), compare_eigenvectors);
}

void print_eigen_vectors_array(Eigenvector **eigen_vectors_array, int n) {
    int i;
    for (i=0; i<n; i++) {
        print_point(eigen_vector_get_point(eigen_vectors_array[i]));
    }
}

void free_eigen_vector(Eigenvector *eigen_vector) {
    free(eigen_vector_get_point(eigen_vector));
    free(eigen_vector);
}

void free_eigen_vectors_array(Eigenvector **eigen_vectors_array, int n) {
    int i;
    for (i=0; i<n; i++) {
        free_eigen_vector(eigen_vectors_array[i]);
    }
}