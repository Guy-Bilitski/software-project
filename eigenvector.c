#include "spkmeans.h"
#include <stdlib.h>
#ifndef EIGENVECTOR_IS_DEFINED
#define EIGENVECTOR_IS_DEFINED
typedef struct Eigenvector {
    Point *point;
    double eigen_value;
} Eigenvector;
#endif

/* Eigenvector API */
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

Eigenvector *create_eigen_vectors_array(int eigenvectors_num) {
    int i;
    Eigenvector *eigen_vectors_array;
    eigen_vectors_array = (Eigenvector *)malloc(sizeof(Eigenvector)*eigenvectors_num);
    for (i=0; i<eigenvectors_num; i++) 
        eigen_vectors_array[i].point = create_empty_point();
    return eigen_vectors_array;
}


/* Getters */
Point *eigen_vector_get_point(Eigenvector *eigen_vector) {
    return eigen_vector->point;
}

double eigen_vector_get_eigen_value(Eigenvector *eigen_vector) {
    return eigen_vector->eigen_value;
}

/* Utils */
int compare_eigenvectors(const void *p1, const void *p2) {
const Eigenvector *v1 = p1, *v2 = p2;
double diff = (v1->eigen_value) - (v2->eigen_value); /* if diff > 0 then v1 should be first, return -1 */
if (diff < 0) return 1;
if (diff > 0) return -1;
return 0;
}

void sort_eigenvectors_array(Eigenvector *array, size_t n) {
    qsort(array, n, sizeof(Eigenvector), compare_eigenvectors);
}

/* cleanup */
void free_eigen_vector(Eigenvector *eigen_vector) {
    free(eigen_vector_get_point(eigen_vector));
    free(eigen_vector);
}

void free_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n) {
    int i;
    for (i=0; i<n; i++) {
        free(eigen_vectors_array[i]);
    }
    free(eigen_vectors_array);
}

/* debugging function */
void print_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n) {
    int i,j;
    int m = eigen_vectors_array[0].point->dim;
    for (j=0; j<n; j++) {
        for (i=0; i<m; i++){
            /* printf("%.4f, ", point_get_entry(eigen_vectors_array[i].point, j)); */
            printf("%.4f, ", eigen_vectors_array[i].eigen_value);
        }
        printf("\n");
            
    }
}