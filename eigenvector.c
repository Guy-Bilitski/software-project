#include "spkmeans.h"
#include <stdlib.h>
#ifndef EIGENVECTOR_IS_DEFINED
#define EIGENVECTOR_IS_DEFINED
typedef struct Eigenvector {
    Point point;
    double eigenvalue;
} Eigenvector;
#endif



int compare_eigenvectors(const void *p1, const void *p2)
{
const Eigenvector *v1 = p1, *v2 = p2;
double diff = (v1->eigenvalue) - (v2->eigenvalue); /* if diff > 0 then v1 should be first, return -1 */
if (diff < 0) return 1;
if (diff > 0) return -1;
return 0;
}

void sort_eigenvectors_array(Eigenvector *array, size_t n) {
    qsort(array, n, sizeof(Eigenvector), compare_eigenvectors);
}