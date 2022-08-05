#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix.c"

double gaussian_RBF(double *x1, double *x2);
double euclidean_distance(double *x1, double *x2);



/*computes w_i in the weighted adjacency matrix*/
double gaussian_RBF(double *x1, double *x2) {
    return 0.0;
}


/*computes euclidean norm*/
double euclidean_distance(double *x1, double *x2) {
    return 0.0;
}


/*int argc, char **argv*/
int main2() {
    Point *point = create_point(3);
    point_set(point, 0, 100);
    point_set(point, 1, 1);
    point_set(point, 2, 8);
    print_point(point);
    free(point);


    return 1;
}