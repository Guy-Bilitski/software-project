#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix.c"
#include "point.c"

double gaussian_RBF(Point *x1, Point *x2);
double euclidean_distance(Point *x1, Point *x2);



/*computes w_i in the weighted adjacency matrix*/
double gaussian_RBF(Point *p1, Point *p2) {
    double distance = euclidean_distance(p1, p2);
    return exp(-(distance / 2));
}


/*computes euclidean distance*/
double euclidean_distance(Point *p1, Point *p2) {
    double sum;
    int i, dim = point_get_dim(p1);
    for (i=0; i<dim; i++) {
        sum += pow(point_get_index(p1, i) - point_get_index(p2, i), 2);
    }
    return sqrt(sum);
}

Matrix *create_weighted_matrix(Matrix *X) {
    int i, j, rows_num = matrix_get_rows_num(X), cols_num = matrix_get_cols_num(X);
    double value;
    Point *p1, *p2;
    Matrix *matrix = create_matrix(rows_num, cols_num, false);
    print_matrix(matrix);
    for (i=0; i<rows_num; i++) {
        for (j=0; j<cols_num; j++) {
            if (i == j) {
                matrix_set(matrix, i, j, 0);
            }
            else {
                p1 = matrix_get_point(X, i);
                p2 = matrix_get_point(X, j);
                value = gaussian_RBF(p1, p2);
                matrix_set(matrix, i, j, value);
            }
        }
    }
    return matrix;
}


/*int argc, char **argv*/
int main2() {
    /*
    Point *point1 = create_empty_point(3);
    Point *point2 = create_empty_point(3);
    point_set_index(point1, 0, 33);
    point_set_index(point1, 1, 34);
    point_set_index(point1, 2, 36);
    point_set_index(point2, 0, 33);
    point_set_index(point2, 1, 34);
    point_set_index(point2, 2, 35);
    print_point(point1);
    printf("\n");
    print_point(point2);
    printf("\n");
    printf("%f", gaussian_RBF(point1, point2));
    free(point1);
    free(point2);*/



    printf("\n");
    Matrix *m = create_matrix(4,4,false);
    Matrix *m2;
    double *new_row = malloc(sizeof(double)*4);
    int i;
    for (i=0; i<4; i++) {
        new_row[i] = i + 100;
    }
    Point *p = create_point(4, new_row);
    matrix_set(m, 1, 2, 10);
    matrix_set(m, 0, 0, 20);
    matrix_set(m, 3, 3, 30);
    matrix_set(m, 3, 3, 5);
    matrix_set_point(m, 2, p);
    print_matrix(m);

    printf("\n");
    m2 = create_weighted_matrix(m);
    printf("\n");
    print_matrix(m2);
    free(m);
    free(m2);


    return 1;
}