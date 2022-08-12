#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix.c"

double gaussian_RBF(Point *x1, Point *x2);  /*computes w_i in the weighted adjacency matrix*/
Matrix *create_weighted_matrix(Matrix *X);  /* creates the weighted matrix */


double gaussian_RBF(Point *p1, Point *p2) {
    double distance = euclidean_distance(p1, p2);
    return exp(-(distance / 2));
}

Matrix *create_weighted_matrix(Matrix *X) {
    int i, j, rows_num = matrix_get_rows_num(X);
    double value;
    Point *p1, *p2;
    Matrix *matrix = create_matrix(rows_num, rows_num, false);
    for (i=0; i<rows_num; i++) {
        for (j=0; j<rows_num; j++) {
            if (i == j) {
                matrix_set_entry(matrix, i, j, 0);
            }
            else {
                p1 = matrix_get_row(X, i);
                p2 = matrix_get_row(X, j);
                value = gaussian_RBF(p1, p2);
                matrix_set_entry(matrix, i, j, value);
            }
        }
    }
    return matrix;
}

Matrix *create_diagonal_degree_matrix(Matrix *matrix) {
    int i, rows_num = matrix_get_rows_num(matrix);
    double value;
    Matrix *diagonal_degree_matrix = create_matrix(rows_num, rows_num, true);
    for (i=0; i<rows_num; i++) {
        value = matrix_get_row_sum(matrix, i);
        matrix_set_entry(diagonal_degree_matrix, i, i, value);
    }
    return diagonal_degree_matrix;
}

void neg_root_to_diag_matrix(Matrix *matrix) {
    assert(_is_matrix_diag(matrix));
    int i, rows_num = matrix_get_rows_num(matrix);
    double value;
    for (i=0; i<rows_num; i++) {
        value = matrix_get_entry(matrix, i, i);
        assert(value != 0);
        value = 1 / sqrt(value);
        matrix_set_entry(matrix, i, i, value);
    }
}


/*int argc, char **argv*/
int main() {
    Matrix *X = generate_matrix(6,3, false);
    Matrix *W = create_weighted_matrix(X);

    print_matrix(X);
    space();
    print_matrix(W);
    space();

    printf("%lx", sizeof(W->data));
    space();

    Matrix *D = create_diagonal_degree_matrix(W);
    print_matrix(D);
    space();
    printf("%lx", sizeof(D->data));
    space();
    neg_root_to_diag_matrix(D);
    print_matrix(D);
    space();
    Point *p = matrix_get_row(D, 3);
    print_point(p);



    _clean_matrix(X);
    _clean_matrix(W);
    _clean_matrix(D);
    


    return 1;
}