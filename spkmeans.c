#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix.c"

typedef struct MaxElement
{
    int i;
    int j;
    double value;
} MaxElement;

typedef struct S_and_C
{
    double s;
    double c;
} S_and_C;


double gaussian_RBF(Point *x1, Point *x2);  /*computes w_i in the weighted adjacency matrix*/
Matrix *create_weighted_matrix(Matrix *X);  /* creates the weighted matrix */
Matrix *normalized_graph_laplacian(Matrix *D_minus_05, Matrix *W);
MaxElement get_off_diagonal_absolute_max(Matrix *matrix);
MaxElement get_off_diagonal_absolute_max(Matrix *matrix);





double gaussian_RBF(Point *p1, Point *p2) {
    double distance = euclidean_distance(p1, p2); /*TODO: CHECK*/
    return exp(-(distance / 2));
}

Matrix *create_weighted_matrix(Matrix *X) {
    int i, j, rows_num = matrix_get_rows_num(X);
    double value;
    Point *p1, *p2;
    Matrix *matrix = create_matrix(rows_num, rows_num);
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
    Matrix *diagonal_degree_matrix = create_diag_matrix(rows_num);
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


Matrix *normalized_graph_laplacian(Matrix *D_minus_05, Matrix *W) {
    int n = matrix_get_rows_num(W);
    Matrix *I = create_identity_matrix(n);
    Matrix *temp = multiply_matrices(D_minus_05, W);
    Matrix *X = multiply_matrices(X, D_minus_05);
    Matrix *Lnorm = sub_matrices(I, X);

    _free_matrix(I);
    _free_matrix(temp);
    _free_matrix(X);
    return Lnorm;
}


/* JACOBI */

MaxElement get_off_diagonal_absolute_max(Matrix *matrix){
    MaxElement max = {-1, -1, 0.};
    if (_is_matrix_diag(matrix))
        return max;

    int rows = matrix->rows;
    int cols = matrix->cols;
    int i,j;
    double current_value;

    for (i=0; i<rows; i++){
        for (j=i+1; j<cols; j++){
            current_value = abs(matrix_get_entry(matrix, i, j));
            if ( max.value <  current_value) {
                max.value = current_value;
                max.i = i;
                max.j = j;
            }
        }
    }
    return max;
}

S_and_C get_s_and_c_for_rotation_matrix(Matrix* A, MaxElement max) {
    double t, theta, s, c, sign;
    S_and_C result;

    theta = ( matrix_get_entry(A, max.j, max.j) - matrix_get_entry(A, max.i, max.i) ) / ( 2 * max.value );
    sign = (theta >= 0) ? 1:-1;
    t = sign / ( abs(theta) + sqrt( theta*theta + 1 ));
    c = 1.0 / sqrt( t*t + 1 );
    s = t*c;
    
    result.s = s;
    result.c = c;
    return result;
}



/*int argc, char **argv*/
int main() {
    Matrix *I5 = create_identity_matrix(5);
    print_matrix(I5);


    return 1;
}