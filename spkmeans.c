#include <float.h>
#include <math.h>
#include <ctype.h>
#include <s_and_c.h>
#include <max_element.h>
#include "matrix.c"


double gaussian_RBF(Point *x1, Point *x2);  /*computes w_i in the weighted adjacency matrix*/
Matrix *create_weighted_matrix(Matrix *X);  /* creates the weighted matrix */
Matrix *create_diagonal_degree_matrix(Matrix *matrix);
void neg_root_to_diag_matrix(Matrix *matrix);
Matrix *normalized_graph_laplacian(Matrix *D_minus_05, Matrix *W);
MaxElement get_off_diagonal_absolute_max(Matrix *matrix);
S_and_C get_s_and_c_for_rotation_matrix(Matrix* A, MaxElement max);
Matrix *build_rotation_matrix(S_and_C s_and_c, MaxElement max_element, int dim); /* returns the rotation matrix p */
void normalize_matrix_rows(Matrix *matrix);


double gaussian_RBF(Point *p1, Point *p2) {
    double distance = euclidean_distance(p1, p2); /*TODO: CHECK*/
    return exp(-(distance / 2));
}

Matrix *create_weighted_matrix(Matrix *X) {
    int i, j, rows_num = matrix_get_rows_num(X);
    double value;
    Point *p1 = (Point *)malloc(sizeof(Point)), *p2 = (Point *)malloc(sizeof(Point));
    Matrix *matrix = create_matrix(rows_num, rows_num);
    for (i=0; i<rows_num; i++) {
        for (j=0; j<rows_num; j++) {
            if (i == j) {
                matrix_set_entry(matrix, i, j, 0);
            }
            else {
                matrix_get_row_to_point(X, p1, i);
                matrix_get_row_to_point(X, p2, j);
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

Matrix *build_rotation_matrix(S_and_C s_and_c, MaxElement max_element, int dim) {
    Matrix *p = create_identity_matrix(dim);
    int s = s_and_c_get_s(s_and_c), c = s_and_c_get_c(s_and_c), i = max_element_get_index1(max_element), j = max_element_get_index2(max_element);
    matrix_set_entry(p, i, i, c);
    matrix_set_entry(p, j, j, c);
    matrix_set_entry(p, i, j, s);
    matrix_set_entry(p, j, i, -s);
    return p;
}

void normalize_matrix_rows(Matrix *matrix) {
    Point *row = (Point *)malloc(sizeof(Point));
    int num_of_rows, num_of_cols;
    double rows_norm, entry;
    int i,j;
    num_of_rows = matrix_get_rows_num(matrix);
    num_of_cols = matrix_get_cols_num(matrix);

    for (i=0; i<num_of_rows; i++) {
        matrix_get_row_to_point(matrix, row, i);
        rows_norm = euclidean_norm(row);

        for (j=0; j<num_of_cols; j++) {
            entry = matrix_get_entry(matrix, i, j) / rows_norm;
            matrix_set_entry(matrix, i, j, entry);
        }
    }
}


/*int argc, char **argv*/
int main() {


    return 1;
}