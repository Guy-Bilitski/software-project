#include <float.h>
#include <math.h>
#include <ctype.h>
#include "spkmeans.h"
#include "s_and_c.c"
#include "max_element.c"
#include "matrix.c"
#include "eigenvector.c"
#include "point.c"
#include "kmeans_io.c"
#include "kmeans.c"
#include "jacobi_output.c"

#define EPSILON 0.00001
#define MAX_NUMBER_OF_ROTATIONS 100


int main (int argc, char **argv) {
    const char *input_filename;
    char *goal;
    Matrix *data_points;

    if (argc != 3) {
        printf("Invalid Input!\n");
        exit(1);
    }
    
    goal = argv[1];
    input_filename = argv[2];

    data_points = input_file_to_matrix(input_filename);

    achieve_goal(data_points, goal);
    free_matrix(data_points);
    return 0;
    
}

/* spkmeans functions */
void achieve_goal(Matrix *data_points, char *goal) {
    JacobiOutput *Jout;
    if (!strcmp(goal, "wam")){
        Matrix *W = wam(data_points);
        print_matrix(W);
        free_matrix(W);
        
        return;
    }
    else if (!strcmp(goal, "ddg")){
        Matrix *D = ddg(data_points);
        print_matrix(D);
        free_matrix(D);
        return;
    }
    else if (!strcmp(goal, "lnorm")){
        Matrix *Lnorm = lnorm(data_points);
        print_matrix(Lnorm);
        free_matrix(Lnorm);
        return;
    }
    else if (!strcmp(goal, "jacobi")){
        Jout = create_empty_jacobi_output();
        jacobi(data_points, Jout);
        print_jacobi_output(Jout);
        free_jacobi_output(Jout);
        return;
    }
    else {
        printf("Invalid Input!\n");
        exit(1);
    }
}

double gaussian_RBF(Point *p1, Point *p2) {
    double distance = euclidean_distance(p1, p2); /*TODO: CHECK*/
    return exp(-(distance / 2));
}

Matrix *create_weighted_matrix(Matrix *X) {
    int i, j, rows_num;
    double value;
    Point *p1, *p2;
    Matrix *matrix;
    rows_num = matrix_get_rows_num(X);
    matrix = create_matrix(rows_num, rows_num);
    p1 = create_empty_point(); p2 = create_empty_point();
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
    free(p1);
    free(p2);
    return matrix;
}

Matrix *create_diagonal_degree_matrix(Matrix *matrix) {
    int i, rows_num = matrix_get_rows_num(matrix);
    double value;
    Matrix *diagonal_degree_matrix = create_matrix(rows_num, rows_num);
    for (i=0; i<rows_num; i++) {
        value = matrix_get_row_sum(matrix, i);
        matrix_set_entry(diagonal_degree_matrix, i, i, value);
    }
    return diagonal_degree_matrix;
}

void neg_root_to_diag_matrix(Matrix *D) {
    int i, rows_num;
    double value;
    assert(_is_matrix_diag(D));
    rows_num = matrix_get_rows_num(D);

    for (i=0; i<rows_num; i++) {
        value = matrix_get_entry(D, i, i);
        assert(value != 0);
        value = 1 / sqrt(value);
        matrix_set_entry(D, i, i, value);
    }
}

Matrix *normalized_graph_laplacian(Matrix *D_minus_05, Matrix *W) {
    int n = matrix_get_rows_num(W);
    Matrix *I = create_identity_matrix(n);
    Matrix *temp = multiply_matrices(D_minus_05, W);
    Matrix *X = multiply_matrices(temp, D_minus_05);
    Matrix *Lnorm = sub_matrices(I, X);

    free_matrix(I);
    free_matrix(temp);
    free_matrix(X);
    return Lnorm;
}

void get_s_and_c_for_rotation_matrix(Matrix* A, MaxElement *max_element, S_and_C *s_and_c) {
    double t, theta, s, c, sign, value = max_element_get_value(max_element);
    int i = max_element_get_index1(max_element), j = max_element_get_index2(max_element);

    theta = ( matrix_get_entry(A, j, j) - matrix_get_entry(A, i, i) ) / ( 2 * value );
    sign = (theta >= 0) ? 1 : -1;
    t = sign / ( fabs(theta) + sqrt( theta*theta + 1 ));
    c = 1.0 / sqrt( t*t + 1 );
    s = t*c;
    
    S_and_C_set_values(s_and_c, s, c);
}

void build_rotation_matrix(S_and_C *s_and_c, MaxElement *max_element, Matrix *identity_matrix) {
    double s = s_and_c_get_s(s_and_c), c = s_and_c_get_c(s_and_c);
    int i = max_element_get_index1(max_element), j = max_element_get_index2(max_element);
    matrix_set_entry(identity_matrix, i, i, c);
    matrix_set_entry(identity_matrix, j, j, c);
    matrix_set_entry(identity_matrix, i, j, s);
    matrix_set_entry(identity_matrix, j, i, -s);
}

void normalize_matrix_rows(Matrix *matrix) { /* TODO: what if the row is 0? */
    Point *row = create_empty_point();
    int num_of_rows, i;
    double row_norm;
    num_of_rows = matrix_get_rows_num(matrix);

    for (i=0; i<num_of_rows; i++) {
        matrix_get_row_to_point(matrix, row, i);
        row_norm = euclidean_norm(row);
        if (row_norm != 0) {
            divide_point_by_value(row, row_norm);
        }
    }
    free(row);
}

double matrix_off(Matrix *matrix) {
    int i, j, rows_num = matrix_get_rows_num(matrix), cols_num = matrix_get_cols_num(matrix);
    double sum = 0;
    for (i=0; i<rows_num; i++) {
        for (j=0; j<cols_num; j++) {
            if (i != j) {
                sum += pow(matrix_get_entry(matrix, i, j), 2);
            }
        }
    }
    return sum;
}

Matrix *transform_matrix(Matrix *matrix, S_and_C *s_and_c, MaxElement *max_element) {
    int row_index, col_index, rows_num = matrix_get_rows_num(matrix), cols_num = matrix_get_cols_num(matrix), i = max_element_get_index1(max_element), j = max_element_get_index2(max_element);
    double matrix_entry, s = s_and_c_get_s(s_and_c), c = s_and_c_get_c(s_and_c);
    Matrix *transformed_matrix = create_matrix(rows_num, rows_num);
    for (row_index=0; row_index<rows_num; row_index++) {
        for (col_index=0; col_index<cols_num; col_index++) {
            matrix_entry = get_value_for_transformed_matrix(matrix, s, c, i, j, row_index, col_index);
            matrix_set_entry(transformed_matrix, row_index, col_index, matrix_entry);
        }
    }
    return transformed_matrix;
}

int get_k_from_sorted_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n) {
    int k, i;
    double maxgap, currentgap;

    k = -1;
    maxgap = -1.;
    for (i=0; i<n/2; i++){ /* MIGHT BE AN ERROR, as it downs't work with n/2! */
        currentgap = eigen_vectors_array[i].eigen_value - eigen_vectors_array[i+1].eigen_value;
        assert(currentgap >= 0);
        if (currentgap > maxgap){
            maxgap = currentgap;
            k = i+1;
        }
    }
    return k;
}

void get_eigen_vectors_from_jacobi_output(JacobiOutput *jacobi_output, Eigenvector *eigen_vectors_array) {
    Matrix *A = jacobi_output_get_A(jacobi_output), *V = jacobi_output_get_V(jacobi_output);
    int i, n;
    n = matrix_get_cols_num(V);

    for (i=0; i<n; i++) {
        matrix_get_column_to_point(V, eigen_vectors_array[i].point, i);
        eigen_vectors_array[i].eigen_value = matrix_get_entry(A, i, i);
    }
}

Matrix *getU(JacobiOutput *jacobi_output, int k) { /* k == 0 if needed to be computed by eigengap heuristic */
    int eigenvectors_num, i, j, n;
    Eigenvector *eigen_vectors_array;
    Matrix *U, *V = jacobi_output_get_V(jacobi_output);
    double entry;
    Point *point_j;

    eigenvectors_num = matrix_get_cols_num(V);
    eigen_vectors_array = create_eigen_vectors_array(eigenvectors_num);
    get_eigen_vectors_from_jacobi_output(jacobi_output, eigen_vectors_array);
    sort_eigenvectors_array(eigen_vectors_array, eigenvectors_num);
    if (k == 0) {
        k = get_k_from_sorted_eigen_vectors_array(eigen_vectors_array, eigenvectors_num);
        /*
        print_matrix_diag(jacobi_output->A);
        printf("cols: %d, rows: %d\n", jacobi_output->A->cols, jacobi_output->A->rows);
        printf("K is %d\n",k);
        */
    }
    n = matrix_get_rows_num(V);
    U = create_matrix(n, k);
    for (j=0; j<k; j++){
        point_j = eigen_vectors_array[j].point; /*the jth eigen vector, jth column*/
        for (i=0; i<n; i++) {
            entry = point_get_entry(point_j, i);
            matrix_set_entry(U, i, j, entry);
        }
    }

    for (i=0; i<n; i++) {
        free(eigen_vectors_array[i].point);
    }
    free(eigen_vectors_array);
    return U;
}

double get_value_for_transformed_matrix(Matrix *old_matrix, double s, double c, int i, int j, int row_index, int col_index) {
     if (col_index == i) {
        if (row_index == i) {
            return pow(c, 2)*matrix_get_entry(old_matrix, i, i) + pow(s, 2)*matrix_get_entry(old_matrix, j, j) - 2*s*c*matrix_get_entry(old_matrix, i ,j);
        } else if (row_index == j) {
            return 0;
        } else {
            return c*matrix_get_entry(old_matrix, row_index, i) - s*matrix_get_entry(old_matrix, row_index, j);
        }
    } else if (col_index == j) {
        if (row_index == j) {
            return pow(s, 2)*matrix_get_entry(old_matrix, i, i) + pow(c, 2)*matrix_get_entry(old_matrix, j, j) + 2*s*c*matrix_get_entry(old_matrix, i, j);
        } else if (row_index == i) {
            return 0;
        } else {
            return c*matrix_get_entry(old_matrix, row_index, j) + s*matrix_get_entry(old_matrix, row_index, i);
        }
    } else if (row_index == i) { /* edge cases have already been validated */
        return c*matrix_get_entry(old_matrix, i, col_index) - s*matrix_get_entry(old_matrix, j, col_index);
    } else if (row_index == j) {
        return c*matrix_get_entry(old_matrix, j, col_index) + s*matrix_get_entry(old_matrix, i, col_index);
    } else {
        return matrix_get_entry(old_matrix, row_index, col_index); 
    }
}

int matrix_converge(double A_off, Matrix *A) {
    return A_off - matrix_off(A) <= EPSILON;
}

/* SPKMEANS API */
Matrix *wam(Matrix* data_points) {
    return create_weighted_matrix(data_points);
}

Matrix *ddg(Matrix* data_points) {
    Matrix *W, *D;
    W = create_weighted_matrix(data_points);
    D = create_diagonal_degree_matrix(W);
    free_matrix(W);
    return D;
}

Matrix *lnorm(Matrix* data_points) {
    Matrix *W, *D, *Lnorm;
    W = create_weighted_matrix(data_points);
    D = create_diagonal_degree_matrix(W);
    neg_root_to_diag_matrix(D);
    Lnorm = normalized_graph_laplacian(D, W);
    free_matrix(W);
    free_matrix(D);
    return Lnorm;
}

JacobiOutput *jacobi(Matrix *A, JacobiOutput *jacobi_output) { /* notice that jacobi function frees Matrix A! */
    int dim;
    int rotation_num;
    double recent_off;
    Matrix *A_tmp, *V_tmp, *V;
    S_and_C *s_and_c;
    MaxElement *max_element;
    Matrix *P;

    dim = matrix_get_rows_num(A);
    V = create_identity_matrix(dim);
    if(_is_matrix_diag(A)) {
        set_jacobi_output_values(jacobi_output, A, V);
        return jacobi_output;
    }

    rotation_num = 0;
    s_and_c = create_empty_S_and_C();
    max_element = create_empty_max_element();
    
    while (rotation_num < MAX_NUMBER_OF_ROTATIONS) {
        recent_off = matrix_off(A);
        matrix_get_non_diagonal_max_absolute_value(A, max_element);
        get_s_and_c_for_rotation_matrix(A, max_element, s_and_c);
        P = create_identity_matrix(dim);
        build_rotation_matrix(s_and_c, max_element, P); rotation_num ++;
        A_tmp = transform_matrix(A, s_and_c, max_element); free_matrix(A); A=A_tmp;
        V_tmp = multiply_matrices(V, P); free(V); V = V_tmp; free_matrix(P); 
        if (matrix_converge(recent_off, A)) {
            break;
        }
    }
    set_jacobi_output_values(jacobi_output, A, V);
    free(max_element);
    free(s_and_c);
    return jacobi_output;
}
