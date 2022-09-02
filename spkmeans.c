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
#include "yacobi_output.c"


#define EPSILON 0.00001
#define MAX_NUMBER_OF_ROTATIONS 100



int main(int argc, char **argv) {
    srand((int) time(NULL)); /* important for random */
    /* this is matrix 1 for yacobi */
    Matrix *A = create_matrix(5,5);
    matrix_set_entry(A, 0, 0, 0.9244561740566206);
    matrix_set_entry(A, 0, 1, 0.8742308098340758);
    matrix_set_entry(A, 0, 2, 0.8957630183913475);
    matrix_set_entry(A, 0, 3, 0.9448353953575863);
    matrix_set_entry(A, 0, 4, 0.2525962896380215);
    matrix_set_entry(A, 1, 0, 0.8742308098340758);
    matrix_set_entry(A, 1, 1, 0.47428491745836876);
    matrix_set_entry(A, 1, 2, 0.02753243942768835);
    matrix_set_entry(A, 1, 3, 0.30301127246390347);
    matrix_set_entry(A, 1, 4, 0.038962315397181646);
    matrix_set_entry(A, 2, 0, 0.8957630183913475);
    matrix_set_entry(A, 2, 1, 0.02753243942768835);
    matrix_set_entry(A, 2, 2, 0.27992567873413254);
    matrix_set_entry(A, 2, 3, 0.7802355856161365);
    matrix_set_entry(A, 2, 4, 0.7596530627838091);
    matrix_set_entry(A, 3, 0, 0.9448353953575863);
    matrix_set_entry(A, 3, 1, 0.30301127246390347);
    matrix_set_entry(A, 3, 2, 0.7802355856161365);
    matrix_set_entry(A, 3, 3, 0.985130302648835);
    matrix_set_entry(A, 3, 4, 0.13884538443846695);
    matrix_set_entry(A, 4, 0, 0.2525962896380215);
    matrix_set_entry(A, 4, 1, 0.038962315397181646);
    matrix_set_entry(A, 4, 2, 0.7596530627838091);
    matrix_set_entry(A, 4, 3, 0.13884538443846695);
    matrix_set_entry(A, 4, 4, 0.3869426064309628);
    print_matrix(A);
    YacobiOutput *yacobi_output = create_empty_yacobi_output();
    Jacobi(A, 3, yacobi_output);
    print_matrix(yacobi_output->V);
}


/* spkmeans functions */
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
    free(p1);
    free(p2);
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

    free_matrix(I);
    free_matrix(temp);
    free_matrix(X);
    return Lnorm;
}


/* JACOBI */
YacobiOutput *Jacobi(Matrix *A, int k, YacobiOutput *yacobi_output) {
    int dim = matrix_get_rows_num(A);
    Matrix *V = create_identity_matrix(dim);
    if(check_if_matrix_is_diagonal(A)) { /* edge case when is already diagonal */
        set_yacobi_output_values(yacobi_output, A, V, k);
        return yacobi_output;
    }

    int rotation_num = 0, need_to_stop = 0;
    double recent_off;
    S_and_C *s_and_c = create_empty_S_and_C();
    MaxElement *max_element = create_empty_max_element();

    while (rotation_num <= MAX_NUMBER_OF_ROTATIONS) {
        recent_off = matrix_off(A);
        matrix_get_non_diagonal_max_absolute_value(A, max_element);
        get_s_and_c_for_rotation_matrix(A, max_element, s_and_c);
        Matrix *P = create_identity_matrix(dim); build_rotation_matrix(s_and_c, max_element, dim, P); rotation_num ++;
        A = transform_matrix(A, s_and_c, max_element);
        V = multiply_matrices(V, P); free_matrix(P);
        if (matrix_converge(recent_off, A)) {
            break;
        }
    }
    set_yacobi_output_values(yacobi_output, A, V, k);
    free(max_element);
    free(s_and_c);
    return yacobi_output;
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

void build_rotation_matrix(S_and_C *s_and_c, MaxElement *max_element, int dim, Matrix *identity_matrix) {
    double s = s_and_c_get_s(s_and_c), c = s_and_c_get_c(s_and_c);
    int i = max_element_get_index1(max_element), j = max_element_get_index2(max_element);
    matrix_set_entry(identity_matrix, i, i, c);
    matrix_set_entry(identity_matrix, j, j, c);
    matrix_set_entry(identity_matrix, i, j, s);
    matrix_set_entry(identity_matrix, j, i, -s);
}

void normalize_matrix_rows(Matrix *matrix) {
    Point *row = (Point *)malloc(sizeof(Point));
    int num_of_rows, num_of_cols;
    double row_norm, entry;
    int i,j;
    num_of_rows = matrix_get_rows_num(matrix);
    num_of_cols = matrix_get_cols_num(matrix);

    for (i=0; i<num_of_rows; i++) {
        matrix_get_row_to_point(matrix, row, i);
        row_norm = euclidean_norm(row);
        divide_point_by_value(row, row_norm);
    }
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

Matrix *getU(Matrix *V, Matrix *A, int k) {
    int i, j, n = A->cols;
    double entry;
    Eigenvector *eigen_vectors_array = (Eigenvector *)malloc(sizeof(Eigenvector)*n);
    Eigenvector *current_eigen_vector;
    Matrix *U;
    Point *p;
    for (i=0; i<n; i++){
        p = &(eigen_vectors_array[i].point);
        matrix_get_column_to_point(V, p, i);
        eigen_vectors_array[i].eigenvalue = matrix_get_entry(A, i, i);
    }
    sort_eigenvectors_array(eigen_vectors_array, n);

    if (k == -1)
        k = get_k_from_sorted_eigenvectors_array(eigen_vectors_array, n);

    U = create_matrix(n, k);

    for (j=0; j<k; j++){
        p = &(eigen_vectors_array[j].point); /* the jth eigen vector, jth column */
        for (i=0; i<n; i++) {
            entry = point_get_entry(p, i);
            matrix_set_entry(U, i, j, entry);
        }
    }

    return U;
}

int get_k_from_sorted_eigenvectors_array(Eigenvector *eigen_vectors_array, int n) {
    double maxgap, currentgap;
    int k, i;
    maxgap = -1.;

    for (i=0; i<n/2; i++){
        currentgap = eigen_vectors_array[i].eigenvalue - eigen_vectors_array[i+1].eigenvalue;
        assert(currentgap >= 0);
        if (currentgap > maxgap){
            maxgap = currentgap;
            k = i;
        }
    }
    return k;
    
}

/* utilities */
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

