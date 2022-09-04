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
#include "yacobi_output.c"


#define EPSILON 0.00001
#define MAX_NUMBER_OF_ROTATIONS 100

int main() {
    int j, i, n;
    YacobiOutput *yacobi_output = create_empty_yacobi_output();
    Eigenvector *eigen_vectors_array;
    Matrix *X;
    char path[] = "testfiles/jacobi_0.txt";
    for (j=0; j<10; j++) {
        path[17] = j + '0';
        X = input_file_to_matrix(path);
        jacobi(X, yacobi_output);
        n = matrix_get_cols_num(yacobi_output->V);
        eigen_vectors_array = create_eigen_vectors_array(n);
        get_eigen_vectors_from_yacobi_output(yacobi_output, eigen_vectors_array);
        sort_eigenvectors_array(eigen_vectors_array, n);
        printf("%d \n", get_k_from_sorted_eigen_vectors_array(eigen_vectors_array, n));
    }
  
}

/*
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
    return 0;
    
}*/

void achieve_goal(Matrix *data_points, char *goal) {
    YacobiOutput *Jout;
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
        Jout = create_empty_yacobi_output();
        jacobi(data_points, Jout);
        print_jacobi_output(Jout);
        free_yacobi_output(Jout);
        return;
    }
    else {
        printf("Invalid Input!\n");
        exit(1);
    }
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
    int num_of_rows;
    double row_norm;
    int i;
    num_of_rows = matrix_get_rows_num(matrix);

    for (i=0; i<num_of_rows; i++) {
        matrix_get_row_to_point(matrix, row, i);
        row_norm = euclidean_norm(row);
        divide_point_by_value(row, row_norm);
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

Matrix *getU(YacobiOutput *yacobi_output, int k) { /* k == 0 if needed to be computed by eigengap heuristic */
    int eigenvectors_num, i, j, n;
    Eigenvector *eigen_vectors_array;
    Matrix *U;
    double entry;
    Point *point_j;

    eigenvectors_num = matrix_get_cols_num(yacobi_output->V);
    eigen_vectors_array = create_eigen_vectors_array(eigenvectors_num);
    get_eigen_vectors_from_yacobi_output(yacobi_output, eigen_vectors_array);
    sort_eigenvectors_array(eigen_vectors_array, eigenvectors_num);
    if (k == 0) 
        k = get_k_from_sorted_eigen_vectors_array(eigen_vectors_array, eigenvectors_num);

    n = matrix_get_rows_num(yacobi_output->V);
    U = create_matrix(n, k);
    for (j=0; j<k; j++){
        point_j = eigen_vectors_array[j].point; /*the jth eigen vector, jth column*/
        for (i=0; i<n; i++) {
            entry = point_get_entry(point_j, i);
            matrix_set_entry(U, i, j, entry);
        }
    }

    for (i=0; i<n; i++)
        free(eigen_vectors_array[i].point);
    free(eigen_vectors_array);
    return U;
}

int get_k_from_sorted_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n) {
    int k, i;
    double maxgap, currentgap;

    k = -1;
    maxgap = -1.;
    print_eigen_vectors_array(eigen_vectors_array, n);
    for (i=0; i<n-1; i++){ /* MIGHT BE AN ERROR, as it downs't work with n/2! */
        currentgap = eigen_vectors_array[i].eigen_value - eigen_vectors_array[i+1].eigen_value;
        assert(currentgap >= 0);
        printf("current gap: %f ", currentgap);
        if (currentgap > maxgap){
            maxgap = currentgap;
            k = i+1;
        }
    }
    return k;
}

void get_eigen_vectors_from_yacobi_output(YacobiOutput *yacobi_output, Eigenvector *eigen_vectors_array) {
    Matrix *A = yacobi_output_get_A(yacobi_output), *V = yacobi_output_get_V(yacobi_output);
    int i, n;
    n = matrix_get_cols_num(V);

    for (i=0; i<n; i++) {
        matrix_get_column_to_point(V, eigen_vectors_array[i].point, i);
        eigen_vectors_array[i].eigen_value = matrix_get_entry(A, i, i);
    }
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



/* SPKMEANS API */

Matrix *wam(Matrix* data_points){
    return create_weighted_matrix(data_points);
}

Matrix *ddg(Matrix* data_points){
    Matrix *W, *D;
    W = create_weighted_matrix(data_points);
    D = create_diagonal_degree_matrix(W);
    free_matrix(W);
    return D;
}

Matrix *lnorm(Matrix* data_points){
    Matrix *W, *D, *Lnorm;
    W = create_weighted_matrix(data_points);
    D = create_diagonal_degree_matrix(W);
    neg_root_to_diag_matrix(D);
    Lnorm = normalized_graph_laplacian(D, W);
    free_matrix(W);
    free_matrix(D);
    return Lnorm;
}

Matrix *matrix_tr(Matrix *m){
    Matrix *t = create_matrix(m->cols, m->rows);
    int i,j;
    for (i=0; i<t->rows; i++){
        for (j=0; j<t->cols; j++){
            matrix_set_entry(t, i, j, matrix_get_entry(m, j, i));
        }
    }
    return t;
}

YacobiOutput *jacobi(Matrix *A, YacobiOutput *yacobi_output) {
    int dim = matrix_get_rows_num(A);
    Matrix *V = create_identity_matrix(dim);
    if(_is_matrix_diag(A)) {
        set_yacobi_output_values(yacobi_output, A, V);
        return yacobi_output;
    }

    int rotation_num = 0;
    double recent_off;
    S_and_C *s_and_c = create_empty_S_and_C();
    MaxElement *max_element = create_empty_max_element();

    while (rotation_num < MAX_NUMBER_OF_ROTATIONS) {
        recent_off = matrix_off(A);
        matrix_get_non_diagonal_max_absolute_value(A, max_element);
        get_s_and_c_for_rotation_matrix(A, max_element, s_and_c);
        Matrix *P = create_identity_matrix(dim);
        build_rotation_matrix(s_and_c, max_element, dim, P); rotation_num ++;
        A = transform_matrix(A, s_and_c, max_element);
        V = multiply_matrices(V, P); free_matrix(P); /* TODO: here we lose the matrix reference */
        if (matrix_converge(recent_off, A)) {
            break;
        }
    }
    set_yacobi_output_values(yacobi_output, A, V);
    free(max_element);
    free(s_and_c);
    return yacobi_output;
}
