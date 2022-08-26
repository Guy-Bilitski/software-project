#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> /* delete */
#include <assert.h>
#include "spkmeans.h"

#define true 1
#define false 0

#ifndef MATRIX_IS_DEFINED
#define MATRIX_IS_DEFINED
typedef struct Matrix {
    int rows;
    int cols;
    int is_diag;
    double *data;
} Matrix;
#endif

/* TODO: Check callocs/mallocs on failure */




/* Matrix API */
/* create */
Matrix *create_matrix(int rows, int cols) {
    int size_of_data;
    size_of_data = rows*cols;
    assert(rows>0 && cols>0);

    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->is_diag = false;
    matrix->data = (double *)calloc(sizeof(double), size_of_data);
    return matrix;
}

Matrix *create_diag_matrix(int n) {
    assert(n>0);
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    matrix->rows = n;
    matrix->cols = n;
    matrix->is_diag = true;
    matrix->data = (double *)calloc(sizeof(double), n);
    return matrix;
}

Matrix *create_identity_matrix(int n) {
    assert(n>0);
    int i;
    Matrix *I = create_matrix(n, n);
    for (i=0; i<n; i++)
        matrix_set_entry(I, i, i, 1.);
    return I;
}

/* getters */
int matrix_get_rows_num(Matrix *matrix) {
    return matrix->rows;
}

int matrix_get_cols_num(Matrix *matrix) {
    return matrix->cols;
}

double *matrix_get_data(Matrix *matrix) {
    return matrix->data;
}

double matrix_get_entry(Matrix *matrix, int row, int col){
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));
    double *matrix_data = matrix_get_data(matrix);

    if (_is_matrix_diag(matrix)){ /* if matrix is diagonal */
        if (row != col){
            return 0;
        } else {
            return matrix_data[col];
        }
        
    } else { /* matrix is NOT diagonal */
        return matrix_data[_get_matrix_index(matrix, row, col)];
    }
}

void matrix_get_row_to_point(Matrix *matrix, Point *point, int row_index) { /* TODO: what about diag matrix? */
    assert(matrix_get_rows_num(matrix) > row_index);
    int cols_num = matrix_get_cols_num(matrix);
    point->dim = cols_num;
    point->offset = 1;
    point->data = matrix_get_data(matrix) + row_index*cols_num;
}

void matrix_get_column_to_point(Matrix *matrix, Point *point, int column_index) {
    int cols_num = matrix_get_cols_num(matrix);
    assert(cols_num > column_index);
    point->dim = matrix_get_rows_num(matrix);
    point->data = matrix_get_data(matrix) + column_index;
    point->offset = cols_num;
}

void matrix_get_non_diagonal_max_absolute_value(Matrix *matrix, MaxElement *max_element) {
    max_element_set_new_values(max_element, matrix_get_entry(matrix, 0, 1), 0, 1);
    int i, j, cols_num = matrix_get_cols_num(matrix), rows_num = matrix_get_rows_num(matrix);
    double current_value;
    for (i=0; i<rows_num; i++) {
        for (j=0; j<cols_num; j++) {
            if (i != j) {
                current_value = abs(matrix_get_entry(matrix, i, j));
                if (current_value > max_element_get_value(max_element)) {
                    max_element_set_new_values(max_element, current_value, i, j);
                }
            }
        }
    }
}

/* setters */
void matrix_set_entry(Matrix *matrix, int row, int col, double value) {
    /* sets an entry in the matrix */
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));
    double *matrix_data = matrix_get_data(matrix);
    int matrix_index = _get_matrix_index(matrix, row, col);

    if (_is_matrix_diag(matrix)) {
        if (row != col && value != 0.) {
            _diag_to_square_matrix(matrix);
             matrix_data[matrix_index] = value;
        } else if (row != col && value == 0){ /* nothing to change */
            return;
        } else {
            matrix_data[row] = value;
        }
    } 
    else { /* matrix is NOT diagonal */
         matrix_data[matrix_index] = value;
    }
}

void matrix_set_row(Matrix *matrix, int row_index, Point *point) {
    /* TODO: what if is diag? */
    /* sets a whole point in the matrix */
    int i, rows_num = matrix_get_rows_num(matrix);
    if (!_is_matrix_diag(matrix)) {
        for (i=0; i<rows_num; i++) {
            matrix_set_entry(matrix, row_index, i, point_get_entry(point, i));
        }
    }
}

/* utilities */
int check_if_matrix_is_diagonal(Matrix *matrix) {
    int i, j;
    for (i=0; i<matrix_get_rows_num(matrix); i++) {
        for (j=0; j<matrix_get_cols_num(matrix); j++) {
            if (i != j && matrix_get_entry(matrix, i, j) != 0) {
                return false;
            }
        }
    }
    return true;
}

double matrix_get_row_sum(Matrix *matrix, int row_index) {
    int col_index, cols_num = matrix_get_cols_num(matrix);
    double sum;
    for (col_index=0; col_index<cols_num; col_index++) {
        sum += matrix_get_entry(matrix, row_index, col_index);
    }
    return sum;
}

void free_matrix(Matrix *matrix) {
    free(matrix_get_data(matrix));
    free(matrix);
}

void matrix_add_point_to_row(Matrix *matrix, int row_index, Point *point){
    int j;
    double entry;
    int dim = point_get_dim(point);
    assert(dim == matrix_get_cols_num(matrix));

    for (j=0; j<dim; j++){
        entry = matrix_get_entry(matrix, row_index, j) + point_get_entry(point, j);
        matrix_set_entry(matrix, row_index, j, entry);
    }
}

void reset_matrix_entries_to_zero(Matrix *matrix){
    int i,j;
    int rows = matrix_get_rows_num(matrix);
    int cols = matrix_get_cols_num(matrix);
    for (i=0; i<rows; i++)
        for (j=0; j<cols; j++)
            matrix_set_entry(matrix, i,j,0);
}

Matrix *multiply_matrices(Matrix *m1, Matrix *m2) {
    assert(matrix_get_cols_num(m1) == matrix_get_rows_num(m2));

    if (_is_matrix_diag(m1) && _is_matrix_diag(m2))
        return _multiply_matrices_diag_with_diag(m1, m2);

    else if (!_is_matrix_diag(m1) && !_is_matrix_diag(m2))
        return _multiply_matrices_nondiag_with_nondiag(m1, m2);
    
    else
        return _multiply_matrices_diag_with_nondiag(m1, m2);
}

/* Matrix inner functions */
int _is_matrix_diag(Matrix *matrix) {
    return matrix->is_diag;
}

int _get_matrix_index(Matrix *matrix, int row, int col) {
    if (!_is_matrix_diag(matrix))
        return row*(matrix->cols) + col;

}

void _diag_to_square_matrix(Matrix *matrix) { 
    /* Changing matrix->data from n sized array (diagonal matrix) to n*n sized array (square matrix) */
    int i,n = matrix_get_rows_num(matrix);
    double *old_data = matrix_get_data(matrix);
    matrix->data = (double *)calloc(sizeof(double), n*n);  /* isn't it safer to create a new matrix? */
    for (i=0; i<n; i++){
        (matrix->data)[i*n + i] = old_data[i];
    }
    matrix->is_diag = 0;
    free(old_data);
}

Matrix *_multiply_matrices_diag_with_diag(Matrix *m1, Matrix *m2) {
    int n = matrix_get_rows_num(m1);
    Matrix *new_matrix = create_diag_matrix(n);
    int i;
    double entry_value;

    for (i=0; i<n; i++){
        entry_value = matrix_get_entry(m1, i, i) * matrix_get_entry(m2, i, i);
        matrix_set_entry(new_matrix, i, i, entry_value);
    }
    return new_matrix;
}

Matrix *_multiply_matrices_diag_with_nondiag(Matrix *m1, Matrix *m2) {
    int i, j;
    double current_entry;
    int rows_num = matrix_get_rows_num(m1);
    int cols_num = matrix_get_cols_num(m2);
    Matrix *new_matrix = create_matrix(rows_num, cols_num);
    int first_is_diag = _is_matrix_diag(m1);

    for (i=0; i<rows_num; i++) {
        for (j=0; j<cols_num; j++) {
            if (first_is_diag)
                current_entry = matrix_get_entry(m1, i, i) * matrix_get_entry(m2, i, j);
            else
                current_entry = matrix_get_entry(m1, i, j) * matrix_get_entry(m2, j, j);
            matrix_set_entry(new_matrix, i, j, current_entry);
        }
    }
    return new_matrix;
}

Matrix *_multiply_matrices_nondiag_with_nondiag(Matrix *m1, Matrix *m2) {
    int i, j;
    int rows_num = matrix_get_rows_num(m1);
    int cols_num = matrix_get_cols_num(m2);
    Matrix *new_matrix = create_matrix(rows_num, cols_num);
    Point *row_point = (Point *)malloc(sizeof(Point));
    Point *col_point = (Point *)malloc(sizeof(Point));
    
    for (i=0; i<rows_num; i++) {
        matrix_get_row_to_point(m1, row_point, i);
        for (j=0; j<cols_num; j++) {
            matrix_get_column_to_point(m2, col_point, j);
            matrix_set_entry(new_matrix, i, j, inner_product(row_point, col_point));
        }
    }
    return new_matrix;
}

Matrix *sub_matrices(Matrix *A, Matrix *B) {
    assert(matrix_get_cols_num(A) == matrix_get_cols_num(B));
    assert(matrix_get_rows_num(A) == matrix_get_rows_num(B));
    int i,j;
    int cols = matrix_get_cols_num(A);
    int rows = matrix_get_rows_num(A);
    double value;
    Matrix *result = (Matrix *)malloc(sizeof(Matrix));

    
    for (i=0; i<rows; i++){
        for (j=0; j<cols; j++) {
            value = matrix_get_entry(A, i, j)-matrix_get_entry(B, i, j);
            matrix_set_entry(result, i, j, value);
        }
    }
    return result;
}


/* debugging function */
void print_matrix(Matrix *matrix) {
    space();
    int i, j;
    double val;
    for (i=0; i<matrix->rows; i++) {
        for (j=0; j<matrix->cols; j++) {
            val = matrix_get_entry(matrix, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
    space();
}

void print_matrix2(Matrix *matrix) {
    int i, j;
    double val;
    printf("[");
    for (i=0; i<matrix->rows; i++) {
        printf("[");
        for (j=0; j<matrix->cols; j++) {
            val = matrix_get_entry(matrix, i, j);
            printf("%lf", val);
            if (j != matrix->cols - 1)
                printf(", ");
        }
        if (i != matrix->rows - 1)
            printf("],\n");
        else printf("]]\n");
    }
    printf("\n");
}

void print_matrix_rows(Matrix *matrix) {
    int i;
    Point *point = (Point *)malloc(sizeof(Point));
    for (i=0; i<matrix_get_rows_num(matrix); i++) {
        matrix_get_row_to_point(matrix, point, i);
        print_point(point);
        printf("\n");
    }
    free(point);
}

void print_matrix_cols(Matrix *matrix) {
    int i;
    Point *point = (Point *)malloc(sizeof(Point));
    for (i=0; i<matrix_get_cols_num(matrix); i++) {
        matrix_get_column_to_point(matrix, point, i);
        print_point(point);
        printf("\n");
    }
    free(point);
}

double RandomReal(double low, double high) {
  double d;

  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}


Matrix *generate_matrix(int rows, int cols, int is_diag) {
    int i,j;
    Matrix *m;

    if (is_diag) {
        m = create_diag_matrix(rows);
        for (i=0; i<rows; i++) {
            matrix_set_entry(m, i, i, RandomReal(0, 10));
        }
    }
    else {
        m = create_matrix(rows, cols);
        for (i=0; i<rows; i++) {
            for (j=0; j<cols; j++) {
                matrix_set_entry(m, i, j, RandomReal(0, 10));
            }
        }
    }

    return m; 
}

Matrix *generate_symmetric_matrix(int n) {
    Matrix *new_matrix = generate_matrix(n, n, false);
    int i, j;
    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            matrix_set_entry(new_matrix, j, i, matrix_get_entry(new_matrix, i, j));
        }
    }
    return new_matrix;
}

void space() {
    printf("\n");
    printf("\n");
}


int main9() {
    srand((int) time(NULL)); /* important for random */
    Matrix *m1 = generate_matrix(5, 5, true);
    Matrix *m2 = generate_matrix(5, 5, true);
    Matrix *m3 = multiply_matrices(m1, m2);
    print_matrix2(m1);
    print_matrix2(m2);
    print_matrix2(m3);
    space();

    // Point *p = (Point *)malloc(sizeof(Point));
    // matrix_get_column_to_point(m1, p, 2);
    // print_point(p);
    // space();
    // matrix_get_row_to_point(m1, p, 1);
    // print_point(p);
    // space();
    // free(p);

    // print_matrix(m2);
    // space();
    
    // Matrix *product = multiply_matrices(m1, m2);

    // print_matrix_rows(product);
    // print_matrix_cols(product);
    
    // _clean_matrix(m1);
    // _clean_matrix(m2);
    // _clean_matrix(product);    

    return 1;
}
/*
gcc matrix.c
a.out
*/