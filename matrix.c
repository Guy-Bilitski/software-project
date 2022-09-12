#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> /* delete */
#include "spkmeans.h"

#define true 1
#define false 0


/* TODO: Check callocs/mallocs on failure */


/* Matrix API */
Matrix *create_matrix(int rows, int cols) {
    int size_of_data, i;
    Matrix *matrix;
    size_of_data = rows*cols;
    assert(rows>0 && cols>0);

    matrix = (Matrix *)malloc(sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->is_not_diag = 0;
    matrix->data = (double *)calloc(sizeof(double), size_of_data);
    for (i=0; i<size_of_data; i++)
        (matrix->data)[i] = 0.;
    return matrix;
}

Matrix *create_identity_matrix(int n) {
    int i;
    Matrix *I;

    assert(n>0);
    I = create_matrix(n,n);
    for (i=0; i<n; i++)
        matrix_set_entry(I, i, i, 1.);
    return I;
}

/* Getters */
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
    double *matrix_data;
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));
    
    matrix_data = matrix_get_data(matrix);
    return matrix_data[_get_matrix_index(matrix, row, col)];
}

void matrix_get_row_to_point(Matrix *matrix, Point *point, int row_index) {
    int cols_num;
    assert(matrix_get_rows_num(matrix) > row_index);
    
    cols_num = matrix_get_cols_num(matrix);
    point_set_values(point, cols_num, 1, matrix_get_data(matrix) + row_index*cols_num);
}

void matrix_get_column_to_point(Matrix *matrix, Point *point, int column_index) {
    int cols_num = matrix_get_cols_num(matrix);
    assert(cols_num > column_index);

    point_set_values(point, matrix_get_rows_num(matrix), cols_num, matrix_get_data(matrix) + column_index);
}

/* Setters */
void matrix_set_entry(Matrix *matrix, int row, int col, double value) {
    int matrix_index;
    double *matrix_data;
    double old_value;
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));

    matrix_data = matrix_get_data(matrix);
    matrix_index = _get_matrix_index(matrix, row, col);
    old_value = matrix_data[matrix_index];
    matrix_data[matrix_index] = value;

    if (row != col){
        if ((old_value == 0.) && (value != 0.))
            (matrix->is_not_diag)++;
        if ((old_value != 0.) && (value == 0.))
            (matrix->is_not_diag)--;
    }
}

/* Utils */
void matrix_get_non_diagonal_max_absolute_value(Matrix *matrix, MaxElement *max_element) {
    int i, j, n;
    double current_value;
    max_element_set_new_values(max_element, matrix_get_entry(matrix, 0, 1), 0, 1);
    n = matrix_get_cols_num(matrix);
    
    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            current_value = matrix_get_entry(matrix, i, j);
            if (fabs(current_value) > fabs(max_element_get_value(max_element))) {
                max_element_set_new_values(max_element, current_value, i, j);
            }
        }
    }
}

double matrix_get_row_sum(Matrix *matrix, int row_index) {
    int col_index, cols_num = matrix_get_cols_num(matrix);
    double sum;
    sum = 0;
    for (col_index=0; col_index<cols_num; col_index++) {
        sum += matrix_get_entry(matrix, row_index, col_index);
    }
    return sum;
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

Matrix *sub_matrices(Matrix *A, Matrix *B) {
    int i,j;
    int cols = matrix_get_cols_num(A);
    int rows = matrix_get_rows_num(A);
    double value;
    Matrix *result;
    assert(matrix_get_cols_num(A) == matrix_get_cols_num(B));
    assert(matrix_get_rows_num(A) == matrix_get_rows_num(B));
    result = create_matrix(rows, cols);

    
    for (i=0; i<rows; i++){
        for (j=0; j<cols; j++) {
            value = matrix_get_entry(A, i, j)-matrix_get_entry(B, i, j);
            matrix_set_entry(result, i, j, value);
        }
    }
    return result;
}

/* cleanup */
void free_matrix(Matrix *matrix) {
    free(matrix_get_data(matrix));
    free(matrix);
}

/* Matrix inner functions */
int _is_matrix_diag(Matrix *matrix) {
    return (!(matrix->is_not_diag)) && (matrix->cols == matrix->rows);
}

int _get_matrix_index(Matrix *matrix, int row, int col) {
    return row*(matrix->cols) + col;
}

Matrix *_multiply_matrices_diag_with_diag(Matrix *m1, Matrix *m2) {
    int n = matrix_get_rows_num(m1);
    Matrix *new_matrix = create_matrix(n,n);
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
    free(row_point);
    free(col_point);
    return new_matrix;
}

/* debugging function */
void print_matrix(Matrix *matrix) {
    int i, j;
    double val;
    int cols_num = matrix->cols;
    int rows_num = matrix->rows;

    for (i=0; i<rows_num; i++) {
        for (j=0; j<cols_num; j++) {
            val = matrix_get_entry(matrix, i, j);
            printf("%.4f", val);
            if (j < cols_num-1){
                printf(",");
            } else {
                printf("\n");
            }
        }

    }
}

void print_matrix2(Matrix *matrix) {
    int i, j;
    double val;
    printf("[");
    for (i=0; i<matrix->rows; i++) {
        printf("[");
        for (j=0; j<matrix->cols; j++) {
            val = matrix_get_entry(matrix, i, j);
            printf("%.4f,", val);
            if (j != matrix->cols - 1)
                printf(", ");
        }
        if (i != matrix->rows - 1)
            printf("],\n");
        else printf("]]\n");
    }
    printf("\n");
}

void print_matrix_diag(Matrix *matrix) {
    int i;
    double val;
    int n;
    assert(matrix_get_cols_num(matrix) == matrix_get_cols_num(matrix));
    n = matrix->cols;

    for (i=0; i<n; i++) {
        val = matrix_get_entry(matrix, i, i);
        printf("%.4f", val);
        if (i < n-1) printf(",");
    }
    printf("\n");
}

void print_matrix_rows(Matrix *matrix) {
    int i;
    Point *point = create_empty_point();
    for (i=0; i<matrix_get_rows_num(matrix); i++) {
        matrix_get_row_to_point(matrix, point, i);
        print_point(point);
        printf("\n");
    }
    free(point);
}

void print_matrix_cols(Matrix *matrix) {
    int i;
    Point *point = create_empty_point();
    for (i=0; i<matrix_get_cols_num(matrix); i++) {
        matrix_get_column_to_point(matrix, point, i);
        print_point(point);
        printf("\n");
    }
    free(point);
}

double RandomReal(double low, double high) /* TODO: delete */ {
  double d;

  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}

Matrix *generate_matrix(int rows, int cols, int is_diag) {
    int i,j;
    Matrix *m;

    if (is_diag) {
        m = create_matrix(rows, rows);
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
