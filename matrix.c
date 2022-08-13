#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> /* delete */
#include <assert.h>
#include "point.c"

#define true 1
#define false 0


typedef struct Matrix {
    int rows;
    int cols;
    int is_diag;
    double *data;
} Matrix;


/* TODO: Check callocs/mallocs on failure */

/* Matrix API */
Matrix *create_matrix(int rows, int cols);
Matrix *create_diag_matrix(int n);  /* creates a matrix with dimensions rows X columns all zeros */

int matrix_get_rows_num(Matrix *matrix);  /* returns the number of rows in the matrix */
int matrix_get_cols_num(Matrix *matrix);  /* returns the number of columns in the matrix */
double *matrix_get_data(Matrix *matrix);  /* returns the matrix data */
double matrix_get_row_sum(Matrix *matrix, int row_index);  /* returns <row_index> row sum of values */

double matrix_get_entry(Matrix *matrix, int row, int col);  /* returns the (row, col) entry */
Point *matrix_get_row(Matrix *matrix, int row_index);  /* returns the <row_index> point (row) of the matrix */
Point *matrix_get_column(Matrix *matrix, int column_index); /* returns the <column_index> column of the matrix */

void matrix_get_row_to_point(Matrix *matrix, Point *point, int row_index);
void matrix_get_column_to_point(Matrix *matrix, Point *point, int column_index);
void matrix_set_entry(Matrix *matrix, int row, int col, double value);  /* sets <value> in (row, col) entry */
void matrix_set_row(Matrix *matrix, int row_index, Point *point);  /* sets a point (row) in the matrix in <row_index> */
void matrix_add_point_to_row(Matrix *matrix, int row_index, Point *point);
Matrix *multiply_matrices(Matrix *m1, Matrix *m2);  /* multiply m1 X m2 and returns the new matrix */
void reset_matrix_entries_to_zero(Matrix *matrix);

/* Matrix inner functions */
int _is_matrix_diag(Matrix *matrix);
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);
void _clean_matrix(Matrix *matrix);
Matrix *_multiply_matrices_diag_with_diag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_diag_with_nondiag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_nondiag_with_nondiag(Matrix *m1, Matrix *m2);

/* debugging functions */
void print_matrix(Matrix *matrix);
void space();


/* Matrix API */
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

int matrix_get_rows_num(Matrix *matrix) {
    return matrix->rows;
}

int matrix_get_cols_num(Matrix *matrix) {
    return matrix->cols;
}

double *matrix_get_data(Matrix *matrix) {
    return matrix->data;
}

double matrix_get_row_sum(Matrix *matrix, int row_index) {
    Point *point = matrix_get_row(matrix, row_index);
    return sum_point_values(point);
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

/* TODO: FIX DATA LEAK */
Point *matrix_get_row(Matrix *matrix, int row_index) { /* TODO: what about diag matrix? */
    assert(matrix_get_rows_num(matrix) > row_index);
    int cols_num = matrix_get_cols_num(matrix);
    double *data = matrix_get_data(matrix) + row_index*cols_num;
    Point *point = create_point(data, cols_num, 0);
    return point;
}

/* TODO: FIX DATA LEAK */
Point *matrix_get_column(Matrix *matrix, int column_index) {
    int cols_num = matrix_get_cols_num(matrix);
    assert(cols_num > column_index);
    int rows_num = matrix_get_rows_num(matrix);
    double *data = matrix_get_data(matrix) + column_index;
    Point *point = create_point(data, rows_num, cols_num);
    return point;
}

void matrix_get_row_to_point(Matrix *matrix, Point *point, int row_index) { /* TODO: what about diag matrix? */
    assert(matrix_get_rows_num(matrix) > row_index);
    int cols_num = matrix_get_cols_num(matrix);
    point->dim = cols_num;
    point->offset = 0;
    point->data = matrix_get_data(matrix) + row_index*cols_num;
}

void matrix_get_column_to_point(Matrix *matrix, Point *point, int column_index) {
    int cols_num = matrix_get_cols_num(matrix);
    assert(cols_num > column_index);
    point->dim = matrix_get_rows_num(matrix);
    point->data = matrix_get_data(matrix) + column_index;
    point->offset = cols_num;
}

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

Matrix *multiply_matrices(Matrix *m1, Matrix *m2) {
    assert(matrix_get_cols_num(m1) == matrix_get_rows_num(m2));

    if (_is_matrix_diag(m1) && _is_matrix_diag(m2))
        return _multiply_matrices_diag_with_diag(m1, m2);

    else if (!_is_matrix_diag(m1) && !_is_matrix_diag(m2))
        return _multiply_matrices_nondiag_with_nondiag(m1, m2);
    
    else
        return _multiply_matrices_diag_with_nondiag(m1, m2);
}


void reset_matrix_entries_to_zero(Matrix *matrix){
    int i,j;
    int rows = matrix_get_rows_num(matrix);
    int cols = matrix_get_cols_num(matrix);
    for (i=0; i<rows; i++)
        for (j=0; j<cols; j++)
            matrix_set_entry(matrix, i,j,0);
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

void _clean_matrix(Matrix *matrix) {
    free(matrix_get_data(matrix));
    free(matrix);
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





/* debugging function */
void print_matrix(Matrix *matrix) {
    int i, j;
    double val;
    for (i=0; i<matrix->rows; i++) {
        for (j=0; j<matrix->cols; j++) {
            val = matrix_get_entry(matrix, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
}

void print_matrix_rows(Matrix *matrix) {
    int i;
    for (i=0; i<matrix_get_rows_num(matrix); i++) {
        Point *p = matrix_get_row(matrix, i);
        print_point(p);
        printf("\n");
    }
}

void print_matrix_cols(Matrix *matrix) {
    int i;
    for (i=0; i<matrix_get_cols_num(matrix); i++) {
        Point *p = matrix_get_column(matrix, i);
        print_point(p);
        printf("\n");
    }
}

double RandomReal(double low, double high)
{
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

void space() {
    printf("\n");
    printf("\n");
}


int main() {
    srand((int) time(NULL)); /* important for random */
    Matrix *m1 = generate_matrix(3, 5, false);
    Matrix *m2 = generate_matrix(5, 2, false);
    print_matrix(m1);
    space();

    Point *p = (Point *)malloc(sizeof(Point));
    matrix_get_column_to_point(m1, p, 2);
    print_point(p);
    space();
    matrix_get_row_to_point(m1, p, 1);
    print_point(p);
    space();
    free(p);

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