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

typedef struct Centroids { /*IMPORTANT: column oriented, unlike matrix*/
    int k;
    int dim;
    int *vectors_in_cluster;
    double *data;
} Centroids;


/* TODO: Check callocs/mallocs on failure */

/* Matrix API */
Matrix *create_matrix(int rows, int cols, int is_diag);  /* creates a matrix with dimensions rows X columns all zeros */

int matrix_get_rows_num(Matrix *matrix);  /* returns the number of rows in the matrix */
int matrix_get_cols_num(Matrix *matrix);  /* returns the number of columns in the matrix */
double *matrix_get_data(Matrix *matrix);  /* returns the matrix data */

double matrix_get_entry(Matrix *matrix, int row, int col);  /* returns the (row, col) entry */
Point *matrix_get_row(Matrix *matrix, int row_index);  /* returns the <row_index> point (row) of the matrix */
Point *matrix_get_column(Matrix *matrix, int column_index); /* returns the <column_index> column of the matrix */

void matrix_set_entry(Matrix *matrix, int row, int col, double value);  /* sets <value> in (row, col) entry */
void matrix_set_point(Matrix *matrix, int row_index, Point *point);  /* sets a point (row) in the matrix in <row_index> */
/*Matrix *multiply_matrices(Matrix *m1, Matrix *m2);  /* multiply m1 X m2 and returns the new matrix */

/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);
void _clean_matrix(Matrix *matrix);


/* Centroids API */
Centroids *init_centroids(int dim, int k);

/* debugging functions */
void print_matrix(Matrix *matrix);
void space();


/* Matrix API */
Matrix *create_matrix(int rows, int cols, int is_diag) {
    if (is_diag) {
        assert(rows==cols);
    }
    assert(rows>0 && cols>0);
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->is_diag = is_diag;
    matrix->data = (double *)calloc(sizeof(double), rows*cols);
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

double matrix_get_entry(Matrix *matrix, int row, int col){
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));

    if (matrix->is_diag){ /* if matrix is diagonal */
        if (row != col){
            return 0;
        } else {
            return (matrix->data)[col];
        }
        
    } else { /* matrix is NOT diagonal */
        return matrix->data[_get_matrix_index(matrix, row, col)];
    }
}


Point *matrix_get_row(Matrix *matrix, int row_index) {
    assert(matrix_get_rows_num(matrix) > row_index);
    int cols_num = matrix_get_cols_num(matrix);
    double *data = matrix_get_data(matrix) + row_index*cols_num;
    Point *point = create_point(data, cols_num, 0);
    return point;
}

Point *matrix_get_column(Matrix *matrix, int column_index) {
    int cols_num = matrix_get_cols_num(matrix);
    assert(cols_num > column_index);
    int rows_num = matrix_get_rows_num(matrix);
    double *data = matrix_get_data(matrix) + column_index;
    Point *point = create_point(data, rows_num, cols_num);
    return point;
}

void matrix_set_entry(Matrix *matrix, int row, int col, double value) {
    /* sets an entry in the matrix */
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));

    if (matrix->is_diag) {
        if (row != col && value != 0) {
            _diag_to_square_matrix(matrix);
            (matrix->data)[_get_matrix_index(matrix, row, col)] = value;
        } else if (row != col && value == 0){ /* nothing to change */
            return;
        }
        else {
            (matrix->data)[row] = value;
        }
    } 
    else { /* matrix is NOT diagonal */
        (matrix->data)[_get_matrix_index(matrix, row, col)] = value;
    }
}

void matrix_set_point(Matrix *matrix, int row_index, Point *point) {
    /* sets a whole point in the matrix */
    int i;
    if (matrix->is_diag == false) {
        for (i=0; i<matrix->rows; i++) {
            matrix_set_entry(matrix, row_index, i, point_get_entry(point, i));
        }
    }
}

Matrix *multiply_matrices(Matrix *m1, Matrix *m2) {
    assert(matrix_get_cols_num(m1) == matrix_get_rows_num(m2));
    int i, j, rows_num = matrix_get_rows_num(m1), columns_num = matrix_get_cols_num(m2);
    Point *row, *column;
    Matrix *new_matrix = create_matrix(rows_num, columns_num, false);
    for (i=0; i<rows_num; i++) {
        row = matrix_get_row(m1, i);
        for (j=0; j<columns_num; j++) {
            column = matrix_get_column(m2, j);
            matrix_set_entry(new_matrix, i, j, multiply_points(row, column));
        }
    }
    free(row);
    free(column);
    return new_matrix;
}


/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col) {
    return row*(matrix->cols) + col;
}

void _diag_to_square_matrix(Matrix *matrix) { 
    /* Changing matrix->data from n sized array (diagonal matrix) to n*n sized array (square matrix) */
    int i,n;
    double *old_data = matrix->data;
    n = matrix->rows;
    matrix->data = (double *)calloc(sizeof(double), n*n);
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


/* Centroids API */
double get_centroid_entry(Centroids *centroids, int centroid_idx, int entry){
    assert(0 <= centroid_idx && centroid_idx < (centroids->k)); /*TODO: Delete*/
    assert(0 <= entry && entry < (centroids->dim)); /*TODO: Delete*/
    assert(0 < (centroids->vectors_in_cluster)[centroid_idx]); /* to keep assertion */
    int num_of_vectors_in_cluster = (centroids->vectors_in_cluster)[centroid_idx];
    return (centroids->data)[centroid_idx*(centroids->dim)+entry] / num_of_vectors_in_cluster;
}

Centroids *init_centroids(int dim, int k) {
    assert(dim>0 && k>0);
    int i;
    Centroids *centroids = (Centroids *)malloc(sizeof(Centroids));
    centroids->dim = dim;
    centroids->k = k;
    centroids->vectors_in_cluster = (int *)calloc(sizeof(int), k);
    centroids->data = (double *)calloc(sizeof(double), dim*k);
    for (i=0; i<k; i++)
        (centroids->vectors_in_cluster)[i] = 0;
    for (i=0; i<dim*k; i++)
        (centroids->data)[i] = 0.;
    return centroids;
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
    Matrix *m = create_matrix(rows, cols, false);
    if (!is_diag) {
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
    print_matrix(m2);
    space();
    
    Matrix *product = multiply_matrices(m1, m2);

    print_matrix_rows(product);
    print_matrix_cols(product);
    
    _clean_matrix(m1);
    _clean_matrix(m2);
    _clean_matrix(product);    

    return 1;
}
/*
gcc matrix.c
a.out
*/
