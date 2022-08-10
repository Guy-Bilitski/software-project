#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
double matrix_get_entry(Matrix *matrix, int row, int col);  /* returns the (row, col) entry */
Point *matrix_get_point(Matrix *matrix, int row_index);  /* returns the <row_index> point (row) of the matrix */
Point *matrix_get_column(Matrix *matrix, int column_index); /* returns the <column_index> column of the matrix */
int matrix_get_rows_num(Matrix *matrix);  /* returns the number of rows in the matrix */
int matrix_get_cols_num(Matrix *matrix);  /* returns the number of columns in the matrix */
void matrix_set_entry(Matrix *matrix, int row, int col, double value);  /* sets <value> in (row, col) entry */
void matrix_set_point(Matrix *matrix, int row_index, Point *point);  /* sets a point (row) in the matrix in <row_index> */
Matrix *multiply_matrices(Matrix *m1, Matrix *m2);  /* multiply m1 X m2 and returns the new matrix */

/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);


/* Centroids API */
Centroids *init_centroids(int dim, int k);

/* debugging functions */
void print_matrix(Matrix *matrix);


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

Point *matrix_get_point(Matrix *matrix, int row_index) {
    int i, cols_num = matrix_get_cols_num(matrix);
    double value;
    Point *point = create_empty_point(cols_num);

    for (i=0; i<cols_num; i++) {
        value = matrix_get_entry(matrix, row_index, i);
        point_set_index(point, i, value);
    }

    return point;
}

Point *matrix_get_column(Matrix *matrix, int column_index) {
    int i, column_dim = matrix_get_rows_num(matrix);
    double value;
    Point *new_point = create_empty_point(column_dim);
    for (i=0; i<column_dim; i++) {
        value = matrix_get_entry(matrix, i, column_index);
        point_set_index(new_point, i, value);
    }
    return new_point;
}

int matrix_get_rows_num(Matrix *matrix) {
    return matrix->rows;
}

int matrix_get_cols_num(Matrix *matrix) {
    return matrix->cols;
}

void matrix_set_entry(Matrix *matrix, int row, int col, double value) {
    /* sets an entry in the matrix */
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));

    if (matrix->is_diag) {
        printf("diag");
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
            matrix_set_entry(matrix, row_index, i, point_get_index(point, i));
        }
    }
}

Matrix *multiply_matrices(Matrix *m1, Matrix *m2) {
    assert(matrix_get_cols_num(m1) == matrix_get_rows_num(m2));
    int i, j, rows_num = matrix_get_rows_num(m1), columns_num = matrix_get_cols_num(m2);
    Point *row, *column;
    Matrix *new_matrix = create_matrix(rows_num, columns_num, false);
    for (i=0; i<rows_num; i++) {
        row = matrix_get_point(m1, i);
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


int main() {
    /* non diagonal matrix */
    /*
    Matrix *m = create_matrix(4,4,false);
    printf("\n");
    printf("%d", matrix_get_rows_num(m));
    printf("\n");
    double *new_row = malloc(sizeof(double)*4);
    int i;
    for (i=0; i<4; i++) {
        new_row[i] = i + 100;
    }
    Point *p = create_point(4, new_row);
    matrix_set_entry(m, 1, 2, 10);
    matrix_set_entry(m, 0, 0, 20);
    matrix_set_entry(m, 3, 3, 30);
    matrix_set_entry(m, 3, 3, 5);
    matrix_set_point(m, 2, p);
    print_matrix(m);

    printf("\n");
    print_point(matrix_get_point(m, 0));
    printf("\n");
    


    printf("\n");


    Matrix *m2 = create_matrix(4,4,false);
    matrix_set_entry(m2, 1, 2, 4);
    matrix_set_entry(m2, 0, 0, 3);
    matrix_set_entry(m2, 0, 2, 5.5);
    matrix_set_entry(m2, 3, 3, 4);
    matrix_set_entry(m2, 2, 2, 3);
    matrix_set_entry(m2, 3, 1, 2);
    matrix_set_entry(m2, 2, 1, 1);
    matrix_set_entry(m2, 1, 1, 1);
    matrix_set_entry(m2, 3, 2, 1.7);

    print_matrix(m2);
    
    printf("\n");

    Matrix *product = multiply_matrices(m, m2);
    print_matrix(product);
    free(product);
    free(m2);
    free(m);
    free(new_row);

*/
    
    printf("\n");
    printf("\n");
    Matrix *m11 = create_matrix(3,4,false);
    Matrix *m22 = create_matrix(4,5,false);

    matrix_set_entry(m11, 1, 2, 4);
    matrix_set_entry(m11, 0, 0, 3);
    matrix_set_entry(m11, 0, 2, 5.5);
    matrix_set_entry(m11, 2, 3, 4);
    matrix_set_entry(m11, 2, 2, 3);
    matrix_set_entry(m11, 2, 1, 2);
    matrix_set_entry(m11, 0, 1, 1);
    matrix_set_entry(m11, 1, 1, 1);
    matrix_set_entry(m11, 0, 3, 1.7);

    print_matrix(m11);


    matrix_set_entry(m22, 1, 2, 4);
    matrix_set_entry(m22, 0, 0, 3);
    matrix_set_entry(m22, 0, 2, 5.5);
    matrix_set_entry(m22, 3, 3, 4);
    matrix_set_entry(m22, 2, 4, 3);
    matrix_set_entry(m22, 3, 1, 2);
    matrix_set_entry(m22, 2, 1, 1);
    matrix_set_entry(m22, 1, 1, 1);
    matrix_set_entry(m22, 3, 4, 1.7);

    printf("\n");
    print_matrix(m22);
    printf("\n");

    Matrix *product1 = multiply_matrices(m11, m22);
    print_matrix(product1);

    free(m22);
    free(m11);    
    free(product1);


    /* diagonal matrix */
    /*
    Matrix *d = create_matrix(4,4,true);
    matrix_set_entry(d, 0, 0, 1);
    matrix_set_entry(d, 1, 1, 2);
    matrix_set_entry(d, 2, 2, 3);
    matrix_set_entry(d, 3, 3, 4);
    matrix_set_entry(d, 3, 3, 5);
    printf("%d%d%d\n",d->cols,d->rows,d->is_diag);
    matrix_set_entry(d, 0, 1, 5);
    print_matrix(m);
    printf("%d%d%d\n",d->cols,d->rows,d->is_diag);
    free(d);
    */

    /* point */
    /*
    Point *point = create_empty_point(3);
    point_set_index(point, 0, 100);
    point_set_index(point, 1, 1);
    point_set_index(point, 2, 8);
    print_point(point);
    free(point);
    */
    


    return 1;
}
/*
gcc matrix.c
a.out
*/
