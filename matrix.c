#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define true 1
#define false 0


typedef struct Matrix {
    int rows;
    int cols;
    int is_diag;
    double *data;
} Matrix;


typedef struct Point {
    int dim;
    double *data;
} Point;

/* TODO: Check callocs/mallocs on failure */

/* Matrix API */
Matrix *create_matrix(int rows, int cols, int is_diag);
double matrix_get(Matrix *matrix, int row, int col);
void matrix_set(Matrix *matrix, int row, int col, double value);
void matrix_set_row(Matrix *matrix, int row, double *new_row);

/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);

/* Point API */
Point *create_point(int dim);
void point_set(Point *point, int index, double value);
double point_get(Point *point, int index);

/* debugging functions */
void *print_point(Point *point);
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

double matrix_get(Matrix *matrix, int row, int col){
    assert(!(matrix->rows <= row || matrix->cols <= col || row < 0 || col < 0));

    if (matrix->is_diag){ /* if matrix is diagonal */
        if (row != col){
            return 0;
        } else {
            return matrix->data[col];
        }
        
    } else { /* matrix is NOT diagonal */
        return matrix->data[_get_matrix_index(matrix, row, col)];
    }
}

void matrix_set(Matrix *matrix, int row, int col, double value) {
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

void matrix_set_row(Matrix *matrix, int row, double *new_row) {
    /* sets a whole row in the matrix */
    int i;
    if (matrix->is_diag == false) {
        for (i=0; i<matrix->rows; i++) {
            matrix_set(matrix, row, i, new_row[i]);
        }
    }
}


/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col) {
    return row*(matrix->rows) + col;
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


/* Point API */
double point_get(Point *point, int index) {
    assert(index >= 0 && index < point->dim);
    return point->data[index];
}

void point_set(Point *point, int index, double value) {
    assert(index >= 0 && index < point->dim);
    point->data[index] = value;
}

Point *create_point(int dim) {
    Point *point = (Point *)malloc(sizeof(Point));
    point->data = (double *)calloc(sizeof(double), dim); 
    point->dim = dim;
    return point;
}


/* debugging function */
void print_matrix(Matrix *matrix) {
    int i, j;
    double val;
    for (i=0; i<matrix->rows; i++) {
        for (j=0; j<matrix->cols; j++) {
            val = matrix_get(matrix, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
}

void *print_point(Point *point) {
    int i;
    for (i=0; i<(point->dim); i++) {
        printf("%f ", (point->data)[i]);
    }
}


int main() {
    /* non diagonal matrix */
    Matrix *m = create_matrix(4,4,false);
    double *new_row = malloc(sizeof(double)*4);
    int i;
    for (i=0; i<4; i++) {
        new_row[i] = i + 100;
    }
    matrix_set(m, 1, 2, 10);
    matrix_set(m, 0, 0, 20);
    matrix_set(m, 3, 3, 30);
    matrix_set(m, 3, 3, 5);
    matrix_set_row(m, 3, new_row);
    print_matrix(m);
    free(m);
    free(new_row);


    printf("\n");

    /* diagonal matrix */
    Matrix *d = create_matrix(4,4,true);
    matrix_set(d, 0, 0, 1);
    matrix_set(d, 1, 1, 2);
    matrix_set(d, 2, 2, 3);
    matrix_set(d, 3, 3, 4);
    matrix_set(d, 3, 3, 5);
    printf("%d%d%d\n",d->cols,d->rows,d->is_diag);
    matrix_set(d, 0, 1, 5);
    print_matrix(m);
    printf("%d%d%d\n",d->cols,d->rows,d->is_diag);
    free(d);

    /* point */
    Point *point = create_point(3);
    point_set(point, 0, 100);
    point_set(point, 1, 1);
    point_set(point, 2, 8);
    print_point(point);
    free(point);


    return 1;
}
/*
gcc matrix.c
a.out
*/
