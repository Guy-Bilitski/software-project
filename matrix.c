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

typedef struct Centroids { /*IMPORTANT: column oriented, unlike matrix*/
    int k;
    int dim;
    int *vectors_in_cluster;
    double *data;
} Centroids;


typedef struct Point {
    int dim;
    double *data;
} Point;

/* TODO: Check callocs/mallocs on failure */

/* Matrix API */
Matrix *create_matrix(int rows, int cols, int is_diag);
double matrix_get(Matrix *matrix, int row, int col);
Point *matrix_get_point(Matrix *matrix, int row_num);
int matrix_get_rows_num(Matrix *matrix);
int matrix_get_cols_num(Matrix *matrix);
void matrix_set(Matrix *matrix, int row, int col, double value);
void matrix_set_point(Matrix *matrix, int row, Point *point);

/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);

/* Point API */
Point *create_empty_point(int dim);
Point *create_point(int dim, double *data);
double point_get_index(Point *point, int index);
int point_get_dim(Point *point);
void point_set_index(Point *point, int index, double value);

/* Centroids API */
Centroids *init_centroids(int dim, int k);

/* debugging functions */
void print_point(Point *point);
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
            return (matrix->data)[col];
        }
        
    } else { /* matrix is NOT diagonal */
        return matrix->data[_get_matrix_index(matrix, row, col)];
    }
}

Point *matrix_get_point(Matrix *matrix, int row_index) {
    int i, row_num = matrix_get_rows_num(matrix);
    double val;
    Point *point = create_empty_point(matrix_get_rows_num(matrix));

    for (i=0; i<row_num; i++) {
        val = matrix_get(matrix, row_index, i);
        point_set_index(point, i, val);
    }

    return point;
}

int matrix_get_rows_num(Matrix *matrix) {
    return matrix->rows;
}

int matrix_get_cols_num(Matrix *matrix) {
    return matrix->cols;
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

void matrix_set_point(Matrix *matrix, int row_index, Point *point) {
    /* sets a whole point in the matrix */
    int i;
    if (matrix->is_diag == false) {
        for (i=0; i<matrix->rows; i++) {
            matrix_set(matrix, row_index, i, point_get_index(point, i));
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
Point *create_empty_point(int dim) {
    Point *point = (Point *)malloc(sizeof(Point));
    point->data = (double *)calloc(sizeof(double), dim); 
    point->dim = dim;
    return point;
}

Point *create_point(int dim, double *data) {
    Point *new_point = create_empty_point(dim);
    new_point->data = data;
    return new_point;
}

double point_get_index(Point *point, int index) {
    assert(index >= 0 && index < point->dim);
    return point->data[index];
}

int point_get_dim(Point *point) {
    return point->dim;
}

void point_set_index(Point *point, int index, double value) {
    /* sets value in a specific index */
    assert(index >= 0 && index < point->dim);
    point->data[index] = value;
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
            val = matrix_get(matrix, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
}

void print_point(Point *point) {
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
    Point *p = create_point(4, new_row);
    matrix_set(m, 1, 2, 10);
    matrix_set(m, 0, 0, 20);
    matrix_set(m, 3, 3, 30);
    matrix_set(m, 3, 3, 5);
    matrix_set_point(m, 2, p);
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
    Point *point = create_empty_point(3);
    point_set_index(point, 0, 100);
    point_set_index(point, 1, 1);
    point_set_index(point, 2, 8);
    print_point(point);
    free(point);


    return 1;
}
/*
gcc matrix.c
a.out
*/
