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
double get(Matrix *matrix, int row, int col);
void set(Matrix *matrix, int row, int col, double value);
Matrix *create_matrix(int rows, int cols, int is_diag);


/* Matrix inner functions */
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);
void print_matrix(Matrix *matrix);





/* Matrix inner functions */

int _get_matrix_index(Matrix *matrix, int row, int col) {
    return row*(matrix->rows) + col;
}

void print_matrix(Matrix *matrix) {
    int i, j;
    double val;
    for (i=0; i<matrix->rows; i++) {
        for (j=0; j<matrix->cols; j++) {
            val = get(matrix, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
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



/* Matrix API */

double get(Matrix *matrix, int row, int col){
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


void set(Matrix *matrix, int row, int col, double value) {
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


Matrix *create_matrix(int rows, int cols, int is_diag) {
    if (is_diag)
        assert(rows==cols);
    assert(rows>0 && cols>0);
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->is_diag = is_diag;
    matrix->data = (double *)calloc(sizeof(double), rows*cols);
    return matrix;
}



int main() {
    Matrix *m = create_matrix(4,4,false);
    set(m, 1, 2, 10);
    set(m, 0, 0, 20);
    set(m, 3, 3, 30);
    set(m, 3, 3, 5);
    print_matrix(m);
    free(m);


    printf("\n");
    Matrix *d = create_matrix(4,4,true);
    set(d, 0, 0, 1);
    set(d, 1, 1, 2);
    set(d, 2, 2, 3);
    set(d, 3, 3, 4);
    set(d, 3, 3, 5);
    printf("%d%d%d\n",d->cols,d->rows,d->is_diag);
    set(d, 0, 1, 5);
    print_matrix(m);
    printf("%d%d%d\n",d->cols,d->rows,d->is_diag);

    return 1;
}
/*
gcc matrix.c
a.out
*/
