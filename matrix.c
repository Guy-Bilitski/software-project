#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Matrix {
    int rows;
    int col;
    double *data;
} Matrix;

typedef struct Diagonal_Matrix {
    int rows;
    int col;
    double *data; /* only for diagonal */
} Diagonal_Matrix;


/* General matrix internal function */
int _get_matrix_index(Matrix *matrix, int row, int col);
void print_matrix(Matrix *matrix);
int get_diag_size(Matrix *matrix);

/* regular matrix functions */
double get_matrix_entry(Matrix *matrix, int row, int col);
void set_matrix_entry(Matrix *matrix, int row, int col, double value);

/* diagonal matrix functions */
double get_diagonal_matrix_entry(Matrix *matrix, int row, int col);
void set_diagonal_matrix_entry(Matrix *matrix, int row, int col, double value);


/* General matrix internal function */

int _get_matrix_index(Matrix *matrix, int row, int col) {
    return row*(matrix->rows) + col;
}

void print_matrix(Matrix *matrix) {
    int i, j;
    double val;
    for (i=0; i<matrix->rows; i++) {
        for (j=0; j<matrix->col; j++) {
            val = get_matrix_entry(matrix, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
}

int get_diag_size(Matrix *matrix) {
    double rows_double = pow(matrix->rows, 2);
    double columns_double = pow(matrix->col, 2);
    return (int)sqrt(rows_double + columns_double);
}


/* regular matrix functions */

double get_matrix_entry(Matrix *matrix, int row, int col) {
    int matrix_index = _get_matrix_index(matrix, row, col);
    return (matrix->data)[matrix_index];
}

void set_matrix_entry(Matrix *matrix, int row, int col, double value) {
    int matrix_index = _get_matrix_index(matrix, row, col);
    (matrix->data)[matrix_index] = value;
}

Matrix *create_matrix(rows, col) {
    Matrix *matrix = malloc(sizeof(Matrix));
    matrix->rows = rows;
    matrix->col = col;
    matrix->data = calloc(sizeof(double), rows*col);
    return matrix;
}


/* diagonal matrix functions */

double get_diagonal_matrix_entry(Matrix *matrix, int row, int col) {
    if (row != col) { /* non diagonal entry value is zero */
        return 0;
    }
    int matrix_index = _get_matrix_index(matrix, row, col);
    return (matrix->data)[matrix_index];
}

void set_diagonal_matrix_entry(Matrix *matrix, int row, int col, double value) {
    if (row == col) { /* it is wrong to change non diagonal entry */
        int matrix_index = _get_matrix_index(matrix, row, col);
        (matrix->data)[matrix_index] = value;
    }   
}



int main() {
    Matrix *m = create_matrix(4,4);
    set_matrix_entry(m, 1, 2, 2);
    print_matrix(m);
    return 1;
}