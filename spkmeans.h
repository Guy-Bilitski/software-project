#ifndef spkmeansh
#define spkmeansh
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

/* ---------------- STRUCTS ---------------- */

#ifndef POINT_IS_DEFINED
#define POINT_IS_DEFINED
typedef struct Point {
    double *data;
    int dim;
    int offset; /* for column representation */
} Point;
#endif

#ifndef MATRIX_IS_DEFINED
#define MATRIX_IS_DEFINED
typedef struct Matrix {
    int rows;
    int cols;
    int is_not_diag;
    double *data;
} Matrix;
#endif

#ifndef EIGENVECTOR_IS_DEFINED
#define EIGENVECTOR_IS_DEFINED
typedef struct Eigenvector {
    Point *point;
    double eigen_value;
} Eigenvector;
#endif

#ifndef S_AND_C_IS_DEFINED
#define S_AND_C_IS_DEFINED
typedef struct S_and_C
{
    double s;
    double c;
} S_and_C;
#endif

#ifndef MAX_ELEMENT_IS_DEFINED
#define MAX_ELEMENT_IS_DEFINED
typedef struct MaxElement
{
    int i;
    int j;
    double value;
} MaxElement;
#endif


#ifndef JACOBI_OUTPUT_IS_DEFINED
#define JACOBI_OUTPUT_IS_DEFINED
typedef struct JacobiOutput
{
    Matrix *A;
    Matrix *V;
} JacobiOutput;
#endif

/* -------------------- POINT PROTOTYPES -------------------- */
/* Struct point is an interface of data stored as a list of doubles, like a mathematical vector */
/* the Point gives a clean, easy to work with way of treating a matrix as rows and columns */
/* using offset field you can create a point the represents a column of a matrix */

/* Point API */
Point *create_empty_point(); /* creates an empty point in the memory */
Point *create_point(double *data, int dim, int offset);  /* creates a point from list */

/* Getters */
double point_get_entry(Point *point, int entry);  /* returns the value in index <index> of point */
int point_get_dim(Point *point);  /* returns point's dim */
int point_get_offset(Point *point);  /* returns point's offset */
double *point_get_data(Point *point);  /* returns point's data */

/* Setters */
void point_set_values(Point *point, int dim, int offset, double *data);  /* sets point new values */

/* Utils */
double inner_product(Point *row_point, Point *column_point);  /* returns row X column scalar */
double euclidean_distance(Point *p1, Point *p2);  /* returns the euclidian distance between two points (vectors) */
double euclidean_norm(Point *p);  /* returns the eucalidian norm of a point (vector) */
void divide_point_by_value(Point *p, double value);  /* divides each point entry by a given value */

/* Point inner functions */
int _convert_point_index(Point *point, int index);  /* converts given index to the "real" index considering the offset */

/* debugging functions */
void print_point(Point *point);

/* -------------------- MATRIX PROTOTYPES -------------------- */
/* Struct matrix is an interface of a list of doubles representing a matrix */
/* a matrix can be diagnal or not, consists rows and column dimensions, and data */
/* upon the matrix interface there are many known matrix functionalities */


/* Matrix API */
Matrix *create_matrix(int rows, int cols);  /* creates a matrix in the memory containing 0 in all entries */
Matrix *create_identity_matrix(int n);  /* crate the I matrix with 1-diagonal 0-non diagonal */

/* Getters */
int matrix_get_rows_num(Matrix *matrix);  /* returns matrix's rows dim */
int matrix_get_cols_num(Matrix *matrix);  /* returns matrix's columns dim */
double *matrix_get_data(Matrix *matrix);  /* returns matrix's data */
double matrix_get_entry(Matrix *matrix, int row, int col);  /* returns the (row, col) entry value */
void matrix_get_row_to_point(Matrix *matrix, Point *point, int row_index); /* inserts a row (represented as a point) into a given point */
void matrix_get_column_to_point(Matrix *matrix, Point *point, int column_index);  /* inserts a column (represented as a point) into a given point */

/* Setters */
void matrix_set_entry(Matrix *matrix, int row, int col, double value);  /* sets <value> in (row, col) entry */

/* Utils */
void matrix_get_non_diagonal_max_absolute_value(Matrix *matrix, MaxElement *max_element);  /* returns the max element of the matrix that is not on the diagonal */
double matrix_get_row_sum(Matrix *matrix, int row_index);  /* returns <row_index> row sum of values */
void matrix_add_point_to_row(Matrix *matrix, int row_index, Point *point); /* adds a row to the matrix value by value */
void reset_matrix_entries_to_zero(Matrix *matrix);  /* resets all metrix entries to zero */
Matrix *multiply_matrices(Matrix *m1, Matrix *m2);  /* multiply m1 X m2 and returns the new matrix */
Matrix *sub_matrices(Matrix *A, Matrix *B); /* sub A - B */

/* Cleanup */
void free_matrix(Matrix *matrix); /* free matrix object and data */

/* Matrix inner functions */
int _is_matrix_diag(Matrix *matrix);
int _get_matrix_index(Matrix *matrix, int row, int col);
Matrix *_multiply_matrices_diag_with_diag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_diag_with_nondiag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_nondiag_with_nondiag(Matrix *m1, Matrix *m2);

/* debugging functions */
void print_matrix(Matrix *matrix);
void print_matrix2(Matrix *matrix);
void print_matrix_diag(Matrix *matrix);
void print_matrix_rows(Matrix *matrix);
void print_matrix_cols(Matrix *matrix);
double RandomReal(double low, double high);
Matrix *generate_matrix(int rows, int cols, int is_diag);
Matrix *generate_symmetric_matrix(int n);
void space();

/* -------------------- EIGENVECTOR PROTOTYPES -------------------- */

/* EIGENVECTOR API */
Eigenvector *create_empty_eigen_vector();
Eigenvector *create_eigen_vector(Point *point, double eigen_value);
Eigenvector *create_eigen_vectors_array(int eigenvectors_num);

/* Getters */
Point *eigen_vector_get_point(Eigenvector *eigen_vector);
double eigen_vector_get_eigen_value(Eigenvector *eigen_vector);

/* Utils */
int compare_eigenvectors(const void *p1, const void *p2);
void sort_eigenvectors_array(Eigenvector *array, size_t n);

/* Cleanup */
void free_eigen_vector(Eigenvector *eigen_vector);

/* debugging functions */
void print_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n);

/* -------------------- S AND C PROTOTYPES -------------------- */

/* S AND C API */
S_and_C *create_empty_S_and_C();

/* Getters */
double s_and_c_get_s(S_and_C *s_and_c);
double s_and_c_get_c(S_and_C *s_and_c);

/* Setters */
void S_and_C_set_values(S_and_C *s_and_c, double s, double c);

/* debugging functions */
void print_s_and_c(S_and_C *s_and_c);

/* -------------------- MAX ELEMENT PROTOTYPES -------------------- */

/* MaxElemnt API */
MaxElement *create_empty_max_element();
MaxElement *create_max_element(double value, int i, int j);

/* Getters */
double max_element_get_value(MaxElement *max_element);
int max_element_get_index1(MaxElement *max_element);
int max_element_get_index2(MaxElement *max_element);

/* Setters */
void max_element_set_new_values(MaxElement *max_element, double value, int i, int j);
void max_element_set_value(MaxElement *max_element, double value);
void max_element_set_index1(MaxElement *max_element, int i);
void max_element_set_index2(MaxElement *max_element, int j);

/* debugging */
void print_max_element(MaxElement *max_element);

/* -------------------- JACOBI OUTPUT PROTOTYPES -------------------- */

/* Jacobi Output API */
JacobiOutput *create_empty_jacobi_output();

/* Getters */
Matrix *jacobi_output_get_A(JacobiOutput *jacobi_output);
Matrix *jacobi_output_get_V(JacobiOutput *jacobi_output);

/* Setters */
void set_jacobi_output_values(JacobiOutput *jacobi_output, Matrix *A, Matrix *V);

/* Cleanup */
void free_jacobi_output(JacobiOutput *jacobi_output);

/* debugging */
void print_jacobi_output(JacobiOutput *J);

/* -------------------- KMEANS-IO PROTOTYPES -------------------- */

/* Kmeans_io API */
int get_dimension(const char *input_file);
int get_n(const char *input_file);
Matrix *input_file_to_matrix(const char *input_file);

/* -------------------- KMEANS PROTOTYPES -------------------- */

/* Kmeans API */
Matrix * kmeans(Matrix *data_points, Matrix *centroids);
double max_distance_between_centroids(Matrix *old_centroids, Matrix *new_centroids);
void kmeans_iteration(Matrix *data_points , Matrix *centroids, Matrix *new_centroids);
int find_closest_centroid(Point *vector, Matrix *centroids);

/* -------------------- SPKMEANS PROTOTYPES -------------------- */

/* spkmeans functions */
void achieve_goal(Matrix *data_points, char *goal);
double gaussian_RBF(Point *x1, Point *x2);  /*computes w_i in the weighted adjacency matrix*/
Matrix *create_weighted_matrix(Matrix *X);  /* creates the weighted matrix */
Matrix *create_diagonal_degree_matrix(Matrix *matrix); /* retruns the I matrix */
void neg_root_to_diag_matrix(Matrix *matrix); /* performs pow of -0.5 for all the diagonal entries */
Matrix *normalized_graph_laplacian(Matrix *D_minus_05, Matrix *W);
void get_s_and_c_for_rotation_matrix(Matrix* A, MaxElement *max_element, S_and_C *s_and_c);
void build_rotation_matrix(S_and_C *s_and_c, MaxElement *max_element, Matrix *identity_matrix); /* inserts the rotation matrix to a given identity matrix */
void normalize_matrix_rows(Matrix *matrix);
double matrix_off(Matrix *matrix); /* returns the value of "off" function on a given matrix */
Matrix *transform_matrix(Matrix *matrix, S_and_C *s_and_c, MaxElement *max_element);  /* permorms matrix transformation */
int get_k_from_sorted_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n);
void get_eigen_vectors_from_jacobi_output(JacobiOutput *jacobi_output, Eigenvector *eigen_vectors_array);
Matrix *getU(JacobiOutput *jacobi_output, int k);
double get_value_for_transformed_matrix(Matrix *old_matrix, double s, double c, int i, int j, int row_index, int col_index); /* returns the expected value of the transformed matrix at (row_index, col_index) based on the rules described at 6. Relations betweeb A and A'*/
int matrix_converge(double A_off, Matrix *A);

/* SPKMEANS API */
Matrix *wam(Matrix* data_points);
Matrix *ddg(Matrix* data_points);
Matrix *lnorm(Matrix* data_points);
JacobiOutput *jacobi(Matrix *A, JacobiOutput *jacobi_output);

#endif