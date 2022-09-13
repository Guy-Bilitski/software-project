#ifndef spkmeansh
#define spkmeansh
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

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
Matrix *multiply_matrices(Matrix *m1, Matrix *m2);  /* multiply m1 X m2 and retumalrns the new matrix */
Matrix *sub_matrices(Matrix *A, Matrix *B); /* sub A - B */

/* Cleanup */
void free_matrix(Matrix *matrix); /* free matrix object and data */

/* Matrix inner functions */
int _is_matrix_diag(Matrix *matrix);
int _get_matrix_index(Matrix *matrix, int row, int col);
Matrix *_multiply_matrices_diag_with_diag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_diag_with_nondiag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_nondiag_with_nondiag(Matrix *m1, Matrix *m2);

void print_matrix(Matrix *matrix);
void print_matrix_diag(Matrix *matrix);

/* -------------------- EIGENVECTOR PROTOTYPES -------------------- */
/* Eigenvector struct represents an eigen vector (point) that has an eigen value */

/* EIGENVECTOR API */
Eigenvector *create_empty_eigen_vector();  /* creates an empty eigenvector in the memory  */
Eigenvector *create_eigen_vector(Point *point, double eigen_value);  /* creates an eigenvector in the memory with given values */
Eigenvector *create_eigen_vectors_array(int eigenvectors_num);  /* creates an eigenvector_num size array of eigenvectors */

/* Getters */
Point *eigen_vector_get_point(Eigenvector *eigen_vector);  /* returns the vector (point) of the eigenvector */
double eigen_vector_get_eigen_value(Eigenvector *eigen_vector);  /* returns the eigen value of the eigenvector */

/* Utils */
int compare_eigenvectors(const void *p1, const void *p2);  /* Eigenvectors comperator based on eigenvector's eigen values */
void sort_eigenvectors_array(Eigenvector *array, size_t n);  /* sorts eigenvectors array based on compare_eigenvectors comperator */

/* Cleanup */
void free_eigen_vector(Eigenvector *eigen_vector);  /* frees Eigen vector struct and its point property */


/* -------------------- S AND C PROTOTYPES -------------------- */
/* S_AND_C struct represents the values of s and c when building a rotation matrix */

/* S AND C API */
S_and_C *create_empty_S_and_C();  /* creates an empty s_and_c in the memory  */

/* Getters */
double s_and_c_get_s(S_and_C *s_and_c);  /* returns s_and_c value of s */
double s_and_c_get_c(S_and_C *s_and_c);  /* returns s_and_c value of c */

/* Setters */
void S_and_C_set_values(S_and_C *s_and_c, double s, double c);  /* sets s,c values to s_and_c struct */

/* -------------------- MAX ELEMENT PROTOTYPES -------------------- */
/* MaxElement struct represents a value with <i,j> coordinated (e.g of a matrix) */

/* MaxElemnt API */
MaxElement *create_empty_max_element();  /* creates an empty max_element in the memory  */
MaxElement *create_max_element(double value, int i, int j);  /* creates max element with given values */

/* Getters */
double max_element_get_value(MaxElement *max_element);  /* returns max_element value */
int max_element_get_index1(MaxElement *max_element);  /* returns max_element first index */
int max_element_get_index2(MaxElement *max_element);  /* returns max_element second index */

/* Setters */
void max_element_set_new_values(MaxElement *max_element, double value, int i, int j);  /* sets max_element values */
void max_element_set_value(MaxElement *max_element, double value);  /* sets max_element value */
void max_element_set_index1(MaxElement *max_element, int i);  /* sets max_element first index */
void max_element_set_index2(MaxElement *max_element, int j);  /* sets max_element second index */


/* -------------------- JACOBI OUTPUT PROTOTYPES -------------------- */
/* JacobiOutput struct represents the output of jacobi process containing */
/* eigen vectors matrix V and diagonal eigen values matrix A */

/* Jacobi Output API */
JacobiOutput *create_empty_jacobi_output();  /* creates an empty jacobi_output in the memory  */

/* Getters */
Matrix *jacobi_output_get_A(JacobiOutput *jacobi_output);  /* returns matrix A of Jacobi Output */
Matrix *jacobi_output_get_V(JacobiOutput *jacobi_output);  /* returns matrix V of Jacobi Output */

/* Setters */
void set_jacobi_output_values(JacobiOutput *jacobi_output, Matrix *A, Matrix *V);  /* sets A,V matrices to Jacoib Output */

/* Cleanup */
void free_jacobi_output(JacobiOutput *jacobi_output);  /* free Jacobi Output struct and A, V matrices */

/* Print */
void print_jacobi_output(JacobiOutput *J); /* prints the eigenvectors and eigenvalues as required */

/* -------------------- KMEANS-IO PROTOTYPES -------------------- */

/* Kmeans_io API */
int get_dimension(const char *input_file);  /* returns the dimension of points given in input_file */
int get_n(const char *input_file);  /* returns the number of data points given in input_file */
Matrix *input_file_to_matrix(const char *input_file);  /* returns data points from input_file as a matrix */

/* -------------------- KMEANS PROTOTYPES -------------------- */

/* Kmeans API */
Matrix * kmeans(Matrix *data_points, Matrix *centroids);  /* performs kmeans algorithm on data_points matrix using the given centroids */
double max_distance_between_centroids(Matrix *old_centroids, Matrix *new_centroids);  /* returns the maximum distance between given centroids calculated one by one */
void kmeans_iteration(Matrix *data_points , Matrix *centroids, Matrix *new_centroids);  /* performs kmeans algorithm single iteration (centroids improvement) */
int find_closest_centroid(Point *vector, Matrix *centroids);  /* returns the index of the closet centroid to a given vector */

/* -------------------- SPKMEANS PROTOTYPES -------------------- */

/* spkmeans functions */
void achieve_goal(Matrix *data_points, char *goal);
double gaussian_RBF(Point *x1, Point *x2);  /*computes w_i in the weighted adjacency matrix*/
Matrix *create_weighted_matrix(Matrix *X);  /* creates the weighted matrix */
Matrix *create_diagonal_degree_matrix(Matrix *matrix); /* retruns the I matrix */
void neg_root_to_diag_matrix(Matrix *matrix); /* performs pow of -0.5 for all the diagonal entries */
Matrix *normalized_graph_laplacian(Matrix *D_minus_05, Matrix *W); /* Construct Lnorm, inner function */
void get_s_and_c_for_rotation_matrix(Matrix* A, MaxElement *max_element, S_and_C *s_and_c); /* inner function, getting values for the rotation matrix construction */
void build_rotation_matrix(S_and_C *s_and_c, MaxElement *max_element, Matrix *identity_matrix); /* inserts the rotation matrix to a given identity matrix */
void normalize_matrix_rows(Matrix *matrix); /* normalize matrix rows for them to be unit vectors. assumes each row != 0 */
double matrix_off(Matrix *matrix); /* returns the value of "off" function on a given matrix */
Matrix *transform_matrix(Matrix *matrix, S_and_C *s_and_c, MaxElement *max_element);  /* permorms matrix transformation */
int get_k_from_sorted_eigen_vectors_array(Eigenvector *eigen_vectors_array, int n); /* The Eigengap Heuristic */
void get_eigen_vectors_from_jacobi_output(JacobiOutput *jacobi_output, Eigenvector *eigen_vectors_array); /* Turning the matrix into eigenvecotrs array, for sorting purposes */
Matrix *getU(JacobiOutput *jacobi_output, int k); /* getting U as in the main algorithm */
double get_value_for_transformed_matrix(Matrix *old_matrix, double s, double c, int i, int j, int row_index, int col_index); /* returns the expected value of the transformed matrix at (row_index, col_index) based on the rules described at 6. Relations betweeb A and A'*/
int matrix_converge(double A_off, Matrix *A); /* stop condition for jacobi */
void retrieve_identity_from_rotation_matrix(MaxElement *max_element, Matrix *P); /* making the existed rotation matrix identity matrix again */
void update_eigenvectors_matrix_V(Matrix *V, Matrix *P, MaxElement *max_element); /* efficient matrix multiplication for that specific use case (updating V <- VxP) */

/* SPKMEANS API */
/* API for the CPython wrapper. each function's wrapper calls a single function */
Matrix *wam(Matrix* data_points);
Matrix *ddg(Matrix* data_points);
Matrix *lnorm(Matrix* data_points);
JacobiOutput *jacobi(Matrix *A, JacobiOutput *jacobi_output);

#endif