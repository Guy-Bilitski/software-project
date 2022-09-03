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
    int is_diag;
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


#ifndef YACOBI_OUTPUT_IS_DEFINED
#define YACOBI_OUTPUT_IS_DEFINED
typedef struct YacobiOutput
{
    Matrix *A;
    Matrix *V;
    int k;
} YacobiOutput;
#endif


/* -------------------- POINT PROTOTYPES -------------------- */

/* Point API */
Point *create_empty_point();
Point *create_point(double *data, int dim, int offset);  /* creates a point from list */

int _convert_point_index(Point *point, int index);  /* converts given index to the real one considering the offset */

/* getters */
double point_get_entry(Point *point, int entry);  /* returns the value in index <index> of point */
int point_get_dim(Point *point);  /* returs the point dimension */
int point_get_offset(Point *point);  /* returs the point offset */
double *point_get_data(Point *point);  /* returs the point data */
void divide_point_by_value(Point *p, double value);

/* setters */
void point_set_entry(Point *point, int index, double value);  /* sets value in index <index> */

double inner_product(Point *row_point, Point *column_point);  /* returns row X column scalar */
double euclidean_distance(Point *p1, Point *p2);  /* returns the euclidian distance between two points */
double sum_point_values(Point *point); /* returns the sum of the points values */
double euclidean_norm(Point *p);

/* debugging functions */
void print_point(Point *point);





/* -------------------- MATRIX PROTOTYPES -------------------- */

/* Matrix API */
/* create */
Matrix *create_matrix(int rows, int cols);
Matrix *create_diag_matrix(int n);  /* creates a matrix with dimensions rows X columns all zeros */
Matrix *create_identity_matrix(int n);

/* getters */
int matrix_get_rows_num(Matrix *matrix);
int matrix_get_cols_num(Matrix *matrix);
double *matrix_get_data(Matrix *matrix);
double matrix_get_entry(Matrix *matrix, int row, int col);  /* returns the (row, col) entry */
void matrix_get_row_to_point(Matrix *matrix, Point *point, int row_index); /* inserts a row into a given point */
void matrix_get_column_to_point(Matrix *matrix, Point *point, int column_index);  /* inserts a column into a given point */
void matrix_get_non_diagonal_max_absolute_value(Matrix *matrix, MaxElement *max_element);  /* returns the max element of the matrix */

/* setters */
void matrix_set_entry(Matrix *matrix, int row, int col, double value);  /* sets <value> in (row, col) entry */
void matrix_set_row(Matrix *matrix, int row_index, Point *point);  /* sets a point (row) in the matrix in <row_index> */

/* utilities */
int check_if_matrix_is_diagonal(Matrix *matrix);  /* goes through matrix non diagonal and check if one of them is not 0 */
double matrix_get_row_sum(Matrix *matrix, int row_index);  /* returns <row_index> row sum of values */
void free_matrix(Matrix *matrix); /* cleanup matrix object and sub-objects */
void matrix_add_point_to_row(Matrix *matrix, int row_index, Point *point); /* TODO: check if relevant */
void reset_matrix_entries_to_zero(Matrix *matrix);  /* resets all metrix entries to zero */
Matrix *multiply_matrices(Matrix *m1, Matrix *m2);  /* multiply m1 X m2 and returns the new matrix */
Matrix *sub_matrices(Matrix *A, Matrix *B); /* sub A - B */

/* Matrix inner functions */
int _is_matrix_diag(Matrix *matrix);
int _get_matrix_index(Matrix *matrix, int row, int col);
void _diag_to_square_matrix(Matrix *matrix);
Matrix *_multiply_matrices_diag_with_diag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_diag_with_nondiag(Matrix *m1, Matrix *m2);
Matrix *_multiply_matrices_nondiag_with_nondiag(Matrix *m1, Matrix *m2);

/* debugging functions */
void print_matrix(Matrix *matrix);
void print_matrix_diag(Matrix *matrix);
Matrix *generate_symmetric_matrix(int n);
void space();




/* -------------------- EIGENVECTOR PROTOTYPES -------------------- */

Eigenvector *create_empty_eigen_vector();
Point *eigen_vector_get_point(Eigenvector *eigen_vector);
double eigen_vector_get_eigen_value(Eigenvector *eigen_vector);
Eigenvector *create_eigen_vector(Point *point, double eigen_value);
int compare_eigenvectors(const void *p1, const void *p2);
void sort_eigenvectors_array(Eigenvector *array, size_t n);
void print_eigen_vectors_array(Eigenvector **eigen_vectors_array, int n);
void free_eigen_vector(Eigenvector *eigen_vector);
void free_eigen_vectors_array(Eigenvector **eigen_vectors_array, int n);



/* -------------------- S AND C PROTOTYPES -------------------- */

/* S AND C API */
S_and_C *create_empty_S_and_C();

/* getters */
double s_and_c_get_s(S_and_C *s_and_c);
double s_and_c_get_c(S_and_C *s_and_c);

/* setters */
void S_and_C_set_values(S_and_C *s_and_c, double s, double c);



/* -------------------- MAX ELEMENT PROTOTYPES -------------------- */

/* MaxElemnt API */
MaxElement *create_empty_max_element();
MaxElement *create_max_element(double value, int i, int j);

/* getters */
double max_element_get_value(MaxElement *max_element);
int max_element_get_index1(MaxElement *max_element);
int max_element_get_index2(MaxElement *max_element);

/* setters */
void max_element_set_new_values(MaxElement *max_element, double value, int i, int j);
void max_element_set_value(MaxElement *max_element, double value);
void max_element_set_index1(MaxElement *max_element, int i);
void max_element_set_index2(MaxElement *max_element, int j);

/* debugging */
void print_max_element(MaxElement *max_element);


/* -------------------- YACOBI OUTPUT PROTOTYPES -------------------- */

YacobiOutput *create_empty_yacobi_output();
Matrix *yacobi_output_get_A(YacobiOutput *yacobi_output);
Matrix *yacobi_output_get_V(YacobiOutput *yacobi_output);
void set_yacobi_output_values(YacobiOutput *yacobi_output, Matrix *A, Matrix *V);
void free_yacobi_output(YacobiOutput *yacobi_output);


/* -------------------- KMEANS-IO PROTOTYPES -------------------- */

int get_dimension(const char *input_file);
int get_n(const char *input_file);
Matrix *input_file_to_matrix(const char *input_file);

/* -------------------- KMEANS PROTOTYPES -------------------- */

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

/* JACOBI */
MaxElement *get_off_diagonal_absolute_max(Matrix *matrix);
void get_s_and_c_for_rotation_matrix(Matrix* A, MaxElement *max_element, S_and_C *s_and_c);
void build_rotation_matrix(S_and_C *s_and_c, MaxElement *max_element, int dim, Matrix *identity_matrix); /* inserts the rotation matrix to a given identity matrix */
void normalize_matrix_rows(Matrix *matrix);
double matrix_off(Matrix *matrix); /* returns the value of "off" function on a given matrix */
Matrix *transform_matrix(Matrix *matrix, S_and_C *s_and_c, MaxElement *max_element);  /* permorms matrix transformation */
void normalize_matrix_rows(Matrix *matrix);
int get_k_from_sorted_eigenvectors_array(Eigenvector **eigen_vectors_array, int n);
int get_k_from_yacobi_output(YacobiOutput *yacobi_output);
Eigenvector **get_eigen_vectors_from_yacobi_output(YacobiOutput *yacobi_output);
Matrix *getU(YacobiOutput *yacobi_output, int k);

/* utilities */
double get_value_for_transformed_matrix(Matrix *old_matrix, double s, double c, int i, int j, int row_index, int col_index); /* returns the expected value of the transformed matrix at (row_index, col_index) based on the rules described at 6. Relations betweeb A and A'*/
int matrix_converge(double A_off, Matrix *A);


/* SPKMEANS API */
Matrix *wam(Matrix* data_points);
Matrix *ddg(Matrix* data_points);
Matrix *lnorm(Matrix* data_points);
YacobiOutput *jacobi(Matrix *A, YacobiOutput *yacobi_output);



/* JACOBI */
void print_jacobi_output(YacobiOutput *J);
void free_yacobi_output(YacobiOutput *yacobi_output);
void set_yacobi_output_values(YacobiOutput *yacobi_output, Matrix *A, Matrix *V, int k);
YacobiOutput *create_empty_yacobi_output();



#endif