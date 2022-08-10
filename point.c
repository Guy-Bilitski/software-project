#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define true 1
#define false 0

typedef struct Point {
    int dim;
    double *data;
} Point;


/* Point API */
Point *create_empty_point(int dim);  /* creates a point with all zero values */
Point *create_point(int dim, double *data);  /* creates a point from list */
double point_get_index(Point *point, int index);  /* returns the value in index <index> of point */
int point_get_dim(Point *point);  /* returs the point dimension */
void point_set_index(Point *point, int index, double value);  /* sets value in index <index> */
double multiply_points(Point *row_point, Point *column_point);  /* returns row X column scalar */
double euclidean_distance(Point *p1, Point *p2);  /* returns the euclidian distance between two points */

/* debugging functions */
void print_point(Point *point);


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
    assert(index >= 0 && index < point->dim);
    point->data[index] = value;
}

double multiply_points(Point *row_point, Point *column_point) {
    int i, points_dim = point_get_dim(row_point);
    double sum = 0;
    for (i=0; i<points_dim; i++) {
        sum += (point_get_index(row_point, i) * point_get_index(column_point, i));
    }
    return sum;
}

double euclidean_distance(Point *p1, Point *p2) {
    double sum;
    int i, dim = point_get_dim(p1);
    for (i=0; i<dim; i++) {
        sum += pow(point_get_index(p1, i) - point_get_index(p2, i), 2);
    }
    return sqrt(sum);
}


/* debugging function */
void print_point(Point *point) {
    int i;
    for (i=0; i<(point->dim); i++) {
        printf("%f ", (point->data)[i]);
    }
}



int main3() {
    return 1;
}