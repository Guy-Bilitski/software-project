#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.h"

#define true 1
#define false 0

#ifndef POINT_IS_DEFINED
#define POINT_IS_DEFINED
typedef struct Point {
    double *data;
    int dim;
    int offset; /* for column representation */
} Point;
#endif


/* Point API */
Point *create_empty_point() {
    Point *new_point = (Point *)malloc(sizeof(Point));
    return new_point;
}

Point *create_point(double *data, int dim, int offset) {
    Point *new_point = create_empty_point();
    new_point->data = data;
    new_point->dim = dim;
    new_point->offset = offset;
    return new_point;
}

/* Getters */
double point_get_entry(Point *point, int index) {
    double *data;
    int real_index;
    assert(index >= 0 && index < point->dim);
    real_index = _convert_point_index(point, index);
    data = point_get_data(point);
    return data[real_index];
}

int point_get_dim(Point *point) {
    return point->dim;
}

int point_get_offset(Point *point) {
    return point->offset;
}

double *point_get_data(Point *point) {
    return point->data;
}

/* Utils */
double inner_product(Point *row_point, Point *column_point) {
    int i, points_dim = point_get_dim(row_point);
    double sum = 0;
    for (i=0; i<points_dim; i++) {
        sum += (point_get_entry(row_point, i) * point_get_entry(column_point, i));
    }
    return sum;
}

double euclidean_norm(Point *p) {
    double sum = 0.;
    int i, dim = point_get_dim(p);
    for (i=0; i<dim; i++) {
        sum += pow(point_get_entry(p, i), 2);
    }
    return sqrt(sum);
}

double euclidean_distance(Point *p1, Point *p2) {
    double sum = 0.;
    int i, dim = point_get_dim(p1);
    for (i=0; i<dim; i++) {
        sum += pow(point_get_entry(p1, i) - point_get_entry(p2, i), 2);
    }
    return sqrt(sum);
}

void divide_point_by_value(Point *p, double value) {
    int i, n, index;
    n = point_get_dim(p); 

    for (i=0; i<n; i++) {
        index = _convert_point_index(p, i);
        point_get_data(p)[index] /= value;
    }
}

/* Point inner functions */
int _convert_point_index(Point *point, int index) {
    int point_offset = point_get_offset(point);
    return point_offset*index;
}

/* debugging function */
void print_point(Point *point) {
    int i;
    for (i=0; i<(point->dim); i++) {
        printf("%f ", point_get_entry(point, i));
    }
}
