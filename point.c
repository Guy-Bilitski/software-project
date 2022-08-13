#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define true 1
#define false 0

typedef struct Point {
    double *data;
    int dim;
    int offset; /* for column representation */
} Point;


/* Point API */
Point *create_point(double *data, int dim, int offset);  /* creates a point from list */

int _convert_point_entry(Point *point, int entry);  /* converts given entry to the real one considering the offset */

double point_get_entry(Point *point, int entry);  /* returns the value in index <index> of point */
int point_get_dim(Point *point);  /* returs the point dimension */
int point_get_offset(Point *point);  /* returs the point offset */
double *point_get_data(Point *point);  /* returs the point data */

void point_set_entry(Point *point, int index, double value);  /* sets value in index <index> */

double inner_product(Point *row_point, Point *column_point);  /* returns row X column scalar */
double euclidean_distance(Point *p1, Point *p2);  /* returns the euclidian distance between two points */
double sum_point_values(Point *point); /* returns the sum of the points values */

/* debugging functions */
void print_point(Point *point);


/* Point API */
Point *create_point(double *data, int dim, int offset) {
    Point *new_point = (Point *)malloc(sizeof(Point));
    new_point->data = data;
    new_point->dim = dim;
    new_point->offset = offset;
    return new_point;
}

int _convert_point_entry(Point *point, int entry) {
    int point_offset = point_get_offset(point);
    if (point_offset == 0) {
        return entry;
    } else {
        return point_offset*entry;
    }
}

double point_get_entry(Point *point, int entry) {
    assert(entry >= 0 && entry < point->dim);
    int real_entry = _convert_point_entry(point, entry);
    double *data = point_get_data(point);
    return data[real_entry];
    
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

void point_set_entry(Point *point, int entry, double value) {
    assert(entry >= 0 && entry < point->dim);
    int real_entry = _convert_point_entry(point, entry);
    double *data = point_get_data(point);
    data[real_entry] = value;
}

double inner_product(Point *row_point, Point *column_point) {
    int i, points_dim = point_get_dim(row_point);
    double sum = 0;
    for (i=0; i<points_dim; i++) {
        sum += (point_get_entry(row_point, i) * point_get_entry(column_point, i));
    }
    return sum;
}

double euclidean_distance(Point *p1, Point *p2) {
    double sum = 0.;
    int i, dim = point_get_dim(p1);
    for (i=0; i<dim; i++) {
        sum += pow(point_get_entry(p1, i) - point_get_entry(p2, i), 2);
    }
    return sqrt(sum);
}

double sum_point_values(Point *point) {
    int i;
    double sum = 0;
    for (i=0; i<point_get_dim(point); i++) {
        sum += point_get_entry(point, i);
    }
    return sum;
}


/* debugging function */
void print_point(Point *point) {
    int i;
    for (i=0; i<(point->dim); i++) {
        printf("%f ", point_get_entry(point, i));
    }
}


int main3() {
    int i;
    int dim = 10;
    double *data = (double *)calloc(sizeof(double), dim);
    for (i=0; i<dim; i++) {
        data[i] = i*2 + 1;
    }
    Point *p = create_point(data, dim, 0);
    print_point(p);
    printf("\n");
    printf("\n");
    Point *p2 = create_point(data, 5, 2);
    print_point(p2);
}