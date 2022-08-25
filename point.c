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


/* Point API */
Point *create_point(double *data, int dim, int offset) {
    Point *new_point = (Point *)malloc(sizeof(Point));
    new_point->data = data;
    new_point->dim = dim;
    new_point->offset = offset;
    return new_point;
}

int _convert_point_index(Point *point, int index) {
    int point_offset = point_get_offset(point);
    return point_offset*index;
}

double point_get_entry(Point *point, int index) {
    assert(index >= 0 && index < point->dim);
    int real_index = _convert_point_index(point, index);
    double *data = point_get_data(point);
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

void point_set_entry(Point *point, int index, double value) {
    assert(index >= 0 && index < point->dim);
    int real_index = _convert_point_index(point, index);
    double *data = point_get_data(point);
    data[real_index] = value;
}

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

double sum_point_values(Point *point) {
    int i;
    double sum = 0;
    for (i=0; i<point_get_dim(point); i++) {
        sum += point_get_entry(point, i);
    }
    return sum;
}

void divide_point_by_value(Point *p, double value) {
    int i, n, index;
    n = point_get_dim(p); 

    for (i=0; i<n; i++) {
        index = _convert_point_index(p, i);
        (p->data)[index] /= value;
    }
}


/* debugging function */
void print_point(Point *point) {
    int i;
    for (i=0; i<(point->dim); i++) {
        printf("%f ", point_get_entry(point, i));
    }
}


int main() {
    int i;
    int dim = 9;
    double *data = (double *)calloc(sizeof(double), dim);
    for (i=0; i<dim; i++) {
        data[i] = i*2;
    }
    Point *p = create_point(data, dim/3, 3);
    print_point(p);
    printf("\n");
    printf("\n");
    divide_point_by_value(p, 2);
    print_point(p);
    printf("\n");
    printf("\n");
}