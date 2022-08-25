#include <float.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include "spkmeans.h"


Matrix * kmeans(Matrix *data_points, Matrix *centroids, int maxiter, double epsilon);
double max_distance_between_centroids(Matrix *old_centroids, Matrix *new_centroids);
void kmeans_iteration(Matrix *data_points , Matrix *centroids, Matrix *new_centroids);
int find_closest_centroid(Point *vector, Matrix *centroids);



Matrix * kmeans(Matrix *data_points, Matrix *centroids, int maxiter, double epsilon)
{
    Matrix *new_centroids, *temp;
    int dim, k;
    int iter;
    double max_distance;

    //Setting variables
    dim = matrix_get_cols_num(centroids);
    k = matrix_get_rows_num(centroids);
    maxiter = maxiter == -1 ? INT_MAX: maxiter;
    new_centroids = create_matrix(k, dim);
    reset_matrix_entries_to_zero(new_centroids);


    for (iter=0; iter < maxiter; iter++) {
        kmeans_iteration(data_points, centroids, new_centroids);
        max_distance = max_distance_between_centroids(centroids, new_centroids);

        temp = centroids;
        centroids = new_centroids;
        new_centroids = temp;
        
        if (max_distance < epsilon) {
            break;
        }
        reset_matrix_entries_to_zero(new_centroids);
    }
    return centroids;
}




double max_distance_between_centroids(Matrix *old_centroids, Matrix *new_centroids) {
    int dim, k;
    int r, c;

    double max_distance = DBL_MIN;
    double current_distance;

    dim = matrix_get_cols_num(new_centroids);
    k = matrix_get_rows_num(new_centroids);

    Point *old_centroid = (Point *)malloc(sizeof(Point));
    Point *new_centroid = (Point *)malloc(sizeof(Point));


    for (r=0; r < k; r++) {
        matrix_get_row_to_point(old_centroids, old_centroid, r);
        matrix_get_row_to_point(new_centroids, new_centroid, r);
        current_distance = euclidean_distance(old_centroid, new_centroid);

        if (current_distance > max_distance) {
            max_distance = current_distance;
        }
    }

    free(old_centroid);
    free(new_centroid);
    return max_distance;
}


void kmeans_iteration(Matrix *data_points , Matrix *centroids, Matrix *new_centroids) {
    int r, c;
    int closet_centroid_index;
    double entry_value;
    int dim = matrix_get_cols_num(data_points);
    int n = matrix_get_rows_num(data_points);
    int k = matrix_get_rows_num(centroids);
    Point *current_vector = (Point *)malloc(sizeof(Point));
    int *num_of_points_in_cluster = (int *)calloc(k, sizeof(int));

    
    for (r=0; r<n; r++) {
        matrix_get_row_to_point(data_points, current_vector, r);
        closet_centroid_index = find_closest_centroid(current_vector, centroids);
        matrix_add_point_to_row(new_centroids, closet_centroid_index, current_vector);
        num_of_points_in_cluster[closet_centroid_index]++;
    }

    for (r=0; r<k; r++){
        assert(num_of_points_in_cluster[r] > 0);
        for (c=0; c<dim; c++){
            entry_value = matrix_get_entry(new_centroids, r, c) / num_of_points_in_cluster[r];
            matrix_set_entry(new_centroids, r, c, entry_value);
        }
    }
    free(current_vector);
    free(num_of_points_in_cluster);
}


int find_closest_centroid(Point *vector, Matrix *centroids) {
    int dim, k;
    int centroid_idx, entry;
    dim = matrix_get_cols_num(centroids);
    k = matrix_get_rows_num(centroids);

    double min_distance = DBL_MAX;
    int min_index = -1;
    double current_distance;
    Point *current_centroid = (Point *)malloc(sizeof(Point));

    for (centroid_idx = 0; centroid_idx < k; centroid_idx++) {
        matrix_get_row_to_point(centroids, current_centroid, centroid_idx);
        current_distance = euclidean_distance(current_centroid, vector);
        if (current_distance < min_distance) {
            min_distance = current_distance;
            min_index = centroid_idx;
        }
    }
    free(current_centroid);
    return min_index;
}
