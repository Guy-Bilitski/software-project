#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>

int initialize_centroids(int k, int dim, char *input_file, double**);
int get_dimension(char *input_file);
double max_distance_between_centroids(int k, int dim, double **old_centroids, double **new_centroids);
int kmeans_iteration(int k, int dim, char *input_file, double **centroids, double **new_centroids);
int find_closest_centroid(int k, int dim, double *vector, double **centroids);
void initarray(int dim, int k, double **arr);
int write_result(int k, int dim, char *outname, double **data);
int checkForZeros(int k, int dim, double **centroids);
void printcentroids(int k, int d, double **c);


int main(int argc, char **argv)
{
    char *output_filename, *input_filename;
    int k, maxiter, dim, initialize_error, error;
    double **centroids, **new_centroids, **temp;
    double maxd;
    int i;

    k = atoi(argv[1]);
    input_filename = argv[argc-2];
    output_filename = argv[argc-1];
    maxiter = argc==5 ? atoi(argv[2]) : INT_MAX;
    if (maxiter <= 0 || k <= 0 || argc < 4 || argc > 5) {
        printf("Invalid Input!");
        return 1;
    }

    dim = get_dimension(input_filename);
    if (dim == -1){
        printf("An Error Has Occurred");
        return 1;
    }

    centroids = calloc(k, sizeof(double *));
    new_centroids = calloc(k, sizeof(double *));
    if (centroids == NULL || new_centroids == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    
    for (i=0; i < k; i++) {
        centroids[i] = calloc(dim + 1, sizeof(double));
        new_centroids[i] = calloc(dim + 1, sizeof(double));
        if (centroids[i] == NULL || new_centroids[i] == NULL){
            printf("An Error Has Occurred");
            return 1;
        }
    }

    initialize_error = initialize_centroids(k, dim, input_filename, centroids);
    if (initialize_error == 1){
        printf("An Error Has Occurred");
        return 1;
    } else if (initialize_error == 2){
        printf("Invalid Input!");
        return 1;
    }

    for (i=0; i < maxiter; i++) {
        error = kmeans_iteration(k, dim, input_filename, centroids, new_centroids);
        if (error || checkForZeros(k, dim, new_centroids)){
            printf("An Error Has Occurred");
            return 1;
        }
        maxd = max_distance_between_centroids(k, dim, centroids, new_centroids);
        temp = &centroids[0];
        centroids = &new_centroids[0];
        new_centroids = &temp[0];
        if (maxd < 0.001) {
            break;
        }
        initarray(dim, k, new_centroids);
    }

    error = write_result(k, dim, output_filename, centroids);
    if (error){
        printf("An Error Has Occurred");
        return 1;
    }
    
    for (i=0; i < k; i++) {
        free(centroids[i]);
        free(new_centroids[i]);
    }
    free(centroids);
    free(new_centroids);
    return 0;
}


double max_distance_between_centroids(int k, int dim, double **old_centroids, double **new_centroids) {
    double max_value = DBL_MIN;
    double current_value;

    int i,j;
    for (i=0; i < k; i++) {
        current_value = 0;
        for (j=0; j < dim; j++) {
            current_value += pow((old_centroids[i][j] / old_centroids[i][dim]) - (new_centroids[i][j] / new_centroids[i][dim]), 2);
        }
        if (current_value > max_value) {
            max_value = current_value;
        }
    } 

    return sqrt(max_value);
}


int kmeans_iteration(int k, int dim, char *input_file, double **centroids, double **new_centroids) {
    FILE *ifp;
    int closet_centroid_index, end;
    double *vector;
    char c;
    vector = calloc(dim, sizeof(double));
    if (vector == NULL){
        return 1;
    }

    ifp = fopen(input_file, "r");
    if (ifp == NULL) {
        return 1;
    }

    while (!feof(ifp)) {
        int j;
        for (j=0; j < dim; j++) {
            end = fscanf(ifp, "%lf%c", &vector[j], &c);
        }
        if (end != 2){
            break;
        }

        closet_centroid_index = find_closest_centroid(k, dim, vector, centroids);

        for (j = 0; j < dim; j++) {
            new_centroids[closet_centroid_index][j] += vector[j];
        }

        new_centroids[closet_centroid_index][dim] ++;
    }
    free(vector);
    fclose(ifp);
    return 0;
}

int checkForZeros(int k, int dim, double **centroids){
    int i;
    for (i = 0; i < k; i++) {
        if (centroids[i][dim] == 0){
            return 1;
        }
    }
    return 0;
}

int find_closest_centroid(int k, int dim, double *vector, double **centroids) {
    double closest_value = DBL_MAX;
    double current_value = 0;
    int closest_index = -1;

    int i,j;
    for (i = 0; i < k; i++) {
        current_value = 0;
        for (j=0; j < dim; j++) {
            current_value += pow((vector[j] - (centroids[i][j] / centroids[i][dim])), 2);
        }
        if (current_value < closest_value) {
            closest_value = current_value;
            closest_index = i;
        }
    }

    return closest_index;
}


int initialize_centroids(int k, int dim, char *input_file, double **datapoints) {
    FILE *ifp;
    int i,j;
    char c;
    double check;
    int fsc;
    ifp = fopen(input_file, "r");
    if (ifp == NULL) {
        return 1;
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            fscanf(ifp, "%lf%c", &datapoints[i][j],&c);
        }
        datapoints[i][dim] = 1;
    }

    rewind(ifp);
    for (i = k; i < INT_MAX; i++) {
        
        for (j = 0; j < dim; j++) {
            fsc = fscanf(ifp, "%lf%c", &check,&c);
            if (fsc < 2){
                break;
            }
            if (j < dim-1 && c != ','){
                fclose(ifp);
                return 2;
            }
            else if (j == dim-1 && !isspace(c)){
                fclose(ifp);
                return 2;
            }
        }
        if (fsc < 2){
            break;
        }
    }

    fclose(ifp);
    return 0;
}

void initarray(int dim, int k, double **arr){
    int i,j;
    for (i=0; i<k; i++){
        for (j=0; j<dim+1; j++){
            arr[i][j]=0;
        }
    }
}

int get_dimension(char *input_file) {
    FILE *ifp = NULL;
    char c;
    int d=1;
    ifp = fopen(input_file, "r");
    if (ifp == NULL) {
        return -1;
    }
    while ((c = fgetc(ifp)) != EOF){
        if (c == ','){
            d += 1;
        } else if (c == '\n')
        {
            break;
        }
    }
    fclose(ifp);

    if (c == EOF){
        return -1;
    }
    return d;
}


int write_result(int k, int dim, char *outname, double **data){
    FILE *ofp;
    int i,j;
    ofp = fopen(outname, "w");
    if (ofp == NULL) {
        return 1;
    }
    for (i = 0; i<k; i++){
        for (j=0; j<dim; j++){
            fprintf(ofp, "%.4f",data[i][j]/data[i][dim]);
            if (j < dim-1){
                fprintf(ofp, ",");
            } else {
                fprintf(ofp, "\n");
            }
        }
    }
    fclose(ofp);
    return 0;
}

void printcentroids(int k, int d, double **c){
    int i,j;
    for (i=0; i<k; i++){
        for (j=0; j<d; j++){
            printf("%.4f,", c[i][j]);
        }
        printf(" numOfVecs=%.4f\n", c[i][d]);
    }
}