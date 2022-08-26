#include "spkmeans.h"


int main1 (int argc, char **argv) {
    const char *goal_strings[invalid - spk + 1] = { "spk", "wam", "ddg", "lnorm", "jacobi", "invalid" };
    int dim, n;
    char *input_filename;
    enum goal goal;
    Matrix *data_points, *output;

    if (argc != 3) {
        printf("Invalid Input!\n");
        exit(1);
    }
    
    goal = get_goal(argv[1]);
    input_filename = argv[2];

    dim = get_dimension(input_filename);
    n = get_n(input_filename);
    data_points = input_file_to_matrix(input_filename, dim, n);

    output = achieve_goal(data_points, goal);
    print_matrix(output);
    return 1;
}



int get_dimension(char *input_file) {
    FILE *ifp = NULL;
    char c;
    int dim=1;
    ifp = fopen(input_file, "r");
    if (ifp == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    while ((c = fgetc(ifp)) != EOF){
        if (c == ','){
            dim += 1;
        } else if (c == '\n')
        {
            break;
        }
    }
    fclose(ifp);

    if (c == EOF){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return dim;
}


int get_n(char *input_file) {
    FILE *ifp;
    int i,j;
    char delimiter;
    double temp;
    int n = 0;

    ifp = fopen(input_file, "r");
    if (ifp == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    while (1) {
        if (fscanf(ifp, "%lf%c", &temp,&delimiter) < 2){
            break;
        }
        if (delimiter == '\n') {
            n += 1;
        }
    }

    fclose(ifp);
    return n;
}


enum goal get_goal(char *goal) {
    if (!strcmp(goal, "wam")){
        return wam;
    }
    else if (!strcmp(goal, "ddg")){
        return ddg;
    }
    else if (!strcmp(goal, "lnorm")){
        return lnorm;
    }
    else if (!strcmp(goal, "jacobi")){
        return jacobi;
    }
    else {
        printf("Invalid Input!\n");
        exit(1);
    }
}



Matrix *input_file_to_matrix(char *input_file, int dim, int n) {
    FILE *ifp;
    int i,j;
    char delimiter;
    double element;
    int elements_count = 0;
    Matrix *data_points = create_matrix(n, dim);

    ifp = fopen(input_file, "r");
    if (ifp == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i=0; i<n; i++) {
        for (j=0; j<dim; j++) {
            if (fscanf(ifp, "%lf%c", &element,&delimiter) < 2){
                break;
            }
            matrix_set_entry(data_points, i, j, element);
            elements_count++;
        }
    }

    if (elements_count != dim*n) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    fclose(ifp);
    return data_points;
}


Matrix *achieve_goal(Matrix *data_points, enum goal goal) {
    Matrix *W, *D, *Lnorm;

    switch (goal)
    {
        case wam:
            W = create_weighted_matrix(data_points);
            return W;
        
        case ddg:
            W = create_weighted_matrix(data_points);
            D = create_diagonal_degree_matrix(W);
            free_matrix(W);
            return D;

        case lnorm:
            W = create_weighted_matrix(data_points);
            D = create_diagonal_degree_matrix(W);
            neg_root_to_diag_matrix(D);
            Lnorm = normalized_graph_laplacian(D, W);
            free_matrix(W);
            free_matrix(D);
            return Lnorm;
        
    }
}