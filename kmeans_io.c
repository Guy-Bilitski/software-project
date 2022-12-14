#include "spkmeans.h"
#include "ctype.h"


int get_dimension(const char *input_file) {
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


int get_n(const char *input_file) {
    FILE *ifp;
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
        if (isspace(delimiter)) {
            n += 1;
        }
    }

    fclose(ifp);
    return n;
}


Matrix *input_file_to_matrix(const char *input_file) {
    FILE *ifp;
    int i,j;
    char delimiter;
    double element;
    int elements_count = 0;
    int dim, n;
    Matrix *data_points;

    dim = get_dimension(input_file);
    n = get_n(input_file);
    data_points = create_matrix(n, dim);

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