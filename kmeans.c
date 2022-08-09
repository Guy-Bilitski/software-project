#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix.c"
#include <limits.h>


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
    dim = centroids->cols;
    k = centroids->rows;
    maxiter = maxiter == -1 ? INT_MAX: maxiter;
    new_centroids = create_matrix(dim, k, 0);
    reset_matrix_entries(new_centroids);


    for (iter=0; iter < maxiter; iter++) {
        kmeans_iteration(data_points, centroids, new_centroids);
        max_distance = max_distance_between_centroids(centroids, new_centroids);

        temp = centroids;
        centroids = new_centroids;
        new_centroids = temp;
        
        if (max_distance < epsilon) {
            break;
        }
        reset_matrix_entries(new_centroids);
    }
    return centroids;
}




double max_distance_between_centroids(Matrix *old_centroids, Matrix *new_centroids) {
    int dim, k;
    int r, c;

    double max_value = DBL_MIN;
    double current_value;
    double old_denomin, new_denomin;

    dim = new_centroids->cols;
    k = new_centroids->rows;

    for (r=0; r < k; r++) {
        current_value = 0.;
        for (c=0; c < dim; c++) {
            current_value += pow(
                matrix_get(old_centroids, r, c) - matrix_get(new_centroids, r, c), 2
                );
        }
        if (current_value > max_value) {
            max_value = current_value;
        }
    }

    return sqrt(max_value);
}


void kmeans_iteration(Matrix *data_points , Matrix *centroids, Matrix *new_centroids) {
    int dim, n;
    int r, c;
    int closet_centroid_index;
    double entry_value;
    
    dim = data_points->cols;
    n = data_points->rows;


    for (i=0; i<n; i++) {
        current_vector = PyList_GetItem(data_points, i);
        closet_centroid_index = find_closest_centroid(current_vector, centroids);
        closest_centroid = PyList_GetItem(new_centroids, closet_centroid_index);

        for (j = 0; j < dim; j++) {
            entry_value = PyFloat_AsDouble(PyList_GetItem(closest_centroid, j)) + PyFloat_AsDouble(PyList_GetItem(current_vector, j));
            if (PyList_SetItem(closest_centroid, j, PyFloat_FromDouble(entry_value))){
                printf("An Error Has Occurred\n");
                exit(1);
            }
        }
        entry_value = PyFloat_AsDouble(PyList_GetItem(closest_centroid, dim)) + 1.;
        PyList_SetItem(closest_centroid, dim, PyFloat_FromDouble(entry_value));
    }
}


int find_closest_centroid(Point *vector, Matrix *centroids) {
    int dim, k;
    int centroid_idx, entry;
    dim = centroids->dim;
    k = centroids->k;

    double min_distance = DBL_MAX;
    int min_index = -1;
    double current_distance = 0;

    for (centroid_idx = 0; centroid_idx < k; centroid_idx++) {
        current_distance = 0;
        for (entry=0; entry < dim; entry++) {
            current_value += pow(
                point_get_index(vector, entry) - 
                get_centroid_entry(centroids, centroid_idx, entry), 2);
        }
        if (current_distance < min_distance) {
            min_distance = current_distance;
            min_index = centroid_idx;
        }
    }

    return closest_index;
}




static PyObject* kmeans_capi(PyObject *self, PyObject *args){
    PyObject *data_points;
    PyObject *init_centroids;
    int maxiter;
    double epsilon;

    if (!(PyArg_ParseTuple(args, "OOid", &data_points, &init_centroids, &maxiter, &epsilon))){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    return Py_BuildValue("O", kmeans(data_points, init_centroids, maxiter, epsilon));
}

static PyMethodDef capiMethods[] = {
    {
        "fit",
        (PyCFunction) kmeans_capi,
        METH_VARARGS,
        PyDoc_STR("Args:\nData-Points: ndarray,\nCentroids: list[list]\nmaxiter: int\nepsilon: float")
    },
    {
        NULL, NULL, 0, NULL
    }
};

static struct PyModuleDef moduledef =
{
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    "My Kmeans Module",
    -1,
    capiMethods
};


PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;
    m=PyModule_Create(&moduledef);
    if (!m) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return m;
}