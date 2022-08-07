#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix.c"


PyObject * kmeans(PyObject *data_points, Centroids *centroids, int maxiter, double epsilon);
double max_distance_between_centroids(PyObject *old_centroids, PyObject *new_centroids);
void kmeans_iteration(PyObject *data_points , PyObject *centroids, PyObject *new_centroids);
void initarray(PyObject *matrix);
int checkForZeros(int k, int dim, double **centroids);
void print_pymatrix(PyObject *matrix);
void print_centroids(PyObject *matrix);
void fix_final_centroids_matrix(PyObject *matrix);


/* Already changed:*/
int find_closest_centroid(Point *vector, Centroids *centroids);


PyObject * kmeans(PyObject *data_points, Centroids *centroids, int maxiter, double epsilon)
{
    Centroids *new_centroids, *temp;
    int dim, k;
    int iter;
    double maxd;

    //Setting variables
    dim = centroids->dim;
    k = centroids->k;
    maxiter = maxiter == -1 ? INT_MAX: maxiter;
    new_centroids = init_centroids(dim, k)


    for (iter=0; iter < maxiter; iter++) {
        kmeans_iteration(data_points, centroids, new_centroids);
        maxd = max_distance_between_centroids(centroids, new_centroids);

        temp = centroids;
        centroids = new_centroids;
        new_centroids = temp;
        
        if (maxd < epsilon) {
            break;
        }
        initarray(new_centroids);
    }
    fix_final_centroids_matrix(centroids);

    return centroids;
}




double max_distance_between_centroids(PyObject *old_centroids, PyObject *new_centroids) {
    Py_ssize_t dim, k;
    Py_ssize_t i,j;
    PyObject *old_c, *new_c;

    double max_value = DBL_MIN;
    double current_value;
    double old_denomin, new_denomin;

    dim = PyList_Size(PyList_GetItem(new_centroids, 0))-1;
    k = PyList_Size(new_centroids);

    for (i=0; i < k; i++) {
        current_value = 0.;
        old_c = PyList_GetItem(old_centroids, i);
        new_c = PyList_GetItem(new_centroids, i);
        for (j=0; j < dim; j++) {
            old_denomin = PyFloat_AsDouble(PyList_GetItem(old_c, dim));
            new_denomin = PyFloat_AsDouble(PyList_GetItem(new_c, dim));
            if (old_denomin == 0 || new_denomin == 0){
                printf("An Error Has Occurred\n");
                exit(1);
            }
            current_value += pow(
                (PyFloat_AsDouble(PyList_GetItem(old_c, j)) / old_denomin) - 
                (PyFloat_AsDouble(PyList_GetItem(new_c, j)) / new_denomin), 2
                );
        }
        if (current_value > max_value) {
            max_value = current_value;
        }
    }

    return sqrt(max_value);
}


void kmeans_iteration(PyObject *data_points , PyObject *centroids, PyObject *new_centroids) {
    Py_ssize_t dim, n;
    Py_ssize_t i,j;
    Py_ssize_t closet_centroid_index;
    PyObject *closest_centroid, *current_vector;
    double entry_value;
    
    dim = PyList_Size(PyList_GetItem(data_points, 0));
    n = PyList_Size(data_points);


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


int find_closest_centroid(Point *vector, Centroids *centroids) {
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


void initarray(PyObject *junk_centroids){
    Py_ssize_t i,j;
    Py_ssize_t dim, k;
    PyObject * centroid_i;
    dim = PyList_Size(PyList_GetItem(junk_centroids, 0))-1;
    k = PyList_Size(junk_centroids);

    for (i=0; i<k; i++){
        centroid_i = PyList_GetItem(junk_centroids, i);
        for (j=0; j<dim+1; j++){
            PyList_SetItem(centroid_i, j, PyFloat_FromDouble(0.));
        }
    }
}



void print_pymatrix(PyObject *matrix){
    Py_ssize_t i, j;
    Py_ssize_t dim;
    Py_ssize_t k;
    PyObject *pyfloat;

    k = PyList_Size(matrix);
    dim = PyList_Size(PyList_GetItem(matrix, 0));
    for (i=0; i < k; i++) {
        for (j=0; j < dim; j++){
            pyfloat = PyList_GetItem(PyList_GetItem(matrix, i), j);
            printf("%f,",PyFloat_AsDouble(pyfloat));
        }
        printf("\n");
    }
}

void print_centroids(PyObject *matrix){
    Py_ssize_t i, j;
    Py_ssize_t dim;
    Py_ssize_t k;
    PyObject *pyfloat;
    double num;

    k = PyList_Size(matrix);
    dim = PyList_Size(PyList_GetItem(matrix, 0));
    for (i=0; i < k; i++) {
        num = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(matrix, i), dim-1));
        for (j=0; j < dim-1; j++){
            pyfloat = PyList_GetItem(PyList_GetItem(matrix, i), j);
            printf("%f,",PyFloat_AsDouble(pyfloat)/num);
        }
        printf("\n");
    }
}

void fix_final_centroids_matrix(PyObject *matrix){
    Py_ssize_t i, j;
    Py_ssize_t dim;
    Py_ssize_t k;
    PyObject *pyfloat, *sublist;
    double num, entry;

    k = PyList_Size(matrix);
    dim = PyList_Size(PyList_GetItem(matrix, 0));

    for (i=0; i < k; i++) { //dividing each entry by sum
        num = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(matrix, i), dim-1));
        for (j=0; j < dim-1; j++){
            pyfloat = PyList_GetItem(PyList_GetItem(matrix, i), j);
            entry = PyFloat_AsDouble(pyfloat)/num;
            PyList_SetItem(PyList_GetItem(matrix, i), j, PyFloat_FromDouble(entry));
        }
    }

    for (i=0; i < k; i++) { //popping the sum entry
        sublist = PyList_GetSlice(PyList_GetItem(matrix, i), 0, dim-1);
        PyList_SetItem(matrix, i, sublist);
    }
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