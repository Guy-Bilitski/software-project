#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.c"

PyObject *matrix_to_pylist(Matrix *matrix){
    int cols, rows;
    PyObject *py_matrix, *temp;
    int i, j;
    double value;

    cols = matrix_get_cols_num(matrix);
    rows = matrix_get_rows_num(matrix);
    py_matrix = PyList_New(rows);

    if (py_matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i=0; i<rows; i++) {

        temp = PyList_New(cols);
        if (temp == NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }

        PyList_SetItem(py_matrix, i, temp);
        for (j=0; j<cols; j++){
            value = matrix_get_entry(matrix, i, j);
            PyList_SetItem(temp, j, PyFloat_FromDouble(value));
        }
    }
    return py_matrix;
}

PyObject *diagonal_matrix_to_pylist(Matrix *matrix){
    int n;
    PyObject *py_matrix;
    int i;
    double value;

    n = matrix->cols;
    py_matrix = PyList_New(n);

    if (py_matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i=0; i<n; i++) {
        value = matrix_get_entry(matrix, i, i);
        PyList_SetItem(py_matrix, i, PyFloat_FromDouble(value));
    }
    return py_matrix;
}


Matrix *pylist_to_matrix(PyObject *pymatrix) {
    int cols, rows;
    Matrix *matrix;
    int i, j;
    double value;
    PyObject *current_row;

    if (!PyList_Check(pymatrix) || !PyList_Check(PyList_GetItem(pymatrix, 0))){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    rows = PyList_Size(pymatrix);
    cols = PyList_Size(PyList_GetItem(pymatrix, 0));
    matrix = create_matrix(rows, cols);
    
    for (i=0; i<rows; i++){
        current_row = PyList_GetItem(pymatrix, i);
        for (j=0; j<cols; j++){
            value = PyFloat_AsDouble(PyList_GetItem(current_row, j));
            matrix_set_entry(matrix, i, j, value);
        }
    }
    return matrix;
}



static PyObject* wam_capi(PyObject *self, PyObject *args){
    PyObject *data_points_as_pylist;
    Matrix *W, *data_points;
    PyObject *output;

    if (!(PyArg_ParseTuple(args, "O", &data_points_as_pylist))){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    data_points = pylist_to_matrix(data_points_as_pylist);
    W = wam(data_points);
    output = matrix_to_pylist(W);
    
    free_matrix(W);
    free_matrix(data_points);
    return Py_BuildValue("O", output);
}


static PyObject* ddg_capi(PyObject *self, PyObject *args){
    PyObject *data_points_as_pylist;
    Matrix *D, *data_points;
    PyObject *output;

    if (!(PyArg_ParseTuple(args, "O", &data_points_as_pylist))){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    data_points = pylist_to_matrix(data_points_as_pylist);
    D = ddg(data_points);
    output = matrix_to_pylist(D);
    free_matrix(D);
    free_matrix(data_points);
    return Py_BuildValue("O", output);
}

static PyObject* lnorm_capi(PyObject *self, PyObject *args){
    PyObject *data_points_as_pylist;
    Matrix *Lnorm, *data_points;
    PyObject *output;

    if (!(PyArg_ParseTuple(args, "O", &data_points_as_pylist))){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    data_points = pylist_to_matrix(data_points_as_pylist);
    Lnorm = lnorm(data_points);
    output = matrix_to_pylist(Lnorm);
    free_matrix(Lnorm);
    free_matrix(data_points);
    return Py_BuildValue("O", output);
}


static PyObject* jacobi_capi(PyObject *self, PyObject *args) {
    PyObject *data_points_as_pylist;
    JacobiOutput *Jout;
    PyObject *V, *A;
    Matrix *sym_matrix;

    if (!(PyArg_ParseTuple(args, "O", &data_points_as_pylist))){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    sym_matrix = pylist_to_matrix(data_points_as_pylist);
    Jout = create_empty_jacobi_output();
    jacobi(sym_matrix, Jout);
    V = matrix_to_pylist(Jout->V);
    A = diagonal_matrix_to_pylist(Jout->A);
    free_jacobi_output(Jout);
    return Py_BuildValue("OO", V, A);
}


static PyObject* transform_data_points_capi(PyObject *self, PyObject *args) {
    PyObject *data_points_as_pylist;
    int k;
    JacobiOutput *Jout;
    PyObject *pylist_U;
    Matrix *data_points_matrix, *laplacian, *U;
    

    if (!(PyArg_ParseTuple(args, "Oi", &data_points_as_pylist, &k))){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    data_points_matrix = pylist_to_matrix(data_points_as_pylist);
    
    laplacian = lnorm(data_points_matrix);
    Jout = create_empty_jacobi_output();
    jacobi(laplacian, Jout);

    U = getU(Jout, k);  
    normalize_matrix_rows(U);
    pylist_U = matrix_to_pylist(U);
    
    free_jacobi_output(Jout);
    free_matrix(data_points_matrix);
    free_matrix(U);
    
    
    return Py_BuildValue("O", pylist_U);
}

static PyObject* kmeans_capi(PyObject *self, PyObject *args){
    PyObject *py_data_points, *py_init_centroids, *output;
    Matrix *init_centroids, *data_points, *final_centroids;

    if (!(PyArg_ParseTuple(args, "OO", &py_data_points, &py_init_centroids))){
        printf("An Error Has Occurred\n");
        exit(1);
    }

    data_points = pylist_to_matrix(py_data_points);
    init_centroids = pylist_to_matrix(py_init_centroids);

    final_centroids = kmeans(data_points, init_centroids);

    output = matrix_to_pylist(final_centroids);
    
    free_matrix(data_points);
    free_matrix(final_centroids);
    
    return Py_BuildValue("O", output);
}


static PyMethodDef capiMethods[] = {
    {
        "wam",
        (PyCFunction) wam_capi,
        METH_VARARGS,
        PyDoc_STR("wam docs...")
    },
    {
        "ddg",
        (PyCFunction) ddg_capi,
        METH_VARARGS,
        PyDoc_STR("ddg docs...")
    },
    {
        "lnorm",
        (PyCFunction) lnorm_capi,
        METH_VARARGS,
        PyDoc_STR("lnorm docs...")
    },
    {
        "jacobi",
        (PyCFunction) jacobi_capi,
        METH_VARARGS,
        PyDoc_STR("jacobi docs...")
    },
    {
        "transform_data_points",
        (PyCFunction) transform_data_points_capi,
        METH_VARARGS,
        PyDoc_STR("kmeans docs...")
    },
    {
        "kmeans",
        (PyCFunction) kmeans_capi,
        METH_VARARGS,
        PyDoc_STR("kmeans docs...")
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