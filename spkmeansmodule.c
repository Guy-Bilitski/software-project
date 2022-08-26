#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


//TODO: Check if prototypes declaration is needed

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


