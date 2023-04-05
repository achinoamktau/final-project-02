#include <Python.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "spkmeans_functions.h"

/* C Api Declarations: */
static double **create2DArrayFromPyObject(PyObject *data, int n, int d);
static PyObject *create2DPyObject(double **data, int n, int d);
#define verifyNotNULL(var) if((var)==NULL) {printf("An Error Has Occured"); exit(-1);}

/*--------------------------------- c-api -----------------------------------------------------*/

static PyObject *wam_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double** wam_matrix = wam(data_points, n, d);
    PyObject* result = create2DPyObject(wam_matrix, n, n);
    // we put here a param that tells if we want to print the matrix or not but i dont like it
    free2DMalloc(data_points, n);
    free_contiguous_mat(wam_matrix);
    return result;
}

static PyObject *ddg_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    double **ddg_matrix = ddg(wam_matrix, n);
    PyObject *result = create2DPyObject(ddg_matrix, n, n);
    free2DMalloc(data_points, n);
    free_contiguous_mat(wam_matrix);
    free_contiguous_mat(ddg_matrix);
    return result;
}
static PyObject *gl_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    double **ddg_matrix = ddg(wam_matrix, n);
    double **gl_matrix = gl(ddg_matrix, wam_matrix, n);   // this what gl_py but also dont like it
    PyObject *result = create2DPyObject(gl_matrix, n, n);
    free2DMalloc(data_points, n);
    free_contiguous_mat(wam_matrix);
    free_contiguous_mat(ddg_matrix);
    free_contiguous_mat(gl_matrix);
    return result;
}

static PyObject *jacobi_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);    // this is a n on n symmetric matrix
    Jacobi_output *jacobi_matrix = jacobi(data_points, n);
    double **jacobi_res = create2DfromJacobi(jacobi_matrix, n);
    PyObject *result = create2DPyObject(jacobi_res, n+1, n);
    free2DMalloc(data_points, n);
    free_contiguous_mat(jacobi_res);
    free(jacobi_matrix->eigenValues);
    free_contiguous_mat(jacobi_matrix->V);
    free(jacobi_matrix);
    return result;
}


// static PyObject *getT_Capi(PyObject *self, PyObject *args) {
//     PyObject *pyPoints;
//     PyObject *pyT;
//     int k, n, d;
//     double **T, **points;
//     if (!PyArg_ParseTuple(args, "iiiO", &k, &n, &d, &pyPoints)) {
//         printf("An Error Has Occured");
//         Py_RETURN_NONE;
//     }

//     points = create2DArrayFromPyObject(pyPoints, n, d);

//     T=getT(points,d,n,&k);
//     pyT = create2DPyObject(T,n,k);

//     free2DMalloc(points, n);
//     free_contiguous_mat(T);
//     return Py_BuildValue("[Oi]", pyT, k);
// }



static PyObject *kmeans_api(PyObject *self, PyObject *args) {
    PyObject *pyInitialCentroids, *pyPoints;
    int k, n, d, max_iter;
    double **initialCentroids, **points;
    if (!PyArg_ParseTuple(args, "iiiiOO", &k, &n, &d, &max_iter, &pyInitialCentroids, &pyPoints)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    initialCentroids = create2DArrayFromPyObject(pyInitialCentroids, k, d);
    points = create2DArrayFromPyObject(pyPoints, n, d);

    double **finalCentroids = kmeans(points, initialCentroids, k, d, n, max_iter);

    print_2d_array(finalCentroids, k, d);

    free2DMalloc(points, n);
    free2DMalloc(initialCentroids, k);

    return Py_True;
}
static double **create2DArrayFromPyObject(PyObject *data, int n, int d) {
    int i, j;
    double **points;
    PyObject *temp_point,*inner_item;

    points = (double **) malloc(n * sizeof(double *));
    verifyNotNULL(points)

    for (i = 0; i < n; i++) {
        double *vector = malloc(d * sizeof(double));
        verifyNotNULL(vector)


        temp_point = PyList_GetItem(data, i);
        for (j = 0; j < d; j++) {
            inner_item = PyList_GetItem(temp_point, j);
            vector[j] = PyFloat_AsDouble(inner_item);
        }
        points[i] = vector;
    }

    return points;
}
static PyObject *create2DPyObject(double** matrix, int n, int d) {
    int i,j;
    PyObject *currentVector, *pyMatrix,*num;
    pyMatrix = PyList_New(n);
    verifyNotNULL(pyMatrix)
    for (i = 0; i < n; i++) {
        currentVector = PyList_New(d);
        verifyNotNULL(currentVector)
        for (j = 0; j < d; j++) {
            num = PyFloat_FromDouble(matrix[i][j]);
            verifyNotNULL(num)
            PyList_SET_ITEM(currentVector, j, num);
        }
    PyList_SET_ITEM(pyMatrix, i, currentVector);
    }
    return pyMatrix;
}



static PyMethodDef methods[] = {
    {"spk",(PyCFunction)kmeans_api, METH_VARARGS, PyDoc_STR("takes 2 python lists, max iteration value, Convergence value")},
    {"wam",(PyCFunction)wam_api,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the weight adjacency matrix")},
    {"ddg",(PyCFunction)ddg_api,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the diagonal degree matrix")},
    {"gl",(PyCFunction)gl_api,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the graph laplacian matrix")},
    {"jacobi",(PyCFunction)jacobi_api,METH_VARARGS,PyDoc_STR("jacobis on a symmetric matrix,second argument should be \"sorted\" for spk purposes\n, returns(values,vectors matrix,k)")},
      /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     

};


static struct PyModuleDef mykmeanssp = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    methods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject* m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}