#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <assert.h>
#include <math.h>

double **T;
int K;
int d;
int N;

static double** c_means(double** DataPts,double** centroids);

/*
 * function that implemented kmeans++ algorithm from HW2
 * the function take an input the DataPts and the initial centroids and return the final centroids
 *
 * in this function we use the function that we declare in the spkmeans.h file and implemented in
 * the spkmeans.c file.
 *
 */
static double** c_means(double** DataPts,double** centroids) {
    int i;
    int t;


    double **copy_clusters;
    double **group_clusters;


    copy_clusters = (double**) malloc(K*sizeof(double*));
    assert1(copy_clusters);
    assert(copy_clusters!=NULL);
    for (i = 0; i< K; i++){
        copy_clusters[i] = (double*) malloc(d*sizeof(double));
        assert2(copy_clusters[i]);
        assert(copy_clusters[i]!=NULL);
    }

    group_clusters = (double**) malloc(K*sizeof(double*));
    assert1(group_clusters);
    assert(group_clusters!=NULL);
    for (i = 0; i< K; i++){
        group_clusters[i] = (double*) malloc((d+1)*sizeof(double));
        assert2(group_clusters[i]);
        assert(group_clusters[i]!=NULL);
    }

    i = 0;
    t = 1;
    while ((i < 300)&&(t == 1)){
        copy_clus(centroids, copy_clusters);
        calc_new_clus(centroids, DataPts, group_clusters);
        new_clus(group_clusters, centroids);
        t = update(centroids, copy_clusters);
        i++;
    }


    for(i = 0 ; i < N  ; i++){
        free(DataPts[i]);
    }
    free(DataPts);

    for(i = 0 ; i < K ; i++){
        free(copy_clusters[i]);
        free(group_clusters[i]);
    }

    free(group_clusters);
    free(copy_clusters);

    return centroids;

}

/*
 * API functions
 * this function is the function that connects between python and C, this is also a wrapping function for cmeans
 * which is the function above (where we implemented k means)
 * in this function we use the function that we declare in the spkmeans.h file and implemented in
 * the spkmeans.c file that implementes the first section of the algorithm of spkmeans
 *
 */
static PyObject* spk_means_capi(PyObject *self, PyObject *args){
    double **DataPts;
    char *goal;

    Py_ssize_t n_py, d_py, i, j;
    PyObject  *_DataPts, *point, *coordinate;
    PyObject *py_T, *py_T_LINE;

    if(!PyArg_ParseTuple(args, "iOs", &K ,&_DataPts, &goal)) {
        return NULL;
    }

    if (!PyList_Check(_DataPts))
        return NULL;

    N = (long) PyObject_Length(_DataPts);
    d = (long) PyObject_Length(PyList_GetItem(_DataPts,0));

    n_py = PyList_Size(_DataPts);
    d_py = PyList_Size(PyList_GetItem(_DataPts,0));

    DataPts= malloc(sizeof(double *) * N);
    assert1(DataPts);
    assert(DataPts != NULL);
    for (i = 0; i < N; ++i) {
        DataPts[i] = malloc(sizeof(double) * d);
        assert2(DataPts[i]);
        assert(DataPts[i] != NULL);
    }

    for (i = 0; i < n_py; i++) {
        point = PyList_GetItem(_DataPts, i);
        for(j = 0; j < d_py; j++){
            coordinate = PyList_GetItem(point, j);
            DataPts[i][j] = PyFloat_AsDouble(coordinate);
        }
    }

    if(!strcmp(goal,"spk")){
        spk_case(DataPts);
    }
    else if(!strcmp(goal,"jacobi")){
        jacobi_case(DataPts);
        Py_RETURN_NONE;
    }
    else if(!strcmp(goal,"wam")){
        wam_case(DataPts);
        Py_RETURN_NONE;
    } else if(!strcmp(goal,"ddg")){
        ddg_case(DataPts);
        Py_RETURN_NONE;
    } else if(!strcmp(goal,"lnorm")){
        lnorm_case(DataPts);
        Py_RETURN_NONE;
    } else {
        printf("Invalid Input!\n");
        exit(0);
    }

    py_T = PyList_New(N);
    if (!PyList_Check(py_T))
        return NULL;

    for (i=0; i<N; i++){
        py_T_LINE = PyList_New(d);
        if (!PyList_Check(py_T_LINE))
            return NULL;
        for (j=0; j<d; j++){
            PyList_SET_ITEM(py_T_LINE, j, Py_BuildValue("d",T[i][j]));
        }
        PyList_SetItem(py_T, i, Py_BuildValue("O", py_T_LINE));
    }

    return Py_BuildValue("O", py_T);

}

/*
 * API functions
 * this function is the function that connects between python and C, this is also a wrapping function for cmeans
 * which is the function above (where we implemented k means)
 * in this function we use the function that we declare in the spkmeans.h file and implemented in
 * the spkmeans.c file that implementes the second section of the algorithm of spkmeans that is to find the
 * K clusters.
 *
 */
static PyObject* c_means_capi(PyObject *self, PyObject *args){
    double **DataPts, **centroids;

    PyObject  *_DataPts, *_centroids, *point, *centroid, *coordinate;
    Py_ssize_t n, k, d_py, i, j;


    if(!PyArg_ParseTuple(args, "OO" ,&_DataPts, &_centroids)) {
        return NULL;
    }

    if (!PyList_Check(_DataPts) || !PyList_Check(_centroids))
        return NULL;


    N = (long) PyObject_Length(_DataPts);
    K = (long) PyObject_Length(_centroids);
    d = (long) PyObject_Length(PyList_GetItem(_DataPts,0));


    n = PyList_Size(_DataPts);
    k = PyList_Size(_centroids);
    d_py = PyList_Size(PyList_GetItem(_DataPts,0));

    DataPts= malloc(sizeof(double *) * N);
    assert1(DataPts);
    assert(DataPts != NULL);
    for (i = 0; i < N; ++i) {
        DataPts[i] = malloc(sizeof(double) * d);
        assert2(DataPts[i]);
        assert(DataPts[i] != NULL);
    }
    centroids = malloc(sizeof(double *) * K);
    assert1(centroids);
    assert(centroids != NULL);

    for (i = 0; i < K; ++i) {
        centroids[i] = malloc(sizeof(double) * d);
        assert2(centroids[i]);
        assert(centroids[i] != NULL);
    }

    for (i = 0; i < k; i++) {
        centroid = PyList_GetItem(_centroids, i);
        for(j = 0; j < d_py; j++){
            coordinate = PyList_GetItem(centroid, j);
            centroids[i][j] = PyFloat_AsDouble(coordinate);
        }
    }

    for (i = 0; i < n; i++) {
        point = PyList_GetItem(_DataPts, i);
        for(j = 0; j < d_py; j++){
            coordinate = PyList_GetItem(point, j);
            DataPts[i][j] = PyFloat_AsDouble(coordinate);
        }
    }

    centroids = c_means(DataPts,centroids);

    for(i = 0; i < K; i++){
        for(j = 0; j < d; j++){
            if (centroids[i][j] < 0.0 && centroids[i][j] > -0.00005){
                printf("%.4f",0.0000);
            } else {
                printf("%.4f",centroids[i][j]);
            }
            if(j < K-1){
                printf("%s",",");
            }
        }
        if(i < K-1){
            printf("\n");
        }
    }printf("\n");


    for(i = 0 ; i < K ; i++){
        free(centroids[i]);
    }
    free(centroids);

    Py_RETURN_NONE;
}

/* defining the name of the function c_means_capi to be kmeans, in order to use it python with that name */
static PyMethodDef _methods[] = {
        {"fit1", (PyCFunction)spk_means_capi, METH_VARARGS, PyDoc_STR("Please enter the fields : K and arrays of : Observations and Goal")},
        {"fit2", (PyCFunction)c_means_capi, METH_VARARGS, PyDoc_STR("Please enter the fields : MAX_ITER and arrays of : Observations and Centroids")},
        {NULL, NULL, 0, NULL}
};

/* defining the name of the module */
static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
        NULL,
        -1,
        _methods
};

PyMODINIT_FUNC
/* initializing the CAPI module */
PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


