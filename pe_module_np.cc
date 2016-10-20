#include <Python.h>
#include <numpy/arrayobject.h>
#include <cstring>

#include "permutation_entropy.cc"

static PyObject *
permutation_entropy_wrapper(PyObject *self, PyObject *args)
{
    
    int pe_n;
    //int ts_ndims;
    double pe_res;

    long int ts_len;

    char * pe_search_method;
    long int * ts_dims;
    double * ts_data_ptr;

    // default arg values
    PyObject * np_in_raw_pyobj = NULL;
    PyObject * np_ts_obj = NULL;
    PyObject * return_value_obj = NULL;
    pe_search_method = NULL;  

    
    // args: PyArrayObject (numpy array) ts
    //       int pe_n
    // np_in_raw_pyobj gets a borrowed reference (no IN/DECREFs needed)
    if (!PyArg_ParseTuple(args, "Oi|z", &np_in_raw_pyobj, 
        &pe_n, &pe_search_method)) {
        goto fail;
    }
    
    // parse to a C-contiguous order and aligned numpy array
    np_ts_obj = PyArray_FROM_OTF(np_in_raw_pyobj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (np_ts_obj == NULL) {
        goto fail;
    }
    
    // get the input array length
    //ts_ndims = PyArray_NDIM(np_ts_obj);
    ts_dims  = PyArray_DIMS(np_ts_obj);

    ts_len = ts_dims[0]; // always use first dimension TODO AXIS OPTION?
    
    ts_data_ptr = (double *) PyArray_DATA(np_ts_obj);
    if (ts_data_ptr == NULL) {
        goto fail;
    }

    // ACTUAL c++ PE FUNCTIONS
    if (pe_search_method == NULL) {
        // TODO AUTO METHOD by default
        pe_res = 666.0;
    } else if (!strcmp(pe_search_method, "full_array_stat")) {
    
        pe_res = permutation_entropy::permutation_entropy_array_stats(ts_data_ptr, (int) ts_len, pe_n);
        //cout << "FULL";        
    } else if (!strcmp(pe_search_method, "dict_array_stat")) {
    
        pe_res = permutation_entropy::permutation_entropy_dictionary_stats(ts_data_ptr, (int) ts_len, pe_n);
        //cout  << "DICTIONARY";
    } else {
        // TODO AUTO METHOD by default
        pe_res = 666.0;
    }

    // build retrun obj (new, reference already INCREF'd)
    return_value_obj = Py_BuildValue("d", pe_res);
    if (return_value_obj == NULL) {
        goto fail;
    }

    // decref auxiliary objs
    Py_DECREF(np_ts_obj);

    return return_value_obj;

    fail:
        // DECREF everything except NULLs
        Py_XDECREF(np_ts_obj);
        Py_XDECREF(return_value_obj);
        return NULL;
}


// define functions in module
static PyMethodDef pe_module_np_methods[] =
{
     {"permutation_entropy",
      permutation_entropy_wrapper,
      METH_VARARGS,
      "evaluates the permutation entropy a numpy array"},
      {NULL, NULL, 0, NULL}
};


// module initialization
PyMODINIT_FUNC
initpe_module_np(void)
{
     (void) Py_InitModule("pe_module_np", pe_module_np_methods);
     import_array();
}


















