#ifdef __cplusplus
extern "C"
#endif

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "okada_dc3d.h"
#include "okada_disloc3d.h"


static PyObject *okada_rect(PyObject *self, PyObject *args) {
    // Initialize NumPy array object
    PyArrayObject *obs;          // ndarray of [n x 3], xyz coordinates of observation stations
    PyArrayObject *models;       // ndarray of [n x 10]  
    double mu; 		     // shear modulus
    double nu; 		     // Poisson's ratio
    
    // Parse arguments from Pyhton 
    if (!PyArg_ParseTuple(args, "O!O!dd", &PyArray_Type, &obs, &PyArray_Type, &models, &mu, &nu)) {
    	PyErr_SetString(PyExc_TypeError, "Expected NumPy arrays and Python floats as input.");
    	return NULL;
    }
    /* if (!PyArray_Check(obs)) { return NULL;} */
    
    // obs and models must be a 2-D ndarray array, 1-D ndarray list is not allowed
    if ((PyArray_NDIM(obs) != 2) || (PyArray_NDIM(models) != 2)) {
        PyErr_SetString(PyExc_ValueError, "The observations and fault models must be 2-D ndarrays of [n x 3] and [n x 10]! If it is a 1-D ndarray list (usually 1 stations or 1 fault model), reshape it to 2-D!\n");
        return NULL;
    }
    
    // Check if obs has 3 columns and models has 10 columns
    if ((PyArray_SIZE(models) % 10 != 0) || PyArray_SIZE(obs) % 3 != 0) {
        PyErr_SetString(PyExc_ValueError, "The observations should be an array of [n x 3] and models should be [n x 10]!\n");
        return NULL;
    }

    // Get the numbers of models and stations
    npy_intp nmodels = PyArray_SIZE(models) / 10;
    npy_intp nobs    = PyArray_SIZE(obs) / 3;
    /*printf("nmodels:  %ld\n", nmodels); 
     printf("nobs:  %ld\n", nobs); */
    
    // Initialize C array (2-D)
    double **c_models = NULL;
    double **c_obs    = NULL;

    // convert to 2-D array if models or obs is a 1-D list
    //if (nmodels == 1) {c_models = PyArray_Newshaple}
    PyArray_AsCArray((PyObject **)&models, &c_models, (npy_intp []){nmodels, 10}, 2, PyArray_DescrFromType(NPY_DOUBLE));
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to convert NumPy array of Fault Models to C array.");
        return NULL;
    }

    PyArray_AsCArray((PyObject **)&obs, &c_obs, (npy_intp[]){nobs, 3}, 2, PyArray_DescrFromType(NPY_DOUBLE));
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to convert NumPy array of Observations Coordinates to C array.");
        return NULL;
    }

    /* printf("\n");
    for (npy_intp i =0; i<nmodels; i++) {
        for (int j=0; j< 10; j++) {
    	printf("%f  ", c_models[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    for (npy_intp i =0; i<nobs; i++) {
        for (int j=0; j< 3; j++) {
    	printf("%f  ", c_obs[i][j]);
        }
        printf("\n");
    } */

    //double 
    //// call disloc3d.c
    //disloc3d(c_models, nmodels, c_obs, nobs, mu, nu, U, D, S, flags);

    
    // free memory
    PyArray_Free((PyObject *)models, c_models);
    PyArray_Free((PyObject *)obs, c_obs);
    
    // return a Python Object Pointer
    Py_RETURN_NONE;
    //return PyFloat_FromDouble(x + y);
}

// ---------------------------------------------------------------------- 
static PyObject *add(PyObject *self, PyObject *args) {
    double x; 
    double y;
    PyArg_ParseTuple(args, "dd", &x, &y);
    return PyFloat_FromDouble(x + y);
}

// ---------------------------------------------------------------------- 

PyDoc_STRVAR(
    okada_rect_doc,
    "Okada_rect calculates displacement, stress and strain using Okada rectangle dislocation source.\n"
    "\n"
    "okada_rect(observations, models, mu, nu)\n"
    
    "- Arguments:\n"
    "  - observations: An 2-D NumPy array of [nobs x 3], xyz Cartesian coordinates of stations, \n"
    " 		       in which z values should <= 0. \n"
    "  - models	     : An 2-D Numpy array of [nmodles x 10], in which: \n"
    "                [x_uc, y_uc, depth, length, width, strike, dip, str_slip, dip_slip, opening]\n "
    "                 x_uc     : x Cartesian coordinates of the fault reference point (center of fault upper edge), \n"
    "    	      y_uc     : y Cartesian coordinates of the fault reference point (center of fault upper edge), \n"
    "  		      depth    : depth of the fault reference point (center of fault upper edge, always positive vales), \n"
    "  		      length   : length of the fault patches, \n"
    "  		      width    : width of the fault patches, \n"
    " 		      strike   : strike angles of the fault patches, \n"
    "  		      dip      : dip angles of the fault patches, \n"
    "  		      str-slip : strike-slip components of the fault patches, \n"
    "  		      dip-slip : dip-slip components of the fault patches, \n"
    "  		      opening  : tensile components of the fault patches, \n"
    "  - mu 	     : Shear modulus, \n"
    "  - nu 	     : Poisson's ratio, \n"
    "\n"
    "- Output:\n"
    "  - u           : displacements,\n"
    "  - d           : 9 spatial derivatives of the displacements,\n"
    "  - s           : 9 stress tensor components, 6 of them are independent,\n"
    "  - flags       : flags [nobs x nmodels].\n"
    "                  flags = 0 : normal,\n "
    "                  flags = 1: the Z value of the obs > 0,\n"
    "	               flags = 10: the depth of the fault upper center point < 0, \n"
    "	               flags = 100: singular point, observation is on fault edges,\n"
    "                  flags could also be 11, 101, 110, 111, indicating the sum of multiple cases.\n "
    "\n"
    "NOTE: The units of displacements are same to dislocation slip (str-slip ...);\n"
    "      Strains are dimensionless, Cartesian coordinates and dislocation slip are better kept as 'm', so that there is no need to make conversion to the strain; \n"
    "      The units of stress depends on shear modulus mu, Pa is recommended. \n"

);

// 4. 模块函数列表
static PyMethodDef method_funcs[] = {
    // function name, function pointer, argument flag, function docs
    {"add", add, METH_VARARGS, "Add two numbers together."},
    {"okada_rect", okada_rect, METH_VARARGS, okada_rect_doc},
    {NULL, NULL, 0, NULL}
};

// ---------------------------------------------------------------------- 
// 3. module definition
static struct PyModuleDef dislocation = {
    PyModuleDef_HEAD_INIT,
    "dislocation", 			  // module name
    "This is a module named dislocation.", // module docs, could be called by help(module_name)
    -1, 
    method_funcs 			  // module functions list 
};

// ---------------------------------------------------------------------- 
// 1. C extension entrance, PyInit_+module name
PyMODINIT_FUNC PyInit_dislocation() {
    printf("Now the module has been imported!\n");
    // 2. create module, the argument type is PyModuleDef
    // return PyModule_Create(&dislocation);
    
    // Initialize Numpy
    PyObject *module = PyModule_Create(&dislocation);
    import_array();
    return module;
}
