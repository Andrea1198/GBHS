
#if defined __PYTHON
#  include <Python.h>
#endif

#include "c_defs.h"

const char *pycode =
#include "pole_fitting.str"
   " ";

   
int F77_FUNC(pole_fitting_pywrap,POLE_FITTING_PYWRAP)(int *nx, double *xdata, double *ydata, 
             int *npoles_fit, double *params)
{
    /* This function needs to be called after Py_Initilize() and before Py_Finalize() */
#if defined __PYTHON
    PyObject *pX, *pY, *pRes, *pValue;
    PyObject *main_module, *main_dict;
    PyObject *np_obj, *p0_obj;
    int i, n, ndim;

    /* Setup the __main__ module for us to use */
    main_module = PyImport_ImportModule("__main__");
    main_dict   = PyModule_GetDict(main_module);

    n=*nx;
    pX = PyList_New(n);
    for (i=0;i<n;i++) {
      pValue = PyFloat_FromDouble(xdata[i]);
      PyList_SetItem(pX,(Py_ssize_t)i,pValue);
    }
    pY = PyList_New(n);
    for (i=0;i<n;i++) {
      pValue = PyFloat_FromDouble(ydata[i]);
      PyList_SetItem(pY,(Py_ssize_t)i,pValue);
    }
    ndim=3*(*npoles_fit)+1;
    p0_obj = PyList_New(ndim);
    for (i=0;i<ndim;i++) {
      pValue = PyFloat_FromDouble(params[i]);
      PyList_SetItem(p0_obj,(Py_ssize_t)i,pValue);
    }
    //
    i=*npoles_fit;
    np_obj = PyLong_FromSsize_t((Py_ssize_t)i);

    /* Inject variables into __main__ */
    PyDict_SetItemString(main_dict, "xdata", pX);
    PyDict_SetItemString(main_dict, "ydata", pY);
    PyDict_SetItemString(main_dict, "npol", np_obj);
    PyDict_SetItemString(main_dict, "p0", p0_obj);

    /* Run the code snippet above in the current environment */
    PyRun_SimpleString(pycode);

    /* Extract the resultant variable, k */
    pRes = PyMapping_GetItemString(main_dict, "params");

    // get values from Py
    for (i=0;i<ndim+1;i++) {
      params[i] = PyFloat_AsDouble( PyList_GetItem(pRes,(Py_ssize_t)i) );
    }

    /* cleanup */
    Py_XDECREF(pValue);
    Py_XDECREF(pX);
    Py_XDECREF(pY);
    Py_XDECREF(np_obj);
    Py_XDECREF(pRes);

    return 0;
#else
    return 1;
#endif
}

