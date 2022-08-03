
#if defined __PYTHON
#  include <Python.h>
#endif

#include "c_defs.h"

int F77_FUNC(python_initialize,PYTHON_INITIALIZE)()
{
#if defined __PYTHON
   Py_Initialize();
#endif
   return 0;
}

int F77_FUNC(python_finalize,PYTHON_FINALIZE)()
{
#if defined __PYTHON
   Py_Finalize();
#endif
   return 0;
}
