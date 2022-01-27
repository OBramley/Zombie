// #define PY_SSIZE_T_CLEAN
#include <Python/Python.h>
#include <string.h>

static PyObject * helloworld(PyObject *self, PyObject *args){
    char *name;
    char greeting[255]="Hello";
    if(!PyArg_ParseTuple(args,"s",&name)){
        return NULL;
    }
    strcat(greeting,name);
    return Py_BuildValue("s", greeting);
}

static char greeting_doc[]= 
    "name( ): Greets with your name\n";

static PyMethodDef greetModule[] = {
    {"helloworld",(PyCFunction)helloworld, METH_VARARGS, greeting_doc}, {NULL}
};

void initname(void){
    Py_InitModule3("helloworld", greetModule, "Extension module example!");
}
