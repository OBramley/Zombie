from distutils.core import setup, Extension

module= Extension("myModule", sources=["hello.c"])

setup( name='Package', 
    version='1.0', 
    description='This is a package module',
    ext_modules=[module]
)