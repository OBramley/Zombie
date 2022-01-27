from setuptools import setup, Extension

setup( name='greetmodule', version='1.0', \
    ext_modules=[Extension('greet',sources=['hello.c'])]
)