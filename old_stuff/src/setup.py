from distutils.core import setup, Extension 

module1 = Extension('c_ham', sources=['ham.c'])

setup(name='ham_c',version='1.0', description='Hope',ext_modules=[module1])