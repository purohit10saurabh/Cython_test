import sys
import os
from setuptools import setup, Extension

setup(
    # Information
    name = "testcython",
    version = "1.0.0",
    # Build instructions
    ext_modules = [Extension("testcython", ["testcython.pyx", os.path.join( "..","..","c++","test.cpp" ) ], include_dirs=[ os.path.join( "..", "..", "c++" ) ], language="c++", extra_compile_args=['-std=c++11','-O3'])]
)
