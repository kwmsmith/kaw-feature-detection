#!/bin/bash

FC=$(which ifort) CFLAGS="-m32" LINKFLAGS="-m32" LIB="fftw3f fftw3" LIBPATH="/home/ksmith/opt/lib" INCLUDES="/home/ksmith/opt/include" python fwc.py configure --name=test build -v *.f90
