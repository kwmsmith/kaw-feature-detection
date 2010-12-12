#!/bin/bash

unlink *.so

FC=$(which ifort) \
CFLAGS="-m32" \
LINKFLAGS="-m32" \
STLIB="fftw3f fftw3" \
STLIBPATH="/home/ksmith/opt/fftw32/lib/" \
INCLUDES="/home/ksmith/opt/fftw32/include/" \
python fwrapc.py configure --name=fwderiv build *.f90

ln -s fwproj/build/fwderiv.so
