FC = /opt/intel/Compiler/11.1/064/bin/intel64/ifort 
FFLAGS = -fPIC

FFT_INC = -I/home/ksmith/opt/include/
FFT_LIBS = -L/home/ksmith/opt/lib/ -lfftw3f

FWRAP = ~/Devel/fwrap/ve-python25/bin/python ~/Devel/fwrap/fwrap-git/fwrapc.py

FWNAME = fwderiv

all:
	$(FC) $(FFLAGS) $(FFT_INC) -c deriv_fft.f90
	$(FC) $(FFLAGS) $(FFT_INC) -c deriv_wrap.f90
	$(FWRAP) deriv_wrap.f90 --build --name=$(FWNAME) --objects deriv_fft.o --override  --fcompiler=intelem $(FFT_LIBS)
	touch $(FWNAME)/__init__.py

clean:
	-rm -r *.o *.mod
	-rm -rf fwderiv
