
FC = gfortran
FCLAGS = -g -O3 -fbounds-check #-ffpe-trap=zero,overflow,underflow 
FFLAGS = -llapack

SOURCE1 = fortran/Modules.f90 fortran/QRBMModule.f90 fortran/QRBMProgram.f90

SOURCE2 = fortran/Modules.f90 fortran/symmQRBMModule.f90 fortran/symmQRBMProgram.f90

SOURCE3 = fortran/Modules.f90 fortran/symmQRBMModule.f90 fortran/QRBMModule.f90 fortran/performanceProgram.f90

all: QRBM symmQRBM performance



symmQRBM: $(SOURCE2)
	$(FC) $(FCLAGS) -o $@ $(SOURCE2) $(FFLAGS)

QRBM: $(SOURCE1)
	$(FC) $(FCLAGS) -o $@ $(SOURCE1) $(FFLAGS) 

performance :  $(SOURCE3)
	$(FC) $(FCLAGS) -o $@ $(SOURCE3) $(FFLAGS) 

clean:
	-rm -f *.mod QRBM symmQRBM performance

