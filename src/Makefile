#intel or fortran compiler
FC=ifort
# FC=gfortran

# GFLAGS=-Wall -Wextra -pedantic -std=f2008 -g -O2 -fno-second-underscore -Wno-compare-reals -fcheck=all -fimplicit-none -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function -fcheck=bounds -I/usr/local/include -ffree-line-length-512

# GFLAGS= -O2 -fast -I/usr/local/include #-ffree-line-length-512 #-fbounds-check   ##INTEL
GFLAGS= -O2 -I/usr/local/include -ffree-line-length-512 #-fbounds-check            ##FORTRAN
CFLAGS= -lm -g -O2 -lcfitsio

OUT=CAIRO

all: blas.o lbfgsb.o linpack.o timer.o mod_constants.o mod_convert.o mod_fits.o mod_inout.o mod_array.o mod_model.o mod_optimize.o mod_minimize.o mod_functions.o mod_start.o mod_rohsa.o main.o
	$(FC) $^ $(CFLAGS) -o $(OUT)

main.o: main.f90
	$(FC) $(GFLAGS) -c $^

# L-BFGS-B package
blas.o: blas.f
	$(FC) $(GFLAGS) -c $^
lbfgsb.o: lbfgsb.f
	$(FC) $(GFLAGS) -c $^
linpack.o: linpack.f
	$(FC) $(GFLAGS) -c $^
timer.o: timer.f
	$(FC) $(GFLAGS) -c $^

# My files
mod_constants.o: mod_constants.f90
	$(FC) $(GFLAGS) -c $^

mod_convert.o: mod_convert.f90
	$(FC) $(GFLAGS) -c $^

mod_fits.o: mod_fits.f90
	$(FC) $(GFLAGS) -c $^

mod_inout.o: mod_inout.f90
	$(FC) $(GFLAGS) -c $^

mod_array.o: mod_array.f90
	$(FC) $(GFLAGS) -c $^

mod_model.o: mod_model.f90
	$(FC) $(GFLAGS) -c $^

mod_optimize.o: mod_optimize.f90
	$(FC) $(GFLAGS) -c $^

mod_minimize.o: mod_minimize.f90
	$(FC) $(GFLAGS) -c $^

mod_functions.o: mod_functions.f90
	$(FC) $(GFLAGS) -c $^

mod_start.o: mod_start.f90
	$(FC) $(GFLAGS) -c $^

mod_rohsa.o: mod_rohsa.f90
	$(FC) $(GFLAGS) -c $^

clean:
	\rm -rf *.o $(OUT) *.mod *~

run: all
	./$(OUT)

rerun: clean all run
