# define the LIBs
LIB = -L/usr/local/Cellar/lapack/3.11 -llapack #one should adjust this path to link your local lapack
LIB1 = -L./adyo_v1_0 -l:libcwf_cpp.a  

# define the compilers
FC = gfortran
F90 = gfortran

BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)

# define the compiling options
FFLAGS = -O2 -Wtabs -fopenmp -Wconversion -ffixed-line-length-0 -cpp -fPIC\
-DBASE=\"$(BASE)\" -DVERDATE='"$(VERDATE)"' -DVERREV='"$(VERREV)"' -DCOMPDATE='"$(COMPDATE)"' 

# define the object files
OBJF = precision.o constants.o system.o mesh.o pot_class.o \
clebsch.o spharm.o gauss_mesh.o channels.o \
gauss.o coulcc.o coul90.o\
generate_laguerre.o rot_potential.o \
solve_eigen.o input.o \
matrix_element.o bound.o npcc.o yamaguchi.o scatt.o  

# define the static library
STLIB = liballmodules.a

.SUFFIXES: .F90 .f90 .f95 .F 

all: subdir $(STLIB) COLOSS

subdir:
	make -C ./adyo_v1_0

$(STLIB): $(OBJF) 
	$(AR) rv $@ $(OBJF)
	ranlib $@

COLOSS:  $(STLIB) COLOSS.o
	$(FC) $(FFLAGS) -o COLOSS COLOSS.o liballmodules.a $(LIB) $(LIB1) -lstdc++

COLOSS.o: COLOSS.F
	$(FC) -c COLOSS.F $(FFLAGS)


.f90.o:
	$(F90) $(FFLAGS) -c $<

.F90.o:
	$(F90) $(FFLAGS) -c $<

.f95.o:
	$(F90) $(FFLAGS) -c $<

.F.o:
	$(F90) $(FFLAGS) -c $<

.f.o:
	$(F90) $(FFLAGS) -c $<

generate_laguerre.o: generate_laguerre.f
	$(F90) $(FFLAGS) -c $< -I./adyo_v1_0 

.PHONY: clean 
clean:
	rm -f COLOSS $(STLIB) *.mod core *.o *.a
