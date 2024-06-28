# compiler
FC = gfortran
CC = g++

# compile flags
FCFLAGS = -c -fpic -ffree-line-length-none
CCFLAGS = -c -fpic -std=c++11

# source files
F90_SRC = cwf_cpp.f90
CPP_SRC = coulomb.cpp

# object files
F90_OBJ = $(F90_SRC:.f90=.o)
CPP_OBJ = $(CPP_SRC:.cpp=.o)

# output library
LIB = libcwf_cpp.a

all: $(LIB)

$(LIB): $(F90_OBJ) $(CPP_OBJ)
	ar rcs $@ $^
	ranlib $@

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

%.o: %.cpp
	$(CC) $(CCFLAGS) -o $@ $<

clean:
	rm -f $(F90_OBJ) $(CPP_OBJ) $(LIB)
