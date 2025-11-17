# define the LIBs
LIB = -L/opt/homebrew/opt/lapack/lib -llapack -lblas
LIB1 = ./adyo_v1_0/libcwf_cpp.a

# define the compilers
FC = gfortran
F90 = gfortran

# define source directory
SRCDIR = src
BUILDDIR = build
VPATH = $(SRCDIR)

BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)

# define the compiling options
FFLAGS = -g -Wtabs -fopenmp -Wconversion -ffixed-line-length-0 -cpp -fPIC\
-DBASE=\"$(BASE)\" -DVERDATE='"$(VERDATE)"' -DVERREV='"$(VERREV)"' -DCOMPDATE='"$(COMPDATE)"' -J$(BUILDDIR)

# define the object files
OBJF = $(BUILDDIR)/precision.o $(BUILDDIR)/constants.o $(BUILDDIR)/system.o $(BUILDDIR)/mesh.o $(BUILDDIR)/pot_class.o \
$(BUILDDIR)/clebsch.o $(BUILDDIR)/spharm.o $(BUILDDIR)/gauss_mesh.o $(BUILDDIR)/channels.o $(BUILDDIR)/readpot.o\
$(BUILDDIR)/gauss.o $(BUILDDIR)/coulcc.o $(BUILDDIR)/coul90.o $(BUILDDIR)/interpolation.o\
$(BUILDDIR)/generate_laguerre.o $(BUILDDIR)/rot_potential.o \
$(BUILDDIR)/solve_eigen.o $(BUILDDIR)/input.o \
$(BUILDDIR)/matrix_element.o $(BUILDDIR)/bound.o $(BUILDDIR)/npcc.o $(BUILDDIR)/yamaguchi.o $(BUILDDIR)/scatt.o

# define the static library
STLIB = $(BUILDDIR)/liballmodules.a

.SUFFIXES: .F90 .f90 .f95 .F

all: $(BUILDDIR) subdir $(STLIB) COLOSS

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

subdir:
	make -C ./adyo_v1_0

$(STLIB): $(OBJF)
	$(AR) rv $@ $(OBJF)
	ranlib $@

COLOSS: $(STLIB) $(BUILDDIR)/COLOSS.o
	$(FC) $(FFLAGS) -o COLOSS $(BUILDDIR)/COLOSS.o $(STLIB) $(LIB) $(LIB1) -lc++

$(BUILDDIR)/COLOSS.o: $(SRCDIR)/COLOSS.F | $(BUILDDIR)
	$(FC) -c $< $(FFLAGS) -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.F90 | $(BUILDDIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.f95 | $(BUILDDIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.F | $(BUILDDIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.f | $(BUILDDIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/generate_laguerre.o: $(SRCDIR)/generate_laguerre.f | $(BUILDDIR)
	$(F90) $(FFLAGS) -c $< -o $@ -I./adyo_v1_0 

.PHONY: clean COLOSS
clean:
	rm -f COLOSS
	rm -rf $(BUILDDIR)
	cd adyo_v1_0 && $(MAKE) clean
