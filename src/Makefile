# Simple Makefile for HV-DFA
# Build targets:
#   make hv_orig       - original FORTRAN program 
#   make python        - Python wrapper (replaces hv_api.exe)
#   make all           - both
#   make clean         - remove build artifacts

FC=gfortran
OPT ?= -O2
USE_OPENMP ?= 1

MODDIR = hvswdpy/_fortran
BINDIR = ../bin

ifeq ($(USE_OPENMP),1)
  OMPFLAG = -fopenmp
else
  OMPFLAG =
endif

FFLAGS ?= $(OPT) $(OMPFLAG) -J$(MODDIR)
LDFLAGS ?= $(OMPFLAG)

# Core sources (shared by API build)
SRCS_COMMON = \
  modules.f90 \
  aux_procedures.f90 \
  Dispersion.f90 \
  root_solver.f90 \
  dispersion_solver.f90 \
  dispersion_equation.f90 \
  greens_rayleigh.f90 \
  greens_love.f90 \
  rayleigh_mode.f90 \
  love_mode.f90 \
  bodywave_integrals.f90

OBJS_COMMON = $(SRCS_COMMON:.f90=.o)

# Python wrapper configuration
PYTHON ?= python3
PYTHON_SETUP = setup.py

.PHONY: all clean python

# Build both targets
all: hv_orig python

# Original FORTRAN CLI: sources are under src/cli/, includes modules from hvswdpy/_fortran
hv_orig: cli/HV.f90
	$(FC) $(FFLAGS) -ffree-form -I$(MODDIR) $< -o $(BINDIR)/hv_orig $(LDFLAGS)

# Python wrapper: build the hvdfa extension using setup.py here
python:
	$(PYTHON) setup.py build_ext --inplace

# Python wrapper:  install into the active environment
python-dev:
	$(PYTHON) -m pip install -e .

# Generic compile rule for manual Fortran compilation
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod $(BINDIR)/hv_orig
	rm -rf build/ ../build/
	rm -f HVSWDpy*.so HVSWDpy-f2pywrappers2.f90 HVSWDpymodule.c
	rm -f HVSWDpy*.cpython-*.so

