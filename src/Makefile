# Simple Makefile for HV-DFA
# Build targets:
#   make hv_orig       - original FORTRAN program 
#   make python        - Python wrapper (replaces hv_api.exe)
#   make all           - both
#   make clean         - remove build artifacts

FC=gfortran
OPT ?= -O2
USE_OPENMP ?= 1
MODDIR = src/hvswdpy/_fortran

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

# Original FORTRAN CLI: sources moved under src/cli/, includes modules from src/hvswdpy/_fortran
hv_orig: src/cli/HV.f90
	$(FC) $(FFLAGS) -ffree-form -I$(MODDIR) $< -o bin/hv_orig $(LDFLAGS)

# Python wrapper: build the hvdfa extension
python: $(PYTHON_SETUP)
	$(PYTHON) $(PYTHON_SETUP) build_ext --inplace

# Generic compile rule (for manual Fortran compilation if needed)
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod bin/hv_orig
	rm -rf build/
	rm -f HVSWDpy*.so HVSWDpy-f2pywrappers2.f90 HVSWDpymodule.c
	rm -f HVSWDpy*.cpython-*.so
	rm -r __pycache__

