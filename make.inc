# make.inc.  Generated from make.inc.in by configure.

# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

# most fortran compilers can directly preprocess c-like directives: use
# 	$(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
# 	$(CPP) $(CPPFLAGS) $< -o $*.F90
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

# .f.o and .c.o: do not modify

.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS)  -c $<



# Top QE directory, useful for locating libraries,  linking QE with plugins
# The following syntax should always point to TOPDIR:
TOPDIR = $(dir $(abspath $(filter %make.inc,$(MAKEFILE_LIST))))
# if it doesn't work, uncomment the following line (edit if needed):

# TOPDIR = /e/coding/GBHS

# DFLAGS  = precompilation options (possible arguments to -D and -U)
#           used by the C compiler and preprocessor
# FDFLAGS = as DFLAGS, for the f90 compiler
# See include/defs.h.README for a list of options and their meaning
# With the exception of IBM xlf, FDFLAGS = $(DFLAGS)
# For IBM xlf, FDFLAGS is the same as DFLAGS with separating commas

# MANUAL_DFLAGS  = additional precompilation option(s), if desired
#                  BEWARE: it does not work for IBM xlf! Manually edit FDFLAGS
MANUAL_DFLAGS  =
DFLAGS         =  -D__FFTW
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
#          typically -I../include -I/some/other/directory/
#          the latter contains .e.g. files needed by FFT libraries

IFLAGS         = -I$(TOPDIR)/include -I../../include/ $(MOD_FLAG)../baselib $(IOTK_INCS) $(PYTHON_INCS)

# MOD_FLAGS = flag used by f90 compiler to locate modules
# Each Makefile defines the list of needed modules in MODFLAGS

MOD_FLAG      = -I

# Compilers: fortran-90, fortran-77, C
# If a parallel compilation is desired, MPIF90 should be a fortran-90
# compiler that produces executables for parallel execution using MPI
# (such as for instance mpif90, mpf90, mpxlf90,...);
# otherwise, an ordinary fortran-90 compiler (f90, g95, xlf90, ifort,...)
# If you have a parallel machine but no suitable candidate for MPIF90,
# try to specify the directory containing "mpif.h" in IFLAGS
# and to specify the location of MPI libraries in MPI_LIBS

MPIF90         = gfortran
F90           = gfortran
CC             = gcc
F77            = gfortran

# C preprocessor and preprocessing flags - for explicit preprocessing,
# if needed (see the compilation rules above)
# preprocessing flags must include DFLAGS and IFLAGS

CPP            = cpp
CPPFLAGS       = -P -traditional $(DFLAGS) $(IFLAGS)

# compiler flags: C, F90, F77
# C flags must include DFLAGS and IFLAGS
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax

CFLAGS         = -O3 $(DFLAGS) $(IFLAGS)
F90FLAGS       = $(FFLAGS) -x f95-cpp-input $(FDFLAGS) $(IFLAGS) $(MODFLAGS)
FFLAGS         = -O3 -g

# compiler flags without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack

FFLAGS_NOOPT   = -O0 -g

# compiler flag needed by some compilers when the main program is not fortran
# Currently used for Yambo

FFLAGS_NOMAIN   = 

# Linker, linker-specific flags (if any)
# Typically LD coincides with F90 or MPIF90, LD_LIBS is empty

LD             = gfortran
LDFLAGS        = -g -pthread
LD_LIBS        = 

# External Libraries (if any) : blas, lapack, fft, MPI

# If you have nothing better, use the local copy via "--with-netlib" :
# BLAS_LIBS = /your/path/to/espresso/LAPACK/blas.a
# BLAS_LIBS_SWITCH = internal

BLAS_LIBS      = $(TOPDIR)/extlibs/blas/lib/libblas.a 
BLAS_LIBS_SWITCH = internal

# If you have nothing better, use the local copy via "--with-netlib" :
# LAPACK_LIBS = /your/path/to/espresso/LAPACK/lapack.a
# LAPACK_LIBS_SWITCH = internal
# For IBM machines with essl (-D__ESSL): load essl BEFORE lapack !
# remember that LAPACK_LIBS precedes BLAS_LIBS in loading order

LAPACK_LIBS    = $(TOPDIR)/extlibs/lapack/lib/liblapack.a
LAPACK_LIBS_SWITCH = internal

SCALAPACK_LIBS = 

IOTK_LIBS      = -L../../extlibs/iotk/lib -liotk
IOTK_LIBS_SWITCH = internal
IOTK_INCS      = $(MOD_FLAG)../../extlibs/iotk/include

# nothing needed here if the the internal copy of FFTW is compiled
# (needs -D__FFTW in DFLAGS)

FFT_LIBS       = 

# HDF5
HDF5_LIB = 
#FOX_LIB  = -L$(TOPDIR)/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys 

# For parallel execution, the correct path to MPI libraries must
# be specified in MPI_LIBS (except for IBM if you use mpxlf)

MPI_LIBS       = 

# IBM-specific: MASS libraries, if available and if -D__MASS is defined in FDFLAGS

MASS_LIBS      = 

# python

PYTHON_LIBS         = 
PYTHON_INCS         = 


# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv

AR             = ar
ARFLAGS        = ruv

# ranlib command. If ranlib is not needed (it isn't in most cases) use
# RANLIB = echo

RANLIB         = ranlib

# all internal and external libraries - do not modify

FLIB_TARGETS   = all

LIBOBJS        = $(TOPDIR)/src/ctools/libagwx_ctools.a  $(TOPDIR)/extlibs/iotk/lib/libiotk.a
LIBS           = $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FOX_LIB) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(HDF5_LIB) $(PYTHON_LIBS) $(LD_LIBS)

# wget or curl - useful to download from network
WGET = curl -o

# Install directory - not currently used
PREFIX = /usr/local
