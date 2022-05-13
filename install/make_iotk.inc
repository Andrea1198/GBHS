# install/make_iotk.inc.  Generated from make_iotk.inc.in by configure.

# -----------------------------
# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

# -----------------------------
# compiler definitions


CC                  = cc
CFLAGS              = -O3 $(DFLAGS) $(IFLAGS)
CPP                 = cpp
CPPFLAGS            = -P -traditional $(DFLAGS) $(IFLAGS)
F90                 = gfortran
MPIF90              = mpif90
F90FLAGS            = $(FFLAGS) -x f95-cpp-input $(FDFLAGS) $(IFLAGS) 
F77                 = f77
FFLAGS              = -O3 -g
FFLAGS_NOOPT        = -O0 -g
LD                  = mpif90
LDFLAGS             = -g -pthread
AR                  = ar
ARFLAGS             = ruv
RANLIB              = ranlib


#------------------------------
# specific definitions

DFLAGS         =  -D__FFTW -D__MPI
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)
IFLAGS         = $(AUTO_INCLUDE)
LIBOBJS        = $(AUTO_LIBS)


# -----------------------------
# libraries

TOPDIR              = /home/apintus/GBHS
#
LIB_MPI             = 
LIB_MASS            = 
#
LIBS                = $(LIB_MPI) $(LIB_MASS)

