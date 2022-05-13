#
#===============================
# AGWX package
#===============================
#
include ../../make.inc

DIR=baselib
LIBOBJS_LOC=$(LIBOBJS) 

OBJS= \
dawson.o \
d1mach.o\
dh12.o\
dhfti.o\
dlpdp.o\
dlsei.o\
dlsi.o\
dwnlit.o\
dwnlsm.o\
dwnlt1.o\
dwnlt2.o\
dwnlt3.o\
dwnnls.o\
fdump.o\
i1mach.o\
j4save.o\
xercnt.o\
xerhlt.o\
xermsg.o\
xerprn.o\
xersve.o\
xgetua.o\
anglib.o\
bspline_module.o \
cleandir.o \
compute_efermi.o \
constants.o \
crystal_io_base.o \
crystal_io.o \
date_and_tim.o \
efermig.o \
erf.o \
error.o \
etsf_error.o \
expint.o \
fft_scalar.o \
files.o \
fitting_get_peaks.o \
fitting_get_trial.o \
flush_unit.o \
funct.o \
functionals.o \
gauss_module.o \
math_common.o \
mod_LIVE_t.o \
mod_stderr.o \
mod_openmp.o \
mod_wrappers.o \
mod_poly.o \
more_functionals.o \
mp_tools.o \
newton.o \
lsda_functionals.o \
wgauss.o \
input_from_file.o \
io_global.o \
kind.o \
locate.o \
log_module.o \
mp.o \
mp_global.o \
mod_nonneg_leastsq.o \
mod_func4fit.o \
mod_fitting.o \
mod_fitting_base.o \
mod_fitting_genetic.o \
parallel_include.o \
parameters.o \
parser_base.o \
parser.o \
qexml.o \
qexpt.o \
radial_gradients.o \
random_numbers.o \
random_orthovect.o \
shmem_include.o \
shutdown.o \
simpsn.o \
smearing_base.o \
smearing.o \
sort.o \
sph_bes.o \
sph_har.o \
startup.o \
sumkg.o \
timing.o \
util.o \
version.o \
warning.o \
headers.o \
ylmr2.o

#poles_base.o \
#dynam_opr.o \
#dsbevx_omp.o \

#
# main targets
#
exec= libagwx_base.a \
  test_legendre.x \
  test_fitting_simple.x \
  test_fitting_gau.x \
  test_bspline.x \
  test_poly.x

all : $(exec)

libagwx_base.a : $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

test_legendre.x: test_legendre.o $(OBJS) ../ctools/libagwx_ctools.a
	$(LD) $(LDFLAGS) -o $@ test_legendre.o \
	       $(OBJS) $(LIBOBJS_LOC) $(LIBS)
	- ( cd ../../bin ; ln -sf ../src/$(DIR)/$@ . )

test_fitting_simple.x: test_fitting_simple.o $(OBJS) ../ctools/libagwx_ctools.a
	$(LD) $(LDFLAGS) -o $@ test_fitting_simple.o \
	       $(OBJS) $(LIBOBJS_LOC) $(LIBS)
	- ( cd ../../bin ; ln -sf ../src/$(DIR)/$@ . )

test_fitting_gau.x: test_fitting_gau.o $(OBJS) ../ctools/libagwx_ctools.a
	$(LD) $(LDFLAGS) -o $@ test_fitting_gau.o \
	       $(OBJS) $(LIBOBJS_LOC) $(LIBS)
	- ( cd ../../bin ; ln -sf ../src/$(DIR)/$@ . )

test_bspline.x: test_bspline.o $(OBJS) ../ctools/libagwx_ctools.a
	$(LD) $(LDFLAGS) -o $@ test_bspline.o \
	       $(OBJS) $(LIBOBJS_LOC) $(LIBS)
	- ( cd ../../bin ; ln -sf ../src/$(DIR)/$@ . )

test_poly.x: test_poly.o $(OBJS) ../ctools/libagwx_ctools.a
	$(LD) $(LDFLAGS) -o $@ test_poly.o \
	       $(OBJS) $(LIBOBJS_LOC) $(LIBS)
	- ( cd ../../bin ; ln -sf ../src/$(DIR)/$@ . )

dsbevx_omp.o:

#
# cleaning
# 
clean: 
	- /bin/rm -f *.x *.o *.a *core *.F90 *.mod work.pc work.pcl
wash:   clean
	        

# dependencies
include make.depend

