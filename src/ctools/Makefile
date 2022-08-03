#
#===============================
# AGWX package
#===============================
#
include ../../make.inc

OBJS= \
c_mkdir.o \
c_printing.o \
cptimer.o \
fft_stick.o \
io.o \
pole_fitting_pywrap.o \
gcube2plt.o \
python_init.o



#
# main targets
#
all : libagwx_ctools.a


libagwx_ctools.a : $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

pole_fitting_pywrap.o: pole_fitting.str
pole_fitting.str:
	./py2str -o pole_fitting.str pole_fitting.py

#
# cleaning
# 
clean: 
	- /bin/rm -f *.o *.a *core *.F90 *.mod work.pc work.pcl *.str

