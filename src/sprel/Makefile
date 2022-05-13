#
#===============================
# AGWX package
#===============================
#
include ../../make.inc

DIR=sprel
LIBOBJS_LOC= ../baselib/libagwx_base.a $(LIBOBJS)

OBJS=   compute_pr_time.o, \
	error.o, \
	lin_alg.o, \
	sprel_data.o, \
	system_init.o, \
 	test_omp.o	



#
# main targets
#
all : test_omp.x 

#spr.x : $(OBJ) 
#	$(LD) $(LDFLAGS) -o $@ spr.o $(OBJS) $(LIBOBJS_LOC) $(LIBS)
#	- ( cd ../../bin ; ln -sf ../src/$(DIR)/$@ . )


test_omp.x : test_omp.o
	$(LD) $(LDFLAGS) -o $@ test_omp.o \
	$(OBJS) $(LIBOBJS_LOC) $(LIBS) - (cd ../../bin ; ln -sf ..src/$(DIR)/$@ . )

#
# cleaning
# 
clean: 
	- /bin/rm -f *.x *.o *.a *core *.F90 *.mod work.pc work.pcl
wash:   clean
	        

# dependencies
include make.depend

