#
#===============================
# AGWX package
#===============================
#
include ../make.inc


#
# MAIN target
#
all: lib_ctools lib_base lib_atomic \
     lib_onedim lib_lattice lib_heg lib_heg_sigma

#
# drv
#
lib_ctools : 
	 if test -d ctools ;     then ( cd ctools;    $(MAKE) ) ; fi

lib_base : ctools
	 if test -d baselib ;    then ( cd baselib;   $(MAKE) ) ; fi

lib_atomic : lib_base
	 if test -d atomic ;     then ( cd atomic;    $(MAKE) ) ; fi

lib_onedim : lib_base
	 if test -d onedim ;     then ( cd onedim;    $(MAKE) ) ; fi

lib_heg : lib_base
	 if test -d heg ;        then ( cd heg;       $(MAKE) ) ; fi

lib_heg_sigma : lib_base lib_lattice
	 if test -d heg_sigma ;  then ( cd heg_sigma; $(MAKE) ) ; fi

lib_lattice : lib_base
	 if test -d lattice ;    then ( cd lattice;   $(MAKE) ) ; fi

#
# cleaning
#
clean wash:
	 if test -d ctools ;     then ( cd ctools;     $(MAKE) $@ ) ; fi
	 if test -d baselib ;    then ( cd baselib;    $(MAKE) $@ ) ; fi
	 if test -d atomic ;     then ( cd atomic;     $(MAKE) $@ ) ; fi
	 if test -d onedim ;     then ( cd onedim;     $(MAKE) $@ ) ; fi
	 if test -d heg ;        then ( cd heg;        $(MAKE) $@ ) ; fi
	 if test -d heg_sigma ;  then ( cd heg_sigma;  $(MAKE) $@ ) ; fi
	 if test -d lattice ;    then ( cd lattice;    $(MAKE) $@ ) ; fi


