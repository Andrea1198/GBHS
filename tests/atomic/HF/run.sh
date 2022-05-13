#! /bin/bash 
#
# HF
#
#================================================================
#
# Input flags for this script (./run.sh FLAG): 
#
MANUAL=" Usage
   run.sh [FLAG]

 where FLAG is one of the following 
 (no FLAG will print this manual page) :
 
 H
 He, He_spin
 Li, B, Be, C
 Ne, Ne_spin
 Na, Mg, Ar, 
 K,  Ca, Kr

 H_epsilon
 He_epsilon, He_spin_epsilon
 Li_epsilon, Be_epsilon, F_epsilon, Ne_epsilon
 Na_epsilon, Mg_epsilon, Ar_epsilon
 K_epsilon,  Ca_epsilon, Kr_epsilon

 all             perform all the above described steps

 check           check results with the reference outputs
 clean           delete all output files and the temporary directory
"
#
#================================================================
#

#
# source common enviroment
. ../environment.conf
#

#
# evaluate the starting choice about what is to run 


run_H=
run_He=
run_He_spin=

run_Li=
run_Be=
run_B=
run_C=
run_F=
run_Ne=
run_Ne_spin=

run_Na=
run_Mg=
run_Ar=

run_K=
run_Ca=
run_Kr=
#
run_H_eps=
run_He_eps=
run_He_spin_eps=

run_Li_eps=
run_Be_eps=
run_F_eps=
run_Ne_eps=

run_Na_eps=
run_Mg_eps=
run_Ar_eps=

run_K_eps=
run_Ca_eps=
run_Kr_eps=

#
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
#INPUT=`echo $1 | tr [:upper:] [:lower:]`
INPUT=`echo $1`

case $INPUT in 
   (H)              run_H=yes ;;
   (He)             run_He=yes ;;
   (He_spin)        run_He_spin=yes ;;

   (Li)             run_Li=yes ;;
   (Be)             run_Be=yes ;;
   (B)              run_B=yes ;;
   (C)              run_C=yes ;;
   (F)              run_F=yes ;;
   (Ne)             run_Ne=yes ;;
   (Ne_spin)        run_Ne_spin=yes ;;

   (Na)             run_Na=yes ;;
   (Mg)             run_Mg=yes ;;
   (Ar)             run_Ar=yes ;;

   (K)              run_K=yes ;;
   (Ca)             run_Ca=yes ;;
   (Kr)             run_Kr=yes ;;
   #
   (H_epsilon)      run_H_eps=yes ;;
   (He_epsilon)     run_He_eps=yes ;;
   (He_spin_epsilon)run_He_spin_eps=yes ;;

   (Li_epsilon)     run_Li_eps=yes ;;
   (Be_epsilon)     run_Be_eps=yes ;;
   (F_epsilon)      run_F_eps=yes  ;;
   (Ne_epsilon)     run_Ne_eps=yes ;;

   (Na_epsilon)     run_Na_eps=yes ;;
   (Mg_epsilon)     run_Mg_eps=yes ;;
   (Ar_epsilon)     run_Ar_eps=yes ;;

   (K_epsilon)      run_K_eps=yes  ;;
   (Ca_epsilon)     run_Ca_eps=yes ;;
   (Kr_epsilon)     run_Kr_eps=yes ;;
   #
   (all)            run_H=yes ;  run_He=yes ; run_He_spin=yes ;
                    run_Li=yes ; run_Be=yes ; run_B=yes  ; run_C=yes ; run_F=yes ; run_Ne=yes ;
                    run_Na=yes ; run_Mg=yes ; run_Ar=yes ;
                    run_K=yes ;  run_Ca=yes ; run_Kr=yes ;
                    run_H_eps=yes  ; run_He_eps=yes ; run_He_spin_eps=yes ;
                    run_Li_eps=yes ; run_Be_eps=yes ; run_F_eps=yes ; run_Ne_eps=yes ;
                    run_Na_eps=yes ; run_Mg_eps=yes ; run_Ar_eps=yes ;
                    run_K_eps=yes  ; run_Ca_eps=yes ; run_Kr_eps=yes ;;
   (check)          CHECK=yes ;;
   (clean)          CLEAN=yes ;;
   (*)              echo " Invalid input FLAG, type ./run.sh for help" ; exit 1 ;;
esac

#
# initialize 
#
if [ -z "$CLEAN" ] ; then
   if [ ! -d "$TMPDIR" ] ; then echo "ERROR: no tmp dir"; exit 1 ; fi
   if [ ! -e ./SCRATCH ] ; then ln -s $TMPDIR ./SCRATCH ; fi
fi
#


#-----------------------------------------------------------------------------


if [ "$run_H" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_H.in > run_H.out
fi
#
if [ "$run_He" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_He.in > run_He.out
fi
if [ "$run_He_spin" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_He_spin.in > run_He_spin.out
fi
#
if [ "$run_Li" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Li.in > run_Li.out
fi
#
if [ "$run_Be" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Be.in > run_Be.out
fi
#
if [ "$run_B" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_B.in > run_B.out
fi
#
if [ "$run_C" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_C.in > run_C.out
fi
#
if [ "$run_F" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_F.in > run_F.out
fi
#
if [ "$run_Ne" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ne.in > run_Ne.out
fi
#
if [ "$run_Ne_spin" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ne_spin.in > run_Ne_spin.out
fi
#
if [ "$run_Na" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Na.in > run_Na.out
fi
#
if [ "$run_Mg" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Mg.in > run_Mg.out
fi
#
if [ "$run_Ar" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ar.in > run_Ar.out
fi
#
if [ "$run_K" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_K.in > run_K.out
fi
#
if [ "$run_Ca" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ca.in > run_Ca.out
fi
#
if [ "$run_Kr" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Kr.in > run_Kr.out
fi

#
# EPSILON CALCs
#
if [ "$run_H_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_H_epsilon.in > run_H_epsilon.out
fi
#
if [ "$run_He_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_He_epsilon.in > run_He_epsilon.out
fi
#
if [ "$run_He_spin_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_He_spin_epsilon.in > run_He_spin_epsilon.out
fi
#
if [ "$run_Li_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Li_epsilon.in > run_Li_epsilon.out
fi
#
if [ "$run_Be_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Be_epsilon.in > run_Be_epsilon.out
fi
#
if [ "$run_F_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_F_epsilon.in > run_F_epsilon.out
fi
#
if [ "$run_Ne_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ne_epsilon.in > run_Ne_epsilon.out
fi
#
if [ "$run_Na_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Na_epsilon.in > run_Na_epsilon.out
fi
#
if [ "$run_Mg_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Mg_epsilon.in > run_Mg_epsilon.out
fi
#
if [ "$run_Ar_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ar_epsilon.in > run_Ar_epsilon.out
fi
#
if [ "$run_K_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_K_epsilon.in > run_K_epsilon.out
fi
#
if [ "$run_Ca_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ca_epsilon.in > run_Ca_epsilon.out
fi
#
if [ "$run_Kr_eps" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Kr_epsilon.in > run_Kr_epsilon.out
fi


#
# running CHECK
#
if [ "$CHECK" = yes ] ; then
   echo "running CHECK... "
   #
   cd $TEST_HOME
   list="disentangle$SUFFIX.out wannier$SUFFIX.out"
   #
   for file in $list
   do
      ../../script/check.sh $file
   done
fi


##
## eventually clean
##
#run_clean  RUN=$CLEAN


#
# exiting
exit 0


