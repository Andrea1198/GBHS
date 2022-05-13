#! /bin/bash 
#
# LDA
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
 Li
 Be
 Ne, Ne_spin 
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
run_Ne=
run_Ne_spin=
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
   (Ne)             run_Ne=yes ;;
   (Ne_spin)        run_Ne_spin=yes ;;
   (all)            run_H=yes ;  run_He=yes ; run_He_spin=yes ; 
                    run_Li=yes ; run_Be=yes ; run_Ne=yes ; run_Ne_spin=yes ;;
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
if [ "$run_Ne" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ne.in > run_Ne.out
fi
#
if [ "$run_Ne_spin" = "yes" ] ; then
    $AGWX_BIN/agf.x < run_Ne_spin.in > run_Ne_spin.out
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


