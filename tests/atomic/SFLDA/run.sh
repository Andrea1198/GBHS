#! /bin/bash 
#
# SFLDA
#
#================================================================
#
# Input flags for this script (./run.sh FLAG): 
#
MANUAL=" Usage
   run.sh [FLAG]

 where FLAG is one of the following 
 (no FLAG will print this manual page) :
 
 He, 
 Be, Ne
 Mg, Ar, Kr

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


run_He=
run_Be=
run_Ne=

run_Mg=
run_Ar=
#
run_Ca=
run_Kr=
#
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
#INPUT=`echo $1 | tr [:upper:] [:lower:]`
INPUT=`echo $1`

case $INPUT in 
   (He)             run_He=yes ;;
   (Be)             run_Be=yes ;;
   (Ne)             run_Ne=yes ;;
   #
   (Mg)             run_Mg=yes ;;
   (Ar)             run_Ar=yes ;;
   #
   (Ca)             run_Ca=yes ;;
   (Kr)             run_Kr=yes ;;
   #
   (all)            run_He=yes ;
                    run_Be=yes ; run_Ne=yes ;
                    run_Mg=yes ; run_Ar=yes ; 
                    run_Ca=yes ; run_Kr=yes ;;
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


if [ "$run_He" = "yes" ] ; then
   echo "running SFLDA He... "
   $AGWX_BIN/agf.x < run_He_lda.in > run_He_lda.out
   $AGWX_BIN/agf.x < run_He.in > run_He.out
fi
#
if [ "$run_Be" = "yes" ] ; then
   echo "running SFLDA Be... "
   $AGWX_BIN/agf.x < run_Be_lda.in > run_Be_lda.out
   $AGWX_BIN/agf.x < run_Be.in > run_Be.out
fi
if [ "$run_Ne" = "yes" ] ; then
   echo "running SFLDA Ne... "
   $AGWX_BIN/agf.x < run_Ne_lda.in > run_Ne_lda.out
   $AGWX_BIN/agf.x < run_Ne.in > run_Ne.out
fi
#
if [ "$run_Mg" = "yes" ] ; then
   echo "running SFLDA Mg... "
   $AGWX_BIN/agf.x < run_Mg_lda.in > run_Mg_lda.out
   $AGWX_BIN/agf.x < run_Mg.in > run_Mg.out
fi
if [ "$run_Ar" = "yes" ] ; then
   echo "running SFLDA Ar... "
   $AGWX_BIN/agf.x < run_Ar_lda.in > run_Ar_lda.out
   $AGWX_BIN/agf.x < run_Ar.in > run_Ar.out
fi
#
if [ "$run_Ca" = "yes" ] ; then
   echo "running SFLDA Ca... "
   $AGWX_BIN/agf.x < run_Ca_lda.in > run_Ca_lda.out
   $AGWX_BIN/agf.x < run_Ca.in > run_Ca.out
fi
if [ "$run_Kr" = "yes" ] ; then
   echo "running SFLDA Kr... "
   $AGWX_BIN/agf.x < run_Kr_lda.in > run_Kr_lda.out
   $AGWX_BIN/agf.x < run_Kr.in > run_Kr.out
fi

#
# running CHECK
#
if [ "$CHECK" = yes ] ; then
   echo "running CHECK... "
   #
   cd $TEST_HOME
   list=" "
   #
   for file in $list
   do
      ../../script/check.sh $file
   done
fi


#
# eventually clean
#
if [ "$CLEAN" = yes ] ; then
  echo "running CLEAN... "
  rm -rf *.dat *.out core* CRASH
fi

#
# exiting
exit 0


