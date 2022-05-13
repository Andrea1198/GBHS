#! /bin/bash 
#
# given a directory in input, 
# assumes a number of files omegaX.dat are inside
# and overimposes all the plots with a vertical shift
#

dir=$1
step_v=$2
step_h=$3
ecut=$4
iprint=$5


if [ ! -d "$dir" ]  ; then echo "ERROR: $dir does not exist"; exit 1; fi
if [ -z "$step_v" ] ; then echo "ERROR: invalid step";        exit 1; fi
if [ -z "$step_h" ] ; then step_h=0.0; fi
if [ -z "$ecut" ]   ; then ecut=1000000.0; fi
if [ -z "$iprint" ] ; then iprint=1; fi

index=1
found=yes
data=""
#
while [ "$found" = "yes" ] ; do
   #
   # index is updated at the end 
   #
   file=$dir/omega$index.dat
   if [ ! -e "$file" ] ; then found=no ; fi
   #
   shift_v=`echo $index $step_v $iprint | awk '{print ($1-1) * $2 / $3}'`
   shift_h=`echo $index $step_h $iprint | awk '{print ($1-1) * $2 / $3}'`
   #
   if [ "$found" = "yes" ] ; then
      #
      data="

      `cat $file | awk -v vshift_v=$shift_v -v vshift_h=$shift_h \
                       -v ecut=$ecut '
                      {
                        if ( match($1, "#") ) next 
                        if ( NF > 1 ) {
                           printf "%15.9f ",$1 + vshift_h
                           for ( i=2; i<= NF; i++) {
                               vpot=$(i)
                               if ( vpot > ecut ) vpot=ecut;
                               if ( vpot < -ecut ) vpot=-ecut;
                               printf "%15.9f ", vpot + vshift_v
                           }
                           printf "\n"
                        } 
                      }'`

      $data
      "
      #
   fi
   #
   let index=index+$iprint
   #
done
#
echo "$data" > $dir.dat

