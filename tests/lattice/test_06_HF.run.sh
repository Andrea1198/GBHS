#! /bin/bash

prefix=06_HF
filetmpl=test_06_HF.in.in
bindir=./

#
# parameters
#
U_list="1.0 2.0 4.0 6.0 8.0 10.0 12.0"
V_list="1.0"

#
if [ "$1" == "-c" ] ; then clean=yes ; fi
if [ "$1" == "-d" ] ; then dry_run=yes ; fi

#
# main loop
#
data_rho=
data_energy=
#
for U in $U_list ; do
for V in $V_list ; do
  #
  label="U${U}_V${V}"
  filein=test_${prefix}.${label}.in
  fileout=test_${prefix}.${label}.out
  savedir=run${prefix}_${label}.save
  #
  # input generation
  if [ -e "$filein" ] ;  then rm $filein ; fi
  if [ -e "$fileout" ] ; then rm $fileout ; fi
  if [ -d "$savedir" ] ; then rm -rf $savedir ; fi
  #
  if [ "$clean" != "yes" ] ; then
     #
     sed "s/@U@/$U/g
          s/@V@/$V/g
          s/@label@/$label/" $filetmpl > $filein
     #
     # running
     if [ "$dry_run" != "yes" ] ; then
        echo "Running $label"
        $bindir/lattice.x < $filein > $fileout
     fi
     #
  fi
  #
  # parse data
  #
  if [ -e "$fileout" ] ; then 
     rhodat=`grep rho $fileout | tail -1`
     etotdat=`grep e_tot $fileout | tail -1`
     esgmdat=`grep e_sgm $fileout | tail -1`
     #
     data_rho=`echo "$data_rho";       echo "$U $V   $rhodat"`
     data_energy=`echo "$data_energy"; echo "$U $V   $etotdat   $esgmdat"`
  fi
  #
done
done

if [ "$clean" != "yes" -a "$dry_run" != "yes" ] ; then
  echo "Data Rho"
  echo "$data_rho"
  echo
  echo "Data ENERGY"
  echo "$data_energy"
fi

exit 0
