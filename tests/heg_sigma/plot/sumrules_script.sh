#!/bin/bash

RS=${1}
unitsList=${2}
KMAX=${4}
LIST_delta='0'
OUT_DIR='GW'
LIST_niter=${3}
HEGSIGMAX='../heg_sgm.x'
filein='sumrules.in'
 
cd rs_${1}
for UNITS in $unitsList; do
for delta in $LIST_delta; do
  for iiter in $LIST_niter; do 
    filename="${iiter}_xgreenfunc_sigma.dat"
    fileout="$OUT_DIR.save/${iiter}_sumrules_delta_${delta}_units_${UNITS}.out"
    filein=${iiter}_sumrules_delta_${delta}_units_${UNITS}.in
    #
cat > ${filein} << EOF
&CONTROL
   prefix = "$OUT_DIR"
calculation = "sumrules"
  units_in = "$UNITS" 
     rs_min = $RS 
        nrs = 1
 /

&SUM_RULES
  sumrules_calc = .true.
  sumrules_smeardos_fact = ${delta}
  sumrules_smeardos_type = 'gaussian'
  sumrules_filein =   '$filename'
  sumrules_fileoutA=  '${iiter}_spectrfunc_delta_${delta}_units_${UNITS}.dat'
  sumrules_fileoutV=  '${iiter}_spectrpot_delta_${delta}_units_${UNITS}.dat'
  sumrules_fileoutnk= '${iiter}_occ_number_delta_${delta}_units_${UNITS}.dat'
  sumrules_fileoutwgt='${iiter}_spectralwgt_delta_${delta}_units_${UNITS}.dat'
  sumrules_kmax      =${KMAX}
  /
EOF
        #
        $HEGSIGMAX "-in" $filein > $fileout & 
    done
done
done
cd ..

