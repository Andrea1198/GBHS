#!/bin/bash

RS="4"
UNITS='auxfunits'
LIST_delta='0.01'
OUT_DIR='RPA_test'
LIST_plotfunc='sigma'
# LIST_k='0.5 1.0 1.4'
LIST_k='1.0'
# LIST_niter='1 2 3 4 5 6 7 8'
LIST_niter='1'
HEGSIGMAX='./heg_sgm.x'
filein='plot.in'
 
for delta in $LIST_delta; do
  for iiter in $LIST_niter; do 
    for k in $LIST_k; do
      for plotfunc in $LIST_plotfunc; do
        if [ "$plotfunc" = "sigma" ]; then
          filename="${iiter}_xsigma.dat"
          typeplot="sigma_vs_w"
        elif [ "$plotfunc" = "greenfunc" ]; then
          filename="${iiter}_xgreenfunc.dat"
          typeplot="gf_vs_w"
        elif [ "$plotfunc" = "polarizability" ]; then
          filename="${iiter}_polarizability.dat"
          typeplot="pol_vs_w"
        elif [ "$plotfunc" = "wscreen" ]; then
          filename="${iiter}_wscreen.dat"
          typeplot="wscreen_vs_w"
        fi
        fileoutplot="${plotfunc}_plot_vs_w_k_${k}_iter_${iiter}_delta_${delta}.dat"
        fileout="$OUT_DIR.save/${plotfunc}_plot_vs_w_k_${k}_iter_${iiter}_delta_${delta}.out"
        #
cat > $filein << EOF
&CONTROL
   prefix = "$OUT_DIR"
calculation = "plot"
  units_in = "$UNITS" 
     rs_min = $RS 
        nrs = 1
 /

&PLOT
   what = "$typeplot"
filein  = "$filename" 
fileout = "$fileoutplot"
plot_vals(1) = $k 
  /
EOF
        #
        $HEGSIGMAX "-in" $filein > $fileout
      done
    done
  done
done
rm $filein
