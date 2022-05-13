#!/bin/bash

LISTdelta='0.01'
outputdir='RPA_test'
# LIST_niter='1 2 3 4 5'
LIST_niter='1'
 
for delta in $LISTdelta; do

for iiter in $LIST_niter; do 

sed -i "/prefix/c\     prefix = '$outputdir'" test_scan_ctrl.in
sed -i "/sctrl_what/c\ sctrl_what= 'gf'" test_scan_ctrl.in
sed -i "/sctrl_filein/c\ sctrl_filein = '$iiter\_xgreenfunc.dat'" test_scan_ctrl.in
sed -i "/sctrl_fileoutg/c\ sctrl_fileoutg= '$iiter\_xgreenfunc_scan_ctrl_grid.dat'" test_scan_ctrl.in
sed -i "/sctrl_fileoutp/c\ sctrl_fileoutp= '$iiter\_xgreenfunc_scan_ctrl_poles.dat'" test_scan_ctrl.in
sed -i "/sctrl_verb/c\ sctrl_verb= .true." test_scan_ctrl.in
sed -i "/sctrl_thr/c\ sctrl_thr= 1.0d-2" test_scan_ctrl.in

./heg_sgm.x < test_scan_ctrl.in > test_scan_ctrl_iter\_$iiter\_delta\_$delta.out
mv test_scan_ctrl_iter\_$iiter\_delta\_$delta.out $outputdir.save

done

done


