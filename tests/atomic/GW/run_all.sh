#!/bin/bash
bindir=../../../bin

# run scf LDA calculation
# $bindir/agf.x < run_H.in > run_H.out
$bindir/agf.x < run_He.in > run_He.out
# $bindir/agf.x < run_Li.in > run_Li.out
# $bindir/agf.x < run_Be.in > run_Be.out
# $bindir/agf.x < run_Ne.in > run_Ne.out
# $bindir/agf.x < run_Na.in > run_Na.out
# $bindir/agf.x < run_Mg.in > run_Mg.out
# $bindir/agf.x < run_Ar.in > run_Ar.out
# $bindir/agf.x < run_K.in > run_K.out
# $bindir/agf.x < run_Ca.in > run_Ca.out
# $bindir/agf.x < run_Kr.in > run_Kr.out

# run epsilon calculations
# $bindir/agf.x < run_H_epsilon.in > run_H_epsilon.out
# $bindir/agf.x < run_He_epsilon.in > run_He_epsilon.out
# $bindir/agf.x < run_Li_epsilon.in > run_Li_epsilon.out
# $bindir/agf.x < run_Be_epsilon.in > run_Be_epsilon.out
# $bindir/agf.x < run_Ne_epsilon.in > run_Ne_epsilon.out
# $bindir/agf.x < run_Na_epsilon.in > run_Na_epsilon.out
# $bindir/agf.x < run_Mg_epsilon.in > run_Mg_epsilon.out
# $bindir/agf.x < run_Ar_epsilon.in > run_Ar_epsilon.out
# $bindir/agf.x < run_K_epsilon.in > run_K_epsilon.out
# $bindir/agf.x < run_Ca_epsilon.in > run_Ca_epsilon.out
# $bindir/agf.x < run_Kr_epsilon.in > run_Kr_epsilon.out

# run GW calculations
# $bindir/agf.x < run_H_GW.in > run_H_GW.out
$bindir/agf.x < run_He_GW.in > run_He_GW.out
# $bindir/agf.x < run_Li_GW.in > run_Li_GW.out
# $bindir/agf.x < run_Be_GW.in > run_Be_GW.out
# $bindir/agf.x < run_Ne_GW.in > run_Ne_GW.out
# $bindir/agf.x < run_Na_GW.in > run_Na_GW.out
# $bindir/agf.x < run_Mg_GW.in > run_Mg_GW.out
# $bindir/agf.x < run_Ar_GW.in > run_Ar_GW.out
# $bindir/agf.x < run_K_GW.in > run_K_GW.out
# $bindir/agf.x < run_Ca_GW.in > run_Ca_GW.out
# $bindir/agf.x < run_Kr_GW.in > run_Kr_GW.out

