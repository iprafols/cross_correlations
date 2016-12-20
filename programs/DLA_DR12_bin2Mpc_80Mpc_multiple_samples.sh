1;2c#!/bin/bash

START=$(date)

# run sample A
#time ./DLA_DR12_bin2Mpc_80Mpc_sampleA.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleA.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleA.err

# run sample C2 (DR9)
#time ./DLA_DR9_bin2Mpc_80Mpc_sampleC2.sh > ../logs/DLA_DR9_bin2Mpc_80Mpc_sampleC2.log 2> ../logs/DLA_DR9_bin2Mpc_80Mpc_sampleC2.err 

# run sample N1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN1_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN1_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN1_fromC1.err

# run sample N2
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.err

# run sample N3
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN3_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN3_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN3_fromC1.err

# run sample Z1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ1_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ1_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ1_fromC1.err

# run sample Z2
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ2_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ2_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ2_fromC1.err

# run sample Z3
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ3_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ3_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ3_fromC1.err

echo 'start: '+$START
echo 'end: '+$(date)
