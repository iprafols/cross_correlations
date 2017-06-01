#!/bin/bash

START=$(date)

# run sample A
time ./DLA_DR12_bin2Mpc_80Mpc_sampleA.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleA.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleA.err

# run sample C1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleC1.err

# run sample C2
time ./DLA_DR12_bin2Mpc_80Mpc_sampleC2.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleC2.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleC2.err

# run sample N1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN1.err

# run sample N2
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN2.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN2.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN2.err

# run sample N3
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN3.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN3.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN3.err

# run sample Z1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ1.err

# run sample Z2
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ2.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ2.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ2.err

# run sample Z3
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ3.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ3.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ3.err

# run sample N1_fromC1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN1_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN1_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN1_fromC1.err

# run sample N2_fromC1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.err

# run sample N3_fromC1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleN3_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN3_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleN3_fromC1.err

# run sample Z1_fromC1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ1_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ1_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ1_fromC1.err

# run sample Z2_fromC1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ2_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ2_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ2_fromC1.err

# run sample Z3_fromC1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleZ3_fromC1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ3_fromC1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleZ3_fromC1.err

# run sample S1
time ./DLA_DR12_bin2Mpc_80Mpc_sampleS1.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS1.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS1.err

# run sample S2
time ./DLA_DR12_bin2Mpc_80Mpc_sampleS2.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS2.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS2.err

# run sample S3
time ./DLA_DR12_bin2Mpc_80Mpc_sampleS3.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS3.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS3.err

# run sample S4
time ./DLA_DR12_bin2Mpc_80Mpc_sampleS4.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS4.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS4.err

# run sample S1_LSN
time ./DLA_DR12_bin2Mpc_80Mpc_sampleS1_LSN.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS1_LSN.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleS1_LSN.err

# run sample SA
time ./DLA_DR12_bin2Mpc_80Mpc_sampleSA.sh > ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleSA.log 2> ../logs/DLA_DR12_bin2Mpc_80Mpc_sampleSA.err


echo 'start: '+$START
echo 'end: '+$(date)
