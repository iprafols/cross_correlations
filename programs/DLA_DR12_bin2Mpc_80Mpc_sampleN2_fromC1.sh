#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./compute_projection_correction.run ../config/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.ini
time ./compute_lya1d.run ../config/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.ini
time ./correlation.run ../config/DLA_DR12_bin2Mpc_80Mpc_sampleN2_fromC1.ini

echo 'start: '+$START
echo 'end: '+$(date)
