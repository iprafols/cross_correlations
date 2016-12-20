#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./compute_projection_correction.run ../config/DLA_DR9_bin2Mpc_80Mpc_sampleC2_cosmoAndreu.ini
time ./compute_lya1d.run ../config/DLA_DR9_bin2Mpc_80Mpc_sampleC2_cosmoAndreu.ini
time ./correlation.run ../config/DLA_DR9_bin2Mpc_80Mpc_sampleC2_cosmoAndreu.ini

echo 'start: '+$START
echo 'end: '+$(date)
