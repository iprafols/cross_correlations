#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=4

time ./compute_projection_correction.run ../config/QSO_DR12_bin4Mpc_200Mpc.ini
time ./compute_lya1d.run ../config/QSO_DR12_bin4Mpc_200Mpc.ini
time ./correlation.run ../config/QSO_DR12_bin4Mpc_200Mpc.ini

echo 'start: '+$START
echo 'end: '+$(date)
