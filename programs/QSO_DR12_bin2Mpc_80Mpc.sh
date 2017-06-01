#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./compute_projection_correction.run ../config/QSO_DR12_bin2Mpc_80Mpc.ini
time ./compute_lya1d.run ../config/QSO_DR12_bin2Mpc_80Mpc.ini
time ./correlation.run ../config/QSO_DR12_bin2Mpc_80Mpc.ini

echo 'start: '+$START
echo 'end: '+$(date)
