#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./compute_lya1d.run ../config/DLA_DR12_bin2Mpc_80Mpc_sampleC2.ini

echo 'start: '+$START
echo 'end: '+$(date)
