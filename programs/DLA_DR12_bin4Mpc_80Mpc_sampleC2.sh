#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./correlation.run ../config/DLA_DR12_bin4Mpc_80Mpc_sampleC2.ini

echo 'start: '+$START
echo 'end: '+$(date)
