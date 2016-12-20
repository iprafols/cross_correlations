#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./test_lya_deltas.run ../config/test_lya_deltas.ini

echo 'start: '+$START
echo 'end: '+$(date)
