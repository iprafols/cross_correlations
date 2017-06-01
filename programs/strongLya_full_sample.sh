#!/bin/bash

START=$(date)

export OMP_NUM_THREADS=8

time ./compute_projection_correction.run ../config/strongLya_full_sample.ini
time ./compute_lya1d.run ../config/strongLya_full_sample.ini
time ./correlation.run ../config/strongLya_full_sample.ini

echo 'start: '+$START
echo 'end: '+$(date)