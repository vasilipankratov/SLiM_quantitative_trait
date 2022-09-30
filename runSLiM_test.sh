#!/bin/bash

n_chr=$1
Ne=$2

l_chr=20000
msprime_seed=999
##################################################################
# msprime
##################################################################

START=$(date +%s)

python3 neutral_burn-in.py -N $Ne -l $l_chr -n $n_chr -s $msprime_seed

END=$(date +%s)

RUNTIME=$(echo "$END - $START" | bc -l)

echo "msprime run finished"
printf "%10s %5d\n %10s %5d\n %15s %5d\n" "msprime run finished \nN of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME


##################################################################
# SLiM
##################################################################

START=$(date +%s)

slim -d N=$Ne -d u=1.25e-8 -d r=1e-8 -d sd_beta=0.1 -d h2=0.8 -d w=INF -d opt=0 -d l_chr=$l_chr -d n_chr=$n_chr -d msprime_seed=$msprime_seed take_msprime_out.slim 

END=$(date +%s)

RUNTIME=$(echo "$END - $START" | bc -l)

echo "SLiM run finished"
printf "%10s %5d\n %10s %5d\n %15s %5d\n" "N of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME
