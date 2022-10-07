#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem=60g
#SBATCH --job-name=slim
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vasilipankratov@gmail.com

module load any/slim/4.0


n_chr=$1
Ne=$2
w=$3

l_chr=20000
msprime_seed=9999
SLiM_seed=7777777
##################################################################
# msprime
##################################################################

#START=$(date +%s)

#python3 neutral_burn-in.py -N $Ne -l $l_chr -n $n_chr -s $msprime_seed

#END=$(date +%s)

#RUNTIME=$(echo "$END - $START" | bc -l)

#echo "msprime run finished"
#printf "%10s %5d\n %10s %5d\n %15s %5d\n" "N of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME


##################################################################
# SLiM burn-in
##################################################################


START=$(date +%s)

slim -s $SLiM_seed -d N=$Ne -d u=1.25e-8 -d r=1e-8 -d sd_beta=0.1 -d V_E=0.4 -d w=$w -d opt=0 -d l_chr=$l_chr -d n_chr=$n_chr -d msprime_seed=$msprime_seed take_msprime_out.slim 

END=$(date +%s)

RUNTIME=$(echo "$END - $START" | bc -l)

echo "SLiM run finished"
printf "%10s %5d\n %10s %5d\n %15s %5d\n" "N of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME


##################################################################
# SLiM demography
##################################################################


START=$(date +%s)

slim -d N=$Ne -d u=1.25e-8 -d r=1e-8 -d sd_beta=0.1 -d V_E=0.4 -d w=$w -d opt=0 -d l_chr=$l_chr -d n_chr=$n_chr -d burn_in_ID=$SLiM_seed demography.slim 

END=$(date +%s)

RUNTIME=$(echo "$END - $START" | bc -l)

echo "SLiM run finished"
printf "%10s %5d\n %10s %5d\n %15s %5d\n" "N of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME
