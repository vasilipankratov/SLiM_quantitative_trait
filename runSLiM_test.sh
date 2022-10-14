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

V_E=0.5
opt=0.0
sd_beta=0.1
l_chr=20000
msprime_seed=9999
SLiM_burnin_seed=7777777
SLiM_demo_seed=1234567891011
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

slim -s $SLiM_burnin_seed -d N=$Ne -d u=1.25e-8 -d r=1e-8 -d sd_beta=$sd_beta -d V_E=$V_E -d w=$w -d opt=$opt -d l_chr=$l_chr -d n_chr=$n_chr -d msprime_seed=$msprime_seed read_msprime_run_burnin.slim 

END=$(date +%s)

RUNTIME=$(echo "$END - $START" | bc -l)

echo "SLiM run finished"
printf "%10s %5d\n %10s %5d\n %15s %5d\n" "N of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME


##################################################################
# SLiM demography
##################################################################


START=$(date +%s)

slim -s $SLiM_demo_seed -d N=$Ne -d u=1.25e-8 -d r=1e-8 -d sd_beta=$sd_beta -d V_E=$V_E -d w=$w -d opt=$opt -d l_chr=$l_chr -d n_chr=$n_chr -d burn_in_ID=$SLiM_burnin_seed demography.slim 

END=$(date +%s)

RUNTIME=$(echo "$END - $START" | bc -l)

echo "SLiM run finished"
printf "%10s %5d\n %10s %5d\n %15s %5d\n" "N of chromosomes: " $n_chr "Pop size: " $Ne "Run time (in sec): " $RUNTIME



##################################################################
# vcf managing
##################################################################


module load tabix
module load bcftools/1.9

file=N-${Ne}_l-${l_chr}_n-${n_chr}_opt-${opt}_w-${w}_sd_beta-${sd_beta}_VE-${V_E}_SLiM_burnin_seed-${SLiM_burnin_seed}_SLiM.demography-${SLiM_demo_seed}


for pop in p31 p41 p51 p7
do
bgzip ${file}_${pop}.vcf
dir=demography_output/run-${SLiM_demo_seed}
mkdir $dir
mv ${file}_${pop}.vcf.gz ${dir}/
cd $dir
bcftools query -l ${file}_${pop}.vcf.gz | sed 's/_//g' | awk -F'\t' -v pop=$pop 'BEGIN {OFS = FS} {print $1, pop $1}' > samples_${pop}.txt
bcftools reheader -s samples_${pop}.txt ${file}_${pop}.vcf.gz | bcftools view -Oz -i 'MT=1' -o m1_${file}_${pop}.vcf.gz
tabix -p vcf m1_${file}_${pop}.vcf.gz
done

bcftools merge -Oz --missing-to-ref -o ${file}_all.vcf.gz m1_{p31,p41,p51,p7}.vcf.gz


