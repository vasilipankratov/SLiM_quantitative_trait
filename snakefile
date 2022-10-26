

import numpy as np
import pandas as pd



configfile: "config.yaml"




## Read in the file with meta data on the SNPs (gwas summary stats + some other info like daf)            
variable_parameters = pd.read_table("variable_parameters.csv").set_index("setup", drop=False)
variable_parameters.index = variable_parameters.index.map(str)

ext_vcf = [".vcf.gz",".vcf.gz.tbi", ".samples.txt"]

output_path = "results/msprime_" + str(config["msprime_seed"]) + \
        "/SLiM_burnin_{slim_burnin_ID}" + \
        "/demography_output/setup_{setup}" + \
        "_run-{runID}" + \
        "/N-" + str(config["initial_Ne"]) + \
        "_l-" + str(config["len_chr"]) + \
        "_n-" + str(config["n_chr"]) + \
        "_opt-" + str(config["trait_opt"]) + \
        "_w-{fitness_peak_w}" + \
        "_sd_beta-" + str(config["sd_beta"]) + \
        "_VE-" + str(config["env_var"]) + \
        "_SLiM_burnin_seed-{slim_burnin_ID}" + \
        "_SLiM.demography-{runID}" + \
        "_all"
        
        
output_vcfs = expand(output_path,
        zip,
        slim_burnin_ID = variable_parameters.loc[:,"slim_burnin_seed"],
        fitness_peak_w = variable_parameters.loc[:,"w"],
        setup = variable_parameters.loc[:,"setup"],
        runID = variable_parameters.loc[:,"slim_demography_seed"])

output_vcfs = [path + "{ext}" for path in output_vcfs]



rule all:
    input:
    	expand(output_vcfs,
        ext = ext_vcf
        )



rule run_msprime:
    output:
    	"results/msprime_{msprime_seed}/msprime_burn-in_Ne_{initial_Ne}_length_{n_chr}x{len_chr}_seed_{msprime_seed}.sim"
    resources:
        mem = 5000,
        mem_mb = 5000,
        disk_mb = 5000,
        time = 3*60
    shell:
        '''
        cd results/msprime_{wildcards.msprime_seed}
        python3 ../../neutral_burn-in.py \
        -N {wildcards.initial_Ne} \
        -l {wildcards.len_chr} \
        -n {wildcards.n_chr} \
        -s {wildcards.msprime_seed}
        '''

rule run_SLiM_burnin:
    input:
    	"results/msprime_{msprime_seed}/msprime_burn-in_Ne_{initial_Ne}_length_{n_chr}x{len_chr}_seed_" + str(config["msprime_seed"]) + ".sim"
    output:
    	multiext("results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.burn-in",
    	".sim",".stats.tsv",".QTLs.freq.tsv",".QTLs.list.tsv")
    params:
        mu = config["mu"],
        rec = config["rec"],
        dir = "results/msprime_" + str(config["msprime_seed"])
    resources:
        mem = 15000,
        time = 5*60
    shell:
    	'''
    	cd {params.dir}
    	slim \
    	-s {wildcards.slim_burnin_seed} \
    	-d N={wildcards.initial_Ne} \
    	-d u={params.mu} \
    	-d r={params.rec} \
    	-d sd_beta={wildcards.sd_beta} \
    	-d V_E={wildcards.env_var} \
    	-d w={wildcards.fitness_w} \
    	-d opt={wildcards.trait_opt} \
    	-d l_chr={wildcards.len_chr} \
    	-d n_chr={wildcards.n_chr} \
    	-d msprime_seed={wildcards.msprime_seed} \
    	../../read_msprime_run_burnin.slim 
    	mv *SLiM_burnin_seed-{wildcards.slim_burnin_seed}_SLiM.burn-in* SLiM_burnin_{wildcards.slim_burnin_seed}/
    	'''



    	
rule run_SLiM_demography:
    input:
        "results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.burn-in.sim"
    output:
    	multiext("results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/demography_output/setup_{setup}_run-{slim_demography_seed}/N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.demography-{slim_demography_seed}",
    	".stats.tsv",".QTLs.freq.tsv",".QTLs.list.tsv",".asPGS.tsv",".param.list.tsv",".p7.vcf",".p31.vcf",".p41.vcf",".p51.vcf")
    params:
        mu = config["mu"],
        rec = config["rec"],
        deviating_pop = lambda wildcards: variable_parameters.loc[wildcards.setup, "deviating_pop"],
        opt_shifted = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_dev_pop"],
        opt_UK = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_UK"],
        w_shifted = lambda wildcards: variable_parameters.loc[wildcards.setup, "fitness_w_dev_pop"],
        w_UK = lambda wildcards: variable_parameters.loc[wildcards.setup, "fitness_w_UK"],
        dir = lambda wildcards: variable_parameters.loc[wildcards.setup, "slim_burnin_seed"]
    resources:
        mem = 30000,
        time = 10*60
    shell:
        '''
        cd "results/msprime_{wildcards.msprime_seed}/SLiM_burnin_{params.dir}"
        slim \
        -s {wildcards.slim_demography_seed} \
        -d setup={wildcards.setup} \
        -d N={wildcards.initial_Ne} \
        -d u={params.mu} \
        -d r={params.rec} \
        -d sd_beta={wildcards.sd_beta} \
        -d V_E={wildcards.env_var} \
        -d w={wildcards.fitness_w} \
        -d opt={wildcards.trait_opt} \
        -d deviating_pop={params.deviating_pop} \
        -d opt_shifted={params.opt_shifted} \
        -d opt_UK={params.opt_UK} \
        -d w_shifted={params.w_shifted} \
        -d w_UK={params.w_UK} \
        -d l_chr={wildcards.len_chr} \
        -d n_chr={wildcards.n_chr} \
        -d burn_in_ID={wildcards.slim_burnin_seed} \
        ../../../demography.slim 
        '''



rule handle_vcfs:
    input:
    	multiext("results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/demography_output/setup_{setup}_run-{slim_demography_seed}/N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.demography-{slim_demography_seed}",
    	".p7.vcf",".p31.vcf",".p41.vcf",".p51.vcf")
    output:
        multiext("results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/demography_output/setup_{setup}_run-{slim_demography_seed}/N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.demography-{slim_demography_seed}_all",
    	".vcf.gz",".vcf.gz.tbi",".samples.txt")
    params:
        demodir = "results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/demography_output/setup_{setup}_run-{slim_demography_seed}/",
        file = "N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.demography-{slim_demography_seed}"
    resources:
        mem = 5000,
        time = 1*60
    shell:
        '''
        module load tabix
        module load bcftools/1.9
        cd {params.demodir}
        for pop in p31 p41 p51 p7
        do
            bgzip {params.file}.${{pop}}.vcf
            bcftools query -l {params.file}.${{pop}}.vcf.gz | sed 's/_//g' | awk -F'\t' -v pop=$pop 'BEGIN {{OFS = FS}} {{print $1, pop "_" $1}}' > samples.${{pop}}.txt
            bcftools reheader -s samples.${{pop}}.txt {params.file}.${{pop}}.vcf.gz | bcftools view -Oz -i 'MT=1' -o m1_{params.file}.${{pop}}.vcf.gz
            tabix -p vcf m1_{params.file}.${{pop}}.vcf.gz
        done
        bcftools merge -Oz --missing-to-ref -o {params.file}_all.vcf.gz m1_{params.file}.{{p31,p41,p51,p7}}.vcf.gz
        tabix -p vcf {params.file}_all.vcf.gz
        rm m1_{params.file}.{{p31,p41,p51,p7}}.vcf.gz*
        mkdir vcf_backup
        mv {params.file}.{{p31,p41,p51,p7}}.vcf.gz* vcf_backup
        mv samples.{{p31,p41,p51,p7}}.txt vcf_backup
        bcftools query -l {params.file}_all.vcf.gz | awk '{print $1, $1}' | awk -v OFS='\t' '{sub(/_i[0-9]*/, "", $2)} 1' > {params.file}_all.samples.txt
        '''
    	




