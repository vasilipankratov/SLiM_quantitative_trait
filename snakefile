#############################################################################
# This snakemake file describes a simulation pipeline used to test contribution of 3 european ancestral populations (WHG, Anatolian Farmers and Yamnaya) to a polygenic trait.
# The simulation consists of 3 steps corresponding to the following rules:
# 1. run_msprime - running a neutral burn-in in msprime. The msprime model describes a single constant Ne population
# 2. run_SLiM_burnin - running a burn-in in SLiM. This step is actually required if one simulates stabilizing selection: this is again a constant Ne single populaiton model, however, a polygenic trati is simulated here and the model accomodates selection on that trait which can be described by the optimal trait value and the width of the fitness peak. Such a burn-in allows to reach equilibrium between drift and selection (no new QTL mutations are currently simulated at this stage - all QTls arise during the msprime stage).
# 3. run_SLiM_demography - this rule actually describes a demographic model reflecting the origin of and subsequent admixture of the 3 componets. It allows to vary the degree of stabilzing selection as well as to shift the trait optimum for one of the ancestral populations as well as for the resulting UK population (allowing to simulate post-admixture selection).


#############################################################################
# load python modules
import numpy as np
import pandas as pd

#############################################################################
# defining the config file. It contains basic parameters of the simualtion which are the same between multiple runs
#initial_Ne: 14000 - the Ne of the original pop during msprime and slim burn-in
#n_chr: 1000 = number of "chromosomes" to simulate
#len_chr: 20000  - length of each chromosome
#mu: 1.25e-8 - mutation rate
#rec: 1e-8 - recombination rate within each chromosome
#msprime_seed: 1162971105 - random seed used in msprime
#trait_opt: 0.0 - default optimal value of the trait in all 3 stages. Can be overriden for one ancestral pop and UK if desired 
#sd_beta: 0.1 - standard deviation of the 0-centered distribution from which per snp betas 
#env_var: 0.9 - environmental variance - when modeling the phenotype the environmental noise is sampled from a 0-centered normal distribution with SD of env_var^0.5

configfile: "config.yaml"


#############################################################################
## Read in the file with parameters for each run         
# I need to do some strange moves (reading some columns as float and then turning them into str and converting to upper case) because python uses "inf" while SliM uses "INF"

# a dictionary with data types used when reading the file with parameters
datatypes = {'setup' : int, 'w' : float, 'slim_burnin_seed' : int, 'slim_demography_seed' : int, 'deviating_pop' : int, 'opt_of_dev_pop' : float, 'opt_of_UK' : float, 'fitness_w_dev_pop' : float, 'fitness_w_UK' : float}

# this file contains variable simulation parameters which wwe want to test.
# each line is a combination of parameter values specific for a given run (given "setup")
# here is a desciption of each column:

# setup - setup index
# w - width of the fitness function peak for modeling stabilizing selection. INF results in neutrality
# slim_burnin_seed - random seed for slim burn-in
# slim_demography_seed - random seed for slim demography
# deviating_pop - an ancestral population with differences of the trait optimum and/or strength of stabilizing selection
# opt_of_dev_pop - trait optimum value in that population
# opt_of_UK - optimum trait in UK pop - provides a way to simulate post-admixture selection
# fitness_w_dev_pop - fitness function width in the deviating pop
# fitness_w_UK - fitness function width in the UK population (one can do for instance relaxation of stabilizing selection)

variable_parameters = pd.read_table("variable_parameters.csv", dtype = datatypes)

# a new dictionary to change 'w-columns' to str
datatypes_w = {'w' : str, 'fitness_w_dev_pop' : str, 'fitness_w_UK' : str}

variable_parameters = variable_parameters.astype(datatypes_w)

# converting to uppercase (inf -> INF)
w_columns = ['w', 'fitness_w_dev_pop', 'fitness_w_UK']
variable_parameters[w_columns] = variable_parameters[w_columns].apply(lambda x: x.str.upper())

# setting the setup column as index and converting the index to string
variable_parameters = variable_parameters.set_index("setup", drop=False)
variable_parameters.index = variable_parameters.index.map(str)



# creating the master prefix with common variables from the config file and creating wild cards for variable parameters
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
        "_SLiM_setup_{setup}" + \
        "_seed-{runID}"

# populating the wild cards in the above <output_path> with values from the <variable_parameters> df        
output_files = expand(output_path,
        zip,
        slim_burnin_ID = variable_parameters.loc[:,"slim_burnin_seed"],
        fitness_peak_w = variable_parameters.loc[:,"w"],
        setup = variable_parameters.loc[:,"setup"],
        runID = variable_parameters.loc[:,"slim_demography_seed"])


# creating the output files names for the "all" rule

output_files = [path + "{ext}" for path in output_files]

# here one can specify the extension of the ultimate output; this should allow to easily extend the pipeline if needed
# output_ext = [".no_rel.vcf.gz"]
output_ext = expand(".h2_" + str(config["h2_cutoff"]) + ".{filetype}",
        filetype = ["asPGS+covma.csv", "pdf"])


rule all:
    input:
        expand(output_files, 
        ext = output_ext)

# this rule runs msprime simulation job which provides neutral genetic diversity at mutation-drift equlibrium to pick QTLs 

rule run_msprime:
    output:
    	"results/msprime_{msprime_seed}/msprime_burn-in_Ne_{initial_Ne}_length_{n_chr}x{len_chr}_seed_{msprime_seed}.sim"
    resources:
        mem = 30*1000,
        time = 3*60
    shell:
        '''
        wd=$(pwd)
        cd results/msprime_{wildcards.msprime_seed}
        python3 ${wd}/scripts/neutral_burn-in.py \
        -N {wildcards.initial_Ne} \
        -l {wildcards.len_chr} \
        -n {wildcards.n_chr} \
        -s {wildcards.msprime_seed}
        '''



# defining some prefices to siplify specifying I/O files

burnin_prefix = "results/msprime_{msprime_seed}/SLiM_burnin_{slim_burnin_seed}/"
demography_prefix = burnin_prefix + "demography_output/setup_{setup}_run-{slim_demography_seed}/"
file_prefix = "N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM_setup_{setup}_seed-{slim_demography_seed}"
master_prefix = demography_prefix + file_prefix


demography_prefix_for_expand = demography_prefix.replace("{", "{{").replace("}", "}}" )
file_prefix_for_expand = file_prefix.replace("{", "{{").replace("}", "}}" )


# this is a SLiM burn-in step. When we introduce stabilizing selection acting on a trait it would be nice to reach equlibrium before simulating any demography. This is what this step does. Note that it is also kept in completely neutral setups for consistency.

rule run_SLiM_burnin:
    input:
    	"results/msprime_{msprime_seed}/msprime_burn-in_Ne_{initial_Ne}_length_{n_chr}x{len_chr}_seed_" + str(config["msprime_seed"]) + ".sim"
    output:
    	multiext(burnin_prefix + "N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.burn-in",
    	".sim",".stats.tsv",".QTLs.freq.tsv",".QTLs.list.tsv")
    params:
        mu = config["mu"],
        rec = config["rec"],
        dir = "results/msprime_" + str(config["msprime_seed"])
    resources:
        mem = 30*1000,
        time = 5*60
    shell:
    	'''
    	wd=$(pwd)
    	module load any/slim/4.0
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
    	${{wd}}/scripts/read_msprime_run_burnin.slim 
    	mv *SLiM_burnin_seed-{wildcards.slim_burnin_seed}_SLiM.burn-in* SLiM_burnin_{wildcards.slim_burnin_seed}/
    	'''




    	
rule run_SLiM_demography:
    input:
        burnin_prefix + "N-{initial_Ne}_l-{len_chr}_n-{n_chr}_opt-{trait_opt}_w-{fitness_w}_sd_beta-{sd_beta}_VE-{env_var}_SLiM_burnin_seed-{slim_burnin_seed}_SLiM.burn-in.sim"
    output:
    	multiext(master_prefix,
    	".stats.tsv",".QTLs.freq.tsv",".QTLs.list.tsv",".asPGS.tsv",".param.list.tsv",".p7.vcf",".p31.vcf",".p41.vcf",".p51.vcf")
    params:
        mu = config["mu"],
        rec = config["rec"],
        deviating_pop = lambda wildcards: variable_parameters.loc[wildcards.setup, "deviating_pop"],
        opt_shifted = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_dev_pop"],
        opt_UK = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_UK"],
        opt_p6 = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_p6"],
        w_shifted = lambda wildcards: variable_parameters.loc[wildcards.setup, "fitness_w_dev_pop"],
        w_UK = lambda wildcards: variable_parameters.loc[wildcards.setup, "fitness_w_UK"],
        dir = lambda wildcards: variable_parameters.loc[wildcards.setup, "slim_burnin_seed"]
    resources:
        mem = 60*1000,
        time = 24*60
    shell:
        '''
        wd=$(pwd)
        module load any/slim/4.0
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
        -d opt_p6={params.opt_p6} \
        -d w_shifted={params.w_shifted} \
        -d w_UK={params.w_UK} \
        -d l_chr={wildcards.len_chr} \
        -d n_chr={wildcards.n_chr} \
        -d burn_in_ID={wildcards.slim_burnin_seed} \
        ${{wd}}/scripts/demography.slim 
        '''



rule handle_vcfs:
    input:
    	multiext(master_prefix,
    	".p7.vcf",".p31.vcf",".p41.vcf",".p51.vcf")
    output:
        multiext(master_prefix + "_all",
    	".vcf.gz",".vcf.gz.tbi",".samples.txt"),
        demography_prefix + "vcf_backup/" + file_prefix + ".p7.vcf.gz"
    params:
        demodir = demography_prefix,
        file = file_prefix
    resources:
        mem = 300,
        time = 1*60
    shell:
        '''
        module load tabix
        module load bcftools/1.9
        cd {params.demodir}
        for pop in p31 p41 p51 p7
        do
            bgzip {params.file}.${{pop}}.vcf
            tabix -p vcf {params.file}.${{pop}}.vcf.gz
            bcftools query -l {params.file}.${{pop}}.vcf.gz | sed 's/_//g' | awk -F'\\t' -v pop=$pop 'BEGIN {{OFS = FS}} {{print $1, pop "_" $1}}' > samples.${{pop}}.txt
            bcftools reheader -s samples.${{pop}}.txt {params.file}.${{pop}}.vcf.gz | bcftools view -Oz -i 'MT=1' -o m1_{params.file}.${{pop}}.vcf.gz
            tabix -p vcf m1_{params.file}.${{pop}}.vcf.gz
        done
        bcftools merge -Oz --missing-to-ref -o {params.file}_all.vcf.gz m1_{params.file}.{{p31,p41,p51,p7}}.vcf.gz
        tabix -p vcf {params.file}_all.vcf.gz
        rm m1_{params.file}.{{p31,p41,p51,p7}}.vcf.gz*
        mkdir -p vcf_backup
        mv {params.file}.{{p31,p41,p51,p7}}.vcf.gz* vcf_backup
        mv samples.{{p31,p41,p51,p7}}.txt vcf_backup
        bcftools query -l {params.file}_all.vcf.gz | awk '{{print $1, $1}}' | awk '{{sub(/_i[0-9]*/, "", $2)}} 1' | sed 's/ /\\t/g' > {params.file}_all.samples.txt
        '''
    	


rule convert_to_plink:
    input:
        master_prefix + "_all.vcf.gz"
    output:
        temp(multiext(master_prefix,
    	".bed",".bim",".fam",".nosex"))
    params:
        file = master_prefix
    resources:
        mem = 8*1000,
        time = 30,
    threads: 8
    shell:
        '''
        module load any/plink/1.9
        plink --vcf {input} --make-bed --out {params.file}
        '''


rule run_king:
    input:
        multiext(master_prefix,
    	".bed",".bim",".fam",".nosex")
    output:
        kin = temp(multiext(master_prefix,
        ".kin",".kin0")),
        rel = master_prefix + ".relatives"
    params:
    	file = master_prefix
    resources:
        mem = 4*1000,
        time = 60,
    threads: 10
    shell:
        '''
        module load any/king/2.2.7
        king -b {params.file}.bed --kinship --prefix {params.file}
	less {params.file}.kin | awk -v OFS='\\t'  '$9 >= 0.0442 {{print $1 "_" $2, $1 "_" $3, $9}}' | tail -n +2 > {output.rel}
        '''


rule filter_relatives:
    input: 
    	master_prefix + ".relatives",
    	master_prefix + "_all.samples.txt"
    output:
        master_prefix + ".no_rel.samples"
    params:
        file = master_prefix
    resources:
        mem = 1000,
        time = 30
    shell:
        '''
        wd=$(pwd)
        python ${{wd}}/scripts/filter_relatives.py -p {params.file}
        '''



rule filter_vcf:
    input:
        vcf = master_prefix + "_all.vcf.gz",
        keep = master_prefix + ".no_rel.samples"
    output:
        vcf = master_prefix + ".no_rel.vcf.gz",
        index = master_prefix + ".no_rel.vcf.gz.tbi"
    resources:
        mem = 500,
        time = 30
    shell:
        '''
        module load bcftools/1.9
        module load tabix
        cut -f1 {input.keep} > {input.keep}.tmp
        bcftools view -Oz -S {input.keep}.tmp -o {output.vcf} {input.vcf}
        rm {input.keep}.tmp
        tabix -p vcf {output.vcf}
        '''


# this defines the regions to be used in covMA calculation. This command should be modified if we want to include only regions meeting a given condition (will do that later on by reading the file with QTL betas and dafs in )
rule create_bed:
	input: 
	    multiext(master_prefix, ".QTLs.freq.tsv",".QTLs.list.tsv")
	output: 
	    multiext(master_prefix + ".h2_" + str(config["h2_cutoff"]), ".bed", ".pos")
	params:
	    h2 = config["h2_cutoff"],
	    file = master_prefix
	shell:
	    '''
	    wd=$(pwd)
	    module load any/R/4.0.3
	    Rscript ${{wd}}/scripts/create_bed.R -l {wildcards.len_chr} -t {params.h2} -p {params.file}
	    '''




fstat_params={'cov_ma':'-a p7 -b p31 -b p41 -b p51 cov_ma'}


# this rule runs for about 2h on ~7K (relatives filtered) individuals, 1000 x 20K genome
rule fstat_single:
    input: 
        vcf = master_prefix + ".no_rel.vcf.gz", 
        bed = master_prefix + ".h2_{h2_cutoff}.bed", 
        sets = master_prefix + ".no_rel.samples"
    output: 
        master_prefix + ".h2_{h2_cutoff}.{stat,f3|f4|fst|cov_ma|f3M}"
    params: 
        opts = lambda w: fstat_params[w['stat']]
    threads: 10
    resources: 
        mem=4072, 
        time= lambda w: 1800
    shell: 
        '''
        wd=$(pwd)
        python ${{wd}}/scripts/Fstat_bed_vcf.py -p {threads} --single {params.opts} {input.vcf} {input.bed} {input.sets} > {output}
        '''

rule prepare_dosage_files:
    input:
        vcf = demography_prefix + "vcf_backup/" + file_prefix + ".p7.vcf.gz",
        pos = master_prefix + ".h2_{h2_cutoff}.pos",
        samples = master_prefix + ".no_rel.samples"
    output:
        expand(demography_prefix_for_expand + "gwas_hits/" + file_prefix_for_expand + ".m{m}.h2_{{h2_cutoff}}.dosage.gz", m = ["3", "4", "5"])
    params: 
        out_dir = demography_prefix + "gwas_hits",
        file_prefix = file_prefix
    resources: 
        mem = 1000, 
        time = 30
    shell: 
        '''
        module load bcftools/1.9
        mkdir -p {params.out_dir}
        less {input.samples} | grep p7 | cut -f1 | sed 's/p7_//g' | sort -V > {params.out_dir}/samples_no_rel.txt 
        for m in {{3..5}}
        do
        bcftools view -S {params.out_dir}/samples_no_rel.txt -R {input.pos} -i "MT = ${{m}}" {input.vcf} | 
        bcftools norm -m+ | bcftools +dosage | cut -f5- | tail -n +2 | sed 's/\.0//g' | gzip > {params.out_dir}/{params.file_prefix}.m${{m}}.h2_{wildcards.h2_cutoff}.dosage.gz
        done
        '''
        
        

# this rule runs for about 2h on ~7K (relatives filtered) individuals, 1000 x 20K genome
rule combine_covma_local_anc:
    input: 
        covma = master_prefix + ".h2_{h2_cutoff}.cov_ma", 
        asPGS = master_prefix + ".asPGS.tsv", 
        dosage = expand(demography_prefix_for_expand + "gwas_hits/" + file_prefix_for_expand + ".m{m}.h2_{{h2_cutoff}}.dosage.gz", m = ["3", "4", "5"]),
        sets = master_prefix + ".no_rel.samples"
    output: 
        multiext(master_prefix + ".h2_{h2_cutoff}", ".pdf", ".asPGS+covma.csv")
    params:
        deviating_pop = lambda wildcards: variable_parameters.loc[wildcards.setup, "deviating_pop"],
        opt_shifted = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_dev_pop"],
        opt_UK = lambda wildcards: variable_parameters.loc[wildcards.setup, "opt_of_UK"],
        w_shifted = lambda wildcards: variable_parameters.loc[wildcards.setup, "fitness_w_dev_pop"],
        w_UK = lambda wildcards: variable_parameters.loc[wildcards.setup, "fitness_w_UK"],
        path = demography_prefix,
        file = file_prefix
    resources: 
        mem = 1000, 
        time = 30
    shell: 
        '''
        wd=$(pwd)
        module load any/R/4.0.3
        Rscript ${{wd}}/scripts/covma_vs_la_server.R \
        -s {wildcards.setup} \
        -w {wildcards.fitness_w} \
        -d {params.deviating_pop} \
        -o {params.opt_shifted} \
        --w_pop {params.w_shifted} \
        --opt_uk {params.opt_UK} \
        --w_uk {params.w_UK} \
        -t {wildcards.h2_cutoff} \
        -f {params.file} \
        -p {params.path}
        '''


