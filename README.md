# SLiM_quantitative_trait
## SLiM simulation for testing out the behaviour of the Marnetto's Cov_MA

This is a hybrid simulation with the neutral burn-in run in msprime and then the last part with potential selection on the olygenic phenotype run in SLiM.

The msprime part requires the following python libraries: msprime, sys, pyslim, tskit, argparse.

Running the simulations requires specifying various parameters in the `config.yaml` and `variable_parameters.csv` files.

More general parameters are described in `config.yaml`:

```
initial_Ne: 14000 # Ne for the burn-in phases
n_chr: 1000 # number of cromosomes
len_chr: 20000 # length of an individual chromosome
mu: 1.25e-8 # mutation rate per base per generation
rec: 1e-8 # recombination rate per base per generation
trait_opt: 0.0 # trait optimum value
h2_cutoff: 0.00015 # per-SNP heritability used to assertain gwas hits
strong_frac: 0 # fraction of causal loci taken from a wider effect size distribution
frac_h2: 0 # fraction of the trait heritabiilty explained by causal loci coming from a wider effect size distribution
beta2: 1.5 # variance of the distribution the snp effects are drawn from
env_var: 0.5 # variance of the distribution the environmental noise is drawn from
```

The above combination of `strong_frac`, `frac_h2`, `beta2` and `env_var` correspond to a signle distribution to draw causal effects from and heritability around 0.5. 

The second file allows to launch multiple simulations at once by specifying the following paramenters for each simulation run.

`setup	msprime_seed	slim_burnin_seed	slim_demography_seed	deviating_pop	opt_of_dev_pop	opt_of_p6	opt_of_UK	w	fitness_w_dev_pop	fitness_w_UK`


