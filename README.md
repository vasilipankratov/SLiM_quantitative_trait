# SLiM_quantitative_trait
## SLiM simulation for testing out the behaviour of the Marnetto's Cov_MA

This is a hybrid simulation with the neutral burn-in run in msprime and then the last part with potential selection on the olygenic phenotype run in SLiM.

For now for illustration purposes I have a .sh file that combines the two part. Later on I'll do a snakemake pipeline for it.
The msprime part requires the following python libraries: msprime, sys, pyslim, tskit, argparse. To run this on rocket one needs a python environment with all those packages.


