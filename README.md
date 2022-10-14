# SLiM_quantitative_trait
## SLiM simulation for testing out the behaviour of the Marnetto's Cov_MA

This is a hybrid simulation with the neutral burn-in run in msprime and then the last part with potential selection on the olygenic phenotype run in SLiM.

For now for illustration purposes I have a .sh file that combines the two part. Later on I'll do a snakemake pipeline for it.
The msprime part requires the following python libraries: msprime, sys, pyslim, tskit, argparse. To run this on rocket one needs a python environment with all those packages.


## vcf files for Massimo
The file m1_p0+p7+p31+p41+p51+p61.vcf.gz contains genotypes for individuals from 6 populations:

p0 - East Asia

p7 - UKBB

p31 - Anatolian Neolithic farmers (a sister group of the admixing ones)

p41 - WHG (a sister group of the admixing ones)

p51 - Yamnaya (a sister group of the admixing ones)

p61 - European Early farmers (a sister group of the admixing ones)

Note: in the simulation itself the admixing populations are labeled p3, p4, p5 etc while the sampled sister group is p31, p41, p51 etc

This vcf file was produced with the simulation code presented here, then I've renamed the samples to follow a p{p}i{i} pattern where p is the population index and j is the individual index (0 to 99 as there are 100 samples for each population) and filtered out only mutations of type m1 (well actually, this means all mutations but avoiding multiallelics in the resulting vcf).

I need to do a couple more checks that the simulation worked as intended but I will have time for this only at the end of the week. However, if any bugs exist, I expect them only to affect the phenotype, not the demography.
