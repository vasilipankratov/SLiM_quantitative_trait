#!/usr/bin/python
import msprime, sys, pyslim, tskit, argparse
from math import log
from math import exp

# currently I run it like
# python3 neutral_burn-in.py 10000 20000 5 999
# passing the script, Ne, chromosome length, number of chromosome and the random seed
# I'll perhaps us

parser = argparse.ArgumentParser(description="")
parser.add_argument('-N', '--Ne', type=int, help='effective pop size', required=True)
parser.add_argument('-n', '--nchr', type=int, help='number of chromosomes to simulate', required=True)
parser.add_argument('-l', '--lchr', type=int, help='chromsome length in bp', required=True)
parser.add_argument('-s', '--seed', type=int, help='random seed', required=True)
args = parser.parse_args()

Ne = args.Ne
l = args.lchr
n = args.nchr
seed = args.seed

print("Ne = " + str(Ne))
print("chromosome lenght = " + str(l))
print("number of chromosomes = " + str(n))
print("random seed = " + str(seed))

mu=1.25e-8 # mutation rate per bp

#output file name
outfile = "neutralBI_Ne_" + str(Ne) + "_length_" + str(n) + "x" + str(l) + "_seed_" + str(seed)

#creating a vector of breaks for the recombination map 
# for a 20kb chromosome it is 0, 19999, 20000, 39999, 40000, ...
a = list(range(l-1, l*n, l))
b = list(range(0, l*n+1, l))
breaks = a + b

breaks.sort()

# and a recombination rate vector
rate = [1e-8, 0.5]
rate = rate*n

#creating the rec map
recomb_map = msprime.RateMap(
  position = breaks,
  rate = rate)

#initializing a demographic model and setting the effective population size
demog_model = msprime.Demography()
demog_model.add_population(initial_size=Ne)

#running the simulation 
# I use the discrete time wright-fisher model because the normal hudson coalescence shouldn't be used due to sample size == Ne
ts = msprime.sim_ancestry(
        samples=Ne,
        demography=demog_model,
        random_seed=seed,
        recombination_rate=recomb_map,
        model = "dtwf")

# annotating the obtained tree sequence with SLiM metadata
mutated = pyslim.annotate(ts, model_type="WF", tick=1, stage="late")

# adding neutral mutations on top
mutated = msprime.sim_mutations(mutated, rate = mu, random_seed = seed,
                                model = msprime.SLiMMutationModel(type = 1))
# simplifying the tree
mutated.simplify()

# writing output
mutated.dump(outfile)


