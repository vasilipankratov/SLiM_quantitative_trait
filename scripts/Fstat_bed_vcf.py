#!/usr/bin/env python

from re import split, search, match,sub
from sys import stderr,stdout,stdin,exit,argv
from os.path import isfile
from itertools import product
import gzip
from os import stat
import argparse
import pysam
from multiprocessing import Pool
#import warnings


def safe_open(filename):
	if search('\.gz$', filename):
		f=gzip.open(filename,'rt')
	else:
		f=stdin if filename=="-" else open(filename,'r')
	return(f)

def floatOrNone(x):
	try:
		return(float(x))
	except ValueError:
		return(None)

def get1stvalue(x):
	return(list(x.values())[0])

####### defined in Patterson et al 2012 appendix A

def h(f,n):
	c = f*n
	h = (c*(n-c)) / (n*(n-1))
	return(h)

#computes frequency of a set from a list of indexes and a list containing fields of a ff line
def freq_comp(ids,vcfl):
	x=[hap for i in ids for hap in vcfl.samples[i]['GT'] if hap is not None]
	try:
		return(sum(x)/len(x))
	except ZeroDivisionError: #when size of mypop is 0
		return(None)

def size_comp(ids,vcfl):
	x=[1 for i in ids for hap in vcfl.samples[i]['GT'] if hap is not None]
	return(len(x))

#fstat computation in a single chr independently
def F_collector(filename,bed,samples,args):
	global F
	assert isfile(filename+'.tbi') or isfile(filename+'.csi'), filename+' has no tabix index'
	vcf=pysam.VariantFile(filename)

	joint_sets=[x for x in args.a_pop+args.b_pop+[args.c_pop]+[args.d_pop] if x is not None and search(r'\+',x)]
	for js in joint_sets:
		samples[js]=[sample for set_id in js.split('+') for sample in samples[set_id]]


	#a_ids contains ids for the individuals who belong to populations a
	if args.single: #fstat computation for each individual in a independently
		a_ids={i:[i] for x in args.a_pop for i in samples[x]}
	else:
		a_ids={x:samples[x] for x in args.a_pop}

	#b_ids contains ids for populations b 
	b_ids={x:samples[x] for x in args.b_pop}
	#c_ids contains ids for population c
	c_ids={args.c_pop:samples[args.c_pop]} if args.c_pop else None
	#d_ids either contains an index for population d or a list of indexes for the individuals who belong to it
	d_ids={args.d_pop:samples[args.d_pop]} if args.d_pop else None

	myfstat=FF_collector(args.stat,a_ids,b_ids,c_ids,d_ids)
	for bedline in bed:
		for vcfline in vcf.fetch(bedline[0], int(bedline[1]), int(bedline[2])):
			myfstat.addline(vcfline)
	return(myfstat.tot) #tot is a dictionary with key= a tuple (a,b) and value=  a list with [tot numerator, tot denominator, snp count, a size, d name]


class FF_collector:
	def __init__(self,f,a,b,c,d):
		available_stats={	'f3':(self.f3,self.prepf3,self.formulaf3),
					'f4':(self.f4,self.prepf4,self.formulaf4),
					'fst':(self.fst,self.prepfst,self.formulafst),
					'cov_ma':(self.cov_ma,self.prepcov_ma,self.formulacov_ma)}
		self.F,self.prep,self.formula=available_stats[f]
		self.a_ids=a
		self.b_ids=b
		self.c_ids=c
		self.d_ids=d
		self.tot={}

	def addline(self,vcfl):
		#this prepares other stuff
		self.prep(vcfl)
		#this computes fstat for all combinations
		for x,y in product(self.a_ids.keys(), self.b_ids.keys()):
			try:
				fstat=self.F(x,y)
			except (ZeroDivisionError, TypeError):
				#warnings.warn(self.formula(x,y)+' could not be computed at position '+str(vcfl.chrom)+'\t'+str(vcfl.pos),stacklevel=2)
				#print('[Warning] '+self.formula(x,y)+' could not be computed at position '+str(vcfl.chrom)+'\t'+str(vcfl.pos),file=stderr)
				continue
			try:
				self.tot[(x,y)][0:3]=[x+y for x,y in zip(self.tot[(x,y)][0:3],fstat+[1])]
			except KeyError:
				self.tot[(x,y)]=fstat+[1,len(self.a_ids[x]),self.formula(x,y)]

	def computeab(self,vcfl):
		#this computes a
		a={k:freq_comp(v,vcfl) for k,v in self.a_ids.items()}
		#this computes b
		b={k:freq_comp(v,vcfl) for k,v in self.b_ids.items()}
		return([a,b])

	#fst(test,anc)
	def fst(s,x,y):
		num=(s.a[x]-s.b[y])**2 - h(s.a[x],s.n_a[x])/(s.n_a[x]) - h(s.b[y],s.n_b[y])/(s.n_b[y])
		den=num + h(s.a[x],s.n_a[x]) + h(s.b[y],s.n_b[y])
		return([num,den])

	def prepfst(self,vcfl):
		self.a,self.b=self.computeab(vcfl)
		self.n_a={k:size_comp(v,vcfl) for k,v in self.a_ids.items()}
		self.n_b={k:size_comp(v,vcfl) for k,v in self.b_ids.items()}

	def formulafst(self,x,y):
		return('FST('+x+','+y+')')

	#f3(test,anc;out)
	def f3(s,x,y):
		num = (s.c-s.a[x])*(s.c-s.b[y]) - h(s.c,s.n_c)/(s.n_c)
		den = (2*s.c*(1-s.c))
		return([num,den])

	def prepf3(self,vcfl):
		self.a,self.b=self.computeab(vcfl)
		ck=list(self.c_ids.keys())[0]
		self.c=freq_comp(self.c_ids[ck],vcfl)
		self.n_c=size_comp(self.c_ids[ck],vcfl)

	def formulaf3(self,x,y):
		ck=list(self.c_ids.keys())[0]
		return('F3('+x+','+y+';'+ck+')')

	#f4(test,mod;anc,out)
	def f4(s,x,y):
		#z=s.a_d_assoc[x]
		num=(s.a[x]-s.d)*(s.b[y]-s.c)
		den=1
		return([num,den])

	def prepf4(self,vcfl):
		self.a,self.b=self.computeab(vcfl)
		ck=list(self.c_ids.keys())[0]
		self.c=freq_comp(self.c_ids[ck],vcfl)
		dk=list(self.d_ids.keys())[0]
		self.d=freq_comp(self.d_ids[dk],vcfl)

	def formulaf4(self,x,y):
		ck=list(self.c_ids.keys())[0]
		dk=list(self.d_ids.keys())[0]
		return('F4('+x+','+dk+';'+y+','+ck+')')

	#covariance_modern_vs_ancient(mod,anc)
	def cov_ma(s,x,y):
		num=(s.a[x]-s.avg_a)*(s.b[y]-s.avg_b)
		den=1
		return([num,den])

	def prepcov_ma(self,vcfl):
		self.a,self.b=self.computeab(vcfl)
		this_a=[x for x in self.a.values() if x is not None]
		try:
			self.avg_a=sum(this_a)/len(this_a)
		except ZeroDivisionError:
			self.avg_a=None
		this_b=[x for x in self.b.values() if x is not None]
		try:
			self.avg_b=sum(this_b)/len(this_b)
		except ZeroDivisionError:
			self.avg_b=None

	def formulacov_ma(self,x,y):
		return('cov_MA('+x+','+y+')')

def main(arguments):

	parser = argparse.ArgumentParser(
		description='''
compute f-statistics according to Patterson 2012 appendix A on a bed of regions, for groups of individuals present in input VCF and defined by a .sets file.
f-statistics:
	f3(c;a,b) as (c-a)(c-b) corrected for sample size of c
	f4(a,d;b,c) as (a-d)(b-c)
	fst(a,b) as f2/(f2 - h_a - h_b) corrected for sample size of a and b

other statistics:
	cov_ma(a,b) as (a-avg(a))(b-avg(b))

More a and b populations are accepted as in -a <pop1> -a <pop2> -b <pop3> -b <pop4>..., to do many statistics in one run. Plus-separated lists of populations will be considered as one population, as in -a <pop1>+<pop2>. at least one a and one b population are required.
On the other hand only one c an one d populations are accepted and required for the statistics that need them.
''',
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('stat', help="statistic to compute among f3, f4, fst, cov_ma")
	parser.add_argument('vcffile', help="VCF file")
	parser.add_argument('bedfile', help="Candidate regions bed file")
	parser.add_argument('setsfile', help="Individual sets file")
	parser.add_argument('--single', help="compute the statistics separately for each individual in 'a', i.e. each individual is a set. When more a are defined each individual in any set is considered independently",action='store_true')
	parser.add_argument('-a','--a_pop', help="Populations to use as a in f3(c;a,b), f4(a,d;b,c), fst(a,b), cov_ma(a,b)",type=str,required=True,action='append')
	parser.add_argument('-b','--b_pop', help="Populations to use as b in f3(c;a,b), f4(a,d;b,c), fst(a,b), cov_ma(a,b)",type=str,required=True,action='append')
	parser.add_argument('-c','--c_pop', help="Reference population to use as c in f3(c;a,b), f4(a,d;b,c)",type=str,default=None)
	parser.add_argument('-d','--d_pop', help="Reference population to use as d in f4(a,d;b,c), if -d is not supplied, by default is composed by all individuals in Sets file",type=str,default=None)
	#parser.add_argument('--rm_a_from_d', help="From the reference population to use as d in f4(a,d;b,c) are removed the individuals already present in a, so each test is 'a' against 'non-a' taken from the individuals in Sets file. only valid with default d",action='store_true')
	parser.add_argument('-p', help="N of parallel processes",type=int, default=1)
	args = parser.parse_args(arguments)
	
	if args.stat=='f3':
		assert args.c_pop, "population c needed to compute f3; supply it with -c"  
	if args.stat=='f4':
		assert args.c_pop, "population c needed to compute f4; supply it with -c"  
		assert args.d_pop, "population d needed to compute f4; supply it with -d"  
	#if args.rm_a_from_d:
	#	assert not args.d_pop, "samples of 'a' can be removed from 'd' only if 'd' are all samples in Sets file, remove option -d"
	#	if args.single:
	#		warn("option --rm_a_from_d is not applied with single mode")
	#		args.rm_a_from_d=False
	#for f in [args.setsfile,args.bedfile]:
	#	assert stat(f).st_size != 0, f+" is empty"

	#open bed with candidate regions
	wg_bed={}
	with safe_open(args.bedfile) as bed:
		for l in bed:
			if l.startswith('#'):
				continue
			i=l.strip().split("\t")
			try:
				wg_bed[i[0]].append(i)
			except KeyError:
				wg_bed[i[0]]=[i]
		
	#open file with trait classification
	samples={}
	with safe_open(args.setsfile) as tr:
	#	samples=dict([l.strip().split('\t') for l in tr.readlines()])
		for l in tr:
			if l.startswith('#'):
				continue
			i=l.strip().split("\t")
			try:
				samples[i[1]].append(i[0])
			except KeyError:
				samples[i[1]]=[i[0]]
	for pop in args.a_pop+args.b_pop+[args.c_pop,args.d_pop]:
		assert pop==None or len(samples[pop])>0 , 'No samples found in setsfile for set'+pop

	#run F statistic for each chromosome independently
	pool=Pool(processes=args.p) # start p worker processes
	vcfs=[sub('#CHR',x,args.vcffile) for x in wg_bed.keys()]
	samples_repeated=[samples for x in wg_bed.keys()]
	args_repeated=[args for x in wg_bed.keys()]
	F_inputs=list(zip(vcfs,wg_bed.values(),samples_repeated,args_repeated))
	results=pool.starmap(F_collector,F_inputs)
	chr_F={k:v for k,v in zip(wg_bed.keys(),results)}
	
	# check that individuals in sets are the same in different chrs removed

	#sum chromosome results
	tot_F={}
	for thischr in chr_F.values():
		for k,v in thischr.items():
			try:
				tot_F[k][0:3]=[x+y for x,y in zip(tot_F[k][0:3],v[0:3])]
			except KeyError:
				tot_F[k]=v
			else:
				assert tot_F[k][3]==v[3] and tot_F[k][4]==v[4]
	#average results and print output
	#d_pop='NA'
	#if args.stat=='f4':
	#	d_pop=args.d_pop
	#	if args.rm_a_from_d: 
	#		d_pop='all_non_a'
	for k,v in tot_F.items():
		#v[0] is the numerator of Fstatistic
		#v[1] is the denominator of Fstatistic
		#v[2] is the snp count
		try:
			Fstatistic=format((v[0]/v[2])/(v[1]/v[2]),'.6g') #fraction of averages
		except ZeroDivisionError:
			Fstatistic='NA'
		a=k[0]
		b=k[1]
		snpsize=str(v[2])
		setsize=str(v[3])
		formula=v[4]
		print("\t".join([a,b,formula,Fstatistic,snpsize,setsize]))

if __name__ == '__main__':
    exit(main(argv[1:]))

