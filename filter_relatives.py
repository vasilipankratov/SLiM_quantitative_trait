#!/usr/bin/env python

import pandas as pd
import argparse
#df is the pandas dataframe containing the data to be filtered
#kinshipfile is the filename of the file containing idA, idB and relatednes
#idcolname is a string, the column name in df where to filter, 
#minkin is a float, the minimal relatedness to consider 2 ids as relatives

parser = argparse.ArgumentParser(description="")
parser.add_argument('-p', '--prefix', type=str, help='define input file prefix', required=True)
args = parser.parse_args()
prefix = args.prefix

samples = pd.read_table(prefix + "_all.samples.txt", header=None)
samples.columns = ['id', 'pop']


def pyrelout(df, kinshipfile, idcolname, minkin=0.0442):
        mysamples=set(df.loc[:,idcolname])
        with open(kinshipfile) as f:
            kin=[x[:2] for x in [line.strip().split('\t') for line in f.readlines()] if float(x[2])>minkin and x[0] in mysamples and x[1] in mysamples and x [0]!=x[1]]
        kind={}
        for x in kin:
                for i in range(2):
                        try:
                                kind[x[i]].add(x[1-i])

                        except KeyError:
                                kind[x[i]]=set([x[1-i]])
              
        keep=[]
        kindkeyset=set(kind.keys())
        df['nrel']=[len(kind[myid]) if myid in kindkeyset else 0 for myid in df[idcolname]]
        df.sort_values('nrel',inplace=True)
        df.drop(['nrel'],axis=1,inplace=True)
        kins=set()
        for i,r in df.iterrows():
                if r[idcolname] in kins:
                        kins.remove(r[idcolname])
                else:
                        keep.append(i)
                        try:
                                kins=kins.union(kind[r[idcolname]])
                        except KeyError:
                                pass
        return(df.loc[keep])

    
keep_rel = pyrelout(samples, prefix + ".relatives", "id")
keep_rel.sort_values('id', inplace = True)
keep_rel.to_csv(prefix + '.no_rel.samples', index=False, sep = "\t", header=False)



