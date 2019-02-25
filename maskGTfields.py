#!/usr/local/bin/python

import os, sys, glob, subprocess
import numpy as np
#import cyvcf2
import pysam
#import scikit-allel
import allel


class vcfInfo(object):
    Chr=""
    Pos=-1
    Ref=""
    Alt=""
    
    def __init__(self,vcfrow):
        self.Chr=vcfrow[0]
        self.Pos=int(vcfrow[1])
        self.Ref=vcfrow[2]
        self.Alt=vcfrow[3]

def main():
    callset=np.genfromtxt(sys.argv[1], delimiter='\t', dtype=None, skip_header=120,usecols=(0,1,3,4))
    variants=[]
    for row in callset:
        v=vcfInfo(row)
        #print (v.Chr,v.Pos,v.Ref,v.Alt)
        variants.append(v)


if(__name__ == "__main__"):
    main()






