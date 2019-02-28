#!/usr/local/bin/python

import os, sys, glob, subprocess
import numpy as np
#import cyvcf2
import pysam
from pysam import VariantFile
#import scikit-allel
import allel


class masking(object):
    Bampath=""
    variantfile=None
    Vcfpath=""
    Variants=[]
    

    def __init__(self,bampath,vcfpath):
        ### populates class variables 
        ### input : string(bampath), string(vcfpath)
        ### output : none
        self.Bampath = bampath
        self.Bamfile = pysam.AlignmentFile(bampath, "rb")
        self.Vcfpath = vcfpath
        self.variantfile = VariantFile(self.Vcfpath)
        self.populateVarlist()
        

    def populateVarlist(self):
        ### populate varaints from vcf file 
        ### input : none
        ### output : none
        for rec in self.variantfile.fetch():
            self.Variants.append(rec)

    def printMaskVars(self):
        ### prints class variables 
        ### input : none
        ### output : none 
        for v in self.Variants:
            print (v.chrom, v.pos, v.ref, v.alts)
        
        print ("Bampath :", self.Bampath)
        print ("Vcfpath:" , self.Vcfpath)

    def maskVariant(self,varRec,bamAlign):
        ### masking the alt base to ref
        ### input : VariantRecord(varRec), AlignedSegment
        ### output : AlignedSegment (modified)
        
        if not self.doesOverlap(varRec,bamAlign):
            return bamAlign
        else:
            print ("found the overlap with variant")
            AlIndex=varRec.pos - (bamAlign.reference_end - bamAlign.reference_length + 1)

            queryBases=bamAlign.query_sequence
            if queryBases[AlIndex] == varRec.ref :
                print (" sequences matches with ref ...Wooohoooo!!")
                return bamAlign
            elif queryBases[AlIndex] == varRec.alts[0] :
                queryBases == self.replaceChar(queryBases,varRec.ref,AlIndex)
                print ("unmasked",bamAlign.query_sequence)
                bamAlign.query_sequence = queryBases
                #print ("ismasked",bamAlign.query_sequence)
                print ("ismasked", queryBases)
                return bamAlign
            else :
                print ("all hell broken, get a break ")
                return bamAlign 

    def replaceChar(self,bamSeq,ref,index):
        bamSeq=bamSeq[:index] + ref + bamSeq[index+len(ref):]
        return bamSeq

    def doesOverlap(self,varRec,bamAlign):
        ### check if vcf overlaps bamalignment
        ### input : VariantRecord(varRec), AlignedSegment
        ### output : boolean
        pos=varRec.pos
        bamEndPos=bamAlign.reference_end
        bamStartPos= bamEndPos - bamAlign.reference_length
        
        return (pos >= bamStartPos and pos <= bamEndPos)

    def maskAllVariants (self):
        iter= self.Bamfile.fetch()
        for x in iter:
            for v in self.Variants :
                read = self.maskVariant (v,x)
        return 
def main():
    M=masking(sys.argv[1],sys.argv[2])
    M.printMaskVars()
    M.maskAllVariants()
    
if(__name__ == "__main__"):
    main()






