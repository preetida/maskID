#!/usr/local/bin/python

import os, sys, glob, subprocess
import numpy as np
#import cyvcf2
import pysam
from pysam import VariantFile
#import scikit-allel
#import allel


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
            return (bamAlign,False)
        elif len(varRec.alts) != 1 :
            return (bamAlign, False)
        else:
            print ("found the overlap with variant")
            AlIndex=(varRec.pos - 1) - (bamAlign.reference_end - bamAlign.reference_length + 1)  
            queryBases=bamAlign.query_sequence
            if queryBases[AlIndex] == varRec.ref :
                return (bamAlign,False)
            elif queryBases[AlIndex] == varRec.alts[0] :
                queryBases = self.replaceChar(queryBases,varRec.ref,AlIndex)
                print ("queryBas",queryBases)
                bamAlign.query_sequence = queryBases
                return (bamAlign,True)
            else :
                print ("Unhandle case for maskvariant")
                return (bamAlign,False) 

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
        masked_bam=pysam.AlignmentFile("masked.norm.bam", "wb", template=self.Bamfile)
        for x in iter:
            read = x
            for v in self.Variants :
                ret = self.maskVariant (v,x)
                b = ret[1]
                if b :
                    read = ret[0]
            masked_bam.write(read)
        masked_bam.close()
        return 

def main():
    M=masking(sys.argv[1],sys.argv[2])
    M.printMaskVars()
    M.maskAllVariants()
    
if(__name__ == "__main__"):
    main()






