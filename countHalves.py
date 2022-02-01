#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 10:55:55 2021

@author: rodri
"""
#%%
######################### READ FASTA ####################################
"""
path - path to the fasta file
asDict - if false it returns a list that respects fasta line order
returns: a dictionary with the description of sequences as keys 
           (removing anything in parenthesis and the '>' symbol)
           and the sequences as values (str), in uppercase
requires: re for text substitution
"""
def readFasta(path=""):
    import re
    f=open(path)
    reader=f.readlines()
    f=open(path)
    ff={} 
    k=re.sub("\(.*\)", "", re.sub(">", "", reader[0])).strip()
    i=0
    seq=""
    while i<len(reader):
        while(i<len(reader)):
            line=next(f)
            i+=1
            if(line.startswith(">")):
                ff[k]=seq
                seq=""
                k=re.sub("\(.*\)", "", re.sub(">", "", line)).strip()
                break
            else:
                seq+=line.replace("\n","").upper()
    if(len(seq)==0):
        seq+=line.replace("\n","").upper()
    ff[k]=seq
    return ff

#%% #################### COUNT HALVES ############################
"""
fastaPath - path to a fasta file with one or more nucleosome sequences
k - list of k-mers to investigate. By default [1,2] (single nucleotides and dinucleotides)
section - position by which to split halves (default 75)
outPath - folder where files halvesCount-[k]mer.csv will be saved
"""

def countHalves(fastaPath, outPath, k=[1,2], section=75):
    import numpy as np
    import os
    import time
    import itertools
    seqs=readFasta(fastaPath)
    
    #0) Determine all the kmers
    for kk in k:
        print("-------------- "+str(kk)+" ----------------------")
        keys=["".join(x) for x in list(itertools.product(["A","C","G","T"],repeat=kk))]
      
        #1) Compute frequencies
        #t0=time.clock()
        fr={}
        for key in keys:
            fr[key]=np.zeros(2*kk)
        if(not os.path.isdir(outPath)):
            os.mkdir(outPath)
        f=open(outPath+"halvesCount-"+str(kk)+"mer.cvs", "w")
        f.write(str(kk)+"mer")
        for frame in range(0,kk):
            print("------------------frame"+str(frame)+"-----------------")
            cont=frame*2+0
            print(str(cont)+"\t"+str(frame)+"\t"+str(kk))
            for seq in seqs.values():       #First half
                for i in range(frame,section,kk):
                   kmer=seq[i:i+kk]
                   if(("N" in kmer)==False):
                       fr[kmer][cont]+=1
            cont=frame*2+1
            print(cont)
            for seq in seqs.values():       #Second half
                for i in range(section+frame,len(seq)-kk+1,kk):
                   kmer=seq[i:i+kk]
                   #if(kmer=="AA"):
                   #    print(seq[i-4:i+4])
                   if(("N" in kmer)==False):
                       fr[kmer][cont]+=1
            if(kk>1):
                f.write("\tframe"+str(frame)+"\t\t")
        if(kk>1):
           f.write("total\n")
        for i in range(kk):
            f.write("\tsection1\tsection2")
        if(kk>1):
            f.write("\tsection1\tsection2")
        f.write("\ttotalDiff\n")
        for kmer in fr.keys():
            f.write(kmer+"\t")
            for i in range(len(fr[kmer])):
                f.write(str(fr[kmer][i])+"\t")
            if(kk>1):
                sum1=sum(fr[kmer][range(0,len(fr[kmer]),2)])
                sum2=sum(fr[kmer][range(1,len(fr[kmer]),2)])
                f.write(str(sum1)+"\t"+str(sum2)+"\t"+str(sum1-sum2)+"\t")
            else:
                f.write(str(fr[kmer][-2]-fr[kmer][-1]))
            f.write("\n")
        f.close()

#%% RUN COUNT HALVES ######################################

countHalves(fastaPath="/home/rodri/data/remaster/plusOne/His7nucs.fasta",
            k=[1,2], section=75,
            outPath="/home/rodri/data/remaster/halves/his7/")

countHalves(fastaPath="/home/rodri/data/remaster/halves/ura4_pombe_WTnuc2-5.fasta",
            k=[1,2], section=75,
            outPath="/home/rodri/data/remaster/halves/ura4_pombe_WTnuc2-5/")
