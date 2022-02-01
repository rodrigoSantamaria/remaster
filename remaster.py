# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 10:38:53 2017

TODO: Add import argparse and the proper arg parsing thing
@author: rodri
"""
#%%
import sys

import frecuencias
import bioio

import numpy as np
import re
import os


bedPath=None
revbedPath=None
fastaPath=None
ifreqPath=None
outPath=None
ofreqPath=None
genomePath=None
fastaSeqsPath=None
fastaRSesqsPath=None
imgFolder=None
outputCSV=None

fpombeSR={}

seqs=None

log2=False
nonorm=False
offset=0
extendLeft=0
extendRight=0
    
smoothed=True
window=20
step=1
klist=[2,1,3]
vline=[]
#%%
#isize=150
offset=0
genomePath="frecuencias/S_pombe_mio.fasta"
bedPath=None
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\nremaster version 0.2')
#    print('For help information for each function, try:\npython remaster.py <function> -h'
    print('\nUsage: ')
#    print('\tpython remaster.py -bed|-fasta|-freqs file [-seqs file] [-ofreqs file] [-rseqs file] [-img folder] [...]'
    print('\tpython remaster.py -bed|-fasta|-freqs file [-seqs file] -out folder [options]')
    print('\n\t--- Input arguments ---')
    print('\t-bed file:\tBed file with same length intervals to compute their corresponding sequences\' frequencies.')
    print('\t\t entries with a strand (5th) column as \'-\' will be taken as inverse complementary sequences.')
    print('\t-genome file:\tFasta file with genome sequences related to the BED file. (only required for -bed)')
    print('\t-fasta file:\tFasta file with same length sequences to compute their frequencies (substitutes -bed).')
    print('\t-freqs file:\tTSV file generated with this program (no bed or fasta file required, skips frequencies computation and uses these ones).')
    print('\t-seqs file:\tFasta file with sequences to remaster.')
    print('')
    print('\t--- Ouput ---')
    print('\t-out folder:\tName of the folder were output will be saved (if it does not exist it will be created.)')
    print('\t\t\tContains:')
    print('\t\t\t1) A fasta file with the remastered sequences corresponding to -seqs sequences')
    print('\t\t\t2) Tab separated (.tsv) files with computed frequencies for nucleotides (suffix -1) and codons (-3).')
    print('\t\t\t3) Figure files (.tiff) with computed frequencies for nucleotides, AT content and codons.')
#    print('\t\t\t-- Tab files are not generated i '
    
#    print('\t-ofreqs file:\tName for the CSV file with computed frequencies.'
#    print('\t-img folder:\tName of the folder to store frequency-related figures.'
#    print('\t-rseqs file:\tName for the fasta file with remastered sequences.'
    print('')
    print('\t--- Additional options ---')
    print('\t-nonorm:\tComputes percentage frequencies without normalizing by genome frequencies.')
    print('\t-log2:\t\tComputes log2 normalization of frequencies (only if -genome is set).')
    print('\t-k klist:\t\tComputes frequencies for the kmers specified in a comma separated integer list. The first number is the kmer whose freqs are used for remastering  (default 2,1,3).')
    print('')
    print('\t-nosmooth:\tDoes not average frequencies by window (see -window and -step).')
    print('\t-vline:\tDraws vertical lines on plots at the given positions (list separated by commas, default none).')
    print('\t-window n:\t Number of nucleotides to use for averaging (smoothing) frequencies (default 20).')
    print('\t-step n: \tNucleotide step for window scrolling when averaging (smoothing) frequencies (default 1).')
    print('')
  
    print('\t-offset n:\tOffset that must be applied to sequence frequencies for remastering (default 0)')
    print('\t-xl n:\t\tNucleotide extension to the left from bed starts (only for -bed, default 0).')
    print('\t-xr n:\t\ttNucleotide extension to the left from bed starts (only for -bed, default 0).')
    print('')
    
    
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
#    print("ARGV", sys.argv)
    
    if len(sys.argv)==1:
        print('\t\tremaster version 0.2')
        print('\t\tFor a list of functions in remaster, please try: python remaster.py -h')
        sys.exit()

    #Checking for required parameters
    if ('-h' in sys.argv)==False:
        if ('-out' in sys.argv)==False:
            print("ERROR: an output folder must be provided with -out (use -h for help).")
            sys.exit()
        if ('-bed' in sys.argv)==False and ('-fasta' in sys.argv)==False and ('-freqs' in sys.argv)==False:
            print("ERROR: either genomic positions (-bed) or sequences (-fasta) for frequency computation are required; or must be provided with -freqs (use -h for help).")
            sys.exit()
            
            
        
    #input retrieval
    for i in range(1, len(sys.argv)):
        
        if sys.argv[i]=='-h' or len(sys.argv)==1:
            printHelp()
            sys.exit()

        #input params
        if sys.argv[i]=="-bed": bedPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-revbed": revbedPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-genome": genomePath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-fasta": fastaPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-freqs": ifreqPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-seqs": 
            fastaSeqsPath=getValue(sys.argv[i], sys.argv[i+1])
            
            
        if sys.argv[i]=="-out":
            outPath=getValue(sys.argv[i], sys.argv[i+1])
            if os.path.isdir(outPath)==False:
                os.makedirs(outPath)
        

        #other configurations
        if sys.argv[i]=="-log2":
            log2=True
        if sys.argv[i]=="-nonorm":
            nonorm=True
        if sys.argv[i]=="-k":
            klist=getValue(sys.argv[i], sys.argv[i+1])
            klist=[int(x) for x in klist.split(",")]
        if sys.argv[i]=="-vline":
            vline=getValue(sys.argv[i], sys.argv[i+1])
            vline=[int(x) for x in vline.split(",")]
        if sys.argv[i]=="-xl":
            extendLeft=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-xr":
            extendRight=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-offset":
            offset=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-window":
            window=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-step":
            step=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-nosmooth":
            smoothed=False
        

    #setting up output file names
    imgFolder=outPath
    if bedPath!=None:
        ofreqPath=outPath+"/"+re.sub("\\..*", "", re.sub(".*/", "",bedPath))+"-freqs"
    elif fastaPath!=None:
        ofreqPath=outPath+"/"+re.sub("\\..*", "", re.sub(".*/", "",fastaPath))+"-freqs"
    
    if fastaPath!=None and bedPath!=None:
        print("Error: only one parameter, -bed or -fasta, should be used as input for frequencies (use -h for more help)")
        sys.exit()
        
    if fastaSeqsPath!=None:
        fastaRSeqsPath=outPath+"/"+re.sub("\\..*", "", re.sub(".*/", "",fastaSeqsPath))+"-rem.fasta"
            
    #processing
        
    # Read  nucleosome sequences and pombe genome, and compute frequencies
    if(ifreqPath==None):
        fpombeSR=frecuencias.pipelineFreq(positionFile=bedPath, seqFile=fastaPath,
                     window=window, step=step, vline=vline,
                     extendRight=extendRight, extendLeft=extendLeft,
                     genomeFile=genomePath,
                     positionRevFile=revbedPath,
                     outFile=ofreqPath,
                     imgFolder=imgFolder,
                     nonorm=nonorm, kmers=klist,
                     smoothed=smoothed, log2=log2
                     )
    else:
        print("Frequencies provided by user")
        fpombeSR["freq"]={}
        fr=frecuencias.readFreqs(ifreqPath)
        klist=[len(list(fr.keys())[0])]
        fpombeSR["freq"][klist[0]]=fr
        if(len(list(fr.keys())[0])==2):
            frecuencias.drawFreqs(fr, list(range(len(list(fr.values())[0]))), label="Frequency", imgFolder=imgFolder)#frecuencias.drawNucleotides(fr,fr, label="frequency", vline=[], title="", imgFolder=imgFolder) 

            
                 
    #%% Remasterize
    #if(fastaSeqsPath!=None and (3 in fpombeSR["freq"].keys())):
    if(fastaSeqsPath!=None):
        
        print("Reading sequences to remaster...")
        seqs=bioio.readFasta(fastaSeqsPath)
        
        print("Original seqs", [len(x) for x in seqs.values()])
        print("Remastering...")

        if(smoothed==False or ifreqPath!=None):
            freq=fpombeSR["freq"][klist[0]] #make it flexible
        else:
            freq=fpombeSR["smoothedFreq"][klist[0]]["freq"] #make it flexible
        
        rem0=frecuencias.remaster(seqs,freq,offset)
        #freq3=fpombeSR["freq"][2]
        #rem0=remaster2(seqs,freq3,offset) #according to the information provided, there's an ATG in nucleosome 1 at pos 59
        if(fastaRSeqsPath!=None):
            print("Writing to",fastaRSeqsPath)
            bioio.writeFasta(rem0, fastaRSeqsPath)
        else:
            print("The remastered sequences are:")
            print(rem0)
   
#%%
