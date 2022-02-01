#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:22:29 2018

@author: rodri
"""

#%%
#NOTE: in python 3 dtype such as S40 converts to b"", not directly to string,
#so you must either to a .decode when decoding or declare dtype as str.
def searchGene(text0, dataGFF, types=["gene"], exact=True, returnVerbose=False):
    import numpy as np
    #import string
    import re
    data=dataGFF
    result=[]
    text=text0.lower()
    #print ("Search ",text," in ", types)
    if(types[0]!="any"):
        wanted_set = set(types)  #Much faster look up than with lists, for larger lists
        @np.vectorize
        def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
        data=dataGFF[selected(dataGFF["type"])]
    if exact:
        for x in data:
            #PY3: string does not have find method now
            #if(string.find(x["id"].lower(),text)>=0 or string.find(x["name"].lower(),text)>=0):
            if(x["id"].lower()==text or x["name"].lower()==text):
              if(returnVerbose):
                  result.append({"start":x["start"], "end":x["end"], "name":x["name"], "chromosome":x["chromosome"], "sense":x["sense"]})
              else:
                  result.append({"start":x["start"], "end":x["end"]})
    else:
        for x in data:
            #PY3 change: re now returns None if nothing found
            #if(re.search(text,str(x["id"]).lower())>=0 or re.search(text,str(x["name"]).lower())>=0):
            if(re.search(text,str(x["id"]).lower())!=None or re.search(text,str(x["name"]).lower())!=None):
              if(returnVerbose):
                  result.append(x)
              else:
                  result.append({"start":x["start"], "end":x["end"]})
        
    return result
#searchGene("cos", gffData, exact=False, returnVerbose=True)    
#%%
def gene2seq(gene, gffData, fastaData):
    import os
    os.chdir("/home/rodri/workspace/python/IBFG/remaster")
    import bioseq

    #1) read gtf and get positions
    gi=searchGene(gene, gffData, types=["gene"], exact=False, returnVerbose=True)
    print(len(gi), "genes found")
    
    
    #4) Associate sequences (fastaData and fastaData2 to the genes of interest (gi)
    seqs={}
    for gene in gi:
        #Get gene sequence
        s=gene["start"]
        e=gene["end"]
        n=gene["name"]
        gid=gene["id"]
        print("Gene: "+n+" "+gid)
        #ch="chromosome="+gene["chromosome"].replace("chr", "")+"]"
        ch=gene["chromosome"]
        at=0
        for k in fastaData.keys():
            if ch in k:
                at=k
                break
        if(gene["sense"]=="+"):
            seqs[n]=fastaData[at][s:e]
        else:
            fastaDataR=bioseq.complement(fastaData[at])
            print(n,len(fastaData[at]), len(fastaDataR),s,e,at)
            seqs[n]=fastaDataR[s:e]
        
    return seqs
    
#gene2seq("YDL039C", ann, genome)    
#%%
"""
 Get the sequences from genome corresponding to positions (starts) with length windowSize
starts - positions dict fom the genome
genome - dict as returned by bioio.readFasta
windowSize - length of the sequences to return (default 150)
comp - if true, it returns the reverse complement sequences (default False)
"""
def start2seq(starts, genome, windowSize=150, comp=False):
    seqs={}
    for track in starts.keys():
        positions=starts[track]
        seqs[track]=[]
        for i in range(len(positions)):
            start=int(positions[i])
            end=start+windowSize
            ss=genome[track][start:end].upper()
            if(len(ss)==windowSize): #to avoid adding small sequences
                if(comp):
                    ss=bioseq.complement(ss)
                seqs[track].append(ss)

    #convert to array (could have done it directly into a list, though)
    seqlist=[]
    for k in seqs.keys():
        for x in seqs[k]:
            seqlist.append(x)
    return seqlist
#%%
#Searches a given text in go descriptions, then searches for genes annotated
#with the corresponding go terms and returns their locations as start-end
def searchGO(text, dataGOA, dataGO, dataGFF):
    import numpy as np
    import string
    result=[]
    #print("Searching GO terms for: |",text,"|")
    for k in dataGO.keys():
      if(string.find(dataGO[k].lower(),text)>=0):
          result.append(k)
    result=set(result)
    res2=[]
    for x in dataGOA:
        if x["go_id"] in result:
            res2.append(x["gene_name"])
    res2=set(res2)
    
    wanted_set = set(["gene", "ORF"])  # Much faster look up than with lists, for larger lists
    @np.vectorize
    def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
    data=dataGFF[selected(dataGFF["type"])]

    res3=[]
    for x in data:
      if x["name"] in res2:
         res3.append({"start":x["start"], "end":x["end"]})
    return res3
    
#%% -------------------------------- SEQ FIND ----------------------------
"""
Given a set of large sequences (as fasta file) and a set of motifs as a file with
motif names in column 1 and contributing sequences, it returns a bed file with
the occurrences of the motifs on the fasta.
"""
def seqFind(fastaPath=None, motifPath=None, bedPath=None, inc=0):
   
    import bioio
    import bioseq
    import re
    
    #0) Read input
    seqs=bioio.readFasta(fastaPath)
    
    f=open(motifPath)
    motifs={}
    for line in f.readlines():
        line=line.strip().split("\t")
        mot=line[0]
        seq=line[1]
        if not mot in motifs.keys():
            motifs[mot]={}
  
        if not seq in motifs[mot].keys():
            motifs[mot][seq]=1
        else:
            motifs[mot][seq]+=1
            
    #1) Search for occurrences
    ocs=[]
    for k in seqs.keys():
        seq=seqs[k]
        #ndrk=k.replace("_rev", "")
        ndrk=k
        kk=ndrk.split(":")
        ch=kk[0]
        kk=kk[1].split("-")
        start=int(kk[0])
        
        #end=int(kk[1])
               
        for mot in motifs.keys():
            mers=motifs[mot].keys()
            for mer in mers:
                imer=bioseq.complement(mer)
                if(mer in seq):
                    strand="+"
                    s=[x.start() for x in re.finditer(mer, seq)]
                    for si in s:
                        #ocs.append({"ch":ch ,"start":(start+si+1), "end":(start+si+1+len(mer)), "strand":strand, "ndr":k.replace("_rev", ""), "motif":mot, "seq":mer}) 
                        ocs.append({"ch":ch ,"start":(start+si+inc), "end":(start+si+len(mer)-1+inc), "strand":strand, "ndr":k.replace("_rev", ""), "motif":mot, "seq":mer}) 
                if(imer in seq):
                    strand="-"
                    s=[x.start() for x in re.finditer(imer, seq)]
                    for si in s:
                        ocs.append({"ch":ch ,"start":(start+si-inc), "end":(start+si+len(mer)-1-inc), "strand":strand, "ndr":k.replace("_rev", ""), "motif":mot, "seq":mer}) 
                        #ocs.append({"ch":ch ,"start":(start+si+1), "end":(start+si+len(mer)+1), "strand":strand, "ndr":k.replace("_rev", ""), "motif":mot, "seq":mer}) 
                    
                
    fw=open(bedPath, "w")
    for oc in ocs:
        fw.write(oc["ch"]+"\t"+str(oc["start"])+"\t"+str(oc["end"])+"\t"+oc["motif"]+"\t"+oc["seq"]+"\t"+oc["strand"]+"\n")
    fw.close()       
    return ocs

#m=seqFind("/home/rodri/data/ndr/ndrSP/version5/NDRs_total_Sp.fasta", "/home/rodri/data/ndr/ndrSP/version5/motivosTest2.csv", "/home/rodri/data/ndr/ndrSP/version5/test.bed")
m=seqFind("/home/rodri/data/ndr/ndrSP/version5/test.fasta", "/home/rodri/data/ndr/ndrSP/version5/motivosTest2.csv", "/home/rodri/data/ndr/ndrSP/version5/testInc.bed", inc=1)

#%%
import sys
import re
import time

import bioio
import bioseq


bedPath=None
fastaPath=None
motifPath=None
inc=0

#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\nbiosearch version 0.2')
    print('\nUsage (python3): ')
    print('\tpython3 biosearch.py -fasta path -seqs path -out path')
    print('\n\t--- Input arguments ---')
    print('\t-fasta file:\tFasta file with nucleotide sequences. Sequence names must have pattern "chromosomeX:start-end".')
    print('\t-seqs file:\tTab delimited file with motifs (column 1) and sequences of the motif (column 2). For example:')
    print('\t\t\tm1\tATCG')
    print('\t\t\tm1\tATCT')
    print('\t\t\tm2\tTTTTGGGG')
    print('\t\t\tm2\tTTTTGCGG')
    print('\t\t\tm2\tATTTGGGC')
    print('\t-zero:\tIf positions at the fasta file have 0-start and output must be 1-start, use this option')
    print('\n\t--- Output arguments ---')
    print('\t-out file:\tBed file with genomic coordinates of occurrences of the motifs in the fasta sequences. Scores are set to the motif sequence found')
    print('')
  
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
    
    if len(sys.argv)==1:
        print('\t\tbiosearch version 0.2')
        print('\t\tFor a list of functions in biosearch, please try: python biosearch.py -h')
        sys.exit()

    #Checking for required parameters
    if ('-h' in sys.argv)==False:
        if ('-out' in sys.argv)==False:
            print("ERROR: an output folder must be provided with -out (use -h for help).")
            sys.exit()
        if ('-out' in sys.argv)==False or ('-fasta' in sys.argv)==False or ('-seqs' in sys.argv)==False:
            print("ERROR: All the parameters -fasta, -seqs and -motif must be defined and pointing to valid file paths.")
            sys.exit()
            
            
        
    #input retrieval
    for i in range(1, len(sys.argv)):
        
        if sys.argv[i]=='-h' or len(sys.argv)==1:
            printHelp()
            sys.exit()

        #input params
        if sys.argv[i]=="-out": outPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-fasta": fastaPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-seqs":  seqsPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-zero":  inc=1
            

    #processing
    print("Searching for sequences...")
    t0=time.clock()
    seqFind(fastaPath, seqsPath, outPath, increase=inc)
    print("... finished in ", round((time.clock()-t0)/60., 2), " minutes")
        
