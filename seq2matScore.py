#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Returns the score of a given sequence into a matrix. The higher the score,
the better the sequence adapts to the matrix
@author: rodri
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:22:29 2018

@author: rodri
"""
#%%
import bioio
import bioseq
import frecuencias
#%%
"""
Given a dinucleotide frequency matrix and a dict of sequences, returns the
score of each sequence in the matrix. All the sequences must have the same
length that the frequency matrix. If not, they are truncated and the exceeding
nucleotides are not accounted.
"""
def seq2matScore(seqs,freqs,k=2, outFile=None):
    import bioseq
    import numpy as np
    ct=bioseq.codonTable()
    ict=bioseq.inverseCodonTable()
    loss={}
    bases = ['T', 'C', 'A', 'G']
    combos = [a+b for a in bases for b in bases] #extend to k
    computedStats=False
    maxLoss=0
    meanLoss=0
    sdLoss=0
    flen=len(freqs["AA"])
    for key in seqs.keys():
        seq=seqs[key]
        if(len(seq)!=flen+k-1):
            print("WARNING: the sequence ",k," length is not equal to the size of the frequency matrix")
        temp=[]
        lossSeq=0
        for i in range(0,len(seq)-k+1):
            combo=seq[i:i+k].upper()
            freqc={}
            for c in combos:
                freqc[c]=freqs[c][i]
            if not computedStats:
                maxLoss+= min(list(freqc.values()))-max(list(freqc.values()))
                meanLoss+= np.mean(list(freqc.values()))-max(list(freqc.values()))
                sdLoss+= np.std(list(freqc.values()))-max(list(freqc.values()))
                
            lossSeq+=freqs[combo][i]-max(list(freqc.values()))
        loss[key]=lossSeq
        computedStats=True
    if(outFile!=None):
        fw=open(outFile, "w")
        fw.write("minLoss\t 0\n")
        fw.write("maxLoss\t"+str(maxLoss)+"\n")
        fw.write("meanLoss\t"+str(meanLoss)+"\n")
        fw.write("sdLoss\t"+str(sdLoss/len(combos))+"\n")
        for key in loss.keys():
            fw.write(key+"\t"+str(loss[key])+"\n")
        fw.close()       
    
    return {"loss":loss, "minLoss":0,"maxLoss":maxLoss, "meanLoss":meanLoss, "sdLoss":sdLoss/len(combos)}
freqs=frecuencias.readFreqs("/home/rodri/data/nucleosomas/remasterAssym/resultsAssym/Nuc_SP_woNuc1NucT_150bp_extended-freqs-2.tsv")
seqs=bioio.readFasta("/home/rodri/data/nucleosomas/Nuc10_SP.fasta")
score=seq2matScore(seqs,freqs, outFile="/home/rodri/data/nucleosomas/scoreTest.txt")
print(score)
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\nbiosearch version 0.1')
    print('\nUsage (python3): ')
    print('\tpython biosearch.py -fasta path -seqs path -out path')
    print('\n\t--- Input arguments ---')
    print('\t-fasta file:\tFasta file with nucleotide sequences. Sequence names must have pattern "chromosomeX:start-end".')
    print('\t-seqs file:\tTab delimited file with motifs (column 1) and sequences of the motif (column 2). For example:')
    print('\t\t\tm1\tATCG')
    print('\t\t\tm1\tATCT')
    print('\t\t\tm2\tTTTTGGGG')
    print('\t\t\tm2\tTTTTGCGG')
    print('\t\t\tm2\tATTTGGGC')
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
            

    #processing
    print("Searching for sequences...")
    t0=time.clock()
    seqFind(fastaPath, seqsPath, outPath)
    print("... finished in ", round((time.clock()-t0)/60., 2), " minutes")
        
