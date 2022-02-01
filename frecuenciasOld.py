# -*- coding: utf-8 -*-
"""
Obtener secuencias de FASTA desde BED y calcular frecuencias de nucleótidos

@author: rodri
"""
#%% 1) Read BED
#Programmed for this task. bedPath is the path to the BED file
#BED file must be headerless and have the same length (3rd-2nd column) in all the rows
#sl is the sequence length, only for checking that are all the same length
def readStarts(bedPath, sl=2000):
    import csv
    csvfile=open(bedPath)
    reader=csv.DictReader(csvfile, fieldnames=['chr','start','end'], delimiter="\t")
    
    seqs={}
    for row in reader:
        ch=row["chr"]
        if(not ch in seqs.keys()):
            seqs[ch]=[]
        seqs[ch].append(int(row["start"]))   
        dif=int(row["end"])-int(row["start"])
        if(dif !=sl):
            print("ERROR: length is ",dif, "at",ch,row["start"])
            break
    print([len(seqs[x]) for x in seqs.keys()])
    return seqs
    
#%% 2) Get Sequences from full-genome fasta files
#Modified from seqview fasta() method
def fasta(path=""):
    import re
    f=open(path)
    reader=f.readlines()
    f=open(path)
    ff={}
    k=re.sub("\(.*\)", "", re.sub(">", "", reader[0])).strip()
    i=0
    while i<len(reader):
        seq=""
        while(i<len(reader)):
            line=f.next()
            i+=1
            if(line.startswith(">")):
                ff[k]=seq.replace("\n","")
                k=re.sub("\(.*\)", "", re.sub(">", "", line)).strip()
                break
            else:
                seq+=line
    ff[k]=seq.replace("\n","")
    return ff
    
 
#%% 3) Get the sequences from genome (fasta) corresponding to positions (search) with length windowSize
def getSeqsFromPos(search, fasta, windowSize=2000, comp=False):
    seqs={}
    for track in search.keys():
        positions=search[track]
        seqs[track]=[]
        for i in range(len(positions)):
            start=int(positions[i])
            end=start+windowSize
            ss=fasta[track][start:end].upper()
            if(len(ss)==windowSize): #to avoid adding small sequences
                if(comp):
                    ss=complement(ss)
                seqs[track].append(ss)
            else:
                print("Size is", len(ss), start, end, end-start, track, windowSize)
    #convert to array (could have done it directly into a list, though)
    seqlist=[]
    for k in seqs.keys():
        for x in seqs[k]:
            seqlist.append(x)
    return seqlist



#%% 4) Compute frequencies
# (modified from sesion1methods.py)
# seqs - an array of nucleotide sequences of the same length (uppper case)
# returns a dictionary with nucleotides as keys and arrays of frequencies for 
#  the key, with a length equal to the length of each sequence in seqs
def freq(seqs):
    import numpy as np
    
    sl=len(seqs[0])
    fr={}
    for k in ["A","C","T","G", "N"]:
        fr[k]=np.zeros(sl)
    for i in range(sl): #for each nucleotide
        for s in seqs:      #for each sequence
            fr[s[i]][i]+=1
    for k in fr.keys():
        for i in range(len(fr[k])):
            fr[k][i]=(float)(fr[k][i])/len(seqs)
    return fr

#%% 4) Compute frequencies (version 2)
"""
 (modified from sesion1methods.py)
 seqs - an array of nucleotide sequences of the same length (uppercase)
 k    - the length of the sequences to compute frecuencies (usually k<4)
           k=1 (nucleotides, default), k=2 (dinucleotides), k=3 (codons)
 step - if k>1, how to determine k-mers (default 1)
       For example if k=3 and we have 6 positions and step=k, we will get
        codons 123 and 456. If for the same example step=1 we will get codons
        123 234 345 456.
 genome - If a dictionary of str is provided (default None) kmer frequencies on
     seqs are normalized by the kmer frequency on such dict (taken as a whole single seq)
 log2 - if normalization is required, whether to provide it in log2 (True, default) or not
 
 returns a dictionary with kmers as keys and arrays of frequencies for 
  the key, with a length equal to the length of each sequence in seqs or
  smaller (proportional to the step)
        """
def freq2(seqs, k=1, step=1, genome=None, log2=True):
    import numpy as np
    import bioio
    import time
    t00=time.clock()

    
    #0) Determine all the kmers
    import itertools
    keys=["".join(x) for x in list(itertools.product(["A","C","G","T"],repeat=k))]
    sl=len(seqs[0])
      
    #1) Compute frequencies
    print("Computing frequencies...")
    t0=time.clock()
    fr={}
    for key in keys:
        fr[key]=np.zeros(len(range(0,sl-k+1,step)))
    cont=0
    for i in range(0,sl-k+1,step):
        for s in seqs:
            fr[s[i:i+k]][cont]+=1
        cont+=1
    
    for key in keys:
        for i in range(len(fr[key])):
            fr[key][i]=(float)(fr[key][i])/len(seqs)
    print("\t... in ", (time.clock()-t0), "s")
    #2) Genome frequency normalization if genomic sequence is provided
    if (genome!=None):
        print("Normalizing by genome...")
        t0=time.clock()
        #2a) Compute genomic frequencies
        fgenome={}
        for key in keys:
            fgenome[key]=0
        for ch in genome.keys():
            gs=genome[ch].upper()
            for i in range(0,len(gs)-k+1,step):
                kmer=gs[i:i+k]
                if(("N" in kmer)==False):
                    fgenome[kmer]+=1
        genlen=sum(fgenome.values())
        for key in keys:
            fgenome[key]/=genlen*1.
        print("\t... in ", (time.clock()-t0), "s")
        t0=time.clock()
        #2b) Normalize by genomic frequencies
        for key in keys:
            for i in range(len(fr[key])):
                if(log2):
                    fr[key][i]=np.log2(fr[key][i]/fgenome[key])
                else:
                    fr[key][i]/=fgenome[key]
                    
                
                
        print("\t... in ", (time.clock()-t0), "s")
    
    print("Frequencies computations takes ", (time.clock()-t00), "s")       
    return fr
#freq2(fpombe["seqs"])
#ff=freq2(["ACGGT","TTTTT"], k=3,step=1, genome=fpombe["genome"])
freq3=freq2(fpombe["seqs"], k=3,step=1, genome=fpombe["genome"])

#%% 4b) Smoothing frequencies 
"""
freqs - a dictionary as returned by freq() above
window - size of the window to be smoothed (by average)
step - amount of nucleotides to advance the window
returns a dictionary with smoothed frequencies in key 'freq' and the starting
           positions of each window used in key 'xpos'
"""
def smooth(freqs, window=2, step=2):
    import numpy as np
    freqs2={}
    for k in freqs.keys():
        xpos=[]
        freqs2[k]=[]
        for i in range(0, len(freqs[k])-window,step): #recorre de 0 al final en saltos de step
            freqs2[k].append(np.mean(freqs[k][i:i+window])) #hace la media para los valores en el intervalo i-i+window
            xpos.append(int(round(i+window*.5))+1)#+1 to avoid start from 0
    return {"freq":freqs2, "xpos":xpos}
        
#%% 5) Write AT content to file in path, given the frequencies (freqs)
# (coded for this purpose)
def writeFreqs(freqs, path):
    #Writing  single nucleotides
    f=open(path+".tsv","w")
    f.write("".join(["\t".join(freqs.keys()),"\t", "AT", "\n"]))
     
    for i in range(len(freqs["A"])):
        atc=freqs["A"][i]+freqs["T"][i]
        ff=""
        for k in freqs.keys():
            ff+=str(round(freqs[k][i],3))+"\t"
        f.write("".join([ff,str(round(atc,3)),"\n"]))
            
    f.close()
    
    #Writing dinucleotide freqs
    difreqs=kmerFreq(freqs,2)
    f=open(path+"2.tsv","w")
    f.write("\t".join(difreqs.keys())+"\n")
     
    for i in range(len(difreqs["AA"])):
        ff=""
        for k in difreqs.keys():
            ff+=str(round(difreqs[k][i],4))+"\t"
        f.write("".join([ff,"\n"]))
            
    f.close()
    
    #Writing trinucleotide freqs
    trifreqs=kmerFreq(freqs,3)
    f=open(path+"3.tsv","w")
    f.write("\t".join(trifreqs.keys())+"\n")
     
    for i in range(len(trifreqs["AAA"])):
        ff=""
        for k in trifreqs.keys():
            ff+=str(round(trifreqs[k][i],4))+"\t"
        f.write("".join([ff,"\n"]))
            
    f.close()

    
#%% 5b) Computes AT content and returns it, given the frequencies (freqs)
def atContent(freqs):
    at=[]
    for i in range(len(freqs["A"])):
        at.append(freqs["A"][i]+freqs["T"][i])
    return at
    
#%%    
"""
 Computes k-mer frequencies (recommended k<4)
 freqs - dictionary of frequencies of single nucleotides as returned by freq()
 k -  size of the kmers (default dinucleotides, k=2)
 returns a dictionary with the kmer sequences as keys and an frequency array as
 values. The array length will be the length of the arrays in freqs divided by k
"""
def kmerFreq(freqs,k=2):
    #Compute all the k-mer combinations
    import itertools
    kmers={}
    keys=["".join(x) for x in list(itertools.product(["A","C","G","T"],repeat=k))]
    for key in keys:
        kmers[key]=[]
        for i in range(0,len(freqs["A"])-k+1,k):
            f=1
            for x in range(i,i+k):
                f*=freqs[key[x-i]][x] #add position frequency
            kmers[key].append(f)
        
    return kmers
#tal2=kmerFreq(tal["freq"],3) 
#sum([tal2[x][1] for x in tal2.keys()])
    
#%%  6) Draw content line 
def drawAT(freqs, xpos, title="", vline=False):
    at=atContent(freqs)
      
    #Mostrando de manera gráfica el skew
    import matplotlib.pyplot as plt
    plt.plot(xpos, at, label='skew')
    plt.ylabel('AT content (%)')
    plt.xlabel('position')
    plt.title(title)
    if(vline!=False):
        plt.plot([vline,vline],[min(at),max(at)], ls="--", color="#aaaaaa")
    plt.show()


    
#%%
def drawNucleotides(freqs, xpos, nucleotides=["A","T","C","G"], title="Nucleotide content", vline=[]):
    import matplotlib.pyplot as plt
    
    colors={"A":"g", "C":"b", "T":"r", "G":"k"}
    maxy=0
    miny=100
    for nuc in nucleotides:
        at=freqs[nuc]
        yvals=map(lambda x: x*100, at)
        plt.plot(xpos, yvals, label=nuc, color=colors[nuc])
        if(max(at)>maxy):
            maxy=max(at)
        if(min(at)<miny):
            miny=min(at)
    plt.legend(loc='center left', frameon=False, ncol=1, fontsize=9, handletextpad=0)

    for v in vline:
        plt.plot([v,v],[miny*100,maxy*100], ls="--", color="#aaaaaa")
       
    plt.ylabel('content (%)')
    plt.xlabel('position')
    plt.title(title)
    plt.show()

#drawNucleotides(fs2["freq"],fs2["xpos"], vline=100, title="S cerevisiae nucleotide content")    

#%%---------------------- INVERSA COMPLEMENTARIA ----------------------
#Retorna la cadena (inversa) complementaria a seq
def complement(seq):
    comp=[]
    complementos={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}#tabla de conversión
    for s in seq.upper()[::-1]:                 #recorrer en sentido inverso
        comp.append(complementos[s])
    return ''.join(comp)




#%% PIPELINE PATA CALCULAR FRECUENCIAS
#positionFile - bed file with positions. All intervals must have the same size and chromosome names must coincide with the ones on the fasta
#intervalSize - size of the intervals on positionFile for error-checking
#seqFile      - fasta file with sequences
#window       - smoothing window for drawing
#step         - smoothing step for drawing
#extendRight/Left - extends the interval to the Right (towards 3') or to the Left (towards 5') this number of nucleotides
#vline        - list with positions at which to draw vertical lines
#outFile      - file with nucleotide frequencies for the sequence intervals in positionFile (default name "frequencies"). 
#                   No file extension should be included in the name,'.tsv' will be automatically added
#                   Additionals files 'outfile2.tsv' and 'outfile3.tsv' will be generated for di and tri-nucleotide frequencies
#positionFile - optional additional bed file with positions on the reverse strand
#               reverse strand will be computed from seqFile and these intervals will be added to frequencies
def pipelineFreq(positionFile=None, intervalSize=None, seqFile=None, window=20, step=2, vline=[], outFile="frequencies", positionRevFile=None, extendRight=0, extendLeft=0):
    import numpy as np    
    import bioio
    dataFASTA=bioio.readFasta(seqFile) #cargar las secuencias
        
    # Cargar la cadena normal
    seqs=[]
    if(positionFile!=None):
        starts=readStarts(positionFile, intervalSize)
        for k in starts.keys():
            starts[k]=np.array(starts[k])-extendLeft
        seqs=getSeqsFromPos(starts,dataFASTA,intervalSize+extendRight)
    
    # Cargar la cadena complementaria
    seqsi=[]
    if(positionRevFile!=None):
        startsi=readStarts(positionRevFile, intervalSize)
        for k in list(startsi.keys()):
            startsi[k]=np.array(startsi[k])-extendRight
        seqsi=getSeqsFromPos(startsi,dataFASTA,intervalSize+extendRight+extendLeft, comp=True)
    
    # Calcular las frecuencias:
    seqst=seqs+seqsi
    fst=freq(seqst)

    # Dibujar las frecuencias:
    fs2t=smooth(fst,window,step)
    drawAT(fs2t["freq"], fs2t["xpos"], "", vline=vline)
    drawNucleotides(fs2t["freq"],fs2t["xpos"], vline=vline, title="")    
    
    # Guardar archivos
    writeFreqs(fs2t["freq"], outFile)
    return {"genome": dataFASTA, "seqs":seqst, "freq":fst, "smoothedFreq":fs2t}    
#%%
#Returns a dictionary where keys are the 3-nucleotide sequence and value the corresponding aminoacid
def codonTable():    
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table

#Returns a dictionary where keys are the aminoacid and values the list of all 3-nuc corresponding sequences
def inverseCodonTable():
    import numpy as np
    ct=codonTable()
    ict={}
    for aa in np.unique(ct.values()):
        ict[aa]=filter(lambda k:ct[k]==aa, ct.keys())
    return ict
    
#%% REMASTER
"""
Returns the sequence that has the largest frequency according to freqs but is 
equivalent to seq in codons
freqs: dictionary with 64 keys (all possible codons) and 
        frequency arrays as values
seq:   nucleotide sequence to remasterize. If seq length is larger than freqs
        values length, freqs is reused from the beginning once fully iterated. 
        E.g. if seq length is 100 and we have freqs for 50 nucleotides, 
        positions 0 and 50 will be remastered with freq elements at postion 0.
offset: in the case that the remastering has to start with a freq position 
        different than the first one, it can be specified with this parameter. 
        (default 0)
"""
def remaster(seq,freqs,offset=0):
    ct=codonTable()
    ict=inverseCodonTable()
    mseq=[]
    gain=0
    flen=len(freqs["AAA"])
    for i in range(0,len(seq)-2,3):
        codon=seq[i:i+3].upper()
        codons=ict[ct[codon]]
        freqc={}
        ixf=(i+offset)%flen
        for c in codons:
            freqc[c]=freqs[c][ixf]
        selcod=max(freqc, key=freqc.get)
        mseq+=selcod
        gain+=max(freqc.values())-freqs[codon][ixf]
    return "".join(mseq)
        
#################### APPLICATION OF METHODS ##################
#%% FREQUENCIES AND GRAPHICS
#Get frequencies and freq graphics    
tal=pipelineFreq(positionFile="/home/rodri/Documentos/investigacion/IBFG/remaster/test/Nuc_SP_total_153bp.bed",
             intervalSize=500, window=20, step=1, vline=[100,250,400],
             extendRight=0, extendLeft=0,
             seqFile="/home/rodri/workspace/python/IBFG/frecuencias/S_pombe_mio.fasta",
             #positionRevFile="/home/rodri/Documentos/workspace/python/IBFG/frecuencias/Triadas_Nuc1,2,3_SP_m_500bp.bed",
             outFile="/home/rodri/Documentos/workspace/python/IBFG/frecuencias/ATspt"             
             )

#Remastering tests
freq3=kmerFreq(tal["freq"],3)
#%%
s="acgtga"
rs={}
k="remastered sequence for "+s
rs[k]=remaster(s,freq3)

import os
os.chdir("/home/rodri/Documentos/workspace/python/IBFG")
import bioio
bioio.writeFASTA(rs,"remastered.fa")


#%% Checking pombe to octosporus
#Get frequencies and freq graphics    
fpombe=pipelineFreq(positionFile="/home/rodri/Documentos/investigacion/IBFG/remaster/Nuc_SP.bed",
             intervalSize=153, window=20, step=1, vline=[75],
             extendRight=0, extendLeft=0,
             seqFile="/home/rodri/workspace/python/IBFG/frecuencias/S_pombe_mio.fasta",
             positionRevFile=None,
             outFile="/home/rodri/Documentos/investigacion/IBFG/remaster/ATspt"             
             )
#%%
#Remastering tests
#freq3=kmerFreq(fpombe["freq"],3)
freq3=freq2(fpombe["seqs"], k=3,step=1, genome=fpombe["genome"])

#%%
import os
os.chdir("/home/rodri/workspace/python/IBFG/remaster")
import bioio
#seqs=bioio.readFasta("/home/rodri/Documentos/investigacion/IBFG/remaster/secuencias_ura4.fasta")             
seqs=readFasta("/home/rodri/Documentos/investigacion/IBFG/remaster/secuencias_ura4.fasta")             
#%%
s=seqs["ura4_japonicus"]
rs={}
k="remastered sequence for ura4_japonicus"
rsl=seqs["ura4japonicus_remasterizado"].upper() #la obtenida por Luis
maxmatch=0
maxcon=""
maxoff=0
rem0=""
freq3=freq2(fpombe["seqs"], k=3,step=1, genome=fpombe["genome"])

#for off in range(50):
for off in range(len(freq3["AAA"])):
#for off in [58]:
    print(off)
    rs[k]=remaster(s,freq3,off)
    rs[k]   #La obtenida por nosotros
    con=""
    match=0
    for i in range(len(rsl)):
        if(rsl[i]==rs[k][i]):
            con+=rsl[i]
            match+=1
        else:
            con+="-"
    if(maxmatch<match):
        maxcon=con
        maxoff=off
        maxmatch=match
        rem0=rs[k]
    print("Match:",match, 1.*match/len(rsl))
#%%
print("Coincidence in nucleotides: ",1.*maxmatch/len(rsl))
print("Coincidence in codons: ",
      sum([1 if rsl[i:i+3]==rem0[i:i+3] else 0 for i in range(0,len(rsl),3)])/(len(rsl)/3.))
print("Coincidence in amino acids: ",
      sum([1 if (rem0[i:i+3] in ict[ct[rsl[i:i+3]]])==True else 0 for i in range(0,len(rsl),3)])/(len(rsl)/3.))
print(maxcon)        
print(maxoff)
for i in range(0,len(rsl)/10,3):
    if(rem0[i:i+3]!=rsl[i:i+3]):
        print(rem0[i:i+3]," --> ",freq3[rem0[i:i+3]][i+maxoff])
        print(rsl[i:i+3]," --> ",freq3[rsl[i:i+3]][i+maxoff])
        print()
        
    
#    print ("con for offset", off, "is", con)
#%%
rem0=remaster(s,freq3,maxoff)
#%%           
import matplotlib.pyplot as plt
ct=codonTable()
ict=inverseCodonTable()

maxy=0
miny=100
for k in ict.keys():
    for nuc in ict[k]:
        at=freq3[nuc]
        #yvals=map(lambda x: x*100, at)
        yvals=map(lambda x: x, at)
        plt.plot(range(len(freq3["AAA"])), yvals, label=nuc)#, color=colors[nuc])
        if(max(at)>maxy):
            maxy=max(at)
        if(min(at)<miny):
            miny=min(at)
    plt.legend(loc='center left', frameon=False, ncol=1, fontsize=9, handletextpad=0)
       
    #plt.ylabel('freq (%)')
    plt.ylabel('log2fc')
    plt.xlabel('position')
    #plt.plot([25,25],[min(at),max(at)], ls="--", color="#aaaaaa")
    plt.title("Codon freqs for aa="+k)
    plt.savefig("aafreqs/"+k+'.tif') 
    plt.show()

#drawNucleotides(fs2["freq"],fs2["xpos"], vline=100, title="S cerevisiae nucleotide content")    





#%%  REMASTER TESTS
#Read sequence to remasterize
import os
os.chdir("/home/rodri/Documentos/workspace/python/IBFG")
import bioio
seqs=bioio.readFasta("/home/rodri/Documentos/investigacion/IBFG/remaster/secuencias_ura4.fasta")             
s=seqs["ura4_japonicus"]

rsl=seqs["ura4japonicus_remasterizado"].upper() #la obtenida por Luis

# Read  nucleosome sequences and pombe genome
fpombeSR=pipelineFreq(positionFile="/home/rodri/Documentos/investigacion/IBFG/remaster/NucleosomasPos2.bed",
             intervalSize=150, window=20, step=1, vline=[75],
             extendRight=0, extendLeft=0,
             seqFile="/home/rodri/Documentos/workspace/python/IBFG/frecuencias/S_pombe_mio.fasta",
             positionRevFile=None,
             outFile=None             
             )
             
#%% Compute frequencies
req3=SRfreq2(fpombeSR["seqs"], k=3,step=1, genome=fpombeSR["genome"])
#%% Remasterize
rem0=remaster(s,freq3,59) #according to the information provided, there's an ATG in nucleosome 1 at pos 59

#%%

