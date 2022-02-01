# -*- coding: utf-8 -*-
"""
Obtener secuencias de FASTA desde BED y calcular frecuencias de nucleótidos

@author: rodri
"""
#%% 1) Read BED
#Programmed for this task. bedPath is the path to the BED file
#BED file must be headerless and have the same length (3rd-2nd column) in all the rows
#sl is the sequence length, only for checking that are all the same length
def readStarts(bedPath, sl=None):
    import csv
    csvfile=open(bedPath)
    reader=csv.DictReader(csvfile, fieldnames=['chr','start','end'], delimiter="\t")
    
    seqs={}
    for row in reader:
        ch=row["chr"]
        if(not ch in seqs.keys()):
            seqs[ch]=[]
        seqs[ch].append(int(row["start"]))   
        dif=int(row["end"])-int(row["start"])+1
#        if(sl!=None and dif !=sl):
#            print("ERROR: length is ",dif, "at",ch,row["start"])
#            break
    #print [len(seqs[x]) for x in seqs.keys()]
    return (seqs,dif)
    
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
            #line=f.next()
            line=next(f) #python 3 syntax
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
             #else:
             #   print("Size is", len(ss), start, end, end-start, track, windowSize)
    #convert to array (could have done it directly into a list, though)
    seqlist=[]
    for k in seqs.keys():
        for x in seqs[k]:
            seqlist.append(x)
    return seqlist

#%% 3) Get the sequences from genome (fasta) corresponding to positions given on a bed file with length windowSize
"""
bed - path to bed file that must contain 5 columns: 
       1st with chromosome names equals to the ones in the fasta dictionary
       2nd and 3rd with start and end that matches with window length
       5th indicating forward (+) or reverse (-) strand
fasta - dictionary of genome chromosomes as returned by bioio.readFasta
windowSize - a number that must match with all the interval sizes (end-start)

Returns a list of sequences at such positions in the genome 
(reverse complementary sequences if the position is in the reverse strand) 
"""
def getSeqsFromBed(bed, fasta, windowSize=2000):
    import csv
    csvfile=open(bed)
    reader=csv.DictReader(csvfile, fieldnames=['chr','start','end', "xx", "strand"], delimiter="\t")
    
    seqs=[]
    for row in reader:
        start=int(row["start"])
        end=int(row["end"])
        ss=fasta[row["chr"]][start:end].upper()
        #if(len(ss)==windowSize): #to check bed sizes
        if(row["strand"]=='-'):
            ss=complement(ss)
        seqs.append(ss)
        #else:
        #    print("Error: size is", len(ss), start, end, end-start, row["chr"], windowSize)
    return seqs
#tal=getSeqsFromBed("/home/rodri/data/nucleosomas/plusOne/Nuc+1_ok.bed",bioio.readFasta("/home/rodri/data/annotations/sp/Spombe.fasta"), 150)

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
def freq2(seqs, k=1, step=1, genome=None, log2=False, nonorm=False):
    import numpy as np
    import time
    t00=time.clock()
    print("Computing frequencies for k-mers of size ",k,"...")

    
    #0) Determine all the kmers
    import itertools
    keys=["".join(x) for x in list(itertools.product(["A","C","G","T"],repeat=k))]
    sl=len(seqs[0])
      
    #1) Compute frequencies
    #t0=time.clock()
    fr={}
    contFr={}
    for key in keys:
        fr[key]=np.zeros(len(range(0,sl-k+1,step)))
        contFr[key]=0
    cont=0
    for i in range(0,sl-k+1,step):
        for s in seqs:
           kmer=s[i:i+k]
           if(("N" in kmer)==False):
               fr[kmer][cont]+=1
               contFr[kmer]+=1
        cont+=1
    
    for key in keys:
        #print(key+":", contFr[key])
        contFr[key]/=len(seqs)*(len(seqs[0])-1)
        for i in range(len(fr[key])):
            fr[key][i]=(float)(fr[key][i])/len(seqs)
            #print(key+":"+fr[key])
    #print(contFr)
#    print("\t... in ", (time.clock()-t0), "s")
    
    #2) Genome frequency normalization if genomic sequence is provided
    if (genome!=None and nonorm==False):
        print("\tNormalizing by genome...")
        #t0=time.clock()
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
            #print(key+":", fgenome[key])
            fgenome[key]/=genlen*1.
            #print(key+":", fgenome[key])
 #       print("\t... in ", (time.clock()-t0), "s")
        #t0=time.clock()
        #2b) Normalize by genomic frequencies
        for key in keys:
            for i in range(len(fr[key])):
                if(log2):
                    if(fr[key][i]!=0):
                        #fr[key][i]=np.log2(fr[key][i]/fgenome[key]) #lo que consideraba yo: log2 y den. genómico
                        #fr[key][i]=np.log(fr[key][i]/fgenome[key]) #con ln
                        #fr[key][i]=np.log(fr[key][i]/contFr[key]) #ln y den. de secuencias
                        fr[key][i]=np.log2(fr[key][i]/contFr[key]) #log2 y den. de secuencias
                else:
                    fr[key][i]/=fgenome[key]
                    
                
                
#        print("\t... in ", (time.clock()-t0), "s")
    
    print("Frequencies computations takes ", (time.clock()-t00), "s")       
    return fr

#freq2(fpombe["seqs"])
#ff=freq2(["ACGGT","TTTTT"], k=3,step=1, genome=fpombe["genome"])
#freq3=freq2(fpombe["seqs"], k=3,step=1, genome=fpombe["genome"])

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
#        for i in range(0, len(freqs[k])-window,step): #recorre de 0 al final en saltos de step
#            freqs2[k].append(np.mean(freqs[k][i:i+window])) #hace la media para los valores en el intervalo i-i+window
#            xpos.append(int(round(i+window*.5))+1)#+1 to avoid start from 0
        for i in range(0, len(freqs[k]),step): #recorre de 0 al final en saltos de step
            freqs2[k].append(np.mean(freqs[k][i:min(i+window, len(freqs[k]))])) #hace la media para los valores en el intervalo i-i+window
            xpos.append(int(round(i+window*.5))+1)#+1 to avoid start from 0
    return {"freq":freqs2, "xpos":xpos}
        
#%% 5) Write AT content to file in path, given the frequencies (freqs)
# (coded for this purpose)
def writeFreqs(freqs0, path):
    #Writing  single nucleotides
    for kk in freqs0.keys():
        if("freq" in freqs0[kk].keys()):
            freqs=freqs0[kk]["freq"]
        else:
            freqs=freqs0[kk]
        f=open(path+"-"+str(kk)+".tsv","w")
        f.write("".join(["\t".join(freqs.keys()), "\n"]))
         
        for i in range(len(freqs[list(freqs.keys())[0]])):
            ff=""
            for k in freqs.keys():
                ff+=str(round(freqs[k][i],3))+"\t"
            f.write("".join([ff,"\n"]))
        f.close()
        
def writeFreqsDict(freqs, path):
        f=open(path,"w")
        f.write("".join(["\t".join(freqs.keys()), "\n"]))
         
        for i in range(len(freqs[list(freqs.keys())[0]])):
            ff=""
            for k in freqs.keys():
                ff+=str(round(freqs[k][i],3))+"\t"
            f.write("".join([ff,"\n"]))
        f.close()
    
    
def writeFreqs0(freqs, path):
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

    
#%% 5b) Computes AT content and returns it, given the frequencies of single nucleotides (freqs)
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
def drawAT(freqs, xpos, title="", vline=False, imgFolder=None):
    at=atContent(freqs)
      
    #Mostrando de manera gráfica el skew
    import matplotlib.pyplot as plt
    plt.ylabel('AT content') #set (%) if not normalized
    plt.xlabel('position')
    plt.title(title)
    #print(len(at),len(xpos))
    plt.plot(xpos, at)
    if(vline!=False):
        plt.plot([vline,vline],[min(at),max(at)], ls="--", color="#aaaaaa")
    if imgFolder!=None:
        import os
        if os.path.isdir(imgFolder)==False:
            os.makedirs(imgFolder)
        plt.savefig(imgFolder+'/ATcontent.tif')
        plt.close()
    else:
          plt.show()
#%%
def drawNucleotides(freqs, xpos, label="content", nucleotides=["A","T","C","G"], 
                    title="Nucleotide content", vline=[], imgFolder=None, imgName=None, backend="agg"):
    import matplotlib.pyplot as plt
    import os
    
    if(backend!=None):
        plt.switch_backend(backend)
    
    colors={"A":"g", "C":"b", "T":"r", "G":"k"}
    maxy=max(freqs[list(freqs.keys())[0]])
    miny=min(freqs[list(freqs.keys())[0]])
    for nuc in nucleotides:
        at=freqs[nuc]
        #yvals=map(lambda x: x*100, at)
        yvals=at
        plt.plot(xpos, yvals, label=nuc, color=colors[nuc])
        if(max(at)>maxy):
            maxy=max(at)
        if(min(at)<miny):
            miny=min(at)
    plt.legend(loc='center left', frameon=False, ncol=1, fontsize=9, handletextpad=0)

    for v in vline:
       #plt.plot([v,v],[miny*100,maxy*100], ls="--", color="#aaaaaa")
       plt.plot([v,v],[miny,maxy], ls="--", color="#aaaaaa")
       
    plt.ylabel(label)
    plt.xlabel('position')
    plt.title(title)
    #plt.show()
    
    if imgFolder!=None:
        if os.path.isdir(imgFolder)==False:
            os.makedirs(imgFolder)
        if imgName!=None:
            plt.savefig(imgFolder+'/'+imgName+'.tif')             
        else:
            plt.savefig(imgFolder+'/ntFreqs.tif') 
        plt.close()

#%%
def drawFreqs(freqs, xpos, label="content", vline=[], imgFolder=None):
    import matplotlib.pyplot as plt
    import os
    
    plt.switch_backend("agg")

    #colors={"A":"g", "C":"b", "T":"r", "G":"k"}
#    maxy=max(freqs[freqs.keys()[0]])
#    miny=min(freqs[freqs.keys()[0]])
    for k in freqs.keys():
        at=freqs[k]
        #yvals=map(lambda x: x*100, at)
        yvals=at
        plt.plot(xpos, yvals, label=k, color="blue")
#        if(max(at)>maxy):
#            maxy=max(at)
#        if(min(at)<miny):
#            miny=min(at)
        #plt.legend(loc='center left', frameon=False, ncol=1, fontsize=9, handletextpad=0)
    
        for v in vline:
           #plt.plot([v,v],[miny,maxy], ls="--", color="#aaaaaa")
           plt.plot([v,v], ls="--", color="#aaaaaa")
           
        plt.ylabel(label) #set (%) if not normalized
        plt.xlabel('position')
        #plt.title(title)
        #plt.show()
        
        if imgFolder!=None:
            if os.path.isdir(imgFolder)==False:
                os.makedirs(imgFolder)
            plt.savefig(imgFolder+"/"+k+'.tif') 
            plt.close()

#drawNucleotides(fs2["freq"],fs2["xpos"], vline=100, title="S cerevisiae nucleotide content")    
#%%
def drawCodons(freq, xpos, label="content", imgFolder=None):       
    import matplotlib.pyplot as plt
    import bioseq
    
    plt.switch_backend("agg")

    ict=bioseq.inverseCodonTable()
    
    maxy=max(freq[list(freq.keys())[0]])
    miny=min(freq[list(freq.keys())[0]])
    for k in list(ict.keys()):
        for nuc in ict[k]:
            at=freq[nuc]
            #yvals=map(lambda x: x*100, at)
            yvals=list(map(lambda x: x, at))
            plt.plot(xpos, yvals, label=nuc)#, color=colors[nuc])
            if(max(at)>maxy):
                maxy=max(at)
            if(min(at)<miny):
                miny=min(at)
        plt.legend(loc='center left', frameon=False, ncol=1, fontsize=9, handletextpad=0)
           
        #plt.ylabel('freq (%)')
        plt.ylabel(label)
        plt.xlabel('position')
        #plt.plot([25,25],[min(at),max(at)], ls="--", color="#aaaaaa")
        plt.title("Codon freqs for aa="+str(k))
        if(imgFolder!=None):
            plt.savefig(imgFolder+"/codon"+str(k)+'.tif') 
        else:
            plt.show()
        plt.close()

##%%------------------ DRAW HEATMAP ----------------
#def drawHeatmap(freqs):
#    return
#import csv
#f=open("/Users/rodri/Documents/workspace/python/IBFG/ura4sj-sp/Nuc_SP_total_153bp-freqs-2.tsv")
#reader=csv.DictReader(f, delimiter="\t")
#freqs=[]
#order=["GA", "TC", "CA", "TG", "CT", "AG", "GT", "AC", "GC", "CG", "GG", "CC", "TA", "AT", "TT", "AA"]
#for row in reader:
#     row2=[row[k] for k in order]
#     freqs.append(row2)
#    
#
#freqs=map(lambda *row: list(row), *freqs)
#import plotly    
#import plotly.graph_objs as go
#plotly.offline.init_notebook_mode(connected=True)
#
#trace = go.Heatmap(z=freqs, y=order)
#data=[trace]
#plotly.offline.plot(data, filename='/Users/rodri/Documents/workspace/python/IBFG/basic-heatmap')

#%%
#trace = go.Heatmap(z=[[1,2,3,4],[5,6,7,8]])
#data=[trace]
#plotly.offline.iplot(data, filename="/file")

#%%---------------------- INVERSA COMPLEMENTARIA ----------------------
#Retorna la cadena (inversa) complementaria a seq
def complement(seq):
    comp=[]
    complementos={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}#tabla de conversión
    for s in seq.upper()[::-1]:                 #recorrer en sentido inverso
        comp.append(complementos[s])
    return ''.join(comp)


#%%
def readFreqs(ifreqPath=""):
    f=open(ifreqPath)
    fr={}
    kmers=f.readline().replace("\t\n", "").replace("\n","").split("\t")
    for k in kmers:
        fr[k]=[]
    for line in f.readlines():
        freqs=line.replace("\n","").split("\t")
        for i in range(len(kmers)):
            fr[kmers[i]].append(float(freqs[i]))
    return fr   

#%%
"""
Inverts the frequency values, so the position 1 in freq will have the 
frequency of the position n, 2 that of n-1, etc.
freq - dictionary of <str, [float]> with frequencies are returned by e.g. readFreqs
returns - another dictionary as frec, but with values in the array inverted
"""
def invertFreqs(freq):
    ifreq={}
    for k in freq.keys():
        ifreq[k]=[]
        for i in range(len(freq[k])-1,-1,-1):
            ifreq[k].append(freq[k][i])
    return ifreq

#%% PIPELINE PATA CALCULAR FRECUENCIAS
"""
positionFile - bed file with positions. All intervals must have the same size and chromosome names must coincide with the ones on the fasta
sqsFile - alternatively to positionFile, a fasta file path can be used with the sequences
intervalSize - size of the intervals on positionFile for error-checking
genomeFile      - fasta file with sequences
window       - smoothing window for drawing
step         - smoothing step for drawing
extendRight/Left - extends the interval to the Right (towards 3') or to the Left (towards 5') this number of nucleotides
vline        - list with positions at which to draw vertical lines
outFile      - file with nucleotide frequencies for the sequence intervals in positionFile (default name "frequencies"). 
                   No file extension should be included in the name,'.tsv' will be automatically added
                   Additionals files 'outfile2.tsv' and 'outfile3.tsv' will be generated for di and tri-nucleotide frequencies
positionRevFile - optional additional bed file with positions on the reverse strand
               reverse strand will be computed from genomeFile and these intervals will be added to frequencies
seqOutFile   - optional fasta file to store sequences found at positions in positionFile in genomeFile 
"""
def pipelineFreq(positionFile=None, seqFile=None, intervalSize=None, genomeFile=None, window=20,
                 step=2, vline=[], outFile=None, imgFolder=None, kmers=[1,3],
                 positionRevFile=None, extendRight=0, extendLeft=0, 
                 smoothed=True, nonorm=False, log2=True, seqOutFile=None):
    import numpy as np
    import bioio
    
    dataFASTA=None
    if genomeFile!=None:
        print("Loading genome...")
        dataFASTA=bioio.readFasta(genomeFile) #cargar las secuencias
        
    seqst=getSeqsFromBed(positionFile, dataFASTA, intervalSize)
        
    # Cargar la cadena normal
#    seqs=[]
#    if(positionFile!=None):
#        
#        print("Loading intervals...")
#        ret=readStarts(positionFile, intervalSize)
#        starts=ret[0]
#        if(intervalSize==None):
#            intervalSize=ret[1]
#        for k in starts.keys():
#            starts[k]=np.array(starts[k])-extendLeft
#        seqs=getSeqsFromPos(starts,dataFASTA,intervalSize+extendRight)
#    elif(seqFile!=None):
#        seqs=[s.upper() for s in fasta(seqFile).values()]
#    else:
#        print("Error: either seqFile or positionFile should be provided")
#        return
#    
#    # Cargar la cadena complementaria
#    seqsi=[]
#    if(positionRevFile!=None):
#        print("Loading reverse intervals...")
#        ret=readStarts(positionRevFile, intervalSize)
#        startsi=ret[0]
#        for k in startsi.keys():
#            startsi[k]=np.array(startsi[k])-extendRight
#        seqsi=getSeqsFromPos(startsi,dataFASTA,intervalSize+extendRight+extendLeft, comp=True)
#    
#    
    # Calcular las frecuencias:
#    seqst=seqs+seqsi
    f={}
    fs={}
    for kk in kmers:
        f[kk]=freq2(seqst, k=kk,genome=dataFASTA, nonorm=nonorm, log2=log2)
        fs[kk]=smooth(f[kk],window,step)
    
    # Dibujar las frecuencias:
    if(imgFolder!=None):
        label="frequency (%)"
        if(nonorm==False):
            if(log2):
                label="log2 fold change"
            else:
                label="fold change"
            
        print("Saving figures...")
        if(smoothed):
            if(1 in fs.keys()):    
                print("\t... for single nucleotide frequencies")
                nf=freq2(seqst, k=1,genome=dataFASTA, nonorm=True)
                drawAT(nf, fs[1]["xpos"], "AT content", vline=vline, imgFolder=imgFolder)
                drawNucleotides(fs[1]["freq"],fs[1]["xpos"], label=label, vline=vline, title="", imgFolder=imgFolder) 
            if(3 in fs.keys()):
                print("\t... for codon frequencies")
                drawCodons(fs[3]["freq"], fs[3]["xpos"], label=label, imgFolder=imgFolder)
            if(2 in fs.keys()):
                drawFreqs(fs[2]["freq"], fs[2]["xpos"], label=label, imgFolder=imgFolder)
        else:
            if(1 in fs.keys()):
                xpos=range(len(f[1][list(f[1].keys())[0]]))
                drawAT(f[1], xpos, vline=vline, imgFolder=imgFolder)
                drawNucleotides(f[1],xpos, label=label, vline=vline, title="", imgFolder=imgFolder) 
            if(3 in fs.keys()):
                xpos=range(len(f[3][list(f[3].keys())[0]]))
                drawCodons(f[3], xpos, label=label, imgFolder=imgFolder)
            if(2 in fs.keys()):
                drawFreqs(f[2]["freq"], f[2]["xpos"], label=label, imgFolder=imgFolder)
        print("Figures saved")
            
    
    # Guardar archivos
    if(outFile!=None):
        print("Saving frequency files...")
        #writeFreqs(fs2t["freq"], outFile)
        if(smoothed):
            writeFreqs(fs, outFile)
        else:
            writeFreqs(f, outFile) 
        print("Files saved")
    if(seqOutFile!=None):
        print("Saving sequences at positions in bed...")
        bioio.writeFasta(seqst, seqOutFile)
    return {"genome": dataFASTA, "seqs":seqst, "freq":f, "smoothedFreq":fs}    

#res=pipelineFreq(positionFile="/home/rodri/Documentos/investigacion/IBFG/remaster/test/Nuc_SP_total_153bp.bed", 
#             genomeFile="/home/rodri/data/annotations/sp/Spombe.fasta",
#             intervalSize=153, nonorm=True)
   
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
def remaster3(seqs,freqs,offset=0):
    import bioseq
    ct=bioseq.codonTable()
    ict=bioseq.inverseCodonTable()
    mseq={}
    flen=len(freqs["AAA"])
    for k in seqs.keys():
        gain=0
        seq=seqs[k]
        temp=[]
        for i in range(0,len(seq)-2,3):
            codon=seq[i:i+3].upper()
            codons=ict[ct[codon]]
            freqc={}
            ixf=(i+offset)%(flen+2) #---> la buena!
            for c in codons:
                freqc[c]=freqs[c][ixf]
            selcod=max(freqc, key=freqc.get)
            temp+=selcod
            gain+=max(list(freqc.values()))-freqs[codon][ixf]
        mseq[k+"_remastered (gain="+str(gain)+")"]="".join(temp)
        
    return mseq

#res3=remaster3(seqs, tal["freq"][3], 0)
"""Remastering from dinucleotides instead of trinucleotides"""
def remaster2nolink(seqs,freqs,offset=0):
    import bioseq
    ct=bioseq.codonTable()  
    ict=bioseq.inverseCodonTable()
    mseq={}
    flen=len(freqs["AA"])
    for k in seqs.keys():
        seq=seqs[k]
        temp=[]
        gain=0
        for i in range(0,len(seq)-2,3):
            codon=seq[i:i+3].upper()
            codons=ict[ct[codon]]
            freqc={}
            ixf=(i+offset)%(flen+1)
            #print(i, ixf, len(freqs["AA"]), offset)
            for c in codons:
                d1c=c[0:2]
                d2c=c[1:3]
                freqc[c]=freqs[d1c][ixf]+freqs[d2c][ixf+1]
            selcod=max(freqc, key=freqc.get)
            temp+=selcod
            freqorig=freqs[codon[0:2]][ixf]+freqs[codon[1:3]][ixf+1]
            gain+=freqc[selcod]-freqorig
        mseq[str(k)+"_remastered (gain="+str(gain)+")"]="".join(temp)
    return mseq
#res2=remaster2(seqs, tal["freq"][2], 0)

#%%
#improve by taking into account the linkage dinucleotide
#With all combinations it explodes computationally
def remaster2allCombs(seq,freqs, offset=0, nn=""):
    import bioseq
    import itertools
    ct=bioseq.codonTable()  
    ict=bioseq.inverseCodonTable()
    
    codonCombo=[]
    for i in range(0,len(seq)-2,3):
        codon=seq[i:i+3].upper()
        codons=ict[ct[codon]]
        codonCombo.append(codons)  #all combinations is too much for 150 nt
    combs = ["".join(x) for x in list(itertools.product(*codonCombo))]
    maxGain=-1000
    maxSeq=""
    for cseq in combs:
        #print(seq,"->",cseq, "(", offset, ")")
        gain=0
        for i in range(0,len(cseq)-1):
            ioff=(i+offset)%len(freqs[list(freqs.keys())[0]])
            gain+=freqs[cseq[i:i+2]][ioff]-freqs[seq[i:i+2]][ioff]
            #print("\t", seq[i:i+2],"->",cseq[i:i+2],"=",freqs[seq[i:i+2]][ioff],freqs[cseq[i:i+2]][ioff], freqs[cseq[i:i+2]][ioff]-freqs[seq[i:i+2]][ioff], "(", ioff,")")
        if(nn!=""):
            ldn=seq[len(seq)-1]+nn #last dinucleotide o link dinucleotide (the one hooking to the next codon)
            cldn=cseq[len(cseq)-1]+nn
            ioff=(i+1+offset)%len(freqs[list(freqs.keys())[0]])
            gain+=freqs[cldn][ioff]-freqs[ldn][ioff]
            #print("\t LAST NUC: ", ldn,  cldn, freqs[ldn][ioff], freqs[cldn][ioff], freqs[cldn][ioff]-freqs[ldn][ioff], "(", ioff, ")")
        #if(gain>0.2):
        #print(seq,"->",cseq, " gain: ", gain)
        if(gain>maxGain):
            maxSeq=cseq
            maxGain=gain
    return (maxSeq,maxGain)


#%% remaster2allCombs generic for any k-mer frequencies (k=1-3)
def remasterAllCombs(seq,freqs, offset=0, nn="", use=None):
    import bioseq
    import itertools
    import numpy as np
    ct=bioseq.codonTable()  
    ict=bioseq.inverseCodonTable()
    klen=len(list(freqs.keys())[0])
    
    #print("-----------------",seq,"--------------------------")
    
    codonCombo=[]
    for i in range(0,len(seq)-2,3):
        codon=seq[i:i+3].upper()
        codons=ict[ct[codon]]
        codonCombo.append(codons)  #all combinations is too much for 150 nt
    combs = ["".join(x) for x in list(itertools.product(*codonCombo))]
    maxGain=-1000
    maxSeq=""
    debug=False
        
    if(use==None):
        for cseq in combs:
            #print(seq,"->",cseq, "(", offset, ")")
            gain=0
            for i in range(0,len(cseq)-klen+1):
                ioff=(i+offset)%len(freqs[list(freqs.keys())[0]])
                gain+=freqs[cseq[i:i+klen]][ioff]-freqs[seq[i:i+klen]][ioff]
                #print("\t", seq[i:i+2],"->",cseq[i:i+2],"=",freqs[seq[i:i+2]][ioff],freqs[cseq[i:i+2]][ioff], freqs[cseq[i:i+2]][ioff]-freqs[seq[i:i+2]][ioff], "(", ioff,")")
            
            for j in range(len(nn)): #linkage dinucleotides or trinucleotidas
                ldn=seq[len(seq)-1]+nn[:len(nn)+j] #last dinucleotide o link dinucleotide (the one hooking to the next codon)
                cldn=cseq[len(cseq)-1]+nn[:len(nn)+j]
                ioff=(i+j+1+offset)%len(freqs[list(freqs.keys())[0]])
                gain+=freqs[cldn][ioff]-freqs[ldn][ioff]
                #print("\t LAST NUC: ", ldn,  cldn, freqs[ldn][ioff], freqs[cldn][ioff], freqs[cldn][ioff]-freqs[ldn][ioff], "(", ioff, ")")
            #print(seq,"->",cseq, " gain: ", gain)
            if(gain>maxGain):
                maxSeq=cseq
                maxGain=gain
    else: #if codon usage is taken into account
        for cseq in combs:
            debug=False
        
            gain=0
            den=1.
                    
            for i in range(0,len(cseq)-klen+1):
                ioff=(i+offset)%len(freqs[list(freqs.keys())[0]])
                if((i+1)%3==0): #LINK DINUCLEOTIDES
                    #pgain=freqs[cseq[i:i+klen]][ioff]-freqs[seq[i:i+klen]][ioff]
                    pgain=0
                    
                    cods1=i-2
                    codon1=seq[cods1:cods1+3]
                    ccodon1=cseq[cods1:cods1+3]
                    #cods2=i+3 #why +3?
                    cods2=i+1
                    codon2=seq[cods2:cods2+3]
                    ccodon2=cseq[cods2:cods2+3]
                    
                    fsigc=np.sign(freqs[cseq[i:i+klen]][ioff])
                    fsig=np.sign(freqs[seq[i:i+klen]][ioff])
                    if(len(codon1)==len(codon2)==3):
                        if(debug):
                            print(codon1+"-"+codon2)
                        pgain=(freqs[cseq[i:i+klen]][ioff]*((use[ccodon1]/den+use[ccodon2]/den)/2))-(freqs[seq[i:i+klen]][ioff]*((use[codon1]/den+use[codon2]/den)/2))
                        
                        pgain=(freqs[cseq[i:i+klen]][ioff]*(pow(use[ccodon1],fsigc)+pow(use[ccodon2],fsigc))/2)-(freqs[seq[i:i+klen]][ioff]*(pow(use[codon1],fsigc)+pow(use[codon2],fsigc))/2)
                    #pgain=0
                    gain+=pgain
                else:       #CODON DINUCLEOTIDES
                    cods=i-i%3
                    codon=seq[cods:cods+3]
                    ccodon=cseq[cods:cods+3]
                    #if debug:
                    #    print(i,cods,"\t",ccodon, cseq[i:i+klen], "\t", codon, seq[i:i+klen])
                    fsigc=np.sign(freqs[cseq[i:i+klen]][ioff])
                    fsig=np.sign(freqs[seq[i:i+klen]][ioff])
                    pgain=(freqs[cseq[i:i+klen]][ioff]*pow(use[ccodon],fsigc))-(freqs[seq[i:i+klen]][ioff]*pow(use[codon], fsig))
                    #pgain=(freqs[cseq[i:i+klen]][ioff]*use[ccodon]/den)-(freqs[seq[i:i+klen]][ioff]*use[codon]/den)
                    #pgain=(freqs[cseq[i:i+klen]][ioff]-freqs[seq[i:i+klen]][ioff])+(use[ccodon]/den-use[codon]/den)
                    #pgain=(freqs[cseq[i:i+klen]][ioff]-freqs[seq[i:i+klen]][ioff])*(use[ccodon]/den-use[codon]/den)
                    gain+=pgain
                if debug:
                    #print("\t", seq[i:i+2],"->",cseq[i:i+2],"=(",freqs[seq[i:i+2]][ioff],"+",use[codon]/1.,") - (",freqs[cseq[i:i+2]][ioff],"+",use[ccodon]/1.,")", "(", ioff,")")
                    print("\t", seq[i:i+klen],"->",cseq[i:i+klen],": ",freqs[seq[i:i+klen]][ioff],"->",freqs[cseq[i:i+klen]][ioff],"\n\t", codon,"->",ccodon,":",use[codon]/den,"->",use[ccodon]/den,"\n\t\t Dn gain: ",pgain, "[", ioff,"]")
            
            for j in range(len(nn)): #linkage dinucleotides or trinucleotides
                ldn=seq[len(seq)-1]+nn[:len(nn)+j] #last dinucleotide o link dinucleotide (the one hooking to the next codon)
                cldn=cseq[len(cseq)-1]+nn[:len(nn)+j]
                ioff=(i+j+1+offset)%len(freqs[list(freqs.keys())[0]])
                #gain+=freqs[cldn][ioff]-freqs[ldn][ioff]
                gain+=0
                if debug:
                    print("\t LAST NUC: ", ldn,  cldn, freqs[ldn][ioff], freqs[cldn][ioff], freqs[cldn][ioff]-freqs[ldn][ioff], "(", ioff, ")")
            if debug:
                print(seq,"->",cseq, " gain: ", gain, "(max=", maxGain,")")
            if(gain>maxGain):
                maxSeq=cseq
                maxGain=gain
    if debug: 
        print("--- Remaster to", seq,"->",maxSeq, " gain: ", maxGain )
    return (maxSeq,maxGain)

#tal=remaster2sal(nucs, luis, 0, use=cu)

#%% We perform the remastering with all combinations but just in patches to avoid combinatory explosion
def remaster2patches(seqs,freqs,patchLength=18):
    mseq={}
    for k in seqs.keys():
        seq=seqs[k]
        maxSeq=""
        maxGain=0
        for i in range(0,len(seq),patchLength): #30 nt patches
            miniseq=seq[i:min(i+patchLength,len(seq))]
            (rminiseq,rminigain)=remaster2allCombs(miniseq,freqs)
            maxSeq+=rminiseq
            maxGain+=rminigain
        mseq[str(k)+"_remastered (gain="+str(maxGain)+")"]=maxSeq
    return mseq

#s="ATTGATTGTCGTCGTAACGTTGTAACGTGCTTTAAATGTC"
##s="ATTGAT"
#print(len(s))
#print(remaster2patches({1:s}, fr))
#print(remaster2patches({1:s}, fr,30))
#print(remaster2nolink({1:s}, fr))
#%% 
"""We perform the remastering having into account link dinucleotides only in
 the case that Ser/Arg/Leu are the following aa"""
def remaster2sal(seqs,freqs,offset=0, aminos=["S", "R", "L"], use=None):
#def remaster2sal(seqs,freqs,offset=0, aminos=[], use=None):
#def remaster2sal(seqs,freqs,offset=0, aminos=["S", "R", "L",  "A", "G", "P", "T", "V"]):
    import bioseq
    mseq={}
    ct=bioseq.codonTable()
    for k in seqs.keys():
        #print("-----------------------",k,"---------------------")
        seq=seqs[k]
        maxSeq=""
        maxGain=0
        i=0
        while(i<len(seq)-2):
            codon=seq[i:i+3].upper()
            miniseq=codon
            start=i
            if(i<len(seq)-5):
                nextCodon=seq[i+3:i+6]
                while len(nextCodon)!=0 and (ct[nextCodon] in aminos)==True and (len(miniseq)<33) and i<(len(seq)-5):
                    i+=3
                    miniseq+=nextCodon
                    nextCodon=seq[i+3:i+6]
            nextNuc=""
            if(i+4<len(seq)):
                    nextNuc=seq[i+3] #new
            #        print("p", i+3)
            #print("NEXT NUC: ", nextNuc,"at pos", i,"+")
            off=offset+start%len(freqs[list(freqs.keys())[0]])
            #print("To remaster ", miniseq,"|",nextNuc,"at", start, off)        
            #(rminiseq, rminigain)=remaster2allCombs(miniseq,freqs,off, nextNuc) 
            (rminiseq, rminigain)=remasterAllCombs(miniseq,freqs,off, nextNuc, use) 
            #print(miniseq,"-->",rminiseq, "\t", rminigain, off)
            maxSeq+=rminiseq
            maxGain+=rminigain
            i+=3
        #print("Seq final", len(maxSeq))
        mseq[str(k)+"_remastered (gain="+str(maxGain)+")"]=maxSeq
    return mseq

#%% Equivalen to remaster2sal for generic di or trinucleotide frequency matrices
def remaster3syn(seqs,freqs,offset=0, aminos=["S", "R", "L"]):
    import bioseq
    mseq={}
    ct=bioseq.codonTable()
    for k in seqs.keys():
        seq=seqs[k]
        maxSeq=""
        maxGain=0
        i=0
        while(i<len(seq)-2):
            codon=seq[i:i+3].upper()
            miniseq=codon
            start=i
            if(i<len(seq)-5):
                nextCodon=seq[i+3:i+6]
                while len(nextCodon)!=0 and (ct[nextCodon] in aminos)==True and (len(miniseq)<33) and i<(len(seq)-5):
                    i+=3
                    miniseq+=nextCodon
                    nextCodon=seq[i+3:i+6]
            nextNuc=""
            if(i+4<len(seq)):
                    nextNuc=seq[i+3] #new
                    if(len(list(freqs.keys())[0])==3):
                        nextNuc+=seq[i+4]
            #        print("p", i+3)
            #print("NEXT NUC: ", nextNuc,"at pos", i,"+")
            off=offset+start%len(freqs[list(freqs.keys())[0]])
            #print("To remaster ", miniseq,"|",nextNuc,"at", start, off)        
            (rminiseq, rminigain)=remasterAllCombs(miniseq,freqs,off, nextNuc) 
            #print(miniseq,"-->",rminiseq, "\t", rminigain, off)
            maxSeq+=rminiseq
            maxGain+=rminigain
            i+=3
        #print("Seq final", len(maxSeq))
        mseq[str(k)+"_remastered (gain="+str(maxGain)+")"]=maxSeq
    return mseq


#%%
def remaster2(seqs,freqs,patchLength=18):
    mseq={}
    for k in seqs.keys():
        seq=seqs[k]
        maxSeq=""
        maxGain=0
        for i in range(0,len(seq),patchLength): #30 nt patches
            miniseq=seq[i:min(i+patchLength,len(seq))]
            (rminiseq,rminigain)=remaster2allCombs(miniseq,freqs)
            maxSeq+=rminiseq
            maxGain+=rminigain
    return mseq
    
#c=remaster2({1:"ATTGATTGTCGTCGTAACGT"}, fr)
#seq="TACCAAGAACCTCTTTTTTGCTTGGATCGAAATTAAAGGTTTAAAAGCAAAGTTATGTCGGACATCGCCTTGAAAACGTACACGGAGCGCGCCAATGTGCATCCTAATGCAGTCGCCAAGAAGTTGCTGCGTTTGATGGACGAGAAGAAGTCAAACCTCTCGGTCGCTGTGGACTTGACGAAGAAGAACCAAGTCCTGGAACTCGTTGACAAGATTGGTCCCTCGATCTGCTTGTTGAAGACTCACATTGACATTGTGGAGGACTTTGACGCAGACATGGTCCAACAGCTGGTTGCTCTGGCTGAGAAGCACAAGTTTTTGATTTTTGAGGATCGCAAGTTTGCCGACATTGGCAACACTGTCAAGCTTCAGTACTCTGCTGGTGTGTACAAGATTGCTTCTTGGGCTGACATCACCAACTGCCACACTGTTCCCGGTGAGGGTATTATTAGCGGTTTGAAGGAAGTTGGTCTTCCCTTGGGTCGTGGTTTGTTGCTTTTGGCTGAGATGTCTTCGAAGGGAACCCTGGCCACTGGAAGCTACACACAAGCCACATTGGAATTGGCTGAGAAGCACAACGACTTCTGTATGGGTTTTATCGCCAGACGCCGCTTCCCCGGTCTCAAGAGCGACTTTATTCACATGACACCCGGTGTCGGTTTGGACGTTAAGGGCGATGGCCTTGGTCAACAATACCGTACACCAGAGGAAGTGATCTGTGAGAGCCAGAGCGATATCATTATCGTCGGTCGCGGTGTCTACGGCAGCGGCCGTGATGCTGCTCAAGAAGCTGAGCGCTACAGAAAGGCCGGCTGGGAGGCCTACCAGCGTCGCATTTCCAAGCAGTAAAAAAAGACTAATGTAAAATTTTTTTGGTTGGTTATTGAAAAAGTCGATGCCTTGTTTGCGTTTGTTTTC"
#seqs={"1":seq}
#ccnl=remaster2nolink(seqs, fr)
#import time
#t0=time.clock()
#cc18=remaster2(seqs, fr,18)
#print("remaster took", time.clock()-t0)
#t0=time.clock()
#cc30=remaster2(seqs, fr,30)
#print("remaster took", time.clock()-t0)

#%%
"""Remastering from nucleotide frequencies"""
def remaster1(seqs,freqs,offset=0):
    import bioseq
    ct=bioseq.codonTable()
    ict=bioseq.inverseCodonTable()
    mseq={}
    flen=len(freqs["A"])
    for k in seqs.keys():
        gain=0
        seq=seqs[k]
        temp=[]
        for i in range(0,len(seq)-2,3):
            codon=seq[i:i+3].upper()
            codons=ict[ct[codon]]
            freqc={}
            ixf=(i+offset)%(flen)
            #print(i, ixf)
            #print(int(ixf/3)+1, codon)
            for c in codons:
                freqc[c]=freqs[c[0]][ixf]+freqs[c[1]][ixf+1]+freqs[c[2]][ixf+2]
                #print("\t", c,freqc[c])
            selcod=max(freqc, key=freqc.get)
            temp+=selcod
            
            freqorig=freqs[codon[0]][ixf]+freqs[codon[1]][ixf+1]+freqs[codon[2]][ixf+2]
            gain+=freqc[selcod]-freqorig
        mseq[k+"_remastered (gain="+str(gain)+")"]="".join(temp)
        
    return mseq

#seq2={"s70":list(seq.values())[0][459:459+153]}
#res1=remaster1(seq2, tal["freq"][1], 0)
#%%
""" Remaster multiplexer
"""
def remaster(seqs, freqs, offset=0, usePath=None):
    k0=list(freqs.keys())[0]
    if(len(k0)==3):
        #return remaster3(seqs, freqs, offset)
        return remaster3syn(seqs, freqs, offset)
    elif(len(k0)==2):
        print("Remastering by dinucleotides")
        use=None
        if(usePath!=None):
            import bioio
            use=bioio.readCodonUsage(usePath, nakamuraFormat=False)
        #return remaster2nolink(seqs,freqs, offset)
        return remaster2sal(seqs,freqs, offset, use=use)
        #return remaster2(seqs,freqs, 9)
    elif(len(k0)==1):
        return remaster1(seqs,freqs, offset)
    else:
        print("Remaster: error, invalid frequency keys")
    return None
#seq="TACCAAGAACCTCTTTTTTGCTTGGATCGAAATTAAAGGTTTAAAAGCAAAGTTATGTCGGACATCGCCTTGAAAACGTACACGGAGCGCGCCAATGTGCATCCTAATGCAGTCGCCAAGAAGTTGCTGCGTTTGATGGACGAGAAGAAGTCAAACCTCTCGGTCGCTGTGGACTTGACGAAGAAGAACCAAGTCCTGGAACTCGTTGACAAGATTGGTCCCTCGATCTGCTTGTTGAAGACTCACATTGACATTGTGGAGGACTTTGACGCAGACATGGTCCAACAGCTGGTTGCTCTGGCTGAGAAGCACAAGTTTTTGATTTTTGAGGATCGCAAGTTTGCCGACATTGGCAACACTGTCAAGCTTCAGTACTCTGCTGGTGTGTACAAGATTGCTTCTTGGGCTGACATCACCAACTGCCACACTGTTCCCGGTGAGGGTATTATTAGCGGTTTGAAGGAAGTTGGTCTTCCCTTGGGTCGTGGTTTGTTGCTTTTGGCTGAGATGTCTTCGAAGGGAACCCTGGCCACTGGAAGCTACACACAAGCCACATTGGAATTGGCTGAGAAGCACAACGACTTCTGTATGGGTTTTATCGCCAGACGCCGCTTCCCCGGTCTCAAGAGCGACTTTATTCACATGACACCCGGTGTCGGTTTGGACGTTAAGGGCGATGGCCTTGGTCAACAATACCGTACACCAGAGGAAGTGATCTGTGAGAGCCAGAGCGATATCATTATCGTCGGTCGCGGTGTCTACGGCAGCGGCCGTGATGCTGCTCAAGAAGCTGAGCGCTACAGAAAGGCCGGCTGGGAGGCCTACCAGCGTCGCATTTCCAAGCAGTAAAAAAAGACTAATGTAAAATTTTTTTGGTTGGTTATTGAAAAAGTCGATGCCTTGTTTGCGTTTGTTTTC"
#seqs={"1":seq}
##print("REMASTER 3")
#res3=remaster(seqs, tal["freq"][3], 0)
#print("REMASTER 2")
#res2=remaster(seqs, tal["freq"][2], 0)
#print("REMASTER 1")
#res1=remaster(seqs, tal["freq"][1], 0)
#%%
""" Returns the dinucleotide gain in freqs for seq1 respect to seq0
"""
def gain(seq1, seq0, freqs):
    if(len(seq1)!=len(seq0)):
        return "Error: sequences of the same length required"
    gain=0
    for i in range(len(seq1)-1):
        d1=seq1[i:i+2]
        d0=seq0[i:i+2]
        fi=i%len(list(freqs.keys())[0])
        gain+=freqs[d1][fi]-freqs[d0][fi]
    return gain
#gain("ATC", "ACC", fr)
    
#%%
def randomChars(length):
    import random
    s=""
    for i in range(length):
        r=-1
        while(r<48 or (r>58 and r<65) or (r>90 and r<97)):
            r=random.randint(48,122)
        s+=chr(r)
    return s
randomChars(20)
    
    