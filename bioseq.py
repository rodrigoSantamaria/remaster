# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 09:58:56 2017

@author: rodri
"""

#%%---------------------- INVERSA COMPLEMENTARIA ----------------------
#Retorna la cadena (inversa) complementaria a seq
def complement(seq):
    comp=""
    complementos={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}#tabla de conversión
    for s in seq.upper()[::-1]:                 #recorrer en sentido inverso
        comp=comp+complementos[s]
    return comp

#%%--------------------- INVERSA ------------------------------
def inverted(seq):
    inv=""
    for s in seq.upper()[::-1]:
        inv=inv+s
    return inv

#------------------------------- ALIGNMENTS -------------------------------
#%% Get the consensus sequence for a dictionary of sequences
#The consensus letters are shown in lowercase if their frequency is in the range [lower,upper)
#and in upper case if are above the upper parameter. If below the 'lower' parameter noconsensus character is given
def consensus(seqs, lower=.5, upper=.99, noconsensus="-"):
    match=[]
    lseq=len(list(seqs.values())[0])
    for j in range(lseq):
        freq=dict(A=0,G=0,T=0,C=0,N=0)
        for s in seqs.values():
            freq[s[j]]+=1
       # print(freq)
       # add=max(freq.iteritems(), key=operator.itemgetter(1)) #not working in python 3
        add=max(freq, key=(lambda k:freq[k]))
        maxf=max(freq.values())/len(seqs.values())
        if(maxf>upper):
            match.append(add.upper())
        elif(maxf>lower):
            match.append(add.lower())
        else:
            match.append(noconsensus)
    return ''.join(match) #return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk

#motifs must be a dictionary with sequences as values
#Returns the number of nucleotides that differ from the consensus
def score(motifs):
    s=0
    con=consensus(motifs)
    for m in motifs.values():
        for j in range(len(m)):
            if(con[j].upper()!=m[j].upper()):
                s+=1
    return s
    
#score(alnSeqs)
#cs=consensus(remseqs)        
"""Given an dict of aligned sequences (aln), returns the pairwise identities
of the sequences, without counting gaps (-) or not-known nucleotides (N) on
any of the sequences. 
"""
def identity(aln):
    ret={}
    for i in range(len(aln)):
        k1=list(aln.keys())[i]
        s1=list(aln.values())[i]
        for j in range(i+1, len(aln)):
            k2=list(aln.keys())[j]
            s2=list(aln.values())[j]
            s=0
            slen=0
            for p in range(len(s2)):
                if( (s1[p].upper()!="N" and s1[p]!="-") and (s2[p].upper()!="N" and s2[p]!="-") ):
                    slen=slen+1
                    if(s1[p].upper()==s2[p].upper()):
                        s+=1
            #print(k1+", "+k2+"\t"+str(round(100*s/slen,2)))
            if(slen>0):
                ret[k1+","+k2]=s/slen
            else:
                ret[k1+","+k2]=0
    return ret
    
#tal=identity(alnSeqs)
    


#%%------------------------ CODON MANIPULATION -------------------------
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
    for aa in np.unique(list(ct.values())):
        ict[aa]=list(filter(lambda k:ct[k]==aa, list(ct.keys())))
    return ict

#Takes a nucleotide sequence and alters it so it contains the same aminoacid
#sequence but mutates randomly across the possible codons
def randomCodons(seq):
    import random
    ct=codonTable()
    ict=inverseCodonTable()
    rseq=""
    for i in range(0,len(seq),3):
        cod=seq[i:i+3]
        rcod=random.sample(ict[ct[cod]],1)[0]
        rseq+=rcod
    return rseq

#%% ---------------------- COMPARISON ----------------------------
#Returns the  ofhomology between two nucleotide sequences s1 and s2 as 
#a dictionary:
#   aa - amino acid homology (considering a 0-frame start)
#   codon - codon homology (considering a 0-frame start)
#   nt - nucleotide homology
def coincidence(s1,s2):
    co={}
    ct=codonTable()
    ict=inverseCodonTable()
    co["codon"]=sum([1 if s1[i:i+3]==s2[i:i+3] else 0 for i in range(0,len(s1),3)])/(len(s1)/3.)
    co["nt"]=sum([1 if s1[i]==s2[i] else 0 for i in range(len(s1))])/(1.*len(s1))
    co["aa"]=sum([1 if (s1[i:i+3] in ict[ct[s2[i:i+3]]])==True else 0 for i in range(0,len(s1),3)])/(len(s1)/3.)
    return co
    
    
#------------------- IUPAC AMBIGUOUS CODES -------------------------
#%% Given an IUPAC sequence, it unfolds it on all its possible nucleotide combinations
def unfold(seq):
   from Bio import Seq
   import itertools
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   r = []
   for i in itertools.product(*[d[j] for j in seq]):
      r.append("".join(i))
   return r
#%%
"""
Given an expression with degenerate nucleotides, it provides a regular expression
for searches
seq - sequence with nucleotides (including degenerate ones (B-Y)
returns a regular expression where each nucleotide has been substituted by
        all the nucleotides it comprises (e.g. AWR --> A[ATW][AGR])
"""
def nuc2exp(seq):
    #mapper according to http://opisthokonta.net/?p=549
    nmap={"B":"[CGTBSKY]",
         "D":"[AGTDRWK]",
         "H":"[ACTHMYW]",
         "K":"[GTK]",
         "M":"[ACM]",
         "N":"[ACGTBDHKMNRSVWY]",
         "R":"[AGR]",
         "S":"[CGS]",
         "V":"[ACGVMSR]",
         "W":"[ATW]",
         "Y":"[CTY]",
         "A":"A",
         "C":"C",
         "G":"G",
         "T":"T",
         }
    exp=""
    for s in seq:
        exp+=nmap[s]
    return exp



#%% ------------------------- MUTATIONS ---------------------
# Dada una secuencia (word), un número de mutaciones puntuales (num_mismatches)
# y una cadena con los posibles componentes (letters), devuelve una variable 
# generator con todas las posibles mutaciones. 
# OJO: el resultado hay queconvertirlo a list para su porterior uso
#Código extraido de: 
#http://stackoverflow.com/questions/11679855/introducing-mutations-in-a-dna-string-in-python
def mutations(word, num_mismatches, letters="ACGT"):
    import itertools
    for locs in itertools.combinations(range(len(word)), num_mismatches):
        this_word = [[char] for char in word]
        for loc in locs:
            orig_char = word[loc]
            this_word[loc] = [l for l in letters if l != orig_char]
        for poss in itertools.product(*this_word):
            yield ''.join(poss)
            
            
#%% ------------------------- MUTACIONES (II) ---------------------
# Como la función anterior, pero ahora retorna una lista con las mutaciones 
#posibles con num_mismathes SNPs *o menos*
def mutationsEqualOrLess(word, num_mismatches, letters="ACGT"):
   matches=set()
   for dd in range(num_mismatches,-1,-1): 
       matches.update(list(mutations(word, dd, letters)))
   return matches

#%% Pruebas, cuántas mutaciones de 2/3 SNPs hay en un 9-mer?
# kmer='CACAGTAGGCG'
# mu=list(mutations(kmer, 2))
# print len(mu)
# mu=list(mutations(kmer, 3))
# print len(mu)
# mu=mutationsEqualOrLess(kmer, 3)
# len(mu)

#%% ------------------------- MUTACIONES (III) --------------------- S2E5
# Retorna como lista todos los posibles k-mers en seq y todas sus mutaciones 
#en hasta 2 SNPs
def allMutations(seq, k=9, d=2):
    kmers=set()
    for i in range(len(seq)-k+1):
        kmer=seq[i:i+k]
        kmers.update(mutationsEqualOrLess(kmer, d))
    return kmers


#%% Mutates a nucleotide sequence seq in k positions.
# Only A, T,C, G allowed as characters in the sequence.
# If repeatPos is True, some of these positions may be mutated several times (default False). 
#If repeatNuc is True, a nucleotide may be mutated into itself (default False)
# w sets the weights or probability (as number adding up to 100) of choosing each
# nucleotide.
def mutate(seq, k=1, repeatPos=False, repeatNuc=False, w={"A":25, "C":25, "G":25, "T":25}):
    import random 
    pos=list(range(len(seq)))
    nuc=["A","T","C","G"]
    ww=[]
    for n in nuc:
        ww.append(w[n])
    rseq=list(seq.upper())
    for i in range(k):
        #select a random position
        rpos=random.choice(pos)
        if not repeatPos:
            pos.remove(rpos)
        #select a random (maybe weighted) nucleotide
        ant_nuc=rseq[rpos]
        rseq[rpos]=random.choices(nuc, weights=ww, k=1)[0]
        if not repeatNuc:
            while(ant_nuc==rseq[rpos]):
                rseq[rpos]=random.choices(nuc, weights=ww, k=1)[0]
                
            
    return "".join(rseq)        
#tal=mutate("TTT", k=3, repeatNuc=True)   
#print(tal)     
    

#%%
def atcontent(seq):
    return seq.count("A")+seq.count("T")