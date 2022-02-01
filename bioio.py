# -*- coding: utf-8 -*-
"""
Methods for I/O of common bioinformatics formats

@author: rodri
"""
#%%
######################### FASTA ####################################
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

#tal=readFasta(seqFile)
#tal=readFasta("/home/rodri/data/remaster/ura4/ura4sj-6nuc/ura4sj-6nuc.fasta")
#[len(x) for x in tal.values()]
#tal=readFasta("/home/rodri/data/annotations/sp/Spombe.fasta")
#[len(x) for x in tal.values()]
#%%

"""
path - path for the new fasta file
seqs - dictionary with sequences. keys will be used for description
        it can also be a list of sequences, in which case the sequence ids will be numbers
"""  
def writeFasta(seqs, path):
    f=open(path,"w")
    if(type(seqs)==dict):
        for k in seqs:
            f.write(">"+k+"\n")
            f.write(seqs[k]+"\n")
    else:
        for k in range(len(seqs)):
            f.write(">"+str(k)+"\n")
            f.write(seqs[k]+"\n")
    f.close()
    

#%%  
######################### BED ####################################
"""
path - path to the bed file
returns: a dictionary with the chromosomes as keys and arrays of starting 
        positions as values.
"""
def readBedStarts(path):
    import csv
    csvfile=open(path)
    reader=csv.DictReader(csvfile, fieldnames=['chr','start','end'], delimiter="\t")
    
    seqs={}
    for row in reader:
        ch=row["chr"]
        if(not ch in seqs.keys()):
            seqs[ch]=[]
        seqs[ch].append(int(row["start"]))   
    return seqs


def writeBedFromLists(path, track, start, end, strand=None):
    if(len(end)!=len(start) or len(track)!=len(start)):
        return "ERROR: start, end and track lists must have the same length"
    f=open(path,"w")
    if(strand==None):
        for i in range(len(track)):
            f.write(str(track[i])+"\t"+str(start[i])+"\t"+str(end[i])+"\n")
    else:
        for i in range(len(track)):
            f.write(str(track[i])+"\t"+str(start[i])+"\t"+str(end[i])+"\t"+str(strand[i])+"\n")
    f.close()
    return

#save to file:
# - path for the file to be saved
# - tracks a dictionary where keys are chromosomes/tracks and for each entry
#    has a list of genomic positions for that track with the following keys:
#       start, end, name, score, strand
#       The start and end keys must be numeric and are mandatory (int)
#       The score key may be missing and if present must be numeric (int or float)
#       The name and strand keys may be missing and must be str   
def writeBed(path, tracks):
    #if(len(end)!=len(start) or len(track)!=len(start)):
    #    return "ERROR: start, end and track lists must have the same length"
    f=open(path,"w")
    for ch in tracks.keys():
        for p in tracks[ch].values():
          #  f.write(str(ch)+"\t"+str(p["start"])+"\t"+str(p["start"]+150)+
             f.write(str(ch)+"\t"+str(p["start"])+"\t"+str(p["end"])+
                    "\t"+(str(p["name"]) if "name" in p.keys() else "")+
                    "\t"+(str(p["score"]) if "score" in p.keys() else "")+
                    "\t"+(str(p["strand"]) if "strand" in p.keys() else "")+"\n")
    f.close()
    return

def writeBedFromDict(path, track):
    f=open(path,"w")
    for v in track.values():
        f.write(str(v["chromosome"])+"\t"+str(v["start"])+"\t"+str(v["end"])+"\t"+str(v["name"])+"\t"+str(v["score"])+"\t"+str(v["strand"])+"\n")
    f.close()
    return

 #%%  
######################### BED ####################################
"""
path - path to the bed file
checkNames - if set, it checks for names in column 4 as for the canon format (default True)
header - if set, it skips the first line (default False)

returns: a dictionary with the names as keys and arrays of dictionaries as values (start and end as ints, value as float, rest as string)
"""
def readBed(path, checkNames=True, header=False):
    import csv
    csvfile=open(path)
    reader=csv.DictReader(csvfile, fieldnames=['chromosome','start','end', 'name', 'score','strand'], delimiter="\t")
    cont=0
    seqs={}
    for row in reader:
        if(cont==0 and header==True):
            cont+=1
            continue
        if(checkNames==False or (row["name"]=="." or row["name"]=="" or row["name"]==None)):
            if(cont==0 and checkNames==True):
                print("WARNING: row names with no value found. Such values will be converted to numbers starting in 0 (0,1,2,3...)")
            row["name"]=str(cont)
            cont+=1
        elif(row["name"] in seqs.keys()):
            print("WARNING: entry for: ", row["name"], "is duplicated, last one in the file will be returned")
        seqs[row["name"]]={"start":int(row["start"]), "end":int(row["end"]), "chromosome":row["chromosome"], "strand":row["strand"], "score": (row["score"])}   
    return seqs

 
#%%
######################### WIG ####################################
"""
path - path to the wig file (str)
method - in the case of variableStep wigs, which method to use to interpolate values (default "slope")
        'step' simply gives  value seq[i]=v to each point i+n up to next step
        'slope' assigns a value in the line between the values of a step and the next
returns: a dictionary with the chromosome names inf the fireld chrom as keys 
           and the count levels as values (numpy.array of numpy.float32)

requires: re for chromosome name parsing
          numpy for array building
          time for performance messages
COMMENTS: changes to float16 if memory issues, to float64 if oveflows (inf values)
TODO: variableStep requires mandatory track lines (but the format says they are optional)
""" 
def readWig(path,method="slope"):
    import numpy
    import time
    import re
    
    t0=time.clock()
    print(path)
    f=open(path)
    seq=f.readlines()
    print ((time.clock()-t0),' s in reading') 
    t0=time.clock()
    print(seq[1])
    if("fixedStep" in seq[1] or "fixedStep" in seq[0]):
        print ("Fixed Step")
        start=[]    #starting and ending lines per track
        end=[]
        names=[]
        
        tline=0 #marks if there is a track line
        for i in range(len(seq)):
            s=seq[i]
            if s[0]=='f':#new track
             start.append(i+1)
             name=re.sub("\n", "", re.sub(" .*$", "", re.sub("^.*chrom=", "", s)))
             names.append(name)
             if(i>2):
                 end.append(i-tline)
             tline=0
            elif s[0]=="t":#track optional line
             tline=1
        end.append(len(seq))
        print((time.clock()-t0),' s in computing sizes', start,end,names)
        
        t0=time.clock()
        ch={}
        for i in range(len(start)):
            ch[names[i]]=numpy.array(seq[start[i]:end[i]], dtype=numpy.float32)
            #float16 gives issues with jsonify but is more memory efficient
            #float32, 64 is better to avoid overflows but less mem efficient
        print ((time.clock()-t0),' s in formatting')
        
        return ch
    else:   #TODO: requires track line in all the entries
        print("Variable step")
        chsize=[]
        for i in range(2,len(seq)):
            s=seq[i]
            if "variableStep" in seq[i]:
             chsize.append((int)(seq[i-2].split("\t")[0]))
        chsize.append((int)(seq[i-2].split("\t")[0].strip()))
        print((time.clock()-t0),' s in computing sizes')
        
        
        t0=time.clock()
        ch=interpolate(seq,chsize, method)
      
        print((time.clock()-t0),' s in <<interpolating>> seqs')
        return ch

#%%writeWig
"""
path - file to store the data in wig format
tracks - dictionary with chromosomes as keys.
NOTE: by now, this is a simple func
"""
def writeWig(path, tracks, fixed=True, start=1, step=1, span=1):
    f=open(path, "w")
    for k in tracks.keys():
        header="fixedStep chrom="+k+" start="+str(start)+" step="+str(step)+" span="+str(step)+"\n"
        f.write(header)
        for x in tracks[k]:
            f.write(str(x)+"\n")
    f.close()        
        
#%% - Interpolation
"""
seq - sequence lines as read from a wig variabale step, withouth initial comments
chsize - wig track sizes
method - interpolation method. Default 'step'
        'step' simply imputes value seq[i]=v to each point i+n up to n=next variableStep value -1
        'slope' imputes in the same ranges than 'step' values in the line between v1 and v2
"""    
def interpolate(seq, chsize, method="step"):
    import numpy, re
    ch={}
    cont=0
    if(method=="step"):
        for i in chsize:
            cont+=1
            print(seq[cont])
            name=re.sub("\n", "", re.sub(" .*$", "", re.sub("^.*chrom=", "", seq[cont])))
            cont+=1
            print(i, name)
            chi=numpy.empty(i, dtype=numpy.float32) #gives issues with jsonify but is much more memory efficient
            index=0
            while((cont+2< len(seq)) and (("variable" in seq[cont+2]) == False)):
                s=seq[cont].split("\t")
                level=(float)(s[1])
                if("variable" in seq[cont-1]):
                    nuc=0
                    nuc2=(float)(s[0])
                else:
                    nuc=(int)(s[0])
                    nuc2=(int)(seq[cont+1].strip().split("\t")[0])
                step=nuc2-nuc
                chi[nuc:nuc2]=[level]*step
                index+=step
                cont+=1
            cont+=1
            ch[name]=chi
    elif(method=="slope"):
        for i in chsize:
            cont+=1
            print(seq[cont])
            name=re.sub("\n", "", re.sub(" .*$", "", re.sub("^.*chrom=", "", seq[cont])))
            cont+=1
            print(i, name)
            chi=numpy.empty(i, dtype=numpy.float32) #gives issues with jsonify but is much more memory efficient

            y0=0
            x0=0
            while((cont+2< len(seq)) and (("variable" in seq[cont+2]) == False)):
                s=seq[cont].split("\t")
                x1=(int)(s[0])
                y1=(float)(s[1])
                
                f=(y1-y0)/(x1-x0)
                for x in range(x0,x1):
                    chi[x]=f*(x-x0)+y0
                    
                cont+=1
                x0=x1
                y0=y1
            cont+=1
            ch[name]=chi        
    return ch;
    
#%%    
############################# GO ###################
"""
path - path to a GO obo file 
returns : a dict where keys are go_ids and values are just go names by now
"""
def readGO(path):
    f=open(path)
    lines=f.readlines()
    data={}
    i=0
    while i<len(lines):
        if(lines[i].startswith("[Term]")):
            go_id=lines[i+1].replace("id:", "").strip()
            go_name=lines[i+2].replace("name:", "").strip()
            data[go_id]=go_name
        i+=1
    return data
    
    
#%%    
############################# GOA ###################
"""
path - path to a GOA  file 
returns : an array of dicts. 
  Each dict contains "gene_id", "gene_name", "go_id", "evidence", "go_type" and "gene_desc"
   which corresponds to columns 1,2,4,6,8 and 9 from the GOA file
""" 
def readGOA(path):
    f=open(path)    
    lines=f.readlines()
    data=[]
    for l in lines:
        if(l.startswith("!")==False):
            vals=l.split("\t")
            data.append({"gene_id":vals[1],"gene_name":vals[2], "go_id":vals[4], "evidence":vals[6], "go_type":vals[8], "gene_desc":vals[9]})
    return data
#NOTE: there's a bunch of utils for GO/GOA in seqview/annotations.py which might be useful
#%%
############################# GFF ###################
"""
path - path to a GFF(3)  file 
structure - Leave to 'single' (legacy or recoverable method to read several gffs at ones for large organisms. See NOTE) 
organism - Legacy parameter to fix deviations from standards in some organisms
Returns a np table ("single") or a dictionary with np tables (structure="multiple") 
with genome functional annotations.
NOTE: for larger GFF files (e.g. dmelanogaster, 600MB) this is not performing at all)
"""
def readGFF(path, structure="single", organism=""):
    import csv
    import time
    
    #print("reading GFF...")
    f=open(path)
    regions=["gene", "transcript", "exon","ncRNA_gene", "tRNA_gene", "snRNA_gene", "snoRNA_gene", "rRNA_gene", "three_prime_UTR", "five_prime_UTR", "CDS", "ORF"]
    fieldnames=["seqid", "source", "type", "start", "end", "score", "sense", "phase", "attributes"]
    reader=csv.DictReader(f, fieldnames=fieldnames, delimiter="\t")
    
    t0=time.clock()
    entries=[]
    
    #print("filtering out comments...")
    for row in reader:
        if(row['seqid'].startswith("#") or not (row["type"] in regions)):   #case with comments between entries. NOTE: tam will be miscalculated in these cases
                continue
        entries.append(row)
        
    #print("time in filtering out entries", (time.clock()-t0))
    t0=time.clock()
 
    import numpy as np
    import re
    #print("populating annotations...")
    sc=("cerevisiae" in organism)
    ce=("elegans" in organism)
    
    #data=np.empty(len(entries),dtype=[("chromosome", "a40"),("type", "a30"), ("start", "i8"), ("end", "i8"), ("sense", "a1"), ("id", "a40"), ("name", "a40")])
    #PY3 change: aX no longer retrieves basic strings, but byte strings
    data=np.empty(len(entries),dtype=[("chromosome", str,40),("type", str,30), ("start", "i8"), ("end", "i8"), ("sense",str,1), ("id", str,40), ("name", str,40)])
                
    for i in range(len(entries)):
        row=entries.pop()
        
        data[i]["start"]=(int)(row["start"])
        data[i]["end"]=(int)(row["end"])
        if(data[i]["start"]>data[i]["end"]):#some GFF erroneusly switch start and end in antisens
            temp=data[i]["start"]
            data[i]["start"]=data[i]["end"]
            data[i]["end"]=temp
        
        seqid=row["seqid"]
        data[i]["chromosome"]=seqid
        
        data[i]["type"]=row["type"]
        data[i]["sense"]=row["sense"]
        if(sc==True):#SGD are 'flexible' about standards... grrr
            gid=re.sub("gene:", "", re.sub(";Name.*$","",re.sub("^.*ID=", "", row["attributes"])))
            name=re.sub("^.*ID=", "", re.sub(";.*$","",re.sub("^.*gene=", "", row["attributes"])))
        elif(ce==True):#WB are 'flexible' about standards... grrr
            gid=re.sub("\"","",re.sub("gene_id", "", re.sub(";.*$", "", row["attributes"]))).strip()
            name=re.sub("\"","",re.sub("gene_id", "", re.sub(";.*$", "", row["attributes"]))).strip()
        else:
            gid=re.sub("gene:", "", re.sub(";Name.*$","",re.sub("^.*ID=", "", row["attributes"])))
            name=re.sub("^.*ID=", "",re.sub(";.*$","",re.sub("^.*ID=.*;Name=", "", row["attributes"])))
        
        data[i]["id"]=gid
        data[i]["name"]=name            
        
    if(structure=="single"):
        return data
    else:
        d={}
        wanted_set = set(data["chromosome"])  # Much faster look up than with lists, for larger lists
        for k in wanted_set:
            @np.vectorize
            def selected(elmt): return elmt in k  # Or: selected = numpy.vectorize(wanted_set.__contains__)
            d[k]=data[selected(data["chromosome"])]

        return d
        
#tal=readGFF(gffPath, organism="Saccharomyces cerevisiae")
#%%Reads a tab delimited tab with the following columns:
#codon usage
#Usage as percentage. Can be used Nakamura 'format', which has column 4 as usage, in occurrences per 1000 
def readCodonUsage(path, nakamuraFormat=False):
    f=open(path)
    cu={}
    if(nakamuraFormat):
        for line in f.readlines():
            line=line.split("\t")
            cu[line[0]]=float(line[3])/10
    else:
        for line in f.readlines():
            line=line.split("\t")
            cu[line[0]]=float(line[1])
    return cu

#%% Read .bic files as returned by BicOverlapper
def readBiclusters(path):
    f=open(path)
    grupos={}
    lines=[x.strip() for x in f.readlines()[2:]]
    for i in range(0,len(lines),3):
        gname=lines[i].split(":")[0]
        rows=lines[i+1].split("\t")
        cols=lines[i+2].split("\t")
        
        grupos[gname]={'rows':rows, 'cols':cols}    
    f.close()
    return grupos


#%% 
"""
Read meme matrix files. The algorithm is a bit ad-hoc for the meme files we 
use, might not be useful for others.
"""
def readMeme(path):
    f=open(path)
    fm={}
    import re
    for l in f.readlines():
        if(l.startswith("MOTIF")):
            title=re.sub(" ", "_", re.sub(":.*","", l)).strip()
            fm[title]=[]
        elif((l.startswith("0") or l.startswith("1"))):
            fm[title].append(l.strip())
    return fm

#%%
""" Reads a XXmotif output sequence file, returning a dict with motifs as keys
and sequences (column site) as values in a list
"""
def readMotifSeqs(path):
    f=open(path)
    mot={}
    key=""
    import re
    for l in f.readlines():
        if(l.startswith("Motif")):
            title=l.split("\t")[0]
            title=title.split(":")
            consensus="".join(title[1].strip().split(" "))
            key=title[0].strip()+" "+consensus
            mot[key]=[]
        elif(l.startswith("     ")):
            mot[key].append(l.split("\t")[1].strip())
    return mot

"""----------------------------------------- CONVERTERS -----------------"""
#%%
def bamm2meme(bammFolder="", memeFile=""):
    import re
    if(memeFile==""):
        memeFile=bammFolder+".meme"
    
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(bammFolder) if isfile(join(bammFolder, f))]
    ihbcps=sorted([bammFolder+"/"+f for f in onlyfiles if f.endswith(".ihbcp")])
    print(ihbcps)
    f=open(memeFile, "w")
    f.write("MEME version 4\n\n")
    f.write("ALPHABET=ACGT\n\n")
    f.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
    
    
    for file in ihbcps:
        motname=re.sub("\\..*$", "", re.sub("^.*/", "", file))
        f.write("MOTIF "+motname+"\n")
        f.write("unknown, retrieved from BaMM file "+file+"\n")
        f.write("letter-probability matrix:")
        #f.write("letter-probability matrix: alength=4, w=")
        fm=open(file)
        motmat=[]
        for line in fm.readlines():
            if(line.count(" ")==4):
                motmat.append(line)
        #f.write(str(len(motmat))+"\n")
        f.writelines(motmat)
        f.write("\n")
    f.close()    
    return "meme sucessfuly built at"+memeFile