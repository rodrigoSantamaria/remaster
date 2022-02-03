


REMASTER
========
_Remaster_ is a python3 script to:
  1. Compute k-mer frequencies by location on a set of same length sequences 
  2. remaster genomic -coding- sequences by maximizing the frequencies at each location, but respecting aminoacids.

We provide it with some sample data.

Usage example
-------------

Let's say we want to compute the dinucleotide frequencies on the positions along 150bps from nucleosome+1 sequences (_Nuc+1.bed_) coming from _Schizosaccaromyces pombe_ genome (_Spombe.fasta_)

<code>python remaster.py -bed Nuc+1_ok.bed -genome Spombe.fasta -log2 -k 2 -out folder</code>

Figures and the dinucleotide frequency matrix (aka signature, a tsv file) are generated within _folder_

Now, we can incorporate this Nuc+1 signature into any sequence. For example, let's take the 6 identified nucleosomes on gene _ura4_

<code>python remaster.py -freqs folder/Nuc+1_ok-freqs-2.tsv -seqs ura4sp_nucs_147bp_inverted.fasta -out folderURA4</code>

More info about the algorithms underneath at XXX (under review paper) and on the code. More info on the usage: 
<code>python remaster.py -h</code> 
