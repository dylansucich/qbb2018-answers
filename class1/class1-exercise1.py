#!/usr/bin/env python3

"""

blastn -db nr -evalue 0.0001 -remote -outfmt "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qseq sseq" -query week1_query.fa -out week1_query.out
input:week1_query.out

awk '{print ">"$1"\n"$4}' week1_query_short.out > seqs.fa

transeq seqs.fa > seqs.pep   

"""
import sys
import os
import pandas as pd
import fasta
import itertools as it


query_file = open(sys.argv[1])
mafft_file = open(sys.argv[2])

unaligned_dna_seq = fasta.FASTAReader(query_file)
aligned_peptide_seq = fasta.FASTAReader(mafft_file)
    
nt_list = []
aa_list = []

for (dna_id, dna), (aa_id, aa) in zip(unaligned_dna_seq, aligned_peptide_seq): # zip parses 2 files at once
      nuclist = []
      aalist = []
      j=0
      aa_list.append(aa)
    
      for i in range(len(aa)):
          a = aa[i]
          nt = dna[j:j+3]
        
          if a == "-":
              nt = "---"
              nuclist.append(nt)
          else:
              j+= 3
              nuclist.append(nt)
        
      nt_list.append(nuclist)
      aa_list.append(aalist)

# print(nt_list)
# print(aa_list)


dS = {}
dN = {}
ds_count = 0
dn_count = 0
query = 0

for codon, aa in zip(nt_list, aa_list):
    for i in range(len(aa)):
        a = aa[i]
        nuc = dna[j:j+3]
        for nuc == "---":
            dN[i] += 1
        elif no != "---":
            dS[i] += 1

print(nt_list)


    
    
# for (query_ident, query_seq), (mafft_ident, mafft_seq) in zip(unaligned_dna_seq, aligned_peptide_seq):
#     new_nseq = ''
#     nid = 0
#     for aa in mafft_seq:
#         if aa == '-':
#             new_nseq += '---'
#         else:
#             new_nseq += query_seq[nid: nid + 3]
#         nid += 3
#     print('>' + query_ident)
#     print(new_nseq)
# AA to NT
# for "-" in AA
#     add "---"
# print(">" + df.iloc[:,2] +"\n" +df.iloc[:,-1])
