#!/usr/bin/env python3

"""

blastn -db nr -evalue 0.0001 -remote -outfmt "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qseq sseq" -query week1_query.fa -out week1_query.out
input:week1_query.out

awk '{print ">"$1"\n"$4}' week1_query_short.out > seqs.fa

transeq seqs.fa > seqs.pep   

USAGE: ./class1-exercise1.py seqs.fa seqs.mafft

"""
import sys
import os
import pandas as pd
import fasta
import itertools as it
from statsmodels.stats import weightstats as stests
import matplotlib.pyplot as plt


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
    
      for i in range(len(aa)):
          a = aa[i]
          aalist.append(a)
          
          nt = dna[j*3:(j*3)+3]
        
          if a == "-":
              nuclist.append("---")
          else:
              nuclist.append(nt)
              j +=1
        
      nt_list.append(nuclist)
      aa_list.append(aalist)

# print(nt_list)
# print(aa_list)


dS = {}
dN = {}
ds_count = 0
dn_count = 0
query = 0

query_dna = nt_list[0]
query_aa = aa_list[0]

ds_count = [0] * len(query_aa)
dn_count = [0] * len(query_aa)
count = [0] * len(query_aa)



for i in range(1, len(aa_list[1:])):
	align = aa_list[i]
	
	for j in range(len(align)):
		if align[j] != query_aa[j]:
			dn_count[j] += 1
		else:
			ds_count[j] += 1
            
for codon_list, aa_list in zip(nt_list[1:], aa_list[1:]):
     for i in range(0, len(codon_list)):
             if codon_list[i] != query_dna[i]:
                 count[i] += 1
                 if aa_list[i] == query_aa[i]:
                     ds_count[i] += 1
                 else:
                     dn_count[i] += 1
     
ratios = []
for i in range(len(ds_count)):
    ratio = (dn_count[i])/(ds_count[i] + 1)
    ratios.append(ratio)

ztest = stests.ztest(dn_count, ds_count)
print(ztest)

fig, ax = plt.subplots(figsize=(20, 8))
ax.scatter(range(0, len(ratios)), ratios, s = 3)
ax.set_xlabel("AA position")
ax.set_ylabel("dN/dS ratio")
ax.set_title("Ratio of Nonsynonymous and Synonymous Mutations at each Amino Acid Position")
plt.tight_layout()
fig.savefig("Nonsynonymous_to_Synonymous.png")
plt.close(fig)

     

    
    
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
