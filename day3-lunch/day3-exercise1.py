#!/usr/bin/env python3

#location of .gtf file of Drosophila Melanogaster genome
# ~/data/genomes/BDGP6.Ensembl.81.gtf

import sys


num_proteincoding_genes = 0

for i, line in enumerate( open( sys.argv[1] ) ):
    if i <= 5:
        continue
    fields = line.rstrip("\r\n").split("\t")
    if "gene" in fields[2] and "protein_coding" in fields[8]:
        num_proteincoding_genes += 1
        
print(num_proteincoding_genes)
        
