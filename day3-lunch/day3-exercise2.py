#!/usr/bin/env python3

#location of .gtf file of Drosophila Melanogaster genome
# ~/data/genomes/BDGP6.Ensembl.81.gtf

import sys


num_proteincoding_genes = 0
gene_biotype_dict = {}

for i, line in enumerate( open( sys.argv[1] ) ):
    if i <= 5:
        continue
    fields = line.rstrip("\r\n").split("\t")
    read_type = fields[2]
    if read_type == "gene":
        fields = line.rstrip("\r\n").split()
        gene_biotype = fields[17]
        
        if gene_biotype in gene_biotype_dict:
            gene_biotype_dict[gene_biotype] += 1
       
        else:
            gene_biotype_dict[gene_biotype] = 1

for name, value in gene_biotype_dict.items():
    print(name, value)
            
