#!/usr/bin/env python3

#usage ./day3-exercise3.py <gtf file>
import sys


find_pos = int(21378950)
gene_start = 0
gene_end = 0
my_dist = 0
pc_dict = {}
npc_dict = {}



for i, line in enumerate( open( sys.argv[1] ) ):
    if i <= 5:
        continue
    fields = line.rstrip("\r\n").split("\t")
    read_type = fields[2]
    chrom = fields[0]
    gene_start = int(fields[3])
    gene_end = int(fields[4])
    gene_id = fields[8]
   # gene_name = fields[13]
    
    if chrom == "3R" and read_type == "gene":
        if "protein_coding" in fields[8]:
            fields = line.rstrip("\r\n").split()
            gene_biotype = fields[17]
            
            if find_pos < gene_start:
                my_dist = gene_start - find_pos
                
            elif find_pos > gene_end:
                my_dist = find_pos - gene_end
                pc_dict[gene_id] = my_dist 
            
        else:
            #nonprotein coding    
            fields = line.rstrip("\r\n").split()
            gene_biotype = fields[17]
            if find_pos < gene_start:
                my_dist = gene_start - find_pos
            
            elif find_pos > gene_end:
                my_dist = find_pos - gene_end
                npc_dict[gene_id] = my_dist

# for name, value in pc_dict.items():
#     print(name[9:19], value)
# for name, value in npc_dict.items():
#     print(name[9:19], value)

sorted_pc_dict = sorted(pc_dict.items(), key=lambda x:x[1])
min_line, min_dist = sorted_pc_dict[0]
print("Closest protein-coding gene: " + min_line[9:19], min_dist)

sorted_npc_dict = sorted(npc_dict.items(), key=lambda x:x[1])
min_line, min_dist = sorted_npc_dict[0]
print("Closest non-coding gene: " + min_line[9:19], min_dist)
#
# print(min_line)
# print(sorted_pc_dict[0])
# print(sorted_npc_dict[0])


   
    

