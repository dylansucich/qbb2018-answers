#!/usr/bin/env python3

"""
Need: velvet, spades, lastz
conda install <Need>

Commands:
velveth Assem 31 -shortPaired -fastq reads_low_1.fastq reads_low_2.fastq
velvetg Assem

spades.py --only-assembler --1 reads_low_1.fastq --2 reads_low_2.fastq -o spadesdir_low

lastz reference.fasta /Users/cmdb/qbb2018-answers/class2/Assem/contigs.fa --format=general:zstart1,size2,end1,name2 --chain --output=velvet_low_dotplot.out
USAGE: ./class2-homework.py /Users/cmdb/qbb2018-answers/class2/spadesdir_low/contigs.fasta 

"""

import sys
import fasta

contigs = fasta.FASTAReader(open(sys.argv[1]))
num_contigs = 0
contig_lengths = []

for ident, sequence in contigs:
    num_contigs += 1
    contig_lengths.append(len(sequence))
    
contig_lengths.sort(reverse=True)
print("Total # of contigs = " + str(num_contigs))
print("Min contig length = " + str(contig_lengths[-1]))
print("Max contig length = " + str(contig_lengths[0]))
print("Avg contig length = " + str(sum(contig_lengths)/num_contigs))

half_length = sum(contig_lengths)/2
counter = 0

for i, length in enumerate(contig_lengths):
    if counter >= half_length:
        print("n50 = " + str(contig_lengths[i-1]))
        break
    else:
        counter += length
        

