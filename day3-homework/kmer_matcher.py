#!/usr/bin/env python3

#this code is from /day3-afternoon/03-count_kmers.py

import sys
import fasta

target = fasta.FASTAreader( sys.stdin )
query = fasta.FASTAreader( open(sys.argv[1]) )
k = int(sys.argv[2])

kmers = {} #from query
kmer = 0
kmers[kmer] = []

print(k)

for ident, sequence in query:
    for posn, v  in enumerate(range( 0, len (sequence) - k )):
        kmer = sequence[posn:posn+k]
        if kmer not in kmers:
            kmers[kmer] = [posn]
        else:
            kmers[kmer].append(posn)
        
for ident, sequence in target:
    for i, value in enumerate(range( 0, len (sequence) - k )):
        target_kmer = sequence[i:i+k]
        if target_kmer in kmers:
            print(ident, i, kmers[target_kmer], target_kmer) 



# sorted_kmers = sorted(kmers.items(), key=lambda x:x[1])
# for kmerfortarget, occurances in sorted_kmers
# print(sorted_kmers)

# for kmers in target:
#     print(kmers, kmers[kmer])
#
    
    

