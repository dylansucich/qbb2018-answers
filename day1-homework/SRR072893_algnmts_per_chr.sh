#!/bin/bash

grep "SRR072893" SRR072893.sam > SRR072893_alignment.sam
cut -f 3 SRR072893_alignment.sam | grep -v "^211" > SRR072893_alignment_cut.sam
datamash -s -g 1 count 1 < SRR072893_algnmts_per_chr_cut.sam > SRR072893_algnmts_per_chr.sam