#!/bin/bash

GENOME=../genomes/BDGP6
GUIDE=../genomes/BDGP6.Ensembl.81.gtf
SAMPLE=~/Users/cmdb/data/rawdata/

for * in $SAMPLE 
do
  mkdir ~/Users/cmdb/qbb2018-answers/day1-homework/$SAMPLE
  fastqc ~/Users/cmdb/qbb2018-answers/day1-homework/$SAMPLE.fastq
  
done

#hisat2 -p 8 -x BDGP6 -U ../day1-homework/m10. -S m10.sam 
