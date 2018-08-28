#!/bin/bash

GENOME="/Users/cmdb/qbb2018-answers/genomes/BDGP6"
ANNOTATION="/Users/cmdb/qbb2018-answers/genomes/BDGP6.Ensembl.81.gtf"


for SAMPLE in  SRR072893 SRR072903 SRR072905 SRR072915
do
  mkdir $SAMPLE
  cd $SAMPLE
  fastqc ../${SAMPLE}.fastq
  hisat2 -p 8 -x ${GENOME} -U ../${SAMPLE}.fastq -S ${SAMPLE}.sam
  samtools view ${SAMPLE}.sam -b -o ${SAMPLE}.bam
  samtools sort -O BAM -o ${SAMPLE}.sorted.bam ${SAMPLE}.bam
  samtools index ${SAMPLE}.sorted.bam
  stringtie -p 8 -e ${SAMPLE}.sorted.bam -G $ANNOTATION -o ${SAMPLE}_stringtie.gtf -B
  cd ../
done