#!/usr/bin/env python3

# 09 11 23 24 27 31 35 39 62 63
#
# samtools view -bS sequence1.sam > sequence1.bam
#
# freebayes -f ref.fa --haplotype-length 500 \ --haplotype-basis-alleles vars.vcf aln.bam >haps.vcf
#
# R64-1-1
#
# snpEFF download R64-1-1.86
#
# id=dp in vcf
#
# Produce the following set of plots for your variants:
#
# the read depth distribution across each called variant
# the genotype quality distribution
# the allele frequency spectrum of your identified variants
# a summary of the predicted effect of each variant as determined by snpEff (barplot?)
#

# for plot 1 x is each line and y is DP
#
# for plot 2 x is ordered least to greatest variants by y of GQ
#
# for plot 3 x is variant y is AF
#
# for plot 4 x is ?

import sys
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#import the call.vcf output file from snpEff
#separate by semicolon

# df = pd.read_csv(sys.argv[1], sep = ";", header=None).iloc[:,5]
# df = pd.read_csv( sys.argv[1], sep=";")

calls = open(sys.argv[1])
data = []
depth = []
allele_frequency = []
depthvallist = []
qual_list = []

for line in calls:
    
    if line.startswith("#"):
        # re.split(";",line)
        continue
    
    fields = line.split()
    qual = qual_list.append(float(fields[5]))

    info = fields[7]
    for id_val in info.split(";"):
        id, val = id_val.split('=')
        if id == "AF":
            allele_frequency.append(float(val.split(",")[0]))
        elif id == "DP":
            depth.append(val)
            
    #alternate is list comprehension
    # name = "Peter"
    # letters = [x for x in name]
    # --> ["P", "e", "t", "e", "r"]
    
    info1 = dict([x.split("=") for x in fields[7].split(";")])
 #    af = int(info.loc["AF"])
 #    dp = int(info["DP"])
    
 
	# af = int(info.loc["AF"])
    
    
    
    
    fordepth = line.rstrip("\r\n").split("\t")
    
    fields = data.append(line.rstrip("\r\n").split("\t"))
  
    depth = fordepth[7].split(";")
    depthisarg2 = int(depth[7].split("=")[1].split(',')[0])
    depthvalues = depthvallist.append(depthisarg2)

    continue
    
df = pd.DataFrame(data)
# df.assign(depthisarg2)
sorteddepth = depthvallist.sort()

fig, ax = plt.subplots()
plt.hist(depthvallist, bins = 1000)
ax.set_title("Depth Frequency Distribution for each Variant")
plt.ylabel("Frequency")
plt.xlim(0, 700)
plt.xlabel("Distribution")
fig.savefig("depth_hist.png")
plt.close()

fig, ax = plt.subplots()
plt.hist(allele_frequency, bins = 10)
ax.set_title("Allele Frequency Distribution for each Called Variant")
plt.ylabel("Frequency")
plt.xlabel("Distribution")
fig.savefig("allele_frequency_hist.png")
plt.close()

fig, ax = plt.subplots()
plt.hist(qual_list, bins = 10000)
plt.xlim(0, 5000)
ax.set_title("Quality Distribution for each Called Variant")
plt.ylabel("Frequency")
plt.xlabel("Distribution")
fig.savefig("quality_hist.png")
plt.close()

effects = pd.read_table(sys.argv[2], header=0)
effects.columns = ['GeneName', 'GeneId', 'TranscriptId', 'BioType', 'HIGH impact', 'LOW impact', 'MODERATE impact', 'MODIFIER impact', 'Conservative Inframe Deletion', 'Conservative Inframe Insertion', 'Disruptive Inframe Deletion', 'Disruptive Inframe Insertion', 'Downstream Gene Variant', 'Frameshift Variant', 'Initiator Codon Variant', 'Intron Variant', 'Missense Variant', 'Non-coding Transcript Exon Variant', 'Splice Acceptor Variant', 'Splice Donor Variant', 'Splice Region Variant', 'Start lost', 'Stop gained', 'Stop lost', 'Stop Retained Variant', 'Synonymous Variant', 'Upstream Gene Variant']

effects_list = effects.columns.tolist()

totals = 0

effect_dict = {}

for i in range(len(effects_list) - 4 ):
	f = effects.loc[:, effects_list[i + 4]]
	total = sum(f)
	effect_dict[effects_list[i+4]] = total
	totals += total
	continue

effects_list1 = []
effects_list1values = []

for key, value in effect_dict.items():
	print(key + ": \t", value)
	effects_list1.append(key)
	effects_list1values.append(value)
	continue

fig, ax = plt.subplots()
plt.bar(effects_list1[1:4], effects_list1values[1:4])
ax.set_title("Variant Impact Effects Barplot")
plt.ylabel("Frequency")
# ax.set_xticks(np.arange(0.5, len(effects_list1) +0.5),)
ax.set_xticklabels(effects_list1, rotation=50)
plt.tight_layout()
fig.savefig("impact_effects_bar.png")
plt.close()


indel = effects_list1[5:9]
indelval = effects_list1values[5:9]

fig, ax = plt.subplots()
plt.bar(indel, indelval)
ax.set_title("Variant In/Del Barplot")
plt.ylabel("Frequency")
# ax.set_xticks(np.arange(0.5, len(effects_list1) +0.5),)
ax.set_xticklabels(effects_list1, rotation=50)
plt.tight_layout()
fig.savefig("InDel_effects_bar.png")
plt.close()

# totalinframe = effects.loc[:, "variants_effect_conservative_inframe_deletion"]

# Type (alphabetical order)	 	Count	Percent
# conservative_inframe_deletion	 	6	0.002%
# conservative_inframe_insertion	 	6	0.002%
# disruptive_inframe_deletion	 	18	0.006%
# disruptive_inframe_insertion	 	9	0.003%
# downstream_gene_variant	 	121,839	42.086%
# frameshift_variant	 	153	0.053%
# initiator_codon_variant	 	1	0%
# intergenic_region	 	15,562	5.375%
# intron_variant	 	392	0.135%
# missense_variant	 	9,358	3.232%
# non_coding_transcript_exon_variant	 	224	0.077%
# splice_acceptor_variant	 	1	0%
# splice_donor_variant	 	1	0%
# splice_region_variant	 	64	0.022%
# start_lost	 	21	0.007%
# stop_gained	 	72	0.025%
# stop_lost	 	32	0.011%
# stop_retained_variant	 	31	0.011%
# synonymous_variant	 	15,676	5.415%
# upstream_gene_variant)


   

# with open(sys.argv[1]) as calls:
#   for line in calls:
#     if line != '#': # strings can be accessed like lists - they're immutable sequences.
#         calls.to_csv()
#         continue





























