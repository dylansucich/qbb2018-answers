#!/usr/bin/env python3
#input is ./class4-hw.py plink.eigenvec \
# BYxRM_segs_saccer3.bam.simplified.vcf \
# BYxRM_PhenoData.txt \ Cadmium_Chloride  Caffeine  Calcium_Chloride    Cisplatin	Cobalt_Chloride	Congo_red	Copper	Cycloheximide	
# Diamide    E6_Berbamine    Ethanol    Formamide    Galactose    Hydrogen_Peroxide    Hydroquinone    Hydroxyurea    Indoleacetic_Acid    Lactate    Lactose    Lithium_Chloride    Magnesium_Chloride    Magnesium_Sulfate    Maltose    Mannose    Menadione    Neomycin    Paraquat    Raffinose    SDS    Sorbitol    Trehalose    Tunicamycin    x4-Hydroxybenzaldehyde    x4NQO    x5-Fluorocytosine    x5-Fluorouracil

import sys
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#Show that all the strains are genetically related
#input is the plink.eigenvec
df = pd.read_csv(sys.argv[1], sep = " ", header= None)

fig, ax = plt.subplots()

ax.scatter(df.iloc[:, 2], df.iloc[:, 3])
ax.set_title("PCA2 vs. PCA1: Genetic Relatedness between strains")
plt.xticks(rotation = 90)
plt.ylabel("PCA2")
plt.xlabel("PCA1")
plt.tight_layout()
fig.savefig("PCA2vsPCA1.png")
plt.close(fig)

# Plot the allele frequency
# input is the simplified vcf
AF = []
vcf = open(sys.argv[2])
for line in vcf:
    if line.startswith("#"):
        continue
    fields = line.strip("\r\n").split("\t")
    AF_Values = fields[7][3:]
    AF_Value = AF_Values.split(",")[0]#might be two
    AF.append(float(AF_Value))

fig, ax = plt.subplots()
ax.hist(AF, bins = 1000)
ax.set_title("Allele Frequency")
plt.ylabel("Frequency")
plt.xlabel("Allele Frequency")
plt.xticks(rotation = 90)
plt.tight_layout()
fig.savefig("Allele_Freq.png")
plt.close(fig)

# Manhattan hint
# Q*
#treatment in arg after plink.qassoc Cadmium_Chloride  Caffeine	Calcium_Chloride	Cisplatin	Cobalt_Chloride	Congo_red	Copper	Cycloheximide	Diamide	E6_Berbamine	Ethanol	Formamide	Galactose	Hydrogen_Peroxide	Hydroquinone	Hydroxyurea	Indoleacetic_Acid	Lactate	Lactose	Lithium_Chloride	Magnesium_Chloride	Magnesium_Sulfate	Maltose	Mannose	Menadione	Neomycin	Paraquat	Raffinose	SDS	Sorbitol	Trehalose	Tunicamycin	x4-Hydroxybenzaldehyde	x4NQO	x5-Fluorocytosine	x5-Fluorouracil	
tx = []
fname = pd.DataFrame(sys.argv[3]) # pull treatment
for name in sys.argv[3:]:
    tx.append(name)
    
print(tx)
    # with open(fname) as f:
    #     for line in f:
    #         #do stuff
    #         pass
    #
    # fig.savefig("manhattan_{}.png".format(treatment))





#
#     for line in calls:
#
#         if line.startswith("#"):
#             # re.split(";",line)
#             continue
#
#         fields = line.split()
#         info = fields[7]
#         for id_val in info.split(";"):
#             id, val = id_val.split('=')
#             if id == "AF":
#                 allele_freq.append(float(val.split(",")[0]))
#             elif id == "DP":
#                 depth.append(float(val))
#
#         #alternate is list comprehension
#         # name = "Peter"
#         # letters = [x for x in name]
#         # --> ["P", "e", "t", "e", "r"]
#
#         info1 = dict([x.split("=") \
#             for x in fields[7].split(";")])
#         af = float(info["AF"])
#         dp = float(info["DP"])
#
#
#
#
# fig, ax = plt.subplots()
# ax.scatter(df.iloc[:, 2], df.iloc[:, 3])
# ax.set_title("Allel Frequency")
# plt.xticks(rotation = 90)
# plt.ylabel("PCA2")
# plt.xlabel("PCA1")
# plt.tight_layout()
# fig.savefig("PCA2vsPCA1.png")
# plt.close(fig)










 

