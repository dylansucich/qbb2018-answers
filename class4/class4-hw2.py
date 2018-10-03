#!/usr/bin/env python3
# ./class4-hw2.py BYxRM_segs_saccer3.bam.simplified.vcf

import sys
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Plot the allele frequency
# input is the simplified vcf
AF = []
vcf = open(sys.argv[1])
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
