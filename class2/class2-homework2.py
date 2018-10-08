#!/usr/bin/env python3

"""
Sort potplot file from lastz first with
sort -k 1 -n velvet_low_dotplot.out > velvet_low_dotplot_sorted.out

Usage: ./class2-homework2.py velvet_low_dotplot_sorted.out velvet_low_dotplot_sorted
"""

import sys
import os
import pandas as pd
import fasta
import itertools as it
from statsmodels.stats import weightstats as stests
import matplotlib.pyplot as plt

count = 0
plt.figure()

for value in open(sys.argv[1]):
    if "zstart1" in value:
        continue
    else:
        fields = value.split("\t")
        plt.plot([int(fields[0]), int(fields[2])], [count, count + int(fields[1])])
        
        count += int(fields[1])

plt.xlim( 0, 100000 )
plt.ylim( 0, 100000 )
plt.xlabel("Reference Position")
plt.ylabel("Contig Position")
plt.title(sys.argv[2])
plt.tight_layout()
plt.savefig(str(sys.argv[2]) + ".png")
plt.close()