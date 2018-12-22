#!/usr/bin/env python3

# sort ER4 peaks.narrowPeak by p-value and take the top 100 sequences
# Usage: class8-homework.py ./memechip_out/fimo_out_1/fimo.gff ./top100_ER4_peaks.narrowPeak

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from collections import Counter


top100_peaks = {}
motifs = []

for line in open(sys.argv[2]):

	fields = line.rstrip('\r\n').split('\t')
	peak_name = fields[3]
	peak_start = int(fields[1])
	peak_stop = int(fields[2])

	top100_peaks[peak_name] = [peak_start, peak_stop]

for line in open(sys.argv[1]):
	if line[0] == '#':
		continue
	fields2 = line.rstrip('\r\n').split('\t')

	motifStart = int(fields2[3])
	motifStop = int(fields2[4])

	peakStart = top100_peaks[peak_name][0]
	peakStop = top100_peaks[peak_name][1]
	
	motifVal = float((motifStart - peakStart) / (peakStop - peakStart))
	motifs.append(motifVal)


fig, ax = plt.subplots()
ax.hist(motifs, bins = 40, color='g')
ax.set_xlabel("Motif Position")
ax.set_ylabel("Frequency")
plt.tight_layout()
fig.savefig("CTCF_Motif.png")
plt.close(fig)






