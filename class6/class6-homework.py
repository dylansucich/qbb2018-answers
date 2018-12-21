#!/usr/bin/env python3

# plot number of peaks gained and number of peaks lost

# Go through start ends in narrow peaks files and find where they align in
# mouse genome and the feature type associated with it

import sys
import pandas as pd
import matplotlib.pyplot as plt

#Usage: class6-homework.py peaksGained peaksLost Mus_musculus.GRCm38.94_features.bed
GainP = 0
GainE = 0
GainI = 0
LostP = 0
LostE = 0
LostI = 0
peaksGainedExon = 0 
peaksGainedpromotor = 0 
peaksGainedintron = 0 
peaksLostExon = 0 
peaksLostpromotor = 0 
peaksLostintron = 0 
peaksGainedFeatures = []
peaksLostFeatures = []
featStartEnd = []
FeaturesType = {}
GainedStartEnd = {}
FeatureStartEnd = {}
LostStartEnd = {}
GainedStarts = []
LostStarts = []

peaksGained = pd.read_table(sys.argv[1])
peaksGained.columns = ['chr', 'start', 'end', "PeakNum", "idk1", "idk2", "idk3", "idk4", "idk5", "idk6"]
peaksLost = pd.read_table(sys.argv[2])
peaksLost.columns = ['chr', 'start', 'end', "PeakNum", "idk1", "idk2", "idk3", "idk4", "idk5", "idk6"]

Features = pd.read_table(sys.argv[3], header=None)
Features.columns = ['chr', 'start', 'end', 'feature', 'idk1', 'Direction']

# for start, end, feature in Features:


for i in Features.index:
	start, end, feature = Features.loc[i, ['start', 'end', 'feature']]
	for position in range(start,end):
		FeatureStartEnd[position] = feature

gainedTotal = []

for i in peaksGained.index:
	Gain_feature = []
	GainStart, GainEnd = peaksGained.loc[i, ["start", "end"]]
	for pos in range(GainStart, GainEnd):
		if pos in FeatureStartEnd:
			GainFeature = FeatureStartEnd[pos]
			if GainFeature not in Gain_feature:
				Gain_feature.append(GainFeature)

	gainedTotal.append(Gain_feature)

for i in gainedTotal:
	if len(i) == 0:
		continue
	for j in i:
		if j == "exon":
			GainE += 1
		if j == "intron":
			GainI += 1
		if j == "promoter":
			GainP += 1

lostTotal = []

for i in peaksLost.index:
	Lost_feature = []
	LostStart, LostEnd = peaksLost.loc[i, ["start", "end"]]
	for pos in range(LostStart, LostEnd):
		if pos in FeatureStartEnd:
			LostFeature = FeatureStartEnd[pos]
			if LostFeature not in Lost_feature:
				Lost_feature.append(LostFeature)
	lostTotal.append(Lost_feature)

for i in lostTotal:
	if len(i) == 0:
		continue
	for j in i:
		if j == "exon":
			LostE += 1
		if j == "intron":
			LostI += 1
		if j == "promoter":
			LostP += 1

print("Gained Promoter", GainP)
print("Gained Intron", GainI)
print("Gained Exon", GainE)

print("")
print("Lost Promoter", LostP)
print("Lost Intron", LostI)
print("Lost Exon", LostE)


Gained = peaksGained.shape[0]
Lost = peaksLost.shape[0]
xvalues = ['Gained', 'Lost']
yvalues = [Gained, Lost]

x2values = ['Promoter(+)', 'Exon(+)', 'Intron(+)', 'Promoter(-)', 'Exon(-)', 'Intron(-)']
y2values = [GainP, GainE, GainI, LostP, LostE, LostI]

fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (12, 4))

ax1.bar(xvalues, yvalues)
ax1.set_xticks(xvalues)
ax1.set_title("CTCF Binding Sites at G1E to ER4 Transition")
ax1.set_ylabel("Counts")
ax1.set_xlabel("CTCF binding sites")

ax2.bar(x2values, y2values)
ax2.set_xticks(x2values)
ax2.set_title("CTCF Binding Site Location at G1E to ER4 Transition")
ax2.set_ylabel("Counts")
ax2.set_xlabel("CTCF binding site locations")

fig.savefig("CTCF_ChIPseq.png")
plt.close(fig)













