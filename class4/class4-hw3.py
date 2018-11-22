#!/usr/bin/env python3
# ./class4-hw2.py BYxRM_PhenoData.txt plink.TREATMENT.qassoc
# ./class4-hw2.py BYxRM_PhenoData.txt *.qassoc

import sys
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Manhattan hint
# Q*
#treatment in arg after plink.Cadmium_Chloride.qassoc Cadmium_Chloride  Caffeine	Calcium_Chloride	Cisplatin	Cobalt_Chloride	Congo_red	Copper	Cycloheximide	Diamide	E6_Berbamine	Ethanol	Formamide	Galactose	Hydrogen_Peroxide	Hydroquinone	Hydroxyurea	Indoleacetic_Acid	Lactate	Lactose	Lithium_Chloride	Magnesium_Chloride	Magnesium_Sulfate	Maltose	Mannose	Menadione	Neomycin	Paraquat	Raffinose	SDS	Sorbitol	Trehalose	Tunicamycin	x4-Hydroxybenzaldehyde	x4NQO	x5-Fluorocytosine	x5-Fluorouracil	x6-Azauracil	Xylose	YNB	YNB:ph3	YNB:ph8	YPD	YPD:15C	YPD:37C	YPD:4C	Zeocin

tx = [] #list of just the names of treatments
txfiles = [] #list of plink.TX.qassoc
numOfTx = 0

txpheno = pd.read_table(sys.argv[1], index_col = None)

# name = pd.DataFrame(sys.argv[1]) # pull treatment
for name in sys.argv[2:]:   
    txfiles.append(name)
    txment = name.split(".")[-2]
    tx.append(txment)

for i in txfiles:
	txtable = pd.read_table(i, delim_whitespace = True)
	name = i.split("/plink.")[1].split(".")[0]
	print(name)
	print(txtable)
	BasePairPosition = txtable.loc[:, "BP"].tolist()
	Pval = -np.log(txtable.loc[:,"P"])
	df2 = txtable.assign(Position= BasePairPosition, Pvalue= Pval)
	print(txtable)
	df2.loc[df2['Pvalue'] > 5, 'Significant Pval'] = df2['Pvalue']
	print(df2)
	fig, ax = plt.subplots(figsize=(18,6))
	ax.scatter(x=df2.loc[:,"Position"], y=df2.loc[:,"Significant Pval"], color="black", marker=".", s= 0.3, label='Significant Values')
	ax.scatter(x=txtable.loc[:,"BP"], y=df2.loc[:,"Pvalue"], color="red", s= 0.1)
	ax.axhline(y= 5, linewidth=0.5, color='r', marker= "_")
	ax.set_title("Manhattan Plot for " + name)
	ax.set_xlabel("SNPs across S.cerevisiae Genome Position")
	ax.set_ylabel("-log(p-value)")
	fig.savefig(name + "_manhattanplot_" + ".png")
	plt.close()
	continue
    # chrP = txdf.loc[:,-1]
    # print(chrP)
   





# with open(fname) as f:
#     for line in f:
#     #do stuff
#         pass
#
#     fig.savefig("manhattan_{}.png".format(treatment))