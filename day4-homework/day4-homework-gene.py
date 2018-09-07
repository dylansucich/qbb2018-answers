#!/usr/bin/env python3

"""
Usage: day4-homework_one_gene.py <gene_name> <samples.csv> <ctab_directory>
Plot mean value for all transcripts of specified gene
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[2])


for index, sample, sex, stage in df.itertuples():
    filename = os.path.join( sys.argv[3], sample, "t_data.ctab" )
    ctab_df = pd.read_table( filename, index_col="t_name" )
    gene_name = sys.argv[1]
    roi = ctab_df.loc[:, "gene_name"] == gene_name
    compiled = ctab_df.loc[roi, "FPKM"]
    compiled_df = pd.DataFrame(compiled)
    

average = compiled_df.mean(axis = 1)


fig, ax = plt.subplots()
ax.scatter( list(compiled_df.index), list(average) )
ax.set_title(gene_name)
plt.xticks(rotation=90)
plt.ylabel("mRNA abundance (FPKM)")
plt.tight_layout()
fig.savefig( "MeanValueofTransciptforOneGene.png" )
plt.close( fig )        
        
