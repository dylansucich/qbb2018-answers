#!/usr/bin/env python3

"""
Usage: day4-homework_one_gene.py <samples.csv> <ctab_directory> <gene_name>
Plot mean value for all transcripts of multiple genes 
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1])

for i in range(2,len(sys.argv)):
    compiled = {}
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join( sys.argv[2], sample, "t_data.ctab" )
        ctab_df = pd.read_table( filename, index_col="t_name" )
        all_genes = ctab_df.gene_name.unique()

    if sys.argv[i] in all_genes:
        gene_name = sys.argv[i]
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
        fig.savefig( "".join(["day4-homework-multiple-genes_", sys.argv[i], ".png" ]))
        plt.close( fig )        
    else:
        print("Provided gene not in data: ", sys.argv[i])
        continue