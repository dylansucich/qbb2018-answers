#!/usr/bin/env python3

"""
Usage: ./day4-homework.py <t_name> <samples.csv> <t_data.ctab directories>

merge FPKM values from all ctab files 
plot similarly to Lott et al., 2011 PLoS Biology 

Transcript ID: FBtr0331261 for Sxl gene
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

transcript_id = sys.argv[1]

def sex_abundance(gender):
    df = pd.read_csv( sys.argv[2] )
    # df2 = pd.read_csv( sys.argv[4] )
 #    df = pd.concat([df1, df2])
    
    soi = df.loc[:,"sex"] == gender
    df = df.loc[soi,:]
    compiled= {}
    fpkms = []
    for index, sample, sex, stage in df.itertuples():
        name = stage
        filename = os.path.join( sys.argv[3], sample, "t_data.ctab" )
        ctab_df = pd.read_table( filename, index_col="t_name" )
        transcript_id = sys.argv[1]
        roi = ctab_df.index == transcript_id
        compiled[name] = ctab_df.loc[roi, "FPKM"]
           
    compiled_df = pd.DataFrame( compiled )
    return(compiled_df)

females = sex_abundance("female")
males = sex_abundance("male")

fig, ax = plt.subplots()
ax.plot( list(females), females.loc[transcript_id,:], color="red")
ax.plot( list(males), males.loc[transcript_id,:], color="blue")
ax.set_title("Sxl")
plt.xticks(rotation=90)
plt.ylabel("mRNA abundance (FPKM)")
plt.xlabel("developmental stage")
plt.legend(bbox_to_anchor=(1.05,0.5), loc= 2, borderaxespad= 0., labels= ['Females', 'Males'])
plt.tight_layout()
fig.savefig( "FBtr0331261_Male_female.png" )
plt.close( fig )