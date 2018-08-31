#!/usr/bin/env python3

"""
Usage: ./day4-homework.py  samples.csv <replicats.csv> <stringtie> 

merge FPKM values from all ctab files 
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def sex_abundance(gender):
    df1 = pd.read_csv( sys.argv[1] )
    df2 = pd.read_csv( sys.argv[2] )
    df = pd.concat([df1, df2])
    
    soi = df.loc[:,"sex"] == gender
    df = df.loc[soi,:]
    compiled= {}
    fpkms = []
    for index, sample, sex, stage in df.itertuples():
        name = stage
        filename = os.path.join( sys.argv[3], sample, "t_data.ctab" )
        ctab_df = pd.read_table( filename, index_col="t_name" )
        compiled[name] = ctab_df.loc[:, "FPKM"]
           
    compiled_df = pd.DataFrame( compiled )
    print(sample)
    return(compiled_df)

females = sex_abundance("female")
subset_females = females.iloc[:,2]
subset_females = subset_females + 1 

males = sex_abundance("male")
subset_males = males.iloc[:,2]
subset_males = subset_males + 1 

y_axis = np.log2(subset_females/subset_males)
x_axis = np.log10(np.sqrt(subset_males * subset_females))


fig, ax = plt.subplots()
ax.scatter(x_axis,y_axis, alpha = 0.1, color="green")
plt.ylabel("Log2 Ratio of Female to Male transcript abundance (FKPM)")
plt.xlabel("Intensity: Log10 Ratio of Female * Male transcript abundance (FKPM)")
ax.set_title("MA-Plot")
fig.savefig("ma_plot.png")
plt.close( fig )

# ax.plot( list(females), females.loc[transcript_id,:], color="red")
# ax.plot( list(males), males.loc[transcript_id,:], color="blue")
# ax.set_title("Sxl")
# plt.xticks(rotation=90)
# plt.ylabel("mRNA abundance (FPKM)")
# plt.xlabel("developmental stage")
# plt.legend(bbox_to_anchor=(1.05,0.5), loc= 2, borderaxespad= 0., labels= ['Females', 'Males'])
# plt.tight_layout()
# fig.savefig( "malefemale.png" )
# plt.close( fig )