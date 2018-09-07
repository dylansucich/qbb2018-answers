#!/usr/bin/env python3

"""
Usage: ./day4-homework.py <t_name> <samples.csv> <t_data.ctab directories> <replicates.csv>

merge FPKM values from all ctab files and replicates  
plot similarly to Lott et al., 2011 PLoS Biology 

"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

transcript_id = sys.argv[1]

def sex_abundance(gender):
    df1 = pd.read_csv( sys.argv[2] )
    df2 = pd.read_csv( sys.argv[4] )
    df = pd.concat([df1, df2])
    
    soi = df.loc[:,"sex"] == gender
    df = df.loc[soi,:]
    compiled= {}

    for index, sample, sex, stage in df.itertuples():
        name = (sample + "_" + stage)
        filename = os.path.join( sys.argv[3], sample, "t_data.ctab" )
        ctab_df = pd.read_table( filename, index_col="t_name" )
        transcript_id = sys.argv[1]
        roi = ctab_df.index == transcript_id
        compiled[name] = ctab_df.loc[roi, "FPKM"]
           
    compiled_df = pd.DataFrame( compiled )
    return(compiled_df)

females = sex_abundance("female")
split_fname = females.columns.str.split("_")
fnames = []
for i in range(0, len(split_fname)):
    new_name = split_fname[i][1]
    fnames.append(new_name)
    
females.columns = fnames

males = sex_abundance("male")
split_mname = males.columns.str.split("_")
mnames = []
for i in range(0, len(split_mname)):
    new_name = split_fname[i][1]
    mnames.append(new_name)

males.columns = mnames

fig, ax = plt.subplots()
ax.plot(list(females)[0:8], females.loc["FBtr0331261", :][0:8], color = "red")
ax.plot(list(males)[0:8], males.loc["FBtr0331261", :][0:8], color = "blue")

ax.plot(list(females)[8:], females.loc["FBtr0331261", :][8:], color = "red", linestyle = "dashed", alpha = 0.5)
ax.plot(list(males)[8:], males.loc["FBtr0331261", :][8:], color = "blue", linestyle = "dashed", alpha = 0.5)

ax.set_title("Sxl")
plt.xticks(rotation = 90)
plt.ylabel("mRNA abundance (FPKM)")
plt.xlabel("developmental stage")
plt.legend(bbox_to_anchor=(1.05,0.5), loc= 2, borderaxespad= 0., labels= ['Females', 'Males', "Female replicates", "Male replicates"])
plt.tight_layout()
fig.savefig("FBtr0331261_with_replicates.png")
plt.close(fig)
