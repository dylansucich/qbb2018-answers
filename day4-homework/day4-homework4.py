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

#
# fig, ax = plt.subplots()
# ax.plot( list(females), females.loc[transcript_id,:], color="red")
# ax.plot( list(males), males.loc[transcript_id,:], color="blue")
# ax.set_title("Sxl")
# plt.xticks(rotation=90)
# plt.ylabel("mRNA abundance (FPKM)")
# plt.xlabel("developmental stage")
# plt.legend(bbox_to_anchor=(1.05,0.5), loc= 2, borderaxespad= 0., labels= ['Females', 'Males'])
# plt.tight_layout()
# fig.savefig( "FBtr0331261_Male_female.png" )
# plt.close( fig )

#!/usr/bin/env python3

# Author: Margaret R. Starostik

# Run with: ./day4-homework_FBtr0331261_replicates.py
#
# """
# Usage: ./day4-homework_FBtr0331261_replicates.py
# Plots FBtr0331261 for male and female samples. Plot is styled similarly to Lott et al., 2011 PLoS Biology.
# Replicates have been added to plots as dotted lines.
# """
#
# import sys
# import os
# import pandas as pd
# import matplotlib.pyplot as plt
#
# def FPKMs_with_replicates(gender, csv_file, stringtie, transcriptID, replicates_file):
#     df1 = pd.read_csv(csv_file)
#     df2 = pd.read_csv(replicates_file)
#
#     df = pd.concat([df1, df2])
#     soi = df.loc[:, "sex"] == gender
#     df = df.loc[soi, :]
#
#     compiled_dict = {}
#     for index, sample, sex, stage in df.itertuples():
#         name = (sample + "_" + stage)
#
#         filename = os.path.join(stringtie, sample, "t_data.ctab")
#          #print(filename)
#
#         ctab_df = pd.read_table(filename, index_col = "t_name")
#
#          # Obtain information only for FPKM
#         transcript_id = transcriptID
#         roi = ctab_df.index == transcript_id
#
#         compiled_dict[name] = ctab_df.loc[roi, "FPKM"]
#
#     compiled_df = pd.DataFrame(compiled_dict)
#     #compiled_df = compiled_df.iloc[:, 4:]
#     return(compiled_df)
#
# females = FPKMs_with_replicates("female", "~/qbb2018/samples.csv", "~/data/results/stringtie", "FBtr0331261", "~/qbb2018/replicates.csv")
# #print(females)
# # Rename the columns so that x-axis of plot is formatted correctly
# split_fname = females.columns.str.split("_")
# fnames = []
# for i in range(0, len(split_fname)):
#     #print(split_fname[i])
#     new_name = split_fname[i][1]
#     fnames.append(new_name)
#     #print(new_name)
# #print(fnames)
# females.columns = fnames
# #print(list(females))
# #print(len(females))
# print(females)
#
# males = FPKMs_with_replicates("male", "~/qbb2018/samples.csv", "~/data/results/stringtie", "FBtr0331261", "~/qbb2018/replicates.csv")
# # Rename the columns so that x-axis of plot is formatted correctly
# split_mname = males.columns.str.split("_")
# mnames = []
# for i in range(0, len(split_mname)):
#     #print(split_mname[i])
#     new_name = split_fname[i][1]
#     mnames.append(new_name)
#     #print(new_name)
# #print(mnames)
# males.columns = mnames
# #print(list(males))
# #print(len(males))
# #print(males)
# #
# # Generate plot
# fig, ax = plt.subplots()
# ax.plot(list(females)[0:8], females.loc["FBtr0331261", :][0:8], color = "red", alpha = 0.6, marker = "o", markersize = 4)
# ax.plot(list(males)[0:8], males.loc["FBtr0331261", :][0:8], color = "blue", alpha = 0.6,  marker = "o", markersize = 4)
#
# # Replicates
# ax.plot(list(females)[8:], females.loc["FBtr0331261", :][8:], color = "purple", alpha = 0.6, linestyle = "dashed",  marker = "o", markersize = 4)
# ax.plot(list(males)[8:], males.loc["FBtr0331261", :][8:], color = "green", alpha = 0.6, linestyle = "dashed",  marker = "o", markersize = 4)
#
# # Adjustments to plots
# fig.suptitle("Sxl", fontsize = 12, style = "italic")
#
# plt.xlabel("developmental stage", fontsize = 12)
# plt.xticks(rotation = 90)
# ax.tick_params(bottom = False, left = False)
#
# plt.ylabel("mRNA abundance (FPKM)", fontsize = 12)
# plt.ylim((-10, 250))
#
# leg = plt.legend(bbox_to_anchor = (1.05, 0.6), loc = 2, borderaxespad = 0., labels = ["female", "male", "female replicate", "male replicate"])
# leg.get_frame().set_linewidth(0.0)
#
# fig.subplots_adjust(top = 0.4)
# plt.tight_layout()
# #plt.gcf().subplots_adjust(top=0.2)
# fig.savefig("day4-homework_FBtr0331261_with_replicates.png")
# plt.close(fig)