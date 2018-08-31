#!/usr/bin/env python3

"""
Usage: ./day5-exercise1.py <t_data.ctab> 

For SRR072893.t_tab in ~/data/results/stringtie/ approximate promoter 
region and output a .bed file with columns named chromosome, start, end, t_name.

no negative
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

df = pd.read_csv( sys.argv[1], sep="\t")
# promoter = {}
# promoter_location = []
#
# for line in open( sys.argv[1] ):
#     fields = line.strip.split()
    
    
# print(df.loc[:,"strand"])
#
#chromosome column, strand start end 

# df["Start > 500"] = np.where(df["start"] > 500, "yes", "no")

pos_roi = df.loc[df["strand"] == "+"]
posover500 = pos_roi.loc[df["start"] > 500]

pos_promoter_start = df.loc[:,"start"] - 500
pos_promoter_end = df.loc[:,"start"] + 500
tsspos_df = posover500.assign( Promoter_Start = pos_promoter_start, Promoter_End = pos_promoter_end )


neg_roi = df.loc[df["strand"] == "-"]
negover500 = neg_roi.loc[df["start"] > 500]

neg_promoter_start = df.loc[:,"start"] + 500
neg_promoter_end = df.loc[:,"start"] - 500
tssneg_df = negover500.assign( Promoter_Start = neg_promoter_start, Promoter_End = neg_promoter_end )
# print(neg_promoter_start)

df2 = pd.concat([tsspos_df, tssneg_df])


newdf = df2.loc[:,["chr", "Promoter_Start", "Promoter_End", "t_name"]] 
newdf.to_csv(sys.stdout, header= False, index=False, sep="\t")


#
#


# print(pos_df)
# promoter_start = df
# positive_strand = df[pos_roi,
#
#
#     promoter_start = over500.loc[:,"start"] - 500
#     promoter_end = over.loc[:,"start"] + 500
#     tss_df = df.assign( PromoterStart = promoter_start, Promoter_end = promoter_end )
#
# promoter_start = df.loc[:,"end"] - 500
# promoter_end = df.loc[:,"end"] + 500
# tss_df = df.assign( PromoterStart = promoter_start, Promoter_end = promoter_end )
#
#
#
#
#
# print(tss_df)
#
#

# test = pd.concat([df, is_strand_positive])
# print(test)

# if df.loc[:,"strand"] == "+" in df.loc[:,:]:
#     promoter_start = df.loc[:,"start"] - 500
#     print(promoter_start)
#     promoter_end = df.loc[:,"start"] + 500
#     print(promoter_start)
# # print()
# print(neg_strand)
#
#     #column start +500
#     #column end -500
#     df.loc[] transcript_start