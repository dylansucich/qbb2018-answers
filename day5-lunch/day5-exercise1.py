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
newdf.to_csv("SRR072893.bed", header= False, index=False, sep="\t")