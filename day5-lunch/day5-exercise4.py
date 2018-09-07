#!/usr/bin/env python3

"""
Usage: ./day5-exercise4.py <sample.ctab> <n_histonemodificationfiles.tab>

Perform ordinary linear regression for all of the Histone markersfour marks to determine how
 predictive each is of gene expression.

   bigWigAverageOverBed in.bw in.bed out.tab
The output columns are:
   name - name field from bed, which should be unique
   size - size of bed (sum of exon sizes
   covered - # bases within exons covered by bigWig
   sum - sum of values over all bases covered
   mean0 - average over bases with non-covered bases counting as zeroes
   mean - average over just covered bases

"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as sm
# from statsmodels.sandbox.regression.predstd import wls_prediction_std

df = pd.read_csv( sys.argv[1], sep="\t")
fpkms = df.loc[:, "FPKM"]

H3K27ac   = pd.read_csv( sys.argv[2], sep="\t", header = None).iloc[:,-1]
H3K27me3  = pd.read_csv( sys.argv[3], sep="\t", header = None).iloc[:,-1]
H3K4me1   = pd.read_csv( sys.argv[4], sep="\t", header = None).iloc[:,-1]
H3K4me3   = pd.read_csv( sys.argv[5], sep="\t", header = None).iloc[:,-1]
H3K9ac    = pd.read_csv( sys.argv[6], sep="\t", header = None).iloc[:,-1]

combined_avg = pd.concat([H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9ac], axis=1)
combined_avg = combined_avg.assign(fpkms = fpkms)
combined_avg.columns = ["H3K27ac_mean", "H3K27me3_mean", "H3K4me1_mean", "H3K4me3_mean", "H3K9ac_mean", "FPKMs"]

model = sm.ols(formula = "FPKMs ~ H3K27ac_mean + H3K27me3_mean + H3K4me1_mean + H3K4me3_mean + H3K9ac_mean", 
    data = combined_avg)
results = model.fit()
print(results.summary())



# come back to this
# for i in range(2, len(sys.argv)):
#     filename = os.path.splitext(sys.argv[i])[0]
#     filename = filename.rsplit("_")[1]
#     filename_df = pd.read_csv( sys.argv[i], sep="\t")
#     filename_mean = filename_df.iloc[:,5]
#     df.assign( {} = filename_mean, axis = 1),




    
    
