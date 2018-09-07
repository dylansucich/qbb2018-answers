#!/usr/bin/env python3

"""
Usage: ./day5-exercise6.py <sample.ctab> <n_histonemodificationfiles.tab>

Transform your response variable (gene expression) into units of log(FPKM + 1)
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as sm
import scipy.stats
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

log_combined_avg_FPKMs = np.log2(combined_avg.loc[:, "FPKMs"] + 1)
combined_avg_with_LogFPKMs = combined_avg.assign(log2FPKMs = log_combined_avg_FPKMs)

model = sm.ols(formula = "log2FPKMs ~ H3K27ac_mean + H3K27me3_mean + H3K4me1_mean + H3K4me3_mean + H3K9ac_mean", 
    data = combined_avg_with_LogFPKMs)
results = model.fit()
print(results.summary())

mu = np.mean(results.resid)
sigma = np.std(results.resid)

fig, ax = plt.subplots()
n, bins, patches = ax.hist(x = results.resid)
y = scipy.stats.norm.pdf(bins, mu, sigma)
ax.plot(bins, y)
ax.set_ylabel = "Frequency"
ax.set_xlabel = "Value"
ax.set_title = "Histogram of Residuals"
fig.tight_layout()
fig.savefig("day5-exercise6_log2residuals.png")
plt.close(fig)