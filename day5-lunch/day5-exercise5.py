#!/usr/bin/env python3

"""
Usage: ./day5-exercise5.py <sample.ctab> <n_histonemodificationfiles.tab>

Plot the residuals of your model and evaluate to evaluate this assumption with regard to your data.

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

model = sm.ols(formula = "FPKMs ~ H3K27ac_mean + H3K27me3_mean + H3K4me1_mean + H3K4me3_mean + H3K9ac_mean", 
    data = combined_avg)
results = model.fit()
results.summary()

mu = np.mean(results.resid)
sigma = np.std(results.resid)

fig, ax = plt.subplots()
n, bins, patches = ax.hist(x = results.resid, bins = 3000)
y = scipy.stats.norm.pdf(bins, mu, sigma)
ax.plot(bins, y)
ax.set_ylabel = "Frequency"
ax.set_xlabel = "Value"
ax.set_title = "Histogram of Residuals"
ax.set_xlim(-30,100)
fig.tight_layout()
fig.savefig("day5-exercise5_residuals.png")
plt.close(fig)