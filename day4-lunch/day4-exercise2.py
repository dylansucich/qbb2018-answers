#!/usr/bin/env python3

#/data/results/stringtie/SRR072893/t_data.ctab
"""
Usage: ./day4-exercise2.py <sample1/t_data.ctab> <sample2/t_data.ctab>

plot FPKM1 on x-axis and FPKM2 on y-axis
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df1 = pd.read_csv( sys.argv[1], sep="\t", index_col="t_name" )
df2 = pd.read_csv( sys.argv[2], sep="\t", index_col="t_name" )

fpkm1 =df1.loc[:,"FPKM"]
fpkm2 =df2.loc[:,"FPKM"]

fit = np.polyfit(fpkm1, fpkm2, 1)
p = np.poly1d(fit)
xplin = np.linspace(0,fpkm1.max())

#for log axis
fig, ax = plt.subplots()
ax.scatter(fpkm1, fpkm2, alpha=0.1)
plt.plot(xplin, p(xplin), '-', color='m')
plt.yscale("log")
plt.xscale("log")
plt.axis([.001, 10000,.001, 10000])
axes = plt.gca()
ax.set_title("FPKM2 vs. FPKM1 Scatter Plot")
plt.ylabel("log(FPKM2)")
plt.xlabel("log(FPKM1)")
plt.savefig("fpkmlog.png")
plt.close()

# #for non-log axis
# fig, ax = plt.subplots()
# ax.scatter(fpkm1, fpkm2, alpha=0.1)
# plt.plot(xplog, p(xplog), '-', color='m')
# # plt.yscale("log")
# # plt.xscale("log")
# ax.set_title("FPKM2 vs. FPKM1 Scatter Plot")
# plt.ylabel("FPKM2")
# plt.xlabel("FPKM1")
# plt.savefig("fpkm.png")
# plt.close()
# ./day4-exercise2.py ~/data/results/stringtie/SRR072893/t_data.ctab ~/data/results/stringtie/SRR072915/t_data.ctab