#!/usr/bin/env python3
# input ./class4-hw1.py plink.eigenvec

import sys
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#Show that all the strains are genetically related
#input is the plink.eigenvec
df = pd.read_csv(sys.argv[1], sep = " ", header= None)

fig, ax = plt.subplots()

ax.scatter(df.iloc[:, 2], df.iloc[:, 3])
ax.set_title("PCA2 vs. PCA1: Genetic Relatedness between strains")
plt.xticks(rotation = 90)
plt.ylabel("PCA2")
plt.xlabel("PCA1")
plt.tight_layout()
fig.savefig("PCA2vsPCA1.png")
plt.close(fig)