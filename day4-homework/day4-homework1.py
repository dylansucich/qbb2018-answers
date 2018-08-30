#!/usr/bin/env python3

"""
Usage: ./day4-homework.py <samples.csv> <t_data.ctab directories> 

merge FPKM values from all ctab files 
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv( sys.argv[1] )
# key = sample stage
# value = FPKM
compiled= {}
# for item in df.itertuples():
#     print( item )
for index, sample, sex, stage in df.itertuples():
    filename = os.path.join( sys.argv[2], sample, "t_data.ctab" )
    
    ctab_df = pd.read_table( filename, index_col="t_name" )
    
    compiled[ sex + "_" + stage ] = ctab_df.loc[:, "FPKM"]

compiled_df = pd.DataFrame( compiled )
compiled_df.to_csv( sys.stdout )

#
# fig, ax = plt.subplots()
# ax.plot( fpkms )
# fig.savefig( "timecourse.png" )
# plt.close( fig )