#!/usr/bin/env python3

#/data/results/stringtie/SRR072893/t_data.ctab
"""
Usage: ./day4-exercise1.py <threshold> <ctab_file1> <ctab_file2> ... <ctab_filen>

create csv file with FPKMs from variable numbers of samples 
- assumes ctab_file in directory with same name
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

th = int( sys.argv[1] )
arg = len(sys.argv)
#print(arg)
samples = {}
sums = {}

for i in range(2, arg):
    #roi_fpkm = df.loc[:,"FPKM"] > float( sys.argv[2] )
    name = sys.argv[i].split( os.sep )[-2]
    fpkm = pd.read_csv( sys.argv[i], sep="\t", index_col="t_name" ).loc[:,"FPKM"]
    samples[name] = fpkm
    #
    # name2= sys.argv[2].split( os.sep )[-2]
    # fpkm2 =pd.read_csv( sys.argv[2], sep="\t", index_col="t_name" ).loc[:,"FPKM"]


fpkms_df = pd.DataFrame(samples)
sum_fpkms = fpkms_df.sum(axis = 1)

#meets_th = sum_fpkms.loc[]

meets_th = sum_fpkms > th

results = sum_fpkms.index[meets_th == True]
print(sum_fpkms, results)
#results = sum_fpkms.index[boolean=meets_th]

#print(results)
#df2 = sum_fpkms.assign( boolean=meets_th )
#df2.loc[:,:].to_csv( sys.stdout, sep="\t" )
# print(df2)
# print(sum_fpkms)
# print(meets_th)
# print(results)
#
# meets_th = sum_fpkms > th
#
# sum_fpkms.loc[meets_th,:].to_csv( sys.stdout, sep="\t", index=False)




#fpkms_df.to_csv( sys.stdout )

# fig, ax = plt.subplots()
# ax[0].hist( np.log(gene_lengths) )
# ax[1].hist( (gene_exons) )
# ax[2].scatter(np.log(gene_lengths), gene_exons)
# fig.savefig("two_plots.png")
# plt.close(fig)


#
#
# sum_fpkms = fpkms_df.sum(axis = 1)
# print(sum_fpkms)
#
# th = int( sys.argv[1] )
#
# print(df.loc[df[:,1] > th)





# for i in range(2, arg):
#     #roi_fpkm = df.loc[:,"FPKM"] > float( sys.argv[2] )
#     name = sys.argv[i].split( os.sep )[-2]
#     total = pd.DataFrame.sum(fpkms_df.loc[:,:])
#     sums[total] = name
#
# print(sums)





# fig, ax = plt.subplots()
# ax[0].hist( np.log(gene_lengths) )
# ax[1].hist( (gene_exons) )
# ax[2].scatter(np.log(gene_lengths), gene_exons)
# fig.savefig("two_plots.png")
# plt.close(fig)


# fpkms = { name1 : fpkm1, name2 : fpkm2 }
#
# fpkms_df = pd.DataFrame(fpkms)
#
# fpkms_df.to_csv( sys.stdout )
#
    
#df.loc[roi,:].to_csv( sys.stdout, sep="\t", index=False)
    

# df = pd.read_csv( sys.argv[1], sep="\t", index_col="t_name" )
#
# gene_lengths =df.loc[:,"length"]
# print( "Mean gene lengths: {0}".format( np.mean( gene_lengths )))
#
# gene_exons = df.loc[:,"num_exons"]
#
# lengths_vs_exons = np.corrcoef( gene_lengths, gene_exons)[0,1]
#
# #print( type( lengths_vs_exons ))
# print( "Pearson's f: {0:0.4f}".format( lengths_vs_exons) )




#plotting
#required
# fig, ax = plt.subplots()
# ax.scatter( gene_lengths, gene_exons)
# ax.set_title("gene_lengths vs gene_exons")
# ax.set_x
# #Your code goes here
# #required
# fig.savefig("scatter.png")
# #required
# plt.close(fig)

#ax as (ax1, ax2) or ax then ax[#]
#the number in subplots is the actual number of the argument 
