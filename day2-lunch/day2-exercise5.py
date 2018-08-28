#!/usr/bin/env python3

import sys 
import pandas as pd

if len(sys.argv) > 1:
    f = open( sys.argv[1] )
else:
    f = sys.stdin


totalreads = 0
mapq = 0

for i, line in enumerate( f ):
    if line[0] =="@":
        continue
    totalreads += 1
    fields = line.strip().split("\t")
        
    if fields[2] != "*":
        mapq += int(fields[4])
        
    averageMapQ = float(mapq / totalreads)

print(averageMapQ)
    
