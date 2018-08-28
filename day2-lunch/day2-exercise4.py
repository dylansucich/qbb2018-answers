#!/usr/bin/env python3

import sys 

if len(sys.argv) > 1:
    f = open( sys.argv[1] )
else:
    f = sys.stdin

count = 0
for i, line in enumerate( f ):
    if line[0] =="@":
        continue
    fields = line.split("\t")
    if fields[2] != "*":
        count += 1
        if count < 11:
            print(fields[2])
            
        else:
            break
       
     
 
#print(count)