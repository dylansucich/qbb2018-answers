#!/usr/bin/env python3

import sys 

if len(sys.argv) > 1:
    f = open( sys.argv[1] )
else:
    f = sys.stdin

count = 0

for i, line in enumerate( f ):
    if line[0] == "@":
        continue
    if "NH:i:1" in line:
        count += 1
        
print(count)