#!/usr/bin/env python3

import sys

Flybase_Uniprot = {}
fly = set()
 

#    For unknown Uniprot_IDs to be ignored, enter "B" as third argument,\n
#   if you would like these filled with "Unknown" enter "U" as the third argument""")
     
for line in open( sys.argv[1] ):
    fly = line.strip("\r\n").split()
    keyFlybase = fly[0]
    valueUni = fly[1]
    Flybase_Uniprot[keyFlybase] = valueUni
   # print(Flybase_Uniprot[keyFlybase])
   # print(Flybase_Uniprot)
   # print(Flybase_Uniprot.keys())

for line in open( sys.argv[2] ):
    fields = line.strip("\r\n").split("\t")
    Flybase_id = fields[8]
    if Flybase_id in Flybase_Uniprot:
        Uniprot_id = Flybase_Uniprot[Flybase_id]
        print(line + "\t" + Uniprot_id)
        #for unknown uniprotID to be filled in with "Unknown" add "U" as third argument
        if sys.argv[3] == "U":
            print(line + "\t" + "Unknown")
        #To ignore these 
        if sys.argv[3] == "B":
            continue
 