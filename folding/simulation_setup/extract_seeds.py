import numpy as np
import sys
import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("--fileprefix",type=str)
args = parser.parse_args()
locals().update(vars(args))

with open(fileprefix) as f:
    lines = f.read().split()

if(len(lines)>5):
    seedlist1=""
    for i in range(int(len(lines)/2)):
        seedlist1=seedlist1+lines[i]+" "

    seedlist2=""
    for i in range(int(len(lines)/2),len(lines)):
        seedlist2=seedlist2+lines[i]+" "

    print(seedlist1)
    print(seedlist2)

else:
    seedlist=""
    for i in range(len(lines)):
        seedlist=seedlist+lines[i]+" "
    print(seedlist)





