import os
import sys
import numpy as np
from gsd import hoomd as gsd
import argparse
import os.path
import glob
import hoomd
import hoomd.md

parser=argparse.ArgumentParser()
parser.add_argument("--fileprefix",type=str)
args = parser.parse_args()
locals().update(vars(args))

outputfile=fileprefix+".allruns.gsd"

hoomd.context.initialize("")
traj = gsd.open(name=outputfile, mode='wb')

k=0
num_files=4

for filename in sorted(glob.glob(fileprefix+".run.*.gsd"),key=lambda x:int(os.path.basename(x).split("_")[12].split('.')[2])):
    if(k>num_files-1):
        break
    trajectory = gsd.open(filename,'rb') # read gsd file
    key=int(os.path.basename(filename).split("_")[12].split('.')[2])
    print(key)
    print(len(trajectory))
    for i in range(len(trajectory)):
        if(k!=0 and i==0):
            continue
        system= trajectory[int(i)]
        traj.append(system) 
    k+=1

