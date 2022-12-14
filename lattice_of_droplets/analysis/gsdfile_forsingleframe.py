import os
import sys
import numpy as np
from gsd import hoomd as gsd
import argparse
import os.path
import hoomd
import hoomd.md

parser=argparse.ArgumentParser()
parser.add_argument("--fileprefix",type=str)
args = parser.parse_args()
locals().update(vars(args))

framenumber=-1

hoomd.context.initialize("")

trajectory_file=fileprefix+'.allruns.gsd'
trajectory = gsd.open(trajectory_file,'rb') # read gsd file
l=len(trajectory)
print("No of frames: ",l)

basename=os.path.basename(fileprefix)
singleframe_file='finalframes/'+str(basename)+'.finalframe.gsd'

system=hoomd.init.read_gsd(trajectory_file,frame=framenumber)
all = hoomd.group.all()
hoomd.dump.gsd(filename=singleframe_file,period=None,group=all, overwrite=True,dynamic=['attribute','momentum','topology'])


