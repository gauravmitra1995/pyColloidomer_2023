#%pylab inline
import seaborn as sns
import os
import numpy as np
from analyze_bonded_dimerortrimer import dynbond_gsd

import matplotlib.pyplot as plt
from functions_cluster_analysis import *

parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str)
args = parser.parse_args()
locals().update(vars(args))

output_file=os.path.splitext(os.path.basename(trajectory_file))[0]+".rotatedtraj.gsd"
modified_traj = gsd.open(name=output_file, mode='wb')

dynbond_pair_trj, Np, rotated_trj = dynbond_gsd(trajectory_file)
print("Number of frames:",len(dynbond_pair_trj))

for i in range(len(rotated_trj)):
    system= rotated_trj[int(i)]
    modified_traj.append(system)


