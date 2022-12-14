import numpy as np
from functions_cluster_analysis import *
import matplotlib.pyplot as plt
import os
import sys
import glob
import argparse
import pickle
from gsd import hoomd as gsd
import seaborn as sns
sns.set_context("poster", font_scale=1.0)


parser=argparse.ArgumentParser()
parser.add_argument("--fileprefix",type=str)
parser.add_argument("--dt",type=float,default=0.0005)
parser.add_argument("--num_cycles",type=int,default=15)
args = parser.parse_args()
locals().update(vars(args))

temperatures=[]
timesteps=[]

#Total 15 cycles of folding/unfolding for this long trajectory 

k=0
for filename in sorted(glob.glob("folding_singletrajectory_forpaper/logfiles_singletrajectory/"+str(fileprefix)+".run.*.log"),key=lambda x:int(os.path.basename(x).split("_")[12].split('.')[2])):
    if(k==int(num_cycles*2)):
        continue
    key=int(os.path.basename(filename).split("_")[12].split('.')[2])
    data = np.genfromtxt(fname=filename, skip_header=True)
    if(k!=0):
        temperatures.append(list(data[:,2][1:]))
        timesteps.append(list(data[:,0][1:]))
    else:
        temperatures.append(list(data[:,2][0:]))
        timesteps.append(list(data[:,0][0:]))
    k+=1


timesteps=np.array(timesteps,dtype=object)
temperatures=np.array(temperatures,dtype=object)

Time=[]
Temperatures=[]

for i in range(temperatures.shape[0]):
    Temperatures.extend(temperatures[i])
    Time.extend(np.array(timesteps[i])*dt)

Time=Time[::10]  #because the log files are written 10 times more frequently than the dump frequency for the .gsd files
Temperatures=Temperatures[::10]

#Get the average valence data 
filename="folding_singletrajectory_forpaper/"+str(fileprefix)+".allruns.averagevalence.data"
data=np.load(filename,allow_pickle=True)
valences=data[1][:]

#Plot temperature and average valence 
fig,ax1=plt.subplots(figsize=(20,15),dpi=100)
ax1.plot(Time[int(len(Time)/2):],Temperatures[int(len(Time)/2):],marker='.',linewidth=7.0,color='darkblue',label=r'$k_{B}T$')
ax2=ax1.twinx()
ax2.plot(Time[int(len(Time)/2):],valences[int(len(Time)/2):],marker='.',linewidth=7.0,color='crimson',label=r'$\langle B_{n} \rangle$')
ax1.legend(loc='lower left',ncol=1,prop={'size': 50})
ax1.set_yticks(np.arange(0.0,1.7,0.1))
ax2.set_yticks(np.arange(0.0,3.50,0.2))
ax2.legend(loc='lower right',ncol=1,prop={'size': 40})
ax1.set_xlabel("Simulation time (in HOOMD units)",fontsize=60,labelpad=20)
ax1.set_ylabel("Temperature (in units of $k_{B} T$)",fontsize=60,labelpad=20)
ax2.set_ylabel(r'$\langle B_{n} \rangle$',fontsize=60,labelpad=20)
ax1.tick_params(axis='y',colors='darkblue',labelsize=50)
ax2.tick_params(axis='y',colors='crimson',labelsize=50)
ax1.tick_params(axis='x',colors='black',labelsize=40)
ax1.yaxis.label.set_color('darkblue')
ax2.yaxis.label.set_color('crimson')
ax1.spines['left'].set_color('darkblue')
ax1.spines['right'].set_color('crimson')
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(5)
fig.tight_layout()
plt.savefig('final_figures/folding_unfolding_cycle.png',bbox_inches='tight')
plt.close()


