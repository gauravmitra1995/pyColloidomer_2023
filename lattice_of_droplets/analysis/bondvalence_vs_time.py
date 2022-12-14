import numpy as np
import os
import pickle 
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob 
import sys
import seaborn as sns
sns.set_context("poster", font_scale=1.0)
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("--Nclusters",type=int,default=81)
parser.add_argument("--R",type=float,default=50.0)
parser.add_argument("--radiusB",type=float,default=1.0)
parser.add_argument("--Np",type=int,default=100)
parser.add_argument("--areafraction",type=float,default=0.3)
parser.add_argument("--epsilon",default=20.7)
parser.add_argument("--gammaA",type=float,default=1.0)
parser.add_argument("--gammapatch",type=float,default=0.0001)
parser.add_argument("--kspring",type=float,default=10.0)
parser.add_argument("--restlength",type=float,default=2.0)

dt=0.001
args = parser.parse_args()
locals().update(vars(args))

dim1=int(np.sqrt(Nclusters))
dim2=dim1

if(epsilon!='infinite'):
    epsilon=float(epsilon)


fraction_allframes_allvalences=[]

filename=str(os.getcwd())+'/plot_data/bondvalencefraction_vs_time_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.txt'

data=np.genfromtxt(filename,comments='#')

for i in range(data.shape[1]):
    if(i==0):
        timelist=data[:,i]
    else:
        fraction_allframes_allvalences.append(data[:,i])

valences=np.arange(6)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


fig,ax=plt.subplots(figsize=(20,15),dpi=100)
for i in range(valences.size):
    plt.plot(timelist,np.array(fraction_allframes_allvalences)[i],color=colors[i],marker='o',markersize=12.0,linewidth=5.0,label=r'$B_{n}  = $'+str(valences[i]))
plt.xlabel('Simulation time in HOOMD units',labelpad=20,fontsize=60)
plt.ylabel('Fraction of droplets',labelpad=20,fontsize=60)
plt.legend(loc='upper right',ncol=2,prop={'size': 40})
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(5)
ax.tick_params(labelsize=50,axis='both', which='major', pad=20)
plt.yticks(fontsize=50)
plt.xticks(fontsize=50)
fig.tight_layout()
plt.savefig('final_figures/bondvalence_vs_time_Np'+str(Np)+'_R'+str(R)+'_areafrac'+str(areafraction)+'_gammaA'+str(gammaA)+'_eps'+str(epsilon)+'.png',bbox_inches='tight')
plt.close()


