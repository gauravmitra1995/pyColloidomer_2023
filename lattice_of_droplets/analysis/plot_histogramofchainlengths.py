import numpy as np
import os
import pickle 
import random
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob 
import sys
import seaborn as sns
sns.set_context("poster", font_scale=1.0)
import argparse
import pandas as pd
from functions_cluster_analysis import *
from scipy.optimize import curve_fit

def func(x,a,b,c):
    return a * (np.exp(-b * (x-3))) + c

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

filename=str(os.getcwd())+'/chainlengths_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.finalframefractionchainsizes.data'

data=np.load(filename,allow_pickle=True)
chainsizes=data[0]
mean_fractions_allchainsizes=data[1]
stddev_fractions_allchainsizes=data[2]

#Load experimental chain length data from .csv file (for concentration of DNA = 0.3 pmol, McMullen et al. PRL 2018, Fig 3(a))
csvfile=str(os.getcwd())+'/expt_data_PRL2018/chainlength_data_PRL2018.csv'
df = pd.read_csv(csvfile,delimiter=' ')
chainsizes_expt=df.iloc[:, 0].tolist()
fractions_expt=df.iloc[:, 1].tolist()

#Plotting the chain length distributions from simulation:
fig,ax=plt.subplots(figsize=(25,18),dpi=100)
plt.bar(chainsizes,mean_fractions_allchainsizes,yerr=stddev_fractions_allchainsizes,width = 1,color='royalblue',linewidth=7.0,error_kw={'elinewidth':8.0},ecolor="black",edgecolor="black",capsize=12)

#fit simulation data to exponential function:
a=mean_fractions_allchainsizes[0]
b=1/4
c=0.0

popt, pcov = curve_fit(func,chainsizes,mean_fractions_allchainsizes)
a=popt[0]
b=popt[1]
c=popt[2]

plt.plot(chainsizes,func(chainsizes,*popt),linestyle='--',color='navy',linewidth=10.0,label='Fit to Simulation')
plt.plot(chainsizes_expt,fractions_expt,marker='s',markersize=40,linestyle='',markerfacecolor='none',markeredgecolor='k',label='Experiment',markeredgewidth=10.0)

ax.set_xlabel('Linear chain lengths '+r'$(N)$',labelpad=20,fontsize=90)
ax.set_ylabel(r'$P(N)$',labelpad=20,fontsize=100)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(8)
ax.tick_params(labelsize=80,axis='both', which='major', length=20, width=10, pad=20)
plt.xticks(np.arange(3,np.max(chainsizes)+1,2),fontsize=80,fontweight='medium')
plt.yticks(np.arange(0,np.max(mean_fractions_allchainsizes)+0.1,0.1),fontsize=80,fontweight='medium')
plt.legend(loc='best',ncol=1,prop={'size': 75},fontsize=75)
fig.tight_layout()
plt.savefig('final_figures/histogram_chainlengths_Np'+str(Np)+'_R'+str(R)+'_areafrac'+str(areafraction)+'_gammaA'+str(gammaA)+'_eps'+str(epsilon)+'.png',bbox_inches='tight')
plt.close()

