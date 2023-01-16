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

filename=str(os.getcwd())+'/histograms_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.histogramvalences.data'

data=np.load(filename,allow_pickle=True)
valences=data[0]
fractions=data[1]

fig,ax=plt.subplots(figsize=(25,18),dpi=100)
plt.bar(valences, fractions, width=0.5,color=['blue','darkorange','green','red','purple','brown'],edgecolor='black',linewidth=7.0)
plt.xlabel(r'$B_{n}$',labelpad=20,fontsize=100)
plt.ylabel(r'$P(B_{n})|^{t=t^{\infty}}$',labelpad=20,fontsize=100)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(8)
ax.tick_params(labelsize=80,axis='both', which='major', length=20, width=10, pad=20)
plt.yticks(np.arange(0,np.max(fractions)+0.1,0.1),fontsize=80,fontweight='medium')
plt.xticks([0,1,2,3,4,5],fontsize=80,fontweight='medium')
fig.tight_layout()
plt.savefig('final_figures/histogram_valences_Np'+str(Np)+'_R'+str(R)+'_areafrac'+str(areafraction)+'_gammaA'+str(gammaA)+'_eps'+str(epsilon)+'.png',bbox_inches='tight')
plt.close()


