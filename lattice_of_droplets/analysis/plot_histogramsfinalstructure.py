import numpy as np
import os
import pickle 
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob 
import sys
import seaborn as sns
import matplotlib.patches as mpatches
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

filename=str(os.getcwd())+'/histograms_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.histogramstructures.data'
data=np.load(filename,allow_pickle=True)
structures=data[0]
fractions=data[1]

structures=['monomer','dimer','chain','loop','other']

fig,ax=plt.subplots(figsize=(25,18),dpi=100)
plt.bar(structures, fractions, width=0.5,color=['darkgreen','crimson','darkmagenta','gold','turquoise'],edgecolor='black',linewidth=7.0)

"""
patch1 = mpatches.Patch(color='darkgreen', label='monomers')
patch2 = mpatches.Patch(color='crimson', label='dimers')
patch3 = mpatches.Patch(color='indigo', label='chains')
patch4 = mpatches.Patch(color='yellow', label='loops')
patch5 = mpatches.Patch(color='cyan',label='other')
"""
plt.xlabel('Structures formed',labelpad=20,fontsize=90)
plt.ylabel('Fraction of droplets',labelpad=20,fontsize=90)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(8)
ax.tick_params(labelsize=80,axis='both', which='major', length=20, width=10, pad=20)
plt.yticks(np.arange(0,np.max(fractions)+0.1,0.1),fontsize=80,fontweight='medium')
plt.xticks(fontsize=75,fontweight='medium')
#text="Key for structures"
#plt.legend(loc='best',title=text,title_fontsize=50,handles=[patch1,patch2,patch3,patch4,patch5],prop={'size':50},ncol=2)
fig.tight_layout()
plt.savefig('final_figures/histogram_structures_Np'+str(Np)+'_R'+str(R)+'_areafrac'+str(areafraction)+'_gammaA'+str(gammaA)+'_eps'+str(epsilon)+'.png',bbox_inches='tight')
plt.close()







