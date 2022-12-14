import numpy as np
import os
import pickle 
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
#parser.add_argument("--areafraction",type=float,default=0.3)
parser.add_argument("--epsilon",default=20.7)
#parser.add_argument("--gammaA",type=float,default=1.0)
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

labellist=[]
meanfractions_success=[]
stddevfractions_success=[]
meanfractions_errors=[]
stddevfractions_errors=[]
meanfractions_valence2=[]
stddevfractions_valence2=[]

#Obtain the avg and SD of fraction of droplets with valence 2 at final frame for the 8 different (phi,gamma) pairs:

for filename in sorted(glob.glob(str(os.getcwd())+'/finalframefractions_valences_averageallseeds/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac*_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA*_gammapatch'+str(gammapatch)+'.averageallseeds.lastframe.fractionvalences.data'),key=lambda x:(float(os.path.basename(x).split("_")[6].replace('areafrac','')),float(os.path.basename(x).split("_")[9].replace('gammaA','')))):
    data=np.load(filename,allow_pickle=True)
    areafraction=float(os.path.basename(filename).split("_")[6].replace('areafrac',''))
    gammaA=float(os.path.basename(filename).split("_")[9].replace('gammaA',''))
    valences=data[0]
    meanfractionvalences=data[1].astype(float)
    stddevfractionvalences=data[2].astype(float)
    label=(areafraction,gammaA)
    labellist.append(label)
    meanfractions_valence2.append(meanfractionvalences[2])
    stddevfractions_valence2.append(stddevfractionvalences[2])

#Obtain the avg and SD of fraction of droplets in "Success" at final frame for the 8 different (phi,gamma) pairs:

for filename in sorted(glob.glob(str(os.getcwd())+'/finalframefractions_structures_averageallseeds/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac*_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA*_gammapatch'+str(gammapatch)+'.averageallseeds.lastframe.fractionbinaryclassifiers.data'),key=lambda x:(float(os.path.basename(x).split("_")[6].replace('areafrac','')),float(os.path.basename(x).split("_")[9].replace('gammaA','')))):
    data=np.load(filename,allow_pickle=True)
    areafraction=float(os.path.basename(filename).split("_")[6].replace('areafrac',''))
    gammaA=float(os.path.basename(filename).split("_")[9].replace('gammaA',''))
    structures=data[0]
    meanfractions=data[1].astype(float)
    stddevfractions=data[2].astype(float)
    label=(areafraction,gammaA)
    meanfractions_errors.append(meanfractions[1])
    stddevfractions_errors.append(stddevfractions[1])
    meanfractions_success.append(meanfractions[0])
    stddevfractions_success.append(stddevfractions[0])

k=0
conditions = ['']*len(labellist)
for l in labellist:
    conditions[k]=str(l[0])+','+'\n'+str(l[1])
    k=k+1

#Obtain the avg and SD of maximum chain lengths at final frame for the 8 different (phi,gamma) pairs:

averagemaxchainlengths=[]
stddevmaxchainlengths=[]

for filename in sorted(glob.glob(str(os.getcwd())+'/chainlengths_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac*_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA*_gammapatch'+str(gammapatch)+'.finalframechainlengthstatistics.data'),key=lambda x:(float(os.path.basename(x).split("_")[6].replace('areafrac','')),float(os.path.basename(x).split("_")[9].replace('gammaA','')))):
    data=np.load(filename,allow_pickle=True)
    areafraction=float(os.path.basename(filename).split("_")[6].replace('areafrac',''))
    gammaA=float(os.path.basename(filename).split("_")[9].replace('gammaA',''))
    label=(areafraction,gammaA)
    averagemaxchainlengths.append(data[0])
    stddevmaxchainlengths.append(data[1])

#Plot all 3 metrics in bar charts for every (phi,gamma) condition:  

for l in range(len(labellist)):
    fig,ax=plt.subplots(figsize=(21,9),dpi=100)
    info1=[meanfractions_success[l],meanfractions_valence2[l],averagemaxchainlengths[l]/np.max(averagemaxchainlengths)]
    errorinfo1=[stddevfractions_success[l],stddevfractions_valence2[l],stddevmaxchainlengths[l]/np.max(averagemaxchainlengths)]
    xlabels1=['Fraction in'+'\n'+'valence 2','Fraction in '+'\n'+'success',r'$N_{max}/N_{max}^{highest}$']
    N=len(info1)
    ind=np.arange(N)
    xbuffer = 3.0
    height=1.0
    plt.xticks(np.arange(0.0,1.1,0.2))
    ax.set_yticks([])
    ax.barh(ind,info1,height=height,xerr=errorinfo1,color=['darkmagenta','green','royalblue'],linewidth=9.0,edgecolor="black",ecolor='black',error_kw={'elinewidth':12.0},capsize=20.0)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.tick_params(labelsize=60,labelleft=False,labelbottom=False,axis='both', which='major', pad=20)
    ax.set_xlim(0.0,1.2)
    ax.invert_yaxis()
    plt.grid(alpha=1.0,color='k',linestyle='--',linewidth=5.0) 

    #patch1 = mpatches.Patch(color='darkmagenta', label='Fraction in success')
    #patch2 = mpatches.Patch(color='green', label='Fraction with '+r'$B_n = 2$')
    #patch3 = mpatches.Patch(color='royalblue', label=r'${\langle N_{max} \rangle} / {\langle N_{max} \rangle}^{highest}$')
    #plt.legend(loc='best',title="Key for metrics",title_fontsize=50,handles=[patch1,patch2,patch3],prop={'size':50})

    fig.tight_layout()
    plt.savefig('final_figures/barcharts_optimizations_Np'+str(Np)+'_R'+str(R)+'_areafrac'+str(labellist[l][0])+'_gammaA'+str(labellist[l][1])+'_eps'+str(epsilon)+'.png',bbox_inches='tight')
    plt.close()
