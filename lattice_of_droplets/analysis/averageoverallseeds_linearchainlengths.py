import numpy as np
import os
import pickle
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob
import sys
import math
import argparse
import pandas as pd
import seaborn as sns
import random
from functions_cluster_analysis import *
sns.set_context("poster", font_scale=1.0)
from scipy.optimize import curve_fit

def flat_list(alist): # [[1,2],[3,4]] --> [1,2,3,4]
    flatlist=list()
    for sublist in alist:
        for i in sublist:
            flatlist.append(i)
    return(flatlist)


def round_to_int(arr):
    result=[]
    for i in range(arr.size):
        true=arr[i]
        nint=int(arr[i])
        if(np.abs(nint-true)<0.5):
            result.append(int(np.floor(true)))
        else:
            result.append(int(np.ceil(true)))
    return np.array(result)

def bootstrap(n_samples,a):
    bootstrap_samples_maxchainlength=[]
    for i in range(n_samples):
        bts=[]
        for elem in a:
            r = random.randint(0,len(a)-1)
            bts.append(a[r])
        bootstrap_samples_maxchainlength.append(bts)

    return bootstrap_samples_maxchainlength


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

args = parser.parse_args()
locals().update(vars(args))

dt=0.001
dumptime=100
dumpperiod=int(dumptime/dt)

if(epsilon!='infinite'):
    epsilon=float(epsilon)

dim1=int(np.sqrt(Nclusters))
dim2=dim1

bondvalencesallframeslist_allseeds={}
bondvalenceslastframe_allseeds={}
timelist_allseeds={}
timelength_allseeds=[]
seedlist=[]

print("*******************************************************************************************************************************************************")
print("Np  R  areafraction  epsilon  restlength  gammaA  gammapatch")
print(Np," ",R," ",areafraction," ",epsilon," ",restlength," ",gammaA," ",gammapatch)



linearchainlengths_allseeds={}
timelength_allseeds=[]
timelist_allseeds={}
seedlist=[]

linearchainlengths_allseedscombined_finalframe=[]
maxlengths_allseeds=[]
meanlengths_allseeds=[]


##### Part 1: Obtain average maximum chain length and SD #####

n_b=500
for filename in sorted(glob.glob(str(os.getcwd())+'/../simulation_setup/lattice/gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'/Nclusters'+str(Nclusters)+'/Np'+str(Np)+'/R'+str(R)+'/radiusB'+str(radiusB)+'/areafraction'+str(areafraction)+'/epsilon'+str(epsilon)+'/kspring'+str(kspring)+'/seed*/restlength'+str(restlength)+'/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'_seed*.allruns.linearchainlengths.data'),key=lambda x:int(((os.path.basename(x).split("_")[11]).split('.')[0]).replace('seed',''))):
    data=np.load(filename,allow_pickle=True)
    seed=int(((os.path.basename(filename).split("_")[11]).split('.')[0]).replace('seed',''))
    seedlist.append(seed)
    if(len(list(data))!=0):
        linearchainlengths_allseeds[seed]=list(data[-1])
        linearchainlengths_allseedscombined_finalframe.append(list(data[-1]))
        maxlength=np.max(data[-1])
        meanlength=np.mean(data[-1])
        maxlengths_allseeds.append(maxlength)
        meanlengths_allseeds.append(meanlength)

#Bootstrapping the maximum chain length data obtained from all the 10 seeds

bootstrap_samples_maxchainlength=bootstrap(n_b,maxlengths_allseeds)
max_allsamples=[]
for b in bootstrap_samples_maxchainlength:
    max_eachsample=np.max(b)
    max_allsamples.append(max_eachsample)

#Calculating the average maximum chain length and the standard deviation after obtaining all the bootstrap samples 
averagemaxlength=np.mean(max_allsamples)
stddevmaxlength=np.std(max_allsamples)
print("Average max chain length and error is: ",averagemaxlength,stddevmaxlength)

output_file_maxchainlengthinfo=str(os.getcwd())+'/chainlengths_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.finalframechainlengthstatistics.data'
y=np.array([averagemaxlength,stddevmaxlength])
y.dump(output_file_maxchainlengthinfo)



##### Part 2: Obtain chain length distributions by bootstrapping and saving the average and SD of fractions of various chain sizes#####

a=linearchainlengths_allseedscombined_finalframe
print("Chain lengths for all seeds combined:")
print(flat_list(a))

#Bootstrapping the chain length data obtained from all the 10 seeds and calculating the fraction of each chain size (starting from 3) 

n_samples=100
lower=min(flat_list(a))
upper=max(flat_list(a))
chainsizes=np.arange(lower,upper+1,1)

fractions_allsamples=[]

for i in range(n_samples):
    bts=[]
    for elem in a:
        r = random.randint(0,len(a)-1)
        bts.append(a[r])
    bootstrap_samples=flat_list(bts)
    counts=np.zeros(len(chainsizes))
    fractions=np.zeros(len(chainsizes))
    for c in range(len(chainsizes)):
        ctr=0
        for b in bootstrap_samples:
            if(b==chainsizes[c]):
                ctr+=1
        counts[c]=ctr
        fractions[c]=np.round(ctr/len(bootstrap_samples),4)
    fractions_allsamples.append(fractions)

fractions_allsamples=np.array(fractions_allsamples)
fractions_allchainsizes=fractions_allsamples.T

mean_fractions_allchainsizes=np.array([np.mean(ele) for ele in fractions_allchainsizes])
stddev_fractions_allchainsizes=np.array([np.std(ele) for ele in fractions_allchainsizes])

output_file_fractionchainsizes=str(os.getcwd())+'/chainlengths_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.finalframefractionchainsizes.data'
z=np.array([chainsizes,mean_fractions_allchainsizes,stddev_fractions_allchainsizes])
z.dump(output_file_fractionchainsizes)
