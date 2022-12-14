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
sns.set_context("poster", font_scale=1.0)


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


#histogram of final structures 

structures=['monomers','dimers','linear chains','loops','other']

seedlist=[]
fraction_structures_allseeds_dict={}
timelist_allseeds={}
timelength_allseeds=[]
for filename in sorted(glob.glob(str(os.getcwd())+'/../simulation_setup/lattice/gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'/Nclusters'+str(Nclusters)+'/Np'+str(Np)+'/R'+str(R)+'/radiusB'+str(radiusB)+'/areafraction'+str(areafraction)+'/epsilon'+str(epsilon)+'/kspring'+str(kspring)+'/seed*/restlength'+str(restlength)+'/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'_seed*.structuralanalysis.csv'),key=lambda x:int(((os.path.basename(x).split("_")[11]).split('.')[0]).replace('seed',''))):
    seed=int(((os.path.basename(filename).split("_")[11]).split('.')[0]).replace('seed',''))
    seedlist.append(seed)
    df = pd.read_csv(filename,delimiter=' ')
    fraction_structures_givenseed_dict=df.to_dict('list')

    l=len(fraction_structures_givenseed_dict[structures[0]])
    frames_list=np.arange(l)
    time_list=frames_list*dumpperiod*dt
    timelist_allseeds[seed]=time_list
    timelength_allseeds.append(l)
    fraction_structures_allseeds_dict[seed]=fraction_structures_givenseed_dict


right_length_array=max(timelength_allseeds)
timelist_allseeds={i:timelist_allseeds[i] for i in timelist_allseeds if len(timelist_allseeds[i])==right_length_array}


fraction_structures_allseeds_dict_noerrorseeds={}

for k in fraction_structures_allseeds_dict.keys():
    fraction_structures_givenseed=fraction_structures_allseeds_dict[k]
    l=len(fraction_structures_givenseed[structures[0]])
    if(l==right_length_array):
        fraction_structures_allseeds_dict_noerrorseeds[k]=fraction_structures_givenseed

fraction_structures_allseeds_dict=fraction_structures_allseeds_dict_noerrorseeds

fractions=[]
for s in structures:
    numberdroplets_givenstructure_allseeds=[]
    for seed in list(fraction_structures_allseeds_dict.keys()):
        fraction=np.array(fraction_structures_allseeds_dict[seed][s])
        nodroplets=round_to_int(fraction*Nclusters)
        numberdroplets_givenstructure_allseeds.append(nodroplets)

    numberdroplets_givenstructure_allseeds=np.array(numberdroplets_givenstructure_allseeds)
    numberdroplets_givenstructure_allseedscombined=numberdroplets_givenstructure_allseeds.T
    fraction_givenstructure_allseedscombined=np.round((np.sum(numberdroplets_givenstructure_allseedscombined[-1])/(Nclusters*numberdroplets_givenstructure_allseedscombined.shape[1])),4)
    fractions.append(fraction_givenstructure_allseedscombined)

print(structures)
print(fractions)

output_file=str(os.getcwd())+'/histograms_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.histogramstructures.data'

y=np.array([structures,fractions],dtype=object)
y.dump(output_file)



