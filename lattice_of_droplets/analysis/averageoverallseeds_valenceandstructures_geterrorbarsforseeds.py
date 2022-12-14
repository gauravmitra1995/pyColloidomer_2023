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


#BOND VALENCE ANALYSIS (Average the fraction of droplets in a given valence data for all the 10 seeds for a given condition and then save the average and SD at the final frame)

for filename in sorted(glob.glob(str(os.getcwd())+'/../simulation_setup/lattice/gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'/Nclusters'+str(Nclusters)+'/Np'+str(Np)+'/R'+str(R)+'/radiusB'+str(radiusB)+'/areafraction'+str(areafraction)+'/epsilon'+str(epsilon)+'/kspring'+str(kspring)+'/seed*/restlength'+str(restlength)+'/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'_seed*.allruns.bondvalencewithtime.data'),key=lambda x:int(((os.path.basename(x).split("_")[11]).split('.')[0]).replace('seed',''))):
    data=np.load(filename,allow_pickle=True)
    seed=int(((os.path.basename(filename).split("_")[11]).split('.')[0]).replace('seed',''))
    seedlist.append(seed)
    timelist=[]
    bondvalencesallframeslist=[]
    for i in range(data.shape[0]):
        timelist.append(data[i][0])
        bondvalencesallframeslist.append(data[i][1])
    timelist_allseeds[seed]=timelist
    timelength_allseeds.append(len(timelist))
    bondvalencesallframeslist_allseeds[seed]=bondvalencesallframeslist
    bondvalenceslastframe_allseeds[seed]=bondvalencesallframeslist[-1]

right_length_array=max(timelength_allseeds)

timelist_allseeds={i:timelist_allseeds[i] for i in timelist_allseeds if len(timelist_allseeds[i])==right_length_array}
bondvalencesallframeslist_allseeds={i:bondvalencesallframeslist_allseeds[i] for i in timelist_allseeds if len(timelist_allseeds[i])==right_length_array}
bondvalenceslastframe_allseeds={i:bondvalenceslastframe_allseeds[i] for i in timelist_allseeds if len(timelist_allseeds[i])==right_length_array}


bondvalencesallframes_allseeds=[]
for s in seedlist:
    bondvalencesallframes_eachseed=np.array(bondvalencesallframeslist_allseeds[s])
    bondvalencesallframes_allseeds.append(bondvalencesallframes_eachseed)

bondvalencesallframes_allseeds=np.array(bondvalencesallframes_allseeds)

valences=np.arange(7)

fractions_lastframe_allseeds=[]

results_dict = {}

for j in range(bondvalencesallframes_allseeds.shape[0]):
    fraction_allframes_allvalences=[]
    for k in valences:
        fraction_allframes=[]
        for i in range(bondvalencesallframes_allseeds.shape[1]):
            frame=i
            bonds_givenframe=[]
            bonds_givenframe=(bondvalencesallframes_allseeds[j][i]).astype(int)
            ctr=0
            for x in range(len(bonds_givenframe)):
                if(bonds_givenframe[x]==valences[k]):
                    ctr=ctr+1
            fraction_allframes.append(ctr/len(bonds_givenframe))    
        fraction_allframes_allvalences.append(fraction_allframes) 
    fraction_allframes_allvalences=np.array(fraction_allframes_allvalences,dtype=object)
    transpose=fraction_allframes_allvalences.T
    fractions_lastframe=transpose[-1]
    fractions_lastframe_allseeds.append(fractions_lastframe)

fractions_lastframe_allseeds=np.array(fractions_lastframe_allseeds,dtype=object).astype(float)
meanfraction_lastframe_allvalences=np.round(np.mean(fractions_lastframe_allseeds,axis=0),3)
stddevfraction_lastframe_allvalences=np.round(np.std(fractions_lastframe_allseeds,axis=0),3)

print("Valences:",valences)
print(fractions_lastframe_allseeds)

output_filevalences_allseeds=str(os.getcwd())+'/finalframefractions_valences_averageallseeds/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.averageallseeds.lastframe.fractionvalences.data'

y=np.array([valences,meanfraction_lastframe_allvalences,stddevfraction_lastframe_allvalences])
y.dump(output_filevalences_allseeds)

print("---------------------------------------------------------------------------------------------------------------------------------------------")

#STRUCTURAL CLASSIFICATION (SUCCESS: Dimers,chains,loops  &   ERRORS: Monomers,other)

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

meanfraction_lastframe_allstructures=[]
stddevfraction_lastframe_allstructures=[]

fractions_lastframe_allseeds=[]
fractions_lastframe_binaryclassifier_allseeds=[]

meanfraction_lastframe_allbinaryclassifiers=[]
stddevfraction_lastframe_allbinaryclassifiers=[]

binaryclassifier=['success','errors']

for seed in list(fraction_structures_allseeds_dict.keys()):
   fraction_allframes_allstructures=[]
   for s in structures:
       fraction_givenstructure_allframes=np.array(fraction_structures_allseeds_dict[seed][s])
       fraction_allframes_allstructures.append(fraction_givenstructure_allframes)
   fraction_allframes_allstructures=np.array(fraction_allframes_allstructures)
   trans=fraction_allframes_allstructures.T
   fractions_lastframe=trans[-1]
   fractions_lastframe_allseeds.append(fractions_lastframe)

   fractions_lastframe_binaryclassifier_allseeds.append([1.0-(fractions_lastframe[0]+fractions_lastframe[-1]),(fractions_lastframe[0]+fractions_lastframe[-1])])

fractions_lastframe_allseeds=np.array(fractions_lastframe_allseeds)
fractions_lastframe_binaryclassifier_allseeds=np.array(fractions_lastframe_binaryclassifier_allseeds)

print("Structures:",structures)
print(fractions_lastframe_allseeds)

meanfraction_lastframe_allstructures=np.mean(fractions_lastframe_allseeds,axis=0)
stddevfraction_lastframe_allstructures=np.std(fractions_lastframe_allseeds,axis=0)

print(meanfraction_lastframe_allstructures)
print(stddevfraction_lastframe_allstructures)

print("---------------------------------------------------------------------------------------------------------------------------------------------")
print("Binary classifier for optimization of structures:",binaryclassifier)
print(fractions_lastframe_binaryclassifier_allseeds)

meanfraction_lastframe_allbinaryclassifiers=np.mean(fractions_lastframe_binaryclassifier_allseeds,axis=0)
stddevfraction_lastframe_allbinaryclassifiers=np.std(fractions_lastframe_binaryclassifier_allseeds,axis=0)

output_filestructures_allseeds=str(os.getcwd())+'/finalframefractions_structures_averageallseeds/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.averageallseeds.lastframe.fractionstructures.data'

z=np.array([structures,meanfraction_lastframe_allstructures,stddevfraction_lastframe_allstructures])
z.dump(output_filestructures_allseeds)

output_filebinarystr_allseeds=str(os.getcwd())+'/finalframefractions_structures_averageallseeds/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.averageallseeds.lastframe.fractionbinaryclassifiers.data'

z=np.array([binaryclassifier,meanfraction_lastframe_allbinaryclassifiers,stddevfraction_lastframe_allbinaryclassifiers])
z.dump(output_filebinarystr_allseeds)

