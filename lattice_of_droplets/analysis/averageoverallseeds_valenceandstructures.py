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

#BOND VALENCE ANALYSIS (Combine valence data from all the 10 seeds for a given condition) 

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


bondvalencesallseeds_allframes=[]

for i in range(len(timelist)):
    bondvalencesallseeds_eachframe=[]
    for j in timelist_allseeds.keys():
        bondvalencesallseeds_eachframe.append(bondvalencesallframeslist_allseeds[j][i])
    bondvalencesallseeds_allframes.append(np.array(flat_list(bondvalencesallseeds_eachframe)))
bondvalencesallseeds_allframes=np.array(bondvalencesallseeds_allframes)


output_file2_allseeds=str(os.getcwd())+'/valence_data/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.bondvalencewithtime.data'

y=bondvalencesallseeds_allframes
y.dump(output_file2_allseeds)

valences=np.arange(7)
fraction_allframes_allvalences=[]

results_dict = {}

for k in valences:
    fraction_allframes=[]
    for i in range(bondvalencesallseeds_allframes.shape[0]):
        frame=i
        bonds_givenframe=[]
        bonds_givenframe=(bondvalencesallseeds_allframes[i]).astype(int)
        ctr=0
        for x in range(len(bonds_givenframe)):
            if(bonds_givenframe[x]==valences[k]):
                ctr=ctr+1
        fraction_allframes.append(ctr/len(bonds_givenframe))    
    fraction_allframes_allvalences.append(fraction_allframes) 


fraction_allframes_allvalences=np.array(fraction_allframes_allvalences,dtype=object)
timelist=np.array(timelist,dtype=object)

res=[]
for i in range((np.array(fraction_allframes_allvalences).T).shape[1]):
    X=(np.array(fraction_allframes_allvalences)[i].T)
    res.append(X)

res=np.array(res).astype(float)
time=(timelist).astype(float)

Y = np.concatenate((time.reshape(time.shape[0],1),res.T),axis=-1)
np.savetxt(fname=str(os.getcwd())+'/plot_data/bondvalencefraction_vs_time_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.txt',X=Y,header="Time     Valence0    Valence1     Valence2     Valence3     Valence4     Valence5     Valence6")


#STRUCTURAL ANALYSIS (Combine structure classification data from all the 10 seeds for a given condition)

structures=['monomers','dimers','linear chains','loops','other']

seedlist=[]
fraction_structures_allseeds_dict={}
timelist_allseeds={}
timelength_allseeds=[]
for filename in sorted(glob.glob(str(os.getcwd())+'/../simulation_setup/lattice/gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'/Nclusters'+str(Nclusters)+'/Np'+str(Np)+'/R'+str(R)+'/radiusB'+str(radiusB)+'/areafraction'+str(areafraction)+'/epsilon'+str(epsilon)+'/kspring'+str(kspring)+'/seed*/restlength'+str(restlength)+'/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'_seed*.structuralanalysis.csv'),key=lambda x:int(((os.path.basename(x).split("_")[11]).split('.')[0]).replace('seed',''))):
    seed=int(((os.path.basename(filename).split("_")[11]).split('.')[0]).replace('seed',''))
    #print("seed:",seed)
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

mean_fraction_structures={}
for s in structures:
    fraction_givenstucture_allseeds=[]
    for seed in list(fraction_structures_allseeds_dict.keys()):
        fraction=np.array(fraction_structures_allseeds_dict[seed][s])
        fraction_givenstucture_allseeds.append(fraction)
    fraction_givenstucture_allseeds=np.array(fraction_givenstucture_allseeds)
    mean_fraction_givenstructure_overallseeds=np.round(np.mean(fraction_givenstucture_allseeds,axis=0),3)
    mean_fraction_structures[s]=list(mean_fraction_givenstructure_overallseeds)

fraction_allframes_allstructures=[]

for s in structures:
    fraction_allframes=mean_fraction_structures[s]
    fraction_allframes_allstructures.append(fraction_allframes)

fraction_allframes_allstructures=np.array(fraction_allframes_allstructures,dtype=object)
timelist=timelist_allseeds[list(fraction_structures_allseeds_dict.keys())[0]]
timelist=np.array(timelist,dtype=object)

res=[]
for i in range((np.array(fraction_allframes_allstructures).T).shape[1]):
    X=(np.array(fraction_allframes_allstructures)[i].T)
    res.append(X)

res=np.array(res).astype(float)
time=(timelist).astype(float)

Y = np.concatenate((time.reshape(time.shape[0],1),res.T),axis=-1)
np.savetxt(fname=str(os.getcwd())+'/structuralanalysis_data/clusterstructures_vs_time_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.txt',X=Y,header="Time     Monomers    Dimers     Chains     Loops     Gels/Branchedchains")


#AVERAGE VALENCE ANALYSIS (Average the already obtained average valences further over all 10 seeds for a given condition and find the standard deviation)

averagevalenceallframeslist_allseeds={}
timelength_allseeds=[]
timelist_allseeds={}
seedlist=[]

for filename in sorted(glob.glob(str(os.getcwd())+'/../simulation_setup/lattice/gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'/Nclusters'+str(Nclusters)+'/Np'+str(Np)+'/R'+str(R)+'/radiusB'+str(radiusB)+'/areafraction'+str(areafraction)+'/epsilon'+str(epsilon)+'/kspring'+str(kspring)+'/seed*/restlength'+str(restlength)+'/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'_seed*.allruns.averagevalence.data'),key=lambda x:int(((os.path.basename(x).split("_")[11]).split('.')[0]).replace('seed',''))):
    data=np.load(filename,allow_pickle=True)
    seed=int(((os.path.basename(filename).split("_")[11]).split('.')[0]).replace('seed',''))
    seedlist.append(seed)
    timelist_allseeds[seed]=data[0]
    timelength_allseeds.append(len(data[0]))
    averagevalenceallframeslist_allseeds[seed]=data[1]

right_length_array=max(timelength_allseeds)

timelist_allseeds={i:timelist_allseeds[i] for i in timelist_allseeds if len(timelist_allseeds[i])==right_length_array}
averagevalenceallframeslist_allseeds={i:averagevalenceallframeslist_allseeds[i] for i in timelist_allseeds if len(timelist_allseeds[i])==right_length_array}


temp_averagevalenceallframes=[]
temp_time=[]
for j in timelist_allseeds.keys():
    temp_time.append(timelist_allseeds[j])
    temp_averagevalenceallframes.append(averagevalenceallframeslist_allseeds[j])

mean_time_allseeds=np.mean(temp_time,axis=0)
mean_averagevalenceallframes_allseeds=np.mean(np.array(temp_averagevalenceallframes,dtype=np.float64),axis=0)
stddev_averagevalenceallframes_allseeds=np.std(np.array(temp_averagevalenceallframes,dtype=np.float64),axis=0)

print("Average valences over all seeds:")
print(mean_averagevalenceallframes_allseeds)


output_file3_allseeds=str(os.getcwd())+'/valence_data/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.averageallseeds.allruns.averagevalence.data'

z=np.array([mean_averagevalenceallframes_allseeds,stddev_averagevalenceallframes_allseeds,mean_time_allseeds],dtype=object)
z.dump(output_file3_allseeds)

