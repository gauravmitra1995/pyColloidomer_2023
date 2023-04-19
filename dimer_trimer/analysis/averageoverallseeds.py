import numpy as np
import os
import pickle
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob
import sys
import argparse

def round_to_int(arr):
    result=[]
    for i in range(arr.size):
        true=arr[i]
        nint=int(arr[i])
        if(np.abs(nint-true)<0.5):
            result.append(np.floor(true))
        else:
            result.append(np.ceil(true))
    return np.array(result)


parser=argparse.ArgumentParser()

#Parameters REQUIRED to be passed as arguments from outside
parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int)
parser.add_argument("--R",type=float)
parser.add_argument("--Np",type=int)

#Parameters not mandatory to be passed as arguments from outside
parser.add_argument("--metropolis",default=1,type=int)
parser.add_argument("--kAB",type=float,default=200.0)
parser.add_argument("--radiusB",type=float,default=1.0)
parser.add_argument("--kspring",type=float,default=10.0)
parser.add_argument("--dimension",type=int,default=2)
parser.add_argument("--gammaA",type=float,default=0.1)
parser.add_argument("--gammapatch",type=float,default=0.0001)
parser.add_argument("--r0",type=float,default=2.0)

args = parser.parse_args()
locals().update(vars(args))

if(epsilon!='infinite'):
    epsilon=float(epsilon)

if(Nclusters==2):
    simulationtype='dimer'
elif(Nclusters==3):
    simulationtype='trimer'


print("*******************************************************************************************************************************************************")
print("Np  R  epsilon  r0  gammaA  gammapatch simulationtype")
print(Np," ",R," ",epsilon," ",r0," ",gammaA," ",gammapatch," ",simulationtype)


timesteplist=[]
unbondedbinderslist_allseeds={}
seedlist=[]

timesteps_seeds=[]

#NUMBER OF BINDERS REMAINING ON DROPLET SURFACE 

for filename in sorted(glob.glob(str(os.getcwd())+"/../simulation_setup/"+str(simulationtype)+"/gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"/selfavoiding/Nclusters"+str(Nclusters)+"/Np"+str(Np)+"/R"+str(R)+"/radiusB"+str(radiusB)+"/kAB"+str(kAB)+"/epsilon"+str(epsilon)+"/dimension"+str(dimension)+"/kspring"+str(kspring)+"/seed*/restlength"+str(r0)+"/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(epsilon)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"_seed*.allruns.num_unbondedpatches.data"),key=lambda x:(int(((os.path.basename(x).split("_")[12]).split(".")[0]).replace('seed','')))):
    seed=int(((os.path.basename(filename).split("_")[12]).split(".")[0]).replace('seed',''))
    data = (np.load(filename,allow_pickle=True))    
    seedlist.append(seed)
    unbondedbinderslist_allseeds[seed]=data[0]
    timesteplist=data[1]
    timesteps_seeds.append(len(timesteplist))

right_length_array=max(timesteps_seeds)

mean_unbondedbinders_allseeds=[]
stddev_unbondedbinders_allseeds=[]

for i in range(len(data[0])):
    temp=[]
    for j in unbondedbinderslist_allseeds.keys():
        if(len(unbondedbinderslist_allseeds[j][i])==right_length_array):
            temp.append(unbondedbinderslist_allseeds[j][i])
    mean_unbondedbinders_allseeds.append(list(np.mean(temp,axis=0)))
    stddev_unbondedbinders_allseeds.append(list(np.std(temp,axis=0)))

output_file_unbondedbinders_allseeds=str(os.getcwd())+"/unbondedpatches_vs_time_data/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(epsilon)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".averageallseeds.allruns.num_unbondedpatches.data"

x=np.array([mean_unbondedbinders_allseeds,stddev_unbondedbinders_allseeds,timesteplist],dtype=object)
x.dump(output_file_unbondedbinders_allseeds)


