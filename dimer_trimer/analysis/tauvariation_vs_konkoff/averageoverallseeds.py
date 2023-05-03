import numpy as np
import os
import pickle
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob
import sys
import random
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
parser.add_argument("--epsilon",default=20.7)
parser.add_argument("--Nclusters",type=int,default=2)
parser.add_argument("--R",type=float,default=50.0)
parser.add_argument("--Np",type=int,default=100)
parser.add_argument("--koninit",type=float)
parser.add_argument("--koffinit",type=float)
parser.add_argument("--gammapatch",type=float,default=0.0001)


#Parameters not mandatory to be passed as arguments from outside
parser.add_argument("--metropolis",default=1,type=int)
parser.add_argument("--kAB",type=float,default=200.0)
parser.add_argument("--radiusB",type=float,default=1.0)
parser.add_argument("--kspring",type=float,default=10.0)
parser.add_argument("--dimension",type=int,default=2)
parser.add_argument("--gammaA",type=float,default=0.1)
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

for filename in sorted(glob.glob(str(os.getcwd())+"/data_averageoverseeds/gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"/epsilon"+str(epsilon)+"/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_kon"+str(koninit)+"_koff"+str(koffinit)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"_seed*.allruns.num_unbondedpatches.data"),key=lambda x:(int(((os.path.basename(x).split("_")[13]).split(".")[0]).replace('seed','')))):
    seed=int(((os.path.basename(filename).split("_")[13]).split(".")[0]).replace('seed',''))
    data = (np.load(filename,allow_pickle=True))
    if(np.any(np.array(data[0][0])>data[0][0][0])):   #discard any erroneous seeds (cases where the two droplets fall apart)
        print("Discarding error seed ",seed)
        continue
    seedlist.append(seed)
    unbondedbinderslist_allseeds[seed]=data[0]
    timesteplist=data[1]
    timesteps_seeds.append(len(timesteplist))

print("total no of seeds chosen: ",len(seedlist))

right_length_array=max(timesteps_seeds)

mean_unbondedbinders_allseeds=[]
stddev_unbondedbinders_allseeds=[]

#Perform bootstrapping of the data for the seeds present, to ensure better statistical averaging and error computation
n_samples=2000

temp=[]
for j in unbondedbinderslist_allseeds.keys():
    if(len(unbondedbinderslist_allseeds[j][0])==right_length_array):
        temp.append(unbondedbinderslist_allseeds[j][0])

bts=[]
for s in range(n_samples):
    r=random.randint(0,len(temp)-1)
    bts.append(np.array(temp[r]))
bts=np.array(bts)

output_file_unbondedbinders_allseeds=str(os.getcwd())+"/data_allseeds/gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"/epsilon"+str(epsilon)+"/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_kon"+str(koninit)+"_koff"+str(koffinit)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".allsamples.allruns.num_unbondedpatches.data"
x=np.array([bts,np.array(timesteplist)],dtype=object)
x.dump(output_file_unbondedbinders_allseeds)
