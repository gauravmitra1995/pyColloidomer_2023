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

print("*******************************************************************************************************************************************************")
print("Np  R  areafraction  epsilon  restlength  gammaA  gammapatch")
print(Np," ",R," ",areafraction," ",epsilon," ",restlength," ",gammaA," ",gammapatch)


#histogram of final valences

fraction_allframes_allvalences=[]

filename=str(os.getcwd())+'/plot_data/bondvalencefraction_vs_time_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.txt'

data=np.genfromtxt(filename,comments='#')

for i in range(data.shape[1]):
    if(i==0):
        timelist=data[:,i]
    elif(i<=6):
        fraction_allframes_allvalences.append(data[:,i])

valences=np.arange(6)


fraction_allframes_allvalences=np.array(fraction_allframes_allvalences)
fraction_allvalences_allframes=fraction_allframes_allvalences.T

fraction_allvalences_steadystate=fraction_allvalences_allframes[-1]

print(valences)
print(fraction_allvalences_steadystate)

output_file=str(os.getcwd())+'/histograms_finalframe/lattice_restl'+str(restlength)+'_Nc'+str(Nclusters)+'_Np'+str(Np)+'_R'+str(R)+'_rB'+str(radiusB)+'_areafrac'+str(areafraction)+'_eps'+str(epsilon)+'_kspring'+str(kspring)+'_gammaA'+str(gammaA)+'_gammapatch'+str(gammapatch)+'.combineallseeds.allruns.histogramvalences.data'

y=np.array([valences,fraction_allvalences_steadystate],dtype=object)
y.dump(output_file)



