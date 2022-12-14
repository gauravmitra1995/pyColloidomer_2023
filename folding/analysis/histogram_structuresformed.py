import os
import sys
import pickle
import numpy as np
from gsd import hoomd as gsd
import argparse
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import glob
sns.set_context("poster", font_scale=1.0)

def countX(lst, x):
    count = 0
    for ele in lst:
        if (ele == x):
            count = count + 1
    return count

def assign_structure(valences,filename):
    no_of_1s=countX(valences,1)
    no_of_2s=countX(valences,2)
    no_of_3s=countX(valences,3)
    no_of_4s=countX(valences,4)
    no_of_5s=countX(valences,5)
    no_of_6s=countX(valences,6)

    structure=""
    count_errors=0
    if(no_of_1s==0 and no_of_2s==2 and no_of_3s==2 and no_of_4s==3 and no_of_5s==0 and no_of_6s==0):
        structure="Ladder"
    elif(no_of_1s==0 and no_of_2s==3 and no_of_3s==1 and no_of_4s==2 and no_of_5s==1 and no_of_6s==0):
        structure="Rocket"
    elif(no_of_1s==0 and no_of_2s==2 and no_of_3s==3 and no_of_4s==1 and no_of_5s==1 and no_of_6s==0):
        structure="Chevron"
    elif(no_of_1s==0 and no_of_2s==0 and no_of_3s==6 and no_of_4s==0 and no_of_5s==0 and no_of_6s==1):
        structure="Flower"
    else:
        structure="Error" #for incompletely folded structures (not matching any of the assigned geometries)

    return structure
    


ensemble=[]
filenames=[]

for filename in sorted(glob.glob("valence_data/*.allruns.valences.data")):
    data=np.load(filename,allow_pickle=True)
    #At frame=100,300,500,700,900 folded structures are present. (since there are 5 folding/unfolding cycles, each cycle is of 200 frames)
    ensemble.append(list(data[100]))
    ensemble.append(list(data[300]))
    ensemble.append(list(data[500]))
    ensemble.append(list(data[700]))
    ensemble.append(list(data[900]))
    filenames.append(filename)
    filenames.append(filename)
    filenames.append(filename)
    filenames.append(filename)
    filenames.append(filename)

structure_ensemble=[]

counts=[]
structures=["Ladder","Chevron","Rocket","Flower"]

k=0
for e in ensemble:
    structure=assign_structure(e,filenames[k])
    structure_ensemble.append(structure)
    k+=1

count_ladder=countX(structure_ensemble,"Ladder")
count_chevron=countX(structure_ensemble,"Chevron")
count_rocket=countX(structure_ensemble,"Rocket")
count_flower=countX(structure_ensemble,"Flower")
count_error=countX(structure_ensemble,"Error")

counts=[count_ladder,count_chevron,count_rocket,count_flower]

fig,ax=plt.subplots(figsize=(20,15),dpi=100)
plt.bar(structures,np.array(counts)/len(structure_ensemble),width=0.5,color=['maroon','lightgreen','teal','deeppink'],edgecolor='black',linewidth=7.0)
plt.xlabel('Possible folded structures',fontsize=60,labelpad=20)
plt.ylabel('Fraction',fontsize=60,labelpad=20)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(5)
ax.tick_params(labelsize=50,axis='both', which='major', pad=20)
plt.xticks(fontsize=50)
plt.yticks(ticks=np.arange(0,0.55,0.05),fontsize=50)
fig.tight_layout()
plt.savefig('final_figures/histogram_foldedstructures.png',bbox_inches='tight')
plt.close()



