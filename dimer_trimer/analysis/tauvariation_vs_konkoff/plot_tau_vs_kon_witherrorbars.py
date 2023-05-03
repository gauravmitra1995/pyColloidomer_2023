from functions_cluster_analysis import *
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import sys
import ast

parser=argparse.ArgumentParser()

#Parameters REQUIRED to be passed as arguments from outside

parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int,default=2)
parser.add_argument("--R",type=float,default=50.0)
parser.add_argument("--Np",type=int,default=100)
parser.add_argument("--gammaA",type=float,default=0.1)
parser.add_argument("--dt",default=0.001,type=float)

args = parser.parse_args()
locals().update(vars(args))

if(epsilon!='infinite'):
    epsilon=float(epsilon)
else:
    epsilon=str(epsilon)
    
fig,ax=plt.subplots(figsize=(25,20),dpi=100)

colors=['blue','red','green']

r=0

file_pattern = os.path.join(os.getcwd(), "data_allseeds", "gammaA"+str(gammaA)+"_gammapatch*", "epsilon" + str(epsilon), "dimer_epsilon" + str(epsilon) + "_gammapatch*_R"+str(R)+".allsamples.recruitmenttimes.txt")
file_list = sorted(glob.glob(file_pattern),key=lambda x: (float(os.path.basename(x).split("_")[2].replace('gammapatch',''))))

for f in file_list:
    gammapatch=float(os.path.basename(f).split("_")[2].replace('gammapatch',''))
    konlist=[]
    meanlogtau_values=[]
    stddevlogtau_values=[]
    with open(f,"r") as file:
        for line in file:
            key, value = line.strip().split(": ")
            kon=ast.literal_eval(key)[0]
            konlist.append(kon)
            tau_values=np.array(ast.literal_eval(value))
            logtau=np.log10(tau_values.astype('float64'))
            meanoflogtau=np.mean(logtau)
            stddevoflogtau=np.std(logtau)
            meanlogtau_values.append(meanoflogtau)
            stddevlogtau_values.append(stddevoflogtau)
    logkon=np.log10(konlist).astype('float64')
    meanlogtau_values=np.array(meanlogtau_values)
    stddevlogtau_values=np.array(stddevlogtau_values)
    model = LinearRegression().fit(logkon.reshape(-1,1),meanlogtau_values.reshape(-1,1))
    slope=model.coef_
    intercept=model.intercept_
    label_str = r'$\gamma_{\mathrm{binder}} = $' + str(gammapatch)
    ax.errorbar(logkon,meanlogtau_values,yerr=stddevlogtau_values,marker='o',linestyle='None',markersize=20,color=colors[r],label=str(label_str),elinewidth=5.0,capsize=10,capthick=5.0)
    ax.plot(logkon,logkon*slope[0]+intercept[0],linestyle='--',marker='o',markersize=0,linewidth=8.0,color=colors[r],label='Fit for '+str(label_str))
    
    r+=1

plt.yticks(fontsize=50)
plt.xticks(fontsize=50)
plt.xlabel(r'$\log_{10} (k_{\mathrm{on}})$',labelpad=20,fontsize=70)
plt.ylabel(r'$\log_{10} (\tau_{\mathrm{rec}})$',labelpad=20,fontsize=70)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(5)
plt.title('Dimer: '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
ax.legend(loc='best', prop={'size':40})
fig.tight_layout()
plt.savefig('final_figures/dimer_tauvskon_eps'+str(epsilon)+'_R'+str(R)+'.png',bbox_inches='tight')
plt.close()

        

    


