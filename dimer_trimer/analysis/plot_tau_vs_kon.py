from functions_cluster_analysis import *
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import sys

parser=argparse.ArgumentParser()

#Parameters REQUIRED to be passed as arguments from outside

parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int,default=2)
parser.add_argument("--R",type=float,default=50.0)
parser.add_argument("--Np",type=int,default=100)
#parser.add_argument("--gammapatch",type=float)
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
gammapatch_list=[]
scatter = []
fit=[]
for filename in sorted(glob.glob(str(os.getcwd())+"/tauvskon/tauvskon_eps"+str(epsilon)+"_gammapatch*_R"+str(R)+".data"),key=lambda x:(float(os.path.basename(x).split("_")[2].replace('gammapatch','')))):
    gammapatch=float(os.path.basename(filename).split("_")[2].replace('gammapatch',''))
    gammapatch_list.append(gammapatch)
    data = (np.load(filename,allow_pickle=True))
    konlist=data[0]
    logkon=np.log10(np.array(konlist).astype('float64'))
    recruitment_time_list=data[1]
    logtaurec=np.log10(np.array(recruitment_time_list).astype('float64'))
    model = LinearRegression().fit(logkon.reshape(-1,1),logtaurec.reshape(-1,1))
    slope=model.coef_
    intercept=model.intercept_
    label_str = r'$\gamma_{\mathrm{binder}} = $' + str(gammapatch)
    scatter_line, =ax.plot(logkon,logtaurec,marker='o',linestyle='None',markersize=30,color=colors[r],label=str(label_str))
    fit_line, =ax.plot(logkon,logkon*slope[0]+intercept[0],linestyle='--',marker='o',markersize=0,linewidth=8.0,color=colors[r],label='Fit for '+str(label_str))
    scatter.append(scatter_line)
    fit.append(fit_line)
    r+=1

plt.yticks(fontsize=50)
plt.xticks(fontsize=50)
plt.xlabel(r'$\log_{10} (k_{\mathrm{on}})$',labelpad=20,fontsize=70)
plt.ylabel(r'$\log_{10} (\tau_{\mathrm{rec}})$',labelpad=20,fontsize=70)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(5)
plt.title('Dimer: '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
ax.legend(handles=scatter+fit,loc='best', prop={'size':40})
fig.tight_layout()
plt.savefig('final_figures/dimer_tauvskon_eps'+str(epsilon)+'_R'+str(R)+'.png',bbox_inches='tight')
plt.close()


#only for high epsilons we have a two time scale patch formation process 
if(epsilon>13):
    fig,ax=plt.subplots(figsize=(25,20),dpi=100)
    colors=['magenta','navy','darkgreen']
    r=0
    gammapatch_list=[]
    scatter = []
    fit=[]
    for filename in sorted(glob.glob(str(os.getcwd())+"/tausatvskon/tausatvskon_eps"+str(epsilon)+"_gammapatch*_R"+str(R)+".data"),key=lambda x:(float(os.path.basename(x).split("_")[2].replace('gammapatch','')))):
        gammapatch=float(os.path.basename(filename).split("_")[2].replace('gammapatch',''))
        gammapatch_list.append(gammapatch)
        data = (np.load(filename,allow_pickle=True))
        konlist=data[0]
        logkon=np.log10(np.array(konlist).astype('float64'))
        patchsaturation_time_list=data[1]
        logtausat=np.log10(np.array(patchsaturation_time_list).astype('float64'))
        model = LinearRegression().fit(logkon.reshape(-1,1),logtausat.reshape(-1,1))
        slope=model.coef_
        intercept=model.intercept_
        label_str = r'$\gamma_{\mathrm{binder}} = $' + str(gammapatch)
        scatter_line, =ax.plot(logkon,logtausat,marker='o',linestyle='None',markersize=30,color=colors[r],label=str(label_str))
        #fit_line, =ax.plot(logkon,logkon*slope[0]+intercept[0],linestyle='--',marker='o',markersize=0,linewidth=8.0,color=colors[r],label='Fit for '+str(label_str))
        scatter.append(scatter_line)
        #fit.append(fit_line)
        r+=1
    plt.yticks(fontsize=50)
    plt.xticks(fontsize=50)
    plt.xlabel(r'$\log_{10} (k_{\mathrm{on}})$',labelpad=20,fontsize=70)
    plt.ylabel(r'$\log_{10} (\tau_{\mathrm{sat}})$',labelpad=20,fontsize=70)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    plt.ylim(0,8)
    plt.title('Dimer: '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
    ax.legend(handles=scatter,loc='best', prop={'size':40})
    fig.tight_layout()
    plt.savefig('final_figures/dimer_tausatvskon_eps'+str(epsilon)+'_R'+str(R)+'.png',bbox_inches='tight')
    plt.close()

