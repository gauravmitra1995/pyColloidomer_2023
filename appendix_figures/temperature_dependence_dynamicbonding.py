import math
import seaborn as sns
sns.set_context("poster", font_scale=1.0)
import numpy as np
import matplotlib.pyplot as plt

#Details about these parameters are described in Appendix under "Temperature dependence of binding/unbinding"

Tm=1.2  #Tmelt
alpha=200.0  #steepness parameter 
kon_init=100.0
kon_melt=0.0
koff_init1=5.0   #epsilon 3.0 
koff_init2=0.00001  #epsilon 16.0

Tlist=np.linspace(1.0,1.4,50)
f=np.zeros(Tlist.size)
f=(np.tanh(alpha*(Tlist-Tm))+1.0)/2.

koff_melt1=kon_init+kon_melt-koff_init1
koff_melt2=kon_init+kon_melt-koff_init2

kon=kon_init*(1-f)+f*kon_melt
koff1=koff_init1*(1-f)+f*koff_melt1
koff2=koff_init2*(1-f)+f*koff_melt2
frac_1=kon/(kon+koff1)
frac_2=kon/(kon+koff2)

fig, ax = plt.subplots(2,1, figsize=(6,7))
ax[0].plot(Tlist,f,color='indigo',linewidth=4.0)
ax[0].set_xlabel(r'$T/T^{*}$',fontsize=20,fontweight='bold')
ax[0].set_ylabel(r'$g(T)$',fontsize=20,fontweight='bold')
ax[0].axvline(x=Tm,linestyle='--',color='black')
ax[0].text(1.2,0.1,r"$\rightarrow T_{\mathrm{melt}}$",fontsize=20)
ax[0].set_title(r"$(a)$", fontsize=20,fontweight='bold')
ax[0].set_ylim(-0.1,1.1)
ax[0].set_yticks([0.0,0.5,1.0])

ax[1].plot(Tlist,frac_2,color='crimson',linewidth=4.0,linestyle='-',label=r'$\varepsilon = $'+str(np.round(np.log(kon_init/koff_init2),0)))
ax[1].plot(Tlist,frac_1,color='crimson',linewidth=4.0,linestyle=':',label=r'$\varepsilon = $'+str(np.round(np.log(kon_init/koff_init1),0)))
ax[1].axvline(x=Tm,linestyle='--',color='black')
ax[1].text(1.2,0.9,r"$\rightarrow T_{\mathrm{melt}}$",fontsize=20)
ax[1].set_xlabel(r'$T/T^{*}$',fontsize=20,fontweight='bold')
ax[1].set_ylabel(r'$f_{\mathrm{bound}}(T)$',fontsize=20,fontweight='bold')
ax[1].axvline(x=Tm,linestyle='--',color='black')
ax[1].set_title(r"$(b)$",fontsize=20,fontweight='bold')
ax[1].set_ylim(-0.1,1.1)
ax[1].set_yticks([0.0,0.5,1.0])
ax[1].legend(loc='lower left',prop={'size': 15})

fig.tight_layout()
plt.savefig('tempdependence.png',bbox_inches='tight')
plt.close()

