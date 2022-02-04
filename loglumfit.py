
"""
    Author: Giulia Stratta
This script

 * Defines model: L(t)=kE(t)/t (from Simone Notebook)
 * Fits model to the data and compute best fit param and cov. matrix


NOTE: if theta_jet is known, a collimation factor "collf" is applied to the luminosities

To run the script:
    ipython --matplotlib
    run loglumfit.py

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
from datetime import datetime
import set_iniparam
import read_XRTdata
import model_ax_all
import fit_model
import define_bestfitparam
import find_DL



print('')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('')
print('   Code: loglumfit_individual_alpha_fixed.py')
print('   Author: G. Stratta')
print('   version: 10/11/2021')

print('   Purpose: compute magnetar parameters that best reproduce X-ray afterglow "plateau" ')
print('')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('')



plt.clf()
plt.close()


grb=input('Type GRB name [e.g. 130603B]:') or '130603B'

print('')
print('         *** GRB',grb,' ***')
print('')

thetaj=float(input('Enter the GRB jet opening angle in deg [e.g. 4.5] or 90.0 if not known: ') or '4.5')
collf=(1-np.cos(thetaj*np.pi/180))

z = float(input('Enter the GRB redshift [e.g. 0.3564]:') or '0.3564')
DL = find_DL.findDL(z)

print('')
print('...reading GRB data')
tx,dtxp,dtxm,fx,dfxp,dfxm,gammaxrt,dgammaxrtp,dgammaxrtm=read_XRTdata.read_data(grb)


print('...reading initial parameters (rest frame luminosity band of the X-ray afterglow and magnetar parameter set (e.g. radius,collimation, etc.))')
Er1,Er2,thetam,alphax,k,omi,spini,B,a1_norm,Li_norm = set_iniparam.setiniparam()


print('...Computing K correction')
beta=gammaxrt-1
Kcorr=(Er1**(1.-beta+1e-15) - Er2**(1.-beta+1e-15)) / ( (0.3*(1+z))**(1.-beta+1e-15)-(10.0*(1+z))**(1.-beta+1e-15))


print('...Computing luminosity in the ',Er1,Er2,' [keV] band')
ltot=(fx*Kcorr*collf)*(4*np.pi*(DL**2))
ltotp=((fx+dfxp)*Kcorr*collf)*(4*np.pi*(DL**2))
ltotm=((fx+dfxm)*Kcorr*collf)*(4*np.pi*(DL**2))
llum=np.log10(ltot)
dllump=np.log10(ltotp)-llum
dllumm=llum - np.log10(ltotm)
# simmetric errors are needed in the fit function
dllum=(dllump+dllumm)/2
# Compute rest frame time
ttot=tx/(1+z)
ltime=np.log10(ttot)


print('...Plotting afterglow luminosity light curve')
plt.title('GRB'+grb+' ['+str(Er1)+'-'+str(Er2)+' keV]')
plt.plot(ltime,llum,'.', label='Swift/XRT data',color='b')
plt.errorbar(ltime, llum, yerr=[dllumm,dllump],fmt='none',ecolor='b')
plt.xlabel('log time from trigger [s]')
plt.ylabel('log Luminosity corr. for beaming [erg s^-1]')
plt.ion()
plt.show()


print('...Defining the fit temporal window')
print('')
lstartTxrt=input('Type plateau start time (in logarithmic scale -see plot-)[e.g. 1.7]:') or '1.7'
lTend=input('Type plateau end time (in logarithmic scale -see plot-)[e.g. 4.0]:') or '4.0'
startTxrt=10**(float(lstartTxrt))
Tend=10**(float(lTend))
print('   Fitting temporal window [s]: '+str(startTxrt)+'-'+str(Tend))

print('')
print('...Printing the fitting temporal window as vertical lines')
plt.axvline(x=np.log10(startTxrt),color='k', linestyle='--')
plt.axvline(x=np.log10(Tend),color='k', linestyle='--')
plt.show()


print('...Defining the data set to be fitted (no steep decay, no post jet break)')
lt0=ltime[np.where(ltime >= np.log10(startTxrt))]
ll0=llum[np.where(ltime >= np.log10(startTxrt))]
dll0=dllum[np.where(ltime >= np.log10(startTxrt))]
dll0p=dllump[np.where(ltime >= np.log10(startTxrt))]
dll0m=dllumm[np.where(ltime >= np.log10(startTxrt))]
lt=lt0[np.where(lt0 <= np.log10(Tend))]
ll=ll0[np.where(lt0 <= np.log10(Tend))]
dll=dll0[np.where(lt0 <= np.log10(Tend))]
dllp=dll0p[np.where(lt0 <= np.log10(Tend))]
dllm=dll0m[np.where(lt0 <= np.log10(Tend))]



print('...Computing luminosity at the plateau start time ')
for i in range(0,len(ttot)-1):
    #if ttot[i] >= startplateau:
    if ttot[i] >= startTxrt:
        #L051=(10**(llum[i])/1e51)
        L051=(10**(llum[i]))
        break
print("        --> Plateau 0.3-30 keV luminosity at start (10^52 erg/s):"+str("%.2e" % L051))


# define the magnetar collimation factor and the factor that takes into account
# the fact that only the portion within the afterglow jet goes to energizes the shock
collfm=(1-np.cos(thetam))
fb = collf/collfm



print('...defining the magnetar model (alpha=0.9)')
def model_ax(logt,k,B,omi):
    return model_ax_all.model_ax_all(logt,k,B,omi,t00=startTxrt,Ka1=a1_norm,KLi=Li_norm,L0=L051,collfact=fb,alphax=alphax)
print('...plotting the initial magnetar model')
t0=np.logspace(-1.,7., num=100, base=10.0)
t1=t0[np.where(t0>startTxrt)]
t=t1[np.where(t1<10**(ltime[-1]))]
logmodel=model_ax(np.log10(t),k,B,omi)


print('...fitting the magnetar model')
popt,pcov=fit_model.fitmodel(model_ax,lt,ll,dll,k,B,omi)
k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,Pms,dPms,E0,dof,mychi,p_value = define_bestfitparam.define(popt,pcov,lt,ll,dll,L0=L051,t0=startTxrt,alphax=alphax,k=k,B=B,omi=omi,model=model_ax)


print('---')
print ('')

save=input('Do you want to write output file and save figure? [e.g. n]') or 'n'

print('')
print('...plotting the best fit magnetar model')
logmodel_bf=model_ax(np.log10(t),k_bf,B_bf,omi_bf)
plt.plot(np.log10(t),logmodel_bf,'r--',label=('alpha='+str(alphax)+';B14='+str(round(B_bf,0))+'+/-'+str(round(dB_bf,0))+';Pms='+str(round(Pms,1))+'+/-'+str(round(dPms,1))+';k='+str(round(k_bf,2))+'+/-'+str(round(dk_bf,2))+';$\chi^2$='+str(round(mychi,2))+';$N_{dof}$='+str(dof)))


print('')
print(' ... defining the name of the output log file:')
thetamstr=str(round(thetam*180/np.pi,0))
output='./output/'+grb+'_'+thetamstr+'deg_'+str(Er1)+'_'+str(Er2)+'keV.log'
print('output file:', output)

if save == 'y':
    print(' ... saving best fit parameters in the output log file:')
    if not os.path.isfile(output):
        os.system('touch '+output)
        out_file = open(output,'a')
        out_file.write("grb,alphax,startTxrt,E0,k,dk,B,dB,Pms,dPms,thetaj,thetam,mychi,dof,p-value,time"+"\n")
        out_file.close()

if save == 'y':
    out_file = open(output,"a")
    out_file.write(grb+","+str(alphax)+","+str(startTxrt)+","+str(E0)+","+str("%.5f" %k_bf)+","+str("%.5f" %dk_bf)+","+str("%.5f" %B_bf)+","+str("%.5f" %dB_bf)+","+str("%.5f" %Pms)+","+str("%.5f" %dPms)+","+str("%.5f" %thetaj)+","+str("%.5f" %thetam)+","+str("%.5f" %mychi)+","+str("%.5f" %dof)+","+str("%.5f" %p_value)+","+str(datetime.now())+"\n")
    out_file.close()


print('...fitting the magnetar model with alpha=0.5')
alphax=0.5
popt,pcov=fit_model.fitmodel(model_ax,lt,ll,dll,k,B,omi)
k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,Pms,dPms,E0,dof,mychi,p_value = define_bestfitparam.define(popt,pcov,lt,ll,dll,L0=L051,t0=startTxrt,alphax=alphax,k=k,B=B,omi=omi,model=model_ax)
print('---')
print ('')
print('...plotting the magnetar model with alpha=0.5')
logmodel_bf_a05=model_ax(np.log10(t),k_bf,B_bf,omi_bf)
plt.plot(np.log10(t),logmodel_bf_a05,'c--',label=('alpha='+str(alphax)+';B14='+str(round(B_bf,0))+'+/-'+str(round(dB_bf,0))+';Pms='+str(round(Pms,1))+'+/-'+str(round(dPms,1))+';k='+str(round(k_bf,2))+'+/-'+str(round(dk_bf,2))+';$\chi^2$='+str(round(mychi,2))+';$N_{dof}$='+str(dof)))

if save == 'y':
    print(' ... saving best fit parameters in the same output log file:')
    out_file = open(output,"a")
    out_file.write(grb+","+str(alphax)+","+str(startTxrt)+","+str(E0)+","+str("%.5f" %k_bf)+","+str("%.5f" %dk_bf)+","+str("%.5f" %B_bf)+","+str("%.5f" %dB_bf)+","+str("%.5f" %Pms)+","+str("%.5f" %dPms)+","+str("%.5f" %thetaj)+","+str("%.5f" %thetam)+","+str("%.5f" %mychi)+","+str("%.5f" %dof)+","+str("%.5f" %p_value)+","+str(datetime.now())+"\n")
    out_file.close()

print('...fitting the magnetar model with alpha=0.1')
alphax=0.1
popt,pcov=fit_model.fitmodel(model_ax,lt,ll,dll,k,B,omi)
k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,Pms,dPms,E0,dof,mychi,p_value = define_bestfitparam.define(popt,pcov,lt,ll,dll,L0=L051,t0=startTxrt,alphax=alphax,k=k,B=B,omi=omi,model=model_ax)
print('---')
print ('')
print('...plotting the magnetar model with alpha=0.1')
logmodel_bf_a01=model_ax(np.log10(t),k_bf,B_bf,omi_bf)
plt.plot(np.log10(t),logmodel_bf_a01,'g--',label=('alpha='+str(alphax)+';B14='+str(round(B_bf,0))+'+/-'+str(round(dB_bf,0))+';Pms='+str(round(Pms,1))+'+/-'+str(round(dPms,1))+';k='+str(round(k_bf,2))+'+/-'+str(round(dk_bf,2))+';$\chi^2$='+str(round(mychi,2))+';$N_{dof}$='+str(dof)))

if save == 'y':
    print(' ... saving best fit parameters in the same output log file:')
    out_file = open(output,"a")
    out_file.write(grb+","+str(alphax)+","+str(startTxrt)+","+str(E0)+","+str("%.5f" %k_bf)+","+str("%.5f" %dk_bf)+","+str("%.5f" %B_bf)+","+str("%.5f" %dB_bf)+","+str("%.5f" %Pms)+","+str("%.5f" %dPms)+","+str("%.5f" %thetaj)+","+str("%.5f" %thetam)+","+str("%.5f" %mychi)+","+str("%.5f" %dof)+","+str("%.5f" %p_value)+","+str(datetime.now())+"\n")
    out_file.close()

plt.legend(fontsize=8)
plt.show()


if save == 'y':
    print('...saving the final plot')
    plt.savefig('./output/'+grb+'_'+thetamstr+'deg_'+str(Er1)+'_'+str(Er2)+'keV.png')

print('')
input('Done! Press enter to plot best fit model and exit')
