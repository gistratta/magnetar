
"""
    Author: Giulia Stratta

    v2: Jan 2024
    Purpose: compute magnetar parameters that best reproduce X-ray afterglow "plateau"
    following formalism described in Dall'Osso et al. 2011 and Stratta et al. 2018

    To run the script:
    >python loglumfit.py

    Input:
    GRB name, redshift, jet half-opening angle

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from datetime import datetime
import set_iniparam
import read_XRTdata
import model_ax_all
import fit_model
import define_bestfitparam
from astropy.cosmology import Planck18 as cosmo
import sys


print('')
print('######################################################################')
print('# 1. Reading and plotting data')
print('#')
print('# INPUT: GRB name, redshift, jet half-opening angle')
print('######################################################################')
print('')


plt.clf()
plt.close()
#plt.ion()

print('')
save = input('Do you want to write output file and save figure? [e.g. n]') or 'n'
#save='n'


# --- Input: 
print('')
grb = input('Type GRB name [e.g. 051221A]:') or '051221A'
z = float(input('Type GRB redshift [e.g. 0.5464]:') or '0.5464')
thetaj = float(input('Type jet half-opening angle in deg or 5.0 if not known [e.g. 5.5]:') or '5.5')
print('')
#---


# --- Defining the jet collimation factor and luminosity distance
collf=(1-np.cos(thetaj*np.pi/180))
DL = cosmo.luminosity_distance(z).to('cm').value

print('...reading initial parameters setting (rest frame luminosity band of the X-ray afterglow and magnetar parameter set (e.g. radius,collimation, etc.))')
print('')
Er1,Er2,thetam,alphax,k,omi,spini,B,a1_norm,Li_norm = set_iniparam.setiniparam()

# --- magnetar collimation factor
collfm=(1-np.cos(thetam))

# --- factor that accounts that only the fraction of magnetar jet within the afterglow jet contributes to energize the shock
fb = collf/collfm

#---
print('...reading GRB data')
tx,dtxp,dtxm,fx,dfxp,dfxm,FnuxmJy,dFnuxp_mJy,dFnuxm_mJy,gammaxrt,dgammaxrtp,dgammaxrtm=read_XRTdata.read_data(grb)

#--- Converting flux to luminosity
print('...Computing K correction')
beta=gammaxrt-1
Kcorr=(Er1**(1.-beta+1e-15) - Er2**(1.-beta+1e-15)) / ( (0.3*(1+z))**(1.-beta+1e-15)-(10.0*(1+z))**(1.-beta+1e-15))

# ---
print('...Computing luminosity in the ',Er1,Er2,' [keV] band')
ltot=(fx*Kcorr*collf)*(4*np.pi*(DL**2))
ltotp=((fx+dfxp)*Kcorr*collf)*(4*np.pi*(DL**2))
ltotm=((fx+dfxm)*Kcorr*collf)*(4*np.pi*(DL**2))
dlp = ltotp-ltot
dlm = ltot-ltotm
dl = (dlp+dlm)/2
llum=np.log10(ltot)
dllump=np.log10(ltotp)-llum
dllumm=llum - np.log10(ltotm)
dllum=(dllump+dllumm)/2          # simmetric errors are needed in the fit function

# Compute rest frame time
ttot=tx/(1+z)
ttotp = (tx+dtxp)/(1+z)
ttotm = (tx+dtxm)/(1+z)
ltime=np.log10(ttot)
dltimep=np.log10(ttotp)-np.log10(ttot)
dltimem=np.log10(ttot)-np.log10(ttotm)

# ---
print('...Plotting afterglow luminosity light curve')
plt.title('GRB'+grb+' '+str(Er1)+'-'+str(Er2)+' keV')
plt.plot(ltime,llum,'.', label='Swift/XRT data',color='b')
plt.errorbar(ltime, llum, yerr=[dllumm,dllump],fmt='none',ecolor='b')
plt.xlabel('log time from trigger [s]')
plt.ylabel('log Luminosity corr. for beaming [erg s^-1]')

plt.savefig('lc_data.png')
#plt.show()



print('######################################################################')
print('# 2. Defining the lightcurve properties and plot first guess model')
print('#')
print('# INPUT: lightcurve morphology keyword, start and end time of model and of plateau, plateau luminosity, steep decay start time, luminosity and decay index')
print('######################################################################')


print('')
print('---------- Lightcurve Morphology Keyword ---')
print(' SP  = steep decay + plateau (+ afterglow)')
print(' PA  = (afterglow +) plateau + afterglow (i.e. no steep decay ')
print('--------------------------------------------')
print('')

print('')
print('--> !!! Please check the plot saved in lc_data.png !!! <--')
print('')

morphology = input('Type morphology keyword  [default: PA]') or 'PA'
print('')


# --- Define empirical functions

def func_brokenpl(x, c0, x0, alpha1, alpha2, omega=3):
    return (c0/(2**(-1/omega)))*((x/x0)**(alpha1*omega) + (x/x0)**(alpha2*omega))**(-(1./omega))

def func_powerlaw(x, x0, m, c):
    return c*((x/x0)**m)


# ---  Plateau duration and luminosity

lT0_aft = float( input('Plateau start time [s] (log scale) [e.g. 3.2]') or '3.2')
T0_aft = 10**(lT0_aft)

ltau0 = float(input('Plateau end time [s] (log scale) [e.g. 4.0]')  or '4.0'  )      # plateau duration in s
tau0 = 10**(ltau0)

# --- Model start and end time (starting early is better to define the plateau morphology)

# Only Data contained within this temporal window will be used in the fit
print('')
lT0 = float(input('Model start time [s] (log scale) [e.g. 1.5]') or '1.5' )
T0 = 10**lT0
lTend = float(input('Model end time [s] (log scale) [e.g. 5.5]') or '5.5' )
Tend = 10**lTend
print('')


#############################################################################

# Fitting the plateau with a broken power law and finding the plateau end time tau and luminosity at tau 
# and use best fit parameters to build first guess magnetar model

# --- fitting the plateau and post plateau with a smoothed borken power law
ttotnew = ttot[np.where(ttot>T0_aft)]
ltotnew = ltot[np.where(ttot>T0_aft)]
popt, pcov = curve_fit(func_brokenpl,ttotnew,ltotnew,sigma = dl[np.where(ttot>T0_aft)],p0=(ltotnew[0],tau0,0.0,1.5),bounds =[(ltotnew[0]/10.0,tau0/10.0,-0.1,1.01),(ltotnew[0]*10,tau0*10,0.75,2.5)], maxfev=10000)
L_tau = popt[0]                         # plateau luminosity at tau
tau = popt[1]                           # plateau duration
alpha2 = popt[3]                        # post plateau decay index
k = alpha2-1                            # radiative efficiency factor k+1 = alpha2


# ---
print('...Building first guess model')

P = np.sqrt(k*3*1e52/((L_tau/collf)*tau))  # NS spin period in ms
B = (np.sqrt(7/tau)) * 100 * P             # NS magnetic field strenght in 10^14 Gauss
omi = 2*np.pi/P                            # NS spin frequency in 1e3 Hz
alphax = 0.1                               # alphax = (3-n)/2 where n is the NS braking index (see Stratta et al. 2018)

#############################################################################
if morphology == 'PA':


    alpha1 = popt[2]                            # plateau decay index
    L0_aft = L_tau*((T0_aft/tau)**(-alpha1))    # plateau initial lum

    # --- Magentar model for PA morphology
    def model_ax(logt,k,B,omi):
        return model_ax_all.model_ax_all_afterglow(logt,k,B,omi,t00=T0,Ka1=a1_norm,KLi=Li_norm,collfact=fb,alphax=alphax,L0_aft=L0_aft,t0_aft=T0_aft)

    # --- Temporal grid for the model plot (only for plot purposes)
    if lT0 < ltime[0]:
        # if model starts before data:
        t=np.logspace(lT0,ltime[-1]+0.5, num=1000, base=10.0)
    else:
        t=np.logspace(ltime[0]-0.3,ltime[-1]+0.5, num=1000, base=10.0)

    # --- First guess magnetar model
    logmodel=model_ax(np.log10(t),k,B,omi)
 

#############################################################################
elif morphology == 'SP':

    L0_aft = ltotnew[0]                     # plateau luminosity at plateau start time  

    # --- computing the steep decay index

    ttotnews = ttot[np.where((ttot<T0_aft) & (ttot>T0) )]
    ltotnews = ltot[np.where((ttot<T0_aft) & (ttot>T0) )]   

    # find midpoint of steep decay:
    ns = len(ttotnews)
    if ns % 2 == 0:
        T0_prompt = ttotnews[int(ns/2.0)]
        L0_prompt = ltotnews[int(ns/2.0)]        
    else:
        nn=ns+1
        T0_prompt = ttotnews[int(nn/2.0)]
        L0_prompt = ltotnews[int(nn/2.0)]    

    popt_s,pcov_s = curve_fit(func_powerlaw,ttotnews,ltotnews,sigma=dl[np.where((ttot<T0_aft) & (ttot>T0))], p0=(T0_prompt,-2.5,L0_prompt),bounds=[(0.5,-10.0,1e40),(1000.0,-1.5,1e52)])
    delta = -popt_s[1]
    
    
    # --- Magentar model for SP morphology
    def model_ax(logt,B,omi,delta):
        return model_ax_all.model_ax_all_steepdecay(logt,B,omi,delta,t00=T0_prompt,Ka1=a1_norm,KLi=Li_norm,collfact=fb,alphax=alphax,L0=L0_prompt,k=k)

    # --- Temporal grid for the model plot (only for plot purposes)
    if lT0 < ltime[0]:
        # if model starts before data:
        t=np.logspace(lT0,ltime[-1]+0.5, num=1000, base=10.0)
    else:
        t=np.logspace(ltime[0]-0.3,ltime[-1]+0.5, num=1000, base=10.0)

    # --- First guess magnetar model
    #logmodel=model_ax(np.log10(t),B,omi,delta)
    logmodel=model_ax(np.log10(t),B,omi,delta)

#############################################################################

# --- Magnetar only
logmodel_plateau_only = model_ax_all.model_plateau(np.log10(t+0.1),k,B,omi,t00=T0,Ka1=a1_norm,KLi=Li_norm,L0=L0_aft,collfact=fb,alphax=alphax)

# ---
print('...plotting the first guess model')
plt.plot(np.log10(t),logmodel,'r',alpha=0.3,label=('First guess model: '+morphology+', alphax='+str(alphax)+'; k='+str(round(k,2))+'; B14='+str(round(B,0))+'; Pms='+str(round(P,1))))
plt.plot(np.log10(t),logmodel_plateau_only,'--r',alpha=0.3,label=('First guess magnetar energy injection only'))

# --- plot the plateau starttime 
yminv=llum.min()-1.0
ymaxv=llum.max()+1.0
plt.vlines(lT0_aft,yminv,ymaxv,colors='b', linestyles='dotted',label='plateau start time')

# --- plot the data used in the fit
if lT0 < ltime[0]:
    plt.vlines(ltime[0],yminv,ymaxv,linestyles='dashed',colors='black',label='used data')
else:
    plt.vlines(lT0,yminv,ymaxv,linestyles='dashed',colors='black',label='used data')

if lTend < ltime[-1]:
    plt.vlines(lTend,yminv,ymaxv,linestyles='dashed',colors='black')
else:
    plt.vlines(ltime[-1],yminv,ymaxv,linestyles='dashed',colors='black')





print('')
print('')
print('#############################################################################')
print('# 3. Fitting the model')
print('#############################################################################')
print('')


# -- defining data to be used for the fit
lt0=ltime[np.where(ltime >= np.log10(T0))]
ll0=llum[np.where(ltime >= np.log10(T0))]
dll0=dllum[np.where(ltime >= np.log10(T0))]
dll0p=dllump[np.where(ltime >= np.log10(T0))]
dll0m=dllumm[np.where(ltime >= np.log10(T0))]

lt=lt0[np.where(lt0 <= np.log10(Tend))]
ll=ll0[np.where(lt0 <= np.log10(Tend))]
dll=dll0[np.where(lt0 <= np.log10(Tend))]
dllp=dll0p[np.where(lt0 <= np.log10(Tend))]
dllm=dll0m[np.where(lt0 <= np.log10(Tend))]

#---
print('...fitting the model ')


######################################################################
if morphology == 'PA':

    # --- Fitting the model 
    popt,pcov=fit_model.fitmodel(model_ax,lt,ll,dll,k,B,omi)
    
    # --- Computing best fit parameters
    k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,Pms,dPms,dof,mychi,p_value = define_bestfitparam.define(popt,pcov,lt,ll,dll,t0=T0,alphax=alphax,k=k,B=B,omi=omi,model=model_ax)

    # --- Best fit model
    logmodel_bf=model_ax(np.log10(t),k_bf,B_bf,omi_bf)

    # --- Plotting best fit model
    plt.plot(np.log10(t),logmodel_bf,color='g',linestyle='-',label=('Best fit with alpha='+str(alphax)+': B14='+str(round(B_bf,0))+'+/-'+str(round(dB_bf,0))+'; Pms='+str(round(Pms,1))+'+/-'+str(round(dPms,1))+'; k='+str(round(k,2)) )) 
 
    # --- Magnetar energy injection only component
    logmodel_plateau_only_bf = model_ax_all.model_plateau(np.log10(t+0.1),k,B_bf,omi_bf,t00=T0,Ka1=a1_norm,KLi=Li_norm,L0=L0_aft,collfact=fb,alphax=alphax)
 
    # --- Plotting magnetar energy injection only component
    #plt.plot(np.log10(t),logmodel_plateau_only_bf,color='g',linestyle=':',label=('Best fit magnetar energy injection only'))


###########################################################################

if morphology == 'SP':

    print('...fixing k to ', round(k,2))

    # --- Fitting the model 
    popt_sp,pcov_sp=fit_model.fitmodel_no_k(model_ax,lt,ll,dll,B,omi,delta)

    # --- Computing best fit parameters
    B_bf,dB_bf,omi_bf,domi_bf,delta_bf,ddelta_bf,Pms,dPms,dof,mychi,p_value = define_bestfitparam.define_no_k(popt_sp,pcov_sp,lt,ll,dll,t0=T0_prompt,alphax=alphax,B=B,omi=omi,delta=delta,model=model_ax)

    # --- Best fit model
    logmodel_bf=model_ax(np.log10(t),B_bf,omi_bf,delta_bf)

    # --- Plotting best fit model
    plt.plot(np.log10(t),logmodel_bf,color='g',linestyle='-',label=('Best fit with alpha='+str(alphax)+': B14='+str(round(B_bf,0))+'+/-'+str(round(dB_bf,0))+'; Pms='+str(round(Pms,1))+'+/-'+str(round(dPms,1))+'; k='+str(round(k,2)) )) 
 
    # --- Magnetar energy injection only component
    logmodel_plateau_only_bf = model_ax_all.model_plateau(np.log10(t+0.1),k,B_bf,omi_bf,t00=T0,Ka1=a1_norm,KLi=Li_norm,L0=L0_aft,collfact=fb,alphax=alphax)
 
    # --- Plotting magnetar energy injection only component
    #plt.plot(np.log10(t),logmodel_plateau_only_bf,color='g',linestyle=':',label=('magnetar energy injection'))



elif morphology == 'SPA':
    print('...fixing k to ',k)
    popt,pcov=fit_model.fitmodel_no_k(model_ax,lt,ll,dll,B,omi,delta)
    B_bf,dB_bf,omi_bf,domi_bf,delta_bf,ddelta_bf,Pms,dPms,dof,mychi,p_value = define_bestfitparam.define_no_k(popt,pcov,lt,ll,dll,t0=T0,alphax=alphax,B=B,omi=omi,delta=delta,model=model_ax)
    logmodel_bf=model_ax(np.log10(t),B_bf,omi_bf,delta_bf)
    plt.plot(np.log10(t),logmodel_bf,'r--',label=('delta='+str(round(delta_bf,1))+'+/-'+str(round(ddelta_bf,1))+' B14='+str(round(B_bf,0))+'+/-'+str(round(dB_bf,0))+' Pms='+str(round(Pms,1))+'+/-'+str(round(dPms,1))+' alpha='+str(alphax)+' k(fix)='+str(k)))

    if save == 'y':
        print(' ... defining the name of the output log file:')
        thetamstr=str(round(thetam*180/np.pi,0))
        output='./output/'+grb+'_thetam_'+thetamstr+'deg_'+str(Er1)+'_'+str(Er2)+'keV'+'_'+morphology+'.log'
        print('output file:', output)

        print(' ... saving best fit parameters in the output log file:')
        if not os.path.isfile(output):
            os.system('touch '+output)
            out_file = open(output,'a')
            out_file.write("grb,alphax,T0,k,B,dB,Pms,dPms,thetaj,thetam,mychi,dof,p-value,time"+"\n")
            out_file.close()

    if save == 'y':
        out_file = open(output,"a")
        out_file.write(grb+","+str(alphax)+","+str(T0)+","+str("%.5f" %k)+","+str("%.5f" %B_bf)+","+str("%.5f" %dB_bf)+","+str("%.5f" %Pms)+","+str("%.5f" %dPms)+","+str("%.5f" %thetaj)+","+str("%.5f" %thetam)+","+str("%.5f" %mychi)+","+str("%.5f" %dof)+","+str("%.5f" %p_value)+","+str(datetime.now())+"\n")
        out_file.close()


plt.title('GRB'+grb+' ['+str(Er1)+'-'+str(Er2)+' keV]')
#plt.plot(ltime,llum,'.', label='Swift/XRT data',color='b')
plt.errorbar(ltime, llum, yerr=[dllumm,dllump],xerr=[dltimem,dltimep],fmt='none',ecolor='b')
#print('...Printing the fitting temporal window as vertical lines')
#plt.axvline(x=np.log10(T0),color='k', linestyle='--')
#plt.axvline(x=np.log10(Tend),color='k', linestyle='--')

#plt.vlines(x=[np.log10(T0),lTend],ymin=42,ymax=45, colors='k',linestyles='dashed')

#plt.ylim(40,45)


plt.legend(fontsize=6,loc='lower left')
plt.xlabel('log time from trigger [s]')
plt.ylabel('log Luminosity corr. for beaming [erg s^-1]')
plt.ylim(llum[-1]-1,llum[0]+1)

if lT0 < ltime[0]:
    plt.xlim(lT0,ltime[-1]+0.5)
if lT0 > ltime[0]:
    plt.xlim(ltime[0]-0.3,ltime[-1]+0.5)





#------------------------------------------------

print('')
print('--------')
print('')
print('GRB',grb,' z=',z, 'theta_jet=', thetaj)
print('')

print('If you are not happy with the fit, here some tips:')
print('')
print('1) try a different morphology keyword')
print('2) try to change model start and end time (i.e. excluding some early and/or late data point)')
print('')
print('(to exit, close the plot window)')


#------------------------------------------------

def saving(save): 
    if save == 'y':
        
        print(' ... defining the name of the output log file:')
        thetamstr=str(round(thetam*180/np.pi,0))
        output='./output/'+grb+'_thetam_'+thetamstr+'deg_'+str(Er1)+'_'+str(Er2)+'keV_'+morphology+'.log'
        
        print('')
        print(' ... saving best fit parameters in the output log file: ',output)

        if morphology == 'PA':

            # --- writing the header
            if not os.path.isfile(output):
                os.system('touch '+output)
                out_file = open(output,'a')
                out_file.write("grb,alphax,k,B14,dB14,Pms,dPms,thetaj,thetam,mychi,dof,p-value,time"+"\n")
                out_file.close()

            out_file = open(output,"a")
            out_file.write(grb+","+str(alphax)+","+str("%.3f" %k_bf)+","+str("%.3f" %dk_bf)+","+str("%.1f" %B_bf)+","+str("%.1f" %dB_bf)+","+str("%.1f" %Pms)+","+str("%.1f" %dPms)+","+str("%.3f" %thetaj)+","+str("%.3f" %thetam)+","+str("%.1f" %mychi)+","+str(dof)+","+str("%.5f" %p_value)+","+str(datetime.now())+"\n")
            out_file.close()


        if morphology == 'SP':
            
            # --- writing the header
            if not os.path.isfile(output):
                os.system('touch '+output)
                out_file = open(output,'a')
                out_file.write("grb,alphax,delta,B14,dB14,Pms,dPms,thetaj,thetam,mychi,dof,p-value,time"+"\n")
                out_file.close() 

            out_file = open(output,"a")
            out_file.write(grb+","+str(alphax)+","+str("%.5f" %delta_bf)+","+str("%.5f" %ddelta_bf)+","+str("%.5f" %B_bf)+","+str("%.5f" %dB_bf)+","+str("%.5f" %Pms)+","+str("%.5f" %dPms)+","+str("%.5f" %thetaj)+","+str("%.5f" %thetam)+","+str("%.5f" %mychi)+","+str("%.5f" %dof)+","+str("%.5f" %p_value)+","+str(datetime.now())+"\n")
            out_file.close()




    if save == 'y':

        print('...saving the final plot (to exit, close the plot window)')
        plt.savefig('./output/'+grb+'_thetam_'+thetamstr+'deg_'+str(Er1)+'_'+str(Er2)+'keV_'+morphology+'.png')

    else: 
        plt.ioff()
        plt.show()
        print('')
        print('Exit without saving')


saving(save)

plt.ioff()
plt.show()

