import scipy.special as sc
from sympy.functions import hyper
import numpy as np


"""
    Magnetar model from Dall'Osso et al. 2011 and Stratta et al. 2018 
    adapted to fit different lightcurve morphologies

"""

def model_ax_all_afterglow(logt,k,B,omi,t00,Ka1,KLi,collfact,alphax,L0_aft,t0_aft):
    """
    Description: 'AP' morphology
        This module should substitute the old "model_ax_all" that where t00 coinceded with starting data time
        while now it is free to be set at any time 

    Parameters:
        k = coefficient that depends on energy electron fraction and shock dynamical evolution (see Dall'Osso et al.2011)
        B = magnetic field (in units of 10^14 Gauss)
        omi = initial spin frequency (in units of 10^3 Hz)
        t00 = start time (--> should be put very early)
        Ka1 = normalization factor of initial spindown timescale a1 (Ka1=2Ic^3/R^6) where a1=Ka1B^2omi^2
        KLi = normalization factor of initial spindown luminosity Li (KL1=R^6/(6c^3)) where Li=KLiB^2omi^4
        t0_aft = time at which the afterglow is dominated by the plateau in s
        L0_aft = afterglow X-ray luminosity at t0_aft in units of erg/s
        collfact = collimation factor ratio of afterglow over magentar (1-cos(theta_jet))/(1-cos(theta_magnetar))
        alphax = (3-n)/2 where n is the breaking index of the magnetar (see Stratta et al. 2018)
    Output:
        luminosity in erg/s
    """
    t = 10**logt
    a1 = Ka1/((B**2.)*(omi**2.))
    Li = KLi*(B**2.)*(omi**4.)
    E0 = L0_aft*( (t00/t0_aft)**(-(1+k)) )*(t00/k)
    hg1_a = collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t00/a1)
    hg2_a = collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t/a1)
    f_ax = (k/t)*(1./(1. + k))*(t**(-k))*( E0*(t00**k)*(1.+k) + Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    return np.log10(f_ax)


def model_ax_all_steepdecay(logt,B,omi,delta,t00,Ka1,KLi,collfact,alphax,L0,k):
    """
    Description: 'SP' morphology
        Magnetar model without the afterglow component, plus a power law to model the steep decay

    Parameters:
        k = coefficient that depends on energy electron fraction and shock dynamical evolution (see Dall'Osso et al.2011)
        B = magnetic field (in units of 10^14 Gauss)
        omi = initial spin frequency (in units of 10^3 Hz)
        t00 = start time of ** steep decay ** 
        Ka1 = normalization factor of initial spindown timescale a1 (Ka1=2Ic^3/R^6) where a1=Ka1B^2omi^2
        KLi = normalization factor of initial spindown luminosity Li (KL1=R^6/(6c^3)) where Li=KLiB^2omi^4
        L0 = steep decay X-ray luminosity at t00 in units of erg/s
        collfact = collimation factor ratio of afterglow over magentar (1-cos(theta_jet))/(1-cos(theta_magnetar))
        alphax = (3-n)/2 where n is the breaking index of the magnetar (see Stratta et al. 2018)
        delta = steep decay power law index (module)
    Output:
        luminosity in erg/s
    """
    t = 10**logt
    a1 = Ka1/((B**2.)*(omi**2.))
    Li = KLi*(B**2.)*(omi**4.)
    hg1_a = collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t00/a1)
    hg2_a = collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t/a1)
    f_ax = (k/t)*(1./(1. + k))*(t**(-k))*( Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    f_ax_sd = f_ax + L0*(t/t00)**(-delta)
    return np.log10(f_ax_sd)


def model_ax_all_afterglow_steepdecay(logt,k,B,omi,delta,t00,Ka1,KLi,L0,collfact,alphax,L0_aft,t0_aft):
    """
    Description: 'SAP'
        Magnetar model plus a power law to model the steep decay
        NOTE: To be used when we have plateau preceded by an afterglow and before by a steep decay

   Parameters:
        k = coefficient that depends on energy electron fraction and shock dynamical evolution (see Dall'Osso et al.2011)
        B = magnetic field (in units of 10^14 Gauss)
        omi = initial spin frequency (in units of 10^3 Hz)
        t00 = model start time 
        t0_aft = time at which afterglow luminosity is defined (it can be the very start of the plateau) 
        L0_aft = afterglow luminosity at t0_aft
        Ka1 = normalization factor of initial spindown timescale a1 (Ka1=2Ic^3/R^6) where a1=Ka1B^2omi^2
        KLi = normalization factor of initial spindown luminosity Li (KL1=R^6/(6c^3)) where Li=KLiB^2omi^4
        L0 = X-ray luminosity at t00 in units of erg/s
        collfact = collimation factor ratio of afterglow over magentar (1-cos(theta_jet))/(1-cos(theta_magnetar))
        alphax = (3-n)/2 where n is the breaking index of the magnetar (see Stratta et al. 2018)
        delta = steep decay power law index (module)
    Output:
        luminosity in erg/s
    """
    t=10**logt
    a1=Ka1/((B**2.)*(omi**2.))
    Li=KLi*(B**2.)*(omi**4.)
    #E0=L0*t00/k
    E0 = L0_aft*( (t00/t0_aft)**(-(1+k)) )*(t00/k)
    hg1_a=collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t00/a1)
    hg2_a=collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t/a1)
    #f_ax=(k/t)*(1./(1. + k))*(t**(-k))*( E0*(t00**k)*(1.+k) + np.heaviside((t-t00),0.5)*Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    #f_ax=(k/t)*(1./(1. + k))*(t**(-k))*( Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    f_ax=(k/t)*(1./(1. + k))*(t**(-k))*( E0*(t00**k)*(1.+k) + Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    f_ax_af_sd = f_ax + L0*(t/t00)**(-delta)
    return np.log10(f_ax_af_sd)


def model_plateau(logt,k,B,omi,t00,Ka1,KLi,L0,collfact,alphax):
    """
    Description:
        Magnetar energy injection only 
        NOTE: To be used only for plotting purposes and it shows the magnetar component only

    Parameters:
        k = coefficient that depends on energy electron fraction and shock dynamical evolution (see Dall'Osso et al.2011)
        B = magnetic field (in units of 10^14 Gauss)
        omi = initial spin frequency (in units of 10^3 Hz)
        t00 = start time of theoretical plateau (very early)
        Ka1 = normalization factor of initial spindown timescale a1 (Ka1=2Ic^3/R^6) where a1=Ka1B^2omi^2
        KLi = normalization factor of initial spindown luminosity Li (KL1=R^6/(6c^3)) where Li=KLiB^2omi^4
        L0 = X-ray luminosity at t00 in units of erg/s
        collfact = collimation factor ratio of afterglow over magentar (1-cos(theta_jet))/(1-cos(theta_magnetar))
        alphax = (3-n)/2 where n is the breaking index of the magnetar (see Stratta et al. 2018)
    Output:
        luminosity in erg/s
    """
    t=10**logt
    a1=Ka1/((B**2.)*(omi**2.))
    Li=KLi*(B**2.)*(omi**4.)
    hg1_a=collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t00/a1)
    hg2_a=collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t/a1)
    f_ax=(k/t)*(1./(1. + k))*(t**(-k))*(Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    return np.log10(f_ax)





