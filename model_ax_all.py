import scipy.special as sc
from sympy.functions import hyper
import numpy as np

def model_ax_all(logt,k,B,omi,t00,Ka1,KLi,L0,collfact,alphax):
    """
    Description:
        Magnetar model (see Stratta et al. 2018 and reference therein)
    Parameters:
        k = coefficient that depends on energy electron fraction and shock dynamical evolution (see Dall'Osso et al.2011)
        B = magnetic field (in units of 10^14 Gauss)
        omi = initial spin frequency (in units of 10^3 Hz)
        t00 = start time of plateau
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
    E0=L0*t00/k
    #E0=KE0*(omi**2)
    hg1_a=collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t00/a1)
    hg2_a=collfact*sc.hyp2f1((alphax-2.)/(alphax-1.), 1.+k, 2.+k, (alphax-1.)*t/a1)
    f_ax=(k/t)*(1./(1. + k))*(t**(-k))*( E0*(t00**k)*(1.+k) + Li*( (t**(1.+k))*hg2_a - (t00**(1.+k))*hg1_a ) )
    return np.log10(f_ax)
