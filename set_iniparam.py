import numpy as np
import astropy.units as u
from astropy import constants as const

"""
# ==> SET MODEL INITIAL PARAMETERS
"""
def setiniparam():

    # --- rest frame energy range at which compute the afterglow luminosity
    Er1=0.3
    Er2=30.0

    # --- Astrophysical Constants
    msun33=(const.M_sun.cgs.value)/1.e33      # units of 10^33 gr
    G=const.G.cgs.value                       # Gravitational constant
    c10=(const.c.cgs.value)/1.e10             # units of 10^10 cm/s

    # --- Magnetar collimation
    thetamdeg=30
    thetam=thetamdeg*np.pi/180

    # --- Neutron star Constants
    r0 =1.2                                                  # units of 10^6 cm (r0=12km)
    M = 1.4
    be = G*msun33*1.e7*M/(r0*(c10**2))                         # compactness
    Ine = (0.247+(0.642*be)+(0.466*(be**2)))*msun33*M*(r0**2)  # Inertia from Lattimer & Prakash in units of 10^45 gr cm^2
    a1_norm = 2.*(Ine*(c10**3)*1.e5)/(r0**6)                 # units of 10^39 where a1=a1_norm/(B^2omi^2) [s]
    theta_i = 90*np.pi/180                                   # angle between NS angular momentum and magnetic field
    Li_norm = (((r0**2)/2)**2)*((1+(np.sin(theta_i))**2)/(c10**3))*1e46  #Li=Li_norm*B^2omi^4
    #E0_norm = 0.5*Ine*1e51                                   # initial spindown energy norm. factor, where E0 = E0_norm*omi^2

    # --- NS Initial values of model parameters
    B = 10.                   # Magnetic field in units of 10^14 Gauss= (gr/cm)^(1/2)*s^-1
    spini = 5.                # initial spin period in units of ms >1ms
    omi = 2.0*np.pi/spini     # initial spin frequency 2pi/spini in units of 10^3 Hz
    k = 0.1                     # coefficient that depends on energy electron fraction and shock dynamical evolution (see Dall'Osso et al.2011)
    alphax = 0.9

    print('   Initial setting:')
    print('   ---')
    print('   X-ray afterglow luminosity rest frame band [keV]: ',Er1,Er2)
    print('   ---')
    print('   Magnetar properties:')
    print('     Jet angle of the magnetar [deg]: ',str(thetamdeg))
    print('     NS mass: M [solar mass] = ', M)
    print('     NS radius: r0 [10^6 cm (10km)] = ', r0)
    print('     NS compactnes: be = ', round(be,2))
    print('     Inertia from Lattimer and Prakash: Ine = (0.247+(0.642be)+(0.466(be^2)) Msun33 M r0^2 [10^45 gr cm^2] = ', round(Ine,2))
    print('     B vs omi inclination angle [deg]:', round(theta_i*180/np.pi,0))
    print('   ---')
    print('   Magnetar model initial values:')
    print('         Magnetic field strenght: B [10^14 G] = ',B)
    print('         Spin: spini [ms] = ', spini)
    print('         k (radiative eff. k=4x$\epsilon_e$):', k)
    print('         alpha:', alphax)
    print('   ---')
    print('')
    return Er1,Er2,thetam,alphax,k,omi,spini,B,a1_norm,Li_norm
