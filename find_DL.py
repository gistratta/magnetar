import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

"""
PURPOSE: 'findDL' compute the luminosity distance by assuming a flat lambda CDM cosmology
PARAMETER: the redshift
OUTPUT: luminosity distance in cm 
"""

def findDL(z):
    cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    DL=(cosmo.luminosity_distance(z)).value*3.08568*1e24
    return DL
