from scipy.optimize import curve_fit
import numpy as np

def fitmodel(model, x, y, dy,k,B,omi):

    # To be used with 'AP' morphology
    # The efficiency factor k is a free variable
    # There is no steep decay index variable

    p0=np.array([k,B,omi])
    popt, pcov = curve_fit(model, x, y, p0, sigma=dy, bounds=([1e-5,1.e-5,1.e-6],[1.0, 950.,12.3],),maxfev=1000000)
    return popt,pcov

def fitmodel_no_k(model, x, y, dy,B,omi,delta):
    
    # To be used with 'SP' morphology
    # The efficiency factor k is a frozen variable (estimated from empirical broken power law model as k=alpha2-1)
    # The steep decay index is a free variable

    p0=np.array([B,omi,delta])
    popt, pcov = curve_fit(model, x, y, p0, sigma=dy, bounds=([1.e-5,1.e-6,1.0],[950.,12.3,10.0]),maxfev=1000000)
    return popt,pcov


def fitmodel_steep(model, x, y, dy,k,B,omi,delta):
    p0=np.array([k,B,omi,delta])
    popt, pcov = curve_fit(model, x, y, p0, sigma=dy, bounds=([1e-5,1.e-5,1.e-6,1.0],[1.0, 950.,12.3,10.0]),maxfev=1000000)
    return popt,pcov
