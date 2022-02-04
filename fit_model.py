from scipy.optimize import curve_fit
import numpy as np

def fitmodel(model, x, y, dy,k,B,omi):
    # Initial guess for the model parameters (length N)
    p0=np.array([k,B,omi])
    #Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized
    # and the estimated covarianza of popt
    popt, pcov = curve_fit(model, x, y, p0, sigma=dy, bounds=([1.e-12,1.e-5,1.e-6],[1.0, 950.,12.3]),maxfev=10000)
    # best fit parameters
    return popt,pcov
