import numpy as np
from scipy.stats import chisquare
import scipy.stats as stats

"""
PURPOSE:    'define_bestfitparam' extract best fit parameters from data fitting

PARAMETERS:
            popt = best fit model parameter vector (see 'fit_model' module)
            pcov = covariance matrix (see 'fit_model' module)
            x,y,dy = logarithm of rest frame time, X-ray luminosity and errors
            L0 = X-ray luminosity at the start of the plateau
            t0 = start of the plateau
            k, B, omi = free model parameters
            model = the model used (see model_ax module)

OUTPUT:     best fit model parameters (k,B,omi) with 1 sigma uncertainties,
            the associated spin period P in ms, the energy at t0,
            and fit statistics (degrees of freedom, chi square, p-value)
"""


#def define_bestfitparam(popt,pcov,L0=L051,t0=startTxrt,alphax=alphax):
def define(popt,pcov,x,y,dy,L0,t0,alphax,k,B,omi,model):
    k_bf=popt[0]
    dk_bf=pcov[0,0]**0.5
    B_bf=popt[1]
    dB_bf=pcov[1,1]**0.5
    omi_bf=popt[2]
    domi_bf=pcov[2,2]**0.5
    Pms=2.0*np.pi/omi_bf
    dPms=Pms*domi_bf/omi_bf
    E0=L0*t0/popt[0]
    E051=E0/1e51
    print(" ")
    print(' --- Best fit parameters [and initial values] for alphax = ', alphax)
    print(" ")
    print(" k ["+str(k)+"]             = "+"%.5f" %k_bf +" +/- "+ "%.5f" %dk_bf)
    print(" B ["+str(B)+"(10^14 G)]    = "+"%.5f" %B_bf + " +/- "+"%.5f" %dB_bf)
    print(" omi ["+"%.2f" % (omi) +"](10^3 Hz)] = "+ "%.5f" %omi_bf +" +/- "+ "%.5f" %domi_bf)
    print('')
    print(" Spin Period (P=2pi/omi) [ms] = "+ "%.2f" %Pms +" +/- "+ "%.2f" %dPms)
    print(" Ekin(Tstart) = E0 = L0*Tstart/k [(10^51 erg)] = "+"%.4f" %E051)

    #a1_bf=a1_norm*((popt[0]**2.)*(popt[1]**2.))
    #print(" a1 = 2(Ine(c10^3)1.e5)/((r0^6)(B^2)(omi^2)) = "+str(a1_bf))
    #print(" L(Ttstart) = "+ "%2.4f" % Ltstart)
    #print(" Ekin=(L(Tstart))*Tstart/k="+str(10**model(np.log10(startTxrt),popt[0],popt[1],popt[2])*startTxrt/popt[0]))
    #E051=(10**(model(np.log10(startTxrt),popt[0],popt[1],popt[2])))*startTxrt/popt[0]
    print(" ")
    ym=model(x,k_bf,B_bf,omi_bf)
    mychi=sum(((y-ym)**2)/(dy**2))
    dof=len(x)-len(popt)
    print(" my chisquare = "+str(mychi))
    print(" dof = "+str(dof))
    print(" reduced chisquare = "+str(mychi/dof))
    p_value = 1.-stats.chi2.cdf(x=mychi,df=dof)
    print(" P value = "+str(p_value))
    return k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,Pms,dPms,E0,dof,mychi,p_value
