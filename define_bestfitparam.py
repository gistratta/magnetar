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


def define(popt,pcov,x,y,dy,t0,alphax,k,B,omi,model):

    # To be used when efficiency k is a free variable and no steep decay

    k_bf=popt[0]
    dk_bf=pcov[0,0]**0.5
    B_bf=popt[1]
    dB_bf=pcov[1,1]**0.5
    omi_bf=popt[2]
    domi_bf=pcov[2,2]**0.5
    Pms=2.0*np.pi/omi_bf
    dPms=Pms*domi_bf/omi_bf
 
    print(" ")
    print(' --- Best fit parameters [and initial values] for alphax = ', alphax)
    print(" ")
    print(" k ["+"%.4f" % (k) +"]             = "+"%.4f" %k_bf +" +/- "+ "%.4f" %dk_bf)
    print(" B ["+"%.1f" % (B) +"(10^14 G)]    = "+"%.1f" %B_bf + " +/- "+"%.1f" %dB_bf)
    print(" omi ["+"%.1f" % (omi) +"](10^3 Hz)] = "+ "%.1f" %omi_bf +" +/- "+ "%.1f" %domi_bf)
    print('')
    print(" Spin Period (P=2pi/omi) [ms] = "+ "%.1f" %Pms +" +/- "+ "%.1f" %dPms)
    print(" ")

    ym=model(x,k_bf,B_bf,omi_bf)
    mychi=sum(((y-ym)**2)/(dy**2))
    dof=len(x)-len(popt)
    print(" my chisquare = "+str(mychi))
    print(" dof = "+str(dof))
    print(" reduced chisquare = "+str(mychi/dof))
    p_value = 1.-stats.chi2.cdf(x=mychi,df=dof)
    print(" P value = "+str(p_value))
    return k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,Pms,dPms,dof,mychi,p_value


def define_no_k(popt,pcov,x,y,dy,t0,alphax,B,omi,delta,model):

    # To be used when efficiency k is frozen variable and the steep decay is a free variable

    B_bf=popt[0]
    dB_bf=pcov[0,0]**0.5
    omi_bf=popt[1]
    domi_bf=pcov[1,1]**0.5
    delta_bf=popt[2]
    ddelta_bf=pcov[2,2]**0.5    
    Pms=2.0*np.pi/omi_bf
    dPms=Pms*domi_bf/omi_bf
 
    print(" ")
    print(' --- Best fit parameters [and initial values] for alphax = ', alphax)
    print(" ")
    print(" B ["+str(B)+"(10^14 G)]    = "+"%.5f" %B_bf + " +/- "+"%.5f" %dB_bf)
    print(" omi ["+"%.2f" % (omi) +"](10^3 Hz)] = "+ "%.5f" %omi_bf +" +/- "+ "%.5f" %domi_bf)
    print(" P (P=2pi/omi) [ms] = "+ "%.2f" %Pms +" +/- "+ "%.2f" %dPms)
    print(" Steep decay PL index ["+str(delta)+"]    = "+"%.5f" %delta_bf + " +/- "+"%.5f" %ddelta_bf)
    print(" ")

    ym=model(x,B_bf,omi_bf,delta_bf)
    mychi=sum(((y-ym)**2)/(dy**2))
    dof=len(x)-len(popt)
    print(" my chisquare = "+str(mychi))
    print(" dof = "+str(dof))
    print(" reduced chisquare = "+str(mychi/dof))
    p_value = 1.-stats.chi2.cdf(x=mychi,df=dof)
    print(" P value = "+str(p_value))
    return B_bf,dB_bf,omi_bf,domi_bf,delta_bf,ddelta_bf,Pms,dPms,dof,mychi,p_value



def define_steep(popt,pcov,x,y,dy,t0,alphax,k,B,omi,delta,model):
    
    k_bf=popt[0]
    dk_bf=pcov[0,0]**0.5
    B_bf=popt[1]
    dB_bf=pcov[1,1]**0.5
    omi_bf=popt[2]
    domi_bf=pcov[2,2]**0.5
    delta_bf=popt[3]
    ddelta_bf=pcov[3,3]**0.5
    Pms=2.0*np.pi/omi_bf
    dPms=Pms*domi_bf/omi_bf
 
    print(" ")
    print(' --- Best fit parameters [and initial values] for alphax = ', alphax)
    print(" ")
    print(" k ["+str(k)+"]             = "+"%.5f" %k_bf +" +/- "+ "%.5f" %dk_bf)
    print(" B ["+str(B)+"(10^14 G)]    = "+"%.5f" %B_bf + " +/- "+"%.5f" %dB_bf)
    print(" omi ["+"%.2f" % (omi) +"](10^3 Hz)] = "+ "%.5f" %omi_bf +" +/- "+ "%.5f" %domi_bf)
    print(" P (2pi/omi) [ms] = "+ "%.2f" %Pms +" +/- "+ "%.2f" %dPms)
    print(" Steep decay PL index ["+str(delta)+"]             = "+"%.5f" %delta_bf +" +/- "+ "%.5f" %ddelta_bf)
    print(" ")

    ym=model(x,k_bf,B_bf,omi_bf,delta_bf)
    mychi=sum(((y-ym)**2)/(dy**2))
    dof=len(x)-len(popt)
    print(" my chisquare = "+str(mychi))
    print(" dof = "+str(dof))
    print(" reduced chisquare = "+str(mychi/dof))
    p_value = 1.-stats.chi2.cdf(x=mychi,df=dof)
    print(" P value = "+str(p_value))
    return k_bf,dk_bf,B_bf,dB_bf,omi_bf,domi_bf,delta_bf,ddelta_bf,Pms,dPms,dof,mychi,p_value

