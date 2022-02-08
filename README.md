## 1. General description

Among the possible interpretations of the "plateau" feature observed in a large
fraction of GRB X-ray afterglows, there is the hypothesis of a newly-born magnetar
pumping energy into the forward shock in the form of spin-down radiation.

Under this assumption, **this code allows to extract magnetar magnetic field and spin period that better describe
the properties of observed plateau in the X-ray afterglow.**

The code is written in Python language and is based on the magnetar model formulated in Stratta et al. 2018, ApJ, 869, 155 

https://iopscience.iop.org/article/10.3847/1538-4357/aadd8f



## 2. To run the code

The code runs under Python 3.8.x and requires the following modules:

> astropy, numpy, matplotlib, scipy,  pandas, sympy

To run the code:

> ipython

> run loglumfit.py

or

> python loglumfit.py

As an example, the code has predefined default settings for the case of the short GRB 130603B (just press "enter" when a question is prompted).


### 2.1 Required input

Most of the input are predefined in "set_iniparam" module (see also section 5). The following ones are asked to user while running the code:

##### GRB name
Once the GRB name is given (e.g. 130603B), the code reads the corresponding data file in the 'data' directory (e.g. ./data/130603B.txt).
The data file must have the Swift/XRT Repository format for the 0.3-10 keV flux light curve including the corresponding photon indexes.

See the guidelines for downloading GRB products provided in the Swift/XRT Repository: https://www.swift.ac.uk/xrt_products/bulk.php


##### GRB jet half-opening angle in deg

If not known, you can set the jet half-opening angle to 90 deg, i.e. isotropic case.

Some useful references where to find GRB  estimated jet opening angles are:
- Wang X-G et al. 2018, Apj, 859, 160 https://doi.org/10.3847%2F1538-4357%2Faabc13
- Wang F. et al. 2020, ApJ, 893, 77 https://iopscience.iop.org/article/10.3847/1538-4357/ab0a86/meta)

##### GRB redshift

Redshift is required to compare the afterglow luminosity with the magnetar one. If not know, you can put it to a fiducial value.

See the useful GRB table by J. Greiner and collaborators:
https://www.mpe.mpg.de/~jcg/grbgen.html


## 3. Output

The code asks if the user wants to save the output. If the answer is yes, then
in the "output" directory two files are generated:

- a log file with the GRB name, the best fit model parameters, the start time of the fit, the GRB and magnetar half-opening angles, and the fit statistics.
The last column shows the time of the fit, so the user can have an overview of multiple runs for the same source.

- a png file with the afterglow light curve and the best fit models for &alpha;=0.9,0.5 and 0.1

The file names indicate
the GRB name, the magnetar half-opening angle and the rest frame energy band
(e.g. 130603B_30.0deg_0.3_30.0keV.log and 130603B_30.0deg_0.3_30.0keV.png)

## 4. The code step by step

1) The code first reads the GRB X-ray afterglow flux and photon indexes vs time and print onto the screen the magnetar predefined setting parameters

2) Then it computes and plots the afterglow beamed corrected luminosity light curve by reading
the rest frame energy band in "set_iniparam" and applying the K-correction,
and by asking for the redshift and the jet opening angle.

3) At this point the code asks to indicate the start and end time if the temporal window we want
to consider for the fit (this step allows to remove the initial steep decay and possible late time steepening, not taken into account by the model)

4) Once the afterglow data set showing the plateau and post-plateau features is defined, the code proceeds in fitting the magnetar model to the data.

5) Results from fit are prompted on the screen and the code asks if the user wants to save the output or not.


## 5. Magnetar default assumptions

There are a number of assumptions on the magnetar properties that are defined in the "set_iniparam.py" module. Each of these assumptions can be changed as desired.

- The magnetar is assumed to radiate all its energy within a cone of half-aperture "thetamdeg=30" [degrees].
- The magnetar mass and radius are set to M = 1.4  [solar mass] and r0 = 1.2  [10 km]
- The magnetar angular momentum direction forms an angle theta_i = 90 [deg] with the magnetic field direction.
- The magnetar bolometric luminosity is approximated with the luminosity computed in the 0.3-30 keV band, defined with the "Er1" and "Er2" variables. The K-correction will takes into account any change of this band.
- The code keeps the magnetar model coefficient k<1 (where k = (4&epsilon;_e dlnt/dlnT), &epsilon;_e is the electron energy fraction and dlnt/dlnT describes the dynamical evolution of the shock, see Dall'Osso et al. 2011)
- The code tests the magnetar model with fixed &alpha; = (3-n)/2 = 0.1, 0.5, 0.9, where n is the braking index (see Stratta et al. 2018)
