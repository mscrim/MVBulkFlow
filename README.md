#Minimum Variance Bulk Flow Estimator
##A C++ implementation by Morag Scrimgeour

A C++ implementation of the Minimum Variance bulk flow estimator, described in the papers:
* Watkins, Feldman & Hudson (2009), MNRAS 392 743 ([Arxiv link](http://arxiv.org/abs/0809.4041)) 
* Feldman, Watkins & Hudson (2010), MNRAS 407 2328 ([Arxiv link](http://arxiv.org/abs/0911.5516))

##The Method
The Minimum Variance estimator is a method for calculating the large-scale bulk flow of a galaxy sample, given an input catalogue of galaxy positions and measured peculiar velocities, with estimated Gaussian uncertainties. The details of the method can be found in the above papers. 

##The C++ code
The main C++ function is MV_bulkflow.cpp. Running make will compile the header files and create the executable MV_bulkflow. Mock galaxy data and ideal survey files, along with a sample matter power spectrum generated using [CAMB](http://camb.info/), are provided in the Data folder.

####MLE.cpp
This file contains functions to implement the Maximum Likelihood Estimate of the bulk flow, from Eq. (2.2) of Kaiser (1988), MNRAS 231, 149 ([ADS link](http://adsabs.harvard.edu/abs/1988MNRAS.231..149K)).

####MV.cpp
This file contains all the necessary functions for implementing the Minimum Variance bulk flow method, in Eqs. (8) to (26) in Feldman, Watkins & Hudson (2010).

####usefulfuncs.cpp
This file contains a variety of useful small functions, for both this pipeline and for other uses. Functions include getRA(x,y) and getDec(x,y,z) - functions for calculating Right Ascension and Declination from Cartesian coordinates - and LinearInterpolate, a simple interpolation function.

####cosmology.cpp
This file contains some useful cosmology functions, including dxx_mpch(z,OM,OL), which calculates comoving distance in Mpc/h as a function of redshift, Omega_m and Omega_L.