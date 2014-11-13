#ifndef MLE_H
#define MLE_H

double BF_err(int n_SNe, double** w_in, double* verr, double* u, 
			  double* uerr, double* RADecerr, double sigmastar2);

double MLE_BF(int n_SNe, double** rihat, double* verr, double sigmastar2,
			  double* vrad, double* uML, double* uMLerr, double& dMLmag, 
			  double& RA_ML, double& Dec_ML, double* dRADecML, 
			  char* MLweightsfile, bool printopt, bool coutopt);

void MLE_weights(int n_SNe, double** rihat, double* verr, double sigmastar2, double** w_in);

double MLE_depth(int n_SNe, double* rnew, double* verr, double sigmastar2);

double simple_sum_BF(int n_SNe, double** rihat, double* vrad, double* verr, double* u, 
					 double* u_err, double& RA_SS, double& Dec_SS, double* dRADec, double& dumag);

double MLE_BF_simple(int ngal, double* vrad,double* verr,double** rihat,  double sigmastar2,
					 double* uML, double* MLerr, double& dMLmag);

#endif