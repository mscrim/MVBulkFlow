/*
 *  Modules for calculating the Maximum Likelihood Estimate bulk flow
 *
 *  Morag Scrimgeour
 *  11 Oct 2014
 *
 */

#ifndef pi
#define pi 3.141592653589793
#endif

#include <math.h>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "usefulfuncs.h"

using namespace std;

/* ------------------------------ BF_err() ------------------------------ */
// Calculate uncertainty on the Bulk Flow magnitude
double BF_err(int n_SNe, double** w_in, double* verr, double* u, 
			  double* uerr, double* RADecerr, double sigmastar2)
/*
 Calculate dumag = J R_ij J^T
 INPUT:
 n_SNe: number of galaxies/supernovae
 w_in: calculated bulk flow weights, w_in[3][n_SNe]
 vrad: measured radial peculiar velocities, vrad[n_SNe]
 verr: uncertainties on peculiar velocities, verr[n_SNe]
 u: calculated bulk flow components, u[3]
 sigmastar2: sigma_star^2
 OUTPUT:
 uerr: uncertainties on bulk flow components, uerr[3]
 RADecerr: uncertainties on RA and Dec, RADecerr[2]
 RETURN:
 dumag: uncertainty on bulk flow magnitude
 */
{
	//-------------------- Magnitude error ------------------------
	
	double umag, dumag, Wsum, result;
	
	double** Re_ij = new double*[3];
	for (int i = 0; i < 3; i++) {
		Re_ij[i] = new double[3];
	}
	
	umag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	
	// Calculate Re_ij
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Wsum = 0.;
			for (int n=0; n<n_SNe; n++){
				Wsum += w_in[i][n] * w_in[j][n] * (verr[n]*verr[n] + sigmastar2);
			}
			Re_ij[i][j] = Wsum;
		}
	}
	
	cout << "Re_ij:" << endl;
	cout << "[ " << Re_ij[0][0] << " , " << Re_ij[0][1] << " , " << Re_ij[0][2] << endl;
	cout << "  " << Re_ij[1][0] << " , " << Re_ij[1][1] << " , " << Re_ij[1][2] << endl;
	cout << "  " << Re_ij[2][0] << " , " << Re_ij[2][1] << " , " << Re_ij[2][2] << " ]" << endl;
	
	uerr[0] = sqrt(Re_ij[0][0]);
	uerr[1] = sqrt(Re_ij[1][1]);
	uerr[2] = sqrt(Re_ij[2][2]);
	
	double J[3]; // Jacobian
	J[0] = u[0]/umag;
	J[1] = u[1]/umag;
	J[2] = u[2]/umag;
	double temp[3];
	
	for (int j=0; j<3; j++){
		temp[j] = 0;
		for (int i=0; i<3; i++){
			temp[j] += J[i]*Re_ij[i][j];
		}
	}
	
	result = 0.;
	for (int j=0; j<3; j++){
		result += temp[j] * J[j];
	}
	
	dumag = sqrt(result);
	
	//--------------------- RA, Dec error -------------------------
	
	double x = u[0];
	double x2 = x*x;
	double y = u[1];
	double y2 = y*y;
	double z = u[2];
	double z2 = z*z;
	double r = umag;
	double r2 = r*r;
	double r3 = r*r*r;
	
	// Jacobian J(RA, Dec)
	double** J_RADec = new double*[2];
	for (int i = 0; i < 2; i++) {
		J_RADec[i] = new double[3];
	}
	J_RADec[0][0] = 180./pi * (-y)/(x2*(1.+y2/x2));
	J_RADec[0][1] = 180./pi * 1./(x*(1.+y2/x2));
	J_RADec[0][2] = 0.;
	J_RADec[1][0] = 180./pi * (-z*x)/(r3*sqrt(1.-z2/r2));
	J_RADec[1][1] = 180./pi * (-z*y)/(r3*sqrt(1.-z2/r2));
	J_RADec[1][2] = 180./pi * (1./sqrt(1.-z2/r2)) * (z2/r3-1./r);
	
	// Uncertainty on RA
	for (int j=0; j<3; j++){
		temp[j] = 0;
		for (int i=0; i<3; i++){
			temp[j] += J_RADec[0][i]*Re_ij[i][j];
		}
	}
	result = 0.;
	for (int j=0; j<3; j++){
		result += temp[j] * J_RADec[0][j];
	}
	RADecerr[0] = sqrt(result);
	
	// Uncertainty on Dec
	for (int j=0; j<3; j++){
		temp[j] = 0;
		for (int i=0; i<3; i++){
			temp[j] += J_RADec[1][i]*Re_ij[i][j];
		}
	}
	result = 0.;
	for (int j=0; j<3; j++){
		result += temp[j] * J_RADec[1][j];
	}
	RADecerr[1] = sqrt(result);
	
	return dumag;
}
/* ------------------------------ MLE_BF() ------------------------------ */
double MLE_BF(int n_SNe, double** rihat, double* verr, double sigmastar2,
			  double* vrad, double* uML, double* uMLerr, double& dMLmag, 
			  double& RA_ML, double& Dec_ML, double* dRADecML, 
			  char* MLweightsfile, bool printopt, bool coutopt)
/*
 This routine calculates the Maximum Likelihood Estimate (MLE) bulk flow
 componenets uML[3], and their errors, MSerr[3]
 n_SNe  -- number SNe or galaxies
 rihat  -- [3][n_SNe] array, rihat=(xi,yi,zi)/sqrt(xi^2+yi^2+zi^2)
 verr   -- root mean square variance velocity error
 sigmastar2  --  sigma_* x sigma_*
 vrad   -- [n_SNe] array, radial peculiar velocity 
 uML    -- [3] array, x,y,z component of MLE bulk flow
 uMLerr  -- [3] array, root mean square variance of uML
 dMLmag -- uncertainty on magnitude of bulk flow
 RA_ML -- RA of bulk flow (or gl if using Galactic coordinates)
 Dec_ML -- Dec of bulk flow (or gb if using Galactic coordinates)
 dRADecML -- [2] array, uncertainty on RA(gl) and Dec(gb)
 MLweightsfile -- output file to print MLE weights
 printopt -- 1 to print, 0 to not print to file
 coutopt -- 1 to cout Aij and Aijinv matrices, 0 to not cout
 Returns MLmag
 */
{
	double MLmag;
	double Aij;
	double weightin;
	int s;
	
	// Create A_ij[3][3]
	gsl_matrix * A_ij = gsl_matrix_alloc (3, 3);
	
	// Create ML_weight[3][n_SNe]
	double** ML_weight = new double*[3];
	for (int i = 0; i < 3; i++) {
		ML_weight[i] = new double[n_SNe];
	}
	
	// Create A_ij
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Aij=0.;
			for (int n=0; n<n_SNe; n++){
				//if (i == 0  && j == 0) cout << n << "  " << (rihat[i][n] * rihat[j][n]) / (verr[n]*verr[n] + sigmastar2) << endl;
				Aij += (rihat[i][n] * rihat[j][n]) / (verr[n]*verr[n] + sigmastar2) ;
				//cout << i << "  " << j << "  " << (rihat[i][n] * rihat[j][n]) / (verr[n]*verr[n] + sigmastar2) << endl;
			}
			gsl_matrix_set(A_ij,i,j,Aij);
		}
	}
	
	if (coutopt) {
		cout << "  " << endl;
		cout << "A_ij:" << endl;
		cout << gsl_matrix_get(A_ij,0,0) << "  " <<gsl_matrix_get(A_ij,0,1) << "  " << gsl_matrix_get(A_ij,0,2) << endl;
		cout << gsl_matrix_get(A_ij,1,0) << "  " <<gsl_matrix_get(A_ij,1,1) << "  " << gsl_matrix_get(A_ij,1,2) << endl;
		cout << gsl_matrix_get(A_ij,2,0) << "  " <<gsl_matrix_get(A_ij,2,1) << "  " << gsl_matrix_get(A_ij,2,2) << endl;
		cout << "" << endl;
	}
	
	// Invert A_ij......................................
	gsl_matrix * A_inv = gsl_matrix_alloc(3, 3);
	gsl_permutation * perm1 = gsl_permutation_alloc(3);
	
	// Make LU decomposition of matrix A_ij
	gsl_linalg_LU_decomp (A_ij, perm1, &s);   // This destroys original A_ij, but don't need this again
	
	// Invert the matrix A_ij
	gsl_linalg_LU_invert (A_ij, perm1, A_inv);
	//..................................................
	
	if (coutopt){
		cout << "Inverse A_ij:" << endl;
		cout << gsl_matrix_get(A_inv,0,0) << "  " <<gsl_matrix_get(A_inv,0,1) << "  " << gsl_matrix_get(A_inv,0,2) << endl;
		cout << gsl_matrix_get(A_inv,1,0) << "  " <<gsl_matrix_get(A_inv,1,1) << "  " << gsl_matrix_get(A_inv,1,2) << endl;
		cout << gsl_matrix_get(A_inv,2,0) << "  " <<gsl_matrix_get(A_inv,2,1) << "  " << gsl_matrix_get(A_inv,2,2) << endl;
		cout << "" << endl;
	}
	
	// Calculate ML_weights
	for (int i=0; i<3; i++){
		uML[i] = 0.;
		for (int n=0; n<n_SNe; n++){
			weightin = 0.;
			for (int j=0; j<3; j++){
				weightin += gsl_matrix_get(A_inv,i,j) * rihat[j][n] / (verr[n]*verr[n] + sigmastar2);
			}
			ML_weight[i][n] = weightin;
			uML[i] += weightin*vrad[n];  // 1./(verr[n]*verr[n] + sigmastar2)*rihat
			//if (i == 0) cout << weightin << "  " << vrad[n] << "  " << uML[i] << endl;
		}
	}
	
	// Calculate Bulk Flow magnitude & direction
	MLmag = sqrt(uML[0]*uML[0] + uML[1]*uML[1] + uML[2]*uML[2]);
	RA_ML = getRA(uML[0], uML[1]);
	Dec_ML = getDec(uML[0], uML[1], uML[2]);
	
	// Bulk flow uncertainties (returns dMLmag, uMLerr and dRADecML)
	dMLmag = BF_err(n_SNe, ML_weight, verr, uML, uMLerr, dRADecML, sigmastar2);
	
	if (coutopt){
		cout << "uMLerr: " << uMLerr[0] << " , " << uMLerr[1] << " , " << uMLerr[2] << endl;
		cout << "dMLmag in function = " << dMLmag << endl;
		cout << "RA = " << RA_ML << " +- " << dRADecML[0] << endl;
		cout << "Dec = " << Dec_ML << " +- " << dRADecML[1] << endl;
		cout << "" << endl;
	}
	
	// Print Maximum Likelihood weights to file
	if (printopt){
		FILE* ML_file;
		ML_file = fopen(MLweightsfile,"w");
		fprintf(ML_file, "# ML_weight[x]  ML_weight[y]  ML_weight[z]\n");
		for (int n=0; n<n_SNe; n++){
			fprintf(ML_file, "%.10f  ",ML_weight[0][n]);  
			fprintf(ML_file, "%.10f  ",ML_weight[1][n]);  
			fprintf(ML_file, "%.10f\n",ML_weight[2][n]);
		}
		fclose(ML_file);
		cout << "Printed to " << MLweightsfile << endl;
	}
	
	return MLmag;
}

/* ------------------------------ MLE_BF() ------------------------------ */
void MLE_weights(int n_SNe, double** rihat, double* verr, double sigmastar2, double** ML_weight)
/*
 This routine calculates the Maximum Likelihood Estimate (MLE) weights
 n_SNe  -- number SNe or galaxies
 rihat  -- [3][n_SNe] array, rihat=(xi,yi,zi)/sqrt(xi^2+yi^2+zi^2)
 verr   -- root mean square variance velocity error
 sigmastar2  --  sigma_* x sigma_*
 w_in -- MLE weights, 3*n array
 */
{
	double Aij;
	double w_in;
	int s;
	
	// Create A_ij[3][3]
	gsl_matrix * A_ij = gsl_matrix_alloc (3, 3);
	
	// Create A_ij
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Aij=0.;
			for (int n=0; n<n_SNe; n++){
				Aij += (rihat[i][n] * rihat[j][n]) / (verr[n]*verr[n] + sigmastar2) ;
			}
			gsl_matrix_set(A_ij,i,j,Aij);
		}
	}
	
	// Invert A_ij......................................
	gsl_matrix * A_inv = gsl_matrix_alloc(3, 3);
	gsl_permutation * perm1 = gsl_permutation_alloc(3);
	
	// Make LU decomposition of matrix A_ij
	gsl_linalg_LU_decomp (A_ij, perm1, &s);   // This destroys original A_ij, but don't need this again
	
	// Invert the matrix A_ij
	gsl_linalg_LU_invert (A_ij, perm1, A_inv);
	//..................................................
	
	// Calculate ML_weights
	for (int i=0; i<3; i++){
		for (int n=0; n<n_SNe; n++){
			w_in = 0.;
			for (int j=0; j<3; j++){
				w_in += gsl_matrix_get(A_inv,i,j) * rihat[j][n] / (verr[n]*verr[n] + sigmastar2);
			}
			ML_weight[i][n] = w_in;
		}
	}
	
	return;
}

/* ------------------------------ MLE_depth() ------------------------------ */
double MLE_depth(int n_SNe, double* rnew, double* verr, double sigmastar2)
{
	double MLEdepth;
	double weightn;
	double sum1, sum2;
	
	sum1 = 0.;
	sum2 = 0.;
	
	for (int i=0; i<n_SNe; i++){
		weightn = 1. / (verr[i]*verr[i] + sigmastar2);
		sum1 += rnew[i] * weightn;
		sum2 += weightn;
	}
	
	MLEdepth = sum1 / sum2;
	
	return MLEdepth;
}
/* ------------------------------ Simple Sum bulk flow ------------------------------ */
double simple_sum_BF(int n_SNe, double** rihat, double* vrad, double* verr, double* u, 
					 double* u_err, double& RA_SS, double& Dec_SS, double* dRADec, double& dumag)
{
	
	double umag, w_0n, w_1n, w_2n;
	
	double** w_in = new double*[3];
	for (int i = 0; i < 3; i++) {
		w_in[i] = new double[n_SNe];
	}
	
	u[0]=0.;
	u[1]=0.;
	u[2]=0.;
	
	for (int i=0; i<n_SNe; i++){
		w_0n = rihat[0][i] / n_SNe;
		w_1n = rihat[1][i] / n_SNe;
		w_2n = rihat[2][i] / n_SNe;
		
		w_in[0][i] = w_0n;
		w_in[1][i] = w_1n;
		w_in[2][i] = w_2n;
		
		u[0] += w_0n * vrad[i];
		u[1] += w_1n * vrad[i];
		u[2] += w_2n * vrad[i];
	}
	
	umag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	RA_SS = getRA(u[0], u[1]);
	Dec_SS = getDec(u[0], u[1], u[2]);
	
	dumag = BF_err(n_SNe, w_in, verr, u, u_err, dRADec, 0.);
	
	return umag;
}

//*****************************************************************************************
double MLE_BF_simple(int ngal, double* vrad,double* verr,double** rihat,  double sigmastar2, 
					 double* uML, double* MLerr, double& dMLmag, double& RA_ML, double& Dec_ML, 
					 double* dRADecML)
{
	double Aij, w_in, MLmag;
	int s;
	
	// Create ML_weight[3][n_SNe]
	double** ML_weight = new double*[3];
	for (int i = 0; i < 3; i++) {
		ML_weight[i] = new double[ngal];
	}
	
	// Create gsl matrix A_ij[3][3]
	gsl_matrix * A_ij = gsl_matrix_alloc(3, 3);
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Aij = 0.;
			for (int n=0; n<ngal; n++){
				Aij += ((rihat[i][n]*rihat[j][n])/(verr[n]*verr[n] + sigmastar2));
			}
			gsl_matrix_set(A_ij,i,j,Aij);
		}
	}
	
	// Invert A_ij......................................
	gsl_matrix * A_inv = gsl_matrix_alloc(3, 3);
	gsl_permutation * perm = gsl_permutation_alloc(3);
	
	// Make LU decomposition of matrix A_ij
	gsl_linalg_LU_decomp (A_ij, perm, &s);   // This destroys original A_ij, but don't need this again
	
	// Invert the matrix A_ij
	gsl_linalg_LU_invert (A_ij, perm, A_inv);
	//..................................................
	
	// Calculate w_in, uML & MLerr 
	for (int i=0; i<3; i++){
		uML[i] = 0.;
		for (int n=0; n<ngal; n++){
			w_in = 0.;
			for (int j=0; j<3; j++){
				w_in += gsl_matrix_get(A_inv,i,j) * rihat[j][n] / (verr[n]*verr[n] + sigmastar2);
			}
			uML[i] += w_in*vrad[n];
			ML_weight[i][n] = w_in;
		}
	}
	
	// MLE bulk flow amplitude
	MLmag = sqrt(uML[0]*uML[0] + uML[1]*uML[1] + uML[2]*uML[2]);
	RA_ML = getRA(uML[0], uML[1]);
	Dec_ML = getDec(uML[0], uML[1], uML[2]);
	
	// Error in mag, dmag
	dMLmag = BF_err(ngal, ML_weight, verr, uML, MLerr, dRADecML, sigmastar2);
	
	return MLmag;
}
