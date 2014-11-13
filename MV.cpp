/*
 *  Functions for calculating the Minimum Variance bulk flow
 *
 *  Morag Scrimgeour
 *  11 Oct 2014
 *
 */

#ifndef pi
#define pi 3.141592653589793
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "usefulfuncs.h"
#include "MLE.h"

using namespace std;

/* ----------------------------- Global Variables --------------------------- */

int s, count;
double k, fmnk, alpha, A, A2, cosalpha, sinalpha, kA;
double cosKA, sinKA, kA2;
double c1,c2;
double Gnm, argument;
double integral;
double OmegaMp;
double H0;
double sigma2_star;
double sum;

/* ----------------------------- Ideal survey ------------------------------- */
void MV_idealsurvey(int NN, double* xN, double* yN, double* zN, double* rN,
					double** w_in, char* idealxyz)
{
	char dummyc[100];
	
	// Read in ideal survey
	ifstream file1(idealxyz);
	if ( !file1 ) { // opened failed
		cout << "cannot open " << idealxyz << " for read-in\n" <<  flush;
		exit( -1 );
	}
	file1 >> dummyc >> dummyc >> dummyc >> dummyc >> dummyc >> dummyc >> dummyc >> dummyc;
	for (int n=0; n<NN; n++){
		file1 >> xN[n] >> yN[n] >> zN[n] >> rN[n] >> w_in[0][n] >> w_in[1][n] >> w_in[2][n];
	}
	file1.close();
}

/* ---------------------------- MV_setparam --------------------------- */

void MV_setparam(double& OmegaM, double& H0_in, double& sigma2_star_in)
{
	OmegaMp = pow(OmegaM,1.1);
	H0 = H0_in;
	sigma2_star = sigma2_star_in;
	
	return;
}

/* ---------------------------- MV_Winfun2 --------------------------- */
void MV_W2_pq(double** W2_pq, int n_SNe, double** rihat, double* r, 
			  double** weight_in, int pkintsize, double* k_arr)
// W2_pq: 9 x pkintsize array
//        Indices correspond to
//         0,0[k]   0,1[k]   0,2[k]            [0][k]   [1][k]   [2][k]
//         1,0[k]   1,1[k]   1,2[k]      =     [3][k]   [4][k]   [5][k]
//         2,0[k]   2,1[k]   2,2[k]            [6][k]   [7][k]   [8][k]
{
	
	int count1;
	double val;
	int neg_count;
	
	// Set W2_pq[9][k] to zero
	for (int i=0; i<9; i++){
		for (int kk=0; kk<pkintsize; kk++){
			W2_pq[i][kk] = 0.;
		}
	}
	
	
	count1 = 0;
	neg_count = 0;
	for (int n=0; n<n_SNe; n++){
		for (int m=0; m<n_SNe; m++){
			
			if (n!=m){  // Setup for fmnk
				cosalpha = rihat[0][n]*rihat[0][m] + rihat[1][n]*rihat[1][m] + rihat[2][n]*rihat[2][m];
				alpha = acos(cosalpha);
				if (cosalpha < -1.0 && cosalpha > -1.001) alpha = acos(-1.0);
				if (cosalpha > 1.0 && cosalpha < 1.001) alpha = acos(1.0);
				sinalpha = sin(alpha);
				
				c1 = cosalpha/3.;
				val = r[n]*r[n] + r[m]*r[m] - 2.*r[n]*r[m]*cosalpha;
				if (val > 0.) {
					A = sqrt(val);
					c2 = 1./(A*A)*r[n]*r[m]*sinalpha*sinalpha;
				}
				if (val == 0.) {
					cout << "A = 0 in fmnk at n = " << n << " , m = " << m << endl;
					A = 0.;
				}
				if (val < 0.) {
					cout << "A = sqrt(-ve) in fmnk at n = " << n << " , m = " << m << endl;
					cout << "Have set A to zero, fmnk to 1/3" << endl;
					A=0.;
					neg_count ++;
				}
				
			}
			
			for (int kk=0; kk<pkintsize; kk++){
				
				if (n==m || A==0) fmnk = 1./3.;
				else {
					kA = k_arr[kk]*A;
					sinKA = sin(kA);
					cosKA = cos(kA);
					kA2 = kA*kA;
					
					fmnk = c1*(sinKA/kA*(3.-6./kA2)+6.*cosKA/kA2)  + c2* ((3./kA2-1.)*sinKA/kA - 3.*cosKA/kA2);
					if (fmnk != fmnk) {
						cout << "Error in MV_W2_pq:" << endl;
						cout << "  fmnk = NAN at n = " << n << " , m = " << m << ", k = " << k_arr[kk] << flush << endl;
						cout << "c1 = " << c1 << " c2 = " << c2 << " cosalpha = " << cosalpha << " A = " << A << " val = " << r[n]*r[n] + r[m]*r[m] - 2.*r[n]*r[m]*cosalpha << flush << endl;
						cout << rihat[0][n] << "  " << rihat[0][m] << flush << endl;
						cout << rihat[1][n] << "  " << rihat[1][m] << flush << endl;
						cout << rihat[2][n] << "  " << rihat[2][m] << flush << endl;
						exit(1);
					}
				}
				
				W2_pq[0][kk] += weight_in[0][n]*weight_in[0][m] * fmnk;
				W2_pq[1][kk] += weight_in[0][n]*weight_in[1][m] * fmnk;
				W2_pq[2][kk] += weight_in[0][n]*weight_in[2][m] * fmnk;
				W2_pq[3][kk] += weight_in[1][n]*weight_in[0][m] * fmnk;
				W2_pq[4][kk] += weight_in[1][n]*weight_in[1][m] * fmnk;
				W2_pq[5][kk] += weight_in[1][n]*weight_in[2][m] * fmnk;
				W2_pq[6][kk] += weight_in[2][n]*weight_in[0][m] * fmnk;
				W2_pq[7][kk] += weight_in[2][n]*weight_in[1][m] * fmnk;
				W2_pq[8][kk] += weight_in[2][n]*weight_in[2][m] * fmnk;
				
			}
			
		}
		count1++;
		if (count1 % 10 == 0) cout << count1 << "  " << flush;
	}
	cout << neg_count << " with A=sqrt(-ve) out of " << n_SNe << " x " << n_SNe << endl;
	cout << endl;
	return;
}

/* ---------------------------- MV_Re --------------------------- */
void MV_Re(double** Re_ij, int n_SNe, double** weight_in, double* verr)
// Eq. 24 in Feldman et al. (2010)
{
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Re_ij[i][j] = 0.;
			for (int n=0; n<n_SNe; n++){
				Re_ij[i][j] += weight_in[i][n]*weight_in[j][n]*(verr[n]*verr[n]+sigma2_star);
			}
		}
	}
	
	return;
}

/* ---------------------------- MV_Rv --------------------------- */
void MV_Rv(double** Rv_ij, double** W2_pq, int n_SNe, int pkintsize,  
		   double dk, double* k_arr, double* Pk_arr, double OM)
{
	int indx;  // For indexing W2_pq
	OmegaMp = pow(OM,1.1);
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			sum = 0.;
			indx = 3*i+j;
			for (int kk=0; kk<pkintsize; kk++){
				sum += Pk_arr[kk]*W2_pq[indx][kk]*dk;
			}
			Rv_ij[i][j] = OmegaMp*100.*100./(2.*pi*pi)*sum;
			
		}
	}
	
	return;
}

/* ---------------------------- MV_Rv_dk_arr --------------------------- */
void MV_Rv_dk_arr(double** Rv_ij, double** W2_pq, int n_SNe, int pkintsize, double* dk, 
				  double* k_arr, double* ps0, double OM, double H_0)
{
	int indx;  // For indexing W2_pq
	
	double f = pow(OM,0.55);
	
	for (int i=0; i<3; i++){
		for (int j=0; j<=i; j++){
			sum = 0.;
			indx = 3*i+j;
			for (int kk=0; kk<pkintsize; kk++){
				sum += ps0[kk]*W2_pq[indx][kk]*dk[kk];
			}
			Rv_ij[i][j] = (f*f*H_0*H_0)/(2.*pi*pi)*sum;
		}
	}
	
	// Fill in second half of matrix
	for (int i=0; i<3; i++){
		for (int j=i+1; j<3; j++){
			Rv_ij[i][j] = Rv_ij[j][i];
		}
	}
	
	
	return;
}

/* ---------------------------- MV_Rpq --------------------------- */
void MV_Rpq(double** R_ij, double** G_nm, int n_SNe, double** weight_in, double* verr)
{
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			R_ij[i][j] = 0.;
			for (int n=0; n<n_SNe; n++){
				for (int m=0; m<n_SNe; m++){
					R_ij[i][j] += weight_in[i][n]*weight_in[j][m]*G_nm[n][m];
				}
			}
		}
	}
	
	return;
}

/* -------------------------------- MV_weights --------------------------- */
void MV_weights(double** weight_in, double* u, double** G_inv, double** Q_in, 
				double** lambda_ij, double** rihat, int n_SNe, double* vrad)
{
	
	double sum, sumq;
	
	for (int i=0; i<3; i++){
		u[i]=0.;
		for (int n=0; n<n_SNe; n++){
			sum=0.;
			for (int m=0; m<n_SNe; m++){
				sumq=0.;
				for (int q=0; q<3; q++){
					sumq += lambda_ij[i][q]*rihat[q][m];
				}
				sum += G_inv[n][m]*(Q_in[i][m] - sumq/2.);
			}
			weight_in[i][n] = sum;
			if (weight_in[i][n] != weight_in[i][n]){
				cout << "weight = NAN in weight loop!!" <<  flush << endl; 
				exit(1);
			}
			u[i] += weight_in[i][n] * vrad[n];
		}
	}
	
	cout << "Bulk flow: ux = " << u[0] << " , uy = " << u[1] << " , uz = " << u[2] << endl;
	
	return;
}

/* -------------------------------- MV_lambda_ij --------------------------- */
void MV_lambda_ij(double** lambda_ij, double** M_ij_inv, double** G_inv, 
				  double** Q_in, double** rihat, int n_SNe)
{
	double sum_mn;
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			sum=0.;
			for (int l=0; l<3; l++){
				sum_mn=0.;
				for (int m=0; m<n_SNe; m++){
					for (int n=0; n<n_SNe; n++){
						sum_mn += G_inv[n][m]*Q_in[l][m]*rihat[j][n];
					}
				}
				sum += (sum_mn - Kroneckerdelta(l,j)) * M_ij_inv[i][l];
			}
			lambda_ij[i][j] = sum;
		}
		
	}
	
	return;
}

/* -------------------------------- MV_Mij --------------------------- */
void MV_Mij(double** M_ij_inv, double** G_inv, int n_SNe, double** rihat)
{
	
	double sum;
	double val;
	
	//...................................................
	gsl_matrix * M_ij = gsl_matrix_alloc(3, 3);
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			sum=0.;
			for (int n=0; n<n_SNe; n++){ 
				for (int m=0; m<n_SNe; m++){                // changed from m<=n   (15/6/12, 1:24pm)
					sum += G_inv[n][m]*rihat[i][n]*rihat[j][m];
				}
			}
			gsl_matrix_set(M_ij,i,j,sum/2.);
		}
	}
	
	// Invert M_ij ......................................
	
	gsl_matrix * M_inv = gsl_matrix_alloc(3,3);
	gsl_permutation * perm2 = gsl_permutation_alloc(3);
	
	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp (M_ij, perm2, &s);   // This destroys original matrix, but don't need this again
	
	// Invert the matrix m
	gsl_linalg_LU_invert (M_ij, perm2, M_inv);
	
	// Output M_ij_inv ......................................
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			val = gsl_matrix_get(M_inv,i,j);
			M_ij_inv[i][j] = val;
		}
	}
	
	return;
}

/* ---------------------------- MV_Qin --------------------------- */

void MV_Qin(double** Q_in, int n_SNe, int NN, double** vnvn, double** w_in)
{
	for (int j=0; j<3; j++){
		for (int n=0; n<n_SNe; n++){
			Q_in[j][n]=0.;
			for (int m=0; m< NN; m++) {
				Q_in[j][n] += w_in[j][m] * vnvn[m][n];
			}
			if (Q_in[j][n] != Q_in[j][n]) { cout << j << "  " << n << Q_in[j][n] <<  flush << endl; exit(1); }
		}
		cout << flush << endl;
	}
	return;
}

/* ---------------------------- MV_vnvn --------------------------- */

void MV_vnvn(double** vnvn, int n_SNe, int NN, double** ri, double* r, double* xN, 
			 double* yN, double* zN, double* rN, int pkintsize, double dk, double* k_arr, double* ps0)
{
	double factor = OmegaMp*H0*H0/(2.*pi*pi);
	cout << "factor in vnvn: " << factor << endl;
	
	cout << "vnvn loop: " << n_SNe << " x " << NN << " calculations..." << flush << endl;
	
	for (int n=0; n<n_SNe; n++){
		if (n>100 && n % 100 == 0 || n<100) cout << n << " " << flush;
		
		for (int m=0; m<NN; m++){
			
			cosalpha = (ri[0][n]*xN[m] + ri[1][n]*yN[m] + ri[2][n]*zN[m])/(r[n]*rN[m]);
			
			if (cosalpha < -1.0 && cosalpha > -1.001) cosalpha = -1.0;
			if (cosalpha > 1.0 && cosalpha < 1.001) cosalpha = 1.0;
			
			alpha = acos(cosalpha);   // [0,pi] radians
			sinalpha = sin(alpha);
			
			A2 = r[n]*r[n] + rN[m]*rN[m] - 2.*r[n]*rN[m]*cosalpha;
			A = sqrt(A2);
			
			c1 = cosalpha/3.;
			c2 = (1./A2)*r[n]*rN[m]*sinalpha*sinalpha;
			
			integral=0.;
			for (int i=0; i<pkintsize; i++){
				
				
				k = k_arr[i];  // kmin + i*dk;
				kA = k*A;
				
				sinKA = sin(kA);
				cosKA = cos(kA);
				kA2 = kA*kA;
				
				fmnk = c1*(sinKA/kA*(3.-6./kA2)+6.*cosKA/kA2)  + c2* ((3./kA2-1.)*sinKA/kA - 3.*cosKA/kA2);
	
				integral += ps0[i]*fmnk*dk;
				
			}
			vnvn[m][n] = integral*factor;
			if (vnvn[m][n] != vnvn[m][n]) {cout << "vnvn= NAN!! " << m << "  " << n << "  " << vnvn[m][n] << flush << endl; exit(1);}
			
		}
	}
	cout << "done" << flush << endl;
	
	return;			
}

/* ---------------------------- MV_Gnm --------------------------- */

void MV_Gnm(double** G_nm, double** P_nm, double** G_inv, int n_SNe, double* r, double** rihat, 
			double* verr, int pkintsize, double dk, double* k_arr, double* ps0)
/*
 *  G_nm  -- a nxn matrix
 *  P_nm  -- a nxn matrix
 *  ps0   -- linearly interpolated power spectrum
 */
{
	time_t time2;
	double val; 
	
	cout << "  In MV_Gnm: n_SNe = " << n_SNe << flush << endl;
	cout << "starting..." << flush << endl;
	
	// --------------------- G_nm ------------------------
	for (int n=0; n<n_SNe; n++){
		if (n % 100 == 0) cout << n << "  " << flush;
		for (int m=0; m<=n; m++){
			
			if (n == m){
				A = 0.;
			}
			else {
				cosalpha = rihat[0][n]*rihat[0][m] + rihat[1][n]*rihat[1][m] + rihat[2][n]*rihat[2][m];   // -1 to 1
				if (cosalpha < -1.0 && cosalpha > -1.001) cosalpha = -1.0;
				if (cosalpha > 1.0 && cosalpha < 1.001) cosalpha = 1.0;
				
				alpha = acos(cosalpha);   // [0, pi]
				
				sinalpha = sin(alpha);			
				argument = r[n]*r[n] + r[m]*r[m] - 2.*r[n]*r[m]*cosalpha;			
				A = sqrt(argument);  // If n=m then argument < 0 but don't use A later
				
				c1 = cosalpha/3.;
				c2 = 1./(A*A)*r[n]*r[m]*sinalpha*sinalpha;
			}
			
			integral=0.;
			
			for (int i=0; i<pkintsize; i++){
				
				if (n == m) fmnk = 1./3.;
				else {
					kA = k_arr[i]*A;
					kA2 = kA*kA;
					sinKA = sin(kA);
					cosKA = cos(kA);
					
					fmnk = c1*(sinKA/kA*(3.-6./kA2)+6.*cosKA/kA2)  + c2* ((3./kA2-1.)*sinKA/kA - 3.*cosKA/kA2);

				}
				
				integral += ps0[i]*fmnk*dk;
				
			} //............. ............................................
			
			Gnm = integral*H0*H0*OmegaMp/(2.*pi*pi);
			G_nm[n][m] = Gnm;
			P_nm[n][m] = c1;
			
			if (Gnm != Gnm ) {
				cout << "nan in G_nm! at n=" << n << " m=" << m << endl;  
				cout << integral << flush << endl;
				exit(1);}
			
		}
		Gnm += (sigma2_star+verr[n]*verr[n]);  // Only add this error for n=m, diagonal terms
		G_nm[n][n]=Gnm;
	}
	cout << endl;
	
	// ----------- Fill in second half of matrix ---------
	for (int n=0; n<n_SNe; n++){
		for (int m=n+1; m<n_SNe; m++){
			G_nm[n][m] = G_nm[m][n];
			P_nm[n][m] = P_nm[m][n];
		}
	}
	
	// ----------------- Invert G_nm ------------------------
	
	// Make a GSL copy of G_nm
	gsl_matrix * G_nm2 = gsl_matrix_alloc (n_SNe, n_SNe);
	
	for (int i=0; i<n_SNe; i++){
		for (int j=0; j<n_SNe; j++){
			val = G_nm[i][j];
			gsl_matrix_set(G_nm2, i, j, val);
		}
	}
	
	time2 = time (NULL);
	// Invert G_nm2 ......................................
	gsl_matrix * G_inv0 = gsl_matrix_alloc(n_SNe, n_SNe);
	gsl_permutation * perm = gsl_permutation_alloc(n_SNe);
	
	// Make LU decomposition of matrix G_nm
	gsl_linalg_LU_decomp (G_nm2, perm, &s);   // This destroys original G_nm2, but don't need this again
	// Invert the matrix G_nm
	gsl_linalg_LU_invert (G_nm2, perm, G_inv0);
	
	double G_inv_sum = 0.;
	
	for (int i=0; i<n_SNe; i++){
		for (int j=0; j<n_SNe; j++){
			val = gsl_matrix_get(G_inv0,i,j);
			G_inv[i][j] = val;
			G_inv_sum += val;
		}
	}
	cout << "  Time to invert Gnm = " << time (NULL) - time2 << " seconds" << flush << endl;
	
	cout << "G_inv_sum = " << G_inv_sum << endl;
	return;
}

/* ------------------------------ MLE_CV_err() ------------------------------ */
// Calculate Cosmic Variance uncertainty on MLE bulk flow
// Need to calculate MV_W2_pq, which takes a long time
double MLE_CV_err(int ngal, double* rad, double** rihat, double* verr, double sigmastar2, int pkintsize,
				  double* k_arr, double* Pk_arr, double dk, double OM, double* uML, double* MLerrCV)
{
	double umag, dumag_CV, result;
	
	double** w_in = new double*[3];
	for (int i=0; i<3; i++){
		w_in[i] = new double[ngal];
	}
	
	MLE_weights(ngal, rihat, verr, sigmastar2, w_in);
	
	double** W2_ij_MLE = new double*[9];
	for (int i=0; i<9; i++){
		W2_ij_MLE[i] = new double[pkintsize];
	}
	
	MV_W2_pq(W2_ij_MLE, ngal, rihat, rad, w_in, pkintsize, k_arr);
	
	double** Rv_ij_MLE = new double*[3];
	for (int i=0; i<3; i++){
		Rv_ij_MLE[i] = new double[3];
	}
	
	MV_Rv(Rv_ij_MLE, W2_ij_MLE, ngal, pkintsize, dk, k_arr, Pk_arr, OM);
	
	cout << "Rv_ij_MLE:" << endl;
	cout << "[ " << Rv_ij_MLE[0][0] << " , " << Rv_ij_MLE[0][1] << " , " << Rv_ij_MLE[0][2] << endl;
	cout << "  " << Rv_ij_MLE[1][0] << " , " << Rv_ij_MLE[1][1] << " , " << Rv_ij_MLE[1][2] << endl;
	cout << "  " << Rv_ij_MLE[2][0] << " , " << Rv_ij_MLE[2][1] << " , " << Rv_ij_MLE[2][2] << " ]" << endl;
	
	umag = sqrt(uML[0]*uML[0] + uML[1]*uML[1] + uML[2]*uML[2]);
	
	MLerrCV[0] = sqrt(Rv_ij_MLE[0][0]);
	MLerrCV[1] = sqrt(Rv_ij_MLE[1][1]);
	MLerrCV[2] = sqrt(Rv_ij_MLE[2][2]);
	
	double J[3]; // Jacobian
	J[0] = uML[0]/umag;
	J[1] = uML[1]/umag;
	J[2] = uML[2]/umag;
	double temp[3];
	
	for (int j=0; j<3; j++){
		temp[j] = 0;
		for (int i=0; i<3; i++){
			temp[j] += J[i]*Rv_ij_MLE[i][j];
		}
	}
	
	result = 0.;
	for (int j=0; j<3; j++){
		result += temp[j] * J[j];
	}
	
	dumag_CV = sqrt(result);
	
	return dumag_CV;
}

/* ----------------------------- DRIVER ------------------------------- */
void MV_bulkflow(int n_SNe, double** ri, double* vrad,
				 double* verr, int NN,double* xN, double* yN, double* zN, double* rN,
				 double** w_in, int pkintsize, double dk, double* k_arr, double* ps0,
				 double* u, double* du, double* due, double* BFresults)				 
/*
 *   First call MV_setparam
 *	 And create all matrices
 *   Then call this as many times as necessary
 *
 *   n_SNe : number supernovae or galaxies
 *   ri[][] : x,y,z coords of galaxy, x=ri[0][i], y=ri[1][i], z=ri[2][i]
 *   vrad[] : radial pecular velocity of each galaxy
 *   u[3];   // best estimates of moments U
 *   du[3];  // error on each moment
 *   due[3]; // noise error on each moment
 */
{
	time_t starttime, time1;
	starttime = time (NULL);
	
	cout << "  Creating arrays... " << flush;
	
	// Create r array
	double* r = new double[n_SNe];
	// Create rihat[][] array
	double** rihat = new double*[3];
	for (int i=0; i<3; i++){
		rihat[i] = new double[n_SNe];
	}
	// Fill r and rihat arrays
	for (int i=0; i<n_SNe; i++){
		r[i] = sqrt(ri[0][i]*ri[0][i] + ri[1][i]*ri[1][i] + ri[2][i]*ri[2][i]);
		rihat[0][i] = ri[0][i]/r[i];
		rihat[1][i] = ri[1][i]/r[i];
		rihat[2][i] = ri[2][i]/r[i];
	}
	// Create G_mn[n_SNe][n_SNe]
	//gsl_matrix * G_nm = gsl_matrix_alloc (n_SNe, n_SNe);
	double** G_nm = new double*[n_SNe];
	for (int i = 0; i < n_SNe; i++) {
		G_nm[i] = new double[n_SNe];
	}	
	// Create G_inv[n_SNe][n_SNe]
	double** G_inv = new double*[n_SNe];
	for (int i = 0; i < n_SNe; i++) {
		G_inv[i] = new double[n_SNe];
	}	
	// Create P_mn[n_SNe][n_SNe]
	double** P_nm = new double*[n_SNe];
	for (int i = 0; i < n_SNe; i++) {
		P_nm[i] = new double[n_SNe];
	}
	// Create Q_in[3][n_SNe]
	double** Q_in = new double*[3];
	for (int i = 0; i < 3; i++) {
		Q_in[i] = new double[n_SNe];
	}
	// Create vnvn[NN][n_Sne]
	double** vnvn = new double*[NN];
	for (int i=0; i<NN; i++){
		vnvn[i] = new double[n_SNe];
	}
	double** M_ij_inv = new double*[3];
	for (int i=0; i<3; i++){
		M_ij_inv[i] = new double[3];
	}
	// Create lambda_ij[3][3]
	double** lambda_ij = new double*[3];
	for (int i = 0; i < 3; i++) {
		lambda_ij[i] = new double[3];
	}
	// Create weight_in[3][n_SNe]
	double** weight_in = new double*[3];
	for (int i = 0; i < 3; i++) {
		weight_in[i] = new double[n_SNe];
	}
	double** R_ij = new double*[3];
	for (int i=0; i<3; i++){
		R_ij[i] = new double[3];
	}
	double** Re_ij = new double*[3];
	for (int i=0; i<3; i++){
		Re_ij[i] = new double[3];
	}
	
	cout << "done" << flush << endl;
	
	//.................... Calculate G_nm & P_nm .....................
	time1 = time (NULL);
	MV_Gnm(G_nm, P_nm, G_inv, n_SNe, r, rihat, verr, pkintsize, dk, k_arr, ps0);
	cout << "   Gnm done in " << (time (NULL) - time1)/60. << " mins..." << endl;
	//...................... Calculate Q_i,n .........................
	time1 = time (NULL);
	cout << "  vnvn: " << NN*n_SNe << " calculations" << endl;
	MV_vnvn(vnvn, n_SNe, NN, ri, r, xN, yN, zN, rN, pkintsize, dk, k_arr, ps0);
	cout << "   vnvn done in " << (time (NULL) - time1)/60. << " mins..." << endl;
	MV_Qin(Q_in, n_SNe, NN, vnvn, w_in);
	cout << "   Qin done in " << (time (NULL) - time1)/60. << " mins..." << endl;
	//........................ Calculate M_ij ........................
	time1 = time (NULL);
	MV_Mij(M_ij_inv, G_inv, n_SNe, rihat);
	cout << "   Mij done in " << (time (NULL) - time1)/60. << " mins..." << endl;
	//........................ Calculate lambda_ij ..................
	time1 = time (NULL);
	MV_lambda_ij(lambda_ij, M_ij_inv, G_inv, Q_in, rihat, n_SNe);
	cout << "   lambda_ij done in " << (time (NULL) - time1)/60. << " mins..." << endl;
	//................. Calculate weights & bulk flow ................
	time1 = time (NULL);
	MV_weights(weight_in, u, G_inv, Q_in,lambda_ij, rihat, n_SNe, vrad);
	cout << "   MV weights done in " << (time (NULL) - time1)/60. << " mins..." << endl;
	//....................... Calculate R_ij .........................
	MV_Rpq(R_ij, G_nm, n_SNe, weight_in, verr);
	//....................... Calculate Re_ij .........................
	MV_Re(Re_ij, n_SNe, weight_in, verr);
	// ..................... Calculate errors on u[] ..................
	du[0] = sqrt(R_ij[0][0]);
	du[1] = sqrt(R_ij[1][1]);
	du[2] = sqrt(R_ij[2][2]);
	
	due[0] = sqrt(Re_ij[0][0]);
	due[1] = sqrt(Re_ij[1][1]);
	due[2] = sqrt(Re_ij[2][2]);
	cout << "   errors done..." << endl;
	// ............................ Clear memory .......................
	delete [] G_nm;
	delete [] G_inv;
	delete [] P_nm;
	delete [] Q_in;
	delete [] vnvn;
	delete [] M_ij_inv;
	delete [] lambda_ij;
	delete [] weight_in;
	delete [] R_ij;
	delete [] Re_ij;
	// .................... Calculate BF results .......................
	
	double BFmag, dBFmag, deBFmag;	
	BFmag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	
	sum=0.;
	for (int i=0; i<3; i++){
		sum += (u[i]/BFmag * du[i])*(u[i]/BFmag * du[i]);
	}
	dBFmag = sqrt(sum);
	
	sum=0.;
	for (int i=0; i<3; i++){
		sum += (u[i]/BFmag * due[i])*(u[i]/BFmag * due[i]);
	}
	deBFmag = sqrt(sum);
	
	BFresults[0] = BFmag;
	BFresults[1] = dBFmag;
	BFresults[2] = deBFmag;
	
	return;
}