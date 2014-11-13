/*  
 *  Code to implement the Minimum Variance bulk flow estimator of
 *  Watkins, Feldman & Hudson (2009) and Feldman, Watkins & Hudson (2010).
 *
 *  Author: Morag Scrimgeour
 *  11 Oct 2014
 *
 */

#define pi 3.141592653589793

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "cosmology.h"
#include "usefulfuncs.h"
#include "MLE.h"
#include "MV.h"

using namespace std;

void printout_Gnm(int ngal, double** G_nm, char* Gnmfile);

int main()
{
	time_t starttime, endtime, time1;
	starttime = time (NULL); 
	
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------- Setup ----------------------------------------
	//-------------------------------------------------------------------------------------------
	
	// ---------------- DATASET -----------------
	char outfolder[] = "Output";
	char datafile[] = "Data/mock_100_RI10.dat";
	char endstring[] = "_mock_100_RI10";
	char idealsurveyfile[] = "Data/mock_ideal_RI10.dat";
	char pk_file[] = "Data/pk_planck_z0_nonlinear_matterpower.dat";
	int ngal = 100; // Number of galaxies in file
	
	// ------------ FOR IDEAL SURVEY ------------
	double RI = 10.;  // radius of ideal Gaussian
	int NN=1000;    // Number points in ideal survey
	double r_cutoff = 160.;  // Max D_z of 6dFGSv, Planck cosmology
	bool cutoff = false;
	
	// ---------------- OPTIONS -----------------
	bool printweights = false;
	bool print_Gnm = false;
	bool idealreadin = false;
	bool printideal = false;
	
	//--------------- PARAMETERS -----------------
	double OM=0.3175; // Planck value
	double OL=1.-OM;
	//double h=0.6711;  // Planck value
	double H0 = 100.;    // units:  [h km s^-1 Mpc^-1]	
	double dk=0.001;     // For integration
	double kmaxx = 2.;   // Accurate but slow
	double sigmastar2 = 250.*250.;  // (km/s)^2
	
	cout << "" << endl;
	cout << "********************* Minimum Variance bulk flow weights ***********************" << endl;
	cout << "R_I = " << RI << " Mpc/h" << endl;
	cout << "sigmastar2 = " << sigmastar2 << endl;
	
	//--------- Check Outfolder exists ------------
	
	char testfile[50];
	sprintf(testfile,"%s/test.txt",outfolder);
	
	ofstream test_file(testfile);
	if ( !testfile ) { // opened failed
		cout << "Outfolder doesn't exist" << flush << endl;
		exit( -1 );
	}
	test_file.close();
	remove(testfile);
	
	cout << endl;
	cout << "Checked outfolder exists: OK" << flush << endl;
	
	//-------------- Output files ------------------
	
	// Bulk flow components
	char* bffile = new char[100];
	sprintf(bffile,"%s/BF_RI%i%s.dat",outfolder,(int)RI,endstring);
	
	// Survey weights
	char* weightsfile = new char[100];
	sprintf(weightsfile,"%s/weights_RI%i%s.dat",outfolder,(int)RI,endstring);
	
	// Maximum Likelihood weights
	char* MLweightsfile = new char[100];
	sprintf(MLweightsfile,"%s/MLweights_RI%i%s.dat",outfolder,(int)RI,endstring);
	
	// Gnm file
	char* Gnmfile = new char[100];
	sprintf(Gnmfile,"%s/Gnm_RI%i%s.dat",outfolder,(int)RI,endstring);
	
	cout << "Set output file names" << flush << endl;
	
	//-------------------------------------------------------------------------------------------
	//---------------------------------- Read in Velocity Data ----------------------------------
	//-------------------------------------------------------------------------------------------
	
	cout << "Reading in datafile... " << flush;
	
	// Create rihat[3][ngal].  rihat=(xi,yi,zi)/sqrt(xi^2+yi^2+zi^2)
	double** rihat = new double*[3]; //[ngal];
	for (int i=0; i<3; i++){
		rihat[i] = new double[ngal];
	}
	
	// Create ri[3][ngal].  ri=(xi,yi,zi)
	double** ri = new double*[3]; //[ngal];
	for (int i=0; i<3; i++){
		ri[i] = new double[ngal];
	}
	
	double* vrad = new double[ngal];
	double* rad = new double[ngal];
	double* verr = new double[ngal];
	double theta, phi, z, ra, dec;
	int n;
	
	ifstream file1(datafile);
	if ( !file1 ) { // opened failed
		cout << "  cannot open " << datafile << " for read-in" << flush << endl;
		exit( -1 );
	}
	n=0;
	while (file1 >> ra >> dec >> z >> vrad[n] >> verr[n]){\
		// Convert ra, dec, z to Cartesian coords
		rad[n] = dxx_mpch(z, OM, OL);
		phi = ra*pi/180.; // radians
		theta = (90.-dec)*pi/180.; // radians
		rihat[0][n] = sin(theta)*cos(phi);
		rihat[1][n] = sin(theta)*sin(phi);
		rihat[2][n] = cos(theta);
		ri[0][n] = rad[n] * rihat[0][n];
		ri[1][n] = rad[n] * rihat[1][n];
		ri[2][n] = rad[n] * rihat[2][n];
		n++;
	}
	file1.close();
	
	cout << "done" << flush << endl;
	
	//-------------------------------------------------------------------------------------------
	//-------------------------------------- Read in ideal survey --------------------------------
	//-------------------------------------------------------------------------------------------
	//      i.e. random, isotropic set of points with density n(r) \propto exp(-r^2/2R_l^2)
	
	double* xN = new double[NN];
	double* yN = new double[NN];
	double* zN = new double[NN];
	double* rN = new double[NN];
	double factor;
	
	// Create w_in[3][NN]
	double** w_in = new double*[3];
	for (int i = 0; i < 3; i++) {
		w_in[i] = new double[NN];
	}
	
	cout << "Reading in Ideal Survey with " << NN << " points... " << flush;
	
	ifstream file2(idealsurveyfile);
	if ( !file2 ) { // opened failed
		cout << endl << "cannot open " << idealsurveyfile <<  " for read-in" << flush << endl;
		exit( -1 );
	}
	for (int n=0; n<NN; n++){
		file2 >> ra >> dec >> z;
		rN[n] = dxx_mpch(z, OM, OL);
		phi = ra*pi/180.; // radians
		theta = (90.-dec)*pi/180.; // radians
		xN[n] = rN[n]*sin(theta)*cos(phi);
		yN[n] = rN[n]*sin(theta)*sin(phi);
		zN[n] = rN[n]*cos(theta);
		
		factor = 3./(rN[n]*NN);
		
		w_in[0][n] = factor*xN[n];
		w_in[1][n] = factor*yN[n];
		w_in[2][n] = factor*zN[n];
	}
	file2.close();	
	cout << "done" << flush << endl;
	
	//-------------------------------------------------------------------------------------------
	//-------------------------------- Read in CAMB P(k) ----------------------------------------
	//-------------------------------------------------------------------------------------------
	
	int pkintsize;
	double* k_arr;
	double* Pk_arr;
	
	readin_CAMB(pk_file, kmaxx, dk, pkintsize, k_arr, Pk_arr);
	
	//-------------------------------------------------------------------------------------------
	//------------------------------------ MLE Depth  -------------------------------------------
	//-------------------------------------------------------------------------------------------
	
	double MLEdepth;
	cout << "Calculating MLE depth..." << flush << endl;
	
	MLEdepth = MLE_depth(ngal, rad, verr, sigmastar2);
	
	cout << "MLE depth: " << MLEdepth << endl;
	
	//-------------------------------------------------------------------------------------------
	//---------------------------------- Run Minimum Variance method ----------------------------------
	//-------------------------------------------------------------------------------------------
	
	cout << "Setting MV parameters... ";
	//--------------------------------------------------------------
	MV_setparam(OM, H0, sigmastar2);
	//--------------------------------------------------------------
	cout << "done" << endl;
	cout << "OmegaM = " << OM << endl;
	cout << "H0 = " << H0 << endl;
	cout << "sigmastar2 = " << sigmastar2 << endl;
	cout << endl;
	
	//-------------------------------------------------------------------------------------------
	//---------------------------------- Calculate G_nm & P_nm ----------------------------------
	//-------------------------------------------------------------------------------------------
	
	// Create G_mn[ngal][ngal]
	//gsl_matrix * G_nm = gsl_matrix_alloc (ngal, ngal);
	double** G_nm = new double*[ngal];
	for (int i = 0; i < ngal; i++) {
		G_nm[i] = new double[ngal];
	}
	
	// Create G_inv[ngal][ngal]
	double** G_inv = new double*[ngal];
	for (int i = 0; i < ngal; i++) {
		G_inv[i] = new double[ngal];
	}
	
	// Create P_mn[ngal][ngal]
	double** P_nm = new double*[ngal];
	for (int i = 0; i < ngal; i++) {
		P_nm[i] = new double[ngal];
	}
	
	cout << "Calculating G_nm ... " << flush;
	time1 = time (NULL);
	
	//--------------------------------------------------------------
	MV_Gnm(G_nm, P_nm, G_inv, ngal, rad, rihat, verr, pkintsize, dk, k_arr, Pk_arr);
	//--------------------------------------------------------------
	
	cout << "done in " << time (NULL) - time1 << " seconds" << flush << endl;
	
	if (print_Gnm) {
		time1 = time(NULL);
		cout << "Printing Gnm..." << flush << endl;
		printout_Gnm(ngal,G_nm,Gnmfile);
		cout << "done in " << time(NULL) - time1 << " seconds" << flush << endl;
	}

	//-------------------------------------------------------------------------------------------
	//---------------------------------- Calculate Q_i,n ----------------------------------------
	//-------------------------------------------------------------------------------------------
	
	// Create Q_in[3][ngal]
	double** Q_in = new double*[3];
	for (int i = 0; i < 3; i++) {
		Q_in[i] = new double[ngal];
	}
	
	// Create vnvn[NN][ngal]
	double** vnvn = new double*[NN];
	for (int i=0; i<NN; i++){
		vnvn[i] = new double[ngal];
	}
	
	time1 = time(NULL);
	cout << "Calculating Q_in... " << flush;
	
	cout << NN*ngal << " calculations" << endl;
	
	time1 = time (NULL);
	cout << "Calculating vnvn ... " << flush;
	
	//--------------------------------------------------------------
	MV_vnvn(vnvn, ngal, NN, ri, rad, xN, yN, zN, rN, pkintsize, dk, k_arr, Pk_arr);
	//--------------------------------------------------------------
	
	cout << "done in " << time (NULL) - time1 << " seconds" << endl;
	
	//-----------------------------------------
	// Calculate Q_in
	
	time1 = time (NULL);
	cout << "Calculating Q_in ... " << flush;
	
	//--------------------------------------------------------------
	MV_Qin(Q_in, ngal, NN, vnvn, w_in);
	//--------------------------------------------------------------
	
	cout << "done in " << time (NULL) - time1 << " seconds" << flush << endl;

	//-------------------------------------------------------------------------------------------
	//---------------------------------- Calculate M_ij -----------------------------------------
	//-------------------------------------------------------------------------------------------
	
	
	double** M_ij_inv = new double*[3];
	for (int i=0; i<3; i++){
		M_ij_inv[i] = new double[3];
	}
	
	time1 = time(NULL);
	cout << "Calculating M_ij and inverse M_ij... " << flush;
	
	//--------------------------------------------------------------
	MV_Mij(M_ij_inv, G_inv, ngal, rihat);
	//--------------------------------------------------------------
	
	cout << "done in " << time (NULL) - time1 << " seconds" << flush << endl;
	
	//-------------------------------------------------------------------------------------------
	//----------------------- Calculate lambda: Feldman 2010 constraint -------------------------
	//-------------------------------------------------------------------------------------------
	
	// Create lambda_ij[3][3]
	double** lambda_ij = new double*[3];
	for (int i = 0; i < 3; i++) {
		lambda_ij[i] = new double[3];
	}

	cout << "Calculating lambda_ij, Feldman (2010) constraint... " << flush << endl;
	time1=time(NULL);
	//--------------------------------------------------------------
	MV_lambda_ij(lambda_ij, M_ij_inv, G_inv, Q_in, rihat, ngal);
	//--------------------------------------------------------------
	cout << "done in " << time(NULL) - time1 << " seconds" << endl;
	
	//-------------------------------------------------------------------------------------------
	//----------- Calculate w_i,n and bulk flow components: Feldman 2010 constraint -------------
	//-------------------------------------------------------------------------------------------
	
	// Create weight_in[3][ngal]
	double** weight_in = new double*[3];
	for (int i = 0; i < 3; i++) {
		weight_in[i] = new double[ngal];
	}
	
	double* u = new double[3];   // best estimates of moments U
	double sumq;
	double w_ave[3];  // average weight_in
	double sumtest2[3];
	double cosalpha;
	
	time1 = time(NULL);
	cout << "Calculating w_in, Feldman (2010) constraint... " << flush << endl;
	//--------------------------------------------------------------
	MV_weights(weight_in, u, G_inv, Q_in,lambda_ij, rihat, ngal, vrad);
	//--------------------------------------------------------------
	cout << "done in " << time(NULL) - time1 << " seconds" << endl;
	cout << endl;
	
	cout << "Performing sumtest 1: Testing constraint...  (should equal Kronecker delta) " << flush << endl;
	
	// Create sumtest[3][3]
	double** sumtest = new double*[3];
	for (int i = 0; i < 3; i++) {
		sumtest[i] = new double[3];
	}
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			sumtest[i][j]=0.;
			for (int n=0; n<ngal; n++){
				sumtest[i][j] += weight_in[i][n]*rihat[j][n];
			}
			cout << i << "  " << j << "  " << sumtest[i][j] << endl; 
		}
	}
	cout << "done" << endl;
	
	cout << "Performing sumtest 2: Testing k=0 of Window Function... " << flush << endl;
	
	for (int i=0; i<3; i++){
		
		sumtest2[i]=0.;
		for (int n=0; n<ngal; n++){
			for (int m=0; m<ngal; m++){
				if (n == m) cosalpha = 1.;
				else cosalpha = rihat[0][n]*rihat[0][m] + rihat[1][n]*rihat[1][m] + rihat[2][n]*rihat[2][m];
				sumtest2[i] += weight_in[i][n]*weight_in[i][m] * cosalpha/3.;
			}
		}
		
	}
	cout << "done" << endl;
	cout << "W^2(k=0) = (" << sumtest2[0] << ", " << sumtest2[1] << ", " << sumtest2[2] << ")     (should all be =1)"<< endl;
	
	//-------------------------------------------------------------------------------------------
	//----------------------- Calculate Bulk Flow Magnitude & Direction, with errors ------------
	//-------------------------------------------------------------------------------------------
	
	double BFmag, dBFmag;
	double* uerr = new double[3];
	double RA_MV, Dec_MV;
	double* dRADecMV = new double[2];
	
	double dBFmagCV;
	double* uerrCV = new double[3];
	
	BFmag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	RA_MV = getRA(u[0], u[1]);
	Dec_MV = getDec(u[0], u[1], u[2]);
	
	// Noise uncertainties
	dBFmag = BF_err(ngal, weight_in, verr, u, uerr, dRADecMV, sigmastar2);
	
	// Cosmic variance uncertainties
	double** R_ij = new double*[3];
	for (int i=0; i<3; i++){
		R_ij[i] = new double[3];
	}
	//--------------------------------------------------------------
	MV_Rpq(R_ij, G_nm, ngal, weight_in, verr);
	//--------------------------------------------------------------
	double** Re_ij = new double*[3];
	for (int i=0; i<3; i++){
		Re_ij[i] = new double[3];
	}
	//--------------------------------------------------------------
	MV_Re(Re_ij, ngal, weight_in, verr);
	//--------------------------------------------------------------
	double** Rv_ij = new double*[3];
	for (int i=0; i<3; i++){
		Rv_ij[i] = new double[3];
	}
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Rv_ij[i][j] = R_ij[i][j] - Re_ij[i][j];
		}
	}
	uerrCV[0] = sqrt(Rv_ij[0][0]);
	uerrCV[1] = sqrt(Rv_ij[1][1]);
	uerrCV[2] = sqrt(Rv_ij[2][2]);
	
	double result;
	double J[3]; // Jacobian
	J[0] = u[0]/BFmag;
	J[1] = u[1]/BFmag;
	J[2] = u[2]/BFmag;
	double temp[3];
	for (int j=0; j<3; j++){
		temp[j] = 0;
		for (int i=0; i<3; i++){
			temp[j] += J[i]*Rv_ij[i][j];
		}
	}
	result = 0.;
	for (int j=0; j<3; j++){
		result += temp[j] * J[j];
	}
	
	dBFmagCV = sqrt(result);
	
	//-------------------------------------------------------------------------------------------
	//------------------------------------------- Output ----------------------------------------
	//-------------------------------------------------------------------------------------------
	
	cout << "----------------------------------------" << endl;
	cout << "Minimum Variance Bulk Flow: " << endl;
	cout << "ux = " << u[0] << ", uy = " << u[1] << ", uz = " << u[2] << endl;
	cout << "dux = " << uerr[0] << ", duy = " << uerr[1] << ", duz = " << uerr[2] << endl;
	cout << "sqrt(sx^2 + sy^2 + sz^2) = " << sqrt(uerr[0]*uerr[0] + uerr[1]*uerr[1] + uerr[2]*uerr[2]) << endl;
	cout << "Mag: " << BFmag << " +- " << dBFmag << endl;
	cout << "RA: " << RA_MV << " +- " << dRADecMV[0] << " , Dec: " << Dec_MV << " +- " << dRADecMV[1] << endl;
	cout << "----------------------------------------" << flush << endl;
	cout << endl;
	
	cout << "----------------------------------------" << endl;
	cout << "Cosmic Variance uncertainties:" << endl;
	cout << "(" << uerrCV[0] << ", " << uerrCV[1] << ", " << uerrCV[2] << ")" << flush << endl;
	cout << "dmag_CV = " << dBFmagCV << flush << endl;
	cout << "----------------------------------------" << endl;
	
	
	cout << "TEX output:" << endl;
	
	cout << "MV ($R_I = " << (int)RI << "h^{-1}\\,{\\rm Mpc}$)  & $" 
	<< round(BFmag) << " \\pm " <<  round(dBFmag) << "(" << round(dBFmagCV) << ")$ & $"
	<< round(u[0]) << " \\pm " << round(uerr[0]) << "(" << round(uerrCV[0]) << ")$ & $"
	<< round(u[1]) << " \\pm " << round(uerr[1]) << "(" << round(uerrCV[1]) << ")$ & $"
	<< round(u[2]) << " \\pm " << round(uerr[2]) << "(" << round(uerrCV[2]) << ")$ & $"  
	<< round(RA_MV) << " \\pm " << round(dRADecMV[0]) << "$ & $"
	<< round(Dec_MV)<< " \\pm " << round(dRADecMV[1]) 
	<< "$ \\\\" << endl;
	
	
	// Print weights to file
	if (printweights){
		ofstream output_file2(weightsfile);
		output_file2 << "#  w_xn  w_yn  w_z" << endl;
		for (int n = 0; n < ngal; n++){
			output_file2 << weight_in[0][n] << "  " << weight_in[1][n] << "  " << weight_in[2][n] << endl;
		}
		output_file2.close();
		cout << "Printed to " << weightsfile << endl;
		cout << endl;
	}
	
	//------------------------------------------------------------------------
	
	endtime = time (NULL);
	cout << "Time taken = " << (endtime-starttime)/60. << " minutes" << endl; 

	return 0;
}

void printout_Gnm(int ngal, double** G_nm, char* Gnmfile)
{
	cout << "  Opening file..." << flush << endl;
	cout << "  Commencing loop..." << flush << endl;
	ofstream outfile(Gnmfile);
	for (int i=0; i<ngal; i++){
		for (int j=0; j<ngal; j++){
			outfile << i << "  " << j << "  " << G_nm[i][j] << endl;
		}
	}
	outfile.close();
	cout << "  done" << flush << endl;
	
	return;
}
