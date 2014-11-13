/*
 *  Functions for Cosmology
 *
 *  Morag I. Scrimgeour
 *  11 Oct 2014
 *
 */

#ifndef pi
#define pi 3.141592653589793
#endif

#include <math.h>
#include <iostream>
#include <fstream>
#include "usefulfuncs.h"

using namespace std;

//********************************************************************************
/* Function DXX - returns comoving distance dxx(z) in Mpc/h */
double dxx_mpch(double zz, double OM, double OL)
{
	double h,OK,ORR,w0,wa,H0,H0y,c_MpcGyr,xx;
	double redshift,a,ezinv;
	double dz = 1e-8;
	
	h = 1;
	
	OK = 1.-OM-OL;
	
	H0 = h*100.; // km/s/Mpc
	H0y = H0 * 1.0225e-3; // Gyr^-1 3.24e-20 * 3.156e16
	
	c_MpcGyr = 3e8 / 3.09e22 * 3.1536e16; // [m/s] / [m/Gpc] * [s/Gyr]
	
	xx=0;
	redshift=0;
	
	// Do integration
	while (redshift <= zz){
		a = 1./(1.+redshift);
		ezinv = 1./sqrt(OM/(a*a*a) + OK/(a*a) + OL);
		
		xx += ezinv*dz;
		
		redshift += dz;
	}
	
	xx = xx * c_MpcGyr / H0y;
	
	return xx;
}

//********************************************************************************
/* Function DXX_MPCH_ARR - returns comoving distance dxx(z) in Mpc/h for an 
   array of redshifts (saves a huge amount of time). */
double* dxx_mpch_arr(int num, double* zz, double OM, double OL)
{
	double h,OK,ORR,w0,wa,H0,H0y,c_MpcGyr,xx;
	double redshift,a,ezinv;
	double dz = 1e-8;
	double zshift, zlast;
	double* dxx = new double[num];
	
	h = 1;
	OK = 1.-OM-OL;
	H0 = h*100.; // km/s/Mpc
	H0y = H0 * 1.0225e-3; // Gyr^-1 3.24e-20 * 3.156e16
	c_MpcGyr = 3e8 / 3.09e22 * 3.1536e16; // [m/s] / [m/Gpc] * [s/Gyr]
	
	// Sort zz
	int* zindx = new int[num];
	sort_arr(zz, num, zindx, 1, 0);  // Doesn't reorder zz
	
	zlast = 0.;
	xx=0.;
	for (int i=0; i<num; i++){
		redshift=zlast;
		zshift = zz[zindx[i]];
		// Do integration
		while (redshift <= zshift){
			a = 1./(1.+redshift);
			ezinv = 1./sqrt(OM/(a*a*a) + OK/(a*a) + OL);
			
			xx += ezinv*dz;
			
			redshift += dz;
		}
		dxx[zindx[i]] = xx * c_MpcGyr / H0y;
		zlast = redshift - dz;
	}
	
	return dxx;
}

//********************************************************************************
void readin_CAMB(char* camb_file, double kmaxx, double dk, int& pkintsize, 
				 double *& k_arr, double *& Pk_arr)
/* Function to read in CAMB file, and linearly interpolate it to return k_arr and 
   Pk_arr, which can be used for integrations over dk.

   camb_file : matter power spectrum file from CAMB
   kmaxx : desired kmax. If kmaxx=0, use kmax from file.
   dk: desired dk
   pkintsize: size of interpolated P(k) array
 
   In main():
      double dk=0.001;
      int pkintsize;
      double* k_arr;
      double* Pk_arr;
      readin_CAMB(camb_file, kmaxx, dk, pkintsize, k_arr, Pk_arr);
 */
{
	cout << "  Reading in CAMB file ..........." << flush << endl;
	cout << "    " << camb_file << flush << endl;
	
	int pksize, nn;
	double dummy;
	
	pksize = 0.;
	// Read in CAMB power spectrum (1) Find size of file
	ifstream file(camb_file);
	if ( !file ) { // opened failed
		cout << "    cannot open CAMB file for read-in" << flush << endl;
		exit( -1 );
	}
    while(file >> dummy >> dummy){
		pksize++;
	}
	file.close();
	
	cout <<"    pksize=" << pksize << endl;
	
	double k0, Pk0;
	double lnk[pksize], lnPk[pksize];
	double kmax = 0., kmin=1000000.;
	
	nn=0;
	ifstream file2(camb_file);
    while(file2 >> k0 >> Pk0){   // k0: h Mpc^{-1}  Pk0: h^3 Mpc^{-3}
		lnk[nn] = log(k0);
		lnPk[nn] = log(Pk0);
		if (k0 > kmax) kmax = k0;
		if (k0 < kmin) kmin = k0;
		nn++;
	}
	file2.close();
	
	cout << "    kmax=" << kmax << ", kmin = " << kmin;	
	if (kmaxx > 0.) {
		kmax = kmaxx;
		cout << ", imposed kmax = " << kmax << endl;
	}
	else cout << ", no imposed kmax" << endl;
	
	pkintsize = (kmax-kmin)/dk;
	
	cout <<"    interpolating with dk = " << dk << endl;
	cout <<"    new pkintsize after kmax cut=" << pkintsize << endl;
	
	k_arr = new double[pkintsize];
	Pk_arr = new double[pkintsize];
	
	cout << "  interpolating P(k)... " << flush;
	// Interpolate P(k) outside loop
	for (int i=0; i<pkintsize; i++){
		k_arr[i] = kmin + (double)i*dk;
		if (k_arr[i] > kmax) {cout << "kmaxx set too high for CAMB read-in" << flush << endl; exit(1);}
		Pk_arr[i] = exp(LinearInterpolate(pksize, lnk, lnPk, log(k_arr[i])));
	}
	
	cout << "done ....." << flush << endl;
	
	return;
}

//********************************************************************************
double H_z(double z, double OM0, double H0, int unit)
/* Returns Hubble constant as a function of redshift, H(z)
   unit = 1: Gyr-1
   unit = 2: km/s/Mpc
   unit = 3: s-1
 */
{
	double OL0 = 1.-OM0;
	double H0s, H0y, H0kmsmpc;
	double h_zz;
	
	H0kmsmpc = H0;
	
	// Get H0 in inverse seconds
	H0s = H0kmsmpc/(3.08567758e19);
	// Convert to inverse Giga-years
	H0y = H0s * 3.156e16 ; // Gyr-1
	
	if (unit == 1) h_zz = H0y* sqrt(OM0*(1.+z)*(1.+z)*(1.+z) + OL0);       // Gyr-1
	if (unit == 2) h_zz = H0kmsmpc * sqrt(OM0*(1.+z)*(1.+z)*(1.+z) + OL0); // km/s/Mpc
	if (unit == 3) h_zz = H0s;
	
	return h_zz;
	
}

//********************************************************************************
double rho_crit(double z, double H0, double OM0)
/* rho_c = 3 H^2 / 8 pi G
   Returns result in Msol Mpc^-3, or h^2 Msol Mpc^-3 for H0=100
   z: redshift
   H0: Hubble constant at z=0
   OM0: Omega_m at z=0
 */
{
	double Hz, rhocrit, G;
	double h = H0/100.;
	
	Hz = H_z(z, OM0, H0, 3); // in s^-1
	G = 6.67384e-11; // m^3 kg^-1 s^-2
	
	rhocrit = 3.*Hz*Hz / (8.*pi*G);  // kg m^-3
	cout << "    rhocrit (h^2 kg m^-3): " << rhocrit / (h*h) << endl;
	
	rhocrit /= 1.99e30; // kg --> Msol
	rhocrit *= (3.08567758e22)*(3.08567758e22)*(3.08567758e22);  // m^-3 --> Mpc^-3
	
	return rhocrit;
}