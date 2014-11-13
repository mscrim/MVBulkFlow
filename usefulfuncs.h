/*
 *  Useful Functions
 *
 *  Morag I. Scrimgeour
 *  11/10/2014
 *
 */

#ifndef USEFULFUNCS_H
#define USEFULFUNCS_H

double truncate(double x, int n); 
double Tophat_WF_Fourier(double k, double R);
bool surveyis(const char* surveytype, const char* surveystring);
bool comp_strings(const char* string1, const char* string2);
double mean_arr(double* arr, int num);
void mean_sigma_arr(double* arr, int num, double* result);
double ang_betw_2_gals(double RA1, double Dec1, double RA2, double Dec2);
double radians(double);
double degrees(double);
int factorial(int);
double getDec(double x, double y, double z);
double getRA(double x, double y);
double* getxyz(double RA, double Dec, double r);
double LinearInterpolate(int n, double x[],double y[], double x0);
double LinearInterpolateRange(int n, double x[],double y[], double x0, double ymin, double ymax);
double Kroneckerdelta(int i, int j);
double absolute(double x1);
void galactic_conv(double& ragal, double& decgal, double& sgl, double& sgb, int optn);
void bprecess(double ra, double dec,double ra_1950, double dec_1950);
void galactic_conv_fromidl(double& ragal, double& decgal, double& sgl, double& sgb, int optn);
double j0(double);
double j1(double);
double j2(double);
void sort_arr(double* arr1, int num, int* indx1, int flag, bool replacearr);

#endif