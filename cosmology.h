#ifndef COSMOLOGY_H
#define COSMOLOGY_H

double dxx_mpch(double, double, double);
double* dxx_mpch_arr(int num, double* zz, double OM, double OL);
void readin_CAMB(char* camb_file, double kmaxx, double dk, int& pkintsize, double *& k_arr, double *& Pk_arr);
double rho_crit(double z, double h, double OM0);
double H_z(double, double, double, int);

#endif