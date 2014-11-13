#ifndef MV_H
#define MV_H

void MV_idealsurvey(int NN, double* xN, double* yN, double* zN, double* rN,
					double** w_in, char* idealxyz);

void MV_setparam(double& OmegaM, double& H0_in, double& sigma2_star_in);

void MV_multi_u(double** multi_u, int nu, int n_SNe, double** weight_in, double* vrad,
				double** ri, double* r);

void MV_W2_pq(double** W2_pq, int n_SNe, double** rihat, double* r, 
			  double** weight_in, int pkintsize, double* k_arr);

void MV_Re(double** Re_ij, int n_SNe, double** weight_in, double* verr);

void MV_Rv(double** Rv_ij, double** W2_pq, int n_SNe, int pkintsize, double dk, double* k_arr, double* ps0, double OM);

void MV_Rv_dk_arr(double** Rv_ij, double** W2_pq, int n_SNe, int pkintsize, double* dk, 
				  double* k_arr, double* ps0, double OM, double H_0);

void MV_Rpq(double** R_ij, double** G_nm, int n_SNe, double** weight_in, double* verr);

void MV_weights(double** weight_in, double* u, double** G_inv, double** Q_in, 
				double** lambda_ij, double** rihat, int n_SNe, double* vrad);

void MV_lambda_ij(double** lambda_ij, double** M_ij_inv, double** G_inv, 
				  double** Q_in, double** rihat, int n_SNe);

void MV_Mij(double** M_ij_inv, double** G_inv, int n_SNe, double** rihat);

void MV_Qin(double** Q_in, int n_SNe, int NN, double** vnvn, double** w_in);

void MV_vnvn(double** vnvn, int n_SNe, int NN, double** ri, double* r, double* xN, 
			 double* yN, double* zN, double* rN, int pkintsize, double dk, double* k_arr, double* ps0);

void MV_Gnm(double** G_nm, double** P_nm, double** G_inv, int n_SNe, double* r, double** rihat, 
			double* verr, int pkintsize, double dk, double* k_arr, double* ps0);

double MLE_CV_err(int ngal, double* rad, double** rihat, double* verr, double sigmastar2, int pkintsize,
				  double* k_arr, double* Pk_arr, double dk, double OM, double* uML, double* MLerrCV);

void MV_bulkflow(int n_SNe, double** ri, double* vrad,
				 double* verr, int NN, double* xN, double* yN, double* zN, double* rN,
				 double** w_in, int pkintsize, double dk, double* k_arr, double* ps0,
				 double* u, double* du, double* due, double* BFresults);

#endif