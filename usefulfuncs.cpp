/*
 *  Useful C++ Functions
 *
 *  Morag Scrimgeour
 *  11 Oct 2014
 *
 */

#ifndef pi
#define pi 3.14159265358979323846264338328 
#endif

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

//********************************************************************************
double truncate(double x, int n) 
// Truncate double x to n decimal places
{
	int indx = pow(10.,(double)n);
	double x2 = floorf(x * indx) / indx; 
	
	return x2;
}
//********************************************************************************
double Tophat_WF_Fourier(double k, double R)
// Tophat window function in Fourier Space, widetilde{W}_{\rm TK}(k,R)
{
	double W_kr, kR;
	
	kR = k*R;
	
	W_kr = 3./kR * j1(kR);  // j_1 = (sin(x)-x*cos(x))/(x*x);
	
	if (R == 0.) W_kr = 1.;
	
	return W_kr;
}

//********************************************************************************
bool surveyis(const char* surveytype, const char* surveystring)
{
	bool equal;
	equal = false;
	
	if (strcmp(surveytype, surveystring) == 0) equal = true;
	
	return equal;
}
//********************************************************************************
bool comp_strings(const char* string1, const char* string2)
{
	bool equal;
	equal = false;
	
	if (strcmp(string1, string2) == 0) equal = true;
	
	return equal;
}
//********************************************************************************
double mean_arr(double* arr, int num)
// Input: arr = array to find mean for
//        num = size of arr
{
	double meanval;
	double sumval = 0.;
	
	for (int i=0; i<num; i++){
		sumval += arr[i];
	}
	meanval = sumval / ((double)num);
	
	return meanval;
}
//********************************************************************************
void mean_sigma_arr(double* arr, int num, double* result)
// Input: arr = array to find standard deviation for
//        num = size of arr
{
	double meanval, var, sig;
	//double* result = new double[2];
	
	meanval = mean_arr(arr, num);
	
	var=0.;
	for (int i=0; i<num; i++){
		var += (arr[i] - meanval)*(arr[i] - meanval);
	}
	var = var / ((double)num-1.);  // Watch out for N-1 factor!
	
	sig = sqrt(var);
	
	result[0] = meanval;
	result[1] = sig;
	
	return;
}
//********************************************************************************
double ang_betw_2_gals(double RA1, double Dec1, double RA2, double Dec2)
// Returns theta in degrees. RA, Dec input in degrees
{
	double theta, costheta;
	double x1, y1, z1, x2, y2, z2;
	
	x1 = sin((90.-Dec1)*pi/180.) * cos(RA1*pi/180.);
	y1 = sin((90.-Dec1)*pi/180.) * sin(RA1*pi/180.);
	z1 = cos((90.-Dec1)*pi/180.);
	
	x2 = sin((90.-Dec2)*pi/180.) * cos(RA2*pi/180.);
	y2 = sin((90.-Dec2)*pi/180.) * sin(RA2*pi/180.);
	z2 = cos((90.-Dec2)*pi/180.);
	
	// Calculate theta via dot product
	costheta = (x1*x2 + y1*y2 + z1*z2);
	theta = acos(costheta)*180./pi;
	
	return theta;
}

//********************************************************************************
// Function radians - convert degrees to radians
double radians(double deg)
{
	double rad = deg * pi / 180.;
	return rad;
}
//********************************************************************************
// Function degrees - convert radians to degrees
double degrees(double rad)
{
	double deg = rad * 180. / pi;
	return deg;
}

//********************************************************************************
// Function FACTORIAL - returns n!
int factorial(int n)
{
	int i=1,fact=1;
	
	while (i <= n){
		fact *= i;
		i++;
	}
	
	return fact;
}

//********************************************************************************
// Function galactic_conv - Convert between celestial and Galactic (or Supergalactic) coordinates
void galactic_conv(double& ra, double& dec, double& gl, double& gb, int optn)
//   optn = 1: RA, Dec --> gl, gb
//   optn = 2: gl, gb --> RA, Dec
//   All in DEGREES, J2000
//
// Convert Galactic to Equatorial coordinates (J2000.0)
// (use at own risk)
// Input: [l,b] in decimal degrees
// Returns: [ra,dec] in decimal degrees
// Source:
// - Book: "Practical astronomy with your calculator" (Peter Duffett-Smith)
// - Wikipedia "Galactic coordinates"
// Tests (examples given on the Wikipedia page):
// >>> ga2equ([0.0, 0.0]).round(3)
// array([ 266.405, -28.936])
// >>> ga2equ([359.9443056, -0.0461944444]).round(3)
// array([ 266.417, -29.008])
{
	double pole_ra, pole_dec, posangle;
	double l, b, ragal, decgal;
	double centre_ra, centre_dec;
	double valj, valk, valq;
	
	// North galactic pole (J2000) -- according to Wikipedia
	pole_ra = radians(192.859508);
	pole_dec = radians(27.128336);
	posangle = radians(122.932-90.0);
	
	// Galactic Centre (J2000)  -- according to Galactic pdf by M Bastion, Macquarie University
	centre_ra = radians(266.4);
	centre_dec = radians(-28.929656275);
	
	// North galactic pole (B1950)
	//pole_ra = radians(192.25)
	//pole_dec = radians(27.4)
	//posangle = radians(123.0-90.0)
	
	if (optn == 1){
		
		ragal = radians(ra);
		decgal = radians(dec);
		
		//centre_dec = atan2((-1.*cos(centre_ra-pole_ra)), tan(pole_dec));
		
		b = asin( sin(decgal)*sin(pole_dec) + cos(decgal)*cos(pole_dec)*cos(ragal-pole_ra) );
		
		valj = ( sin(decgal)*cos(pole_dec)-cos(decgal)*sin(pole_dec)*cos(ragal-pole_ra) ) / cos(b);
		valk = asin( cos(decgal)*sin(ragal-pole_ra) / cos(b));
		valq = acos(sin(centre_dec)/cos(pole_dec));
		
		if (valj < 0) l = valq + valk - pi;
		else l = valq - valk;
		
		if (l < 0) l += 2.*pi;
		
		gb = degrees(b);
		gl = degrees(l);
		
	}
	else if (optn == 2){
		// Convert l,b to radians
		l = radians(gl);
		b = radians(gb);
		
		ra = atan2((cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra;
		dec = asin(cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec));
		
		ra = degrees(ra);
		dec = degrees(dec);
		
	}
	else {
		cout << "Option not available" << flush << endl;
		exit(1);
	}
	
	return;
}

//********************************************************************************
// Function bprecess - Precess positions from J2000.0 (FK5) to B1950.0 (FK4)
void bprecess(double ra, double dec,double ra_1950, double dec_1950)
//
// PURPOSE:
//       Precess positions from J2000.0 (FK5) to B1950.0 (FK4)
// EXPLANATION:
//       Calculates the mean place of a star at B1950.0 on the FK4 system from
//       the mean place at J2000.0 on the FK5 system.
// RA,DEC - Input J2000 right ascension and declination in *degrees*.
{
	double rad_vel;
	double parallax;
	double epoch = 2000.0;
	
	double radeg = 180./pi;
	double sec_to_radian = 1.0/radeg/3600.0;
	
	double M[6][6];
	M[0][0] = +0.9999256795; M[0][1] = -0.0111814828; M[0][2] = -0.0048590040; M[0][3] = -0.000551; M[0][4] = -0.238560; M[0][5] = +0.435730;
	M[1][0] = +0.0111814828; M[1][1] = +0.9999374849; M[1][2] = -0.0000271557; M[1][3] = +0.238509; M[1][4] = -0.002667; M[1][5] = -0.008541;
	M[2][0] = +0.0048590039; M[2][1] = -0.0000271771; M[2][2] = +0.9999881946; M[2][3] = -0.435614; M[2][4] = +0.012254; M[2][5] = +0.002117;
	M[3][0] = -0.00000242389840; M[3][1] = +0.00000002710544; M[3][2] = +0.00000001177742; M[3][3] = +0.99990432; M[3][4] = -0.01118145; M[3][5] = -0.00485852;
	M[4][0] = -0.00000002710544; M[4][1] = -0.00000242392702; M[4][2] = +0.00000000006585; M[4][3] = +0.01118145; M[4][4] = +0.99991613; M[4][5] = -0.00002716;
	M[5][0] = -0.00000001177742; M[5][1] = +0.00000000006585; M[5][2] = -0.00000242404995; M[5][3] = +0.00485852; M[5][4] = -0.00002717; M[5][5] = +0.99996684;
	
	double A_dot[3];
	A_dot[0] = 1e-3*1.244;           //in arc seconds per century
	A_dot[1] = 1e-3*-1.579;
	A_dot[2] = 1e-3*-0.660;
	
	
	double ra_rad, dec_rad, cosra, sinra, cosdec, sindec;
	
	ra_rad = ra/radeg;
	dec_rad = dec/radeg;
	cosra =  cos( ra_rad );
	sinra = sin( ra_rad );
	cosdec = cos( dec_rad );
	sindec = sin( dec_rad );
	
	dec_1950 = dec*0.;
	ra_1950 = ra*0.;
	
	// Following statement moved inside loop in Feb 2000.
	double A[3];
	A[0] = 1e-6*-1.62557;
	A[1] = 1e-6*-0.31919;
	A[2] =  1e-6*-0.13843;        //in radians
	
	double r0[3];
	r0[0] = cosra*cosdec;
	r0[1] = sinra*cosdec;
	r0[2] = sindec;
	
	double r0_dot[3];
	r0_dot[0] = 0.0;
	r0_dot[1] = 0.0;
	r0_dot[2] = 0.0;
	
	
	double R_0[6];
	for (int i=0; i<3; i++){
		R_0[i] = r0[i];
		R_0[i+3] = r0_dot[i];
	}
	
	double R_1[6];
	// Multiply columns of M by rows of R_0
	for (int i=0; i<6; i++){
		R_1[i] = 0.;
		for (int j=0; j<6; j++){
			R_1[i] +=  M[i][j]*R_0[j];
		}
	}
	
	// Include the effects of the E-terms of aberration to form r and r_dot.
	
	double r1[3], r1_dot[3];
	for (int i=0; i<3; i++){
		r1[i] = R_1[i];
		r1_dot[i] = R_1[i+3];
	}
	
	for (int i=0; i<3; i++){
		r1[i] = r1[i] + sec_to_radian * r1_dot[i] * (epoch - 1950.0)/100.;
		A[i] = A[i] + sec_to_radian * A_dot[i] * (epoch - 1950.0)/100.;
	}
	
	double x1, y1, z1, rmag;
	x1 = R_1[0]; y1 = R_1[1]; z1 = R_1[2];
	rmag = sqrt( x1*x1 + y1*y1 + z1*z1 );
	
	double s1[3], s1_dot[3], s[3];
	for (int i=0; i<3; i++){
		s1[i] = r1[i]/rmag;
		s1_dot[i] = r1_dot[i]/rmag;
		s[i] = s1[i];
	}
	
	
	double totsA;
	double r[3];
	
	for (int j=0; j<3; j++){
		// calculate total (s*A)
		totsA = 0;
		for (int i=0; i<3; i++){
			totsA += s[i]*A[i];
		}
		// r = s1 + A - (total(s * A))*s
		for (int i=0; i<3; i++){
			r[i] = s1[i] + A[i] - totsA*s[i];
			s[i] = r[i]/rmag;
		}
	}
	
	double x,y,z, r2;
	x = r[0];
	y = r[1];
	z = r[2];
	r2 = x*x + y*y + z*z;
	rmag = sqrt(r2);
	
	dec_1950 = asin( z / rmag);
	ra_1950 = atan2( y, x);
	
	if (ra_1950 < 0.){
		ra_1950 = ra_1950 + 2.*pi;
	}
	
	ra_1950 = ra_1950*radeg;
	dec_1950 = dec_1950*radeg;
	
	return;
}
//********************************************************************************
// Function galactic_conv - Convert between celestial and Galactic (or Supergalactic) coordinates
void galactic_conv_fromidl(double ragal, double decgal, double sgl, double sgb, int optn)
//
//	ragal, decgal : degrees
// 
//	Convert between celestial and Galactic (or Supergalactic) coordinates
//	
//	Program to convert right ascension (ra) and declination (dec) to 
//	Galactic longitude (gl) and latitude (gb) 
//	
//	optn = 1 for galactic to RA, Dec
//	optn = 2 for supergalactic to RA, Dec
//	optn = 3 for RA, Dec to galactic
//	optn = 4 for RA, Dec to supergalactic
{
	
	// Assumes J2000
	double ras, decs, sdp, cdp;
	double sdec, cdec, singb, cosgb, sine, cose;
	double radeg = 180.0/pi; double radhrs = radeg/15.0;
	double rapol, decpol, dlon;
	double ra2, dec2;
	double gl,gb;
	double cgb;
	
	//GALACTIC
	//if (optn == 1 || optn == 3) {rapol = 12.82; decpol = 27.40; dlon = 123.0;}  // Christina's code
	if (optn == 1 || optn == 3) {rapol = 12.0 + 49.0/60.0; decpol = 27.40; dlon = 123.0;}
	// SUPERGALACTIC
	//if (optn == 2 || optn == 4) {rapol = 18.917; decpol = 15.709; dlon = 26.73153707;}  // Christina's code
	if (optn == 2 || optn == 4) {rapol = 283.18940711/15.0; decpol = 15.64407736; dlon = 26.73153707;}
	
	sdp = sin(decpol/radeg);
	cdp = sqrt(1.0-sdp*sdp);
	
	
	// RA, Dec to galactic / supergalactic
	if (optn == 3 || optn == 4) {
		ras = ragal;
		decs = decgal;
		
		bprecess(ras,decs,ra2,dec2);  // FIX UP
		ras = ra2;
		decs = dec2;
		
		ras = ras/radeg - rapol/radhrs; 
        sdec = sin(decs/radeg);
        cdec = sqrt(1.0-sdec*sdec);
        sgb = sdec*sdp + cdec*cdp*cos(ras);
        gb = radeg * asin(sgb);
        cgb = sqrt(1.0-sgb*sgb);
        sine = cdec * sin(ras) / cgb;
        cose = (sdec-sdp*sgb) / (cdp*cgb);
        gl = dlon - radeg*atan2(sine,cose);
		if (gl < 0.) gl = gl + 360.0;
		
		//ras = ras/raddeg - rapol/radhrs; 
		//sdec = sin(decs/raddeg);
		//cdec = sqrt(1.0-sdec*sdec);
		//singb = sdec*sdp + cdec*cdp*cos(ras);
		//sgb = raddeg*asin(singb);
		
		//cosgb = sqrt(1.0-singb*singb);
		//sine = cdec*sin(ras) / cosgb;
		//cose = (sdec-sdp*singb) / (cdp*cosgb);
		//sgl = dlon - raddeg*atan2(sine,cose);
		//if (sgl < 0) { sgl = sgl + 360.0; }
	}
	else {
		cout << "Option not yet available" << endl;
		exit(1);
	}
	
	return;	
}

//********************************************************************************
// getDec and getRA
double getDec(double x, double y, double z)
{
	// calculate Declination in degrees from Cartesian coords in [Mpc/h]
	double Dec;
	double r;
	
	r = sqrt(x*x + y*y + z*z);
	
	Dec = 90. - acos(z/r) * 180. / pi;
	
	return Dec;
}
double getRA(double x, double y)
{
	// calculate RA in degrees from Cartesian coords in [Mpc/h]
	double RA;
	
	RA = atan2(y,x) * 180. / pi;
	if (RA < 0.) RA += 360.;
	
	return RA;
}
double* getxyz(double RA, double Dec, double r)
{
	double* result = new double[3];
	double x, y, z;
	double theta, phi;
	
	phi = RA*pi/180.;
	theta = (90. - Dec)*pi/180.;
	
	x = r*sin(theta)*cos(phi);
	y = r*sin(theta)*sin(phi);
	z = r*cos(theta);
	
	result[0] = x;
	result[1] = y;
	result[2] = z;
	
	return result;
}

//********************************************************************************
double LinearInterpolate(int n, double x[],double y[], double x0)
{
	double x1,x2,y1,y2,y0;
	int check;
	
	for (int i=0; i<n; i++) {
		if (x[i] <= x0 & x[i+1] > x0) {
			x1 = x[i];
			x2 = x[i+1];
			y1 = y[i];
			y2 = y[i+1];
			check = 3;
		}
	}
	if (check != 3) {
		cout << "Could not interpolate" << endl;
		cout << "    x0 = " << x0 << " , x range: [" << x[0] << " , " << x[n-1] << "]" << endl;
		cout << "    n = " << n << " , y range: [" << y[0] << " , " << y[n-1] << "]" << flush << endl;
		exit( -1 );
	}
	
	y0 = y1 + ((y2-y1)/(x2-x1))*(x0-x1);
	
	return y0;
}
//********************************************************************************
double LinearInterpolateRange(int n, double x[],double y[], double x0, double ymin, double ymax)
// Linear interpolate array, but only between ymin and ymax
// (Can use for a function that has >1 intercept)
{
	double x1,x2,y1,y2,y0;
	int check;
	double minx=x[0],maxx=x[n-1];
	
	
	// For monotonically increasing
	for (int i=0; i<n; i++) {
		if (y[i] > ymin && y[i] < ymax){
			//cout << i << "  " << x[i] << "  " << y[i] << endl;
			if (x[i] < minx) minx = x[i];
			if (x[i] > maxx) maxx = x[i];
			if (x[i] <= x0 & x[i+1] > x0) {
				x1 = x[i];
				x2 = x[i+1];
				y1 = y[i];
				y2 = y[i+1];
				check = 3;
			}
		}
	}
	// For monotonically decreasing
	if (check != 3) {
		for (int i=0; i<n; i++) {
			if (y[i] > ymin && y[i] < ymax){
				if (x[i] > x0 & x[i+1] <= x0) {
					x1 = x[i];
					x2 = x[i+1];
					y1 = y[i];
					y2 = y[i+1];
					check = 3;
				}
			}
		}
	}	
	
	if (check != 3) {
		cout << "Could not interpolate in range [" << ymin << "," << ymax << "]" << endl;
		cout << "    x0 = " << x0 << " , array range: [" << minx << " , " << maxx << "]" << flush << endl;
		exit( -1 );
	}
	
	y0 = y1 + ((y2-y1)/(x2-x1))*(x0-x1);
	
	return y0;
}

//********************************************************************************
// delta(i,j) - returns 1 if i=j and 0 if i != j
double Kroneckerdelta(int i, int j)
{
	int answer;
	if (i == j) answer = 1;
	else answer = 0;
	
	return answer;
}

//********************************************************************************
// Absolute value function
double absolute(double x1)
{	
	double x_abs;
	
	if (x1 >= 0) x_abs = x1;
	else x_abs = -x1;
	
	return x_abs;
}

//********************************************************************************
// j0(x) - returns the 0th-order spherical Bessel function as a function of x (x in radians)
double j0(double x)
{
	double j_0 = sin(x)/x;
	
	return j_0;
}
//********************************************************************************
// j1(x) - returns the 0th-order spherical Bessel function as a function of x (x in radians)
double j1(double x)
{
	double j_1 = (sin(x)-x*cos(x))/(x*x);
	
	return j_1;
}
//********************************************************************************
// j2(x) - returns the 2nd-order spherical Bessel function as a function of x (x in radians)
double j2(double x)
{
	double j_2 = (3./(x*x) - 1.)*sin(x)/x - 3*cos(x)/(x*x);
	
	return j_2;
}

//********************************************************************************
//double rround(double x, int n)
//{	
//	int y1;
//	double y2;

//	y1 = x * pow(10,n) +0.5;

//	y2 = (double)y1 / pow(10,n);

//	return y2;
//}


//********************************************************************************
// SORT ARRAY
typedef std::pair<double,int> indxpair;
// Set up comparator to use in std::sort
bool comparator_asc( const indxpair& l, const indxpair& r)
{
	return l.first < r.first;
}
bool comparator_desc( const indxpair& l, const indxpair& r)
{
	return l.first > r.first;
}
void sort_arr(double* arr1, int num, int* indx1, int flag, bool replacearr)
/*
 * arr1 : array to sort
 * num : size of array
 * indx1: array to hold sorted indexes
 * Flag = 2: Sort descending
 * Flag = 1: Sort ascending
 *
 10... done in 0 s
 100... done in 0 s
 1000... done in 0 s
 10000... done in 0 s
 100000... done in 0 s
 1000000... done in 1 s
 5647390... done in 5 s
 9474168... done in 11 s
 */
{
	vector<indxpair> pairArr (num);
	
	for (int i=0; i<num; i++){
		pairArr[i] = make_pair(arr1[i],i);
	}
	
	// Sort pairArr
	if (flag == 1) std::sort (pairArr.begin(), pairArr.end(), comparator_asc);
	if (flag == 2) std::sort (pairArr.begin(), pairArr.end(), comparator_desc);
	
	for (int i=0; i<num; i++){
		if (replacearr) arr1[i] = pairArr[i].first;
		indx1[i] = pairArr[i].second;
	}
	
	return;		
}
