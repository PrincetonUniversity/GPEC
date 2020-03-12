#include <complex.h>
#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

double complex xj = I;

double xintgrnd(int* neq, double* x, double* y, double* ydot, 
		int* energy_imaxis, double* ximag, int* xnutype, 
		double* energy_we, double* energy_nuk, double* energy_leff,
		double* xnufac, double* energy_wb, int* energy_n, 
		double* energy_wd, int* xf0type, double* energy_wn, 
		double* energy_wt, int* qt) {
	double complex cx, fx, nux, denom;
	double output;
	//Fortran logicals are true if the last bit is 1
	if (*energy_imaxis % 2 == 1) {
                output = 1.0;
		cx = xj * (*x);
	} else {
                output = 0.0;
		cx = (*x) + xj * (*ximag);
	}
	switch (*xnutype) {
		case 1: 
			nux = 0.0;
			break;
		case 2:
			nux = 0.00001 * (*energy_we);
			break;
		case 3:
			nux = (*energy_nuk);
			break;
		default:
			if (x == 0) {
				nux = DBL_MAX;
			} else {
				nux = (*energy_nuk)*(1+0.25*(*energy_leff)*(*energy_leff))*pow(cx,-1.5);
			}
	}
	nux = (*xnufac) * nux;

	denom = xj*((*energy_leff)*(*energy_wb)*sqrt(cx)+(double)(*energy_n)*((*energy_we)+(*energy_wd)*cx))-nux;

	switch (*xf0type) {
		case 0:
			fx = ((*energy_we)+(*energy_wn)+(*energy_wt)*(cx-1.5))* pow(cx,2.5)*exp(-cx) /denom;
			break;
		case 1:
			fx = ((*energy_we)+(*energy_wn)+(*energy_wt)*2)*pow(cx,2.5)*exp(-cx) /denom;
			break;
		default:
			fx = pow(cx,2.5)*exp(-cx)/((double)(*energy_n)) /1.0;
			break;
	}

	if (*qt % 2 == 1) {
		fx *= (cx-2.5);
	}

	ydot[0] = creal(fx);
	ydot[1] = cimag(fx);
	return output;
}
