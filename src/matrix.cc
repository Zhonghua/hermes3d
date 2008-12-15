#include "config.h"
#include "common.h"
#include "matrix.h"

#include <common/error.h>


#define TINY 1.0e-20


// ludcmp, lubksb - LU decomposition and back-substitution routines from
// the book Numerical Recipes in C, adjusted to zero-based indexing

void ludcmp(double **a, int n, int *indx, double *d) {
	int i, imax = 0, j, k;
	double big, dum, sum, temp;
	double *vv = new double[n];
	MEM_CHECK(vv);

	*d = 1.0;
	for (i = 0; i < n; i++) {
		big=0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
		if (big == 0.0) { assert(false); EXIT(ERR_FAILURE, "Singular matrix in routine LUDCMP!"); }
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n-1) {
			dum = 1.0 / (a[j][j]);
			for (i = j+1; i < n; i++) a[i][j] *= dum;
		}
	}
	delete [] vv;
}


void lubksb(double **a, int n, int *indx, double *b) {
	int i, ip, j;
	double sum;

	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
		b[i] = sum;
	}
	for (i = n-1; i >= 0; i--) {
		sum = b[i];
		for (j = i+1; j < n; j++) sum -= a[i][j]*b[j];
		b[i] = sum / a[i][i];
	}
}


// choldc, cholsl - Cholesky decomposition and solution routines from
// the book Numerical Recipes in C, adjusted to zero-based indexing

void choldc(double **a, int n, double p[]) {
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			double sum = a[i][j];
			k = i;
			while (--k >= 0)
				sum -= a[i][k] * a[j][k];

			if (i == j) {
				if (sum <= 0.0) {
					assert(false);
					EXIT(ERR_FAILURE, "CHOLDC failed!");
				}
				else
					p[i] = sqrt(sum);
			}
			else
				a[j][i] = sum / p[i];
		}
	}
}

