// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 Pavel Kus <pavel.kus@gmail.com>
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef __MATRIX_H
#define __MATRIX_H

#include "common.h"
#include <common/error.h>


/// Creates a new (full) matrix with m rows and n columns with entries of the type T.
/// The entries can be accessed by matrix[i][j]. To delete the matrix, just
/// do "delete matrix".
template<typename T>
T **new_matrix(int m, int n = 0) {
	if (!n) n = m;
	T **vec = (T **) new char[sizeof(T *) * m + sizeof(T) * m * n];
	MEM_CHECK(vec);
	memset(vec, 0, sizeof(T *) * m + sizeof(T) * m * n);
	T *row = (T *) (vec + m);
	for (int i = 0; i < m; i++, row += n)
		vec[i] = row;
	return vec;
}

/// Given a matrix a[n][n], this routine replaces it by the LU decomposition of a rowwise
/// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
/// indx[n] is an output vector that records the row permutation effected by the partial
/// pivoting; d is output as +-1 depending on whether the number of row interchanges was even
/// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
/// or invert a matrix.
void ludcmp(double **a, int n, int *indx, double *d);

/// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
/// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
/// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
/// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
/// and can be left in place for successive calls with different right-hand sides b. This routine takes
/// into account the possibility that b will begin with many zero elements, so it is efficient for use
/// in matrix inversion.
void lubksb(double **a, int n, int* indx, double *b);


/// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
/// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
/// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
/// elements which are returned in p[n].
void choldc(double **a, int n, double p[]);

/// Solves the set of n linear equations A*x = b, where a is a positive-definite symmetric matrix.
/// a[n][n] and p[n] are input as the output of the routine choldc. Only the lower
/// subdiagonal portion of a is accessed. b[n] is input as the right-hand side vector. The
/// solution vector is returned in x[n]. a, n, and p are not modified and can be left in place
/// for successive calls with different right-hand sides b. b is not modified unless you identify b and
/// x in the calling sequence, which is allowed. The right-hand side b can be complex, in which case
/// the solution x is also complex.
template<typename T>
void cholsl(double **a, int n, double p[], T b[], T x[]) {
	int i, k;
	T sum;

	for (i = 0; i < n; i++) {
		sum = b[i];
		k = i;
		while (--k >= 0)
			sum -= a[i][k] * x[k];
		x[i] = sum / p[i];
	}

	for (i = n-1; i >= 0; i--) {
		sum = x[i];
		k = i;
		while (++k < n)
			sum -= a[k][i] * x[k];
		x[i] = sum / p[i];
	}
}


#endif
