// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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


#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>
#include <common/array.h>

double A[3][3] = {
	{ -1,   6, -12 },
	{  0, -13,  30 },
	{  0,  -9,  20 }
};

double B[3][3] = {
	{  4,   2,   1 },
	{  2,   1,   3 },
	{  3,   3,   1 }
};


double S[3][3] = {
	{  7,   2,  -4 },
	{ 10,   3,  -6 },
	{  6,   2,  -3 }
};

double M[3][3] = {
	{ -5,  -1,   3 },
	{ 12,   3,  -7 },
	{ -3,  -1,   2 }
};

// helpers ////////////////////////////////////////////////////////////////////

bool testPrint(bool value, const char *msg, bool correct) {
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

void build_matrix(SparseMatrix *a) {
	// matrix
	a->prealloc(3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			a->pre_add_ij(i, j);
	a->alloc();

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			a->update(i, j, A[i][j]);
}

void build_matrix(SparseMatrix *s, SparseMatrix *m) {
	// stiffness matrix
	s->prealloc(3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			s->pre_add_ij(i, j);
	s->alloc();

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			s->update(i, j, A[i][j]);

	// mass matrix
	m->prealloc(3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m->pre_add_ij(i, j);
	m->alloc();

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m->update(i, j, B[i][j]);
}

void results(SlepcEigenSolver *solver, int n) {
	int nc = solver->get_converged();
	for (int i = 0; i < nc; i++) {
		double kr, ki;							// real and imaginary part of the eigenvalue
		double *xr = new double [n + 1];	// real and imaginary part of the eigenvector
		double *xi = new double [n + 1];

		printf("%d-th pair\n", i);

		solver->get_eigen_pair(i, &kr, &ki, xr, xi);

		printf("  Eigenvalue  : % lf + % lfi\n", kr, ki);
		printf("  Eigenvector : ");
		for (int j = 1; j <= n; j++) {
			if (j > 1) printf(", ");
			printf("% lf", xr[j]);
		}
		printf("\n");

		double err = solver->compute_relative_error(i);
		printf("  Rel. error  : % e\n", err);

		delete [] xr;
		delete [] xi;
	}
}

int main(int argc, char *argv[]) {
	int ret = ERR_SUCCESS;

#ifdef WITH_SLEPC
	SlepcInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
#endif

	if (strcasecmp(argv[1], "slepc") == 0) {
#ifdef WITH_SLEPC
		PetscMatrix a;
		build_matrix(&a);
		a.finish();

		SlepcEigenSolver solver(a);
		solver.solve();
		results(&solver, 3);
#endif
	}
	else if (strcasecmp(argv[1], "slepc-gen") == 0) {
#ifdef WITH_SLEPC
		PetscMatrix a, b;
		build_matrix(&a, &b);
		a.finish();
		b.finish();

		SlepcEigenSolver solver(a, b);
		solver.solve();
		results(&solver, 3);
#endif
	}
	else
		ret = ERR_FAILURE;

#ifdef WITH_PETSC
	SlepcFinalize();
#endif

	return ret;
}

