// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
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

#include "../h3dconfig.h"
#include "slepc.h"

SlepcEigenSolver::SlepcEigenSolver(PetscMatrix &a)
	: a(a), b(a)
{
#ifdef WITH_SLEPC
	EPSCreate(PETSC_COMM_WORLD, &eps);
	EPSSetOperators(eps, a.matrix, PETSC_NULL);
	EPSSetProblemType(eps, EPS_NHEP);
#endif
}

SlepcEigenSolver::SlepcEigenSolver(PetscMatrix &a, PetscMatrix &b)
	: a(a), b(b)
{
#ifdef WITH_SLEPC
	EPSCreate(PETSC_COMM_WORLD, &eps);
	EPSSetOperators(eps, a.matrix, b.matrix);
#endif
}

SlepcEigenSolver::~SlepcEigenSolver() {
#ifdef WITH_SLEPC
	EPSDestroy(eps);
#endif
}

bool SlepcEigenSolver::solve() {
#ifdef WITH_SLEPC
	// default initialization
	PetscErrorCode ec;
	ec = EPSSolve(eps);
	return !ec;
#else
	return false;
#endif
}

int SlepcEigenSolver::get_converged() {
#ifdef WITH_SLEPC
	int nconv;
	EPSGetConverged(eps, &nconv);
	return nconv;
#else
	return 0;
#endif
}

void SlepcEigenSolver::get_eigen_pair(int j, scalar *kr, scalar *ki, scalar *xr, scalar *xi) {
#ifdef WITH_SLEPC
	Vec vr, vi;
	MatGetVecs(a.matrix, PETSC_NULL, &vr);
	MatGetVecs(a.matrix, PETSC_NULL, &vi);

	EPSGetEigenpair(eps, (PetscInt) j, (PetscScalar *) kr, (PetscScalar *) ki, vr, vi);

	// get the size of the eigenvector
	int n;
	VecGetSize(vr, &n);
	scalar *a;

	xr[0] = 1.0;				// dirichlet DOF
	VecGetArray(vr, &a);
	memcpy(xr + 1, a, n * sizeof(scalar));
	VecRestoreArray(vr, &a);

	xi[0] = 1.0;				// dirichlet DOF
	VecGetArray(vi, &a);
	memcpy(xi + 1, a, n * sizeof(scalar));
	VecRestoreArray(vi, &a);
#endif
}

double SlepcEigenSolver::compute_relative_error(int j) {
#ifdef WITH_SLEPC
	PetscReal err = 1.0;
	EPSComputeRelativeError(eps, j, &err);
	return (double) err;
#else
	return 0.0;
#endif
}
