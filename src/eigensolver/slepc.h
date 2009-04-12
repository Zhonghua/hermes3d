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


#ifndef _SLEPC_H_
#define _SLEPC_H_

#include "../h3dconfig.h"
#include "../eigensolver.h"
#include "../linsolver/petsc.h"

#ifdef WITH_SLEPC
#include <slepceps.h>
#endif

class SlepcEigenSolver : public EigenSolver {
public:
	SlepcEigenSolver(PetscMatrix &a, PetscMatrix &b);
	virtual ~SlepcEigenSolver();

	/// Solve the eigenvalue problem
	///
	/// @return true if the problem was solved
	virtual bool solve();

	/// Return number of pairs that converged
	///
	/// @return number of eigenvalue/eigenvector pairs that converged
	int get_converged();

	/// Get eigenvalue/eigenvector
	///
	/// @param[in] j - index of eigen value/vector pair (0..nconv - 1, nconv is the value returned by get_converged)
	/// @param[out] kr - real part of the eigenvalue
	/// @param[out] ki - imaginary part of the eigenvalue
	/// @param[out] xr - real part of the eigenvector
	/// @param[out] xi - imaginary part of the eigenvector
	void get_eigen_pair(int j, scalar *kr, scalar *ki, PetscVector *xr, PetscVector *xi);

	/// Get relative error of the computed eigenvalue/eigenvector pair
	///
	/// @param[in] j - index of the eigenvalue/eigenvector pair (0..nconv - 1, nconv is the value returned by get_converged)
	/// @return the relative error of the j-th eigenvalue/eigenvector pair
	double compute_relative_error(int j);

protected:
#ifdef WITH_SLEPC
	EPS eps;			/// eigensolver context
#endif
	PetscMatrix &a;
	PetscMatrix &b;
};

#endif /* _SLEPC_H_ */
