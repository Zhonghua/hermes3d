// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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

#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

#include "../h3dconfig.h"
#include "../matrix.h"
#include "../linsolver.h"

#ifdef WITH_PETSC
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#endif

/// Wrapper of PETSc matrix, to store matrices used with PETSc in its native format
///
class PetscMatrix : public SparseMatrix {
public:
	PetscMatrix();
	virtual ~PetscMatrix();

	virtual void alloc();
	virtual void free();
	virtual void update(int m, int n, scalar v);
	virtual void update(int m, int n, scalar **mat, int *rows, int *cols);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
	virtual int get_matrix_size() const;

	// finish assembling of the matrix
	virtual void finish();

protected:
#ifdef WITH_PETSC
	Mat matrix;
#endif
	bool inited;

	friend class PetscLinearSolver;
};

/// Wrapper of PETSc vector, to store vectors used with PETSc in its native format
///
class PetscVector : public Vector {
public:
	PetscVector();
	virtual ~PetscVector();

	virtual void alloc(int ndofs);
	virtual void free();
	virtual void update(int idx, scalar y);
	virtual void update(int n, int *idx, scalar *y);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

	// finish assembling of the vector
	virtual void finish();

protected:
#ifdef WITH_PETSC
	Vec vec;
#endif
	bool inited;

	friend class PetscLinearSolver;
};

/// Encapsulation of PETSc linear solver
///
/// @ingroup linearsolvers
class PetscLinearSolver : public LinearSolver {
public:
	PetscLinearSolver(PetscMatrix &mat, PetscVector &rhs);
	virtual ~PetscLinearSolver();

	virtual bool solve();

protected:
	PetscMatrix &m;
	PetscVector &rhs;
};

#endif

