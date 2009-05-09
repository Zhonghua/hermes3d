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

#ifndef _UMFPACK_SOLVER_H_
#define _UMFPACK_SOLVER_H_

#include "../linsolver.h"
#include "../matrix.h"

class UMFPackMatrix : public SparseMatrix {
public:
	UMFPackMatrix();
	virtual ~UMFPackMatrix();

	virtual void alloc();
	virtual void free();
	virtual void update(int m, int n, scalar v);
	virtual void update(int m, int n, scalar **mat, int *rows, int *cols);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
	virtual int get_matrix_size() const;

protected:
	// UMFPack specific data structures for storing matrix, rhs
	int *Ap;
	int *Ai;
	scalar *Ax;

	static void insert_value(int *Ai, scalar *Ax, int Alen, int idx, scalar value);

	friend class UMFPackLinearSolver;
};

class UMFPackVector : public Vector {
public:
	UMFPackVector();
	virtual ~UMFPackVector();

	virtual void alloc(int ndofs);
	virtual void free();
	virtual void update(int idx, scalar y);
	virtual void update(int n, int *idx, scalar *y);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

//protected:
	scalar *v;

	friend class UMFPackLinearSolver;
};


/// Encapsulation of UMFPACK linear solver
///
/// @ingroup linearsolvers
class UMFPackLinearSolver : public LinearSolver {
public:
	UMFPackLinearSolver(UMFPackMatrix &m, UMFPackVector &rhs);
	virtual ~UMFPackLinearSolver();

	virtual bool solve();

protected:
	UMFPackMatrix &m;
	UMFPackVector &rhs;
};

#endif
