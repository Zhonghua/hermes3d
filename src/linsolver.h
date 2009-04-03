// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2008 Miroslav Simko <msimko@miners.utep.edu>
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

#ifndef _LINSOLVER_H_
#define _LINSOLVER_H_

#include "common.h"
#include "matrix.h"

/// @defgroup linearsolvers Linear solvers
///
/// TODO: description

enum ESparseMatrixRepresentation {
	SMATRIX_ROW = 0,
	SMATRIX_COLUMN = 1,
	SMATRIX_OTHER
};

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp

/// Abstract class for encapsulation of a linear solver
///
///
/// TODO: Adjust interface to support faster update of matrix and rhs
///
/// @ingroup linearsolvers
class LinearSolver {
public:
	LinearSolver() {
		ndofs = 0;
		pages = NULL;
	}

	virtual ~LinearSolver() {
		delete [] pages;
	}

	/// prepare memory
	///
	/// @param[in] ndofs - number of unknowns
	virtual void prealloc(int ndofs) = 0;

	/// add indices of nonzero matrix element
	///
	/// @param[in] row  - row index
	/// @param[in] col  - column index
	virtual void pre_add_ij(int row, int col) = 0;

	/// allocate the memory for stiffness matrix and right-hand side
	virtual void alloc() = 0;

	/// free the memory associated with stiffness matrix and right-hand side
	virtual void free() = 0;

	/// update the stiffness matrix
	///
	/// @param[in] m    - the row where to update
	/// @param[in] n    - the column where to update
	/// @param[in] v    - value
	virtual void update_matrix(int m, int n, scalar v) = 0;

	/// update the stiffness matrix
	///
	/// @param[in] m         - number of rows of given block
	/// @param[in] n         - number of columns of given block
	/// @param[in] matrix    - block of values
	/// @param[in] rows      - array with row indexes
	/// @param[in] cols      - array with column indexes
	virtual void update_matrix(int m, int n, scalar **mat, int *rows, int *cols) = 0;

	/// update right-hand side
	///
	/// @param[in] idx - indices where to update
	/// @param[in] y   - value
	virtual void update_rhs(int idx, scalar y) = 0;

	/// update right-hand side
	///
	/// @param[in] n   - number of positions to update
	/// @param[in] idx - indices where to update
	/// @param[in] y   - values
	virtual void update_rhs(int n, int *idx, scalar *y) = 0;

	/// called before assembling
	///
	/// override this if you need to perform something special before the assembling
	virtual void begin_assembling() {}

	/// called after assembling
	///
	/// override this if you need to perform something special after the assembling
	virtual void finish_assembling() {}

	/// solve the system
	///
	/// @param[out] sln - array of doubles that will receive the solution (length has to be equal to ndofs)
	/// @return true if the system was solved successfully
	virtual bool solve_system(scalar *sln) = 0;

	/// dumping matrix and right-hand side
	///
	virtual bool dump_matrix(FILE *file, const char *var_name, EMatrixDumpFormat = DF_MATLAB_SPARSE) { return false; }
	virtual bool dump_rhs(FILE *file, const char *var_name, EMatrixDumpFormat = DF_MATLAB_SPARSE) { return false; }

	virtual int get_matrix_size() const { return 0; }

	// row or column representation of sparse matrix
	virtual ESparseMatrixRepresentation matrix_representation() = 0;

protected:
	static const int PAGE_SIZE = 62;

	struct Page {
		int count;
		int idx[PAGE_SIZE];
		Page *next;
	};

	int ndofs; // number of unknowns
	Page **pages;

	static int sort_and_store_indices(Page *page, int *buffer, int *max);
	static int get_num_indices(Page **pages, int ndofs);
	static void insert_value(int *Ai, scalar *Ax, int Alen, int idx, scalar value);
};

#endif
