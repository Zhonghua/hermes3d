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


#include "petsc.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/utils.h>

PetscMatrix::PetscMatrix() {
	inited = false;
}

PetscMatrix::~PetscMatrix() {
	free();
}

void PetscMatrix::alloc() {
#ifdef WITH_PETSC
	assert(pages != NULL);

	// calc nnz
	int *nnz = new int[this->ndofs];
	MEM_CHECK(nnz);

	// fill in nnz
	int aisize = get_num_indices();
	int *ai = new int[aisize];
	MEM_CHECK(ai);

	// sort the indices and remove duplicities, insert into ai
	int pos = 0;
	for (int i = 0; i < ndofs; i++) {
		nnz[i] = sort_and_store_indices(pages[i], ai + pos, ai + aisize);
		pos += nnz[i];
	}
	delete [] pages; pages = NULL;
	delete [] ai;

	//
	MatCreateSeqAIJ(PETSC_COMM_SELF, this->ndofs, this->ndofs, 0, nnz, &matrix);
//	MatSetOption(matrix, MAT_ROW_ORIENTED);
//	MatSetOption(matrix, MAT_ROWS_SORTED);

	delete [] nnz;

	inited = true;
#endif
}

void PetscMatrix::free() {
#ifdef WITH_PETSC
	if (inited) MatDestroy(matrix);
	inited = false;
#endif
}

void PetscMatrix::update(int m, int n, scalar v) {
#ifdef WITH_PETSC
	MatSetValues(matrix, 1, &m, 1, &n, (PetscScalar *) &v, ADD_VALUES);
#endif
}

void PetscMatrix::update(int m, int n, scalar **mat, int *rows, int *cols) {
#ifdef WITH_PETSC
	// TODO: pass in just the block of the matrix without DIRICHLET_DOFs (so that can use MatSetValues directly without checking
	// row and cols for -1)
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)	 {		// cols
			if (mat[i][j] != 0.0 && rows[i] != DIRICHLET_DOF && cols[j] != DIRICHLET_DOF) {		// ignore "dirichlet DOF"
				MatSetValues(matrix, 1, rows + i, 1, cols + j, (PetscScalar *) &(mat[i][j]), ADD_VALUES);
			}
		}
#endif
}

bool PetscMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
#ifdef WITH_PETSC
	switch (fmt) {
		case DF_MATLAB_SPARSE:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		case DF_HERMES_BIN: {
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;
		}

		case DF_PLAIN_ASCII:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
#endif
	return false;
}

int PetscMatrix::get_matrix_size() const {
	return 0;
}

void PetscMatrix::finish() {
#ifdef WITH_PETSC
	MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
#endif
}

// PETSc vector //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscVector::PetscVector() {
#ifdef WITH_PETSC
	inited = false;
#endif
}

PetscVector::~PetscVector() {
	free();
}

void PetscVector::alloc(int n) {
#ifdef WITH_PETSC
	free();
	ndofs = n;
	VecCreateSeq(PETSC_COMM_SELF, ndofs, &vec);
	inited = true;
#endif
}

void PetscVector::free() {
#ifdef WITH_PETSC
	if (inited) VecDestroy(vec);
	inited = false;
#endif
}

void PetscVector::update(int idx, scalar y) {
#ifdef WITH_PETSC
	if (idx >= 0) VecSetValues(vec, 1, &idx, (PetscScalar *) &y, ADD_VALUES);
#endif
}

void PetscVector::update(int n, int *idx, scalar *y) {
#ifdef WITH_PETSC
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0) VecSetValues(vec, 1, idx + i, (PetscScalar *) (y + i), ADD_VALUES);
#endif
}

bool PetscVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
#ifdef WITH_PETSC
	PetscScalar *a;

	switch (fmt) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx1\n%s = [\n", this->ndofs, var_name);
			a = get_vector();
			for (int i = 0; i < this->ndofs; i++)
				fprintf(file, SCALAR_FMT "\n", SCALAR(a[i]));
			fprintf(file, " ];\n");
			restore(a);
			return true;

		case DF_HERMES_BIN: {
			hermes_fwrite("H2DR\001\000\000\000", 1, 8, file);
			a = get_vector();
			int ssize = sizeof(scalar);
			hermes_fwrite(&ssize, sizeof(int), 1, file);
			hermes_fwrite(&ndofs, sizeof(int), 1, file);
			hermes_fwrite(a, sizeof(scalar), ndofs, file);
			restore(a);
			return true;
		}

		case DF_PLAIN_ASCII:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
#endif
	return false;
}

void PetscVector::finish() {
#ifdef WITH_PETSC
	VecAssemblyBegin(vec);
	VecAssemblyEnd(vec);
#endif
}

scalar *PetscVector::get_vector() {
#ifdef WITH_PETSC
	scalar *a;
	VecGetArray(vec, &a);
	return a;
#else
	return NULL;
#endif
}

void PetscVector::restore(scalar *v) {
#ifdef WITH_PETSC
	VecRestoreArray(vec, &v);
#endif
}

// PETSc linear solver ///////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscLinearSolver::PetscLinearSolver(PetscMatrix &mat, PetscVector &rhs)
	: LinearSolver(), m(mat), rhs(rhs)
{
#ifdef WITH_PETSC
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}


PetscLinearSolver::~PetscLinearSolver() {
#ifdef WITH_PETSC
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

bool PetscLinearSolver::solve() {
#ifdef WITH_PETSC
	assert(m.ndofs == rhs.ndofs);

	m.finish();
	rhs.finish();

	PetscErrorCode ec;
	KSP ksp;
	Vec x;

	KSPCreate(PETSC_COMM_WORLD, &ksp);

	KSPSetOperators(ksp, m.matrix, m.matrix, DIFFERENT_NONZERO_PATTERN);
	KSPSetFromOptions(ksp);
	VecDuplicate(rhs.vec, &x);

	ec = KSPSolve(ksp, rhs.vec, x);
	if (ec) return false;

	// allocate memory for solution vector
	delete [] sln;
	sln = new scalar[m.ndofs + 1];
	MEM_CHECK(sln);
	sln[0] = 1.0;					// Dirichlet DOF has always coefficient 1.0

	// copy solution to the output solution vector
	// index map vector (basic serial code uses the map sln[i] = x[i] for all dofs.
	int *idx = new int[m.ndofs];
	MEM_CHECK(idx);
	for (int i = 0; i < m.ndofs; i++) idx[i] = i;
	VecGetValues(x, m.ndofs, idx, (PetscScalar *) (sln + 1));
	delete [] idx;

	KSPDestroy(ksp);
	VecDestroy(x);

	return true;
#else
	return false;
#endif
}
