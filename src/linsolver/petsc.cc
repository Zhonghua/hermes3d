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

#include "../h3dconfig.h"

#ifdef WITH_PETSC
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#endif

#include "petsc.h"
#include <common/trace.h>
#include <common/error.h>

// all declarations from PETSc must be in .cc file, because PETSc is optional
struct msyst {
#ifdef WITH_PETSC
	// PETSc specific data structures for storing matrix, rhs, vector of unknowns
	Mat matrix;
	Vec rhs;
	Vec x;
	KSP ksp;
	PetscErrorCode ec;
#endif
};

PetscLinearSolver::PetscLinearSolver() {
#ifdef WITH_PETSC
	ms = NULL;
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}


PetscLinearSolver::~PetscLinearSolver() {
#ifdef WITH_PETSC
	free();
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::prealloc(int ndofs) {
#ifdef WITH_PETSC
	this->ndofs = ndofs;

	pages = new Page *[ndofs];
	MEM_CHECK(pages);
	memset(pages, 0, ndofs * sizeof(Page *));
#endif
}

void PetscLinearSolver::pre_add_ij(int row, int col) {
#ifdef WITH_PETSC
	if (pages[row] == NULL || pages[row]->count >= PAGE_SIZE) {
		Page *new_page = new Page;
		MEM_CHECK(new_page);
		new_page->count = 0;
		new_page->next = pages[row];
		pages[row] = new_page;
	}
	pages[row]->idx[pages[row]->count++] = col;
#endif
}

void PetscLinearSolver::alloc() {
#ifdef WITH_PETSC
	assert(pages != NULL);
	this->ms = new msyst;
	MEM_CHECK(this->ms);

	// calc nnz
	int *nnz = new int[this->ndofs];
	MEM_CHECK(nnz);

	// fill in nnz
	int aisize = LinearSolver::get_num_indices(pages, ndofs);
	int *ai = new int[aisize];
	MEM_CHECK(ai);

	// sort the indices and remove duplicities, insert into ai
	int pos = 0;
	for (int i = 0; i < ndofs; i++) {
		nnz[i] = LinearSolver::sort_and_store_indices(pages[i], ai + pos, ai + aisize);
		pos += nnz[i];
	}
	delete [] pages; pages = NULL;
	delete [] ai;

	//
	MatCreateSeqAIJ(PETSC_COMM_SELF, this->ndofs, this->ndofs, 0, nnz, &ms->matrix);
//	MatSetOption(ms->matrix, MAT_ROW_ORIENTED);
//	MatSetOption(ms->matrix, MAT_ROWS_SORTED);

	// create rhs vector
	VecCreateSeq(PETSC_COMM_SELF, this->ndofs, &ms->rhs);
//	VecSetOption(ms->rhs, VEC_IGNORE_NEGATIVE_INDICES); // in some newer version
	// create vector of unknowns
	VecDuplicate(ms->rhs, &ms->x);

	delete [] nnz;

	//
	KSPCreate(PETSC_COMM_WORLD, &ms->ksp);
#endif
}

void PetscLinearSolver::free() {
#ifdef WITH_PETSC
	if (ms != NULL) {
		// commented out, becuase it is causing SEGFAULTS
		MatDestroy(ms->matrix);
		VecDestroy(ms->rhs);
		VecDestroy(ms->x);
		KSPDestroy(ms->ksp);

		delete ms;
		ms = NULL;
	}
#endif
}

void PetscLinearSolver::update_matrix(int m, int n, scalar v) {
#ifdef WITH_PETSC
	assert(ms != NULL);
	MatSetValues(ms->matrix, 1, &m, 1, &n, (PetscScalar *) &v, ADD_VALUES);
#endif
}

void PetscLinearSolver::update_matrix(int m, int n, scalar **mat, int *rows, int *cols) {
#ifdef WITH_PETSC
	assert(ms != NULL);
	// TODO: pass in just the block of the matrix without DIRICHLET_DOFs (so that can use MatSetValues directly without checking row and cols for -1)
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)	 {		// cols
			if (mat[i][j] != 0.0 && rows[i] != DIRICHLET_DOF && cols[j] != DIRICHLET_DOF) {		// ignore "dirichlet DOF"
				MatSetValues(ms->matrix, 1, rows + i, 1, cols + j, (PetscScalar *) &(mat[i][j]), ADD_VALUES);
			}
		}
#endif
}

void PetscLinearSolver::update_rhs(int idx, scalar y) {
#ifdef WITH_PETSC
	assert(ms != NULL);
	if (idx >= 0) VecSetValues(ms->rhs, 1, &idx, (PetscScalar *) &y, ADD_VALUES);
#endif
}

void PetscLinearSolver::update_rhs(int n, int *idx, scalar *y) {
#ifdef WITH_PETSC
	assert(ms != NULL);
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0) VecSetValues(ms->rhs, 1, idx + i, (PetscScalar *) (y + i), ADD_VALUES);
#endif
}

bool PetscLinearSolver::solve_system(scalar *sln) {
#ifdef WITH_PETSC
	KSPSetOperators(ms->ksp, ms->matrix, ms->matrix, DIFFERENT_NONZERO_PATTERN);
	KSPSetFromOptions(ms->ksp);
	ms->ec = KSPSolve(ms->ksp, ms->rhs, ms->x);
	if (ms->ec) return false;

	// copy solution to the output solution vector
	// index map vector (basic serial code uses the map sln[i] = x[i] for all dofs.
	int *idx = new int[ndofs];
	MEM_CHECK(idx);
	for (int i = 0; i < ndofs; i++) idx[i] = i;
	VecGetValues(ms->x, ndofs, idx, (PetscScalar *) sln);
	delete[] idx;

	return true;
#else
	return false;
#endif
}

void PetscLinearSolver::begin_assembling() {
}

void PetscLinearSolver::finish_assembling() {
#ifdef WITH_PETSC
	MatAssemblyBegin(ms->matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ms->matrix, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(ms->rhs);
	VecAssemblyEnd(ms->rhs);
#endif
}

bool PetscLinearSolver::dump_matrix(FILE *file, const char *var_name, EMatrixDumpFormat format) {
#ifdef WITH_PETSC
//	PetscViewer v;
//	PetscViewerASCIIOpen(PETSC_COMM_SELF, file_name, &v);
//	PetscViewerSetFormat(v, PETSC_VIEWER_ASCII_MATLAB);
//	MatView(ms->matrix, v);
//	PetscViewerDestroy(v);
#endif
	return false;
}

bool PetscLinearSolver::dump_rhs(FILE *file, const char *var_name, EMatrixDumpFormat format) {
#ifdef WITH_PETSC
//	PetscViewer v;
//	PetscViewerASCIIOpen(PETSC_COMM_SELF, file_name, &v);
//	PetscViewerSetFormat(v, PETSC_VIEWER_ASCII_MATLAB);
//	VecView(ms->rhs, v);
//	PetscViewerDestroy(v);
#endif
	return false;
}

