//
// petsc.cc
//

#include "../config.h"

#ifdef USE_PETSC
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
#ifdef USE_PETSC
	// PETSc specific data structures for storing matrix, rhs, vector of unknowns
	Mat matrix;
	Vec rhs;
	Vec x;
#endif
};


PetscLinearSolver::PetscLinearSolver() {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
//	initialized = false;
	ms = new msyst;
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

PetscLinearSolver::~PetscLinearSolver() {
#ifdef USE_PETSC
	free();
	delete ms;
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::prealloc(int ndofs) {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::pre_add_ij(int row, int col) {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::alloc() {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);

	free();

	this->ndofs = ndofs;
/*
	// stiffness matrix
	int *o_nnz = new int[ndofs];
	if (o_nnz == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(o_nnz, 0, ndofs * sizeof(int));

	MatCreateMPIAIJ(PETSC_COMM_WORLD, ndofs, ndofs, ndofs, ndofs, 0, d_nnz, 0, o_nnz, &(ms->matrix));

	// fill the matrix with zero entries
	double *dbuf = new double[lbuf];
	if (dbuf == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(dbuf, 0, lbuf * sizeof(double));
	for (int i = 0; i < ndofs; i++)
		MatSetValues(ms->matrix, 1, &i, d_nnz[i], idx[i], dbuf, INSERT_VALUES);

	MatAssemblyBegin(ms->matrix, MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(ms->matrix, MAT_FLUSH_ASSEMBLY);

	delete [] dbuf;
	delete [] o_nnz;

	// RHS
	VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ndofs, &(ms->rhs));
	// vector of unknowns
	VecDuplicate(ms->rhs, &(ms->x));

	initialized = true;
*/
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::free() {
#ifdef USE_PETSC
	if (initialized) {
		// commented out, becuase it is causing SEGFAULTS
//		MatDestroy(ms->matrix);
//		VecDestroy(ms->rhs);
//		VecDestroy(ms->x);
		initialized = false;
	}
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::update_matrix(int m, int n, scalar v) {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
//	MatSetValues(ms->matrix, 1, &m, 1, &n, &v, ADD_VALUES);
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::update_rhs(int idx, scalar y) {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
//	VecSetValues(ms->rhs, 1, &idx, &y, ADD_VALUES);
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

bool PetscLinearSolver::solve_system(double *sln) {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
/*	KSP ksp;
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetOperators(ksp, ms->matrix, ms->matrix, DIFFERENT_NONZERO_PATTERN);
	KSPSetFromOptions(ksp);
	PetscErrorCode ec = KSPSolve(ksp, ms->rhs, ms->x);
	if (ec) return false;
	KSPDestroy(ksp);

	// copy solution to the output solution vector
	// index map vector (basic serial code uses the map sln[i] = x[i] for all dofs.
	int *idx = new int[ndofs];
	if (idx == NULL) EXIT(ERR_OUT_OF_MEMORY);
	for (int i = 0; i < ndofs; i++)
		idx[i] = i;
	VecGetValues(ms->x, ndofs, idx, sln);
	delete [] idx;
*/
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif

	return true;
}

void PetscLinearSolver::begin_assembling() {
#ifdef USE_PETSC
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::finish_assembling() {
#ifdef USE_PETSC
	EXIT(ERR_NOT_IMPLEMENTED);
/*	MatAssemblyBegin(ms->matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ms->matrix, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(ms->rhs);
	VecAssemblyEnd(ms->rhs);
*/
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::dump_matrix(FILE *file, EMatrixDumpFormat format/* = DF_MATLAB_SPARSE*/) {
#ifdef USE_PETSC
/*	FIXME: dump to a stream
	PetscViewer v;
	PetscViewerASCIIOpen(PETSC_COMM_SELF, file_name, &v);
	PetscViewerSetFormat(v, PETSC_VIEWER_ASCII_MATLAB);
	MatView(ms->matrix, v);
	PetscViewerDestroy(v);
*/
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

void PetscLinearSolver::dump_rhs(FILE *file, EMatrixDumpFormat format/* = DF_MATLAB_SPARSE*/) {
#ifdef USE_PETSC
/*	FIXME: dump to a stream
	PetscViewer v;
	PetscViewerASCIIOpen(PETSC_COMM_SELF, file_name, &v);
	PetscViewerSetFormat(v, PETSC_VIEWER_ASCII_MATLAB);
	VecView(ms->rhs, v);
	PetscViewerDestroy(v);
*/
#else
	EXIT(ERR_PETSC_NOT_COMPILED);
#endif
}

