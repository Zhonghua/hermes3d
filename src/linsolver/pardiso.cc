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

#include <common/trace.h>
#include <common/error.h>

#include "../h3dconfig.h"
#include "pardiso.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int pardisoinit_(void *, int *, int *);

extern int
    pardiso_(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *);

#define PARDISOINIT pardisoinit_
#define PARDISO pardiso_

#ifdef __cplusplus
}
#endif

PardisoLinearSolver::PardisoLinearSolver() {
#ifdef WITH_PARDISO
	Ap = NULL;
	Ai = NULL;
	Ax = NULL;
	srhs = NULL;
#else
	EXIT(ERR_PARDISO_NOT_COMPILED);
#endif
}

PardisoLinearSolver::~PardisoLinearSolver() {
#ifdef WITH_PARDISO
	free();
#else
	EXIT(ERR_PARDISO_NOT_COMPILED);
#endif
}

void PardisoLinearSolver::prealloc(int ndofs) {
#ifdef WITH_PARDISO
	free();

	this->ndofs = ndofs;

	pages = new Page *[ndofs];
	if (pages == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(pages, 0, ndofs * sizeof(Page *));
#endif
}

void PardisoLinearSolver::pre_add_ij(int row, int col) {
#ifdef WITH_PARDISO
	int tmp = row;
	row = col;
	col = tmp;

	if (pages[col] == NULL || pages[col]->count >= PAGE_SIZE) {
		Page *new_page = new Page;
		new_page->count = 0;
		new_page->next = pages[col];
		pages[col] = new_page;
	}
	pages[col]->idx[pages[col]->count++] = row;
#endif
}

void PardisoLinearSolver::alloc() {
#ifdef WITH_PARDISO
	free();
	assert(pages != NULL);

	// initialize the arrays Ap and Ai
	Ap = new int[ndofs + 1];
	int aisize = get_num_indices(pages, ndofs);
	Ai = new int[aisize];
	if (Ai == NULL) ERROR(ERR_OUT_OF_MEMORY);

	// sort the indices and remove duplicities, insert into Ai
	int i, pos = 0;
	for (i = 0; i < ndofs; i++) {
		Ap[i] = pos;
		pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
	}
	Ap[i] = pos;

	delete[] pages;
	pages = NULL;

	Ax = new scalar[Ap[ndofs]];
	if (Ax == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(Ax, 0, sizeof(scalar) * Ap[ndofs]);

	// RHS
	srhs = new scalar[ndofs];
	if (srhs == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(srhs, 0, ndofs * sizeof(scalar));
#endif
}

void PardisoLinearSolver::free() {
#ifdef WITH_PARDISO
	delete[] Ap; Ap = NULL;
	delete[] Ai; Ai = NULL;
	delete[] Ax; Ax = NULL;
	delete[] srhs; 	srhs = NULL;
#endif
}

void PardisoLinearSolver::update_matrix(int row, int col, scalar v) {
#ifdef WITH_PARDISO
	int tmp = row;
	row = col;
	col = tmp;

	insert_value(Ai + Ap[col], Ax + Ap[col], Ap[col + 1] - Ap[col], row, v);
#endif
}

void PardisoLinearSolver::update_matrix(int m, int n, double **mat, int *rows, int *cols) {
#ifdef WITH_PARDISO
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)	 {		// cols
			if (mat[i][j] != 0.0 && rows[i] != -1 && cols[j] != -1) {		// -1 is a "dirichlet DOF" -> ignore it
				update_matrix(rows[i], cols[j], mat[i][j]);
			}
		}
#endif
}

void PardisoLinearSolver::update_rhs(int idx, scalar y) {
#ifdef WITH_PARDISO
	if (idx >= 0) srhs[idx] += y;
#endif
}

void PardisoLinearSolver::update_rhs(int n, int *idx, scalar *y) {
#ifdef WITH_PARDISO
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0) srhs[idx[i]] += y[i];
#endif
}

bool PardisoLinearSolver::solve_system(double *sln) {
#ifdef WITH_PARDISO
	bool res = true;

	int n = ndofs;
	int *ia = Ap;
	int *ja = Ai;
	double *a = Ax;
	double *b = srhs;
	double *x = sln;

	//TODO: improve error handling

	try {
		// Numbers of processors, value of OMP_NUM_THREADS
		int num_procs;
		char *var = getenv("OMP_NUM_THREADS");
		if (var != NULL) sscanf(var, "%d", &num_procs);
		else num_procs = 1;

		int mtype = 11; // Real unsymmetric matrix
		int nrhs = 1; // Number of right hand sides.
		int nnz = ia[n]; //number of nonzero elements

		// Internal solver memory pointer pt,
		// 32-bit: int pt[64]; 64-bit: long int pt[64]
		// or void *pt[64] should be OK on both architectures
		void *pt[64];
		// Pardiso control parameters.
		int iparm[64];
		int maxfct, mnum, phase, error, msglvl;
		// Auxiliary variables.
		double ddum; // Double dummy
		int idum; // Integer dummy.

		iparm[2] = num_procs;

		maxfct = 1; // Maximum number of numerical factorizations.
		mnum = 1; // Which factorization to use.
		//  msglvl = 1; // Print statistical information in file
		msglvl = 0; // Do not print statistical information
		error = 0; // Initialize error flag

		// Convert matrix from 0-based C-notation to Fortran 1-based notation.
		for (int i = 0; i < n + 1; i++) ia[i] += 1;
		for (int i = 0; i < nnz; i++) ja[i] += 1;

		// Setup Pardiso control parameters.
		PARDISOINIT(pt, &mtype, iparm);

		// .. Reordering and Symbolic Factorization. This step also allocates
		// all memory that is necessary for the factorization.
		phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			ERROR("ERROR during symbolic factorization: ", error);
			throw ERR_FAILURE;
		}

		DEBUG_PRINT("Reordering completed ... ");
		DEBUG_PRINT("Number of nonzeros in factors = ", iparm[17]);
		DEBUG_PRINT("Number of factorization MFLOPS = ", iparm[18]);

		// .. Numerical factorization.
		phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			ERROR("ERROR during numerical factorization: ", error);
			throw ERR_FAILURE;
		}
		DEBUG_PRINT("Factorization completed ... ");

		// .. Back substitution and iterative refinement.
		phase = 33;
		iparm[7] = 1; // Max numbers of iterative refinement steps.
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
		if (error != 0) {
			ERROR("ERROR during solution: ", error);
			throw ERR_FAILURE;
		}
		DEBUG_PRINT("Solving completed ... ");

		//  Convert matrix back to 0-based C-notation.
		for (int i = 0; i < n + 1; i++) ia[i] -= 1;
		for (int i = 0; i < nnz; i++) ja[i] -= 1;

		// .. Termination and release of memory.
		phase = -1; // Release internal memory.
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	}
	catch (int e) {
		res = false;
	}

	return res;
#else
	return false;
#endif
}

bool PardisoLinearSolver::dump_matrix(FILE *file, const char *var_name, EMatrixDumpFormat format/* = DF_MATLAB_SPARSE*/) {
#ifdef WITH_PARDISO
	// TODO: check if OK (use unsymmetric matrix)
	switch (format) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", ndofs, ndofs, Ap[ndofs], Ap[ndofs]);
			for (int j = 0; j < ndofs; j++)
				for (int i = Ap[j]; i < Ap[j + 1]; i++)
					fprintf(file, "%d %d %.18e\n", j + 1, Ai[i] + 1, Ax[i]);
			fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

			return true;

		case DF_HERMES_BIN:
		case DF_PLAIN_ASCII:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
#else
	return false;
#endif
}

bool PardisoLinearSolver::dump_rhs(FILE *file, const char *var_name, EMatrixDumpFormat format/* = DF_MATLAB_SPARSE*/) {
#ifdef WITH_PARDISO
	switch (format) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx1\n%s = [\n", ndofs, var_name);
			for (int i = 0; i < this->ndofs; i++)
				fprintf(file, "%.18e\n", srhs[i]);
			fprintf(file, " ];\n");
			return true;

		case DF_HERMES_BIN:
		case DF_PLAIN_ASCII:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
#else
	return false;
#endif
}

int PardisoLinearSolver::get_matrix_size() const {
#ifdef WITH_PARDISO
	assert(Ap != NULL);
	return (sizeof(int) + sizeof(scalar)) * (Ap[ndofs] + ndofs);
#else
	return -1;
#endif
}

