//
// umfpacksolver.cc
//

#include "../config.h"

#ifdef USE_UMFPACK
extern "C" {
#include <umfpack.h>
}
#endif

#include "umfpack.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/utils.h>
#include <iostream>

#ifndef COMPLEX
  // real case
  #define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I)   umfpack_di_symbolic(m, n, Ap, Ai, Ax, S, C, I)
  #define umfpack_numeric(Ap, Ai, Ax, S, N, C, I)       umfpack_di_numeric(Ap, Ai, Ax, S, N, C, I)
  #define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I) umfpack_di_solve(sys, Ap, Ai, Ax, X, B, N, C, I)
  #define umfpack_free_symbolic                         umfpack_di_free_symbolic
  #define umfpack_free_numeric                          umfpack_di_free_numeric
  #define umfpack_defaults                              umfpack_di_defaults
#else
  // macros for calling complex UMFPACK in packed-complex mode
  #define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I) \
          umfpack_zi_symbolic(m, n, Ap, Ai, (double*) (Ax), NULL, S, C, I)
  #define umfpack_numeric(Ap, Ai, Ax, S, N, C, I) \
          umfpack_zi_numeric(Ap, Ai, (double*) (Ax), NULL, S, N, C, I)
  #define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I) \
          umfpack_zi_solve(sys, Ap, Ai, (double*) (Ax), NULL, (double*) (X), NULL, (double*) (B), NULL, N, C, I)
  #define umfpack_free_symbolic umfpack_di_free_symbolic
  #define umfpack_free_numeric  umfpack_zi_free_numeric
  #define umfpack_defaults      umfpack_zi_defaults
#endif


UMFPackLinearSolver::UMFPackLinearSolver() {
#ifdef USE_UMFPACK
	Ap = NULL;
	Ai = NULL;
	Ax = NULL;
	srhs = NULL;
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

UMFPackLinearSolver::~UMFPackLinearSolver() {
#ifdef USE_UMFPACK
	free();
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::prealloc(int ndofs) {
#ifdef USE_UMFPACK
	free();

	this->ndofs = ndofs;

	pages = new Page *[ndofs];
	if (pages == NULL) EXIT(ERR_OUT_OF_MEMORY, "Out of memory. Error pre-allocating pages.");
	memset(pages, 0, ndofs * sizeof(Page *));
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::pre_add_ij(int row, int col) {
#ifdef USE_UMFPACK
	if (pages[col] == NULL || pages[col]->count >= PAGE_SIZE) {
		Page *new_page = new Page;
		if (new_page == NULL) EXIT(ERR_OUT_OF_MEMORY, "Out of memory. Error allocating a page.");
		new_page->count = 0;
		new_page->next = pages[col];
		pages[col] = new_page;
	}
	pages[col]->idx[pages[col]->count++] = row;
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::alloc() {
#ifdef USE_UMFPACK
	free();
	assert(pages != NULL);

	// initialize the arrays Ap and Ai
	Ap = new int[ndofs + 1];
	if (Ap == NULL) EXIT(ERR_OUT_OF_MEMORY, "Out of memory. Error allocating stiffness matrix (Ap).");
	int aisize = get_num_indices(pages, ndofs);
	Ai = new int[aisize];
	if (Ai == NULL) EXIT(ERR_OUT_OF_MEMORY, "Out of memory. Error allocating stiffness matrix (Ai).");

	// sort the indices and remove duplicities, insert into Ai
	int i, pos = 0;
	for (i = 0; i < ndofs; i++) {
		Ap[i] = pos;
		pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
	}
	Ap[i] = pos;

	delete [] pages; pages = NULL;

	Ax = new scalar[Ap[ndofs]];
	if (Ax == NULL) EXIT(ERR_OUT_OF_MEMORY, "Out of memory. Error allocating stiffness matrix (Ax).");
	memset(Ax, 0, sizeof(scalar) * Ap[ndofs]);

	// RHS
	srhs = new scalar[ndofs];
	if (srhs == NULL) EXIT(ERR_OUT_OF_MEMORY, "Out of memory. Error allocating the RHS vector.");
	memset(srhs, 0, ndofs * sizeof(scalar));
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::free() {
#ifdef USE_UMFPACK
	delete [] Ap; Ap = NULL;
	delete [] Ai; Ai = NULL;
	delete [] Ax; Ax = NULL;
	delete [] srhs; srhs = NULL;
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::update_matrix(int row, int col, scalar v) {
#ifdef USE_UMFPACK
	insert_value(Ai + Ap[col], Ax + Ap[col], Ap[col + 1] - Ap[col], row, v);
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::update_matrix(int m, int n, scalar **mat, int *rows, int *cols) {
#ifdef USE_UMFPACK
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)	 {		// cols
			if (mat[i][j] != 0.0 && rows[i] != -1 && cols[j] != -1) {		// -1 is a "dirichlet DOF" -> ignore it
				update_matrix(rows[i], cols[j], mat[i][j]);
			}
		}
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}


void UMFPackLinearSolver::update_rhs(int idx, scalar y) {
#ifdef USE_UMFPACK
	if (idx >= 0)
		srhs[idx] += y;
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

void UMFPackLinearSolver::update_rhs(int n, int *idx, scalar *y) {
#ifdef USE_UMFPACK
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0)
			srhs[idx[i]] += y[i];
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

#ifdef USE_UMFPACK
static void check_status(const char *fn_name, int status) {
	switch (status) {
		case UMFPACK_OK:							break;
	    case UMFPACK_WARNING_singular_matrix:		ERROR("%s: singular matrix!", fn_name); break;
		case UMFPACK_ERROR_out_of_memory:			ERROR("%s: out of memory!", fn_name); break;
		case UMFPACK_ERROR_argument_missing:		ERROR("%s: argument missing", fn_name); break;
		case UMFPACK_ERROR_invalid_Symbolic_object:	ERROR("%s: invalid Symbolic object", fn_name); break;
		case UMFPACK_ERROR_invalid_Numeric_object:	ERROR("%s: invalid Numeric object", fn_name); break;
		case UMFPACK_ERROR_different_pattern:		ERROR("%s: different pattern", fn_name); break;
		case UMFPACK_ERROR_invalid_system:			ERROR("%s: invalid system", fn_name); break;
		case UMFPACK_ERROR_n_nonpositive:			ERROR("%s: n nonpositive", fn_name); break;
		case UMFPACK_ERROR_invalid_matrix:			ERROR("%s: invalid matrix", fn_name); break;
		case UMFPACK_ERROR_internal_error:			ERROR("%s: internal error", fn_name); break;
		default:									ERROR("%s: unknown error (%d)", fn_name, status); break;
	}
}
#endif

bool UMFPackLinearSolver::solve_system(scalar *sln) {
#ifdef USE_UMFPACK
	bool res = true;

	void *symbolic, *numeric;
	int status;

	status = umfpack_symbolic(ndofs, ndofs, Ap, Ai, Ax, &symbolic, NULL, NULL);
	if (status != UMFPACK_OK) { check_status("umfpack_di_symbolic", status); return false; }
	if (symbolic == NULL) EXIT(ERR_FAILURE, "umfpack_di_symbolic error: symbolic == NULL");

	status = umfpack_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL);
	if (status != UMFPACK_OK) { check_status("umfpack_di_numeric", status); return false; }
	if (numeric == NULL) EXIT(ERR_FAILURE, "umfpack_di_numeric error: numeric == NULL");

	status = umfpack_solve(UMFPACK_A, Ap, Ai, Ax, sln, srhs, numeric, NULL, NULL);
	if (status != UMFPACK_OK) { check_status("umfpack_di_solve", status); return false; }

	umfpack_free_symbolic(&symbolic);
	umfpack_free_numeric(&numeric);

	return true;
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

bool UMFPackLinearSolver::dump_matrix(FILE *file, const char *var_name, EMatrixDumpFormat format/* = DF_MATLAB_SPARSE*/) {
#ifdef USE_UMFPACK
	switch (format) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n",
				ndofs, ndofs, Ap[ndofs], Ap[ndofs]);
			for (int j = 0; j < ndofs; j++)
				for (int i = Ap[j]; i < Ap[j + 1]; i++)
					fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, SCALAR(Ax[i]));
			fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

			return true;

		case DF_HERMES_BIN: {
			hermes_fwrite("H2DX\001\000\000\000", 1, 8, file);
			int ssize = sizeof(scalar);
			int nnz = Ap[ndofs];
			hermes_fwrite(&ssize, sizeof(int), 1, file);
			hermes_fwrite(&ndofs, sizeof(int), 1, file);
			hermes_fwrite(&nnz, sizeof(int), 1, file);
			hermes_fwrite(Ap, sizeof(int), ndofs + 1, file);
			hermes_fwrite(Ai, sizeof(int), nnz, file);
			hermes_fwrite(Ax, sizeof(scalar), nnz, file);
			return true;
		}

		case DF_PLAIN_ASCII:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

bool UMFPackLinearSolver::dump_rhs(FILE *file, const char *var_name, EMatrixDumpFormat format/* = DF_MATLAB_SPARSE*/) {
#ifdef USE_UMFPACK
	switch (format) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx1\n%s = [\n", ndofs, var_name);
			for (int i = 0; i < this->ndofs; i++)
				fprintf(file, SCALAR_FMT "\n", SCALAR(srhs[i]));
			fprintf(file, " ];\n");
			return true;

			case DF_HERMES_BIN: {
				hermes_fwrite("H2DR\001\000\000\000", 1, 8, file);
				int ssize = sizeof(scalar);
				hermes_fwrite(&ssize, sizeof(int), 1, file);
				hermes_fwrite(&ndofs, sizeof(int), 1, file);
				hermes_fwrite(srhs, sizeof(scalar), ndofs, file);
				return true;
			}

		case DF_PLAIN_ASCII:
			EXIT(ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
#else
	EXIT(ERR_UMFPACK_NOT_COMPILED);
#endif
}

int UMFPackLinearSolver::get_matrix_size() const {
	assert(Ap != NULL);
	return (sizeof(int) + sizeof(scalar)) * (Ap[ndofs] + ndofs);
}
