// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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

/*
 * main.cc
 *
 * Test of linear solver
 * Read matrix and RHS from file and print out the result
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>
#include <common/array.h>

// max row length in input file
#define MAX_ROW_LEN									1024

struct MatrixEntry {
	MatrixEntry() { }
	MatrixEntry(int m, int n, double value) {
		this->m = m;
		this->n = n;
		this->value = value;
	}

	int m, n;			// position
	double value;
};

// helpers ////////////////////////////////////////////////////////////////////

bool testPrint(bool value, const char *msg, bool correct) {
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

bool read_n_nums(char *row, int n, double values[]) {
	int i = 0;
	char delims[] = " \t\n\r";
	char *token = strtok(row, delims);
	while (token != NULL && i < n) {
		int n;
		sscanf(token, "%d", &n);
		values[i++] = n;

		token = strtok(NULL, delims);
	}

	return (i == n);
}

int read_matrix_and_rhs(char *file_name, int &n, Array<MatrixEntry> &mat, Array<double> &rhs) {
	FILE *file = fopen(file_name, "r");
	if (file == NULL)
		return ERR_CAN_NOT_OPEN_FILE;

	enum EState {
		STATE_N,
		STATE_MATRIX,
		STATE_RHS,
	} state = STATE_N;

	double buffer[3];
	char row[MAX_ROW_LEN];
	while (fgets(row, MAX_ROW_LEN, file) != NULL) {
		switch (state) {
			case STATE_N:
				if (read_n_nums(row, 1, buffer)) {
					n = (int) buffer[0];
					state = STATE_MATRIX;
				}
				break;

			case STATE_MATRIX:
				if (read_n_nums(row, 3, buffer)) {
					mat.add(MatrixEntry((int) buffer[0], (int) buffer[1], buffer[2]));
				}
				else
					state = STATE_RHS;
				break;

			case STATE_RHS:
				if (read_n_nums(row, 2, buffer)) {
					rhs[(int) buffer[0]] = buffer[1];
				}
				break;
		}
	}

	fclose(file);

	return ERR_SUCCESS;
}

void build_matrix(int n, Array<MatrixEntry> &ar_mat, Array<double> &ar_rhs, SparseMatrix *mat, Vector *rhs) {
	// matrix
	mat->prealloc(n);
	for (Word_t i = ar_mat.first(); i != INVALID_IDX; i = ar_mat.next(i)) {
		MatrixEntry &me = ar_mat[i];
		mat->pre_add_ij(me.m, me.n);
	}
	mat->alloc();

	for (Word_t i = ar_mat.first(); i != INVALID_IDX; i = ar_mat.next(i)) {
		MatrixEntry &me = ar_mat[i];
		mat->update(me.m, me.n, me.value);
	}

	// RHS
	rhs->alloc(n);
	for (Word_t i = ar_rhs.first(); i != INVALID_IDX; i = ar_rhs.next(i)) {
		rhs->update((int) i, ar_rhs[i]);
	}
}

void build_matrix_block(int n, Array<MatrixEntry> &ar_mat, Array<double> &ar_rhs, SparseMatrix *matrix, Vector *rhs) {
	// matrix
	matrix->prealloc(n);
	for (Word_t i = ar_mat.first(); i != INVALID_IDX; i = ar_mat.next(i)) {
		MatrixEntry &me = ar_mat[i];
		matrix->pre_add_ij(me.m, me.n);
	}
	matrix->alloc();

	double **mat = new_matrix<double>(n, n);
	int *cols = new int[n];
	int *rows = new int[n];
	for (int i = 0; i < n; i++) {
		cols[i] = i;
		rows[i] = i;
	}
	for (Word_t i = ar_mat.first(); i != INVALID_IDX; i = ar_mat.next(i)) {
		MatrixEntry &me = ar_mat[i];
		mat[me.m][me.n] = me.value;
	}
	matrix->update(n, n, mat, rows, cols);

	// rhs
	rhs->alloc(n);
	double *rs = new double[n];
	for (Word_t i = ar_rhs.first(); i != INVALID_IDX; i = ar_rhs.next(i)) {
		rs[i] = ar_rhs[i];
	}
	rhs->update(n, rows, rs);
}

//
// tests themselves
//

void solve(LinearSolver &solver, int n) {
	if (solver.solve()) {
		double *sln = solver.get_solution();
		for (int i = 1; i < n + 1; i++) {
			printf("%lf\n", sln[i]);
		}
	}
	else {
		printf("Unable to solve\n");
	}
}

int main(int argc, char *argv[]) {
//	TRACE_START("trace.txt");
	DEBUG_OUTPUT_OFF;
	SET_VERBOSE_LEVEL(0);

	int ret = ERR_SUCCESS;

#ifdef WITH_PETSC
	// do NOT forget to call this when using PETSc solver
	PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
	// disable PETSc error handler
	PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);
#endif

	if (argc < 3)
		return ERR_NOT_ENOUGH_PARAMS;

	int n;
	Array<MatrixEntry> ar_mat;
	Array<double> ar_rhs;

	if (read_matrix_and_rhs(argv[2], n, ar_mat, ar_rhs) != ERR_SUCCESS)
		return ERR_FAILURE;

	if (strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
		PetscMatrix mat;
		PetscVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		PetscLinearSolver solver(mat, rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
		PetscMatrix mat;
		PetscVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		PetscLinearSolver solver(mat, rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
		UMFPackMatrix mat;
		UMFPackVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		UMFPackLinearSolver solver(mat, rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
		UMFPackMatrix mat;
		UMFPackVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		UMFPackLinearSolver solver(mat, rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "pardiso") == 0) {
#ifdef WITH_PARDISO
#endif
	}
	else if (strcasecmp(argv[1], "pardiso-block") == 0) {
#ifdef WITH_PARDISO
#endif
	}
	else
		ret = ERR_FAILURE;

#ifdef WITH_PETSC
	// do NOT forget to call this when using PETSc solver
	PetscFinalize();
#endif

	return ret;
}

