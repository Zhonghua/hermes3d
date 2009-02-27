// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2009 Pavel Kus <pavel.kus@gmail.com>
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

int read_matrix_and_rhs(char *file_name, int &n, Array<MatrixEntry> &matrix, Array<double> &rhs) {
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
					matrix.add(MatrixEntry((int) buffer[0], (int) buffer[1], buffer[2]));
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

//
// tests themselves
//

void solve(LinearSolver &solver, int n, Array<MatrixEntry> &matrix, Array<double> &rhs) {
	// precalculate sparse structure
	solver.prealloc(n);

	for (Word_t i = matrix.first(); i != INVALID_IDX; i = matrix.next(i)) {
		MatrixEntry &me = matrix[i];
		solver.pre_add_ij(me.m, me.n);
	}

	solver.alloc();

	//
	solver.begin_assembling();

	for (Word_t i = matrix.first(); i != INVALID_IDX; i = matrix.next(i)) {
		MatrixEntry &me = matrix[i];
		solver.update_matrix(me.m, me.n, me.value);
	}

	for (Word_t i = rhs.first(); i != INVALID_IDX; i = rhs.next(i)) {
		solver.update_rhs(i, rhs[i]);
	}

	// !!! do not forget to call this !!!
	solver.finish_assembling();

	double *sln = new double [n];
	if (solver.solve_system(sln)) {
		for (int i = 0; i < n; i++) {
			printf("%lf\n", sln[i]);
		}
	}
	else {
		printf("Unable to solve\n");
	}

	delete [] sln;
}

int test_linear_solver(LinearSolver &solver, char *file_name) {
	Array<MatrixEntry> matrix;		// matrix
	Array<double> rhs;				// right-hand side
	int n;							// number of unknowns

	if (read_matrix_and_rhs(file_name, n, matrix, rhs) != ERR_SUCCESS)
		return ERR_FAILURE;

	solve(solver, n, matrix, rhs);

	return ERR_SUCCESS;
}

//
// testing
//

int test_linear_solver_block(LinearSolver &solver, char *file_name) {
	Array<MatrixEntry> matrix;		// matrix
	Array<double> rhs;				// right-hand side
	int n;							// number of unknowns

	if (read_matrix_and_rhs(file_name, n, matrix, rhs) != ERR_SUCCESS)
		return ERR_FAILURE;

	// precalculate sparse structure
	solver.prealloc(n);

	for (Word_t i = matrix.first(); i != INVALID_IDX; i = matrix.next(i)) {
		MatrixEntry &me = matrix[i];
		solver.pre_add_ij(me.m, me.n);
	}

	solver.alloc();

	//
	solver.begin_assembling();

	double **mat = new_matrix<double>(n, n);
	int *cols = new int[n];
	int *rows = new int[n];
	for (int i = 0; i < n; i++) {
		cols[i] = i;
		rows[i] = i;
	}

	for (Word_t i = matrix.first(); i != INVALID_IDX; i = matrix.next(i)) {
		MatrixEntry &me = matrix[i];
		mat[me.m][me.n] = me.value;
	}
	solver.update_matrix(n, n, mat, rows, cols);

	// right-hand side
	double *rs = new double[n];
	for (Word_t i = rhs.first(); i != INVALID_IDX; i = rhs.next(i)) {
		rs[i] = rhs[i];
	}
	solver.update_rhs(n, rows, rs);

	// !!! do not forget to call this !!!
	solver.finish_assembling();

	double *sln = new double [n];
	if (solver.solve_system(sln)) {
		for (int i = 0; i < n; i++) {
			printf("%lf\n", sln[i]);
		}
	}
	else {
		printf("Unable to solve\n");
	}

	delete [] sln;
	delete [] mat;

	return ERR_SUCCESS;
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

	if (strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
		PetscLinearSolver solver;
		ret = test_linear_solver(solver, argv[2]);
#endif
	}
	else if (strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
		PetscLinearSolver solver;
		ret = test_linear_solver_block(solver, argv[2]);
#endif
	}
	else if (strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
		UMFPackLinearSolver solver;
		ret = test_linear_solver(solver, argv[2]);
#endif
	}
	else if (strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
		UMFPackLinearSolver solver;
		ret = test_linear_solver_block(solver, argv[2]);
#endif
	}
	else if (strcasecmp(argv[1], "pardiso") == 0) {
#ifdef WITH_PARDISO
	    PardisoLinearSolver solver;
	    ret = test_linear_solver(solver, argv[2]);
#endif
	}
	else if (strcasecmp(argv[1], "pardiso-block") == 0) {
#ifdef WITH_PARDISO
	    PardisoLinearSolver solver;
		ret = test_linear_solver_block(solver, argv[2]);
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

