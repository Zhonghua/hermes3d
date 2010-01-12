// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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
	MatrixEntry(int m, int n, scalar value) {
		this->m = m;
		this->n = n;
		this->value = value;
	}

	int m, n;			// position
	scalar value;
};

// helpers ////////////////////////////////////////////////////////////////////

bool testPrint(bool value, const char *msg, bool correct)
{
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

bool read_n_nums(char *row, int n, double values[])
{
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

int read_matrix_and_rhs(char *file_name, int &n, Array<MatrixEntry> &mat, Array<scalar> &rhs)
{
#ifndef COMPLEX
	FILE *file = fopen(file_name, "r");
	if (file == NULL) return ERR_FAILURE;

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
#else
	n = 3;
	mat.add(MatrixEntry(0, 0, scalar(1, 2)));
	mat.add(MatrixEntry(1, 1, scalar(1, 4)));
	mat.add(MatrixEntry(2, 2, scalar(1, 6)));

	rhs[0] = scalar(2, 1);
	rhs[1] = scalar(4, 1);
	rhs[2] = scalar(6, 2);
#endif
	return ERR_SUCCESS;
}

void build_matrix(int n, Array<MatrixEntry> &ar_mat, Array<scalar> &ar_rhs, SparseMatrix *mat,
                  Vector *rhs)
{
	// matrix
	mat->prealloc(n);
	for (Word_t i = ar_mat.first(); i != INVALID_IDX; i = ar_mat.next(i)) {
		MatrixEntry &me = ar_mat[i];
		mat->pre_add_ij(me.m, me.n);
	}

	mat->alloc();
	for (Word_t i = ar_mat.first(); i != INVALID_IDX; i = ar_mat.next(i)) {
		MatrixEntry &me = ar_mat[i];
		mat->add(me.m, me.n, me.value);
	}
	mat->finish();

	// RHS
	rhs->alloc(n);
	for (Word_t i = ar_rhs.first(); i != INVALID_IDX; i = ar_rhs.next(i)) {
		rhs->add((int) i, ar_rhs[i]);
	}
	rhs->finish();
}

void build_matrix_block(int n, Array<MatrixEntry> &ar_mat, Array<scalar> &ar_rhs,
                        SparseMatrix *matrix, Vector *rhs)
{
	// matrix
	matrix->prealloc(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			matrix->pre_add_ij(i, j);

	matrix->alloc();
	scalar **mat = new_matrix<scalar>(n, n);
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
	matrix->add(n, n, mat, rows, cols);
	matrix->finish();

	// rhs
	rhs->alloc(n);
	scalar *rs = new scalar[n];
	for (Word_t i = ar_rhs.first(); i != INVALID_IDX; i = ar_rhs.next(i)) {
		rs[i] = ar_rhs[i];
	}
	rhs->add(n, rows, rs);
	rhs->finish();
}

//
// tests themselves
//

void solve(Solver &solver, int n)
{
	if (solver.solve()) {
		scalar *sln = solver.get_solution();
		for (int i = 0; i < n; i++) {
			printf(SCALAR_FMT"\n", SCALAR(sln[i]));
		}
	}
	else {
		printf("Unable to solve\n");
	}
}

int main(int argc, char *argv[])
{
	int ret = ERR_SUCCESS;

#ifdef WITH_PETSC
	// do NOT forget to call this when using PETSc solver
	PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
	// disable PETSc error handler
	PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);
#endif

#ifndef COMPLEX
	if (argc < 3) die("Not enough parameters");
#else
	if (argc < 2) die("Not enough parameters");
#endif

	int n;
	Array<MatrixEntry> ar_mat;
	Array<scalar> ar_rhs;

	if (read_matrix_and_rhs(argv[2], n, ar_mat, ar_rhs) != ERR_SUCCESS)
		die("Failed to read the matrix and rhs.");

	if (strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
		PetscMatrix mat;
		PetscVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		PetscLinearSolver solver(&mat, &rhs);
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

		UMFPackLinearSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
		UMFPackMatrix mat;
		UMFPackVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		UMFPackLinearSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "pardiso") == 0) {
#ifdef WITH_PARDISO
		PardisoMatrix mat;
		PardisoVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		PardisoLinearSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "pardiso-block") == 0) {
#ifdef WITH_PARDISO
		PardisoMatrix mat;
		PardisoVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		PardisoLinearSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
		EpetraMatrix mat;
		EpetraVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		AztecOOSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "aztecoo-block") == 0) {
#ifdef WITH_TRILINOS
		EpetraMatrix mat;
		EpetraVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		AztecOOSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
		EpetraMatrix mat;
		EpetraVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		if (AmesosSolver::is_available("Klu")) {
			AmesosSolver solver("Klu", &mat, &rhs);
			solve(solver, n);
		}
#endif
	}
	else if (strcasecmp(argv[1], "amesos-block") == 0) {
#ifdef WITH_TRILINOS
		EpetraMatrix mat;
		EpetraVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		if (AmesosSolver::is_available("Klu")) {
			AmesosSolver solver("Klu", &mat, &rhs);
			solve(solver, n);
		}
#endif
	}
	else if (strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
		MumpsMatrix mat;
		MumpsVector rhs;
		build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

		MumpsSolver solver(&mat, &rhs);
		solve(solver, n);
#endif
	}
	else if (strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
		MumpsMatrix mat;
		MumpsVector rhs;
		build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

		MumpsSolver solver(&mat, &rhs);
		solve(solver, n);
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
