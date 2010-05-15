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

#ifndef _SOLVER_EPETRA_H_
#define _SOLVER_EPETRA_H_

#include "../h3dconfig.h"
#include "../solver.h"
#include "../matrix.h"

#ifdef HAVE_EPETRA
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#endif


class EpetraMatrix : public SparseMatrix {
public:
	EpetraMatrix();
#ifdef HAVE_EPETRA
	EpetraMatrix(Epetra_RowMatrix &mat);
#endif
	virtual ~EpetraMatrix();

	virtual void prealloc(int n);
	virtual void pre_add_ij(int row, int col);
	virtual void finish();

	virtual void alloc();
	virtual void free();
	virtual scalar get(int m, int n);
	virtual int get_num_row_entries(int row);
	virtual void extract_row_copy(int row, int len, int &n_entries, double *vals, int *idxs);
	virtual void zero();
	virtual void add(int m, int n, scalar v);
	virtual void add(int m, int n, scalar **mat, int *rows, int *cols);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
	virtual int get_matrix_size() const;
	virtual double get_fill_in() const;

protected:
#ifdef HAVE_EPETRA
	Epetra_BlockMap *std_map;
	Epetra_CrsGraph *grph;
	Epetra_CrsMatrix *mat;
#ifdef COMPLEX
	Epetra_CrsMatrix *mat_im;		// imaginary part of the matrix, mat holds the real part
#endif
	bool owner;
#endif

	friend class AmesosSolver;
	friend class AztecOOSolver;
	friend class NoxSolver;
	friend class IfpackPrecond;
	friend class MlPrecond;
};


class EpetraVector : public Vector {
public:
	EpetraVector();
#ifdef HAVE_EPETRA
	EpetraVector(const Epetra_Vector &v);
#endif
	virtual ~EpetraVector();

	virtual void alloc(int ndofs);
	virtual void free();
#ifdef HAVE_EPETRA
	virtual scalar get(int idx) { return (*vec)[idx]; }
#ifndef COMPLEX
	virtual void extract(scalar *v) const { vec->ExtractCopy(v); }
#else
	virtual void extract(scalar *v) const { }
#endif
#else
	virtual scalar get(int idx) { return 0.0; }
	virtual void extract(scalar *v) const { }
#endif
	virtual void zero();
	virtual void set(int idx, scalar y);
	virtual void add(int idx, scalar y);
	virtual void add(int n, int *idx, scalar *y);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
#ifdef HAVE_EPETRA
	Epetra_BlockMap *std_map;
	Epetra_Vector *vec;
#ifdef COMPLEX
	Epetra_Vector *vec_im;		// imaginary part of the vector, vec holds the real part
#endif
	bool owner;
#endif

	friend class AmesosSolver;
	friend class AztecOOSolver;
	friend class NoxSolver;
};

#endif
