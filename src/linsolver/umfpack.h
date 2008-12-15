#ifndef _UMFPACK_SOLVER_H_
#define _UMFPACK_SOLVER_H_

#include "../linsolver.h"

/// Encapsulation of UMFPACK linear solver
///
/// @ingroup linearsolvers
class UMFPackLinearSolver : public LinearSolver {
public:
	UMFPackLinearSolver();
	virtual ~UMFPackLinearSolver();

	virtual void prealloc(int ndofs);
	virtual void pre_add_ij(int row, int col);
	virtual void alloc();
	virtual void free();

	virtual void update_matrix(int m, int n, scalar v);
	virtual void update_matrix(int m, int n, scalar **mat, int *rows, int *cols);
	virtual void update_rhs(int idx, scalar y);
	virtual void update_rhs(int n, int *idx, scalar *y);

	virtual bool solve_system(scalar *sln);

	virtual bool dump_matrix(FILE *file, const char *var_name, EMatrixDumpFormat = DF_MATLAB_SPARSE);
	virtual bool dump_rhs(FILE *file, const char *var_name, EMatrixDumpFormat = DF_MATLAB_SPARSE);

	virtual ESparseMatrixRepresentation matrix_representation() { return SMATRIX_COLUMN; }

	virtual int get_matrix_size() const;

protected:
	// UMFPack specific data structures for storing matrix, rhs
	int *Ap;
	int *Ai;
	scalar *Ax;
	scalar *srhs;

	// mem stat
	int mem_size;
};

#endif
