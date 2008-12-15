#ifndef PARDISOSOLVER_H_
#define PARDISOSOLVER_H_

#include "../linsolver.h"

/// Encapsulation of pardiso linear solver
///
/// @ingroup linearsolvers
class PardisoLinearSolver : public LinearSolver {
public:
	PardisoLinearSolver();
	virtual ~PardisoLinearSolver();

	virtual void prealloc(int ndofs);
	virtual void pre_add_ij(int row, int col);
	virtual void alloc();
	virtual void free();

	virtual void update_matrix(int m, int n, scalar v);
	virtual void update_matrix(int m, int n, double **mat, int *rows, int *cols);
	virtual void update_rhs(int idx, scalar y);
	virtual void update_rhs(int n, int *idx, scalar *y) ;

	virtual bool solve_system(double *sln);

	virtual bool dump_matrix(FILE *file, const char *var_name, EMatrixDumpFormat format = DF_MATLAB_SPARSE);
	virtual bool dump_rhs(FILE *file, const char *var_name, EMatrixDumpFormat format = DF_MATLAB_SPARSE);

	virtual ESparseMatrixRepresentation matrix_representation() { return SMATRIX_ROW; }

	virtual int get_matrix_size() const;

protected:
	// Pardiso specific data structures for storing matrix, rhs
	int *Ap;
	int *Ai;
	scalar *Ax;
	scalar *srhs;
};

#endif /*PARDISOSOLVER_H_*/
