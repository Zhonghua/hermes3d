#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

#include "../linsolver.h"

struct msyst;

/// Encapsulation of PETSc linear solver
///
/// @ingroup linearsolvers
class PetscLinearSolver : public LinearSolver {
public:
	PetscLinearSolver();
	virtual ~PetscLinearSolver();

	virtual void prealloc(int ndofs);
	virtual void pre_add_ij(int row, int col);
	virtual void alloc();
	virtual void free();

	virtual void update_matrix(int m, int n, scalar v);
	virtual void update_matrix(int m, int n, double **mat, int *rows, int *cols);
	virtual void update_rhs(int idx, scalar y);
	virtual void update_rhs(int n, int *idx, scalar *y) ;

	virtual void begin_assembling();
	virtual void finish_assembling();

	virtual bool solve_system(double *sln);

	virtual void dump_matrix(FILE *file, EMatrixDumpFormat format = DF_MATLAB_SPARSE);
	virtual void dump_rhs(FILE *file, EMatrixDumpFormat format = DF_MATLAB_SPARSE);

	virtual ESparseMatrixRepresentation matrix_representation() { return SMATRIX_OTHER; }

protected:
	msyst *ms;

	bool initialized;
};

#endif

