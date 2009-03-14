// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2007 - 2009 Pavel Kus <pavel.kus@gmail.com>
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

#ifndef _DISCRETIZATION_H_
#define _DISCRETIZATION_H_

#include "linsolver.h"
#include "space.h"
#include "precalc.h"
#include "function.h"
#include "refmap.h"
#include "solution.h"

/// @defgroup assembling Assembling
///
///

/// Represents the discretization of the solved problem
///
/// TODO
/// @ingroup assembling
class Discretization {
public:
	Discretization(LinearSolver *lsolver);
	virtual ~Discretization();

	void free();

	//
	void set_num_equations(int neq);
	int get_num_equations() { return neq; }

	Space *get_space(int idx);
	PrecalcShapeset *get_pss(int idx);

	void set_spaces(int num, ...);
	void set_spaces(int num, Space **spaces);

	void set_pss(int num, ...);
	void set_pss(int num, PrecalcShapeset **pss);

	void set_bilinear_form(int i, int j,
		scalar (*bilinear_form_unsym)(RealFunction*, RealFunction*, RefMap*, RefMap*),
		scalar (*bilinear_form_sym)(RealFunction*, RealFunction*, RefMap*, RefMap*) = NULL,
		scalar (*bilinear_form_surf)(RealFunction*, RealFunction*, RefMap*, RefMap*, FacePos *) = NULL);

	void set_linear_form(int i,
		scalar (*linear_form)(RealFunction*, RefMap*),
		scalar (*linear_form_surf)(RealFunction*, RefMap*, FacePos *) = NULL);

	//
	void create_stiffness_matrix();
	void assemble_stiffness_matrix_and_rhs(bool rhsonly = false);
	bool solve_system(int n, ...);

protected:
	LinearSolver *linear_solver;	// linear solver
	int neq;						// number of equations
	int ndofs;						// number of unknowns
	Space **space;					// spaces
	PrecalcShapeset **pss;			// shapeset
	scalar *solution_vector;		// vector of the solution

	struct BiForm {
		scalar (*unsym)(RealFunction *, RealFunction *, RefMap *, RefMap *);
		scalar (*sym)(RealFunction *, RealFunction *, RefMap *, RefMap *);
		scalar (*surf)(RealFunction *, RealFunction *, RefMap *, RefMap *, FacePos *);

		BiForm() {
			unsym = NULL;
			sym = NULL;
			surf = NULL;
		}
	};

	struct LiForm {
		scalar (*lf)(RealFunction *, RefMap *);
		scalar (*surf)(RealFunction *, RefMap *, FacePos *);

		LiForm() {
			lf = NULL;
			surf = NULL;
		}
	};

	BiForm **biform;
	LiForm  *liform;

	void precalculate_sparse_structure(LinearSolver* solver);

	void free_solution_vector();
};

void update_limit_table(EMode3D mode);

#endif
