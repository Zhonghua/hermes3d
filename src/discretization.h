#ifndef _DISCRETIZATION_H_
#define _DISCRETIZATION_H_

#include <common/freqmap.h>
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

//	void set_space(int idx, Space *space);
	void set_spaces(int num, ...);
	void set_spaces(int num, Space **spaces);

//	void set_pss(int idx, PrecalcShapeset *pss);
	void set_pss(int num, ...);
	void set_pss(int num, PrecalcShapeset **pss);

//	void set_space_and_pss(int idx, Space *space, PrecalcShapeset *pss);

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

	//experimental function, to get info about type of each dof
	void static_condensation_info(char* file_name);

	void uncache_lsm(Word_t eid);

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
//	void allocate_matrices();
//	void assemble_matrices(PrecalcShapeset *ref_map_pss);

	void free_solution_vector();

	FreqMap<Word_t, scalar **> lsm_cache;
	FreqMap<Word_t, scalar **> lsm_cache_surf;
};

void update_limit_table(EMode3D mode);

#endif
