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
#include "matrix.h"
#include "forms.h"

#define FORM_CB(a)	a<fn_t, scalar>, a<fn_order_t, forder_t>

struct FcPos {
	int marker;			///< face marker
	int face;			///< element face number (local)

	// for internal use
	Element *base;
	Space *space, *space_u, *space_v;

};

/// @defgroup assembling Assembling
///
///

/// Represents the discretization of the solved problem
///
/// TODO: description
/// @ingroup assembling
class Discretization {
public:
	Discretization(int neq);
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
		scalar (*form)(int, double *, fn_t *, fn_t *, geom_t<double> *),
		forder_t (*order)(int, double *, fn_order_t *, fn_order_t *, geom_t<forder_t> *));
	void set_bilinear_form_surf(int i, int j,
		scalar (*form)(int, double *, fn_t *, fn_t *, FacePos *, geom_t<double> *),
		forder_t (*order)(int, double *, fn_order_t *, fn_order_t *, FacePos *, geom_t<forder_t> *));

	void set_linear_form(int i,
		scalar (*form)(int, double *, fn_t *, geom_t<double> *),
		forder_t (*order)(int, double *, fn_order_t *, geom_t<forder_t> *));
	void set_linear_form_surf(int i,
		scalar (*form)(int, double *, fn_t *, FacePos *, geom_t<double> *),
	    forder_t (*order)(int, double *, fn_order_t *, FacePos *, geom_t<forder_t> *));

	//
	void create(Matrix *matrix, Vector *rhs = NULL);
	void assemble(Matrix *matrix, Vector *rhs = NULL);

protected:
	int neq;						// number of equations
	int ndofs;						// number of unknowns

	struct BiForm {
		struct {
			scalar (*form)(int n, double *wt, fn_t *u, fn_t *v, geom_t<double> *e);
			forder_t (*order)(int n, double *wt, fn_order_t *u, fn_order_t *v, geom_t<forder_t> *e);
		} vol;

		struct {
			scalar (*form)(int n, double *wt, fn_t *u, fn_t *v, FacePos *fp, geom_t<double> *e);
			forder_t (*order)(int n, double *wt, fn_order_t *u, fn_order_t *v, FacePos *fp, geom_t<forder_t> *e);
		} surf;

		BiForm() {
			vol.form = NULL;
			vol.order = NULL;
			surf.form = NULL;
			surf.order = NULL;
		}
	};

	struct LiForm {
		struct {
			scalar (*form)(int n, double *wt, fn_t *v, geom_t<double> *e);
			forder_t (*order)(int n, double *wt, fn_order_t *v, geom_t<forder_t> *e);
		} vol;

		struct {
			scalar (*form)(int n, double *wt, fn_t *v, FacePos *fp, geom_t<double> *e);
			forder_t (*order)(int n, double *wt, fn_order_t *v, FacePos *fp, geom_t<forder_t> *e);
		} surf;

		LiForm() {
			vol.form = NULL;
			vol.order = NULL;
			surf.form = NULL;
			surf.order = NULL;
		}
	};

	scalar eval_bi_form(BiForm *bi, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv);
	scalar eval_bi_form_surf(BiForm *bi, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, FacePos *fp);
	scalar eval_li_form(LiForm *li, PrecalcShapeset *fv, RefMap *rv);
	scalar eval_li_form_surf(LiForm *li, PrecalcShapeset *fv, RefMap *rv, FacePos *fp);

	// Data
	Space **space;					// spaces
	PrecalcShapeset **pss;			// shapeset
	BiForm **biform;
	LiForm *liform;
};

void update_limit_table(EMode3D mode);

#endif
