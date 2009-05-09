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

#include <stdarg.h>
#include "h3dconfig.h"
#include "discretization.h"
#include "traverse.h"
#include "forms.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>
#include "common.h"

void update_limit_table(EMode3D mode) {
	_F_
	// do nothing
}

// Discretization /////////////////////////////////////////////////////////////

Discretization::Discretization(int neq) {
	_F_

	if (neq <= 0) ERROR("Invalid number of equations.");
	if (neq > 10) WARNING("Large number of equations (%d). Is this the intent?", neq);

	// initialize the bilinear form
	biform = new BiForm *[neq];
	MEM_CHECK(biform);
	for (int i = 0; i < neq; i++) {
		biform[i] = new BiForm[neq];
		MEM_CHECK(biform[i]);
	}

	// init the rest of the arrays
	liform = new LiForm[neq];
	MEM_CHECK(liform);

	space = new Space *[neq];
	MEM_CHECK(space);
	pss = new PrecalcShapeset *[neq];
	MEM_CHECK(pss);

	for (int j = 0; j < neq; j++) {
		space[j] = NULL;
		pss[j] = NULL;
	}

	this->neq = neq;
}

Discretization::~Discretization() {
	_F_
	free();
}

void Discretization::free() {
	_F_
	delete [] space; space = NULL;
	delete [] pss; pss = NULL;

	for (int i = 0; i < neq; i++)
		delete [] biform[i];
	delete [] biform; biform = NULL;
	delete [] liform; liform = NULL;
}

void Discretization::set_num_equations(int neq) {
	_F_
	if (neq <= 0) ERROR("Invalid number of equations.");
	if (neq > 10) WARNING("Large number of equations (%d). Is this the intent?", neq);
	free();

	// initialize the bilinear form
	biform = new BiForm *[neq];
	MEM_CHECK(biform);
	for (int i = 0; i < neq; i++) {
		biform[i] = new BiForm[neq];
		MEM_CHECK(biform[i]);
	}

	// init the rest of the arrays
	liform = new LiForm[neq];
	MEM_CHECK(liform);
	space = new Space *[neq];
	MEM_CHECK(space);
	pss = new PrecalcShapeset *[neq];
	MEM_CHECK(pss);

	for (int j = 0; j < neq; j++) {
		space[j] = NULL;
		pss[j] = NULL;
	}

	this->neq = neq;
}

Space *Discretization::get_space(int idx) {
	assert(0 <= idx && idx < this->neq);
	return space[idx];
}

PrecalcShapeset *Discretization::get_pss(int idx) {
	assert(0 <= idx && idx < this->neq);
	return pss[idx];
}

void Discretization::set_spaces(int num, ...) {
	_F_
	assert(0 <= num && num <= neq);
	va_list ap;
	va_start(ap, num);
	for (int i = 0; i < num; i++)
		space[i] = va_arg(ap, Space *);
	va_end(ap);
}

void Discretization::set_spaces(int num, Space **s) {
	_F_
	assert(0 <= num && num <= neq);
	for (int i = 0; i < num; i++)
		space[i] = s[i];
}

void Discretization::set_pss(int num, ...) {
	_F_
	assert(0 <= num && num <= neq);
	va_list ap;
	va_start(ap, num);
	for (int i = 0; i < num; i++)
		this->pss[i] = va_arg(ap, PrecalcShapeset *);
	va_end(ap);
}

void Discretization::set_pss(int num, PrecalcShapeset **pss) {
	_F_
	assert(0 <= num && num <= neq);
	for (int i = 0; i < num; i++)
		this->pss[i] = pss[i];
}

void Discretization::set_bilinear_form(int i, int j,
	scalar (*form)(int, double *, fn_t *, fn_t *, geom_t<double> *),
	forder_t (*order)(int, double *, fn_order_t *, fn_order_t *, geom_t<forder_t> *))
{
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) ERROR("Bad equation number.");

	biform[i][j].vol.form = form;
	biform[i][j].vol.order = order;
}

void Discretization::set_bilinear_form_surf(int i, int j,
	scalar (*form)(int, double *, fn_t *, fn_t *, FacePos *, geom_t<double> *),
	forder_t (*order)(int, double *, fn_order_t *, fn_order_t *, FacePos *, geom_t<forder_t> *))
{
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) ERROR("Bad equation number.");

	biform[i][j].surf.form = form;
	biform[i][j].surf.order = order;
}

void Discretization::set_linear_form(int i,
	scalar (*form)(int, double *, fn_t *, geom_t<double> *),
	forder_t (*order)(int, double *, fn_order_t *, geom_t<forder_t> *))
{
	_F_
	if (i < 0 || i >= neq) ERROR("Bad equation number.");

	liform[i].vol.form = form;
	liform[i].vol.order = order;
}

void Discretization::set_linear_form_surf(int i,
	scalar (*form)(int, double *, fn_t *, FacePos *, geom_t<double> *),
    forder_t (*order)(int, double *, fn_order_t *, FacePos *, geom_t<forder_t> *))
{
	_F_
	if (i < 0 || i >= neq) ERROR("Bad equation number.");

	liform[i].surf.form = form;
	liform[i].surf.order = order;
}

void Discretization::create(Matrix *matrix, Vector *rhs) {
	_F_
	// calculate the total number of DOFs
	ndofs = 0;
	for (int i = 0; i < neq; i++)
		ndofs += space[i]->get_dof_count();
	if (ndofs == 0) return;

	if (matrix != NULL) {
		SparseMatrix *mat = dynamic_cast<SparseMatrix *>(matrix);
		if (mat == NULL) {
			WARNING("Calling Discretization::create() with a pointer to matrix that is not sparse matrix.");
			return;
		}

		mat->free();
		mat->prealloc(ndofs);

		AsmList *al = new AsmList[neq];
		MEM_CHECK(al);

		// init multi-mesh traversal
		Mesh **meshes = new Mesh *[neq];
		MEM_CHECK(meshes);
		for (int i = 0; i < neq; i++)
			meshes[i] = space[i]->get_mesh();
		Traverse trav;
		trav.begin(neq, meshes);

		Element **e;
		while ((e = trav.get_next_state(NULL, NULL)) != NULL) {
			// obtain assembly lists for the element at all spaces
			for (int i = 0; i < neq; i++)
				space[i]->get_element_assembly_list(e[i], al + i);

			// go through all equation-blocks of the local stiffness matrix
			for (int m = 0; m < neq; m++) {
				for (int n = 0; n < neq; n++) {
					BiForm *bf = biform[m] + n;
					if (bf->vol.form == NULL && bf->surf.form == NULL) continue;

					// pretend assembling of the element stiffness matrix
					for (int j = 0; j < al[n].cnt; j++) {
						// skip dirichlet dofs in 'j'
						if (al[n].dof[j] < 0) continue;

						for (int i = 0; i < al[m].cnt; i++) {
							// skip dirichlet dofs in 'i'
							if (al[m].dof[i] < 0) continue;

							// register the corresponding nonzero matrix element
							mat->pre_add_ij(al[m].dof[i], al[n].dof[j]);
						}
					}
				}
			}
		}

		trav.finish();
		delete [] meshes;
		delete [] al;

		mat->alloc();
	}

	if (rhs != NULL) rhs->alloc(ndofs);
}

void Discretization::assemble(Matrix *matrix, Vector *rhs) {
	_F_
	if (ndofs == 0) return;

	bool bnd[10];
	FacePos fp[10];					// FIXME: magic number - number of faces

	bool *nat = new bool[neq];
	MEM_CHECK(nat);

	// initialize assembly lists
	AsmList *al = new AsmList[neq];
	MEM_CHECK(al);

	// initialize refmaps
	RefMap *refmap = new RefMap[neq];
	MEM_CHECK(refmap);
	for (int i = 0; i < neq; i++) {
		refmap[i].set_mesh(space[i]->get_mesh());
	}

	// create slave pss's for test functions, init quadrature points
	PrecalcShapeset **spss = new PrecalcShapeset *[neq];
	MEM_CHECK(spss);
	for (int i = 0; i < neq; i++) {
		spss[i] = new PrecalcShapeset(pss[i]);
		MEM_CHECK(spss[i]);
	}

	// init multi-mesh traversal
	int nm = neq;
	Mesh *meshes[nm];
	Transformable *fn[nm];
	for (int i = 0; i < neq; i++)
		meshes[i] = space[i]->get_mesh();
	memcpy(fn, pss, neq * sizeof(Transformable *));

	// loop through all elements
	Element **e;
	Traverse trav;
	trav.begin(nm, meshes, fn);
	while ((e = trav.get_next_state(bnd, fp)) != NULL) {
		for (int eq = 0; eq < neq; eq++) {
			pss[eq]->set_quad(get_quadrature(e[eq]->get_mode()));
			spss[eq]->set_quad(get_quadrature(e[eq]->get_mode()));

			// obtain assembly lists for the element at all spaces, set appropriate mode for each pss
			space[eq]->get_element_assembly_list(e[eq], al + eq);

			spss[eq]->set_active_element(e[eq]);
			spss[eq]->set_master_transform();

		    refmap[eq].set_active_element(e[eq]);
			refmap[eq].force_transform(pss[eq]->get_transform(), pss[eq]->get_ctm());
		}

		// go through all equation-blocks of the element stiffness matrix, assemble volume integrals
		AsmList *am = al;
		for (int m = 0; m < neq; m++, am++) {
			PrecalcShapeset *fv = spss[m];

			scalar *lrhs = new scalar[am->cnt];			// local right-hand side
			MEM_CHECK(lrhs);
			memset(lrhs, 0, am->cnt * sizeof(scalar));

			if (matrix != NULL) {
				AsmList *an = al;
				for (int n = 0; n < neq; n++, an++) {
					PrecalcShapeset *fu = pss[n];

					BiForm *bf = biform[m] + n;
					if (bf->vol.form == NULL) continue;

					// assemble the (m,n)-block of the stiffness matrix, one column at a time
					scalar **lsm = new_matrix<scalar>(am->cnt, an->cnt);			// local stiffness matrix

					for (int j = 0; j < an->cnt; j++) {
						fu->set_active_shape(an->idx[j]);
						int l = an->dof[j];

						for (int i = 0; i < am->cnt; i++) {
							if (am->dof[i] < 0) continue;
							fv->set_active_shape(am->idx[i]);
							scalar bi = eval_bi_form(bf, fu, fv, refmap + n, refmap + m) * an->coef[j] * am->coef[i];
							if (l >= 0) lsm[i][j] += bi;
							else lrhs[i] -= bi;
						}
					}

					// update the matrix
					matrix->update(am->cnt, an->cnt, lsm, am->dof, an->dof);
					delete [] lsm;
				}
			}

			// assemble rhs (linear form)
			if (rhs != NULL) {
				if (liform[m].vol.form != NULL) {
					for (int i = 0; i < am->cnt; i++) {
						if (am->dof[i] >= 0) {
							fv->set_active_shape(am->idx[i]);
							lrhs[i] += eval_li_form(liform + m, fv, refmap + m) * am->coef[i];
						}
					}

					// update rhs
					rhs->update(am->cnt, am->dof, lrhs);
				}
			}

			delete [] lrhs;
		}

		// assemble surface integrals now: loop through boundary edges of the element
		if (matrix == NULL) continue;

		for (int iface = 0; iface < e[0]->get_num_of_faces(); iface++) {
			// check if the face is an outer face
			if (!bnd[iface]) continue;

			// obtain the list of shape functions which are nonzero on this edge
			for (int eq = 0; eq < neq; eq++)
				if ((nat[eq] = (space[eq]->bc_type_callback(fp[iface].marker) == BC_NATURAL)))
					space[eq]->get_boundary_assembly_list(e[eq], iface, al + eq);

			// loop through the equation-blocks
			AsmList *am = al;
			for (int m = 0; m < neq; m++, am++) {
				if (!nat[m]) continue;

				scalar *lrhs = new scalar[am->cnt];			// local right-hand side
				MEM_CHECK(lrhs);
				memset(lrhs, 0, am->cnt * sizeof(scalar));

				PrecalcShapeset *fv = spss[m];
				fp[iface].base = trav.get_base();
				fp[iface].space_v = space[m];
				AsmList *an = al;
				for (int n = 0; n < neq; n++, an++) {
					if (!nat[n]) continue;
					BiForm *bf = biform[m] + n;
					if (bf->surf.form == NULL) continue;
					PrecalcShapeset *fu = pss[n];
					fp[iface].space_u = space[n];

					// assemble the surface part of the bilinear form
					scalar **lsm = new_matrix<scalar>(am->cnt, an->cnt);			// local stiffness matrix
					for (int j = 0; j < an->cnt; j++) {
						fu->set_active_shape(an->idx[j]);
						int l = an->dof[j];

						for (int i = 0; i < am->cnt; i++) {
							int k = am->dof[i];
							if (k < 0) continue;
							fv->set_active_shape(am->idx[i]);
							scalar bi = eval_bi_form_surf(bf, fu, fv, refmap + n, refmap + m, fp + iface) * an->coef[j] * am->coef[i];
							if (l >= 0) lsm[i][j] += bi;
							else lrhs[i] -= bi;
						}
					}

					// update matrix
					matrix->update(am->cnt, an->cnt, lsm, am->dof, an->dof);

					delete [] lsm;
				}

				// assemble the surface part of the linear form
				if (liform[m].surf.form == NULL) continue;
				if (rhs != NULL) {
					for (int i = 0; i < am->cnt; i++) {
						if (am->dof[i] < 0) continue;
						fv->set_active_shape(am->idx[i]);
						lrhs[i] += eval_li_form_surf(liform + m, fv, refmap + m, fp + iface) * am->coef[i];
					}
					rhs->update(am->cnt, am->dof, lrhs);
				}
				delete [] lrhs;
			}
		}
	}
	trav.finish();

	// free memory
	for (int i = 0; i < neq; i++)
		delete spss[i];
	delete [] spss;
	delete [] refmap;
	delete [] al;
	delete [] nat;
}

//

scalar Discretization::eval_bi_form(BiForm *bi, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv) {
	_F_

	double fake_wt = 1.0;
	geom_t<forder_t> fake_e;
	init_geom(fake_e);

	fn_order_t ou(fu->get_fn_order());
	fn_order_t ov(fv->get_fn_order());
	forder_t o = bi->vol.order(1, &fake_wt, &ou, &ov, &fake_e);

	order3_t order = o.get_order() + ru->get_inv_ref_order();
	order.limit();

	int np;
	double *jwt = NULL;
	geom_t<double> e;
	fn_t u, v;

	init_jwt(fu, ru, order, np, jwt);
	init_geom(e, ru, order);
	init_fn(fu, ru, order, u);
	init_fn(fv, rv, order, v);
	scalar result = bi->vol.form(np, jwt, &u, &v, &e);

	delete [] jwt;

	return result;
}

scalar Discretization::eval_bi_form_surf(BiForm *bi, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	_F_

	double fake_wt = 1.0;
	geom_t<forder_t> fake_e;
	init_geom(fake_e);
//	FcPos<forder_t> ffp(fp);


	fn_order_t ou(fu->get_fn_order());
	fn_order_t ov(fv->get_fn_order());
	forder_t o = bi->surf.order(1, &fake_wt, &ou, &ov, fp, &fake_e);

	order3_t vol_order = o.get_order() + ru->get_inv_ref_order();
	vol_order.limit();
	order2_t order = vol_order.get_face_order(fp->face);

	int np;
	double *jwt = NULL;
	geom_t<double> e;
	fn_t u, v;

	init_jwt(fu, ru, fp->face, order, np, jwt);
	init_geom(e, ru, fp->face, order);
	init_fn(fu, ru, fp->face, order, u);
	init_fn(fv, rv, fp->face, order, v);
	scalar result = bi->surf.form(np, jwt, &u, &v, fp, &e);

	delete [] jwt;

	return result;
}

scalar Discretization::eval_li_form(LiForm *li, PrecalcShapeset *fv, RefMap *rv) {
	_F_

	double fake_wt = 1.0;
	geom_t<forder_t> fake_e;
	init_geom(fake_e);

	fn_order_t ov(fv->get_fn_order());
	forder_t o = li->vol.order(1, &fake_wt, &ov, &fake_e);
	order3_t order = o.get_order() + rv->get_inv_ref_order();
	order.limit();

	int np = 0;
	double *jwt = NULL;
	geom_t<double> e;
	fn_t v;

	init_jwt(fv, rv, order, np, jwt);
	init_geom(e, rv, order);
	init_fn(fv, rv, order, v);
	scalar result = li->vol.form(np, jwt, &v, &e);

	delete [] jwt;

	return result;
}

scalar Discretization::eval_li_form_surf(LiForm *li, PrecalcShapeset *fv, RefMap *rv, FacePos *fp) {
	_F_

	double fake_wt = 1.0;
	geom_t<forder_t> fake_e;
	init_geom(fake_e);

	fn_order_t ov(fv->get_fn_order());
	forder_t o = li->surf.order(1, &fake_wt, &ov, fp, &fake_e);
	order3_t vol_order = o.get_order() + rv->get_inv_ref_order();
	vol_order.limit();
	order2_t order = vol_order.get_face_order(fp->face);

	int np = 0;
	double *jwt = NULL;
	geom_t<double> e;
	fn_t v;

	init_jwt(fv, rv, fp->face, order, np, jwt);
	init_geom(e, rv, fp->face, order);
	init_fn(fv, rv, fp->face, order, v);
	scalar result = li->surf.form(np, jwt, &v, fp, &e);

	delete [] jwt;

	return result;
}
