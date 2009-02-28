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

#include "config.h"
#include "common.h"
#include "norm.h"
#include "quad.h"
#include "discretization.h"
#include "refmap.h"
#include "integrals/h1.h"
#include "traverse.h"


/// Calculates the absolute error between sln1 and sln2 using function fn
double calc_error(double (*fn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*), MeshFunction *sln1, MeshFunction *sln2) {
	Mesh *meshes[2] = { sln1->get_mesh(), sln2->get_mesh() };
	Transformable *tr[2] = { sln1, sln2 };
	Traverse trav;
	trav.begin(2, meshes, tr);

	double error = 0.0;
	Element **ee;
	while ((ee = trav.get_next_state(NULL, NULL)) != NULL) {
		EMode3D mode = ee[0]->get_mode();

		update_limit_table(mode);
		Quad3D *quad = get_quadrature(mode);

		sln1->set_quad(quad);
		sln2->set_quad(quad);

		RefMap *ru = sln1->get_refmap();
		RefMap *rv = sln2->get_refmap();

		error += fn(sln1, sln2, ru, rv);
	}
	trav.finish();
	return error > 0.0 ? sqrt(error) : error;
}

/// Calculates the norm of sln using function fn
double calc_norm(double (*fn)(MeshFunction*, RefMap*), MeshFunction *sln) {
	double norm = 0.0;
	Element *e;
	Mesh *mesh = sln->get_mesh();

	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		Element *e = mesh->elements[eid];
		Quad3D *quad = get_quadrature(e->get_mode());

		sln->set_active_element(e);
		sln->set_quad(quad);
		RefMap *ru = sln->get_refmap();

		norm += fn(sln, ru);
	}

	return norm > 0.0 ? sqrt(norm) : norm;
}

// H1 space /////////////////////////////////////////////////////////////////////////////////////////

// function used to calculate error in H1 norm
double error_fn_h1(MeshFunction *sln1, MeshFunction *sln2, RefMap *ru, RefMap *rv) {
	Quad3D *quad = sln1->get_quad();

	order3_t o = max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
	o.limit();
	qorder_t qord = ELEM_QORDER(o);
	sln1->set_quad_order(qord);
	sln2->set_quad_order(qord);

	scalar *uval, *vval, *dudx, *dudy, *dudz, *dvdx, *dvdy, *dvdz;
	uval = sln1->get_fn_values();
	vval = sln2->get_fn_values();
	sln1->get_dx_dy_dz_values(dudx, dudy, dudz);
	sln2->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	H1_INTEGRATE_EXPRESSION(sqr(uval[i] - vval[i]) + sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]) + sqr(dudz[i] - dvdz[i]));
	return result;
}

// function used to calculate H1 norm of the solution
double norm_fn_h1(MeshFunction *sln, RefMap *ru) {
	Quad3D *quad = sln->get_quad();

	order3_t o = sln->get_fn_order() + ru->get_inv_ref_order();
	o.limit();
	qorder_t qord = ELEM_QORDER(o);
	sln->set_quad_order(qord);

	scalar *uval, *dudx, *dudy, *dudz;
	uval = sln->get_fn_values();
	sln->get_dx_dy_dz_values(dudx, dudy, dudz);

	H1_INTEGRATE_EXPRESSION(sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]) + sqr(dudz[i]));
	return result;
}


double h1_error(MeshFunction *sln1, MeshFunction *sln2) {
	double error = calc_error(error_fn_h1, sln1, sln2);
	double norm = calc_norm(norm_fn_h1, sln2);
	return error / norm;
}

double h1_norm(MeshFunction *sln) {
	return calc_norm(norm_fn_h1, sln);
}

// function used to calculate error in L2 norm
double error_fn_l2(MeshFunction *sln1, MeshFunction *sln2, RefMap *ru, RefMap *rv) {
	Quad3D *quad = sln1->get_quad();

	order3_t o = max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
	o.limit();
	qorder_t qord = ELEM_QORDER(o);
	sln1->set_quad_order(qord);
	sln2->set_quad_order(qord);

	scalar *uval, *vval;
	uval = sln1->get_fn_values();
	vval = sln2->get_fn_values();

	H1_INTEGRATE_EXPRESSION(sqr(uval[i] - vval[i]));
	return result;
}


// function used to calculate L2 norm of the solution
double norm_fn_l2(MeshFunction *sln, RefMap *ru) {
	Quad3D *quad = sln->get_quad();

	order3_t o = sln->get_fn_order() + ru->get_inv_ref_order();
	o.limit();
	qorder_t qord = ELEM_QORDER(o);
	sln->set_quad_order(qord);
	scalar *uval = sln->get_fn_values();

	H1_INTEGRATE_EXPRESSION(sqr(uval[i]));
	return result;
}


double l2_error(MeshFunction *sln1, MeshFunction *sln2) {
	double error = calc_error(error_fn_l2, sln1, sln2);
	double norm = calc_norm(norm_fn_l2, sln2);
	return error / norm;
}

double l2_norm(MeshFunction *sln) {
	return calc_norm(norm_fn_l2, sln);
}

