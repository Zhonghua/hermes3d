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

#include "h3dconfig.h"
#include "common.h"
#include "solution.h"
#include "function.cc" // non-inline template members
#include <common/error.h>
#include <common/callstack.h>

//// MeshFunction //////////////////////////////////////////////////////////////////////////////////

MeshFunction::MeshFunction(Mesh *mesh) :
	ScalarFunction()
{
	_F_
	this->mesh = mesh;
	this->refmap = new RefMap(mesh);
	MEM_CHECK(this->refmap);
	this->element = NULL;		// this comes with Transformable
	this->seq = 0;
	this->noinc = false;
}

MeshFunction::~MeshFunction() {
	_F_
	delete refmap;
}

void MeshFunction::set_quad(Quad3D *quad) {
	_F_
	ScalarFunction::set_quad(quad);
	refmap->set_quad(quad);
}

void MeshFunction::set_active_element(Element *e) {
	_F_
	element = e;
	mode = e->get_mode();
	refmap->set_active_element(e);
	reset_transform();
}

//// Solution //////////////////////////////////////////////////////////////////////////////////////

Solution::Solution(Mesh *mesh) : MeshFunction(mesh) {
	_F_
	memset(tables, 0, sizeof(tables));
	memset(elems, 0, sizeof(elems));
	memset(oldest, 0, sizeof(oldest));
	dir_coef = 1.0;
	vec = NULL;
	owner = false;
	transform = true;
	space = NULL;
	pss = NULL;
	slave_pss = NULL;
}

Solution::~Solution() {
	_F_
	free();
	if (slave_pss != NULL) delete slave_pss;
}

void Solution::free() {
	_F_
	free_tables();

	if (owner && vec) {
		delete[] vec;
		vec = NULL;
		owner = false;
	}
}

void Solution::free_tables() {
	_F_
	for (int i = 0; i < QUAD_COUNT; i++)
		for (int j = 0; j < NUM_ELEMENTS; j++)
			free_sub_tables(&(tables[i][j]));
}

void Solution::set_space_and_pss(Space *space, PrecalcShapeset *pss) {
	_F_
	if (space->get_shapeset() != pss->get_shapeset()) ERROR("'space' and 'pss' must have the same shapesets.");

	this->space = space;
	this->pss = pss;
	this->mesh = space->get_mesh();
	if (slave_pss != NULL) delete slave_pss;
	slave_pss = new PrecalcShapeset(pss);
	MEM_CHECK(slave_pss);
	num_components = pss->num_components;
}

void Solution::set_solution_vector(scalar *vec, bool owner, scalar dir_coef) {
	_F_
	free();
	this->dir_coef = dir_coef;
	this->vec = vec;
	this->owner = owner;
}

void Solution::set_zero_vector() {
	_F_
	free();
	int ndofs = space->get_max_dof() + 1;
	vec = new scalar[ndofs];
	MEM_CHECK(vec);
	memset(vec, 0, sizeof(scalar) * ndofs);
	owner = true;
}

void Solution::set_quad(Quad3D *quad) {
	_F_
	MeshFunction::set_quad(quad);
	slave_pss->set_quad(quad);
}

void Solution::set_active_element(Element *e) {
	_F_
	MeshFunction::set_active_element(e);

	// try finding an existing table for e
	for (cur_elem = 0; cur_elem < NUM_ELEMENTS; cur_elem++)
		if (elems[cur_quad][cur_elem] == e) break;

	// if not found, free the oldest one and use its slot
	if (cur_elem >= NUM_ELEMENTS) {
		if (tables[cur_quad][oldest[cur_quad]] != NULL) free_sub_tables(&(tables[cur_quad][oldest[cur_quad]]));

		cur_elem = oldest[cur_quad];
		if (++oldest[cur_quad] >= NUM_ELEMENTS) oldest[cur_quad] = 0;

		elems[cur_quad][cur_elem] = e;
	}

	mode = e->get_mode();

	assert(space != NULL);
	space->get_element_assembly_list(e, &(al[cur_elem]));
	assert(slave_pss != NULL);
	slave_pss->set_active_element(e);

	sub_tables = &(tables[cur_quad][cur_elem]);
	update_nodes_ptr();

	order = space->get_element_order(e->id);
}

void Solution::precalculate(qorder_t qord, int mask) {
	_F_
	int i, j, k, l;

	// if we are required to transform vectors, we must precalculate their components
	const int GRAD = FN_DX_0 | FN_DY_0 | FN_DZ_0;
	const int GRAD_ALL = FN_DX | FN_DY | FN_DZ;
	if (transform) {
		if (num_components == 1 && ((mask & FN_DX_0) || (mask & FN_DY_0) || (mask & FN_DZ_0)))
			mask |= GRAD;
		else if (num_components == 3 && ((mask & FN_VAL_0) || (mask & FN_VAL_1) || (mask & FN_VAL_2)))
			mask |= FN_VAL;
		if (num_components == 3 && (
				(mask & FN_DX_0) || (mask & FN_DY_0) || (mask & FN_DZ_0) ||
				(mask & FN_DX_1) || (mask & FN_DY_1) || (mask & FN_DZ_1) ||
				(mask & FN_DX_2) || (mask & FN_DY_2) || (mask & FN_DZ_2)))
			mask |= GRAD_ALL;
	}

	int oldmask = (cur_node != NULL) ? cur_node->mask : 0;
	int newmask = mask | oldmask;

	Quad3D *quad = quads[cur_quad];
	assert(quad != NULL);
	QuadPt3D *pt = NULL;
	int np = 0;
	switch (qord.type) {
		case QOT_ELEMENT:
			np = quad->get_num_points(order3_t::from_int(qord.order));
			pt = quad->get_points(order3_t::from_int(qord.order));
			break;

		case QOT_FACE:
			np = quad->get_face_num_points(qord.face, order2_t::from_int(qord.order));
			pt = quad->get_face_points(qord.face, order2_t::from_int(qord.order));
			break;

		case QOT_EDGE:
			np = quad->get_edge_num_points(qord.order);
			pt = quad->get_edge_points(qord.edge, qord.order);
			break;

		case QOT_VERTEX:
			np = quad->get_vertex_num_points();
			pt = quad->get_vertex_points();
			break;

		default: assert(false);
	}

	Node *node = new_node(newmask, np);
	MEM_CHECK(node);

	// transform integration points by the current matrix
	scalar x[np], y[np], z[np];
	for (i = 0; i < np; i++) {
		x[i] = pt[i].x * ctm->m[0] + ctm->t[0];
		y[i] = pt[i].y * ctm->m[1] + ctm->t[1];
		z[i] = pt[i].z * ctm->m[2] + ctm->t[2];
	}

	// reuse old tables, zero new tables
	for (j = 0; j < num_components; j++)
		for (i = 0; i < VALUE_TYPES; i++)
			if (newmask & idx2mask[i][j]) {
				if (oldmask & idx2mask[i][j]) memcpy(node->values[j][i], cur_node->values[j][i], np * sizeof(scalar));
				else memset(node->values[j][i], 0, np * sizeof(scalar));
			}

	// update ctm, force it to the slave pss
	slave_pss->force_transform(sub_idx, ctm);
	// assemble linear combination of the solution
	AsmList *pal = al + cur_elem;
	for (k = 0; k < pal->cnt; k++) {
		slave_pss->set_active_shape(pal->idx[k]);
		slave_pss->set_quad_order(qord, mask);
		scalar coef = pal->coef[k] * (pal->dof[k] >= 0 ? vec[pal->dof[k]] : dir_coef);

		for (j = 0; j < num_components; j++)
			for (l = 0; l < VALUE_TYPES; l++)
				if ((newmask & idx2mask[l][j]) && !(oldmask & idx2mask[l][j])) {
					scalar *table = node->values[j][l];
					double *val = slave_pss->get_values(j, l);
					for (i = 0; i < np; i++) {
						table[i] += val[i] * coef;
					}
				}
	}

	// transform gradient or vector solution
	if (transform) {
		bool trans = false;
		bool trans_hcurl_grad = false; //it has to be transformed as well!
		scalar *tab1, *tab2, *tab3;
		const int GRAD = FN_DX_0 | FN_DY_0 | FN_DZ_0;
		if (num_components == 1 && (newmask & GRAD) == GRAD && (oldmask & GRAD) != GRAD) {
			trans = true;
			tab1 = node->values[0][DX];
			tab2 = node->values[0][DY];
			tab3 = node->values[0][DZ];
		}
		else if (num_components == 3 && (newmask & FN_VAL) == FN_VAL && (oldmask & FN_VAL) != FN_VAL) {
			trans = true;
			tab1 = node->values[0][FN];
			tab2 = node->values[1][FN];
			tab3 = node->values[2][FN];
		}
		if (num_components == 3 && (newmask & GRAD) == GRAD && (oldmask & GRAD) != GRAD) {
			trans_hcurl_grad = true;
		}

		if (trans) {
			double3x3 *m = NULL;
			int mstep;

	        update_refmap();
			if (refmap->is_jacobian_const()) {
				// only one inverse ref. map
				m = refmap->get_const_inv_ref_map();
				mstep = 0;
			}
			else {
				// array of inverse ref. maps
				mstep = 1;
				switch (qord.type) {
					case QOT_ELEMENT: m = refmap->get_inv_ref_map(order3_t::from_int(qord.order)); break;
					case QOT_FACE: m = refmap->get_face_inv_ref_map(qord.face, order2_t::from_int(qord.order)); break;
					case QOT_EDGE: m = refmap->get_edge_inv_ref_map(qord.face, order1_t(qord.order)); break;
					case QOT_VERTEX: m = refmap->get_vertex_inv_ref_map(); break;
					default: assert(false);
				}

			}

			for (i = 0; i < np; i++, m += mstep) {
				scalar vx = tab1[i], vy = tab2[i], vz = tab3[i];
				tab1[i] = (*m)[0][0] * vx + (*m)[0][1] * vy + (*m)[0][2] * vz;
				tab2[i] = (*m)[1][0] * vx + (*m)[1][1] * vy + (*m)[1][2] * vz;
				tab3[i] = (*m)[2][0] * vx + (*m)[2][1] * vy + (*m)[2][2] * vz;
			}
		}

		if (trans_hcurl_grad) {
			double3x3 *m;
			int mstep;
			if (refmap->is_jacobian_const()) {
				// only one inverse ref. map
				m = refmap->get_const_inv_ref_map();
				mstep = 0;
			}
			else {
				// array of inverse ref. maps
				m = refmap->get_inv_ref_map(order3_t::from_int(qord.order));
				mstep = 1;
			}

			scalar vx[3], vy[3], vz[3], vhx[3], vhy[3], vhz[3];

			for (i = 0; i < np; i++, m += mstep) {
				for(int c = 0; c < 3; c++){
					vx[c] = node->values[c][DX][i];
					vy[c] = node->values[c][DY][i];
					vz[c] = node->values[c][DZ][i];
				}
				for(int c = 0; c < 3; c++){
					vhx[c] = (*m)[0][0] * vx[c] + (*m)[0][1] * vy[c] + (*m)[0][2] * vz[c];
					vhy[c] = (*m)[1][0] * vx[c] + (*m)[1][1] * vy[c] + (*m)[1][2] * vz[c];
					vhz[c] = (*m)[2][0] * vx[c] + (*m)[2][1] * vy[c] + (*m)[2][2] * vz[c];
				}
				for(int c = 0; c < 3; c++){
					node->values[c][DX][i] = (*m)[c][0] * vhx[0] + (*m)[c][1] * vhx[1] + (*m)[c][2] * vhx[2];
					node->values[c][DY][i] = (*m)[c][0] * vhy[0] + (*m)[c][1] * vhy[1] + (*m)[c][2] * vhy[2];
					node->values[c][DZ][i] = (*m)[c][0] * vhz[0] + (*m)[c][1] * vhz[1] + (*m)[c][2] * vhz[2];
				}
			}
		}
	}

	// remove the old node and attach the new one
	replace_cur_node(node);
}

scalar Solution::get_sln_value(double x, double y, double z, EValueType which, int component) {
	_F_
	Shapeset *shapeset = slave_pss->get_shapeset();
	AsmList *pal = al + cur_elem;
	scalar result = 0.0;
	for (int k = 0; k < pal->cnt; k++) {
		scalar coef = pal->coef[k] * (pal->dof[k] >= 0 ? vec[pal->dof[k]] : dir_coef);
		result += coef * shapeset->get_value(which, pal->idx[k], x, y, z, component);
	}
	return result;
}

void Solution::save_solution_vector(char *filename, int ndofs) {
	_F_
	FILE *f = fopen(filename, "wb");
	if (f == NULL) ERROR("Cannot open %s for writing.", filename);
	fwrite(vec, sizeof(scalar), ndofs, f);
	fclose(f);
}


void Solution::load_solution_vector(char *filename, int ndofs) {
	_F_
	// TODO: check that set_space_and_pss has been called
	FILE *f = fopen(filename, "rb");
	if (f == NULL) ERROR("Cannot open %s.", filename);
	free();
	vec = new scalar[ndofs];
	MEM_CHECK(vec);
	owner = true;
	fread(vec, sizeof(scalar), ndofs, f);
	fclose(f);
}

void Solution::enable_transform(bool enable) {
	_F_
	if (transform != enable) free_tables();
	transform = enable;
}

//// ExactSolution /////////////////////////////////////////////////////////////////////////////////

ExactSolution::ExactSolution(Mesh *mesh, exact_fn_t fn)
	: MeshFunction(mesh)
{
	_F_
	this->fn = fn;
	this->num_components = 1;
	memset(tables, 0, sizeof(tables));
}

ExactSolution::ExactSolution(Mesh *mesh, exact_vec_fn_t fn)
	: MeshFunction(mesh)
{
	_F_
	this->fn_vec = fn;
	this->num_components = 3;
	memset(tables, 0, sizeof(tables));
}

ExactSolution::~ExactSolution() {
	_F_
	free();
}

void ExactSolution::free() {
	_F_
	for (int i = 0; i < QUAD_COUNT; i++)
		free_sub_tables(&(tables[i]));
}

void ExactSolution::set_active_element(Element *e) {
	_F_
	MeshFunction::set_active_element(e);

	if (tables[cur_quad] != NULL) free_sub_tables(&(tables[cur_quad]));
	sub_tables = &(tables[cur_quad]);
	update_nodes_ptr();

	switch (e->get_mode()) {
		case MODE_TETRAHEDRON: order = order3_t(MAX_QUAD_ORDER_TETRA); break;
		case MODE_HEXAHEDRON: order = order3_t(MAX_QUAD_ORDER, MAX_QUAD_ORDER, MAX_QUAD_ORDER); break;
		default: break;
	}
}


void ExactSolution::precalculate(qorder_t qord, int mask) {
	_F_
	Quad3D *quad = quads[cur_quad];
	assert(quad != NULL);
	QuadPt3D *pt = NULL;
	int np = 0;
	switch (qord.type) {
		case QOT_ELEMENT:
			np = quad->get_num_points(order3_t::from_int(qord.order));
			pt = quad->get_points(order3_t::from_int(qord.order));
			break;

		case QOT_FACE:
			np = quad->get_face_num_points(qord.face, order2_t::from_int(qord.order));
			pt = quad->get_face_points(qord.face, order2_t::from_int(qord.order));
			break;

		case QOT_EDGE:
			np = quad->get_edge_num_points(qord.order);
			pt = quad->get_edge_points(qord.edge, qord.order);
			break;

		case QOT_VERTEX:
			np = quad->get_vertex_num_points();
			pt = quad->get_vertex_points();
			break;

		default: assert(false);
	}

	assert(!(mask & ~FN_DEFAULT));
	mask = FN_DEFAULT;
	Node *node = new_node(mask, np);

    update_refmap();
    double *x, *y, *z;
    switch (qord.type) {
    	case QOT_ELEMENT:
			x = refmap->get_phys_x(order3_t::from_int(qord.order));
			y = refmap->get_phys_y(order3_t::from_int(qord.order));
			z = refmap->get_phys_z(order3_t::from_int(qord.order));
    		break;

    	case QOT_FACE:
    		x = refmap->get_face_phys_x(qord.face, order2_t::from_int(qord.order));
    		y = refmap->get_face_phys_y(qord.face, order2_t::from_int(qord.order));
    		z = refmap->get_face_phys_z(qord.face, order2_t::from_int(qord.order));
    		break;

    	case QOT_EDGE:
			x = refmap->get_edge_phys_x(qord.edge, qord.order);
			y = refmap->get_edge_phys_y(qord.edge, qord.order);
			z = refmap->get_edge_phys_z(qord.edge, qord.order);
    		break;

		case QOT_VERTEX:
			x = refmap->get_vertex_phys_x();
			y = refmap->get_vertex_phys_y();
			z = refmap->get_vertex_phys_z();
    		break;
    }

	// evaluate the exact solution
    if (num_components == 1) {
		for (int i = 0; i < np; i++) {
			scalar val, dx = 0.0, dy = 0.0, dz = 0.0;
			val = fn(x[i], y[i], z[i], dx, dy, dz);
			node->values[0][FN][i] = val;
			node->values[0][DX][i] = dx;
			node->values[0][DY][i] = dy;
			node->values[0][DZ][i] = dz;
		}
    }
    else if (num_components == 3) {
		for (int i = 0; i < np; i++) {
			scalar3 dx = { 0.0, 0.0, 0.0 };
			scalar3 dy = { 0.0, 0.0, 0.0 };
			scalar3 dz = { 0.0, 0.0, 0.0 };
			scalar3 &val  = fn_vec(x[i], y[i], z[i], dx, dy, dz);
			for (int j = 0; j < num_components; j++) {
				node->values[j][FN][i] = val[j];
				node->values[j][DX][i] = dx[j];
				node->values[j][DY][i] = dy[j];
				node->values[j][DZ][i] = dz[j];
			}
		}
    }
    else {
    	EXIT(ERR_FAILURE, "Invalid number of components.");
    }

	// remove the old node and attach the new one
	replace_cur_node(node);
}

//// ConstantSolution //////////////////////////////////////////////////////////////////////////////

void ConstantSolution::precalculate(qorder_t order, int mask) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
}

