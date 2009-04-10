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
#include "mesh.h"
#include "refmap.h"
#include "refdomain.h"
#include <common/error.h>
#include <common/trace.h>

#include "shapeset/common.h"
#include "shapeset/refmapss.h"

//

#ifdef WITH_TETRA
	static RefMapShapesetTetra		ref_map_shapeset_tetra;
	static PrecalcShapeset 			ref_map_pss_tetra(&ref_map_shapeset_tetra);
	#define REFMAP_SHAPESET_TETRA	&ref_map_shapeset_tetra
	#define REFMAP_PSS_TETRA		&ref_map_pss_tetra
#else
	#define REFMAP_SHAPESET_TETRA	NULL
	#define REFMAP_PSS_TETRA		NULL
#endif

#ifdef WITH_HEX
	static RefMapShapesetHex 		ref_map_shapeset_hex;
	static PrecalcShapeset			ref_map_pss_hex(&ref_map_shapeset_hex);
	#define REFMAP_SHAPESET_HEX		&ref_map_shapeset_hex
	#define REFMAP_PSS_HEX			&ref_map_pss_hex
#else
	#define REFMAP_SHAPESET_HEX		NULL
	#define REFMAP_PSS_HEX			NULL
#endif

// TODO: prisms

static PrecalcShapeset *ref_map_pss[] = { REFMAP_PSS_TETRA, REFMAP_PSS_HEX, NULL };


// RefMap /////////////////////////////////////////////////////////////////////////////////////////

RefMap::RefMap() {
	this->mesh = NULL;
	this->quad = NULL;
	this->pss = NULL;

	nodes = NULL;
	cur_node = NULL;
	overflow = NULL;
}

RefMap::RefMap(Mesh *mesh) {
	this->mesh = mesh;
	this->quad = NULL;
	this->pss = NULL;

	nodes = NULL;
	cur_node = NULL;
	overflow = NULL;
}

RefMap::~RefMap() {
	free();
}

void RefMap::set_quad(Quad3D *quad) {
	this->quad = quad;
	assert(this->pss != NULL);		// if this asserts, check that you called set_active_element() before calling set_quad()
	this->pss->set_quad(quad);
}

void RefMap::set_active_element(Element *e) {
	assert(e != NULL);

	if (e != element) free();

	EMode3D mode = e->get_mode();

	pss = ref_map_pss[mode];
	pss->set_active_element(e);

	quad = get_quadrature(mode);
	pss->set_quad(quad);

	if (e == element) return;
	element = e;

	reset_transform();
	update_cur_node();

	is_const = mode == MODE_TETRAHEDRON;

	int nvertices = element->get_num_of_vertices();

	// prepare the shapes and coefficients of the reference map
	Shapeset *shapeset = this->pss->get_shapeset();
	int i, k = 0;
	for (i = 0; i < nvertices; i++)
		indices[k++] = shapeset->get_vertex_index(i);

	// straight element
	for (int iv = 0; iv < nvertices; iv++)
		vertex[iv] = *mesh->vertices[e->get_vertex(iv)];
	coefs = vertex;
	num_coefs = nvertices;

	// calculate the order of the reference map
	switch (mode) {
		case MODE_TETRAHEDRON: ref_order = order3_t(0); break;
		case MODE_HEXAHEDRON:  ref_order = order3_t(1, 1, 1); break;
		case MODE_PRISM: EXIT(ERR_NOT_IMPLEMENTED); break;
	}

	// calculate the order of the inverse reference map
	switch (mode) {
		case MODE_TETRAHEDRON: inv_ref_order = order3_t(0); break;
		case MODE_HEXAHEDRON:  inv_ref_order = order3_t(1, 1, 1); break;
		case MODE_PRISM: EXIT(ERR_NOT_IMPLEMENTED); break;
	}

	// constant inverse reference map
	if (this->is_const) calc_const_inv_ref_map();
	else const_jacobian = 0.0;
}


void RefMap::push_transform(int son) {
	Transformable::push_transform(son);
	update_cur_node();
	const_jacobian *= 0.125;
}

void RefMap::pop_transform() {
	Transformable::pop_transform();
	update_cur_node();
	const_jacobian *= 8;
}

void RefMap::force_transform(uint64 sub_idx, Trf *ctm) {
	this->sub_idx = sub_idx;
	stack[top] = *ctm;
	ctm = stack + top;
	update_cur_node();
	if (is_const)
		calc_const_inv_ref_map();
}

// helpers

void RefMap::calc_inv_ref_map(order3_t order) {
	qorder_t qord = ELEM_QORDER(order);
	int np = quad->get_num_points(order3_t::from_int(qord.order));

	double3x3 *m = new double3x3[np];
	MEM_CHECK(m);
	memset(m, 0, np * sizeof(double3x3));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		double *dx, *dy, *dz;
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(qord);
		pss->get_dx_dy_dz_values(dx, dy, dz);
		for (int j = 0; j < np; j++) {
			m[j][0][0] += coefs[i].x * dx[j];
			m[j][0][1] += coefs[i].x * dy[j];
			m[j][0][2] += coefs[i].x * dz[j];
			m[j][1][0] += coefs[i].y * dx[j];
			m[j][1][1] += coefs[i].y * dy[j];
			m[j][1][2] += coefs[i].y * dz[j];
			m[j][2][0] += coefs[i].z * dx[j];
			m[j][2][1] += coefs[i].z * dy[j];
			m[j][2][2] += coefs[i].z * dz[j];
		}
	}

	double trj = get_transform_jacobian();
	double3x3 *irm = new double3x3[np];
	MEM_CHECK(irm);
	double *jac = new double[np];
	MEM_CHECK(jac);
	for (int i = 0; i < np; i++) {
		jac[i] =
			m[i][0][0] * m[i][1][1] * m[i][2][2] + m[i][0][1] * m[i][1][2] * m[i][2][0] + m[i][0][2] * m[i][1][0] * m[i][2][1] -
			m[i][2][0] * m[i][1][1] * m[i][0][2] - m[i][2][1] * m[i][1][2] * m[i][0][0] - m[i][2][2] * m[i][1][0] * m[i][0][1];

		double ij = 1.0 / jac[i];
		irm[i][0][0] = (m[i][1][1] * m[i][2][2] - m[i][1][2] * m[i][2][1]) * ij;
		irm[i][1][0] = (m[i][0][2] * m[i][2][1] - m[i][0][1] * m[i][2][2]) * ij;
		irm[i][2][0] = (m[i][0][1] * m[i][1][2] - m[i][0][2] * m[i][1][1]) * ij;
		irm[i][0][1] = (m[i][1][2] * m[i][2][0] - m[i][1][0] * m[i][2][2]) * ij;
		irm[i][1][1] = (m[i][0][0] * m[i][2][2] - m[i][0][2] * m[i][2][0]) * ij;
		irm[i][2][1] = (m[i][0][2] * m[i][1][0] - m[i][0][0] * m[i][1][2]) * ij;
		irm[i][0][2] = (m[i][1][0] * m[i][2][1] - m[i][1][1] * m[i][2][0]) * ij;
		irm[i][1][2] = (m[i][0][1] * m[i][2][0] - m[i][0][0] * m[i][2][1]) * ij;
		irm[i][2][2] = (m[i][0][0] * m[i][1][1] - m[i][0][1] * m[i][1][0]) * ij;

	    jac[i] *= trj;
	}

	cur_node->jacobian[order.get_idx()] = jac;
	cur_node->inv_ref_map[order.get_idx()] = irm;
	cur_node->ref_map[order.get_idx()] = m;
}


void RefMap::calc_const_inv_ref_map() {
	// for linear tetrahedra only
	// TODO: does not take in account the transformation (we do not have it for tetras, so this will probably work)

	double3x3 m = {
		{ (vertex[1].x - vertex[0].x) / 2, (vertex[2].x - vertex[0].x) / 2, (vertex[3].x - vertex[0].x) / 2 },
		{ (vertex[1].y - vertex[0].y) / 2, (vertex[2].y - vertex[0].y) / 2, (vertex[3].y - vertex[0].y) / 2 },
		{ (vertex[1].z - vertex[0].z) / 2, (vertex[2].z - vertex[0].z) / 2, (vertex[3].z - vertex[0].z) / 2 }
	};

	const_jacobian =
		m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1] -
		m[2][0] * m[1][1] * m[0][2] - m[2][1] * m[1][2] * m[0][0] - m[2][2] * m[1][0] * m[0][1];

	double ij = 1.0 / const_jacobian;

	const_inv_ref_map[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * ij;
	const_inv_ref_map[1][0] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * ij;
	const_inv_ref_map[2][0] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * ij;
	const_inv_ref_map[0][1] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * ij;
	const_inv_ref_map[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * ij;
	const_inv_ref_map[2][1] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * ij;
	const_inv_ref_map[0][2] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * ij;
	const_inv_ref_map[1][2] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * ij;
	const_inv_ref_map[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * ij;

//	const_jacobian *= get_transform_jacobian();
}


void RefMap::calc_phys_x(order3_t order) {
	// transform all x coordinates of the integration points
	int np = quad->get_num_points(order);
	double *x = cur_node->phys_x[order.get_idx()] = new double[np];
	MEM_CHECK(x);
	memset(x, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(ELEM_QORDER(order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			x[j] += coefs[i].x * fn[j];
	}
}


void RefMap::calc_phys_y(order3_t order) {
	// transform all y coordinates of the integration points
	int np = quad->get_num_points(order);
	double *y = cur_node->phys_y[order.get_idx()] = new double[np];
	MEM_CHECK(y);
	memset(y, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(ELEM_QORDER(order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			y[j] += coefs[i].y * fn[j];
	}
}


void RefMap::calc_phys_z(order3_t order) {
	// transform all z coordinates of the integration points
	int np = quad->get_num_points(order);
	double *z = cur_node->phys_z[order.get_idx()] = new double[np];
	MEM_CHECK(z);
	memset(z, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(ELEM_QORDER(order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			z[j] += coefs[i].z * fn[j];
	}
}

// edge related //

void RefMap::calc_edge_phys_x(int edge, order1_t order) {
	// transform all x coordinates of the integration points
	int np = quad->get_edge_num_points(order);
	double *x = cur_node->edge_phys_x[edge][order] = new double[np];
	MEM_CHECK(x);
	memset(x, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(EDGE_QORDER(edge, order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			x[j] += coefs[i].x * fn[j];
	}
}

void RefMap::calc_edge_phys_y(int edge, order1_t order) {
	// transform all y coordinates of the integration points
	int np = quad->get_edge_num_points(order);
	double *y = cur_node->edge_phys_y[edge][order] = new double[np];
	MEM_CHECK(y);
	memset(y, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(EDGE_QORDER(edge, order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			y[j] += coefs[i].y * fn[j];
	}
}

void RefMap::calc_edge_phys_z(int edge, order1_t order) {
	// transform all z coordinates of the integration points
	int np = quad->get_edge_num_points(order);
	double *z = cur_node->edge_phys_z[edge][order] = new double[np];
	MEM_CHECK(z);
	memset(z, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(EDGE_QORDER(edge, order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			z[j] += coefs[i].z * fn[j];
	}
}

// face related //

// TODO: rewrite in similar way for non-constant functions
void RefMap::calc_face_const_jacobian(int face) {
	assert(cur_node->face_mode[face] == MODE_TRIANGLE);

	// physical triangle
	const int *face_vtx = element->get_face_vertices(face);
	Vertex vtx[Tri::NUM_VERTICES];
	for (int i = 0; i < element->get_num_of_vertices(); i++)
		vtx[i] = vertex[face_vtx[i]];

	double3x3 m = {
		{ (vtx[1].x - vtx[0].x), (vtx[2].x - vtx[0].x), 1.0 },
		{ (vtx[1].y - vtx[0].y), (vtx[2].y - vtx[0].y), 1.0 },
		{ (vtx[1].z - vtx[0].z), (vtx[2].z - vtx[0].z), 1.0 }
	};
	double phys_triangle_surface = 0.5 *
		sqrt(
			sqr(m[1][0] * m[2][1] - m[1][1] * m[2][0]) +
			sqr(m[0][0] * m[2][1] - m[0][1] * m[2][0]) +
			sqr(m[0][0] * m[1][1] - m[0][1] * m[1][0]));

	// reference triangle from a tetrahedron
	const int *ref_idx = RefTetra::get_face_vertices(face);
	const Point3D *ref_vtx = RefTetra::get_vertices();

	double3x3 n = {
		{ ref_vtx[ref_idx[1]].x - ref_vtx[ref_idx[0]].x, ref_vtx[ref_idx[2]].x - ref_vtx[ref_idx[0]].x, 1.0 },
		{ ref_vtx[ref_idx[1]].y - ref_vtx[ref_idx[0]].y, ref_vtx[ref_idx[2]].y - ref_vtx[ref_idx[0]].y, 1.0 },
		{ ref_vtx[ref_idx[1]].z - ref_vtx[ref_idx[0]].z, ref_vtx[ref_idx[2]].z - ref_vtx[ref_idx[0]].z, 1.0 }
	};

	double ref_triangle_surface = 0.5 *
		sqrt(
			sqr(n[1][0] * n[2][1] - n[1][1] * n[2][0]) +
			sqr(n[0][0] * n[2][1] - n[0][1] * n[2][0]) +
			sqr(n[0][0] * n[1][1] - n[0][1] * n[1][0]));


	cur_node->face_const_jacobian[face] = phys_triangle_surface / ref_triangle_surface;
}


void RefMap::calc_face_jacobian(int face, order2_t order) {
	assert(mesh != NULL);

	int np = quad->get_face_num_points(face, order);

	double3x3 *m = new double3x3[np];
	MEM_CHECK(m);
	memset(m, 0, np * sizeof(double3x3));
	const int *face_vertices = RefHex::get_face_vertices(face);

	qorder_t surf_order = FACE_QORDER(face, order);

	for (int i = 0; i < RefHex::get_face_num_of_vertices(face); i++) {
		double *dx, *dy, *dz;
		pss->set_active_shape(indices[face_vertices[i]]);
		pss->set_quad_order(surf_order);
		pss->get_dx_dy_dz_values(dx, dy, dz);
		for (int j = 0; j < np; j++) {
			m[j][0][0] += vertex[face_vertices[i]].x * dx[j];
			m[j][0][1] += vertex[face_vertices[i]].x * dy[j];
			m[j][0][2] += vertex[face_vertices[i]].x * dz[j];
			m[j][1][0] += vertex[face_vertices[i]].y * dx[j];
			m[j][1][1] += vertex[face_vertices[i]].y * dy[j];
			m[j][1][2] += vertex[face_vertices[i]].y * dz[j];
			m[j][2][0] += vertex[face_vertices[i]].z * dx[j];
			m[j][2][1] += vertex[face_vertices[i]].z * dy[j];
			m[j][2][2] += vertex[face_vertices[i]].z * dz[j];
		}
	}

	// what vectors to take to make vector product
	int dir1, dir2;
	switch (face) {
		case 0:
		case 1:
			dir1 = 1; dir2 = 2;
			break;

		case 2:
		case 3:
			dir1 = 0; dir2 = 2;
			break;

		case 4:
		case 5:
			dir1 = 0; dir2 = 1;
			break;
	}

	// in fact, it is not jacobian, but vector product
	double *jac = new double[np];
	MEM_CHECK(jac);
	for (int i = 0; i < np; i++) {
		Point3D vec1 = { m[i][0][dir1], m[i][1][dir1], m[i][2][dir1] };
		Point3D vec2 = { m[i][0][dir2], m[i][1][dir2], m[i][2][dir2] };
		Point3D vec_product = cross_product(vec1, vec2);
		jac[i] = norm(vec_product);
	}

	cur_node->face_jacobian[face][order.get_idx()] = jac;
	delete [] m;
}

// this is in fact identical to calc_inv_ref_map
// the only difference is, that everything is calculated in integration points on given face
void RefMap::calc_face_inv_ref_map(int face, order2_t order) {
	int np = quad->get_face_num_points(face, order);

	double3x3 *m = new double3x3[np];
	MEM_CHECK(m);
	memset(m, 0, np * sizeof(double3x3));

	for (int i = 0; i < num_coefs; i++) {
		double *dx, *dy, *dz;
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(FACE_QORDER(face, order));
		pss->get_dx_dy_dz_values(dx, dy, dz);
		for (int j = 0; j < np; j++) {
			m[j][0][0] += vertex[i].x * dx[j];
			m[j][0][1] += vertex[i].x * dy[j];
			m[j][0][2] += vertex[i].x * dz[j];
			m[j][1][0] += vertex[i].y * dx[j];
			m[j][1][1] += vertex[i].y * dy[j];
			m[j][1][2] += vertex[i].y * dz[j];
			m[j][2][0] += vertex[i].z * dx[j];
			m[j][2][1] += vertex[i].z * dy[j];
			m[j][2][2] += vertex[i].z * dz[j];
		}
	}

	double3x3 *irm = new double3x3[np];
	MEM_CHECK(irm);
	double *jac = new double[np];
	MEM_CHECK(jac);
	for (int i = 0; i < np; i++) {
		jac[i] =
				m[i][0][0] * m[i][1][1] * m[i][2][2] + m[i][0][1] * m[i][1][2] * m[i][2][0] + m[i][0][2] * m[i][1][0] * m[i][2][1] -
				m[i][2][0] * m[i][1][1] * m[i][0][2] - m[i][2][1] * m[i][1][2] * m[i][0][0] - m[i][2][2] * m[i][1][0] * m[i][0][1];

		double ij = 1.0 / jac[i];
		irm[i][0][0] = (m[i][1][1] * m[i][2][2] - m[i][1][2] * m[i][2][1]) * ij;
		irm[i][1][0] = (m[i][0][2] * m[i][2][1] - m[i][0][1] * m[i][2][2]) * ij;
		irm[i][2][0] = (m[i][0][1] * m[i][1][2] - m[i][0][2] * m[i][1][1]) * ij;
		irm[i][0][1] = (m[i][1][2] * m[i][2][0] - m[i][1][0] * m[i][2][2]) * ij;
		irm[i][1][1] = (m[i][0][0] * m[i][2][2] - m[i][0][2] * m[i][2][0]) * ij;
		irm[i][2][1] = (m[i][0][2] * m[i][1][0] - m[i][0][0] * m[i][1][2]) * ij;
		irm[i][0][2] = (m[i][1][0] * m[i][2][1] - m[i][1][1] * m[i][2][0]) * ij;
		irm[i][1][2] = (m[i][0][1] * m[i][2][0] - m[i][0][0] * m[i][2][1]) * ij;
		irm[i][2][2] = (m[i][0][0] * m[i][1][1] - m[i][0][1] * m[i][1][0]) * ij;
	}

	cur_node->face_inv_ref_map[face][order.get_idx()] = irm;
	cur_node->face_ref_map[face][order.get_idx()] = m;

	// jac is jacobian of 3d map, evaluated in face integration points (2d)
	// it is NOT face jacobian, since face jacobian is jacobian of 2d map ( R^2 -> R^2 )
	// therefore face jacobian has to be calculated in other way and jac can be deleted
	delete [] jac;
}

void RefMap::calc_face_phys_x(int face, order2_t order) {
	// transform all x coordinates of the integration points
	int np = quad->get_face_num_points(face, order);
	double *x = cur_node->face_phys_x[face][order.get_idx()] = new double[np];
	MEM_CHECK(x);
	memset(x, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(FACE_QORDER(face, order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			x[j] += coefs[i].x * fn[j];
	}
}

void RefMap::calc_face_phys_y(int face, order2_t order) {
	// transform all y coordinates of the integration points
	int np = quad->get_face_num_points(face, order);
	double *y = cur_node->face_phys_y[face][order.get_idx()] = new double[np];
	MEM_CHECK(y);
	memset(y, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(FACE_QORDER(face, order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			y[j] += coefs[i].y * fn[j];
	}
}

void RefMap::calc_face_phys_z(int face, order2_t order) {
	// transform all z coordinates of the integration points
	int np = quad->get_face_num_points(face, order);
	double *z = cur_node->face_phys_z[face][order.get_idx()] = new double[np];
	MEM_CHECK(z);
	memset(z, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(FACE_QORDER(face, order));
		double *fn = pss->get_fn_values();
		for (int j = 0; j < np; j++)
			z[j] += coefs[i].z * fn[j];
	}
}

void RefMap::calc_face_normal(int face, order2_t order) {
	assert(mesh != NULL);

	int np = quad->get_face_num_points(face, order);
	double3x3 *m = get_face_ref_map(face, order);

	Point3D *normal = new Point3D[np];
	MEM_CHECK(normal);
	int t_dir_1, t_dir_2; //directions of tangents ot the reference face such that t_dir_1 x t_dir_2 = outer normal
	switch (face) {
		case 0: t_dir_1 = 2; t_dir_2 = 1; break;
		case 1: t_dir_1 = 1; t_dir_2 = 2; break;
		case 2: t_dir_1 = 0; t_dir_2 = 2; break;
		case 3: t_dir_1 = 2; t_dir_2 = 0; break;
		case 4: t_dir_1 = 1; t_dir_2 = 0; break;
		case 5: t_dir_1 = 0; t_dir_2 = 1; break;
	}
	for (int i = 0; i < np; i++) {
		Point3D tangent1 = { m[i][0][t_dir_1], m[i][1][t_dir_1], m[i][2][t_dir_1] };
		Point3D tangent2 = { m[i][0][t_dir_2], m[i][1][t_dir_2], m[i][2][t_dir_2] };
		normal[i] = normalize(cross_product(tangent1, tangent2));
	}

	cur_node->face_normal[face][order.get_idx()] = normal;
}

void RefMap::calc_vertex_phys() {
	int nvtx = element->get_num_of_vertices();
	double *x = cur_node->vertex_phys_x = new double [nvtx]; MEM_CHECK(x);
	double *y = cur_node->vertex_phys_y = new double [nvtx]; MEM_CHECK(y);
	double *z = cur_node->vertex_phys_z = new double [nvtx]; MEM_CHECK(z);

	memset(x, 0, nvtx * sizeof(double));
	memset(y, 0, nvtx * sizeof(double));
	memset(z, 0, nvtx * sizeof(double));

	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(VTX_QORDER());
		double *fn = pss->get_fn_values();
		for (int j = 0; j < nvtx; j++) {
			x[j] += coefs[i].x * fn[j];
			y[j] += coefs[i].y * fn[j];
			z[j] += coefs[i].z * fn[j];
		}
	}
}

//

void RefMap::init_node(Node **pp) {
	Node *node = *pp = new Node;
	MEM_CHECK(node);

	assert(element != NULL);
	node->nedges = element->get_num_of_edges();
	node->nfaces = element->get_num_of_faces();

	// edges
	node->edge_inv_ref_map = new ArrayPtr<double3x3>[node->nedges]; MEM_CHECK(node->edge_inv_ref_map);
	node->edge_phys_x = new Array<double *>[node->nedges]; MEM_CHECK(node->edge_phys_x);
	node->edge_phys_y = new Array<double *>[node->nedges]; MEM_CHECK(node->edge_phys_y);
	node->edge_phys_z = new Array<double *>[node->nedges]; MEM_CHECK(node->edge_phys_z);

	// faces
	node->face_mode = new EMode2D[node->nfaces]; MEM_CHECK(node->face_mode);
	node->face_const_jacobian = new double[node->nfaces]; MEM_CHECK(node->face_const_jacobian);
	node->face_const_normal = new Point3D[node->nfaces]; MEM_CHECK(node->face_const_normal);
	for (int iface = 0; iface < node->nfaces; iface++) {
		node->face_mode[iface] = element->get_face_mode(iface);
		node->face_const_jacobian[iface] = 0.0;			// initialize jacobians to zero (meaning the jacobian needs to be computed)

		node->face_const_normal[iface].x = 0.;			// initialize normals to zero (meaning the normals needs to be computed)
		node->face_const_normal[iface].y = 0.;
		node->face_const_normal[iface].z = 0.;
	}
	node->face_jacobian    = new ArrayPtr<double>[node->nfaces];   MEM_CHECK(node->face_jacobian);
	node->face_normal      = new ArrayPtr<Point3D>[node->nfaces];   MEM_CHECK(node->face_normal);
	node->face_ref_map     = new ArrayPtr<double3x3>[node->nfaces]; MEM_CHECK(node->face_ref_map);
	node->face_inv_ref_map = new ArrayPtr<double3x3>[node->nfaces]; MEM_CHECK(node->face_inv_ref_map);

	node->face_phys_x = new Array<double *>[node->nfaces]; MEM_CHECK(node->face_phys_x);
	node->face_phys_y = new Array<double *>[node->nfaces]; MEM_CHECK(node->face_phys_y);
	node->face_phys_z = new Array<double *>[node->nfaces]; MEM_CHECK(node->face_phys_z);

	node->vertex_inv_ref_map = NULL;
	node->vertex_phys_x = NULL;
	node->vertex_phys_y = NULL;
	node->vertex_phys_z = NULL;
}


void RefMap::free_node(Node *node) {
	for (Word_t idx = node->jacobian.first(); idx != INVALID_IDX; idx = node->jacobian.next(idx))
		delete [] node->jacobian[idx];
	node->jacobian.remove_all();
	for (Word_t idx = node->ref_map.first(); idx != INVALID_IDX; idx = node->ref_map.next(idx))
		delete [] node->ref_map[idx];
	node->ref_map.remove_all();
	for (Word_t idx = node->inv_ref_map.first(); idx != INVALID_IDX; idx = node->inv_ref_map.next(idx))
		delete [] node->inv_ref_map[idx];
	node->inv_ref_map.remove_all();

	for (Word_t idx = node->phys_x.first(); idx != INVALID_IDX; idx = node->phys_x.next(idx))
		delete [] node->phys_x[idx];
	node->phys_x.remove_all();
	for (Word_t idx = node->phys_y.first(); idx != INVALID_IDX; idx = node->phys_y.next(idx))
		delete [] node->phys_y[idx];
	node->phys_y.remove_all();
	for (Word_t idx = node->phys_z.first(); idx != INVALID_IDX; idx = node->phys_z.next(idx))
		delete [] node->phys_z[idx];
	node->phys_z.remove_all();

	// edges
	for (int edge = 0; edge < node->nedges; edge++) {
		if (node->edge_inv_ref_map != NULL) {
			for (Word_t idx = node->edge_inv_ref_map[edge].first(); idx != INVALID_IDX; idx = node->edge_inv_ref_map[edge].next(idx))
				delete [] node->edge_inv_ref_map[edge][idx];
			node->edge_inv_ref_map[edge].remove_all();
		}
		if (node->edge_phys_x != NULL)
			for (Word_t idx = node->edge_phys_x[edge].first(); idx != INVALID_IDX; idx = node->edge_phys_x[edge].next(idx))
				delete [] node->edge_phys_x[edge][idx];
		if (node->edge_phys_y != NULL)
			for (Word_t idx = node->edge_phys_y[edge].first(); idx != INVALID_IDX; idx = node->edge_phys_y[edge].next(idx))
				delete [] node->edge_phys_y[edge][idx];
		if (node->edge_phys_z != NULL)
			for (Word_t idx = node->edge_phys_z[edge].first(); idx != INVALID_IDX; idx = node->edge_phys_z[edge].next(idx))
				delete [] node->edge_phys_z[edge][idx];
	}
	delete [] node->edge_phys_x;
	delete [] node->edge_phys_y;
	delete [] node->edge_phys_z;

	// faces
	for (int face = 0; face < node->nfaces; face++) {
		if (node->face_phys_x != NULL)
			for (Word_t idx = node->face_phys_x[face].first(); idx != INVALID_IDX; idx = node->face_phys_x[face].next(idx))
				delete [] node->face_phys_x[face][idx];
		if (node->face_phys_y != NULL)
			for (Word_t idx = node->face_phys_y[face].first(); idx != INVALID_IDX; idx = node->face_phys_y[face].next(idx))
				delete [] node->face_phys_y[face][idx];
		if (node->face_phys_z != NULL)
			for (Word_t idx = node->face_phys_z[face].first(); idx != INVALID_IDX; idx = node->face_phys_z[face].next(idx))
				delete [] node->face_phys_z[face][idx];
		if (node->face_jacobian != NULL)
			for (Word_t idx = node->face_jacobian[face].first(); idx != INVALID_IDX; idx = node->face_jacobian[face].next(idx))
				delete [] node->face_jacobian[face][idx];
		if (node->face_normal != NULL)
			for (Word_t idx = node->face_normal[face].first(); idx != INVALID_IDX; idx = node->face_normal[face].next(idx))
				delete [] node->face_normal[face][idx];
		if (node->face_ref_map != NULL) {
			for (Word_t idx = node->face_ref_map[face].first(); idx != INVALID_IDX; idx = node->face_ref_map[face].next(idx))
				delete [] node->face_ref_map[face][idx];
			node->face_ref_map[face].remove_all();
		}
		if (node->face_inv_ref_map != NULL) {
			for (Word_t idx = node->face_inv_ref_map[face].first(); idx != INVALID_IDX; idx = node->face_inv_ref_map[face].next(idx))
				delete [] node->face_inv_ref_map[face][idx];
			node->face_inv_ref_map[face].remove_all();
		}

	}
	delete [] node->face_phys_x;
	delete [] node->face_phys_y;
	delete [] node->face_phys_z;

	delete node->vertex_inv_ref_map;
	delete [] node->vertex_phys_x;
	delete [] node->vertex_phys_y;
	delete [] node->vertex_phys_z;

	delete [] node->face_mode;
	delete [] node->face_const_jacobian;
	delete [] node->face_jacobian;
	delete [] node->face_const_normal;
	delete [] node->face_normal;
	delete [] node->face_ref_map;
	delete [] node->face_inv_ref_map;

	delete node;
}

void RefMap::update_cur_node() {
	Node **pp = (sub_idx > max_idx) ? handle_overflow() : (Node **) JudyLIns(&nodes, sub_idx, NULL);
	if (*pp == NULL) init_node(pp);
	cur_node = *pp;
}

void RefMap::free() {
	unsigned long idx = 0;
	Node **pp = (Node **) JudyLFirst(nodes, &idx, NULL);
	while (pp != NULL) {
		free_node(*pp);
		pp = (Node **) JudyLNext(nodes, &idx, NULL);
	}
	JudyLFreeArray(&nodes, NULL);

	if (overflow != NULL) { free_node(overflow); overflow = NULL; }
}


RefMap::Node **RefMap::handle_overflow() {
	if (overflow != NULL) free_node(overflow);
	overflow = NULL;
	return &overflow;
}

void RefMap::calc_edge_inv_ref_map(int edge, order1_t order) {
	qorder_t qord = EDGE_QORDER(edge, order);
	int np = quad->get_edge_num_points(order1_t(qord.order));

	double3x3 *m = new double3x3[np];
	MEM_CHECK(m);
	memset(m, 0, np * sizeof(double3x3));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		double *dx, *dy, *dz;
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(qord);
		pss->get_dx_dy_dz_values(dx, dy, dz);
		for (int j = 0; j < np; j++) {
			m[j][0][0] += coefs[i].x * dx[j];
			m[j][0][1] += coefs[i].x * dy[j];
			m[j][0][2] += coefs[i].x * dz[j];
			m[j][1][0] += coefs[i].y * dx[j];
			m[j][1][1] += coefs[i].y * dy[j];
			m[j][1][2] += coefs[i].y * dz[j];
			m[j][2][0] += coefs[i].z * dx[j];
			m[j][2][1] += coefs[i].z * dy[j];
			m[j][2][2] += coefs[i].z * dz[j];
		}
	}

	double trj = get_transform_jacobian();
	double3x3 *irm = new double3x3[np];
	MEM_CHECK(irm);
	double *jac = new double[np];
	MEM_CHECK(jac);
	for (int i = 0; i < np; i++) {
		jac[i] =
			m[i][0][0] * m[i][1][1] * m[i][2][2] + m[i][0][1] * m[i][1][2] * m[i][2][0] + m[i][0][2] * m[i][1][0] * m[i][2][1] -
			m[i][2][0] * m[i][1][1] * m[i][0][2] - m[i][2][1] * m[i][1][2] * m[i][0][0] - m[i][2][2] * m[i][1][0] * m[i][0][1];

		double ij = 1.0 / jac[i];
		irm[i][0][0] = (m[i][1][1] * m[i][2][2] - m[i][1][2] * m[i][2][1]) * ij;
		irm[i][1][0] = (m[i][0][2] * m[i][2][1] - m[i][0][1] * m[i][2][2]) * ij;
		irm[i][2][0] = (m[i][0][1] * m[i][1][2] - m[i][0][2] * m[i][1][1]) * ij;
		irm[i][0][1] = (m[i][1][2] * m[i][2][0] - m[i][1][0] * m[i][2][2]) * ij;
		irm[i][1][1] = (m[i][0][0] * m[i][2][2] - m[i][0][2] * m[i][2][0]) * ij;
		irm[i][2][1] = (m[i][0][2] * m[i][1][0] - m[i][0][0] * m[i][1][2]) * ij;
		irm[i][0][2] = (m[i][1][0] * m[i][2][1] - m[i][1][1] * m[i][2][0]) * ij;
		irm[i][1][2] = (m[i][0][1] * m[i][2][0] - m[i][0][0] * m[i][2][1]) * ij;
		irm[i][2][2] = (m[i][0][0] * m[i][1][1] - m[i][0][1] * m[i][1][0]) * ij;

	    jac[i] *= trj;
	}

	cur_node->edge_inv_ref_map[edge][order] = irm;

	delete [] m;
	delete [] jac;
}

void RefMap::calc_vertex_inv_ref_map() {
	int np = quad->get_vertex_num_points();

	double3x3 *m = new double3x3[np];
	MEM_CHECK(m);
	memset(m, 0, np * sizeof(double3x3));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < num_coefs; i++) {
		double *dx, *dy, *dz;
		pss->set_active_shape(indices[i]);
		pss->set_quad_order(VTX_QORDER());
		pss->get_dx_dy_dz_values(dx, dy, dz);
		for (int j = 0; j < np; j++) {
			m[j][0][0] += coefs[i].x * dx[j];
			m[j][0][1] += coefs[i].x * dy[j];
			m[j][0][2] += coefs[i].x * dz[j];
			m[j][1][0] += coefs[i].y * dx[j];
			m[j][1][1] += coefs[i].y * dy[j];
			m[j][1][2] += coefs[i].y * dz[j];
			m[j][2][0] += coefs[i].z * dx[j];
			m[j][2][1] += coefs[i].z * dy[j];
			m[j][2][2] += coefs[i].z * dz[j];
		}
	}

	double trj = get_transform_jacobian();
	double3x3 *irm = new double3x3[np];
	MEM_CHECK(irm);
	double *jac = new double[np];
	MEM_CHECK(jac);
	for (int i = 0; i < np; i++) {
		jac[i] =
			m[i][0][0] * m[i][1][1] * m[i][2][2] + m[i][0][1] * m[i][1][2] * m[i][2][0] + m[i][0][2] * m[i][1][0] * m[i][2][1] -
			m[i][2][0] * m[i][1][1] * m[i][0][2] - m[i][2][1] * m[i][1][2] * m[i][0][0] - m[i][2][2] * m[i][1][0] * m[i][0][1];

		double ij = 1.0 / jac[i];
		irm[i][0][0] = (m[i][1][1] * m[i][2][2] - m[i][1][2] * m[i][2][1]) * ij;
		irm[i][1][0] = (m[i][0][2] * m[i][2][1] - m[i][0][1] * m[i][2][2]) * ij;
		irm[i][2][0] = (m[i][0][1] * m[i][1][2] - m[i][0][2] * m[i][1][1]) * ij;
		irm[i][0][1] = (m[i][1][2] * m[i][2][0] - m[i][1][0] * m[i][2][2]) * ij;
		irm[i][1][1] = (m[i][0][0] * m[i][2][2] - m[i][0][2] * m[i][2][0]) * ij;
		irm[i][2][1] = (m[i][0][2] * m[i][1][0] - m[i][0][0] * m[i][1][2]) * ij;
		irm[i][0][2] = (m[i][1][0] * m[i][2][1] - m[i][1][1] * m[i][2][0]) * ij;
		irm[i][1][2] = (m[i][0][1] * m[i][2][0] - m[i][0][0] * m[i][2][1]) * ij;
		irm[i][2][2] = (m[i][0][0] * m[i][1][1] - m[i][0][1] * m[i][1][0]) * ij;

	    jac[i] *= trj;
	}

	cur_node->vertex_inv_ref_map = irm;

	delete [] m;
	delete [] jac;
}
