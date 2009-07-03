// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 Pavel Kus <pavel.kus@gmail.com>
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

#include "../h3dconfig.h"
#include "hcurl.h"
#include "../matrix.h"
#include "../refmap.h"
#include <common/bitarray.h>
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

//#define ADD_ASMLIST_THRESHOLD					1e-13
#define ADD_ASMLIST_THRESHOLD					0

HcurlSpace::HcurlSpace(Mesh *mesh, Shapeset *ss) :
		Space(mesh, ss)
{
	_F_
	this->type = Hcurl;
}

HcurlSpace::~HcurlSpace() {
	_F_
}

Space *HcurlSpace::dup(Mesh *mesh) const {
	_F_
	HcurlSpace *space = new HcurlSpace(mesh, shapeset);
	space->copy_callbacks(this);
	return space;
}

// ndofs ////

int HcurlSpace::get_vertex_ndofs() {
	return 0;
}

int HcurlSpace::get_edge_ndofs(order1_t order) {
	return order + 1;
}

int HcurlSpace::get_face_ndofs(order2_t order) {
	switch (order.type) {
		case MODE_QUAD: return (order.x + 1) * order.y + order.x * (order.y + 1);
		case MODE_TRIANGLE: EXIT(ERR_NOT_IMPLEMENTED); return -1;
		default: EXIT(ERR_UNKNOWN_MODE); return -1;
	}
}

int HcurlSpace::get_element_ndofs(order3_t order) {
	switch (order.type) {
		case MODE_HEXAHEDRON: return (order.x + 1) * order.y * order.z + order.x * (order.y + 1) * order.z + order.x * order.y * (order.z + 1);
		case MODE_TETRAHEDRON: EXIT(ERR_NOT_IMPLEMENTED); return -1;
		default: EXIT(ERR_UNKNOWN_MODE); return -1;
	}
}

//

void HcurlSpace::assign_dofs_internal() {
	_F_
	BitArray init_edges;
	BitArray init_faces;

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *e = mesh->elements[idx];
		// edge dofs
		for (int iedge = 0; iedge < e->get_num_of_edges(); iedge++) {
			Word_t eid = mesh->get_edge_id(e, iedge);
			EdgeData *ed = en_data[eid];
			assert(ed != NULL);
			if (!init_edges.is_set(eid) && !ed->ced) {
				assign_edge_dofs(eid);
				init_edges.set(eid);
			}
		}
	}

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *e = mesh->elements[idx];
		// face dofs
		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			FaceData *fd = fn_data[fid];
			assert(fd != NULL);
			if (!init_faces.is_set(fid) && !fd->ced) {
				assign_face_dofs(fid);
				init_faces.set(fid);
			}
		}
	}

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		assign_bubble_dofs(idx);
	}
}

// assembly lists ////

void HcurlSpace::get_element_assembly_list(Element *e, AsmList *al) {
	_F_
	al->clear();
	for (int i = 0; i < e->get_num_of_edges(); i++) get_edge_assembly_list(e, i, al);
	for (int i = 0; i < e->get_num_of_faces(); i++) get_face_assembly_list(e, i, al);
	get_bubble_assembly_list(e, al);
}

void HcurlSpace::get_boundary_assembly_list(Element *e, int face, AsmList *al) {
	_F_
	al->clear();
	const int *face_edges = e->get_face_edges(face);
	for (int i = 0; i < e->get_face_num_of_edges(face); i++) get_edge_assembly_list(e, face_edges[i], al);
	get_face_assembly_list(e, face, al);
}

// boundary projections ////
// we allowe only zero bc (see hcurl.h), so we only check, whether this is true
// and fill projection with zeros
void HcurlSpace::calc_vertex_boundary_projection(Element *elem, int ivertex) {
	_F_
	Word_t vtx = elem->get_vertex(ivertex);
	VertexData *vnode = vn_data[vtx];
	Vertex *v = mesh->vertices[vtx];
	if (vnode->bc_type == BC_ESSENTIAL) {
		// FIXME: use bc_vec_value_callback_by_coord
		if ((vnode->bc_proj = bc_value_callback_by_coord(vnode->marker, v->x, v->y, v->z)) != 0.)
			EXIT(ERR_NOT_IMPLEMENTED);  //projection of nonzero bc not implemented, see comment in .h
	}
}

// we allowe only zero bc (see hcurl.h), so we only check, whether this is true
// and fill projection with zeros
void HcurlSpace::calc_edge_boundary_projection(Element *elem, int iedge) {
	_F_
	Word_t edge = mesh->get_edge_id(elem, iedge);
	EdgeData *enode = en_data[edge];
	if (enode->bc_type != BC_ESSENTIAL) return;			// process only Dirichlet BC
	if (enode->bc_proj != NULL) return;					// projection already calculated

	scalar *proj_rhs = new scalar[enode->n];
	if (proj_rhs == NULL) EXIT(ERR_OUT_OF_MEMORY);
	for(int i = 0; i < enode->n; i++)
		proj_rhs[i] = 0.;

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());
	order1_t order_rhs = quad->get_edge_max_order(iedge);
	int np = quad->get_edge_num_points(iedge, order_rhs);
	QuadPt3D *pt = quad->get_edge_points(iedge, order_rhs);

	double *edge_phys_x = ref_map.get_phys_x(np, pt);
	double *edge_phys_y = ref_map.get_phys_y(np, pt);
	double *edge_phys_z = ref_map.get_phys_z(np, pt);

	for (int k = 0; k < np; k++) {
		// FIXME: use bc_vec_value_callback_by_coord
		if (bc_value_callback_by_coord(enode->marker, edge_phys_x[k], edge_phys_y[k], edge_phys_z[k]) != 0.)
			EXIT(ERR_NOT_IMPLEMENTED);  //projection of nonzero bc not implemented, see comment in .h
	}

	delete [] edge_phys_x;
	delete [] edge_phys_y;
	delete [] edge_phys_z;

	// save vector of zeros as a projection
	enode->bc_proj = proj_rhs;
}

// we allowe only zero bc (see hcurl.h), so we only check, whether this is true
// and fill projection with zeros
void HcurlSpace::calc_face_boundary_projection(Element *elem, int iface) {
	_F_
	Word_t facet_idx = mesh->get_facet_id(elem, iface);
	FaceData *fnode = fn_data[facet_idx];

	if (fnode->bc_type != BC_ESSENTIAL) return;
	if (fnode->bc_proj != NULL) return;

	scalar *proj_rhs = new scalar[fnode->n];
	if (proj_rhs == NULL) EXIT(ERR_OUT_OF_MEMORY);
	for(int i = 0; i < fnode->n; i++)
		proj_rhs[i] = 0.;

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());
	order2_t order_rhs = quad->get_face_max_order(iface);
	int np = quad->get_face_num_points(iface, order_rhs);
	QuadPt3D *pt = quad->get_face_points(iface, order_rhs);

	double *face_phys_x = ref_map.get_phys_x(np, pt);
	double *face_phys_y = ref_map.get_phys_y(np, pt);
	double *face_phys_z = ref_map.get_phys_z(np, pt);

	for (int k = 0; k < quad->get_face_num_points(iface, order_rhs); k++) {
		// FIXME: use bc_vec_value_callback_by_coord
		if (bc_value_callback_by_coord(fnode->marker, face_phys_x[k], face_phys_y[k], face_phys_z[k]) != 0.)
			EXIT(ERR_NOT_IMPLEMENTED);  //projection of nonzero bc not implemented, see comment in .h
	}

	delete [] face_phys_x;
	delete [] face_phys_y;
	delete [] face_phys_z;

	// save vector of zeros as a projection
	fnode->bc_proj = proj_rhs;
}


void HcurlSpace::update_constraints() {
}
