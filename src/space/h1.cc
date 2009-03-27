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

#include "../config.h"
#include "h1.h"
#include "../matrix.h"
#include "../refmap.h"
#include <common/bitarray.h>
#include <common/trace.h>
#include <common/error.h>

//#define ADD_ASMLIST_THRESHOLD					1e-13
#define ADD_ASMLIST_THRESHOLD					0

H1Space::H1Space(Mesh *mesh, Shapeset *ss) :
	Space(mesh, ss)
{
}

H1Space::~H1Space() {
}

Space *H1Space::dup(Mesh *mesh) const {
	H1Space *space = new H1Space(mesh, shapeset);
	space->copy_callbacks(this);
	return space;
}

// ndofs ////

int H1Space::get_vertex_ndofs() {
	return 1;
}

int H1Space::get_edge_ndofs(order1_t order) {
	return order - 1;
}

int H1Space::get_face_ndofs(order2_t order) {
	switch (order.type) {
		case MODE_TRIANGLE: return (order.order - 1) * (order.order - 2) / 2;
		case MODE_QUAD: return (order.x - 1) * (order.y - 1);
		default: EXIT(ERR_UNKNOWN_MODE); return -1;
	}
}

int H1Space::get_element_ndofs(order3_t order) {
	switch (order.type) {
		case MODE_TETRAHEDRON: return (order.order - 1) * (order.order - 2) * (order.order - 3) / 6;
		case MODE_HEXAHEDRON: return (order.x - 1) * (order.y - 1) * (order.z - 1);
		default: EXIT(ERR_UNKNOWN_MODE); return -1;
	}
}

void H1Space::assign_dofs_internal() {
	BitArray init_vertices;
	BitArray init_edges;
	BitArray init_faces;

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *e = mesh->elements[idx];
		// vertex dofs
		for (int ivtx = 0; ivtx < e->get_num_of_vertices(); ivtx++) {
			Word_t vid = e->get_vertex(ivtx);
			VertexData *vd = vn_data[vid];
			assert(vd != NULL);
			if (!init_vertices.is_set(vid) && !vd->ced) {
				assign_vertex_dofs(vid);
				init_vertices.set(vid);
			}
		}
	}

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

void H1Space::get_element_assembly_list(Element *e, AsmList *al) {
	al->clear();
	for (int i = 0; i < e->get_num_of_vertices(); i++) get_vertex_assembly_list(e, i, al);
	for (int i = 0; i < e->get_num_of_edges(); i++) get_edge_assembly_list(e, i, al);
	for (int i = 0; i < e->get_num_of_faces(); i++) get_face_assembly_list(e, i, al);
	get_bubble_assembly_list(e, al);
}

void H1Space::get_boundary_assembly_list(Element *e, int face, AsmList *al) {
	al->clear();
	const int *face_vtcs = e->get_face_vertices(face);
	const int *face_edges = e->get_face_edges(face);
	for (int i = 0; i < e->get_face_num_of_vertices(face); i++) get_vertex_assembly_list(e, face_vtcs[i], al);
	for (int i = 0; i < e->get_face_num_of_edges(face); i++) get_edge_assembly_list(e, face_edges[i], al);
	get_face_assembly_list(e, face, al);
}


// boundary projections ////

void H1Space::calc_vertex_boundary_projection(Element *elem, int ivertex) {
	Word_t vtx = elem->get_vertex(ivertex);
	VertexData *vnode = vn_data[vtx];
	Vertex *v = mesh->vertices[vtx];
	if (vnode->bc_type == BC_ESSENTIAL)
		vnode->bc_proj = bc_value_callback_by_coord(vnode->marker, v->x, v->y, v->z);
}

void H1Space::calc_edge_boundary_projection(Element *elem, int iedge) {
	Word_t edge_id = mesh->get_edge_id(elem, iedge);
	EdgeData *enode = en_data[edge_id];
	if (enode->bc_type != BC_ESSENTIAL) return;			// process only Dirichlet BC
	if (enode->bc_proj != NULL) return;					// projection already calculated

	int num_fns;
	if (enode->ced) {
		assert(enode->edge_ncomponents > 0);
		Word_t edge_id = enode->edge_baselist[0].edge_id;
		num_fns = en_data[edge_id]->n;
	}
	else {
		num_fns = enode->n;
	}
	if (num_fns <= 0) return;

	double **proj_mat = new_matrix<double>(num_fns, num_fns);
	if (proj_mat == NULL) EXIT(ERR_OUT_OF_MEMORY);
	double *proj_rhs = new double[num_fns];
	if (proj_rhs == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(proj_rhs, 0, sizeof(double) * num_fns);

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());
	// local edge vertex numbers
	const int *local_edge_vtx = elem->get_edge_vertices(iedge);
	// edge vertices (global indices)
	Word_t edge_vtx[2] = { elem->get_vertex(local_edge_vtx[0]), elem->get_vertex(local_edge_vtx[1]) };

	int vtx_fn_idx[] = { shapeset->get_vertex_index(local_edge_vtx[0]), shapeset->get_vertex_index(local_edge_vtx[1]) };
	// function values at vertices
	double vtx_fn_coef[] = { vn_data[edge_vtx[0]]->bc_proj, vn_data[edge_vtx[1]]->bc_proj };
	int *edge_fn_idx = new int[num_fns];
	if (enode->ced && enode->edge_ncomponents > 0) {
		BaseEdgeComponent *ecomp = enode->edge_baselist + 0;
		EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

		int *indices = shapeset->get_edge_indices(iedge, ecomp->ori, cng_enode->order);
		for (int j = 0; j < cng_enode->n; j++) {
			int order = shapeset->get_order(indices[j]).get_edge_order(iedge);
			edge_fn_idx[j] = shapeset->get_constrained_edge_index(iedge, ecomp->ori, order, ecomp->part);
		}
	}
	else {
		int ori = elem->get_edge_orientation(iedge);								// edge orientation
		int *idx = shapeset->get_edge_indices(iedge, ori, enode->order);			// indices of edge functions
		for (int m = 0; m < enode->n; m++)
			edge_fn_idx[m] = idx[m];
	}

	for (int i = 0; i < num_fns; i++) {
		int iidx = edge_fn_idx[i];
		for (int j = i; j < num_fns; j++) {
			int jidx = edge_fn_idx[j];

			order3_t order = shapeset->get_order(iidx) + shapeset->get_order(jidx);
			int edge_order = order.get_edge_order(iedge);
			QuadPt3D *pt = quad->get_edge_points(iedge, edge_order);
			double value = 0.0;
			for (int k = 0; k < quad->get_edge_num_points(edge_order); k++)
				value += pt[k].w *
					shapeset->get_fn_value(iidx, pt[k].x, pt[k].y, pt[k].z, 0) *
					shapeset->get_fn_value(jidx, pt[k].x, pt[k].y, pt[k].z, 0);
			proj_mat[i][j] += value;
		}

		int order_rhs = quad->get_edge_max_order(iedge);
		double *edge_phys_x = ref_map.get_edge_phys_x(iedge, order_rhs);
		double *edge_phys_y = ref_map.get_edge_phys_y(iedge, order_rhs);
		double *edge_phys_z = ref_map.get_edge_phys_z(iedge, order_rhs);

		double value = 0.0;
		QuadPt3D *pt = quad->get_edge_points(iedge, order_rhs);
		for (int k = 0; k < quad->get_edge_num_points(order_rhs); k++) {
			scalar g =
				vtx_fn_coef[0] * shapeset->get_fn_value(vtx_fn_idx[0], pt[k].x, pt[k].y, pt[k].z, 0) +
				vtx_fn_coef[1] * shapeset->get_fn_value(vtx_fn_idx[1], pt[k].x, pt[k].y, pt[k].z, 0);
			value += pt[k].w *
				shapeset->get_fn_value(iidx, pt[k].x, pt[k].y, pt[k].z, 0) *
				(bc_value_callback_by_coord(enode->marker, edge_phys_x[k], edge_phys_y[k], edge_phys_z[k]) - g);
		}
		proj_rhs[i] += value;
	}

	double *chol_p = new double[num_fns];
	choldc(proj_mat, num_fns, chol_p);
	cholsl(proj_mat, num_fns, chol_p, proj_rhs, proj_rhs);
	delete [] chol_p;

	enode->bc_proj = proj_rhs;

	delete [] proj_mat;
	delete [] edge_fn_idx;
}

void H1Space::calc_face_boundary_projection(Element *elem, int iface) {
#ifdef COMPLEX
	assert(0); //not implemented
#else
	Word_t facet_idx = mesh->get_facet_id(elem, iface);
	FaceData *fnode = fn_data[facet_idx];

	if (fnode->bc_type != BC_ESSENTIAL) return;
	if (fnode->bc_proj != NULL) return;
	if (fnode->n <= 0) return;

	double **proj_mat = new_matrix<double>(fnode->n, fnode->n);
	if (proj_mat == NULL) EXIT(ERR_OUT_OF_MEMORY);
	double *proj_rhs = new double[fnode->n];
	if (proj_rhs == NULL) EXIT(ERR_OUT_OF_MEMORY);
	memset(proj_rhs, 0, sizeof(double) * fnode->n);

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());

	const int *local_face_vertex = elem->get_face_vertices(iface);
	const int *local_face_edge = elem->get_face_edges(iface);

	// get total number of vertex + edge functions
	int num_fns = elem->get_face_num_of_vertices(iface);
	for (int edge = 0; edge < elem->get_face_num_of_edges(iface); edge++) {
		Word_t edge_idx = mesh->get_edge_id(elem, local_face_edge[edge]);
		EdgeData *enode = en_data[edge_idx];
		if (enode->ced && enode->edge_ncomponents > 0 && enode->edge_baselist != NULL) {
			Word_t eid = enode->edge_baselist[0].edge_id;
			num_fns += en_data[eid]->n;
		}
		else
			num_fns += enode->n;
	}

	double *coef = new double[num_fns];
	if (coef == NULL) EXIT(ERR_OUT_OF_MEMORY);
	int *fn_idx = new int[num_fns];
	if (fn_idx == NULL) EXIT(ERR_OUT_OF_MEMORY);

	int m = 0;
	// vertex projection coefficients
	for (int vtx = 0; vtx < elem->get_face_num_of_vertices(iface); vtx++, m++) {
		VertexData *vnode = vn_data[elem->get_vertex(local_face_vertex[vtx])];
		coef[m] = vnode->bc_proj;
		fn_idx[m] = shapeset->get_vertex_index(local_face_vertex[vtx]);
	}
	// edge projection coefficients
	for (int edge = 0; edge < elem->get_face_num_of_edges(iface); edge++) {
		Word_t edge_idx = mesh->get_edge_id(elem, local_face_edge[edge]);
		EdgeData *enode = en_data[edge_idx];

		if (enode->ced && enode->edge_ncomponents > 0 && enode->edge_baselist != NULL) {
				BaseEdgeComponent *ecomp = enode->edge_baselist + 0;
				EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

				int *indices = shapeset->get_edge_indices(local_face_edge[edge], ecomp->ori, cng_enode->order);
				for (int j = 0; j < cng_enode->n; j++, m++) {
					int order = shapeset->get_order(indices[j]).get_edge_order(local_face_edge[edge]);
					fn_idx[m] = shapeset->get_constrained_edge_index(local_face_edge[edge], ecomp->ori, order, ecomp->part);
					coef[m] = cng_enode->bc_proj[j];
				}
		}
		else {
			if (enode->n > 0) {
				int edge_ori = elem->get_edge_orientation(local_face_edge[edge]);
				int *edge_fn_idx = shapeset->get_edge_indices(local_face_edge[edge], edge_ori, enode->order);

				for (int i = 0; i < enode->n; i++, m++) {
					coef[m] = enode->bc_proj[i];
					fn_idx[m] = edge_fn_idx[i];
				}
			}
		}
	}

	// do it //
	int face_ori = elem->get_face_orientation(iface);
	int *face_fn_idx = shapeset->get_face_indices(iface, face_ori, fnode->order);

	for (int i = 0; i < fnode->n; i++) {
		int iidx = face_fn_idx[i];
		for (int j = i; j < fnode->n; j++) {
			int jidx = face_fn_idx[j];

			order3_t order = shapeset->get_order(iidx) + shapeset->get_order(jidx);
			order2_t face_order = order.get_face_order(iface);

			QuadPt3D *pt = quad->get_face_points(iface, face_order);
			double value = 0.0;
			for (int k = 0; k < quad->get_face_num_points(iface, face_order); k++) {
				value += pt[k].w *
					shapeset->get_fn_value(iidx, pt[k].x, pt[k].y, pt[k].z, 0) *
					shapeset->get_fn_value(jidx, pt[k].x, pt[k].y, pt[k].z, 0);
			}
			proj_mat[i][j] += value;
		}

		order2_t order_rhs = quad->get_face_max_order(iface);
		double *face_phys_x = ref_map.get_face_phys_x(iface, order_rhs);
		double *face_phys_y = ref_map.get_face_phys_y(iface, order_rhs);
		double *face_phys_z = ref_map.get_face_phys_z(iface, order_rhs);

		double value = 0.0;
		QuadPt3D *pt = quad->get_face_points(iface, order_rhs);
		for (int k = 0; k < quad->get_face_num_points(iface, order_rhs); k++) {
			double g = 0.0; // lin. combination of vertex + edge functions
			for (int l = 0; l < num_fns; l++)
				g += coef[l] * shapeset->get_fn_value(fn_idx[l], pt[k].x, pt[k].y, pt[k].z, 0);

			value += pt[k].w *
				shapeset->get_fn_value(iidx, pt[k].x, pt[k].y, pt[k].z, 0) *
				(bc_value_callback_by_coord(fnode->marker, face_phys_x[k], face_phys_y[k], face_phys_z[k]) - g);
		}

		proj_rhs[i] += value;
	}

	// solve the system using a precalculated Cholesky decomposed projection matrix
	double *chol_p = new double[fnode->n];
	choldc(proj_mat, fnode->n, chol_p);
	cholsl(proj_mat, fnode->n, chol_p, proj_rhs, proj_rhs);
	delete [] chol_p;

	fnode->bc_proj = proj_rhs;

	delete [] fn_idx;
	delete [] coef;
	delete [] proj_mat;
#endif
}

// CED stuff ////

void H1Space::update_constrained_nodes(Word_t fid) {
	Facet *facet = mesh->facets.get(fid);
	assert(facet != NULL);

	if (facet->type == Facet::OUTER)
		return;

	if (facet->ractive || facet->lactive) {
		facet->dump();			// WTF?
	}
	else {
		for (int i = 0; i < 4; i++) {
			if (facet->sons[i] != INVALID_IDX)
				update_constrained_nodes(facet->sons[i]);
		}
	}
}

void H1Space::update_constraints() {
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *e = mesh->elements[idx];
		for (int face = 0; face < e->get_num_of_faces(); face++) {
			Word_t fid = mesh->get_facet_id(e, face);
			update_constrained_nodes(fid);
		}
	}
}
