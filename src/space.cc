//
// space.cc
//

#include "config.h"
#include "common.h"
#include "space.h"
#include "matrix.h"
#include <common/error.h>
#include <common/timer.h>

//#define ADD_ASMLIST_THRESHOLD					1e-13
#define ADD_ASMLIST_THRESHOLD					0

#define PRINTF(...)
//#define PRINTF printf

#define INVALID_EDGE_ORDER						-1

#define CHECK_ELEMENT_ID(id) \
	if ((id) < 1 || (id) > mesh->elements.count())\
		EXIT(ERR_FAILURE, "Invalid element id (eid = %d).", id);\
	assert(mesh->elements.exists(id));


Space::Space(Mesh *mesh, Shapeset *shapeset) :
	mesh(mesh), shapeset(shapeset)
{
	set_bc_types(NULL);
	set_bc_values((scalar(*)(int, double, double, double, int)) NULL);
	init_data_tables();
}


Space::~Space() {
	free_data_tables();

	for (Word_t i = fi_data.first(); i != INVALID_IDX; i = fi_data.next(i))
		delete fi_data[i];
	fi_data.remove_all();

	for (Word_t i = ei_data.first(); i != INVALID_IDX; i = ei_data.next(i))
		delete ei_data[i];
	ei_data.remove_all();
}

void Space::init_data_tables() {
	assert(mesh != NULL);

	FOR_ALL_ELEMENTS(idx, mesh) {
		if (mesh->elements[idx]->active) {
			elm_data[idx] = new ElementData;
			MEM_CHECK(elm_data[idx]);
		}
	}
}

void Space::free_data_tables() {
	Word_t i;
	FOR_ALL_VERTEX_NODES(i)
		delete vn_data[i];
	vn_data.remove_all();

	FOR_ALL_EDGE_NODES(i)
		delete en_data[i];
	en_data.remove_all();

	FOR_ALL_FACE_NODES(i)
		delete fn_data[i];
	fn_data.remove_all();

	FOR_ALL_ELEMENT_NODES(i)
		delete elm_data[i];
	elm_data.remove_all();
}

// element orders ///////////////////////////////////////////////////////////////////////////////


void Space::set_element_order(Word_t eid, order3_t order) {
	CHECK_ELEMENT_ID(eid);

	// TODO: check for validity of order
	if (!elm_data.exists(eid)) {
		elm_data[eid] = new ElementData;
		MEM_CHECK(elm_data[eid]);
	}

	elm_data[eid]->order = order;
}

order3_t Space::get_element_order(Word_t eid) const {
	CHECK_ELEMENT_ID(eid);
	assert(elm_data.exists(eid));
	assert(elm_data[eid] != NULL);
	return elm_data[eid]->order;
}

void Space::set_uniform_order(order3_t order) {
	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		assert(elm_data.exists(eid));
		assert(elm_data[eid] != NULL);
		elm_data[eid]->order = order;
	}
}

void Space::set_order_recurrent(Word_t eid, order3_t order) {
	Element *e = mesh->elements[eid];
	if (e->active) {
		assert(elm_data.exists(e->id));
		assert(elm_data[e->id] != NULL);
		elm_data[e->id]->order = order;
	}
	else {
		for (int i = 0; i < e->get_num_of_sons(); i++) {
			Word_t son = e->get_son(i);
			if (son != INVALID_IDX)
				set_order_recurrent(son, order);
		}
	}
}

inline int LIMIT_ELEMENT_ORDER(int a) {
	if (a > MAX_ELEMENT_ORDER) return MAX_ELEMENT_ORDER;
	else return a;
}

void Space::copy_orders(const Space &space, int inc) {
	Mesh *cmesh = space.get_mesh();
	Word_t eid;
	FOR_ALL_ACTIVE_ELEMENTS(eid, cmesh) {
		order3_t oo = space.get_element_order(eid);
		assert(cmesh->elements[eid]->get_mode() == mesh->elements[eid]->get_mode());

		order3_t order;
		switch (cmesh->elements[eid]->get_mode()) {
			case MODE_TETRAHEDRON: order = oo + order3_t(inc); break;
			case MODE_HEXAHEDRON: order = oo + order3_t(inc, inc, inc); break;
			default: EXIT(ERR_NOT_IMPLEMENTED); break;
		}
		order.limit();

		set_order_recurrent(eid, order);
	}
}

void Space::enforce_minimum_rule() {
	// TODO: minimum CED rule

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *elem = mesh->elements[idx];
		ElementData *elem_node = elm_data[idx];
		order3_t elm_order = elem_node->order;

		switch (elem->get_mode()) {
			case MODE_TETRAHEDRON:
				// on faces
				for (int iface = 0; iface < elem->get_num_of_faces(); iface++) {
					Word_t fidx = mesh->get_facet_id(elem, iface);
					assert(fn_data.exists(fidx));
					FaceData *fnode = fn_data[fidx];
					order2_t forder = elem_node->order.get_face_order(iface);
					if (fnode->order.invalid() || forder.order < fnode->order.order)
						fnode->order = forder;
				}

				// on edges
				for (int iedge = 0; iedge < elem->get_num_of_edges(); iedge++) {
					Word_t eidx = mesh->get_edge_id(elem, iedge);
					assert(en_data.exists(eidx));

					EdgeData *enode = en_data[eidx];
					order1_t eorder = elem_node->order.get_edge_order(iedge);
					if (enode->order == INVALID_EDGE_ORDER || eorder < enode->order)
						enode->order = eorder;
				}
				break;

			case MODE_HEXAHEDRON:
				// on faces
				for (int iface = 0; iface < elem->get_num_of_faces(); iface++) {
					Word_t fidx = mesh->get_facet_id(elem, iface);
					FaceData *fnode = fn_data[fidx];

					if (!fnode->ced) {
						Facet *facet = mesh->facets.get(fidx);

						order2_t forder = elem_node->order.get_face_order(iface);
						if (elem->get_face_orientation(iface) >= 4) forder = order2_t(forder.y, forder.x);		// switch h- and v- order

						if (fnode->order.invalid())
							fnode->order = forder;
						else
							fnode->order = order2_t(std::min(fnode->order.x, forder.x), std::min(fnode->order.y, forder.y));
					}
				}

				// on edges
				for (int iedge = 0; iedge < elem->get_num_of_edges(); iedge++) {
					Word_t eidx = mesh->get_edge_id(elem, iedge);
					assert(eidx != INVALID_IDX);
					if (mesh->edges[eidx].is_active()) {
						EdgeData *enode = en_data[eidx];
						if (!enode->ced) {
							order1_t eorder = elem_node->order.get_edge_order(iedge);
							if (enode->order == INVALID_EDGE_ORDER || eorder < enode->order)
								enode->order = eorder;
						}
					}
				}
				break;

			default:
				EXIT(ERR_NOT_IMPLEMENTED);
		}
	}
}

//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void Space::assign_vertex_dofs(Word_t vid) {
	VertexData *node = vn_data[vid];
	int ndofs = get_vertex_ndofs();
	if (node->bc_type == BC_ESSENTIAL) {
		node->dof = DIRICHLET_DOF;
	}
	else {
		node->dof = next_dof;
		next_dof += ndofs * stride;
	}
	node->n = ndofs;
}

void Space::assign_edge_dofs(Word_t idx) {
	EdgeData *node = en_data[idx];
	int ndofs = get_edge_ndofs(node->order);
	if (node->bc_type == BC_ESSENTIAL) {
		node->dof = DIRICHLET_DOF;
	}
	else {
		node->dof = next_dof;
		next_dof += ndofs * stride;
	}
	node->n = ndofs;
}

void Space::assign_face_dofs(Word_t idx) {
	FaceData *node = fn_data[idx];
	int ndofs = get_face_ndofs(node->order);
	if (node->bc_type == BC_ESSENTIAL) {
		node->dof = DIRICHLET_DOF;
	}
	else {
		node->dof = next_dof;
		next_dof += ndofs * stride;
	}
	node->n = ndofs;
}

void Space::assign_bubble_dofs(Word_t idx) {
	ElementData *enode = elm_data[idx];
	int ndofs = get_element_ndofs(enode->order);
	enode->n = ndofs;
	enode->dof = next_dof;
	next_dof += ndofs * stride;
}

// assembly lists ////

void Space::get_vertex_assembly_list(Element *e, int ivertex, AsmList *al) {
	Word_t vtx = e->get_vertex(ivertex);
	VertexData *vnode = vn_data[vtx];
	int index = shapeset->get_vertex_index(ivertex);

	if (vnode->ced) {
//		printf("nc = %d\n", vnode->ncomponents);
		for (int i = 0; i < vnode->ncomponents; i++) {
			if (vnode->baselist[i].coef != 0) {
//				printf(" - dof = %d, coef = %lf\n", vnode->baselist[i].dof, vnode->baselist[i].coef);
				al->add(index, vnode->baselist[i].dof, vnode->baselist[i].coef);
			}
		}
//		printf("--\n");
	}
	else {
		double coef = vnode->dof >= 0 ? 1.0 : vnode->bc_proj;
		assert(vnode->dof >= DIRICHLET_DOF && vnode->dof < get_dof_count());
		al->add(index, vnode->dof, coef);
	}
}

void Space::get_edge_assembly_list(Element *elem, int iedge, AsmList *al) {
	Word_t edge_id = mesh->get_edge_id(elem, iedge);
	EdgeData *enode = en_data[edge_id];
	int ori = elem->get_edge_orientation(iedge);

	if (enode->ced) {
//		printf("  --\n");
		// edge constrained by an edge
		for (int i = 0; i < enode->edge_ncomponents; i++) {
			BaseEdgeComponent *ecomp = enode->edge_baselist + i;
			EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node
			assert(cng_enode->ced == false);

//			printf("edge_id = %d, eid = %d, iedge = %d, ori = %d, order = %d, part = %d\n",
//				edge_id, elem->id, iedge, ecomp->ori, cng_enode->order, ecomp->part.part);

			int *indices = shapeset->get_edge_indices(iedge, 0, cng_enode->order);		// iedge bude 0 (?)
			if (cng_enode->dof >= 0) {
				for (int j = 0, dof = cng_enode->dof; j < cng_enode->n; j++, dof += stride) {
					order1_t order = shapeset->get_order(indices[j]).get_edge_order(iedge);
					int idx = shapeset->get_constrained_edge_index(iedge, ecomp->ori, order, ecomp->part);
					al->add(idx, dof, ecomp->coef);
				}
			}
			else {
				for (int j = 0; j < cng_enode->n; j++) {
					order1_t order = shapeset->get_order(indices[j]).get_edge_order(iedge);
					int idx = shapeset->get_constrained_edge_index(iedge, ecomp->ori, order, ecomp->part);
					al->add(idx, DIRICHLET_DOF, ecomp->coef * cng_enode->bc_proj[j]);
				}
			}
		}
//		printf("  --\n");
		// edge constrained by face
		for (int i = 0; i < enode->face_ncomponents; i++) {
			BaseFaceComponent *fcomp = enode->face_baselist + i;
			FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining edge node
			assert(cng_fnode->ced == false);

//			printf("eid = %d, e = %d, iface = %d, ori = %d, order = %d (%d, %d), part = %d, %d, dir = %d\n", edge_id, elem->id, fcomp->iface, fcomp->ori,
//			       cng_fnode->order, GET_QUAD_ORDER_1(cng_fnode->order), GET_QUAD_ORDER_2(cng_fnode->order),
//			       fcomp->part.horz, fcomp->part.vert, fcomp->dir);

//			int *indices = shapeset->get_face_indices(fcomp->iface, fcomp->ori, cng_fnode->order);
			int *indices = shapeset->get_face_indices(fcomp->iface, 0, cng_fnode->order);
			if (cng_fnode->dof >= 0) {
				for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, dof += stride) {
					order2_t order = shapeset->get_order(indices[j]).get_face_order(fcomp->iface);
					int idx = shapeset->get_constrained_edge_face_index(iedge, fcomp->ori, order, fcomp->part, fcomp->dir);
					al->add(idx, dof, fcomp->coef);
				}
			}
			else {
				for (int j = 0; j < cng_fnode->n; j++) {
					order2_t order = shapeset->get_order(indices[j]).get_face_order(fcomp->iface);
					int idx = shapeset->get_constrained_edge_face_index(iedge, fcomp->ori, order, fcomp->part, fcomp->dir);
					al->add(idx, DIRICHLET_DOF, fcomp->coef * cng_fnode->bc_proj[j]);
				}
			}
		}
	}
	else {
		int *indices = shapeset->get_edge_indices(iedge, ori, enode->order);
		if (enode->dof >= 0) {
			for (int j = 0, dof = enode->dof; j < enode->n; j++, dof += stride) {
				al->add(indices[j], dof, 1.0);
			}
		}
		else if (enode->bc_proj != NULL) {
			for (int j = 0; j < enode->n; j++) {
				double coef = enode->bc_proj[j];
				al->add(indices[j], DIRICHLET_DOF, coef);
			}
		}
	}
}

void Space::get_face_assembly_list(Element *elem, int iface, AsmList *al) {
	Word_t face_id = mesh->get_facet_id(elem, iface);
	FaceData *fnode = fn_data[face_id];

	if (fnode->ced) {
		Facet *facet = mesh->facets[face_id];

		if (fnode->facet_id != -1) {
			FaceData *cng_fnode = fn_data[fnode->facet_id];
//			printf("eid = %d, iface = %d, ori = %d, order = %d (%d, %d), part = (%d, %d)\n", elem->id, iface, fnode->ori, cng_fnode->order,
//			       GET_QUAD_ORDER_1(cng_fnode->order), GET_QUAD_ORDER_2(cng_fnode->order),
//			       fnode->part.horz, fnode->part.vert);

			int *indices = shapeset->get_face_indices(iface, 0, cng_fnode->order);
			if (cng_fnode->dof >= 0) {
				for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, dof += stride) {
					order2_t order = shapeset->get_order(indices[j]).get_face_order(iface);
					int idx = shapeset->get_constrained_face_index(iface, fnode->ori, order, fnode->part);
					assert(dof >= DIRICHLET_DOF && dof < get_dof_count());
					al->add(idx, dof, 1.0);
				}
			}
			else {
				for (int j = 0; j < cng_fnode->n; j++) {
					order2_t order = shapeset->get_order(indices[j]).get_face_order(iface);
					assert(false);
//					int idx = shapeset->get_constrained_face_index(iface, fnode->ori, order, fnode->part);
//					al->add(idx, DIRICHLET_DOF, cng_fnode->bc_proj[j]);
				}
			}
		}
	}
	else {
		int ori = elem->get_face_orientation(iface);

//		printf("elem %d, iface = %d : %d, ori = %d\n", elem->id, iface, face_id, ori);

		int *indices = shapeset->get_face_indices(iface, ori, fnode->order);
		if (fnode->dof >= 0) {
			for (int j = 0, dof = fnode->dof; j < fnode->n; j++, dof += stride)
				al->add(indices[j], dof, 1.0);
		}
		else if (fnode->bc_proj != NULL) {
			for (int j = 0; j < fnode->n; j++) {
				double coef = fnode->bc_proj[j];
				al->add(indices[j], DIRICHLET_DOF, coef);
			}
		}
	}
}

void Space::get_bubble_assembly_list(Element *e, AsmList *al) {
	ElementData *enode = elm_data[e->id];

	int *indices = shapeset->get_bubble_indices(enode->order);
	for (int j = 0, dof = enode->dof; j < enode->n; j++, dof += stride)
		al->add(indices[j], dof, 1.0);
}

// BC ////

void Space::set_bc_info(NodeData *node, EBCType bc, int marker) {
	if (bc == BC_ESSENTIAL || (bc == BC_NATURAL && node->bc_type == BC_NONE)) {
		node->bc_type = bc;
		node->marker = marker;
	}
}

void Space::set_bc_information() {
	FOR_ALL_FACETS(idx, mesh) {
		Facet *facet = mesh->facets.get(idx);
		assert(facet != NULL);

//		if (facet->ractive && facet->lactive && facet->type == Facet::OUTER) {
		if (facet->type == Facet::OUTER) {
			Boundary *bdr = mesh->boundaries[facet->right];
			EBCType bc_type = bc_type_callback(bdr->marker);

			int marker = bdr->marker;
			// set boundary condition for face
			assert(fn_data.exists(idx));
			fn_data[idx]->bc_type = bc_type;
			if (fn_data[idx]->marker == MARKER_UNDEFINED)
				fn_data[idx]->marker = marker;
//			else
//				WARNING("Marker already defined for face #%d.", idx);

			Element *elem = mesh->elements[facet->left];
			int iface = facet->left_face_num;

			// set boundary condition for vertices on the face
			int vtx_num = elem->get_face_num_of_vertices(iface);
			Word_t vtcs[vtx_num];
			elem->get_face_vertices(iface, vtcs);
			for (int i = 0; i < vtx_num; i++) {
				assert(vn_data.exists(vtcs[i]));
				set_bc_info(vn_data[vtcs[i]], bc_type, marker);
			}

			// set boundary condition for edges on the face
			const int *face_edges = elem->get_face_edges(iface);
			for (int i = 0; i < elem->get_face_num_of_edges(iface); i++) {
				int edge_id = mesh->get_edge_id(elem, face_edges[i]);
				if (mesh->edges[edge_id].bnd) {
					assert(en_data.exists(edge_id));
					set_bc_info(en_data[edge_id], bc_type, marker);

//					printf("BC: edgeid = %d > %d\n", edge_id, bc_type);
				}
			}
		}
	}

//	FOR_ALL_EDGES(edge_id, mesh) {
//		if (mesh->edges[edge_id].bnd) {
//			assert(en_data.exists(edge_id));
//			set_bc_info(en_data[edge_id], bc_type, marker);
//
//			printf("BC: edgeid = %d > %d\n", edge_id, bc_type);
//		}
//
//	}
}

// find constraints ///////////////////////////////////////////////////////////

Space::VertexData *Space::create_vertex_node_data(Word_t vid, bool ced) {
	VertexData *vd = vn_data[vid];
	if (vd == NULL) {
		vd = vn_data[vid] = new VertexData;
		MEM_CHECK(vd);
		vd->ced = ced;
		if (ced) {
			vd->baselist = NULL;
			vd->ncomponents = 0;
		}
		else {
			vd->dof = DOF_UNASSIGNED;
//			vd->order = -1;
			vd->n = -1;
		}
	}
	else {
		// not CED and should be => change it (but ont the other way)
		if (!vd->ced && ced) {
			vd->ced = ced;
			vd->baselist = NULL;
			vd->ncomponents = 0;
		}
	}

//	else {
//		if (vd->ced) {
//			free(vd->baselist); vd->baselist = NULL; vd->ncomponents = 0;
//		}
//		vd->ced = ced;
//		if (ced) {
//			vd->baselist = NULL;
//			vd->ncomponents = 0;
//		}
//		else {
//			vd->dof = DOF_UNASSIGNED;
//			vd->order = -1;
//			vd->n = -1;
//		}
//	}

	return vd;
}

Space::EdgeData *Space::create_edge_node_data(Word_t eid, bool ced) {
//	if (eid == 799) {
//		printf("edge %d: %d\n", eid, ced);
//	}

	EdgeData *ed = en_data[eid];
	if (ed == NULL) {
		ed = en_data[eid] = new EdgeData;
		MEM_CHECK(ed);
		ed->ced = ced;
		if (ced) {
			ed->edge_baselist = NULL;
			ed->edge_ncomponents = 0;
			ed->face_baselist = NULL;
			ed->face_ncomponents = 0;
		}
		else {
			ed->order = INVALID_EDGE_ORDER;
			ed->dof = DOF_UNASSIGNED;
			ed->n = -1;
		}
	}
	else {
		// not CED and should be => change it (but ont the other way)
		if (!ed->ced && ced) {
			ed->ced = ced;
			ed->edge_baselist = NULL;
			ed->edge_ncomponents = 0;
			ed->face_baselist = NULL;
			ed->face_ncomponents = 0;
		}

//		if (ed->ced && !ced) {
//			printf("converting ced edge to non-ced one (%d)\n", eid);
//			assert(false);
//		}
	}

//	else {
//		if (ed->ced) {
//			free(ed->edge_baselist); ed->edge_baselist = NULL; ed->edge_ncomponents = 0;
//			free(ed->face_baselist); ed->face_baselist = NULL; ed->face_ncomponents = 0;
//		}
//		ed->ced = ced;
//		if (ced) {
//			ed->edge_baselist = NULL;
//			ed->edge_ncomponents = 0;
//			ed->face_baselist = NULL;
//			ed->face_ncomponents = 0;
//		}
//		else {
//			ed->order = -1;
//			ed->dof = DOF_UNASSIGNED;
//			ed->n = -1;
//		}
//	}

	return ed;
}

Space::FaceData *Space::create_face_node_data(Word_t fid, bool ced) {
//	printf(" - %d => %d\n", fid, ced);

	FaceData *fd = fn_data[fid];
	if (fd == NULL) {
		fd = fn_data[fid] = new FaceData;
		MEM_CHECK(fd);
		fd->ced = ced;
		if (ced) {
			fd->facet_id = INVALID_IDX;
			fd->ori = 0;
			fd->part.horz = 0;
			fd->part.vert = 0;
		}
		else {
//			fd->order = -1;
			fd->dof = DOF_UNASSIGNED;
			fd->n = -1;
		}
	}
	else {
		if (!fd->ced && ced) {
			fd->ced = ced;
			fd->facet_id = INVALID_IDX;
			fd->ori = 0;
			fd->part.horz = 0;
			fd->part.vert = 0;
		}
	}

//	else {
//		if (fd->ced) {
//			fd->facet_id = INVALID_IDX;
//			fd->ori = 0;
//			fd->part.horz = 0;
//			fd->part.vert = 0;
//		}
//		fd->ced = ced;
//		if (ced) {
//			fd->facet_id = INVALID_IDX;
//			fd->ori = 0;
//			fd->part.horz = 0;
//			fd->part.vert = 0;
//		}
//		else {
//			fd->order = -1;
//			fd->dof = DOF_UNASSIGNED;
//			fd->n = -1;
//		}
//	}

	return fd;
}

void Space::fc_face(Word_t eid, int iface, bool ced) {
	Element *elem = mesh->elements[eid];
	// vertices
	int nv = elem->get_face_num_of_vertices(iface);
	Word_t vtcs[nv];
	elem->get_face_vertices(iface, vtcs);

	Word_t fid = mesh->get_facet_id(elem, iface);
	Facet *facet = mesh->facets[fid];

	// set CEDs
	Word_t emp[4], fmp;
	switch (facet->ref_mask) {
		case REFT_QUAD_HORZ:
			emp[0] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[1] = mesh->peek_midpoint(vtcs[3], vtcs[0]);

			// vertices
			create_vertex_node_data(emp[0], ced);
			create_vertex_node_data(emp[1], ced);
			// edges
			create_edge_node_data(mesh->get_edge_id(vtcs[1], emp[0]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[0], vtcs[2]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[3], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[1], vtcs[0]), ced);

			create_edge_node_data(mesh->get_edge_id(vtcs[0], vtcs[1]), false);
			create_edge_node_data(mesh->get_edge_id(emp[0], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[2], vtcs[3]), false);
			break;

		case REFT_QUAD_VERT:
			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[2], vtcs[3]);

			// vertices
			create_vertex_node_data(emp[0], ced);
			create_vertex_node_data(emp[1], ced);
			// edges by edges
			create_edge_node_data(mesh->get_edge_id(vtcs[0], emp[0]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[0], vtcs[1]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[2], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[1], vtcs[3]), ced);

			create_edge_node_data(mesh->get_edge_id(vtcs[0], vtcs[3]), false);
			create_edge_node_data(mesh->get_edge_id(emp[0], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[1], vtcs[2]), false);
			break;

		case REFT_QUAD_BOTH:
			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[2] = mesh->peek_midpoint(vtcs[2], vtcs[3]);
			emp[3] = mesh->peek_midpoint(vtcs[3], vtcs[0]);
			fmp = mesh->peek_midpoint(emp[0], emp[2]);

			// vertices
			for (int iv = 0; iv < Quad::NUM_VERTICES; iv++)
				create_vertex_node_data(emp[iv], ced);
			create_vertex_node_data(fmp, ced);
			// edges
			for (int i = 0; i < Quad::NUM_VERTICES; i++) {
				int i1 = (i + 1) % Quad::NUM_VERTICES;
				create_edge_node_data(mesh->get_edge_id(vtcs[i], emp[i]), ced);
				create_edge_node_data(mesh->get_edge_id(emp[i], vtcs[i1]), ced);
				create_edge_node_data(mesh->get_edge_id(emp[i], fmp), ced);
			}
			break;
	}

	// faces (common for all types of refinements)
	for (int i = 0; i < Facet::MAX_SONS; i++) {
		int sid = facet->sons[i];
		if (sid != INVALID_IDX) create_face_node_data(sid, ced);
	}
}

void Space::fc_face_left(Word_t fid) {
	if (fid == INVALID_IDX) return;

	Facet *facet = mesh->facets[fid];
	fc_face(facet->left, facet->left_face_num, true);
	// recur to sons
	for (int i = 0; i < Facet::MAX_SONS; i++)
		fc_face_left(facet->sons[i]);
}

void Space::fc_face_right(Word_t fid) {
	if (fid == INVALID_IDX) return;

	Facet *facet = mesh->facets[fid];
	fc_face(facet->right, facet->right_face_num, true);
	// recur to sons
	for (int i = 0; i < Facet::MAX_SONS; i++)
		fc_face_right(facet->sons[i]);
}

void Space::fc_element(Word_t idx) {
	if (idx == INVALID_IDX) return;

	Element *elem = mesh->elements[idx];
	for (int iface = 0; iface < elem->get_num_of_faces(); iface++) {
//		printf("elem #%d, iface %d\n", idx, iface);

		Word_t fid = mesh->get_facet_id(elem, iface);
		Facet *facet = mesh->facets[fid];
		assert(facet != NULL);

		// vertices
		int nv = elem->get_face_num_of_vertices(iface);
		Word_t vtcs[nv];
		elem->get_face_vertices(iface, vtcs);
		for (int iv = 0; iv < nv; iv++)
			create_vertex_node_data(vtcs[iv], false);
		// edges
		int ne = elem->get_face_num_of_edges(iface);
		const int *edge_idx = elem->get_face_edges(iface);
		for (int ie = 0; ie < ne; ie++)
			create_edge_node_data(mesh->get_edge_id(elem, edge_idx[ie]), false);

		//
		create_face_node_data(fid, false);

		// handle possible CEDs
		if (facet->type == Facet::INNER) {
			if ((facet->lactive && !facet->ractive) && (facet->right == idx && facet->right_face_num == iface)) {
//				create_face_node_data(fid, true);
				fc_face_right(fid);
			}
			else if ((!facet->lactive && facet->ractive) && (facet->left == idx && facet->left_face_num == iface)) {
//				create_face_node_data(fid, true);
				fc_face_left(fid);
			}
//			else
//				create_face_node_data(fid, false);
		}
//		else
//			create_face_node_data(fid, false);
	}

	// process child elements
	for (int ich = 0; ich < elem->get_num_of_sons(); ich++) {
		Word_t elem_id = elem->get_son(ich);
		fc_element(elem_id);
	}
}

inline void Space::output_component(BaseVertexComponent *&current, BaseVertexComponent *&last, BaseVertexComponent *min, bool add) {
	// if the edge is already in the list, just add half of the other coef
	if (last != NULL && last->dof == min->dof) {
		PRINTF(" * dof already in list (%d, last->coef = % lf, min->coef = % lf.\n", min->dof, last->coef, min->coef);
		if (add) {
//			if (min->dof == DIRICHLET_DOF) last->coef += min->coef; /* do nothing */
			last->coef += min->coef;
		}
		else {
//			if (min->dof == DIRICHLET_DOF) last->coef += 0.5 * min->coef;
//			else
			last->coef += min->coef * 0.5;
		}
		return;
	}

	// output new vertex component
	current->dof = min->dof;
	if (add) {
//		if (min->dof == DIRICHLET_DOF) current->coef = 0.5 * min->coef;
//		else
		current->coef = min->coef;
	}
	else {
//		if (min->dof == DIRICHLET_DOF) current->coef = 0.5 * min->coef;
//		else
		current->coef = min->coef * 0.5;
	}
	last = current++;
}

Space::BaseVertexComponent *Space::merge_baselist(BaseVertexComponent *l1, int n1, BaseVertexComponent *l2, int n2, int &ncomponents, bool add) {
	if (l1 == NULL && l2 == NULL) { ncomponents = 0; return NULL; }
	if (l1 == NULL) { ncomponents = n2; return l2; }
	if (l2 == NULL) { ncomponents = n1; return l1; }

	// estimate the upper bound of the result size
	int max_result = n1 + n2;

	BaseVertexComponent *result = (BaseVertexComponent *) malloc(max_result * sizeof(BaseVertexComponent));
	BaseVertexComponent *current = result;
	BaseVertexComponent *last = NULL;

	// main loop - always output the component with smaller dof so that we get a sorted array
	int i1 = 0, i2 = 0;
	while (i1 < n1 && i2 < n2) {
		if (l1[i1].dof < l2[i2].dof) output_component(current, last, l1 + i1++, add);
		else output_component(current, last, l2 + i2++, add);
	}

	// finish the longer baselist
	while (i1 < n1) output_component(current, last, l1 + i1++, add);
	while (i2 < n2) output_component(current, last, l2 + i2++, add);

	// if we produced less components than we expected, reallocate the resulting array
	// ...this should be OK as we are always shrinking the array so no copying should occur
	ncomponents = current - result;
//	if (ncomponents < max_result) return (BaseVertexComponent *) realloc(result, ncomponents * sizeof(BaseVertexComponent));
//	else
	return result;
}

inline void Space::output_component(BaseEdgeComponent *&current, BaseEdgeComponent *&last, BaseEdgeComponent *min, bool add) {
	// if the edge is already in the list, just add half of the other coef
	if (last != NULL && last->edge_id == min->edge_id) {
		PRINTF(" * edge already in list (%d, last->coef = % lf, min->coef = % lf.\n", min->edge_id, last->coef, min->coef);
		if (add) last->coef += min->coef;
		else last->coef += min->coef * 0.5;
		return;
	}

	// output new edge component
	current->edge_id = min->edge_id;
	current->ori = min->ori;
	current->part = min->part;
	if (add) current->coef = min->coef;
	else current->coef = min->coef * 0.5;
	last = current++;
}

Space::BaseEdgeComponent *Space::merge_baselist(BaseEdgeComponent *l1, int n1, BaseEdgeComponent *l2, int n2, int &ncomponents, bool add) {
	if (l1 == NULL && l2 == NULL) { ncomponents = 0; return NULL; }
	if (l1 == NULL) { ncomponents = n2; return l2; }
	if (l2 == NULL) { ncomponents = n1; return l1; }

	int max_result = n1 + n2;
//	BaseEdgeComponent *result = new BaseEdgeComponent[max_result];
	BaseEdgeComponent *result = (BaseEdgeComponent *) malloc(max_result * sizeof(BaseEdgeComponent));
	BaseEdgeComponent *current = result;
	BaseEdgeComponent *last = NULL;

	// main loop - always output the component with smaller edge_id so that we get a sorted array
	int i1 = 0, i2 = 0;
	while (i1 < n1 && i2 < n2) {
		if (l1[i1].edge_id < l2[i2].edge_id) output_component(current, last, l1 + i1++, add);
		else output_component(current, last, l2 + i2++, add);
	}

	// finish the longer baselist
	while (i1 < n1) output_component(current, last, l1 + i1++, add);
	while (i2 < n2) output_component(current, last, l2 + i2++, add);

	ncomponents = current - result;
//	if (ncomponents < max_result) return (BaseEdgeComponent *) realloc(result, ncomponents * sizeof(BaseEdgeComponent));
//	else
	return result;
}

inline void Space::output_component_over(BaseFaceComponent *&current, BaseFaceComponent *min, BaseFaceComponent *m) {
	if (min != NULL) {
		current->face_id = min->face_id;
		current->ori = min->ori;
		current->iface = min->iface;
		current->part = min->part;
		current->coef = min->coef;
		current->dir = min->dir;
	}
	else if (m != NULL) {
		current->face_id = m->face_id;
		current->ori = m->ori;
		current->iface = m->iface;
		current->part = m->part;
		current->coef = m->coef;
		current->dir = m->dir;
	}
	current++;
}

inline void Space::output_component(BaseFaceComponent *&current, BaseFaceComponent *&last, BaseFaceComponent *min, bool add) {
	if (last != NULL && last->face_id == min->face_id &&
		last->part.vert == min->part.vert && last->part.horz == min->part.horz && last->dir == min->dir) {
		PRINTF(" * face already in list (%d, last->coef = % lf, min->coef = % lf, last->part = %d, min->part = %d)\n",
			min->face_id, last->coef, min->coef, last->part.horz, min->part.horz);
		last->coef += min->coef * 0.5;
		return;
	}

	// output new face component
	current->face_id = min->face_id;
	current->ori = min->ori;
	current->iface = min->iface;
	current->part = min->part;
	current->dir = min->dir;
	if (add) current->coef = min->coef;
	else current->coef = min->coef * 0.5;
	last = current++;
}

Space::BaseFaceComponent *Space::merge_baselist(BaseFaceComponent *l1, int n1, BaseFaceComponent *l2, int n2, int &ncomponents, Word_t fid, bool add) {
	if (l1 == NULL && l2 == NULL) { ncomponents = 0; return NULL; }

	int max_result = n1 + n2;
//	BaseFaceComponent *result = new BaseFaceComponent[max_result];
	BaseFaceComponent *result = (BaseFaceComponent *) malloc(max_result * sizeof(BaseFaceComponent));
	BaseFaceComponent *current = result;
	BaseFaceComponent *last = NULL;

	// main loop - always output the component with smaller face_id so that we get a sorted array
	int i1 = 0, i2 = 0;
	if (add) {
		while (i1 < n1 && i2 < n2) {
			if (l1[i1].face_id < l2[i2].face_id) output_component(current, last, l1 + i1++, add);
			else output_component(current, last, l2 + i2++, add);
		}

		// finish the longer baselist
		while (i1 < n1) output_component(current, last, l1 + i1++, add);
		while (i2 < n2) output_component(current, last, l2 + i2++, add);

		ncomponents = current - result;
	}
	else {
		while (i1 < n1 && i2 < n2) {
			if (l1[i1].face_id == fid) i1++;
			else if (l2[i2].face_id == fid) i2++;
			else {
				if ((l1[i1].face_id < l2[i2].face_id) ||
					(l1[i1].face_id == l2[i2].face_id && l1[i1].part.horz < l2[i2].part.horz) ||
					(l1[i1].face_id == l2[i2].face_id && l1[i1].part.horz == l2[i2].part.horz && l1[i1].part.vert < l2[i2].part.vert) ||
					(l1[i1].face_id == l2[i2].face_id && l1[i1].part.horz == l2[i2].part.horz && l1[i1].part.vert == l2[i2].part.vert && l1[i1].dir < l2[i2].dir))
					output_component(current, last, l1 + i1++, add);
				else
					output_component(current, last, l2 + i2++, add);
			}
		}

		// finish the longer baselist
		while (i1 < n1) {
			if (l1[i1].face_id == fid) i1++;
			else output_component(current, last, l1 + i1++, add);
		}
		while (i2 < n2) {
			if (l2[i2].face_id == fid) i2++;
			else output_component(current, last, l2 + i2++, add);
		}

		ncomponents = current - result;
	}

//	if (ncomponents < max_result) return (BaseFaceComponent *) realloc(result, ncomponents * sizeof(BaseFaceComponent));
//	else return result;
	return result;
}

/// @param[in] vtx1 - vertex 1
/// @param[in] vtx2 - vertex 2
///
void Space::calc_vertex_vertex_ced(Word_t vtx1, Word_t vtx2) {
	assert(vtx1 != INVALID_IDX);
	assert(vtx2 != INVALID_IDX);
	VertexData *vd[] = { vn_data[vtx1], vn_data[vtx2] };
	Word_t mid_pt = mesh->peek_midpoint(vtx1, vtx2);
	assert(mid_pt != INVALID_IDX);

	PRINTF("calc vertex/vertex #%d\n", mid_pt);

	VertexData *vd_mid = vn_data[mid_pt];
	assert(vd_mid != NULL);

    BaseVertexComponent *bl[2], dummy_bl[2];	// base lists of vtx1 and vtx2
	int nc[2] = { 0, 0 }; // number of components of bl[0] and bl[1]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 2; k++) {
		if (vd[k]->ced) {
			bl[k] = vd[k]->baselist;
			nc[k] = vd[k]->ncomponents;
		}
		else {	// make up an artificial baselist
			dummy_bl[k].dof = vd[k]->dof;
			dummy_bl[k].coef = (vd[k]->dof >= 0) ? 1.0 : vd[k]->bc_proj;
			bl[k] = &dummy_bl[k];
			nc[k] = 1;
		}
	}

//	PRINTF("--\n");
//	for (int i = 0; i < nc[0]; i++)
//		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, bl[0][i].dof, bl[0][i].coef);
//	PRINTF("--\n");
//	for (int i = 0; i < nc[1]; i++)
//		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, bl[1][i].dof, bl[1][i].coef);
//	PRINTF("--\n");

	int ncomp = 0;
	vd_mid->baselist = merge_baselist(bl[0], nc[0], bl[1], nc[1], ncomp, false);
	vd_mid->ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++) {
		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd_mid->baselist[i].dof, vd_mid->baselist[i].coef);
	}
}

/// @param[in] vtx1 - vertex 1
/// @param[in] vtx2 - vertex 2
///
void Space::calc_mid_vertex_vertex_ced(Word_t mid, Word_t vtx1, Word_t vtx2, Word_t vtx3, Word_t vtx4) {
	assert(vtx1 != INVALID_IDX);
	assert(vtx2 != INVALID_IDX);
	assert(vtx3 != INVALID_IDX);
	assert(vtx4 != INVALID_IDX);
	VertexData *vd[] = { vn_data[vtx1], vn_data[vtx2], vn_data[vtx3], vn_data[vtx4] };

	PRINTF("calc mid vertex/vertex #%d (%d, %d, %d, %d)\n", mid, vtx1, vtx2, vtx3, vtx4);

	VertexData *vd_mid = vn_data[mid];
	assert(vd_mid != NULL);

    BaseVertexComponent *bl[4], dummy_bl[4];	// base lists of vtx1-4
	int nc[4] = { 0, 0, 0, 0 }; // number of components of bl[0-3]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 4; k++) {
		if (vd[k]->ced) {
			bl[k] = vd[k]->baselist;
			nc[k] = vd[k]->ncomponents;
		}
		else {	// make up an artificial baselist
			dummy_bl[k].dof = vd[k]->dof;
			dummy_bl[k].coef = (vd[k]->dof >= 0) ? 1.0 : vd[k]->bc_proj;
			bl[k] = &dummy_bl[k];
			nc[k] = 1;
		}
	}

	int tmp_nc[] = { 0, 0 };
	BaseVertexComponent *tmp_bl[2];
	tmp_bl[0] = merge_baselist(bl[0], nc[0], bl[2], nc[2], tmp_nc[0], false);
	tmp_bl[1] = merge_baselist(bl[1], nc[1], bl[3], nc[3], tmp_nc[1], false);

//	PRINTF("--\n");
//	for (int i = 0; i < tmp_nc[0]; i++) {
//		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, tmp_bl[0][i].dof, tmp_bl[0][i].coef);
//	}
//	PRINTF("--\n");
//	for (int i = 0; i < tmp_nc[1]; i++) {
//		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, tmp_bl[1][i].dof, tmp_bl[1][i].coef);
//	}
//	PRINTF("--\n");

	int ncomp = 0;
	vd_mid->baselist = merge_baselist(tmp_bl[0], tmp_nc[0], tmp_bl[1], tmp_nc[1], ncomp, false);
	vd_mid->ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++)
		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd_mid->baselist[i].dof, vd_mid->baselist[i].coef);
}

void Space::calc_vertex_edge_ced(Word_t vtx, Word_t eid, int ori, int part) {
	PRINTF("calc vertex/edge #%d\n", vtx);

	PRINTF(" - eid = %d, part = %d, ori = %d\n", eid, part, ori);
//	printf(" - eid = %d, part = %d, ori = %d\n", eid, part, ori);

	assert(eid != INVALID_IDX);
	EdgeData *ed = en_data[eid];
	assert(ed != NULL);

	assert(vtx != INVALID_IDX);
	VertexData *vd = vn_data[vtx];
	assert(vd != NULL);

	double lo, hi;
	int ncomp = 0;
	if (ed->ced) {
//		printf("B\n");

		// count the number of components to merge
		int nc = 0;
		for (int i = 0; i < ed->edge_ncomponents; i++) {
			BaseEdgeComponent *ecomp = ed->edge_baselist + i;
			EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node
			nc += cng_enode->n;
		}
		for (int i = 0; i < ed->face_ncomponents; i++) {
			BaseFaceComponent *fcomp = ed->face_baselist + i;
			FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining face node
			nc += cng_fnode->n;
		}
//		BaseVertexComponent *baselist = new BaseVertexComponent[nc];
		BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(nc * sizeof(BaseVertexComponent));
//		printf("nc = %d\n", nc);

//		get_interval_part(part, lo, hi);
//		double mid = (lo + hi) * 0.5;
//		printf("   - part = %d (% lf, % lf)\n", part, lo, hi);

		int nci = 0;
		// update the edge part
		for (int i = 0; i < ed->edge_ncomponents; i++) {
			BaseEdgeComponent *ecomp = ed->edge_baselist + i;
			EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

//			get_interval_part(ecomp->part.part, lo, hi);
//			double mid = (lo + hi) * 0.5;

//			double edge_coef = get_edge_coef(face_to_edge_part(part));
//			printf("   - ECOMP: edge_id = %d, part = %d, coef = %lf | edge_coef = %lf, mid = %lf\n",
//				ecomp->edge_id, ecomp->part.part, ecomp->coef, edge_coef, mid);
//			printf("   - ECOMP: edge_id = %d, part = %d, coef = %lf\n",
//				ecomp->edge_id, ecomp->part.part, ecomp->coef);

			int *indices = shapeset->get_edge_indices(0, ecomp->ori, cng_enode->order);
			for (int j = 0, dof = cng_enode->dof; j < cng_enode->n; j++, nci++) {
				order1_t order = shapeset->get_order(indices[j]).get_edge_order(0);
				int idx = shapeset->get_constrained_edge_index(0, ecomp->ori, order, ecomp->part);
//				printf("   - dof = %d, coef = % lf, fn = % lf, (order = %d)\n",
//				       dof, ecomp->coef, shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0), order);
				baselist[nci].dof = dof;
				baselist[nci].coef = (ecomp->coef) * shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0);
				if (cng_enode->dof == DIRICHLET_DOF) baselist[nci].coef *= cng_enode->bc_proj[j];
				else dof += stride;
			}
		}
		// update the face part
		for (int i = 0; i < ed->face_ncomponents; i++) {
			BaseFaceComponent *fcomp = ed->face_baselist + i;
			FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining face node

//			double h_mid, v_mid;
//			if (fcomp->part.ori == PART_ORI_VERT) {
//				get_edge_part(fcomp->part.epart, h_mid);
//				get_interval_part(fcomp->part.fpart, lo, hi);
//				v_mid = (lo + hi) * 0.5;
//			}
//			else {
//				get_interval_part(fcomp->part.fpart, lo, hi);
//				h_mid = (lo + hi) * 0.5;
//				get_edge_part(fcomp->part.epart, v_mid);
//			}

//			get_interval_part(fcomp->part.fpart, lo, hi);
//			double mid = (lo + hi) * 0.5;
//
//			double h_mid;
//			get_edge_part(fcomp->part.epart, h_mid);
//			PRINTF("   - face_id = %d, part = (f = %d, e = %d, ori = %d), ori = %d, coef = %lf, pt = (% lf, %lf)\n",
//				fcomp->face_id, fcomp->part.fpart, fcomp->part.epart, fcomp->part.ori, fcomp->ori, fcomp->coef, h_mid, v_mid);
//			PRINTF("   - FCOMP: face_id = %d, part = (f = %d, e = %d, ori = %d), ori = %d, coef = %lf\n",
//				fcomp->face_id, fcomp->part.fpart, fcomp->part.epart, fcomp->part.ori, fcomp->ori, fcomp->coef);
//			printf("   - FCOMP: face_id = %d, part = (%d, %d), dir = %d, ori = %d, coef = %lf\n",
//				fcomp->face_id, fcomp->part.horz, fcomp->part.vert, fcomp->dir, fcomp->ori, fcomp->coef);

			int *indices = shapeset->get_face_indices(2, fcomp->ori, cng_fnode->order);
			for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, nci++) {
				// FIXME: Hex-specific
				order2_t order = shapeset->get_order(indices[j]).get_face_order(2);
				int idx = shapeset->get_constrained_edge_face_index(0, fcomp->ori, order, fcomp->part, fcomp->dir);

//				PRINTF("   - dof = %d, coef = % lf, fn = % lf, (order = %d)\n",
//				       dof, fcomp->coef, shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0), order);
//				printf("   - dof = %d, coef = % lf, fn = % lf, (order = %s)\n",
//					   dof, fcomp->coef, shapeset->get_fn_value(idx, 0.0, -1.0, -1.0, 0), order.str());

				baselist[nci].dof = dof;
				baselist[nci].coef = (fcomp->coef) * shapeset->get_fn_value(idx, 0.0, -1.0, -1.0, 0);
				if (cng_fnode->dof == DIRICHLET_DOF) baselist[nci].coef *= cng_fnode->bc_proj[j];
				else dof += stride;
			}
		}

//		PRINTF("--\n");
//		for (int i = 0; i < vd->ncomponents; i++)
//			PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);
//		PRINTF("--\n");
//		for (int i = 0; i < nc; i++)
//			PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, baselist[i].dof, baselist[i].coef);
//		PRINTF("--\n");

//		printf("--\n");
//		for (int i = 0; i < vd->ncomponents; i++)
//			printf(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);
//		printf("--\n");
//		for (int i = 0; i < nc; i++)
//			printf(" - [%d]: dof = %d, coef = %lf\n", i, baselist[i].dof, baselist[i].coef);
//		printf("--\n");

		BaseVertexComponent *del = vd->baselist;
		vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, nc, ncomp, true);
		vd->ncomponents = ncomp;

//		free(del);
//		free(baselist);
	}
	else {
		get_interval_part(part, lo, hi);
		double mid = (lo + hi) * 0.5;

/*		int *indices = shapeset->get_edge_indices(0, ori, ed->order);
		BaseVertexComponent *baselist;
		if (ed->dof == DIRICHLET_DOF) {
			baselist = new BaseVertexComponent[0];
			baselist[0].coef = 0.0;
			baselist[0].dof = DIRICHLET_DOF;

			for (int j = 0; j < ed->n; j++)
				baselist[0].coef += ed->bc_proj[j] * shapeset->get_fn_value(indices[j], mid, -1.0, -1.0, 0);

			vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, 1, ncomp, true);
			vd->ncomponents = ncomp;
		}
		else {
			baselist = new BaseVertexComponent[ed->n];

			for (int j = 0, dof = ed->dof; j < ed->n; j++, dof += stride) {
				baselist[j].dof = dof;
				baselist[j].coef = shapeset->get_fn_value(indices[j], mid, -1.0, -1.0, 0);
			}

			vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, ed->n, ncomp, true);
			vd->ncomponents = ncomp;
		}
*/
//		BaseVertexComponent *baselist = new BaseVertexComponent[ed->n];
		BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(ed->n * sizeof(BaseVertexComponent));

		int *indices = shapeset->get_edge_indices(0, ori, ed->order);
		for (int j = 0, dof = ed->dof; j < ed->n; j++) {
			baselist[j].dof = dof;
			baselist[j].coef = shapeset->get_fn_value(indices[j], mid, -1.0, -1.0, 0);
			if (ed->dof == DIRICHLET_DOF) baselist[j].coef *= ed->bc_proj[j];
			else dof += stride;
		}

		BaseVertexComponent *del = vd->baselist;
		vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, ed->n, ncomp, true);
		vd->ncomponents = ncomp;

//		free(del);
//		free(baselist);
	}

	for (int i = 0; i < vd->ncomponents; i++) {
//		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);
//		printf(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);

//		if (vd->baselist[i].dof == -1)
//			vd->bc_proj = vd->baselist[i].coef;
	}
}

void Space::calc_mid_vertex_edge_ced(Word_t vtx, Word_t fmp, Word_t eid, int ori, int part) {
	PRINTF("calc mid vertex/edge #%d, [%d | %d]\n", vtx, eid, fmp);

	assert(eid != INVALID_IDX);

	EdgeData *ed = en_data[eid];
    BaseEdgeComponent *ebl, edummy_bl;
	int enc = 0;
	BaseFaceComponent *fbl, fdummy_bl;
	int fnc = 0;

	// get baselists of edge node[0] and edge node[1]
	if (ed->ced) {
		ebl = ed->edge_baselist;
		enc = ed->edge_ncomponents;

		fbl = ed->face_baselist;
		fnc = ed->face_ncomponents;
	}
	else {	// make up an artificial baselist (we care about edge id and coef
		edummy_bl.edge_id = eid;
		edummy_bl.ori = ori;
		edummy_bl.part.part = part;
		edummy_bl.coef = 1.0;

		ebl = &edummy_bl;
		enc = 1;

		fbl = NULL;
		fnc = 0;
	}

	int en_comp = enc;
	BaseEdgeComponent *edge_baselist = ebl;

	int fn_comp = fnc;
	BaseFaceComponent *face_baselist = fbl;

	for (int i = 0; i < en_comp; i++) {
		BaseEdgeComponent ec = edge_baselist[i];
		PRINTF(" - [%d]: edge_id = %d, ori = %d, part = %d, coef = %lf\n",
			i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
	}
	for (int i = 0; i < fn_comp; i++) {
		BaseFaceComponent fc = face_baselist[i];
//		PRINTF(" - [%d]: face_id = %d, part = (f = %d, e = %d, ori = %d), ori = %d, coef = %lf\n",
//			i, fc.face_id, fc.part.fpart, fc.part.epart, fc.part.ori, fc.ori, fc.coef);
		PRINTF(" - [%d]: face_id = %d, part = (%d, %d), dir = %d, ori = %d, coef = %lf\n",
			i, fc.face_id, fc.part.horz, fc.part.vert, fc.dir, fc.ori, fc.coef);
	}

//	printf(" fc0 = %d\n", ed[0]->face_ncomponents);
//	printf(" fc1 = %d\n", ed[1]->face_ncomponents);
//
//	printf("---\n");

	// -- /////////////////

	// count the number of components to update
	int nc = 0;
	for (int i = 0; i < en_comp; i++) {
		BaseEdgeComponent *ecomp = edge_baselist + i;
		EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node
		nc += cng_enode->n;
		PRINTF(" - enc = %d\n", cng_enode->n);
	}
	for (int i = 0; i < fn_comp; i++) {
		BaseFaceComponent *fcomp = face_baselist + i;
		FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining face node
		nc += cng_fnode->n;
		PRINTF(" - fnc = %d\n", cng_fnode->n);
	}

//	BaseVertexComponent *baselist = new BaseVertexComponent[nc];
	BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(nc * sizeof(BaseVertexComponent));
	PRINTF(" - nc = %d\n", nc);

	assert(vtx != INVALID_IDX);
	VertexData *vd = vn_data[vtx];
	assert(vd != NULL);

//	double lo, hi;
//	double mid;

	int nci = 0;
	for (int i = 0; i < en_comp; i++) {
		BaseEdgeComponent *ecomp = edge_baselist + i;
		EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

//		get_interval_part(ecomp->part.part, lo, hi);
//		mid = (lo + hi) * 0.5;

//		PRINTF("   - ECOMP: edge_id = %d, part = %d, coef = %lf, mid = %lf, order = %d\n",
//			ecomp->edge_id, ecomp->part.part, ecomp->coef, mid, cng_enode->order);
		PRINTF("   - ECOMP: edge_id = %d, part = %d, coef = %lf, order = %d\n",
			ecomp->edge_id, ecomp->part.part, ecomp->coef, cng_enode->order);

		int *indices = shapeset->get_edge_indices(0, ecomp->ori, cng_enode->order);		// iedge bude 0 (?)
		for (int j = 0, dof = cng_enode->dof; j < cng_enode->n; j++, nci++) {
			int order = shapeset->get_order(indices[j]).get_edge_order(0);
			int idx = shapeset->get_constrained_edge_index(0, ecomp->ori, order, ecomp->part);
//			PRINTF("   - dof = %d, coef = % lf, fn = % lf, (order = %d)\n",
//			       dof, ecomp->coef, shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0), order);

			baselist[nci].dof = dof;
			baselist[nci].coef = ecomp->coef * shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0);
			if (cng_enode->dof == DIRICHLET_DOF) baselist[nci].coef *= cng_enode->bc_proj[j];
			else dof += stride;
		}
	}
	for (int i = 0; i < fn_comp; i++) {
		BaseFaceComponent *fcomp = face_baselist + i;
		FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining edge node

//		double h_mid, v_mid;
//		if (fcomp->part.ori == PART_ORI_VERT) {
//			get_edge_part(fcomp->part.epart, h_mid);
//			get_interval_part(fcomp->part.fpart, lo, hi);
//			v_mid = (lo + hi) * 0.5;
//		}
//		else {
//			get_interval_part(fcomp->part.fpart, lo, hi);
//			h_mid = (lo + hi) * 0.5;
//			get_edge_part(fcomp->part.epart, v_mid);
//		}

//		PRINTF("   - FCOMP = %d, part = (f = %d, e = %d, ori = %d), ori = %d, coef = %lf, pt = (% lf, %lf)\n",
//			fcomp->face_id, fcomp->part.fpart, fcomp->part.epart, fcomp->part.ori, fcomp->ori, fcomp->coef, h_mid, v_mid);
		PRINTF("   - FCOMP = %d, part = (%d, %d), dir = %d, ori = %d, coef = %lf\n",
			fcomp->face_id, fcomp->part.horz, fcomp->part.vert, fcomp->dir, fcomp->ori, fcomp->coef);

		int *indices = shapeset->get_face_indices(2, fcomp->ori, cng_fnode->order);
		for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, nci++) {
			// FIXME: Hex-specific
			order2_t order = shapeset->get_order(indices[j]).get_face_order(2);
			int idx = shapeset->get_constrained_edge_face_index(0, fcomp->ori, order, fcomp->part, fcomp->dir);

//			PRINTF("   - dof = %d, coef = % lf, fn = % lf, (order = %d)\n",
//			       dof, fcomp->coef, shapeset->get_fn_value(idx, 0.0, -1.0, -1.0, 0), order);

			baselist[nci].dof = dof;
			baselist[nci].coef = fcomp->coef * shapeset->get_fn_value(idx, 0.0, -1.0, -1.0, 0);
			if (cng_fnode->dof == DIRICHLET_DOF) baselist[nci].coef *= cng_fnode->bc_proj[j];
			else dof += stride;
		}
	}

	int ncomp = 0;
	BaseVertexComponent *del = vd->baselist;
	vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, nc, ncomp, true);
	vd->ncomponents = ncomp;
//	free(del);

	for (int i = 0; i < ncomp; i++)
		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);

	{
		assert(fmp != INVALID_IDX);
		VertexData *fmp_vd = vn_data[fmp];
		assert(fmp_vd != NULL);

		//
		for (int i = 0; i < nc; i++)
			baselist[i].coef *= 0.5;

		BaseVertexComponent *bl, dummy_bl;	// base lists of vtx1-4
		int fnc = 0; // number of components of bl[0-3]

		// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
		if (fmp_vd->ced) {
			bl = fmp_vd->baselist;
			fnc = fmp_vd->ncomponents;
		}
		else {	// make up an artificial baselist
			dummy_bl.dof = fmp_vd->dof;
			dummy_bl.coef = (fmp_vd->dof >= 0) ? 1.0 : fmp_vd->bc_proj;
			bl = &dummy_bl;
			fnc = 1;
		}

		fmp_vd->baselist = merge_baselist(bl, fnc, baselist, nc, ncomp, true);
		fmp_vd->ncomponents = ncomp;
	}

//	free(baselist);

//	delete [] edge_baselist;
}


/// @param vtx - id of the vertex node for which we calculate the constrain
/// @param fid - id of the constraining facet
/// @param ori - the orientation of the constraining facet
void Space::calc_vertex_face_ced(Word_t vtx, Word_t fid, int ori, int iface, int hpart, int vpart) {
	PRINTF("calc vertex/face #%d\n", vtx);

	FaceData *fd = fn_data[fid];
	assert(fd != NULL);

	VertexData *vd = vn_data[vtx];
	assert(vd != NULL);

	double h_lo, h_hi, v_lo, v_hi;
	get_interval_part(hpart, h_lo, h_hi);
	get_interval_part(vpart, v_lo, v_hi);
	double h_mid = (h_lo + h_hi) * 0.5;
	double v_mid = (v_lo + v_hi) * 0.5;

	PRINTF(" - fid = %d, ori = %d, iface = %d, part = (%d, %d), comp = %d\n", fid, ori, iface, hpart, vpart, vd->ncomponents);

	if (fd->ced) {
		EXIT(ERR_FAILURE, "Unusual vertex/face CED situation, please report.");
	}
	else {
//		BaseVertexComponent *baselist = new BaseVertexComponent[fd->n];
		BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(fd->n * sizeof(BaseVertexComponent));

		int *indices = shapeset->get_face_indices(2, ori, fd->order);
		for (int j = 0, dof = fd->dof; j < fd->n; j++) {
			// FIXME: hex-specific
			order2_t order = shapeset->get_order(indices[j]).get_face_order(2);
			Part part;
			part.horz = hpart;
			part.vert = vpart;
			int idx = shapeset->get_constrained_face_index(2, ori, order, part);

			baselist[j].dof = dof;
			baselist[j].coef = shapeset->get_fn_value(idx, 0.0, -1.0, 0.0,  0);
			if (fd->dof == DIRICHLET_DOF) baselist[j].coef *= fd->bc_proj[j];
			else dof += stride;
			PRINTF(" - [%d]: dof = %d, coef = %lf\n", j, baselist[j].dof, baselist[j].coef);
		}

		int ncomp = 0;
		BaseVertexComponent *del = vd->baselist;
		vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, fd->n, ncomp, true);
		vd->ncomponents = ncomp;

		PRINTF("--\n");
		for (int i = 0; i < vd->ncomponents; i++)
			PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);

//		free(del);
//		free(baselist);
	}
}

/// @param[in] seid - small edge id
/// @param[in] eid - big edge id
/// @param[in] ori - prt orientation
/// @param[in] epart - part of the edge
/// @param[in] part - part of the edge
void Space::calc_edge_edge_ced(Word_t seid, Word_t eid, int ori, int epart, int part) {
	PRINTF("calc edge/edge #%d, #%d\n", seid, eid);

//	if (seid == 100) printf("********\n");

	assert(eid != INVALID_IDX);
	EdgeData *cng_ed = en_data[eid]; // constraining edge
	assert(cng_ed != NULL);

	assert(seid != INVALID_IDX);
	EdgeData *ed = en_data[seid];
	assert(ed != NULL);

	PRINTF(" *** part = %d, cng_eid = %d\n", part, eid);

	if (cng_ed->ced) {
		PRINTF(" - CED version\n");
		// use the baselist from the "parent" edge and update the part and ori (?)
		int ncomp = cng_ed->edge_ncomponents;
//		BaseEdgeComponent *edge_bl = new BaseEdgeComponent[ncomp];
		BaseEdgeComponent *edge_bl = (BaseEdgeComponent *) malloc(ncomp * sizeof(BaseEdgeComponent));
		for (int i = 0; i < ncomp; i++) {
			edge_bl[i] = cng_ed->edge_baselist[i];
			edge_bl[i].part.part = combine_face_part(edge_bl[i].part.part, epart);
		}
		ed->edge_baselist = edge_bl;
		ed->edge_ncomponents = ncomp;

		// face components
		ncomp = cng_ed->face_ncomponents;
//		BaseFaceComponent *face_bl = new BaseFaceComponent[ncomp];
		BaseFaceComponent *face_bl = (BaseFaceComponent *) malloc(ncomp * sizeof(BaseFaceComponent));
		for (int i = 0; i < ncomp; i++) {
			face_bl[i] = cng_ed->face_baselist[i];

//			BaseFaceComponent fc = face_bl[i];
//			PRINTF(" - [%d]: face_id = %d, ori = %d, part = (%d, %d), dir = %d, coef = %lf\n",
//				i, fc.face_id, fc.ori, fc.part.horz, fc.part.vert, fc.dir, fc.coef);
//			printf("---------------------------------------------------\n");
//			if (face_bl[i].dir == PART_ORI_VERT)
//				printf("comb: %d, %d\n", face_bl[i].part.vert, epart);
//			else
//				printf("comb: %d, %d\n", face_bl[i].part.horz, epart);
//			printf("---------------------------------------------------\n");

			if (face_bl[i].dir == PART_ORI_VERT) face_bl[i].part.vert = combine_face_part(face_bl[i].part.vert, epart);
			else face_bl[i].part.horz = combine_face_part(face_bl[i].part.horz, epart);
//			face_bl[i].part.fpart = combine_face_part(face_bl[i].part.fpart, epart);
//			assert(false);
		}
		ed->face_baselist = face_bl;
		ed->face_ncomponents = ncomp;

		for (int i = 0; i < ed->edge_ncomponents; i++) {
			BaseEdgeComponent ec = ed->edge_baselist[i];
			PRINTF(" - [%d]: edge_id = %d, ori = %d, part = %d, coef = %lf\n",
				i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
		}

		for (int i = 0; i < ed->face_ncomponents; i++) {
			BaseFaceComponent fc = ed->face_baselist[i];
			PRINTF(" - [%d]: face_id = %d, ori = %d, part = (%d, %d), dir = %d, coef = %lf\n",
				i, fc.face_id, fc.ori, fc.part.horz, fc.part.vert, fc.dir, fc.coef);
		}

	}
	else {
		int nc = 1;
//		BaseEdgeComponent *baselist = new BaseEdgeComponent[nc];
		BaseEdgeComponent *baselist = (BaseEdgeComponent *) malloc(nc * sizeof(BaseEdgeComponent));
		baselist[0].edge_id = eid;
		baselist[0].ori = ori;
		baselist[0].part.part = part;
		baselist[0].coef = 1.0;

		BaseEdgeComponent *del = ed->edge_baselist;
		int ncomp = 0;
		ed->edge_baselist = merge_baselist(ed->edge_baselist, ed->edge_ncomponents, baselist, nc, ncomp, false);
		ed->edge_ncomponents = ncomp;

		for (int i = 0; i < ncomp; i++) {
			BaseEdgeComponent ec = ed->edge_baselist[i];
			PRINTF(" - [%d]: edge_id = %d, ori = %d, part = %d, coef = %lf\n",
				i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
		}

//		free(del);
//		free(baselist);
	}
}

void Space::calc_mid_edge_edge_ced(Word_t meid, Word_t eid[], int ori[], int epart, int part) {
	PRINTF("calc mid edge/edge #%d\n", meid);

	assert(eid[0] != INVALID_IDX);
	assert(eid[1] != INVALID_IDX);

	assert(meid != INVALID_IDX);
	EdgeData *mid_ed = en_data[meid];
	assert(mid_ed != NULL);

	EdgeData *ed[] = { en_data[eid[0]], en_data[eid[1]] };
    BaseEdgeComponent *bl[2], dummy_bl[2];		// base lists of eid1 and eid2
	int nc[2] = { 0, 0 };						// number of components of bl[0] and bl[1]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 2; k++) {
		if (ed[k]->ced) {
			PRINTF(" - CED version\n");
			// use the baselist from the "parent" edge and update the part and ori (?)
			int ncomp = ed[k]->edge_ncomponents;
//			BaseEdgeComponent *edge_bl = new BaseEdgeComponent[ncomp];
			BaseEdgeComponent *edge_bl = (BaseEdgeComponent *) malloc(ncomp * sizeof(BaseEdgeComponent));
			for (int i = 0; i < ncomp; i++) {
				edge_bl[i] = ed[k]->edge_baselist[i];
				edge_bl[i].part.part = combine_face_part(edge_bl[i].part.part, epart);
			}

			bl[k] = edge_bl;
			nc[k] = ncomp;

			// FIXME: edge_bl was allocated, but not deallocated
		}
		else {	// make up an artificial baselist
			dummy_bl[k].edge_id = eid[k];
			dummy_bl[k].ori = ori[k];
			dummy_bl[k].part.part = part;
			dummy_bl[k].coef = 1.0;

			bl[k] = &dummy_bl[k];
			nc[k] = 1;
		}
	}

	int ncomp = 0;
	mid_ed->edge_baselist = merge_baselist(bl[0], nc[0], bl[1], nc[1], ncomp, false);
	mid_ed->edge_ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++) {
		BaseEdgeComponent ec = mid_ed->edge_baselist[i];
		PRINTF(" - [%d]: edge_id = %d, ori = %d, part = %d, coef = %lf\n",
			i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
	}
}


/// @param[in] eid - edge id
/// @param[in] fid - constraining face id
/// @param[in] ori - orientation of the constraining face
/// @param[in] part_ori - part ori
/// @param[in] fpart - face part
/// @param[in] epart - edge part
void Space::calc_edge_face_ced(Word_t mid_eid, Word_t eid[], Word_t fid, int ori, int iface, int part_ori, int fpart, int epart) {
	PRINTF("calc edge/face #%d\n", mid_eid);

	assert(fid != INVALID_IDX);
	FaceData *cng_fnode = fn_data[fid]; // constraining face
	assert(cng_fnode != NULL);

	assert(mid_eid != INVALID_IDX);
	EdgeData *mid_ed = en_data[mid_eid];
	assert(mid_ed != NULL);

	PRINTF(" - n = %d | face_id = %d, ori = %d, iface = %d, part (%d, %d), dir = %d\n",
		cng_fnode->n, fid, ori, iface, fpart, epart, part_ori);

	EdgeData *ed[] = { en_data[eid[0]], en_data[eid[1]] };
    BaseFaceComponent *bl[2];					// base lists of eid1 and eid2
	int nc[2] = { 0, 0 };						// number of components of bl[0] and bl[1]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 2; k++) {
		if (ed[k]->ced) {
			bl[k] = ed[k]->face_baselist;
			nc[k] = ed[k]->face_ncomponents;
		}
		else {
			bl[k] = NULL;
			nc[k] = 0;
		}
	}

	int mnc = 0;
	BaseFaceComponent *mbl = merge_baselist(bl[0], nc[0], bl[1], nc[1], mnc, fid, false);

	BaseFaceComponent dummy_bl;
	dummy_bl.face_id = fid;
	dummy_bl.ori = ori;
	dummy_bl.iface = iface;
//	dummy_bl.part.ori = part_ori;
//	dummy_bl.part.fpart = fpart;
//	dummy_bl.part.epart = epart;
	dummy_bl.part.horz = fpart;
	dummy_bl.part.vert = epart;
	dummy_bl.dir = part_ori;
	dummy_bl.coef = 1.0;

	int ncomp = 0;
	mid_ed->face_baselist = merge_baselist(&dummy_bl, 1, mbl, mnc, ncomp, fid, true);
	mid_ed->face_ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++) {
		BaseFaceComponent fc = mid_ed->face_baselist[i];
//		PRINTF(" - [%d]: face_id = %d, ori = %d, iface = %d, part = (ori = %d, fpart = %d, epart = %d), coef = %lf\n",
//			i, fc.face_id, fc.ori, fc.iface, fc.part.ori, fc.part.fpart, fc.part.epart, fc.coef);
		PRINTF(" - [%d]: face_id = %d, ori = %d, iface = %d, part = (%d, %d), dir = %d, coef = %lf\n",
			i, fc.face_id, fc.ori, fc.iface, fc.part.horz, fc.part.vert, fc.dir, fc.coef);
	}
}

/// @param[in] sfid - small facet id (constrained)
/// @param[in] fid - constraining facet id
/// @param[in] ori - orientation of the constraining face
/// @param[in] hpart - horizontal part
/// @param[in] vpart - vertical part
void Space::calc_face_face_ced(Word_t sfid, Word_t fid, int ori, int hpart, int vpart) {
	PRINTF("calc face/face #%d\n", sfid);

	FaceData *fd = fn_data[sfid];
	assert(fd != NULL);
	fd->facet_id = fid;
	fd->ori = ori;
	fd->part.horz = hpart;
	fd->part.vert = vpart;

	PRINTF(" - part = (%d, %d)\n", hpart, vpart);
}

void Space::uc_edge(Word_t eid, int iedge) {
	Element *elem = mesh->elements[eid];

	Word_t edge_id = mesh->get_edge_id(elem, iedge);
	assert(edge_id != INVALID_IDX);
	if (!ei_data.exists(edge_id)) return;

//	printf("*");

	EdgeInfo *ei = ei_data[edge_id];
	assert(ei != NULL);

//	Facet *facet = mesh->facets[fid];
//	assert(facet != NULL);

//	printf("element #%d, face %d, ori = %d\n", eid, iface, elem->get_face_orientation(iface));
//	printf("element #%d, edge %d\n", eid, iedge);

	// vertices
	Word_t vtcs[Edge::NUM_VERTICES];
	elem->get_edge_vertices(iedge, vtcs);

	Word_t seid, sub_eid[2];
	EdgeInfo *sei, *sub_ei[2];

	Word_t emp;
	emp = mesh->peek_midpoint(vtcs[0], vtcs[1]);

	// edges
	sub_eid[0] = mesh->get_edge_id(vtcs[0], emp);
	sub_ei[0] = sei = new EdgeInfo();
	MEM_CHECK(sei);
	sei->part = get_lower_part(ei->part);
	ei_data[sub_eid[0]] = sei;

	sub_eid[1] = mesh->get_edge_id(emp, vtcs[1]);
	sub_ei[1] = sei = new EdgeInfo();
	MEM_CHECK(sei);
	sei->part = get_higher_part(ei->part);
	ei_data[sub_eid[1]] = sei;

	// --- ////

	// vertex
	calc_vertex_vertex_ced(vtcs[0], vtcs[1]);

	// edge by edge
	Word_t cng_edge_id = mesh->get_edge_id(elem, iedge);
	int cng_edge_ori = elem->get_edge_orientation(iedge);

	calc_vertex_edge_ced(emp, cng_edge_id, cng_edge_ori, ei->part);

	calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp), cng_edge_id, cng_edge_ori, 1, sub_ei[0]->part);
	calc_edge_edge_ced(mesh->get_edge_id(emp, vtcs[1]), cng_edge_id, cng_edge_ori, 2, sub_ei[1]->part);
}

void Space::uc_face(Word_t eid, int iface) {
	Element *elem = mesh->elements[eid];

	Word_t fid = mesh->get_facet_id(elem, iface);
	assert(fid != INVALID_IDX);
	if (!fi_data.exists(fid)) return;

//	printf("*");

	FaceInfo *fi = fi_data[fid];
	assert(fi != NULL);

	Facet *facet = mesh->facets[fid];
	assert(facet != NULL);

//	printf("element #%d, face %d, ori = %d\n", eid, iface, elem->get_face_orientation(iface));
//	printf("element #%d, face %d (reft = %d)\n", eid, iface, facet->ref_mask);

	// vertices
	int nv = elem->get_face_num_of_vertices(iface);
	Word_t vtcs[nv];
	elem->get_face_vertices(iface, vtcs);
	// face edges
	const int *face_edge = elem->get_face_edges(iface);

	assert(fi->elem_id != INVALID_IDX);
	Element *big_elem = mesh->elements[fi->elem_id];

	int cng_face_id = mesh->get_facet_id(big_elem, fi->face);
	int cng_face_ori = big_elem->get_face_orientation(iface);
	const int *cng_face_edge = big_elem->get_face_edges(iface);

	PRINTF(" - big_elem = %d, cng_face_ori = %d\n", fi->elem_id, cng_face_ori);

	Element *par_elem = mesh->elements[eid];
	int par_face_ori = par_elem->get_face_orientation(iface);
	PRINTF(" - par_elem = %d, par_face_ori = %d\n", eid, par_face_ori);

	int h_edge_part = face_to_edge_part(fi->h_part);
	int v_edge_part = face_to_edge_part(fi->v_part);
//	double h_edge_coef = get_edge_coef(h_edge_part);
//	double v_edge_coef = get_edge_coef(v_edge_part);


	Word_t emp[4], fmp;		// four edge mid-points, one face mid-point
	Word_t cng_edge_id;		// constraining edge id
	int cng_edge_ori;		// orientation of the constraining edge
	Word_t edge_id[4];		// ID of two edges (left-right | upper-lower)
	int edge_ori[4];		// orientation of two edges

	Word_t seid, sfid, sub_fid[4], mid_edge_id;
	FaceInfo *sfi, *sub_fi[4];
//	EdgeInfo *ei, *sei;
	int part0, part1, pt;

	int part_ori;			// orientation of edge/face constraint
	double bc_coef;

	switch (facet->ref_mask) {
		case REFT_QUAD_HORZ:
			PRINTF("HORZ\n");

			emp[0] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[1] = mesh->peek_midpoint(vtcs[3], vtcs[0]);

			// faces
			sub_fid[0] = mesh->get_facet_id(4, vtcs[0], vtcs[1], emp[0], emp[1]);
			sub_fi[0] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = fi->h_part;
			sfi->v_part = get_lower_part(fi->v_part);
			fi_data[sub_fid[0]] = sfi;

			sub_fid[1] = mesh->get_facet_id(4, emp[1], emp[0], vtcs[2], vtcs[3]);
			sub_fi[1] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = fi->h_part;
			sfi->v_part = get_higher_part(fi->v_part);
			fi_data[sub_fid[1]] = sfi;

			// --- ////

			// vertex
			calc_vertex_vertex_ced(vtcs[1], vtcs[2]);
			calc_vertex_vertex_ced(vtcs[3], vtcs[0]);

			// edge by edge
			cng_edge_id = mesh->get_edge_id(elem, face_edge[1]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[1]);

//			printf("emp = %d, eid = %d, iface = %d\n", emp[0], eid, iface);
//			if (emp[0] == 20 && eid == 6 && iface == 0) {
//				printf("BBB\n");
//				calc_vertex_edge_ced(emp[0], cng_edge_id, cng_edge_ori, fi->h_part);
//			}
//			else
//				calc_vertex_edge_ced(emp[0], cng_edge_id, cng_edge_ori, fi->v_part);
			calc_vertex_edge_ced(emp[0], cng_edge_id, cng_edge_ori, fi->v_part);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[1], emp[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[0], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->v_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[3]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[3]);

			calc_vertex_edge_ced(emp[1], cng_edge_id, cng_edge_ori, fi->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp[1]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[1], vtcs[3]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->v_part);

			// mid edge
			mid_edge_id = mesh->get_edge_id(emp[0], emp[1]);
			edge_id[0] = mesh->get_edge_id(vtcs[0], vtcs[1]);
			edge_id[1] = mesh->get_edge_id(vtcs[2], vtcs[3]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[0]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[2]);

			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 0, fi->h_part);

			// edge by face
			// face orientations 4-7 have switched vertical and horizontal orientation (refer to the Pavel Solin's grey book)
//			part_ori = cng_face_ori < 4 ? PART_ORI_VERT : PART_ORI_HORZ;
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, v_edge_part);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, PART_ORI_VERT, fi->h_part, v_edge_part);
//			cng_face_ori = par_face_ori;

//			part_ori = cng_face_ori < 4 ? PART_ORI_HORZ : PART_ORI_VERT;
			part_ori = PART_ORI_HORZ;
//			printf(" - %d*, %d\n", fi->h_part, fi->v_part);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, v_edge_part);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, fi->v_part);

			// face by face
			calc_face_face_ced(sub_fid[0], cng_face_id, cng_face_ori, sub_fi[0]->h_part, sub_fi[0]->v_part);
			calc_face_face_ced(sub_fid[1], cng_face_id, cng_face_ori, sub_fi[1]->h_part, sub_fi[1]->v_part);
			break;

		case REFT_QUAD_VERT:
			PRINTF("VERT\n");

			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[2], vtcs[3]);

			// faces
			sub_fid[0] = mesh->get_facet_id(4, vtcs[0], emp[0], emp[1], vtcs[3]);
			sub_fi[0] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_lower_part(fi->h_part);
			sfi->v_part = fi->v_part;
			fi_data[sub_fid[0]] = sfi;

			sub_fid[1] = mesh->get_facet_id(4, emp[0], vtcs[1], vtcs[2], emp[1]);
			sub_fi[1] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_higher_part(fi->h_part);
			sfi->v_part = fi->v_part;
			fi_data[sub_fid[1]] = sfi;

			// --- ////

			// vertex
			calc_vertex_vertex_ced(vtcs[0], vtcs[1]);
			calc_vertex_vertex_ced(vtcs[2], vtcs[3]);

			// edge by edge
			cng_edge_id = mesh->get_edge_id(elem, face_edge[0]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[0]);

			calc_vertex_edge_ced(emp[0], cng_edge_id, cng_edge_ori, fi->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[0], vtcs[1]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->h_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[2]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[2]);

			calc_vertex_edge_ced(emp[1], cng_edge_id, cng_edge_ori, fi->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[3], emp[1]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[1], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->h_part);

			// mid edge
			mid_edge_id = mesh->get_edge_id(emp[0], emp[1]);
			edge_id[0] = mesh->get_edge_id(elem, face_edge[3]);
			edge_id[1] = mesh->get_edge_id(elem, face_edge[1]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[3]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[1]);

			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 0, fi->v_part);

			// edge by face
			// face orientations 4-7 have switched vertical and horizontal orientation (refer to the Pavel Solin's grey book)
//			part_ori = cng_face_ori < 4 ? PART_ORI_VERT : PART_ORI_HORZ;
			part_ori = PART_ORI_VERT;
//			printf(" - %d, %d*\n", fi->h_part, fi->v_part);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->v_part, h_edge_part);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, fi->v_part);

			// face by face
			calc_face_face_ced(sub_fid[0], cng_face_id, cng_face_ori, sub_fi[0]->h_part, sub_fi[0]->v_part);
			calc_face_face_ced(sub_fid[1], cng_face_id, cng_face_ori, sub_fi[1]->h_part, sub_fi[1]->v_part);
			break;

		case REFT_QUAD_BOTH:
			PRINTF("BOTH\n");

			// TODO
			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[2] = mesh->peek_midpoint(vtcs[2], vtcs[3]);
			emp[3] = mesh->peek_midpoint(vtcs[3], vtcs[0]);
			fmp = mesh->peek_midpoint(emp[0], emp[2]);

			// faces
			sub_fid[0] = mesh->get_facet_id(4, vtcs[0], emp[0], fmp, emp[3]);
			sub_fi[0] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_lower_part(fi->h_part);
			sfi->v_part = get_lower_part(fi->v_part);
			fi_data[sub_fid[0]] = sfi;

			sub_fid[1] = mesh->get_facet_id(4, emp[0], vtcs[1], emp[1], fmp);
			sub_fi[1] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_higher_part(fi->h_part);
			sfi->v_part = get_lower_part(fi->v_part);
			fi_data[sub_fid[1]] = sfi;

			sub_fid[2] = mesh->get_facet_id(4, fmp, emp[1], vtcs[2], emp[2]);
			sub_fi[2] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_higher_part(fi->h_part);
			sfi->v_part = get_higher_part(fi->v_part);
			fi_data[sub_fid[2]] = sfi;

			sub_fid[3] = mesh->get_facet_id(4, emp[3], fmp, emp[2], vtcs[3]);
			sub_fi[3] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_lower_part(fi->h_part);
			sfi->v_part = get_higher_part(fi->v_part);
			fi_data[sub_fid[3]] = sfi;

			// --- ////

			// vertex by vertex
			calc_vertex_vertex_ced(vtcs[0], vtcs[1]);
			calc_vertex_vertex_ced(vtcs[1], vtcs[2]);
			calc_vertex_vertex_ced(vtcs[2], vtcs[3]);
			calc_vertex_vertex_ced(vtcs[3], vtcs[0]);
			calc_mid_vertex_vertex_ced(fmp, emp[0], emp[1], emp[2], emp[3]);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[0]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[0]);
			calc_mid_vertex_edge_ced(emp[0], fmp, cng_edge_id, cng_edge_ori, fi->h_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[1]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[1]);
			calc_mid_vertex_edge_ced(emp[1], fmp, cng_edge_id, cng_edge_ori, fi->v_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[2]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[2]);
			calc_mid_vertex_edge_ced(emp[2], fmp, cng_edge_id, cng_edge_ori, fi->h_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[3]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[3]);
			calc_mid_vertex_edge_ced(emp[3], fmp, cng_edge_id, cng_edge_ori, fi->v_part);

			// vertex by face
			calc_vertex_face_ced(fmp, cng_face_id, cng_face_ori, iface, fi->h_part, fi->v_part);

			// edge by edge
			cng_edge_id = mesh->get_edge_id(elem, face_edge[0]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[0]);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[0], vtcs[1]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->h_part);

			// egde by edge (right)
			cng_edge_id = mesh->get_edge_id(elem, face_edge[1]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[1]);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[1], emp[1]), cng_edge_id, cng_edge_ori, 1, sub_fi[1]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[1], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[2]->v_part);

			// edge by edge (upper)
			cng_edge_id = mesh->get_edge_id(elem, face_edge[2]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[2]);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[3], emp[2]), cng_edge_id, cng_edge_ori, 1, sub_fi[3]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[2], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[2]->h_part);

			// edge by edge (left)
			cng_edge_id = mesh->get_edge_id(elem, face_edge[3]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[3]);

			calc_edge_edge_ced(mesh->get_edge_id(emp[3], vtcs[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[3], emp[3]), cng_edge_id, cng_edge_ori, 2, sub_fi[3]->v_part);

			// mid egdes (vert)
			edge_id[0] = mesh->get_edge_id(elem, face_edge[3]);
			edge_id[1] = mesh->get_edge_id(elem, face_edge[1]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[3]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[1]);

			mid_edge_id = mesh->get_edge_id(emp[0], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 1, sub_fi[0]->v_part);

			mid_edge_id = mesh->get_edge_id(emp[2], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 2, sub_fi[3]->v_part);

			// mid edges (horz)
			edge_id[0] = mesh->get_edge_id(elem, face_edge[0]);
			edge_id[1] = mesh->get_edge_id(elem, face_edge[2]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[0]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[2]);
			//

			mid_edge_id = mesh->get_edge_id(emp[3], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 1, sub_fi[0]->h_part);
			//
			mid_edge_id = mesh->get_edge_id(emp[1], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 2, sub_fi[1]->h_part);

			// edge by face
//			part_ori = cng_face_ori < 4 ? PART_ORI_VERT : PART_ORI_HORZ;
			part_ori = PART_ORI_VERT;

			edge_id[0] = mesh->get_edge_id(vtcs[0], emp[3]);
			edge_id[1] = mesh->get_edge_id(vtcs[1], emp[1]);
			mid_edge_id = mesh->get_edge_id(emp[0], fmp);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[0]->v_part, h_edge_part);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, sub_fi[0]->v_part);

			edge_id[0] = mesh->get_edge_id(emp[3], vtcs[3]);
			edge_id[1] = mesh->get_edge_id(emp[1], vtcs[2]);
			mid_edge_id = mesh->get_edge_id(emp[2], fmp);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[3]->v_part, h_edge_part);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, sub_fi[3]->v_part);

//			part_ori = cng_face_ori < 4 ? PART_ORI_HORZ : PART_ORI_VERT;
			part_ori = PART_ORI_HORZ;

			edge_id[0] = mesh->get_edge_id(vtcs[0], emp[0]);
			edge_id[1] = mesh->get_edge_id(vtcs[3], emp[2]);
			mid_edge_id = mesh->get_edge_id(emp[3], fmp);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[0]->h_part, v_edge_part);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[0]->h_part, fi->v_part);

			edge_id[0] = mesh->get_edge_id(emp[0], vtcs[1]);
			edge_id[1] = mesh->get_edge_id(emp[2], vtcs[2]);
			mid_edge_id = mesh->get_edge_id(fmp, emp[1]);
//			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[1]->h_part, v_edge_part);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[1]->h_part, fi->v_part);


			// face by face
			calc_face_face_ced(sub_fid[0], cng_face_id, cng_face_ori, sub_fi[0]->h_part, sub_fi[0]->v_part);
			calc_face_face_ced(sub_fid[1], cng_face_id, cng_face_ori, sub_fi[1]->h_part, sub_fi[1]->v_part);
			calc_face_face_ced(sub_fid[2], cng_face_id, cng_face_ori, sub_fi[2]->h_part, sub_fi[2]->v_part);
			calc_face_face_ced(sub_fid[3], cng_face_id, cng_face_ori, sub_fi[3]->h_part, sub_fi[3]->v_part);

			break;
	}
}

void Space::uc_element(Word_t idx) {
	if (idx == INVALID_IDX) return;

	Element *e = mesh->elements[idx];
//	printf("uc element #%d (active = %d)\n", idx, e->active);

//	// update constraints on element edges
	for (int iedge = 0; iedge < e->get_num_of_edges(); iedge++) {
//		uc_edge(idx, iedge);
		Word_t edge_id = mesh->get_edge_id(e, iedge);
		Edge edg = mesh->edges[edge_id];

//		printf("iedge = %d, (active = %d, bnd = %d), edge_id = %d (%p)", iedge, edg.is_active(), edg.bnd, edge_id, en_data[edge_id]->bc_proj);
		if (edg.is_active())
			calc_edge_boundary_projection(e, iedge);
//		printf("\n");
	}
//	printf("\n");

	for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
		Word_t fid = mesh->get_facet_id(e, iface);
		Facet *facet = mesh->facets[fid];

//		printf("iface = %d, type = %d, (l, r) = (%d, %d), fid = %d", iface, facet->type, facet->lactive, facet->ractive, fid);
		if (facet->ractive && facet->lactive && mesh->facets[fid]->type == Facet::OUTER)
			calc_face_boundary_projection(e, iface);

		if (facet->ced(idx, iface)) {
			if (!fi_data.exists(fid)) {
				FaceInfo *fi = NULL;
				switch (facet->mode) {
					case MODE_QUAD:
//						printf("   - creating new for face: %d (iface = %d)\n", fid, iface);
						fi_data[fid] = new FaceInfo(MODE_QUAD, idx, iface);
						MEM_CHECK(fi_data[fid]);
						break;

					case MODE_TRIANGLE:
						EXIT(ERR_NOT_IMPLEMENTED);
						break;

					default:
						EXIT(ERR_UNKNOWN_MODE);
						break;
				}
			}

			// face
//			printf("iface = %d\n", iface);
			uc_face(idx, iface);
//			printf(" ced");
		}

//		printf("\n");
	}
//	printf("\n");

#if 0
	Element *e = mesh->elements[idx];
	printf("uc element #%d (active = %d)\n", idx, e->active);
/*	for (int i = 0; i < e->get_num_of_faces(); i++) {
//		printf(" * face %2d: ori = %d, ", i, e->get_face_orientation(i));

		Word_t fid = mesh->get_facet_id(e, i);
		if (!fi_data.exists(fid)) {
//			printf("NULL\n");
		}
		else {
			FaceInfo *fi = fi_data[fid];
//			printf("elem = %d, face = %d, part = (%d, %d)\n", fi->elem_id, fi->face, fi->h_part, fi->v_part);
		}
	}
//	printf(" ---\n");
*/
	if (e->active) {
		// FIXME: update edge-, face- nodes

		// update boundary projections (only edge and face, vertex ones are already done)
		for (int iedge = 0; iedge < e->get_num_of_edges(); iedge++) {
			Word_t eid = mesh->get_edge_id(e, iedge);
			if (en_data[eid]->bnd)
				calc_edge_boundary_projection(e, iedge);
		}

		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			if (mesh->facets[fid]->type == Facet::OUTER)
				calc_face_boundary_projection(e, iface);
		}
	}
	else {
		// update boundary projections (only edge and face, vertex ones are already done)
		for (int iedge = 0; iedge < e->get_num_of_edges(); iedge++) {
			Word_t eid = mesh->get_edge_id(e, iedge);
			if (en_data[eid]->bnd)
				calc_edge_boundary_projection(e, iedge);
		}

		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			if (mesh->facets[fid]->type == Facet::OUTER)
				calc_face_boundary_projection(e, iface);
		}


		// create EdgeInfo, FaceInfo where we don't have them yet
		for (int iedge = 0; iedge < e->get_num_of_edges(); iedge++) {
			Word_t eid = mesh->get_edge_id(e, iedge);

			if (en_data[eid]->ced) {
				if (!ei_data.exists(eid)) {
					printf("   - creating new for edge: %d (iface = %d)\n", eid, iedge);
					ei_data[eid] = new EdgeInfo();
					MEM_CHECK(ei_data[eid]);
				}
			}
		}

//		printf("Element #%d\n", e->id);
		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			Facet *facet = mesh->facets[fid];

			if (facet->ced(idx, iface)) {
				// face
				if (!fi_data.exists(fid)) {
					FaceInfo *fi = NULL;
					switch (facet->mode) {
						case MODE_QUAD:
							printf("   - creating new for face: %d (iface = %d)\n", fid, iface);
							fi_data[fid] = new FaceInfo(MODE_QUAD, idx, iface);
							MEM_CHECK(fi_data[fid]);
							break;

						case MODE_TRIANGLE:
							EXIT(ERR_NOT_IMPLEMENTED);
							break;

						default:
							EXIT(ERR_UNKNOWN_MODE);
							break;
					}
				}
			}
		}

		// update constraints on element edges
		for (int iedge = 0; iedge < e->get_num_of_edges(); iedge++) {
			uc_edge(idx, iedge);
			printf("iedge = %d, ", iedge);
		}
		printf("\n");
		// update constraints on element faces
		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			uc_face(idx, iface);
			printf("iface = %d, ", iface);
		}
		printf("\n");

		// recursively down
		switch (e->reft) {
			case REFT_HEX_X:
			case REFT_HEX_Y:
			case REFT_HEX_Z:
				uc_element(e->get_son(0));
				uc_element(e->get_son(1));
				break;

			case REFT_HEX_XY:
			case REFT_HEX_YZ:
			case REFT_HEX_XZ:
				uc_element(e->get_son(0));
				uc_element(e->get_son(1));
				uc_element(e->get_son(2));
				uc_element(e->get_son(3));
				break;

			case REFT_HEX_XYZ:
				uc_element(e->get_son(0));
				uc_element(e->get_son(1));
				uc_element(e->get_son(2));
				uc_element(e->get_son(3));
				uc_element(e->get_son(4));
				uc_element(e->get_son(5));
				uc_element(e->get_son(6));
				uc_element(e->get_son(7));
				break;
		}
	}
#endif
}

//

int Space::assign_dofs(int first_dof, int stride) {
	this->first_dof = next_dof = first_dof;
	this->stride = stride;

	// free data
	Word_t i;
	FOR_ALL_VERTEX_NODES(i)
		delete vn_data[i];
	vn_data.remove_all();

	FOR_ALL_EDGE_NODES(i)
		delete en_data[i];
	en_data.remove_all();

	FOR_ALL_FACE_NODES(i)
		delete fn_data[i];
	fn_data.remove_all();

	for (Word_t i = fi_data.first(); i != INVALID_IDX; i = fi_data.next(i))
		delete fi_data[i];
	fi_data.remove_all();

	for (Word_t i = ei_data.first(); i != INVALID_IDX; i = ei_data.next(i))
		delete ei_data[i];
	ei_data.remove_all();

	// find constraints
	Word_t idx;
	FOR_ALL_BASE_ELEMENTS(idx, mesh)
		fc_element(idx);

	enforce_minimum_rule();
	set_bc_information();

	assign_dofs_internal();
//	free_extra_data();
//	calc_boundary_projections();
	update_constraints();

	return get_dof_count();
}

void Space::update_constraints() {
	Word_t idx;

#if 1
	FOR_ALL_ACTIVE_ELEMENTS(elm_idx, mesh) {
		Element *e = mesh->elements[elm_idx];
		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::OUTER) {
				// mark the vertices on the boundary
				const int *vtx = e->get_face_vertices(iface);
				for (int iv = 0; iv < e->get_face_num_of_vertices(iface); iv++) {
					calc_vertex_boundary_projection(e, vtx[iv]);
//					vn_data[e->get_vertex(vtx[iv])]->bnd = 1;
				}

				// mark the edges on the boundary
				const int *edge = e->get_face_edges(iface);
				for (int ie = 0; ie < e->get_face_num_of_edges(iface); ie++) {
					Word_t edge_id = mesh->get_edge_id(e, edge[ie]);
//					en_data[edge_id]->bnd = 1;
					if (mesh->edges[edge_id].bnd == 0) {
						printf("edge #%d should be boundary edge.\n", edge_id);
					}
					assert(mesh->edges[edge_id].bnd == 1);
//					calc_edge_boundary_projection(e, edge[ie]);
				}

/*				for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
					Word_t fid = mesh->get_facet_id(e, iface);
					Facet *facet = mesh->facets[fid];
					if (facet->type == Facet::OUTER) {
						const int *edge = e->get_face_edges(iface);
						for (int ie = 0; ie < e->get_face_num_of_edges(iface); ie++)
							calc_edge_boundary_projection(e, edge[ie]);

						calc_face_boundary_projection(e, iface);
					}
				}
*/
			}
		}
	}
#endif

	// TODO: free the items in the arrays
	for (Word_t i = fi_data.first(); i != INVALID_IDX; i = fi_data.next(i))
		delete fi_data[i];
	fi_data.remove_all();

	for (Word_t i = ei_data.first(); i != INVALID_IDX; i = ei_data.next(i))
		delete ei_data[i];
	ei_data.remove_all();

//	FOR_ALL_BASE_ELEMENTS(idx, mesh)
//		uc_element(idx);

	FOR_ALL_ELEMENTS(idx, mesh)
		uc_element(idx);
}

//// BC stuff /////////////////////////////////////////////////////////////////////////////////////

static EBCType default_bc_type(int marker) {
	return BC_ESSENTIAL;
}

static scalar default_bc_value_by_coord(int marker, double x, double y, double z, int comp) {
	return 0.0;
}

void Space::set_bc_types(EBCType(*bc_type_callback)(int)) {
	if (bc_type_callback == NULL) bc_type_callback = default_bc_type;
	this->bc_type_callback = bc_type_callback;
}

void Space::set_bc_values(scalar(*bc_value_callback_by_coord)(int, double, double, double, int)) {
	if (bc_value_callback_by_coord == NULL) bc_value_callback_by_coord = default_bc_value_by_coord;
	this->bc_value_callback_by_coord = bc_value_callback_by_coord;
}

void Space::copy_callbacks(const Space *space) {
	bc_type_callback = space->bc_type_callback;
	bc_value_callback_by_coord = space->bc_value_callback_by_coord;
}

void Space::calc_boundary_projections() {
/*	Word_t elm_idx;
	FOR_ALL_ACTIVE_ELEMENTS(elm_idx, mesh) {
		Element *e = mesh->elements[elm_idx];
		for (int vtx = 0; vtx < e->get_num_of_vertices(); vtx++)
			calc_vertex_boundary_projection(e, vtx);
		for (int edge = 0; edge < e->get_num_of_edges(); edge++)
			calc_edge_boundary_projection(e, edge);
		for (int face = 0; face < e->get_num_of_faces(); face++)
			calc_face_boundary_projection(e, face);
	}
*/
/*	Element *e = mesh->elements[2];
	FaceData *fd = fn_data[13];
	assert(fd != NULL);
	fd->dof = -1;
	fd->n = 1;
	fd->bc_type = BC_ESSENTIAL;
	fd->bc_proj = NULL;
	fd->order = 52;
	calc_face_boundary_projection(e, 5);
*/
	FOR_ALL_ACTIVE_ELEMENTS(elm_idx, mesh) {
		Element *e = mesh->elements[elm_idx];
//		printf("elem #%d: \n");
		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			Facet *facet = mesh->facets[fid];
//			printf("fid = %d: ", fid);
			if (facet->type == Facet::OUTER) {
//				printf("vtx =");
				const int *vtx = e->get_face_vertices(iface);
				for (int iv = 0; iv < e->get_face_num_of_vertices(iface); iv++) {
					calc_vertex_boundary_projection(e, vtx[iv]);
//					vn_data[e->get_vertex(vtx[iv])]->bnd = 1;		// mark the vertices on the boundary
//					printf(" %d", vtx[iv]);
//					printf(" | ");
				}

//				const int *edge = e->get_face_edges(iface);
//				for (int ie = 0; ie < e->get_face_num_of_edges(iface); ie++) {
//					if (en_data[mesh->get_edge_id(e, edge[ie])]->bnd)
//						calc_edge_boundary_projection(e, edge[ie]);
//				}
				const int *edge = e->get_face_edges(iface);
//				printf("eid =");
				for (int ie = 0; ie < e->get_face_num_of_edges(iface); ie++) {
//					printf(" %d", mesh->get_edge_id(e, edge[ie]));
//					if (en_data[mesh->get_edge_id(e, edge[ie])]->bnd)
						calc_edge_boundary_projection(e, edge[ie]);
//					printf(" | ");
				}

//				if (facet->type == Facet::OUTER) {
					calc_face_boundary_projection(e, iface);
//				}
//					printf("\n");

/*				const int *edge = e->get_face_edges(iface);
				for (int ie = 0; ie < e->get_face_num_of_edges(iface); ie++) {
					if (!en_data[edge[ie]]->ced)
						calc_edge_boundary_projection(e, edge[ie]);
				}

				if (!fn_data[fid]->ced)
					calc_face_boundary_projection(e, iface);
*/
			}
		}
//		printf("\n");
	}

	/*
	FOR_ALL_BASE_ELEMENTS(eid, mesh) {
		Element *e = mesh->elements[eid];
		for (int iface = 0; iface < e->get_num_of_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::OUTER) {
				const int *vtx = e->get_face_vertices(iface);
				for (int iv = 0; iv < e->get_face_num_of_vertices(iface); iv++) {
					calc_vertex_boundary_projection(e, vtx[iv]);
					vn_data[e->get_vertex(vtx[iv])]->bnd = 1;		// mark the vertices on the boundary
				}

				const int *edge = e->get_face_edges(iface);
				for (int ie = 0; ie < e->get_face_num_of_edges(iface); ie++)
					calc_edge_boundary_projection(e, edge[ie]);

				calc_face_boundary_projection(e, iface);
			}
		}
	}
	*/
}

