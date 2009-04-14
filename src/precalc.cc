// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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
#include "precalc.h"
#include "function.cc" // non-inline template members
#include <common/error.h>
#include <common/callstack.h>


PrecalcShapeset::PrecalcShapeset(Shapeset *shapeset)
               : Function<double>()
{
	_F_
	this->shapeset = shapeset;
	master_pss = NULL;
	num_components = shapeset->get_num_components();
	assert(num_components == 1 || num_components == 3);
	tables = NULL;
}


PrecalcShapeset::PrecalcShapeset(PrecalcShapeset *pss)
               : Function<double>()
{
	_F_
	while (pss->is_slave())
		pss = pss->master_pss;
	master_pss = pss;
	shapeset = pss->shapeset;
	num_components = pss->num_components;
	tables = NULL;
}

PrecalcShapeset::~PrecalcShapeset() {
	_F_
	free();
	JudyLFreeArray(&tables, NULL);
}

void PrecalcShapeset::set_quad(Quad3D *quad) {
	_F_
	RealFunction::set_quad(quad);

	set_active_shape(0);
}

void PrecalcShapeset::set_active_shape(int index) {
	_F_
	// Each precalculated table is accessed and uniquely identified by the
	// following seven items:
	//
	//   - cur_quad:  quadrature table selector (0-7)
	//   - mode:      mode of the shape function (tetra, hex, prism) (no)
	//   - index:     shape function index
	//   - sub_idx:   the index of the sub-element (not fow now, this is for multi mesh)
	//   - order:     integration rule order
	//   - component: shape function component (0-1)
	//   - val/d/dd:  values, dx, dy, dz, ddx, ddy, ddz (0-4)
	//
	// The table database is implemented as a three-way chained Judy array.
	// The key to the first Judy array ('tables') is formed by cur_quad,
	// mode and index. This gives a pointer to the second Judy array, which
	// is indexed solely by sub_idx. The last Judy array is the node table,
	// understood by the base class and indexed by order. The component and
	// val/d/dd indices are used directly in the Node structure.

	unsigned key = ((unsigned) (index) << 3) | cur_quad;

	void **tab = (master_pss == NULL) ? &tables : &(master_pss->tables);
	sub_tables = (void **) JudyLIns(tab, key, NULL);
	update_nodes_ptr();

	this->index = index;
	order = shapeset->get_order(index);
}


void PrecalcShapeset::set_active_element(Element *e) {
	_F_
	if (e->get_mode() != shapeset->get_mode())
		EXIT(ERR_FAILURE, "Using element with incorrect shapeset.");

	element = e;
}

void PrecalcShapeset::precalculate(qorder_t qord, int mask) {
	_F_
	// initialization
	Quad3D *quad = get_quad();
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

	int oldmask = (cur_node != NULL) ? cur_node->mask : 0;
	int newmask = mask | oldmask;
	Node *node = new_node(newmask, np);
	MEM_CHECK(node);

	// precalculate all required tables
	for (int ic = 0; ic < num_components; ic++) {
		for (int j = 0; j < VALUE_TYPES; j++) {
			if (newmask & idx2mask[j][ic]) {
				if (oldmask & idx2mask[j][ic])
					memcpy(node->values[ic][j], cur_node->values[ic][j], np * sizeof(double));
				else
					for (int k = 0; k < np; k++)
						node->values[ic][j][k] = shapeset->get_value(j, index,
								ctm->m[0] * pt[k].x + ctm->t[0],
								ctm->m[1] * pt[k].y + ctm->t[1],
								ctm->m[2] * pt[k].z + ctm->t[2], ic);
			}
		}
	}

	// remove the old node and attach the new one to the Judy array
	replace_cur_node(node);
}

void PrecalcShapeset::free() {
	_F_
	if (master_pss != NULL) return;

	// iterate through the primary Judy array
	unsigned long key = 0;
	void **sub = (void **) JudyLFirst(tables, &key, NULL);
	while (sub != NULL) {
		// free the secondary and tertiary Judy arrays
		free_sub_tables(sub);
		JudyLDel(&tables, key, NULL);
		sub = JudyLNext(tables, &key, NULL);
	}
}

void PrecalcShapeset::set_master_transform() {
	_F_
	assert(master_pss != NULL);
	sub_idx = master_pss->sub_idx;
	top = master_pss->top;
	stack[top] = *(master_pss->ctm);
	ctm = stack + top;
}


void PrecalcShapeset::dump_info(int quad, FILE *f) {
	_F_
	unsigned long key = 0, n1 = 0, m1 = 0, n2 = 0, n3 = 0, size = 0;
	void **sub = (void **) JudyLFirst(tables, &key, NULL);
	while (sub != NULL) {
		if ((key & 7) == (unsigned) quad) {
			fprintf(f, "PRIMARY TABLE, index=%ld\n", (key >> 3));
			unsigned long idx = 0;
			void **nodes = (void **) JudyLFirst(*sub, &idx, NULL);
			while (nodes != NULL) {
				fprintf(f, "   SUB TABLE, sub_idx=%ld\n      NODES: ", idx); n2++;
				unsigned long order = 0;
				void **pp = (void **) JudyLFirst(*nodes, &order, NULL);
				while (pp != NULL) {
					fprintf(f, "%ld ", order); n3++;
					size += ((Node*) *pp)->size;
					pp = JudyLNext(*nodes, &order, NULL);
				}
				fprintf(f, "\n");
				nodes = JudyLNext(*sub, &idx, NULL);
			}
			fprintf(f, "\n\n"); n1++;
		}
		sub = JudyLNext(tables, &key, NULL); m1++;
	}

	fprintf(f, "Number of primary tables: %ld (%ld for all quadratures)\n"
	        "Avg. size of sub table:   %g\n"
	        "Avg. number of nodes:     %g\n"
	        "Total number of nodes:    %ld\n"
	        "Total size of all nodes:  %ld bytes\n",
	        n1, m1, (double) n2 / n1, (double) n3 / n2, n3, size);
}
