//
// precalc.cc
//

#include "config.h"
#include "common.h"
#include "precalc.h"
#include "function.cc" // non-inline template members
#include <common/error.h>


PrecalcShapeset::PrecalcShapeset(Shapeset *shapeset)
               : Function<double>()
{
	this->shapeset = shapeset;
	master_pss = NULL;
	num_components = shapeset->get_num_components();
	assert(num_components == 1 || num_components == 3);
	tables = NULL;
	update_max_index();
}


PrecalcShapeset::PrecalcShapeset(PrecalcShapeset *pss)
               : Function<double>()
{
	while (pss->is_slave())
		pss = pss->master_pss;
	master_pss = pss;
	shapeset = pss->shapeset;
	num_components = pss->num_components;
	tables = NULL;
	update_max_index();
}


PrecalcShapeset::~PrecalcShapeset() {
	free();
	JudyLFreeArray(&tables, NULL);
}


void PrecalcShapeset::update_max_index() {
	max_index = shapeset->get_max_index();
}


void PrecalcShapeset::set_quad(Quad3D *quad) {
	RealFunction::set_quad(quad);

	set_active_shape(0);
//	cur_quad = ((master_pss != NULL) ? master_pss : this)->register_quad(quad);
//	max_order = quad->get_max_order();
}

void PrecalcShapeset::set_active_shape(int index) {
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

	assert(max_index >= index);

	unsigned key = ((unsigned) (max_index - index) << 3) | cur_quad;
//	unsigned key = cur_quad | (mode << 3) | ((unsigned) (max_index[mode] - index) << 4);

	void **tab = (master_pss == NULL) ? &tables : &(master_pss->tables);
	sub_tables = (void **) JudyLIns(tab, key, NULL);
	update_nodes_ptr();

	this->index = index;
	order = shapeset->get_order(index);
}


void PrecalcShapeset::set_active_element(Element *e) {
	if (e->get_mode() != shapeset->get_mode())
		EXIT(ERR_FAILURE, "Using element with incorrect shapeset.");

	Quad3D *quad = get_quad();
//	max_order = quad->get_max_order();

	element = e;
}

void PrecalcShapeset::precalculate(Qorder qord, int mask) {
	// initialization
	Quad3D *quad = get_quad();
	assert(quad != NULL);
	QuadPt3D *pt = NULL;
	int np = 0;
	switch (qord.type) {
		case QOT_ELEMENT:
			np = quad->get_num_points(qord.order);
			pt = quad->get_points(qord.order);
			break;

		case QOT_FACE:
			np = quad->get_face_num_points(qord.face, qord.order);
			pt = quad->get_face_points(qord.face, qord.order);
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
	for (int j = 0; j < num_components; j++) {
		for (int k = 0; k < VALUE_TYPES; k++) {
			if (newmask & idx2mask[k][j]) {
				if (oldmask & idx2mask[k][j])
					memcpy(node->values[j][k], cur_node->values[j][k], np * sizeof(double));
				else
					for (int i = 0; i < np; i++)
						node->values[j][k][i] = shapeset->get_value(k, index,
								ctm->m[0] * pt[i].x + ctm->t[0],
								ctm->m[1] * pt[i].y + ctm->t[1],
								ctm->m[2] * pt[i].z + ctm->t[2], j);
			}
		}
	}

	// remove the old node and attach the new one to the Judy array
	replace_cur_node(node);
}

void PrecalcShapeset::free() {
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
	assert(master_pss != NULL);
	sub_idx = master_pss->sub_idx;
	top = master_pss->top;
	stack[top] = *(master_pss->ctm);
	ctm = stack + top;
}


void PrecalcShapeset::dump_info(int quad, FILE *f) {
	unsigned long key = 0, n1 = 0, m1 = 0, n2 = 0, n3 = 0, size = 0;
	void **sub = (void **) JudyLFirst(tables, &key, NULL);
	while (sub != NULL) {
		if ((key & 7) == quad) {
			fprintf(f, "PRIMARY TABLE, index=%d\n", max_index - (key >> 3));
			unsigned long idx = 0;
			void **nodes = (void **) JudyLFirst(*sub, &idx, NULL);
			while (nodes != NULL) {
				fprintf(f, "   SUB TABLE, sub_idx=%d\n      NODES: ", idx); n2++;
				unsigned long order = 0;
				void **pp = (void **) JudyLFirst(*nodes, &order, NULL);
				while (pp != NULL) {
					fprintf(f, "%d ", order); n3++;
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

	fprintf(f, "Number of primary tables: %d (%d for all quadratures)\n"
	        "Avg. size of sub table:   %g\n"
	        "Avg. number of nodes:     %g\n"
	        "Total number of nodes:    %d\n"
	        "Total size of all nodes:  %d bytes\n",
	        n1, m1, (double) n2 / n1, (double) n3 / n2, n3, size);
}
