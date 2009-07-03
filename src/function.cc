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
#include "function.h"
#include <common/callstack.h>


// the order of items must match values of EValueType
template<typename TYPE>
int Function<TYPE>::idx2mask[][COMPONENTS] = {
	{ FN_VAL_0, FN_VAL_1, FN_VAL_2 },
	{ FN_DX_0,  FN_DX_1,  FN_DX_2  },
	{ FN_DY_0,  FN_DY_1,  FN_DY_2  },
	{ FN_DZ_0,  FN_DZ_1,  FN_DZ_2  },
	{ FN_DXX_0, FN_DXX_1, FN_DXX_2 },
	{ FN_DYY_0, FN_DYY_1, FN_DYY_2 },
	{ FN_DZZ_0, FN_DZZ_1, FN_DZZ_2 },
	{ FN_DXY_0, FN_DXY_1, FN_DXY_2 },
	{ FN_DYZ_0, FN_DYZ_1, FN_DYZ_2 },
	{ FN_DXZ_0, FN_DXZ_1, FN_DXZ_2 }
};


template<typename TYPE>
Function<TYPE>::Function() {
	_F_
	order = 0;
	max_mem = total_mem = 0;

	nodes = NULL;
	cur_node = NULL;
	sub_tables = NULL;

	memset(quads, 0, sizeof(quads));
	cur_quad = 0;

	sub_idx = 0;
}


template<typename TYPE>
Function<TYPE>::~Function() {
	_F_
}

template<typename TYPE>
void Function<TYPE>::set_quad(Quad3D *quad) {
	_F_
	int i;

	// check to see if we already have the quadrature
	for (i = 0; i < QUAD_COUNT; i++)
		if (quads[i] == quad) {
			cur_quad = i;
			return;
		}

	// if not, add the quadrature to a free slot
	for (i = 0; i < QUAD_COUNT; i++)
		if (quads[i] == NULL) {
			quads[i] = quad;
			cur_quad = i;
			return;
		}

	ERROR("too many quadratures.");
}

template<typename TYPE>
void Function<TYPE>::push_transform(int son) {
	_F_
	Transformable::push_transform(son);
	if (sub_tables) update_nodes_ptr(); // fixme
}


template<typename TYPE>
void Function<TYPE>::pop_transform() {
	_F_
	Transformable::pop_transform();
	if (sub_tables) update_nodes_ptr(); // fixme
}

template<typename TYPE>
typename Function<TYPE>::Node *Function<TYPE>::new_node(int mask, int num_points) {
	_F_
	// get the number of tables
	int nt = 0, m = mask;
	if (num_components < 3) m &= FN_VAL_0 | FN_DX_0 | FN_DY_0 | FN_DZ_0 | FN_DXX_0 | FN_DYY_0 | FN_DZZ_0 | FN_DXY_0 | FN_DXZ_0 | FN_DYZ_0;
	while (m) {
		nt += m & 1;
		m >>= 1;
	}

	// allocate a node including its data part, init table pointers
	int size = sizeof(Node) + sizeof(TYPE) * num_points * nt;
	Node *node = (Node *) malloc(size);
	node->mask = mask;
	node->size = size;
	memset(node->values, 0, sizeof(node->values));
	TYPE *data = node->data;
	for (int j = 0; j < num_components; j++) {
		for (int i = 0; i < VALUE_TYPES; i++)
			if (mask & idx2mask[i][j]) {
				node->values[j][i] = data;
				data += num_points;
			}
	}

	total_mem += size;
	if (max_mem < total_mem) max_mem = total_mem;

	return node;
}

template<typename TYPE>
void Function<TYPE>::free_nodes(void **nodes) {
	_F_
	// free all nodes stored in the tertiary Judy array
	unsigned long order = 0;
	void **pp = (void **) JudyLFirst(*nodes, &order, NULL);
	while (pp != NULL) {
		// free the concrete Node structure
		total_mem -= ((Node *) *pp)->size;
		::free(*pp);
		pp = JudyLNext(*nodes, &order, NULL);
	}
	JudyLFreeArray(nodes, NULL);
}

template<typename TYPE>
void Function<TYPE>::free_sub_tables(void **sub) {
	_F_
	// iterate through the specified secondary (sub_idx) Judy array
	unsigned long idx = 0;
	void **nodes = (void **) JudyLFirst(*sub, &idx, NULL);
	while (nodes != NULL) {
	    free_nodes(nodes);
		nodes = JudyLNext(*sub, &idx, NULL);
	}
	JudyLFreeArray(sub, NULL);
}

#undef CHECK_PARAMS
#undef CHECK_TABLE
