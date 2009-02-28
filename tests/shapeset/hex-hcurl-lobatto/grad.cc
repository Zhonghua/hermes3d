// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 David Andrs <dandrs@unr.edu>
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

/*
 * grad.cc
 *
 * testing correctness of a gradients
 *
 */

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

void hcurl_int_vol(RealFunction *fu, double3 result) {
	// TODO: PUT YOUR CODE HERE
}

void hcurl_int_surf(RealFunction *fu, double3 result) {
	// TODO: PUT YOUR CODE HERE
}

bool test_grad(int fn_idx, Shapeset *shapeset) {
	PrecalcShapeset pss_u(shapeset);
	pss_u.set_quad(get_quadrature(MODE));

	pss_u.set_active_shape(fn_idx);

//	printf("  fn #%d\n", fn_idx);

	printf(".");
	fflush(stdout);			// prevent caching of output (to see that it did not freeze)

	double3 vol_val = { 0.0, 0.0, 0.0 };
	double3 surf_val = { 0.0, 0.0, 0.0 };

	hcurl_int_vol(&pss_u, vol_val);
	hcurl_int_surf(&pss_u, surf_val);

	if (fabs(vol_val[0] - surf_val[0]) > EPS || fabs(vol_val[1] - surf_val[1]) > EPS || fabs(vol_val[2] - surf_val[2]) > EPS) {
		printf("\n");
		ERROR("Gradient values for fn #%d do not match", fn_idx);
		return false;
	}

	return true;
}

bool test_gradients(Shapeset *shapeset) {
	printf("IV. gradients\n");

	// vertex fns
	printf("* Vertex functions\n");
	for (int i = 0; i < Hex::NUM_VERTICES; i++) {
		int fn_idx = shapeset->get_vertex_index(i);
		if (!test_grad(fn_idx, shapeset))
			return false;
	}

	// edge fns
	printf("\n* Edge functions\n");
	for (int i = 0; i < Hex::NUM_EDGES; i++) {
		int order = MAX_ELEMENT_ORDER;
		int *edge_idx = shapeset->get_edge_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_edge_fns(order); j++) {
			if (!test_grad(edge_idx[j], shapeset))
				return false;
		}
	}

	// face fns
	printf("\n* Face functions\n");
	for (int i = 0; i < Hex::NUM_FACES; i++) {
		int order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
		int *face_idx = shapeset->get_face_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_face_fns(order); j++) {
			if (!test_grad(face_idx[j], shapeset))
				return false;
		}
	}

	// bubble
	printf("\n* Bubble functions\n");
	int order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
	int *bubble_idx = shapeset->get_bubble_indices(order);
	for (int j = 0; j < shapeset->get_num_bubble_fns(order); j++) {
		if (!test_grad(bubble_idx[j], shapeset))
			return false;
	}

	return true;
}

