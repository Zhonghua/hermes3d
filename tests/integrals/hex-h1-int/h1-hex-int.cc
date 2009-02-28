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

/*
 * h1-hex-int.cc
 *
 * Testing H1-integrals
 *
 */

//
// NOTES: everything should be passed into this test via command line
// - mesh
// - indices of shapes functions (two)
// - exact values of integrals (seven)
//

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

// first two Lobatto shape functions
#define l0(x) ((1.0 - (x)) * 0.5)
#define l1(x) ((1.0 + (x)) * 0.5)

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 2) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	// Load mesh
	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

	// init
	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	RefMap ru(&mesh), rv(&mesh);
	PrecalcShapeset fu(pss), fv(pss);

	// FIXME: meshes can have more than 1 element
	Element *elem = mesh.elements[0];
	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);
	pss.set_quad(quad);
	fu.set_quad(quad);
	fv.set_quad(quad);

	fu.set_active_element(elem);
	fv.set_active_element(elem);

    ru.set_active_element(elem);
    rv.set_active_element(elem);

    // pick U and V function
	fu.set_active_shape(shapeset.get_vertex_index(0));
	fv.set_active_shape(shapeset.get_vertex_index(0));

	double result;

	//
	// volume integrals
	//

	result = int_u(&fu, &ru);
	printf("int_u = %lf\n", result);

	result = int_u_v(&fu, &fv, &ru, &rv);
	printf("int_u_v = %lf\n", result);

	// TODO: int_F_v

	result = int_grad_u_grad_v(&fu, &fv, &ru, &rv);
	printf("int_grad_u_grad_v = %lf\n", result);

	//
	// surface integrals
	//

	// surf_int_v
	result = 0.0;
	for (Word_t fid = mesh.facets.first(); fid != INVALID_IDX; fid = mesh.facets.next(fid)) {
		Facet *facet = mesh.facets.get(fid);
		assert(facet != NULL);

		if (facet->type == Facet::OUTER) {
			// FIXME: for meshes with more than 1 element, we need to set the active element for fu, fv, ru, rv
			FacePos fp;
			fp.face = facet->left_face_num;
			result += surf_int_v(&fv, &rv, &fp);
		}

	}
	printf("surf_int_v = %lf\n", result);

	// TODO: surf_int_G_v
	// TODO: surf_int_u_v

	TRACE_END;

	return res;
}

