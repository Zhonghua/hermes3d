//
// main.cc
//
// Testing if the values calculated by PrecalcShapset class are correct.
// The purpose is to check if the stuff using Qorder works okey.
//
//
// TODO:
// - testing specified order
// - testing dx, dy, dx values
// - testing different shapesets
// - test transformations
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

void test_vertex_values(PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(VTX_QORDER());

	double *val = pss.get_fn_values();
	int np = quad->get_vertex_num_points();

	printf(" Vertex values:\n");
	for (int k = 0; k < np; k++) {
		printf(" % lf", val[k]);
	}
	printf("\n");
}

void test_edge_values(int edge, int order, PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(EDGE_QORDER(edge, order));

	double *val = pss.get_fn_values();
	int np = quad->get_edge_num_points(order);

	printf(" Edge #%d values:\n", edge);
	for (int k = 0; k < np; k++) {
		printf(" % lf", val[k]);
	}
	printf("\n");
}

void test_face_values(int face, int order, PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(FACE_QORDER(face, order));

	double *val = pss.get_fn_values();
	int np = quad->get_face_num_points(face, order);

	printf(" Face #%d values:\n", face);
	for (int k = 0; k < np; k++) {
		printf(" % lf", val[k]);
	}
	printf("\n");
}

void test_elem_values(int order, PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(ELEM_QORDER(order));

	double *val = pss.get_fn_values();
	int np = quad->get_num_points(order);

	printf(" Elem values:\n");
	for (int k = 0; k < np; k++) {
		printf(" % lf", val[k]);
	}
	printf("\n");
}

// Main ////

int main(int argc, char *argv[]) {
	int ret = ERROR_SUCCESS;

	if (argc < 2) return ERR_NOT_ENOUGH_PARAMS;

	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	int order = MAKE_HEX_ORDER(10, 10, 10);

	FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
		Element *e = mesh.elements[idx];
		RefMap rm(&mesh);

		Quad3D *quad = get_quadrature(e->get_mode());

		pss.set_quad(quad);
		pss.set_active_element(e);
	    rm.set_active_element(e);

	    int eorder = get_hex_edge_order(0, order);
	    int forder = get_hex_face_order(0, order);

	    // functions to test
	    int indices[] = {
	    	shapeset.get_vertex_index(0),
	    	shapeset.get_edge_indices(0, 0, eorder)[0],
	    	shapeset.get_face_indices(0, 0, forder)[0],
	    	shapeset.get_bubble_indices(order)[0]
	    };

	    for (int i = 0; i < countof(indices); i++) {
	    	int index = indices[i];
	    	printf("Testing function #%d\n", index);

			pss.set_active_shape(index);

			// test values
			test_vertex_values(pss, quad);
			for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
				test_edge_values(iedge, get_hex_edge_order(iedge, order), pss, quad);
			for (int iface = 0; iface < Hex::NUM_FACES; iface++)
				test_face_values(iface, get_hex_face_order(iface, order), pss, quad);
			test_elem_values(order, pss, quad);
	    }
	}



	return ret;
}
