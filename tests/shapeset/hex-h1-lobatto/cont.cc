// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 David Andrs <dandrs@unr.edu>
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

//
// testing continuity of function values over faces
//

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

// Everything is tested on the following geometry
//
//
//      7             6             11
//        +-----------+-----------+
//       /|          /|          /|
//      / |         / |         / |
//     /  |      5 /  |     10 /  |
//  4 +-----------+-----------+   |
//    |   |       |   |       |   |
//    |   +-------|---+-------|---+
//    |  / 3      |  / 2      |  / 9
//    | /         | /         | /
//    |/          |/          |/
//    +-----------+-----------+
//   0            1            8
//
//  ^ z
//  |
//  | / y
//  |/
//  +---> x

// vertices
Point3D vtcs[] = {
	{ -2, -1, -1 },
	{  0, -1, -1 },
	{  0,  1, -1 },
	{ -2,  1, -1 },
	{ -2, -1,  1 },
	{  0, -1,  1 },
	{  0,  1,  1 },
	{ -2,  1,  1 },
	{  2, -1, -1 },
	{  2,  1, -1 },
	{  2, -1,  1 },
	{  2,  1,  1 },
};

// CAUTION: do not change the following unless you know what you are doing
// the order of hexes is crucial (function get_face_perm_ori depends on the ordering)
Word_t hex[2][24][8] = {
	// hex 1
	{
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },
		{  3,  0,  1,  2,  7,  4,  5,  6 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  0,  4,  5,  1,  3,  7,  6,  2 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  7,  6,  5,  4,  3,  2,  1,  0 },
		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  3,  2,  6,  7,  0,  1,  5,  4 }
	},
	// hex 2
	{
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5, 10 },
		{  2,  1,  8,  9,  6,  5, 10, 11 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{  1,  5, 10,  8,  2,  6, 11,  9 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  6, 11, 10,  5,  2,  9,  8,  1 },
		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{  2,  9, 11,  6,  1,  8, 10,  5 }
	}
};


// vertex numbers
// - indexing: [element][position of the element]
int vtx_no[2][24] = {
	{ 1, 0, 3, 2,   0, 3, 2, 1,  0, 3, 2, 1,  4, 5, 6, 7,  6, 7, 4, 5,  4, 7, 6, 5 },
	{ 0, 3, 2, 1,   4, 7, 6, 5,  1, 0, 3, 2,  0, 1, 2, 3,  7, 4, 5, 6,  7, 6, 5, 4 }
};

// edge numbers
// - indexing: [element][position of the element]
// - horizontal edge shared by elements in y-axis
int edge_no[2][24] = {
	{ 1, 0, 3, 2,   3, 2, 1, 0,  4, 7, 6, 5,  8, 9,10,11,  9,10,11, 8,  4, 7, 6, 5 },
	{ 3, 2, 1, 0,  11,10, 9, 8,  5, 4, 7, 6,  0, 1, 2, 3, 11, 8, 9,10,  7, 6, 5, 4 }
};

// - indexing: [element][position of the element]
// - vertical edge shared by elements in z-axis
int edge_v_no[2][24] = {
	{ 5, 4, 7, 6,   0, 3, 2, 1,  3, 2, 1, 0,  11, 8, 9,10,  6, 7, 4, 5,   8,11,10, 9 },
	{ 4, 7, 6, 5,   8,11,10, 9,  1, 0, 3, 2,   3, 0, 1, 2,  7, 4, 5, 6,  10, 9, 8,11 }
};

// face numbers
// - indexing: [element][position of the element]
int face_no[2][24] = {
	{ 1, 2, 0, 3,   4, 4, 4, 4,  0, 3, 1, 2,  5, 5, 5, 5,  1, 3, 0, 2,  2, 0, 3, 1 },
	{ 0, 3, 1, 2,   5, 5, 5, 5,  1, 2, 0, 3,  4, 4, 4, 4,  0, 2, 1, 3,  3, 1, 2, 0 }
};

// sets the orientation according to the edge and its orientation
void get_edge_perm_ori(int edge, int ori, double c[]) {
	switch (edge) {
		// edges in x-axis
		case 0:
		case 8:
		case 10:
		case 2:
			c[0] = ori == 0 ? 1 : -1;
			c[1] = c[2] = 1;
			break;

		// edges in y-axis
		case 1:
		case 9:
		case 11:
		case 3:
			c[0] = c[2] = 1;
			c[1] = ori == 0 ? 1 : -1;
			break;

		// edges in z-axis
		case 4:
		case 5:
		case 6:
		case 7:
			c[0] = c[1] = 1;
			c[2] = ori == 0 ? 1 : -1;
			break;
	}
}

// sets the permutation and the orientation according to the orientation of the hex
void get_face_perm_ori(int pos, int perm[], double c[]) {
	switch (pos) {
		//
		case 0:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = 1; c[2] = 1;
			break;

		case 1:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = 1; c[2] = 1;
			break;

		case 2:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = -1; c[2] = 1;
			break;

		case 3:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = -1; c[1] = 1; c[2] = 1;
			break;

		//
		case 4:
			perm[0] = 1; perm[1] = 0; perm[2] = 2;
			c[0] = 1; c[1] = 1; c[2] = 1;
			break;

		case 5:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = -1; c[2] = 1;
			break;

		case 6:
			perm[0] = 1; perm[1] = 0; perm[2] = 2;
			c[0] = -1; c[1] = -1; c[2] = 1;
			break;

		case 7:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = -1; c[1] = 1; c[2] = 1;
			break;

		//
		case 8:
			perm[0] = 0; perm[1] = 2; perm[2] = 1;
			c[0] = 1; c[1] = 1; c[2] = 1;
			break;

		case 9:
			perm[0] = 2; perm[1] = 1; perm[2] = 0;
			c[0] = 1; c[1] = 1; c[2] = 1;
			break;

		case 10:
			perm[0] = 0; perm[1] = 2; perm[2] = 1;
			c[0] = 1; c[1] = -1; c[2] = 1;
			break;

		case 11:
			perm[0] = 2; perm[1] = 1; perm[2] = 0;
			c[0] = -1; c[1] = 1; c[2] = 1;
			break;

		//
		case 12:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = 1; c[2] = 1;
			break;

		case 13:
			perm[0] = 1; perm[1] = 0; perm[2] = 2;
			c[0] = -1; c[1] = 1; c[2] = 1;
			break;

		case 14:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = -1; c[1] = -1; c[2] = 1;
			break;

		case 15:
			perm[0] = 1; perm[1] = 0; perm[2] = 2;
			c[0] = 1; c[1] = -1; c[2] = 1;
			break;

		//
		case 16:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = -1; c[2] = -1;
			break;

		case 17:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = 1; c[2] = -1;
			break;

		case 18:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = 1; c[1] = 1; c[2] = -1;
			break;

		case 19:
			perm[0] = 0; perm[1] = 1; perm[2] = 2;
			c[0] = -1; c[1] = 1; c[2] = -1;
			break;

		//
		case 20:
			perm[0] = 2; perm[1] = 1; perm[2] = 0;
			c[0] = 1; c[1] = 1; c[2] = -1;
			break;

		case 21:
			perm[0] = 0; perm[1] = 2; perm[2] = 1;
			c[0] = 1; c[1] = -1; c[2] = -1;
			break;

		case 22:
			perm[0] = 2; perm[1] = 1; perm[2] = 0;
			c[0] = -1; c[1] = 1; c[2] = -1;
			break;

		case 23:
			perm[0] = 0; perm[1] = 2; perm[2] = 1;
			c[0] = 1; c[1] = 1; c[2] = -1;
			break;
	}
}

//
// test the continuity of fn. values for vertex shape functions
// - check at the vertex, on the edge and on the face
//
// @return false 	if fn. values do not match
//
bool test_cont_values_of_vertex_fns(Mesh *mesh, int pos0, int pos1, Shapeset *shapeset) {
	Quad3D *quad = get_quadrature(MODE);

	printf("  - vtx_no = %d, %d\n", vtx_no[0][pos0], vtx_no[1][pos1]);

	int perm0[] = { 0, 1, 2 };
	int perm1[] = { 0, 1, 2 };

	double c0[] = { 1, 1, 1 };
	double c1[] = { 1, 1, 1 };

	// elements
	Element *hex[] = {
		mesh->elements[1],
		mesh->elements[2]
	};

	// vertex functions on both faces
	int vtx_fn[] = {
		shapeset->get_vertex_index(vtx_no[0][pos0]),
		shapeset->get_vertex_index(vtx_no[1][pos1]),
	};

	// check the same value at vertex
	printf("    - at vertex\n");

	const Point3D *hex_vtcs = REF_DOMAIN::get_vertices();
	const Point3D *vtx_pt[] = {
		&hex_vtcs[vtx_no[0][pos0]],			// vertex[0]
		&hex_vtcs[vtx_no[1][pos1]],			// vertex[1]
	};

	double x0 = vtx_pt[0]->x;
	double y0 = vtx_pt[0]->y;
	double z0 = vtx_pt[0]->z;

	double x1 = vtx_pt[1]->x;
	double y1 = vtx_pt[1]->y;
	double z1 = vtx_pt[1]->z;

	if (fabs(shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)) > EPS) {
		printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
			vtx_fn[0], vtx_fn[1],
			x0, y0, z0,
			x1, y1, z1,
			shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0),
			shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0),
			shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)
		);

		printf(" - NE!");
		printf("\n");
 //		ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
		return false;
	}

	// on face
	printf("    - on face\n");

	get_face_perm_ori(pos0, perm0, c0);
	get_face_perm_ori(pos1, perm1, c1);

	order2_t face_order = quad->get_face_max_order(face_no[0][pos0]);
	QuadPt3D *face_pts[] = {
		quad->get_face_points(face_no[0][pos0], face_order),
		quad->get_face_points(face_no[1][pos1], face_order)
	};

	for (int k = 0; k < quad->get_face_num_points(face_no[0][pos0], face_order); k++) {
		double x0 = c0[0] * face_pts[0][k][perm0[0]];
		double y0 = c0[1] * face_pts[0][k][perm0[1]];
		double z0 = c0[2] * face_pts[0][k][perm0[2]];

		double x1 = c1[0] * face_pts[1][k][perm1[0]];
		double y1 = c1[1] * face_pts[1][k][perm1[1]];
		double z1 = c1[2] * face_pts[1][k][perm1[2]];

		if (fabs(shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)) > EPS) {
			printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
				vtx_fn[0], vtx_fn[1],
				x0, y0, z0,
				x1, y1, z1,
				shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0),
				shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0),
				shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)
			);
			printf(" - NE!");
			printf("\n");
//			ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
			return false;
		}
	}

	// on edges
	order1_t edge_order = MAX_QUAD_ORDER;

	// on H edge
	int edge_ori[] =  {
		hex[0]->get_edge_orientation(edge_no[0][pos0]),
		hex[1]->get_edge_orientation(edge_no[1][pos1])
	};

	get_edge_perm_ori(edge_no[0][pos0], edge_ori[0], c0);
	get_edge_perm_ori(edge_no[1][pos1], edge_ori[1], c1);

	printf("    - on edge - h (%d, %d), ori = (%d, %d)\n", edge_no[0][pos0], edge_no[1][pos1], edge_ori[0], edge_ori[1]);

	QuadPt3D *edge_pts[] = {
		quad->get_edge_points(edge_no[0][pos0], edge_order),
		quad->get_edge_points(edge_no[1][pos1], edge_order)
	};

	for (int k = 0; k < quad->get_edge_num_points(edge_order); k++) {
		// position on ref. domain for edge[0]
		double x0 = c0[0] * edge_pts[0][k][0];
		double y0 = c0[1] * edge_pts[0][k][1];
		double z0 = c0[2] * edge_pts[0][k][2];

		// position on ref. domain for edge[1]
		double x1 = c1[0] * edge_pts[1][k][0];
		double y1 = c1[1] * edge_pts[1][k][1];
		double z1 = c1[2] * edge_pts[1][k][2];

		if (fabs(shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)) > EPS) {
			printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
				vtx_fn[0], vtx_fn[1],
				x0, y0, z0,
				x1, y1, z1,
				shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0),
				shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0),
				shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)
			);

			printf(" - NE!");
			printf("\n");
//			ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
			return false;
		}
	}

	// on V edge
	int edge_v_ori[] =  {
		hex[0]->get_edge_orientation(edge_v_no[0][pos0]),
		hex[1]->get_edge_orientation(edge_v_no[1][pos1])
	};

	get_edge_perm_ori(edge_v_no[0][pos0], edge_v_ori[0], c0);
	get_edge_perm_ori(edge_v_no[1][pos1], edge_v_ori[1], c1);

	printf("    - on edge - v (%d, %d), ori = (%d, %d)\n", edge_v_no[0][pos0], edge_v_no[1][pos1], edge_v_ori[0], edge_v_ori[1]);

	QuadPt3D *edge_v_pts[] = {
		quad->get_edge_points(edge_v_no[0][pos0], edge_order),
		quad->get_edge_points(edge_v_no[1][pos1], edge_order)
	};

	for (int k = 0; k < quad->get_edge_num_points(edge_order); k++) {
		// position on ref. domain for face 0
		double x0 = c0[0] * edge_v_pts[0][k][0];
		double y0 = c0[1] * edge_v_pts[0][k][1];
		double z0 = c0[2] * edge_v_pts[0][k][2];

		// position on ref. domain for face 1
		double x1 = c1[0] * edge_v_pts[1][k][0];
		double y1 = c1[1] * edge_v_pts[1][k][1];
		double z1 = c1[2] * edge_v_pts[1][k][2];

		if (fabs(shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)) > EPS) {
			printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
				vtx_fn[0], vtx_fn[1],
				x0, y0, z0,
				x1, y1, z1,
				shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0),
				shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0),
				shapeset->get_fn_value(vtx_fn[0], x0, y0, z0, 0) - shapeset->get_fn_value(vtx_fn[1], x1, y1, z1, 0)
			);

			printf(" - NE!");
			printf("\n");
//			ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
			return false;
		}
	}

	return true;
}


//
// test the continuity of fn. values for edge shape functions
// - check on the edge and on the face
//
// @return false 	if fn. values do not match
//
bool test_cont_values_of_edge_fns(Mesh *mesh, int pos0, int pos1, Shapeset *shapeset) {
	Quad3D *quad = get_quadrature(MODE);

	// functions up to the order 'order' are tested
	int order = MAX_ELEMENT_ORDER;

	printf("  - edge = %d, %d\n", edge_no[0][pos0], edge_no[1][pos1]);

	int perm0[] = { 0, 1, 2 };
	double c0[] = { 1, 1, 1 };

	int perm1[] = { 0, 1, 2 };
	double c1[] = { 1, 1, 1 };

	// elements
	Element *hex[] = {
		mesh->elements[1],
		mesh->elements[2]
	};

	// orientations of both edge functions
	int edge_ori[] =  {
		hex[0]->get_edge_orientation(edge_no[0][pos0]),
		hex[1]->get_edge_orientation(edge_no[1][pos1])
	};

	// face functions on both faces
	int *edge_fn[] = {
		shapeset->get_edge_indices(edge_no[0][pos0], edge_ori[0], order),
		shapeset->get_edge_indices(edge_no[1][pos1], edge_ori[1], order),
	};

	get_face_perm_ori(pos0, perm0, c0);
	get_face_perm_ori(pos1, perm1, c1);

	// on face
	order2_t face_order = quad->get_face_max_order(face_no[0][pos0]);
//	int face_order = MAKE_QUAD_ORDER(4, 3);
	QuadPt3D *face_pts[] = {
		quad->get_face_points(face_no[0][pos0], face_order),
		quad->get_face_points(face_no[1][pos1], face_order)
	};

	// on face
	printf("    - on face\n");
	for (int ei = 0; ei < shapeset->get_num_edge_fns(order); ei++) {
		for (int k = 0; k < quad->get_face_num_points(edge_no[0][pos0], face_order); k++) {
			// position on ref. domain for face 0
			double x0 = c0[0] * face_pts[0][k][perm0[0]];
			double y0 = c0[1] * face_pts[0][k][perm0[1]];
			double z0 = c0[2] * face_pts[0][k][perm0[2]];

			// position on ref. domain for face 1
			double x1 = c1[0] * face_pts[1][k][perm1[0]];
			double y1 = c1[1] * face_pts[1][k][perm1[1]];
			double z1 = c1[2] * face_pts[1][k][perm1[2]];

			if (fabs(shapeset->get_fn_value(edge_fn[0][ei], x0, y0, z0, 0) - shapeset->get_fn_value(edge_fn[1][ei], x1, y1, z1, 0)) > EPS) {
				printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
					edge_fn[0][ei], edge_fn[1][ei],
					x0, y0, z0,
					x1, y1, z1,
					shapeset->get_fn_value(edge_fn[0][ei], x0, y0, z0, 0),
					shapeset->get_fn_value(edge_fn[1][ei], x1, y1, z1, 0),
					shapeset->get_fn_value(edge_fn[0][ei], x0, y0, z0, 0) - shapeset->get_fn_value(edge_fn[1][ei], x1, y1, z1, 0)
				);
				printf(" - NE!");
				printf("\n");
//				ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
				return false;
			}
		}
	}

	// on edge
	printf("    - on edge\n");
	int edge_order = MAX_QUAD_ORDER;
//	int edge_order = 4;

	get_edge_perm_ori(edge_no[0][pos0], edge_ori[0], c0);
	get_edge_perm_ori(edge_no[1][pos1], edge_ori[1], c1);

	QuadPt3D *edge_pts[] = {
		quad->get_edge_points(edge_no[0][pos0], edge_order),
		quad->get_edge_points(edge_no[1][pos1], edge_order)
	};

	for (int ei = 0; ei < shapeset->get_num_edge_fns(order); ei++) {
		for (int k = 0; k < quad->get_edge_num_points(edge_order); k++) {
			// position on ref. domain for edge[0]
			double x0 = c0[0] * edge_pts[0][k][0];
			double y0 = c0[1] * edge_pts[0][k][1];
			double z0 = c0[2] * edge_pts[0][k][2];

			// position on ref. domain for edge[1]
			double x1 = c1[0] * edge_pts[1][k][0];
			double y1 = c1[1] * edge_pts[1][k][1];
			double z1 = c1[2] * edge_pts[1][k][2];

			if (fabs(shapeset->get_fn_value(edge_fn[0][ei], x0, y0, z0, 0) - shapeset->get_fn_value(edge_fn[1][ei], x1, y1, z1, 0)) > EPS) {
				printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
					edge_fn[0][ei], edge_fn[1][ei],
					x0, y0, z0,
					x1, y1, z1,
					shapeset->get_fn_value(edge_fn[0][ei], x0, y0, z0, 0),
					shapeset->get_fn_value(edge_fn[1][ei], x1, y1, z1, 0),
					shapeset->get_fn_value(edge_fn[0][ei], x0, y0, z0, 0) - shapeset->get_fn_value(edge_fn[1][ei], x1, y1, z1, 0)
				);
				printf(" - NE!");
				printf("\n");
//				ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
				return false;
			}
		}
	}

	return true;
}

//
// test the continuity of fn. values for face shape functions
// - check on the face
//
// @return false 	if fn. values do not match
//
bool test_cont_values_of_face_fns(Mesh *mesh, int pos0, int pos1, Shapeset *shapeset) {
	Quad3D *quad = get_quadrature(MODE);

	int perm0[] = { 0, 1, 2 };
	double c0[] = { 1, 1, 1 };
	int perm1[] = { 0, 1, 2 };
	double c1[] = { 1, 1, 1 };

	// functions up to the order 'order' are tested
	order2_t order(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
//	int order = MAKE_QUAD_ORDER(4, 3);

	// elements
	Element *hex[] = {
		mesh->elements[1],
		mesh->elements[2]
	};

	// face

	// orientations of both faces
	int face_ori[] =  {
		hex[0]->get_face_orientation(face_no[0][pos0]),
		hex[1]->get_face_orientation(face_no[1][pos1])
	};

	// face functions on both faces
	int *face_fn[] = {
		shapeset->get_face_indices(face_no[0][pos0], face_ori[0], order),
		shapeset->get_face_indices(face_no[1][pos1], face_ori[1], order),
	};

	get_face_perm_ori(pos0, perm0, c0);
	get_face_perm_ori(pos1, perm1, c1);

//	int face_order = MAKE_QUAD_ORDER(4, 3);
	order2_t face_order = quad->get_face_max_order(face_no[0][pos0]);
	QuadPt3D *face_pts[] = {
		quad->get_face_points(face_no[0][pos0], face_order),
		quad->get_face_points(face_no[1][pos1], face_order)
	};

	printf("  - face = %d, %d, face_ori = %d, %d, face_no = %d, %d\n",
		face_no[0][pos0], face_no[1][pos1], face_ori[0], face_ori[1], face_no[0][pos0], face_no[1][pos1]);

	for (int fi = 0; fi < shapeset->get_num_face_fns(order); fi++) {
		for (int k = 0; k < quad->get_face_num_points(face_no[0][pos0], face_order); k++) {
			// position on ref. domain for face 0
			double x0 = c0[0] * face_pts[0][k][perm0[0]];
			double y0 = c0[1] * face_pts[0][k][perm0[1]];
			double z0 = c0[2] * face_pts[0][k][perm0[2]];

			// position on ref. domain for face 1
			double x1 = c1[0] * face_pts[1][k][perm1[0]];
			double y1 = c1[1] * face_pts[1][k][perm1[1]];
			double z1 = c1[2] * face_pts[1][k][perm1[2]];

			if (fabs(shapeset->get_fn_value(face_fn[0][fi], x0, y0, z0, 0) - shapeset->get_fn_value(face_fn[1][fi], x1, y1, z1, 0)) > EPS) {
				printf("%5x, %5x: % lf, % lf, % lf == % lf, % lf, % lf | % lf, % lf (diff = % g)",
					face_fn[0][fi], face_fn[1][fi],
					x0, y0, z0,
					x1, y1, z1,
					shapeset->get_fn_value(face_fn[0][fi], x0, y0, z0, 0),
					shapeset->get_fn_value(face_fn[1][fi], x1, y1, z1, 0),
					shapeset->get_fn_value(face_fn[0][fi], x0, y0, z0, 0) - shapeset->get_fn_value(face_fn[1][fi], x1, y1, z1, 0)
				);
				printf(" - NE!");
				printf("\n");
//				printf("\n");
//				ERROR("Face fn #%d does not match face fn #%d on faces (%d, %d).", face_fn[0][fi], face_fn[1][fi], face0, face1);
				return false;
			}
		}
	}

	return true;
}

// mesh stuff

Element *add_hex(Array<Element *> &elements, int *idx) {
	Element *hex = new Hex(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5], idx[6], idx[7]);
	Word_t pos = elements.add(hex);
	hex->id = pos;

	return hex;
}

//
// main check function
//
bool test_continuity(Shapeset *shapeset) {
	printf("III. continuity\n");

	// combine all possible situations how two hexahedrons can face each other
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			// build the mesh
			Mesh mesh;
			for (unsigned int k = 0; k < countof(vtcs); k++)
				mesh.add_vertex(vtcs[k].x, vtcs[k].y, vtcs[k].z);
			mesh.add_hex(hex[0][i] + 0);
			mesh.add_hex(hex[1][j] + 0);
			mesh.ugh();

			// test
			printf("pos = %d, %d\n", i, j);

			if (!test_cont_values_of_vertex_fns(&mesh, i, j, shapeset))
				return false;

			if (!test_cont_values_of_edge_fns(&mesh, i, j, shapeset))
				return false;

			if (!test_cont_values_of_face_fns(&mesh, i, j, shapeset))
				return false;
		}
	}

	return true;
}

