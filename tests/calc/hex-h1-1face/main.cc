// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2009 Pavel Kus <pavel.kus@gmail.com>
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
 * hang-nodes-continuity.cc
 *
 * usage: $0 <mesh file> <element id> <refinement id> [<element id> <refinement id>...]
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

#define BEGIN_BLOCK						{
#define END_BLOCK							}

//#define DIRICHLET
#define NEWTON

#define X2_Y2_Z2
//#define XM_YN_ZO
#define XM_YN_ZO_2


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

// mesh
Word_t hex[2][48][8] = {
	// hex 1
	{
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },
		{  3,  0,  1,  2,  7,  4,  5,  6 },

		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },
		{  7,  6,  5,  4,  3,  2,  1,  0 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  0,  4,  5,  1,  3,  7,  6,  2 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  7,  6,  5,  4,  3,  2,  1,  0 },
		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },

		{  3,  0,  1,  2,  7,  4,  5,  6 },
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  0,  4,  5,  1,  3,  7,  6,  2 }
	},
	// hex 2
	{
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5, 10 },
		{  2,  1,  8,  9,  6,  5, 10, 11 },

		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },
		{  6, 11, 10,  5,  2,  9,  8,  1 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{  1,  5, 10,  8,  2,  6, 11,  9 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  6, 11, 10,  5,  2,  9,  8,  1 },
		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },

		{  2,  1,  8,  9,  6,  5, 10, 11 },
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5,  10 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{  1,  5, 10,  8,  2,  6, 11,  9 }
	}
};

Word_t bnd[10][5] = {
	{ 0,  3,  7,  4, 1 },
	{ 8,  9, 11, 10, 2 },
	{ 0,  1,  5,  4, 3 },
	{ 1,  8, 10,  5, 3 },
	{ 3,  2,  6,  7, 4 },
	{ 2,  9, 11,  6, 4 },
	{ 0,  1,  2,  3, 5 },
	{ 1,  8,  9,  2, 5 },
	{ 4,  5,  6,  7, 6 },
	{ 5, 10, 11,  6, 6 }
};


int m = 2, n = 2, o = 2;

double fnc(double x, double y, double z) {
#ifdef XM_YN_ZO
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
#elif defined XM_YN_ZO_2
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#elif defined X2_Y2_Z2
	return x*x + y*y + z*z;
#endif
}

double dfnc(double x, double y, double z) {
#ifdef XM_YN_ZO
	double ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	double ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	double ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);

#elif defined XM_YN_ZO_2
	double ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 2 * z;
	double ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
	double ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);
#elif defined X2_Y2_Z2
	return -6.0;
#endif
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
#ifdef XM_YN_ZO
	dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o - 1) - pow(x, 3) + 4 * pow(z, 3);
#elif defined XM_YN_ZO_2
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3);
#elif defined X2_Y2_Z2
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;
#endif

	return fnc(x, y, z);
}

//

EBCType bc_types(int marker) {
#ifdef DIRICHLET
	return BC_ESSENTIAL;
#elif defined NEWTON
	return BC_NATURAL;
#endif
}

double bc_values(int marker, double x, double y, double z, int comp) {
#ifdef DIRICHLET
	return fnc(x, y, z);
#elif defined NEWTON
	switch (marker) {
#ifdef XM_YN_ZO
		case 1: return -(m * pow(x, m - 1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z) + fnc(x, y, z);
		case 2: return   m * pow(x, m - 1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z + fnc(x, y, z);
		case 3: return -(n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2)) + fnc(x, y, z);
		case 4: return   n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2) + fnc(x, y, z);
		case 5: return -(o * pow(x, m) * pow(y, n) * pow(z, o - 1) - pow(x, 3) + 4 * pow(z, 3)) + fnc(x, y, z);
		case 6: return   o * pow(x, m) * pow(y, n) * pow(z, o - 1) - pow(x, 3) + 4 * pow(z, 3) + fnc(x, y, z);
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE); return 0.0;
#elif defined XM_YN_ZO_2
		case 1: return -(m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z) + fnc(x, y, z);
		case 2: return   m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z + fnc(x, y, z);
		case 3: return -(n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2)) + fnc(x, y, z);
		case 4: return   n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2) + fnc(x, y, z);
		case 5: return -(o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3)) + fnc(x, y, z);
		case 6: return   o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3) + fnc(x, y, z);
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE); return 0.0;
#elif defined X2_Y2_Z2
		case 1: return -(2*x) + fnc(x, y, z);
		case 2: return  (2*x) + fnc(x, y, z);
		case 3: return -(2*y) + fnc(x, y, z);
		case 4: return  (2*y) + fnc(x, y, z);
		case 5: return -(2*z) + fnc(x, y, z);
		case 6: return  (2*z) + fnc(x, y, z);
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE); return 0.0;
#endif
	}
#endif
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar bilinear_form_surf(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	return surf_int_u_v(fu, fv, ru, rv, fp);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return int_F_v(dfnc, fv, rv);
}

scalar linear_form_surf(RealFunction *fv, RefMap *rv, FacePos *fp) {
	return surf_int_G_v(fv, rv, fp);
}

// helpers ////////////////////////////////////////////////////////////////////////////////////////

int parse_reft(char *str) {
	if (strcasecmp(str, "x") == 0) return REFT_HEX_X;
	else if (strcasecmp(str, "y") == 0) return REFT_HEX_Y;
	else if (strcasecmp(str, "z") == 0) return REFT_HEX_Z;
	else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return REFT_HEX_XY;
	else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return REFT_HEX_XZ;
	else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return REFT_HEX_YZ;
	else if (strcasecmp(str, "xyz") == 0) return REFT_HEX_XYZ;
	else return REFT_HEX_NONE;
}

///

const double EPS = 10e-14;


// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;


#ifdef WITH_PETSC
	PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	try {
		for (int i = 0; i < 48; i++) {
			for (int j = 0; j < 48; j++) {
//		int i = 5; {
//		int j = 0; {
				printf("Config: %d, %d ", i, j);

				Mesh mesh;

				for (Word_t k = 0; k < countof(vtcs); k++)
					mesh.add_vertex(vtcs[k].x, vtcs[k].y, vtcs[k].z);
				Word_t h1[] = {
						hex[0][i][0] + 1, hex[0][i][1] + 1, hex[0][i][2] + 1, hex[0][i][3] + 1,
						hex[0][i][4] + 1, hex[0][i][5] + 1, hex[0][i][6] + 1, hex[0][i][7] + 1 };
				mesh.add_hex(h1);
				Word_t h2[] = {
						hex[1][j][0] + 1, hex[1][j][1] + 1, hex[1][j][2] + 1, hex[1][j][3] + 1,
						hex[1][j][4] + 1, hex[1][j][5] + 1, hex[1][j][6] + 1, hex[1][j][7] + 1 };
				mesh.add_hex(h2);
				// bc
				for (Word_t k = 0; k < countof(bnd); k++) {
					Word_t facet_idxs[Quad::NUM_VERTICES] = { bnd[k][0] + 1, bnd[k][1] + 1, bnd[k][2] + 1, bnd[k][3] + 1 };
					mesh.add_quad_boundary(facet_idxs, bnd[k][4]);
				}

				mesh.ugh();
//				mesh.dump();

//				Element *hx[] = { mesh.elements[1], mesh.elements[2] };
//				printf("[%d, %d]\n", hx[0]->get_face_orientation(1), hx[1]->get_face_orientation(2));

//				Word_t fidx[4];
//				hx[1]->get_face_vertices(2, fidx);
//				printf("FI: %d, %d, %d, %d\n", fidx[0], fidx[1], fidx[2], fidx[3]);
				printf("\n");

#ifdef OUTPUT_DIR
				BEGIN_BLOCK
					// output the mesh
					const char *of_name = OUTPUT_DIR "/ref.msh";
					FILE *ofile = fopen(of_name, "w");
					if (ofile != NULL) {
						GmshOutputEngine output(ofile);
						output.out(&mesh);
						fclose(ofile);
					}
					else {
						ERROR("Can not not open '%s' for writing.", of_name);
					}
				END_BLOCK
#endif

				H1ShapesetLobattoHex shapeset;
				PrecalcShapeset pss(&shapeset);

//				printf("* Setting the space up\n");
				H1Space space(&mesh, &shapeset);
				space.set_bc_types(bc_types);
				space.set_bc_values(bc_values);

#ifdef XM_YN_ZO
				int dir_x = 4, dir_y = 4, dir_z = 4;
#elif defined XM_YN_ZO_2
//				int dir_x = 2, dir_y = 3, dir_z = 4;
				int dir_x = 4, dir_y = 4, dir_z = 4;
#elif defined X2_Y2_Z2
				int dir_x = 2, dir_y = 2, dir_z = 2;
#endif
				order3_t o(dir_x, dir_y, dir_z);
//				printf("  - Setting uniform order to (%d, %d, %d)\n", dir_x, dir_y, dir_z);
				space.set_uniform_order(o);

				space.assign_dofs();

//				printf("* Calculating a solution\n");

#if defined WITH_UMFPACK
				UMFPackLinearSolver solver;
#elif defined WITH_PARDISO
				PardisoLinearSolver solver;
#elif defined WITH_PETSC
				PetscLinearSolver solver;
#endif

				Discretization d(&solver);
				d.set_num_equations(1);
				d.set_spaces(1, &space);
				d.set_pss(1, &pss);
#ifdef DIRICHLET
				d.set_bilinear_form(0, 0, bilinear_form);
				d.set_linear_form(0, linear_form);
#elif defined NEWTON
				d.set_bilinear_form(0, 0, bilinear_form, NULL, bilinear_form_surf);
				d.set_linear_form(0, linear_form, linear_form_surf);
#endif

				// assemble siffness matrix
				d.create_stiffness_matrix();

				Timer assemble_timer("Assembling stiffness matrix");
				assemble_timer.start();
				d.assemble_stiffness_matrix_and_rhs();
				assemble_timer.stop();

				// solve the stiffness matrix
				Solution sln(&mesh);
				bool solved = d.solve_system(1, &sln);
				if (!solved) throw ERR_FAILURE;

//				{
//					char file_name[1024];
//					sprintf(file_name, "%s/matrix-%d-%d", OUTPUT_DIR, i, j);
//					FILE *file = fopen(file_name, "w");
//					if (file != NULL) {
//						solver.dump_matrix(file, "A");
//						solver.dump_rhs(file, "b");
//
//						fclose(file);
//					}
//				}


//				double *s = sln.get_solution_vector();
//				for (int ii = 0; ii <= ndofs; ii++)
//					printf("x[% 3d] = % lf\n", ii, s[ii]);


				ExactSolution exsln(&mesh, exact_solution);
				// norm
				double h1_sln_norm = h1_norm(&sln);
				double h1_err_norm = h1_error(&sln, &exsln);
				printf(" - H1 solution norm:   % le\n", h1_sln_norm);
				printf(" - H1 error norm:      % le\n", h1_err_norm);

				double l2_sln_norm = l2_norm(&sln);
				double l2_err_norm = l2_error(&sln, &exsln);
				printf(" - L2 solution norm:   % le\n", l2_sln_norm);
				printf(" - L2 error norm:      % le\n", l2_err_norm);

				assert(h1_sln_norm > 0 && h1_err_norm > 0);
				assert(l2_sln_norm > 0 && l2_err_norm > 0);

//				// out fn
//				char fname[4096];
//				sprintf(fname, "%s/cfg-%d-%d.pos", OUTPUT_DIR, i, j);
//				FILE *fnf = fopen(fname, "w");
//				assert(fnf != NULL);
//				GmshOutputEngine out(fnf);
//				char var[64];
//				sprintf(var, "%d_%d", i, j);
//				out.out(&sln, var);
//				fclose(fnf);
//
//				char mfname[4096];
//				sprintf(mfname, "%s/mesh-%d-%d.ref", OUTPUT_DIR, i, j);
//				FILE *mfnf = fopen(mfname, "w");
//				assert(mfnf != NULL);
//				GmshOutputEngine outm(mfnf);
//				outm.out(&mesh);
//				fclose(mfnf);

				if (h1_err_norm > EPS || l2_err_norm > EPS) {
					// calculated solution is not enough precise
					printf("Solution is not precise enough.\n");
					throw ERR_FAILURE;
				}

				printf("Passed\n");
			}
		}
	}
	catch (int e) {
		res = e;
		printf("Failed\n");
	}

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

