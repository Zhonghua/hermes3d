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
//#define NEWTON

#define X2_Y2_Z2
//#define XM_YN_ZO
//#define XM_YN_ZO_2

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

double bc_values(int marker, double x, double y, double z) {
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


////////////////////////////////////////////////////////////////////////////////

// maximal level of refinement considered (= 0 .. # of ref)
#define MAX_LEVEL						4
#define NUM_RULES						((MAX_LEVEL + 1) * (MAX_LEVEL + 1))

// special quadrature used in this test to define set of points, where
// continuity is tested. Quadrature is used in order to use RefMap abilities
// calculate physical coordinates of points from the reference domain.
//
// points are chosen in this way:
// - only face points are used
// - there are more levels of points
// - 0-level consist of 16 points chosen symmetricaly
// - other levels (or orders) simulate division of adjacent elements and map
//   points in such way, that some points from divided and some from nondivided
//   face will still match in the physical domain
// - this is of course only possible up to some level of division, controled by
//   constant LEVEL. for example, for hex1 :
//   LEVEL = 1  divisions 0 x 1 y or 0 x 1 y 3 z are ok, but 0 x 1 y 3 y not
//   LEVEL = 2  division 0 x 1 y 3 y is ok, but 0 x 1 y 3 y 5 y not
// - if LEVEL is not sufficient, there will be some faces, that will not be tested,
//   because no points from the face will match to points from the constraining face
class ContQuad : public Quad3D {
public:
	ContQuad() {
//		max_order = max_edge_order = max_face_order = NUM_RULES;

		int my_np_1d[MAX_LEVEL + 1];
		double my_tables_1d[MAX_LEVEL + 1][100];

		my_np_1d[0] = 4;
		my_tables_1d[0][0] = -0.71;  my_tables_1d[0][1] = -0.59;
		my_tables_1d[0][2] =  0.59;  my_tables_1d[0][3] =  0.71;

		for(int order = 1; order <= MAX_LEVEL; order++) {
			my_np_1d[order] = 2 * my_np_1d[order - 1];
			for(int i = 0; i < my_np_1d[order - 1]; i++) {
				my_tables_1d[order][i] = (my_tables_1d[order - 1][i] - 1.) / 2.;
				my_tables_1d[order][i + my_np_1d[order - 1]] = (my_tables_1d[order - 1][i] + 1.) / 2.;
			}
		}

		face_tables = new Array<QuadPt3D *>[Hex::NUM_FACES];
		for (int order = 0; order < NUM_RULES; order++) {
			int ord1 = order / (MAX_LEVEL + 1);
			int ord2 = order % (MAX_LEVEL + 1);
			np_face[order] = my_np_1d[ord1] * my_np_1d[ord2];
			for (int face = 0; face < Hex::NUM_FACES; face++)
				face_tables[face][order] = new QuadPt3D[np_face[order]];

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[0][order][n].x = -1;
					face_tables[0][order][n].y = my_tables_1d[ord1][k];
					face_tables[0][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[1][order][n].x = 1;
					face_tables[1][order][n].y = my_tables_1d[ord1][k];
					face_tables[1][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[2][order][n].x = my_tables_1d[ord1][k];
					face_tables[2][order][n].y = -1;
					face_tables[2][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[3][order][n].x = my_tables_1d[ord1][k];
					face_tables[3][order][n].y = 1;
					face_tables[3][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[4][order][n].x = my_tables_1d[ord1][k];
					face_tables[4][order][n].y = my_tables_1d[ord2][l];
					face_tables[4][order][n].z = -1;
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[5][order][n].x = my_tables_1d[ord1][k];
					face_tables[5][order][n].y = my_tables_1d[ord2][l];
					face_tables[5][order][n].z = 1;
				}
			}
		}
	}

	~ContQuad() {
		for (int face = 0; face < Hex::NUM_FACES; face++)
			for (int order = 0; order < NUM_RULES; order++)
				delete [] face_tables[face][order];
		delete [] face_tables;
	}
};

struct Point {
	double ref_x, ref_y, ref_z;
	double phys_x, phys_y, phys_z;
	Word_t elm_idx, fac_idx;

	Point(int idx, int f_idx, double rx, double ry, double rz, double px, double py, double pz) {
		elm_idx = idx; fac_idx = f_idx;
		ref_x = rx;  ref_y = ry;  ref_z = rz;
		phys_x = px; phys_y = py; phys_z = pz;
	}
};

typedef
	int (*compfn)(const void*, const void*);

int compare(Point **pt1, Point **pt2) {
	double val1 = 1000000. * (*pt1)->phys_x + 1000. * (*pt1)->phys_y + (*pt1)->phys_z;
	double val2 = 1000000. * (*pt2)->phys_x + 1000. * (*pt2)->phys_y + (*pt2)->phys_z;

	if (val1 < val2) return -1;
	else if (val1 > val2) return 1;
	else return 0;
}

const double EPS = 1e-13;
const double TOLERANCE = 1e-10;

bool equal(Point *pt1, Point *pt2) {
	if (fabs(pt1->phys_x - pt2->phys_x) > TOLERANCE) return false;
	if (fabs(pt1->phys_y - pt2->phys_y) > TOLERANCE) return false;
	if (fabs(pt1->phys_z - pt2->phys_z) > TOLERANCE) return false;
	return true;
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;


#ifdef WITH_PETSC
	PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 2) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

	// apply refinements
	for (int i = 2; i < argc; i += 2) {
		int elem_id, reft_id;
		sscanf(args[i], "%d", &elem_id);
		reft_id = parse_reft(args[i + 1]);
		mesh.refine_element(elem_id + 1, reft_id);
	}
	mesh.dump();


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
//	shapeset.preload_products();

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_bc_values(bc_values);

#ifdef XM_YN_ZO
	int dir_x = 4, dir_y = 4, dir_z = 4;
#elif defined XM_YN_ZO_2
	int dir_x = 2, dir_y = 3, dir_z = 4;
#elif defined X2_Y2_Z2
	int dir_x = 2, dir_y = 2, dir_z = 2;
#endif
	order3_t o(dir_x, dir_y, dir_z);
	printf("  - Setting uniform order to (%d, %d, %d)\n", dir_x, dir_y, dir_z);
	space.set_uniform_order(o);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

#if defined WITH_UMFPACK
	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(mat, rhs);
#elif defined WITH_PARDISO
	PardisoLinearSolver solver;
#elif defined WITH_PETSC
	PetscMatrix mat;
	PetscVector rhs;
	PetscLinearSolver solver(mat, rhs);
#endif

	Discretization d;
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

	// assemble stiffness matrix
	d.create(&mat, &rhs);
	d.assemble(&mat, &rhs);

#ifdef OUTPUT_DIR
	{
		char file_name[1024];
		sprintf(file_name, "%s/matrix", OUTPUT_DIR);
		FILE *file = fopen(file_name, "w");
		if (file != NULL) {
			mat.dump(file, "A");
			rhs.dump(file, "b");

			fclose(file);
		}
	}
#endif

	try {
		// solve the stiffness matrix
		bool solved = solver.solve();

		if (!solved) throw ERR_FAILURE;

		Solution sln(&mesh);
		sln.set_space_and_pss(&space, &pss);
		sln.set_solution_vector(solver.get_solution(), false);

		{
			double *s = sln.get_solution_vector();
			for (int i = 0; i <= ndofs; i++)
				printf("x[% 3d] = % lf\n", i, s[i]);
		}

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

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
			printf("Solution is not precise enough.\n");
		}

		//
		//
		// the main code starts here
		//
		//
#if 1
		ContQuad my_quad;
		RefMap ref_map(&mesh);

		int num_points = 0;
		for (int order = 0; order < NUM_RULES; order++) {
			FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
				for (int iface = 0; iface < Hex::NUM_FACES; iface++)
					num_points += my_quad.get_face_num_points(iface, order);
			}
		}

		Point **points = new Point *[num_points];
		int ipt = 0;

		// find points
		for(int order = 0; order < NUM_RULES; order++) {
			FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
				Element *e = mesh.elements[idx];
				ref_map.set_active_element(e);
				ref_map.set_quad(&my_quad);
				for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
					Word_t fac_idx = mesh.get_facet_id(e, iface);

					QuadPt3D *quad_pts = my_quad.get_face_points(iface, order);
					int np = my_quad.get_face_num_points(iface, order);
					double *phys_x = ref_map.get_face_phys_x(iface, order);
					double *phys_y = ref_map.get_face_phys_y(iface, order);
					double *phys_z = ref_map.get_face_phys_z(iface, order);

					// for each face and each integration point store reference and physical coordinates
					for (int pt = 0; pt < np; pt++) {
						points[ipt++] = new Point(idx, fac_idx,
							quad_pts[pt].x, quad_pts[pt].y, quad_pts[pt].z,
							phys_x[pt], phys_y[pt], phys_z[pt]);
					}
				}
			}
		}

		// sort points according to first phys_x, then phys_y and phys_z
		// it means, that two points, with almost identical physical coordinates
		// (even though from different elements) will be next to each other in the array
		qsort((void *) points, num_points, sizeof(Point*), (compfn)compare);

		int *pairs = new int [num_points];
		int num_pairs = 0;

		// choose those indicies, that correspond to pairs with identical physical coordinates
		// and store them in field pairs
		for (int i = 0; i < num_points - 1; i++) {
			if (equal(points[i], points[i+1])) {
				pairs[num_pairs++] = i;
			}
		}

		// check, whether we tested points from all inner active facets
		// this is done only for testing of correctness of the test itself
		int nonchecked_faces = 0;
		FOR_ALL_FACETS(fid, &mesh) {
			bool ok = false;
			Facet *fac = mesh.facets[fid];
			if (fac->type == Facet::OUTER) continue;
			if (!(fac->ractive || fac->lactive)) continue;
			for (int i = 0; i < num_pairs; i++) {
				if ((points[pairs[i]]->fac_idx == fid) || (points[pairs[i] + 1]->fac_idx == fid)) {
					ok = true;
					break;
				}
			}

			if (!ok) nonchecked_faces++;
		}


		// loop over all basis functions
		for (int dof = 0; dof < ndofs; dof++) {
			printf("processing dof %d...\n", dof);

			// prepare solution correspondig to basis function with dof dof
			double sln_vector[ndofs + 1];
			memset(sln_vector, 0, (ndofs + 1) * sizeof(double));
			sln_vector[dof] = 1.0;
			sln.set_solution_vector(sln_vector, false);

			double max_difference = 0.;
			double max_pt_x, max_pt_y, max_pt_z, max_val_1, max_val_2;
			Word_t max_elm_1, max_elm_2;

			// loop over all pairs of points, that correspond to one point in the physical domain
			for(int pair = 0; pair < num_pairs; pair++) {
				int i = pairs[pair];
				Element *e1 = mesh.elements[points[i]->elm_idx];
				sln.set_active_element(e1);
				double val1 = sln.get_sln_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, FN, 0);

				Element *e2 = mesh.elements[points[i + 1]->elm_idx];
				sln.set_active_element(e2);
				double val2 = sln.get_sln_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, FN, 0);

				//value in this point should be the same, no matter from which element we go
				double difference = fabs(val1 - val2);
				if (difference > max_difference) {
					max_difference = difference;
					max_pt_x = points[i]->phys_x;
					max_pt_y = points[i]->phys_y;
					max_pt_z = points[i]->phys_z;
					max_val_1 = val1;
					max_val_2 = val2;
					max_elm_1 = points[i]->elm_idx;
					max_elm_2 = points[i + 1]->elm_idx;
				}
			}

			if (max_difference > TOLERANCE) {
				printf("base fn %d NOT continuous between elements %ld and %ld @ (% lf, % lf, % lf), max difference %g (%.15g <-> %.15g)\n",
						 dof, max_elm_1, max_elm_2 , max_pt_x, max_pt_y, max_pt_z, max_difference, max_val_1, max_val_2);
				res = ERR_FAILURE;
			}

		}

		delete [] pairs;
		delete [] points;

		printf("continuity tested in %d points and %d inner faces with at least one active adjacat element were not tested\n", num_pairs, nonchecked_faces);
#endif
		if (res != ERR_SUCCESS) throw res;

		printf("Passed\n");
	}
	catch (int e) {
		res = e;
		printf("Failed\n");
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

