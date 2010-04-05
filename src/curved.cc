// This file is part of Hermes3D
//
// Copyright (c) 2010 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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
#include "curved.h"
#include "matrix.h"
#include "quad.h"
#include "shapeset/lobatto.h"
#include "shapeset.h"
#include "shapeset/refmapss.h"
#include "shapefn.h"

static double **edge_proj_matrix = NULL;	// projection matrix for each edge is the same
static double **face_proj_matrix = NULL;	// projection matrix for each hex face is the same
static double **bubble_proj_matrix = NULL;	// projection matrix for hex bubbles

static double *edge_p = NULL;				// diagonal vector in cholesky factorization
static double *face_p = NULL;				// diagonal vector in cholesky factorization
static double *bubble_p = NULL;				// diagonal vector in cholesky factorization

#ifdef WITH_TETRA
	// TODO: implement me!
	#define REFMAP_SHAPESET_TETRA	NULL
	#define REFMAP_PSS_TETRA		NULL
#else
	#define REFMAP_SHAPESET_TETRA	NULL
	#define REFMAP_PSS_TETRA		NULL
#endif

#ifdef WITH_HEX
	static RefMapShapesetHex 		ref_map_ss;
	static ShapeFunction			ref_map_pss(&ref_map_ss);
	#define REFMAP_SHAPESET_HEX		&ref_map_ss
	#define REFMAP_PSS_HEX			&ref_map_pss
#else
	#define REFMAP_SHAPESET_HEX		NULL
	#define REFMAP_PSS_HEX			NULL
#endif

// FIXME: use this for tetrahedral curvilinear elements (plus fix the code below)
//static ShapeFunction *ref_map_pss[] = { REFMAP_PSS_TETRA, REFMAP_PSS_HEX, NULL };

//// projection based interpolation ////////////////////////////////////////////////////////////////

// preparation of projection matrices, Cholesky factorization
static void precalculate_cholesky_projection_matrix_edge()
{
	int n = 0;
	int indices[] = { 1 };

	edge_proj_matrix = new_matrix<double>(n, n);

	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);
	// calculate projection matrix of maximum order
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			int ii = indices[i], ij = indices[j];

			order3_t ord = ref_map_ss.get_order(ii) + ref_map_ss.get_order(ij);
			order1_t eo = ord.get_edge_order(0);
			QuadPt3D *pt = quad->get_edge_points(0, eo);
			int np = quad->get_edge_num_points(0, eo);

			ref_map_pss.set_active_shape(ii);
			ref_map_pss.precalculate(np, pt, FN_DEFAULT);
			double *fni = ref_map_pss.get_fn_values();

			ref_map_pss.set_active_shape(ij);
			ref_map_pss.precalculate(np, pt, FN_DEFAULT);
			double *fnj = ref_map_pss.get_fn_values();

			double val = 0.0;
			for (int k = 0; k < np; k++)
				val += pt[k].w * (fni[k] * fnj[k]);

			edge_proj_matrix[i][j] = edge_proj_matrix[j][i] = val;
		}
	}

	// Cholesky factorization of the matrix
	edge_p = new double[n];
	choldc(edge_proj_matrix, n, edge_p);
}

static void precalculate_cholesky_projection_matrix_face()
{
	int n = 0;
	int indices[] = { 1 };

	face_proj_matrix = new_matrix<double>(n, n);

	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);
	// calculate projection matrix
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			int ii = indices[i], ij = indices[j];

			order3_t ord = ref_map_ss.get_order(ii) + ref_map_ss.get_order(ij);
			order2_t fo = ord.get_face_order(0);
			QuadPt3D *pt = quad->get_face_points(0, fo);
			int np = quad->get_face_num_points(0, fo);

			ref_map_pss.set_active_shape(ii);
			ref_map_pss.precalculate(np, pt, FN_DEFAULT);
			double *fni = ref_map_pss.get_fn_values();

			ref_map_pss.set_active_shape(ij);
			ref_map_pss.precalculate(np, pt, FN_DEFAULT);
			double *fnj = ref_map_pss.get_fn_values();

			double val = 0.0;
			for (int k = 0; k < np; k++)
				val += pt[k].w * (fni[k] * fnj[k]);

			face_proj_matrix[i][j] = face_proj_matrix[j][i] = val;
		}
	}

	// Cholesky factorization of the matrix
	face_p = new double[n];
	choldc(face_proj_matrix, n, face_p);
}

// calculate the L1 products (\phi_i, \phi_j) for all 0 <= i,j < n, n is the number of bubble functions
static void precalculate_cholesky_projection_matrices_bubble()
{
	int n = 0;
	int indices[] = { 1 };

	bubble_proj_matrix = new_matrix<double>(n, n);

	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			int ii = indices[i], ij = indices[j];

			order3_t o = ref_map_ss.get_order(ii) + ref_map_ss.get_order(ij);
			QuadPt3D *pt = quad->get_points(o);
			int np = quad->get_num_points(o);

			ref_map_pss.set_active_shape(ii);
			ref_map_pss.precalculate(np, pt, FN_DEFAULT);
			double *fni = ref_map_pss.get_fn_values();

			ref_map_pss.set_active_shape(ij);
			ref_map_pss.precalculate(np, pt, FN_DEFAULT);
			double *fnj = ref_map_pss.get_fn_values();

			double val = 0.0;
			for (int k = 0; k < np; k++)
				val += pt[k].w * (fni[k] * fnj[k]);

			bubble_proj_matrix[i][j] = bubble_proj_matrix[j][i] = val;
		}
	}

	// cholesky factorization of the matrix
	bubble_p = new double[n];
	choldc(bubble_proj_matrix, n, bubble_p);
}

///

static void calc_edge_projection(Element *e, int edge, Curve **curves, int order, double2 *proj)
{
/*
	ref_map_pss.set_active_element(e);

  int i, j, k;
  int mo1 = quad1d.get_max_order();
  int np = quad1d.get_num_points(mo1);
  int ne = order - 1;
  int mode = e->get_mode();

  assert(np <= 15 && ne <= 10);
  double2 fn[15];
  double rhside[2][10];
  memset(fn, 0, sizeof(double2) * np);
  memset(rhside[0], 0, sizeof(double) * ne);
  memset(rhside[1], 0, sizeof(double) * ne);

  double a_1, a_2, b_1, b_2;
  a_1 = ctm.m[0] * ref_vert[mode][edge][0] + ctm.t[0];
  a_2 = ctm.m[1] * ref_vert[mode][edge][1] + ctm.t[1];
  b_1 = ctm.m[0] * ref_vert[mode][e->next_vert(edge)][0] + ctm.t[0];
  b_2 = ctm.m[1] * ref_vert[mode][e->next_vert(edge)][1] + ctm.t[1];

  // values of nonpolynomial function in two vertices
  double2 fa, fb;
  calc_ref_map(e, nurbs, a_1, a_2, fa);
  calc_ref_map(e, nurbs, b_1, b_2, fb);

  double2* pt = quad1d.get_points(mo1);
  for (j = 0; j < np; j++) // over all integration points
  {
    double2 x, v;
    double t = pt[j][0];
    edge_coord(e, edge, t, x, v);
    calc_ref_map(e, nurbs, x[0], x[1], fn[j]);

    for (k = 0; k < 2; k++)
      fn[j][k] = fn[j][k] - (fa[k] + (t+1)/2.0 * (fb[k] - fa[k]));
  }

  double2* result = proj + e->nvert + edge * (order - 1);
  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < ne; i++)
    {
      for (j = 0; j < np; j++)
      {
        double t = pt[j][0];
        double fi = lob[i+2](t);
        rhside[k][i] += pt[j][1] * (fi * fn[j][k]);
      }
    }
    // solve
    cholsl(edge_proj_matrix, ne, edge_p, rhside[k], rhside[k]);
    for (i = 0; i < ne; i++)
      result[i][k] = rhside[k][i];
  }
*/
}

//// face part of projection based interpolation /////////////////////////////////////////////////

// TODO

//// bubble part of projection based interpolation /////////////////////////////////////////////////

static void old_projection(Element *e, int order, double2 *proj, double *old[2])
{
/*
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);

  for (unsigned int k = 0; k < e->nvert; k++) // loop over vertices
  {
    // vertex basis functions in all integration points
    double* vd;
    int index_v = ref_map_shapeset.get_vertex_index(k);
    ref_map_pss.set_active_shape(index_v);
    ref_map_pss.set_quad_order(mo2);
    vd = ref_map_pss.get_fn_values();

    for (int m = 0; m < 2; m++)   // part 0 or 1
      for (int j = 0; j < np; j++)
        old[m][j] += proj[k][m] * vd[j];

    for (int ii = 0; ii < order - 1; ii++)
    {
      // edge basis functions in all integration points
      double* ed;
      int index_e = ref_map_shapeset.get_edge_index(k,0,ii+2);
      ref_map_pss.set_active_shape(index_e);
      ref_map_pss.set_quad_order(mo2);
      ed = ref_map_pss.get_fn_values();

      for (int m = 0; m < 2; m++)  //part 0 or 1
        for (int j = 0; j < np; j++)
          old[m][j] += proj[e->nvert + k * (order-1) + ii][m] * ed[j];
    }
  }
*/
}


static void calc_bubble_projection(Element *e, Curve **curves, int order, double2 *proj)
{
/*
  ref_map_pss.set_active_element(e);

  int i, j, k;
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);
  int qo = e->is_quad() ? make_quad_order(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);

  AUTOLA_OR(double2, fn, np);
  memset(fn, 0, sizeof(double2) * np);

  double* rhside[2];
  double* old[2];
  for (i = 0; i < 2; i++) {
    rhside[i] = new double[nb];
    old[i] = new double[np];
    memset(rhside[i], 0, sizeof(double) * nb);
    memset(old[i], 0, sizeof(double) * np);
  }

  // compute known part of projection (vertex and edge part)
  old_projection(e, order, proj, old);

  // fn values of both components of nonpolynomial function
  double3* pt = quad2d.get_points(mo2);
  for (j = 0; j < np; j++)  // over all integration points
  {
    double2 a;
    a[0] = ctm.m[0] * pt[j][0] + ctm.t[0];
    a[1] = ctm.m[1] * pt[j][1] + ctm.t[1];
    calc_ref_map(e, nurbs, a[0], a[1], fn[j]);
  }

  double2* result = proj + e->nvert + e->nvert * (order - 1);
  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < nb; i++) // loop over bubble basis functions
    {
      // bubble basis functions in all integration points
      double *bfn;
      int index_i = ref_map_shapeset.get_bubble_indices(qo)[i];
      ref_map_pss.set_active_shape(index_i);
      ref_map_pss.set_quad_order(mo2);
      bfn = ref_map_pss.get_fn_values();

      for (j = 0; j < np; j++) // over all integration points
        rhside[k][i] += pt[j][2] * (bfn[j] * (fn[j][k] - old[k][j]));
    }

    // solve
    if (e->nvert == 3)
      cholsl(bubble_proj_matrix_tri, nb, bubble_tri_p, rhside[k], rhside[k]);
    else
      cholsl(bubble_proj_matrix_quad, nb, bubble_quad_p, rhside[k], rhside[k]);

    for (i = 0; i < nb; i++)
      result[i][k] = rhside[k][i];
  }

  for (i = 0; i < 2; i++) {
    delete [] rhside[i];
    delete [] old[i];
  }
*/
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void ref_map_projection(Mesh *mesh, Element *e, Curve **curves, const order3_t &order, Vertex *proj)
{
	// vertex part
//	Word_t vtcs[10];
//	e->get_vertices(vtcs);
	for (unsigned int i = 0; i < e->get_num_vertices(); i++) {
//		Vertex *v = mesh->vertices[vtcs[i]];
//		proj[i] = *v;
		proj[i] = *mesh->vertices[e->get_vertex(i)];
	}

	proj[8].x = 0;
	proj[8].y = 0.5;
	proj[8].z = 0;

	proj[24].x = 0;
	proj[24].y = 0.5;
	proj[24].z = 0;


/*	if (e->cm->toplevel == false)
		e = e->cm->parent;

	// edge part
	for (int edge = 0; edge < e->get_num_edges(); edge++)
		calc_edge_projection(e, edge, curves, order, proj);

	// face part
	for (int face = 0; face < e->get_num_faces(); face++)
		calc_face_projection(e, face, curves, order, proj);

	// bubble part
	calc_bubble_projection(e, curves, order, proj);
*/
}


CurvMap::CurvMap(CurvMap *cm)
{
	// NOTE: ugly
	memcpy(this, cm, sizeof(CurvMap));

	coefs = new Vertex[nc];
	memcpy(coefs, cm->coefs, sizeof(Vertex) * nc);
}


CurvMap::~CurvMap()
{
	if (coefs != NULL)
		delete [] coefs;
}

void CurvMap::update_refmap_coefs(Mesh *mesh, Element *e)
{
//	ref_map_pss.set_quad_2d(&quad2d);
//	ref_map_pss.set_active_element(e);

	// calculation of projection matrices
//	if (edge_proj_matrix == NULL) precalculate_cholesky_projection_matrix_edge();
//	if (face_proj_matrix == NULL) precalculate_cholesky_projection_matrix_face();
//	if (bubble_proj_matrix == NULL) precalculate_cholesky_projection_matrices_bubble();

	order = order3_t(3, 3, 3);		// cubic curves

	// allocate projection coefficients
	int n_verts = e->get_num_vertices();
	int n_edges = e->get_num_edges();
	int n_faces = e->get_num_faces();
	int n_edge_fns = 2;						// cubic - 1
	int n_face_fns = 2 * 2;					// (cubic - 1) (cubic - 1)
	int n_bubble_fns = 2 * 2 * 2;
//	int ne = order - 1;
//	int qo = e->is_quad() ? make_quad_order(order, order) : order;
//	int n_bubble_fns = ref_map_ss.get_num_bubbles(qo);
	nc = n_verts + n_edges * n_edge_fns + n_faces * n_face_fns + n_bubble_fns;
	if (coefs != NULL) delete [] coefs;
	coefs = new Vertex[nc];

	// WARNING: do not change the format of the array 'coefs'. If it changes,
	// RefMap::set_active_element() has to be changed too.

	for (int i = 0; i < nc; i++) {
		coefs[i].x = 0;
		coefs[i].y = 0;
		coefs[i].z = 0;
	}

//	printf(" nc = %d\n", nc);

//	Nurbs** nurbs;
//	if (toplevel == false) {
//		ref_map_pss.set_active_element(e);
//		ref_map_pss.set_transform(part);
//		nurbs = parent->cm->nurbs;
//	}
//	else {
//		ref_map_pss.reset_transform();
//		nurbs = e->cm->nurbs;
//	}

//	ctm = *(ref_map_pss.get_ctm());
//	ref_map_pss.reset_transform(); // fixme - do we need this?

	// calculation of new projection coefficients
	ref_map_projection(mesh, e, curves, order, coefs);

	for (int i = 0; i < nc; i++) {
//			printf(" %d", indices[i]);
//		printf("[%lf, %lf, %lf]\n", coefs[i].x, coefs[i].y, coefs[i].z);
//		printf("[%lf, %lf, %lf]\n", e->cm->coefs[i].x, e->cm->coefs[i].y, e->cm->coefs[i].z);
	}
//	die("A");
}
