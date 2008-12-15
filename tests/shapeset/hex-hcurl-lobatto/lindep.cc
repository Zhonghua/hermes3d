/*
 * lindep.cc
 * 
 * testing linear indepenedence of a shape functions
 * 
 */

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

// TODO: Mostly the same as for H1 => review the code

// l2 product
double l2_product(RealFunction *fu, RealFunction *fv) {
	Quad3D *quad = fu->get_quad();

	// integrate with maximum order
	int o = quad->get_max_order();

	fu->set_quad_order(o);
	fv->set_quad_order(o);

	double *uval = fu->get_fn_values(); 
	double *vval = fv->get_fn_values(); 

	// integrating over reference brick -> jacobian is 1.0 (we do not have to bother with refmap)
	double result = 0.0;
	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);
	for (int i = 0; i < np; i++)
		result += pt[i].w * (uval[i] * vval[i]);
	
	return result;
}
	

bool test_lin_indep(Shapeset *shapeset) {
	printf("I. linear independency\n");

	UMFPackLinearSolver solver;
	
	PrecalcShapeset pss_u(shapeset), pss_v(shapeset);
	pss_u.set_quad(get_quadrature(MODE));
	pss_v.set_quad(get_quadrature(MODE));

	int n = 
		Hex::NUM_VERTICES * 1 +			// 1 vertex fn 
		Hex::NUM_EDGES * shapeset->get_num_edge_fns(MAX_ELEMENT_ORDER) +
		Hex::NUM_FACES * shapeset->get_num_face_fns(MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER)) +
		shapeset->get_num_bubble_fns(MAKE_HEX_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER));

	printf("number of functions = %d\n", n);
	
	int *fn_idx = new int [n];
	int m = 0;
	// vertex fns
	for (int i = 0; i < Hex::NUM_VERTICES; i++, m++)
		fn_idx[m] = shapeset->get_vertex_index(i);
	// edge fns
	for (int i = 0; i < Hex::NUM_EDGES; i++) {
		int order = MAX_ELEMENT_ORDER;
		int *edge_idx = shapeset->get_edge_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_edge_fns(order); j++, m++)
			fn_idx[m] = edge_idx[j]; 
	}
	// face fns
	for (int i = 0; i < Hex::NUM_FACES; i++) {
		int order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
		int *face_idx = shapeset->get_face_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_face_fns(order); j++, m++)
			fn_idx[m] = face_idx[j]; 
	}
	// bubble
	int order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
	int *bubble_idx = shapeset->get_bubble_indices(order);
	for (int j = 0; j < shapeset->get_num_bubble_fns(order); j++, m++)
		fn_idx[m] = bubble_idx[j]; 
	
		
	// precalc structure
	solver.prealloc(n);	
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			solver.pre_add_ij(i, j);
	solver.alloc();

	printf("assembling matrix ");
	
	for (int i = 0; i < n; i++) {
		pss_u.set_active_shape(fn_idx[i]);

		printf(".");
		fflush(stdout);			// prevent caching of output (to see that it did not freeze)

		for (int j = 0; j < n; j++) {			
			pss_v.set_active_shape(fn_idx[j]);
			
			double value = l2_product(&pss_u, &pss_v); 
			
			solver.update_matrix(i, j, value);
		}
	}
	printf("\n");
	
	for (int i = 0; i < n; i++)
		solver.update_rhs(i, 0.0);

	printf("solving matrix\n");
	
	// solve the system
	double *sln = new double [n];
	if (solver.solve_system(sln)) {
		bool indep = true;
		for (int i = 0; i < n; i++) {
			if (sln[i] >= EPS) {
				indep = false;
				break;
			}
		}
		
		if (indep)
			printf("ok\n");
		else
			printf("Shape functions are not linearly independent\n");
	}
	else {
		printf("Shape functions are not linearly independent\n");
	}

	delete [] sln;
	
	delete [] fn_idx;
	
	return true;
}
