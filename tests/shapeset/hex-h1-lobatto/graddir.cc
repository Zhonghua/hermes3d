#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

double hh = 1e-3;
double tol = hh * hh * 20.;

int np = 5;
double ptx[5] = {-0.1, 0.6, -0.4, 0.5, 0.3};
double pty[5] = {0., -0.94, 0.42, 0.13, 0.7};
double ptz[5] = {0.8, -0.4, 0.2, -0.7, -0.1};

double dif, maxdif;

bool test(Shapeset *ss, int fnum)
{
	scalar val, valph, valder;
	
	bool passed = true;
	
	for(int point = 0; point < np; point++){
		for(int comp = 0; comp < ss->get_num_components(); comp++){
			// dx
			val = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point], comp);
			valph = ss->get_fn_value(fnum, ptx[point] + hh, pty[point], ptz[point], comp);
			valder = ss->get_dx_value(fnum, ptx[point] + hh/2., pty[point], ptz[point], comp);
			if ( (dif = fabs(valph - val - hh*valder)) > tol){
				passed = false;
				printf("GRADTEST ERROR @(%lf, %lf, %lf), function %d, component %d dx is %lf, should be %lf\n",
						ptx[point], pty[point], ptz[point], fnum, comp, valder, (valph-val)/hh);
			}
			if (dif > maxdif)
				maxdif = dif;
					
			// dy
			val = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point], comp);
			valph = ss->get_fn_value(fnum, ptx[point], pty[point] + hh, ptz[point], comp);
			valder = ss->get_dy_value(fnum, ptx[point], pty[point] + hh/2., ptz[point], comp);
			if ( (dif = fabs(valph - val - hh*valder)) > tol){
				passed = false;
				printf("GRADTEST ERROR @(%lf, %lf, %lf), function %d, component %d dy is %lf, should be %lf\n",
						ptx[point], pty[point], ptz[point], fnum, comp, valder, (valph-val)/hh);
			}
			if (dif > maxdif)
				maxdif = dif;
					
			// dz
			val = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point], comp);
			valph = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point] + hh, comp);
			valder = ss->get_dz_value(fnum, ptx[point], pty[point], ptz[point] + hh/2., comp);
			if ( (dif = fabs(valph - val - hh*valder)) > tol){
				passed = false;
				printf("GRADTEST ERROR @(%lf, %lf, %lf), function %d, component %d dz is %lf, should be %lf\n",
						ptx[point], pty[point], ptz[point], fnum, comp, valder, (valph-val)/hh);
			}
			if (dif > maxdif)
				maxdif = dif;
			
		}
	}
	return passed;
}

bool test_gradients_directly(Shapeset *ss)
{
	printf("V. direct check of the gradient values\n");

	maxdif = 0.;
	bool passed = true;
	int index, *indices, ii;
	
	for(int iv = 0; iv < 8; iv++){
		index = ss->get_vertex_index(iv);
		if(! test(ss, index)) 
			passed = false;
	}
	
	int order = MAX_ELEMENT_ORDER;
	for(int ie = 0; ie < 12; ie++){
		for(int ori = 0; ori < 2; ori++){ 
			indices = ss->get_edge_indices(ie, ori, order);
			for(ii = 0; ii < ss->get_num_edge_fns(order); ii++)
				if(! test(ss, indices[ii]))
					passed = false;
		}
	}
	
	order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
	for(int ic = 0; ic < 6; ic++){
		for(int ori = 0; ori < 8; ori++){ 
			indices = ss->get_face_indices(ic, ori, order);
			for(ii = 0; ii < ss->get_num_face_fns(order); ii++)
				if(! test(ss, indices[ii]))
					passed = false;
		}
	}
	
	order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
	indices = ss->get_bubble_indices(order);
	for(ii = 0; ii < ss->get_num_bubble_fns(order); ii++)
		if(! test(ss, indices[ii]))
			passed = false;
	
	printf("maximal difference is %g, which is %g * h^2\n", maxdif, maxdif/hh/hh);
	
	return passed;
}
























