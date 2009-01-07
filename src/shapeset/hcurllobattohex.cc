//
// hcurl_lobatto_hex.cc
//

#include "../config.h"
#include "../refdomain.h"
#include "common.h"
#include "lobatto.h"
#include "hcurllobattohex.h"
#include <common/error.h>
#include "matrix.h"

#include "mesh.h"

static const int base_coding = MAX_ELEMENT_ORDER + 1;		// order: 0..10

int HCurlShapesetLobattoHex::max_fns_per_edge;
int HCurlShapesetLobattoHex::max_fns_per_face;
int HCurlShapesetLobattoHex::max_bubble_fns;

int HCurlShapesetLobattoHex::face_offset;
int HCurlShapesetLobattoHex::bubble_offset;


// t_direction is direction of edge, v_direction is direction of vectors of shapefunction
// for the edge functions they are the same
void HCurlShapesetLobattoHex::check_edge_fn_index(int edge_fn_index, int indices[3], int &edge_index, int &which_legendre, int &t_direction, int &v_direction) {
	//index of edge
	edge_index = edge_fn_index / (max_fns_per_edge * 8);
	//local index of edge function on the edge
	//int n = edge_fn_index - max_fns_per_edge * edge_index;
	int n = (edge_fn_index / 8) % max_fns_per_edge;

	switch(edge_index){

		//edge e1
		case 0: //E[0] = legendre_fn_tab_1d[n](x) * lobatto_fn_tab_1d[0](y) * lobatto_fn_tab_1d[0](z);
			t_direction = 0;
			indices[1] = 0;
			indices[2] = 0;
			break;

		//edge e2
		case 1: //E[1] = lobatto_fn_tab_1d[1](x) * legendre_fn_tab_1d[n](y) * lobatto_fn_tab_1d[0](z);
			t_direction = 1;
			indices[0] = 1;
			indices[2] = 0;
			break;

		//edge e3
		case 2: //E[0] = legendre_fn_tab_1d[n](x) * lobatto_fn_tab_1d[1](y) * lobatto_fn_tab_1d[0](z);
			t_direction = 0;
			indices[1] = 1;
			indices[2] = 0;
			break;

		//edge e4
		case 3: //E[1] = lobatto_fn_tab_1d[0](x) * legendre_fn_tab_1d[n](y) * lobatto_fn_tab_1d[0](z);
			t_direction = 1;
			indices[0] = 0;
			indices[2] = 0;
			break;

		//edge e5
		case 4: //E[2] = lobatto_fn_tab_1d[0](x) * lobatto_fn_tab_1d[0](y) * legendre_fn_tab_1d[n](z);
			t_direction = 2;
			indices[0] = 0;
			indices[1] = 0;
			break;

		//edge e6
		case 5: //E[2] = lobatto_fn_tab_1d[1](x) * lobatto_fn_tab_1d[0](y) * legendre_fn_tab_1d[n](z);
			t_direction = 2;
			indices[0] = 1;
			indices[1] = 0;
			break;

		//edge e7
		case 6: //E[2] = lobatto_fn_tab_1d[1](x) * lobatto_fn_tab_1d[1](y) * legendre_fn_tab_1d[n](z);
			t_direction = 2;
			indices[0] = 1;
			indices[1] = 1;
			break;

		//edge e8
		case 7: //E[2] = lobatto_fn_tab_1d[0](x) * lobatto_fn_tab_1d[1](y) * legendre_fn_tab_1d[n](z);
			t_direction = 2;
			indices[0] = 0;
			indices[1] = 1;
			break;

		//edge e9
		case 8: //E[0] = legendre_fn_tab_1d[n](x) * lobatto_fn_tab_1d[0](y) * lobatto_fn_tab_1d[1](z);
			t_direction = 0;
			indices[1] = 0;
			indices[2] = 1;
			break;

		//edge e10
		case 9: //E[1] = lobatto_fn_tab_1d[1](x) * legendre_fn_tab_1d[n](y) * lobatto_fn_tab_1d[1](z);
			t_direction = 1;
			indices[0] = 1;
			indices[2] = 1;
			break;

		//edge e11
		case 10: //E[0] = legendre_fn_tab_1d[n](x) * lobatto_fn_tab_1d[1](y) * lobatto_fn_tab_1d[1](z);
			t_direction = 0;
			indices[1] = 1;
			indices[2] = 1;
			break;

		//edge e12
		case 11: //E[1] = lobatto_fn_tab_1d[0](x) * legendre_fn_tab_1d[n](y) * lobatto_fn_tab_1d[1](z);
			t_direction = 1;
			indices[0] = 0;
			indices[2] = 1;
			break;
		default: ERROR("Wrong edge index in check_edge_fn_index(). Index %d", edge_fn_index);
   }

	indices[t_direction] = n;
	which_legendre = t_direction; //legendre polynomial allways in the direction of edge
	v_direction = t_direction; //vectors of function in the direction of edge
	assert(which_legendre < 3);
}

// t_direction_1 and _2 are two tangent vectors to the face corresponding to the face function
// v_direction is direction of vectors of the shape function. It is allways either t_dir_1 or t_dir_2
void HCurlShapesetLobattoHex::check_face_fn_index(int face_fn_index, int indices[3],  int &face_index,
		int &which_legendre, int &t_direction_1, int &t_direction_2, int &v_direction) {
	//Face functions are ordered as follows: 0..max_fns_per_face-1 are face functions on the face s1,
	//max_fns_per_face..2*max_fns_per_face-1 face functions corresponding to face s2 etc.
	//Calculate index of face (0 <= face_index <= NUMBER_OF_FACES-1):
	face_index = face_fn_index / (max_fns_per_face * 8);

	//There are two groups of face functions on every face [see relation (2.67) in SoSeDo]
	//  Group 0: n1 = 0..max_face_order, n2 = 2..max_face_order+1.
	//  Group 1: n1 = 2..max_face_order+1, n2 = 0..max_face_order
	//Calculate the relative index of the face function on the face it belongs to
	//int n = face_fn_index - max_fns_per_face * face_index;
	int n = (face_fn_index / 8) % max_fns_per_face;

	//if in first half then group = 0, otherwise group = 1
	int group;
	if (n < max_fns_per_face/2) group = 0;
	else {
		group = 1;
		n -= max_fns_per_face/2; //now n becomes the index of the face function within the second group
	}
	//Next decompose n into the indices of the Legendre polynomials and the Lobatto functions.
	//This decomposition is different for each group.
	//For the first group the table has max_face_order+1 columns (Legendre polynomials starting with L0)
	//and max_face_order rows (Lobatto functions starting with l2). Here, n1 indexes the Legendre polynomials
	//and n2 the Lobatto functions.
	//For the second group the table has max_face_order columns (Lobatto functions starting with l2)
	//and max_face_order+1 rows (Legendre polynomials starting with L0). Here, n1 indexes the Legendre
	//polynomials and n2 the Lobatto functions.
	//The faster-running index in both tables is n2
	int n1, n2;
	if (group == 0) {
		//n1 ... indices of Legendre polynomials, going from 0 to max_face_order
		//n2 ... indices of Lobatto functions, going from 2 to max_face_order + 1
		n1 = n / MAX_ELEMENT_ORDER;
		n2 = n % MAX_ELEMENT_ORDER;
		n2 += 2; //index shift for Lobatto functions
	} else {
		//n1 ... indices of Lobatto functions, going from 2 to max_face_order + 1
		//n2 ... indices of Legendre polynomials, going from 0 to max_face_order
		n1 = n / (MAX_ELEMENT_ORDER + 1);
		n2 = n % (MAX_ELEMENT_ORDER + 1);
		n1 += 2; //index shift for Lobatto functions
	}

	switch(face_index) {
		case 0:
			// (group=0) E[1]=lobatto_fn_tab_1d[0](x)*legendre_fn_tab_1d[n1](y)*lobatto_fn_tab_1d[n2](z); E[0]=E[2]=0
			//or (group=1) E[2]=lobatto_fn_tab_1d[0](x)*lobatto_fn_tab_1d[n1](y)*legendre_fn_tab_1d[n2](z); E[0]=E[1]=0
			t_direction_1 = 1;
			t_direction_2 = 2;
			indices[0] = 0;
			break;
		case 1:
			//E[1] = lobatto_fn_tab_1d[1](x) * legendre_fn_tab_1d[n1](y) * lobatto_fn_tab_1d[n2](z);
			//or E[2] = lobatto_fn_tab_1d[1](x) * lobatto_fn_tab_1d[n1](y) * legendre_fn_tab_1d[n2](z);
			t_direction_1 = 1;
			t_direction_2 = 2;
			indices[0] = 1;
			break;
		case 2:
			//E[0] = legendre_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[0](y) * lobatto_fn_tab_1d[n2](z);
			//or E[2] = lobatto_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[0](y) * legendre_fn_tab_1d[n2](z);
			t_direction_1 = 0;
			t_direction_2 = 2;
			indices[1] = 0;
			break;
		case 3:
			//E[0] = legendre_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[1](y) * lobatto_fn_tab_1d[n2](z);
			//or E[2] = lobatto_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[1](y) * legendre_fn_tab_1d[n2](z);
			t_direction_1 = 0;
			t_direction_2 = 2;
			indices[1] = 1;
			break;
		case 4:
			//E[0] = legendre_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[n2](y) * lobatto_fn_tab_1d[0](z);
			//or E[1] = lobatto_fn_tab_1d[n1](x) * legendre_fn_tab_1d[n2](y) * lobatto_fn_tab_1d[0](z);
			t_direction_1 = 0;
			t_direction_2 = 1;
			indices[2] = 0;
			break;
		case 5:
			//E[0] = legendre_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[n2](y) * lobatto_fn_tab_1d[1](z);
			//or E[1] = lobatto_fn_tab_1d[n1](x) * legendre_fn_tab_1d[n2](y) * lobatto_fn_tab_1d[1](z);
			t_direction_1 = 0;
			t_direction_2 = 1;
			indices[2] = 1;
			break;
		default: ERROR("Wrong face index in check_face_fn_index().");
	}

	indices[t_direction_1] = n1;
	indices[t_direction_2] = n2;
	if(group == 0){
		which_legendre = t_direction_1;
		v_direction = t_direction_1;
	}
	else{
		which_legendre = t_direction_2;
		v_direction = t_direction_2;
	}
}

void HCurlShapesetLobattoHex::check_bubble_fn_index(int bubble_fn_index, int indices[3], int &which_legendre, int &v_direction) {
	//Bubble functions are split into three groups according to the spatial component in which
	//they are nonzero. By 'n' we denote the local index withing the group;
	int group, n;

	bubble_fn_index /= 8;

	if (bubble_fn_index < max_bubble_fns / 3) {
		group = 0;
		n = bubble_fn_index;
	}
	else {
		if (bubble_fn_index < 2 * max_bubble_fns / 3) {
			group = 1;
			n = bubble_fn_index - max_bubble_fns / 3;
		}
		else {
			group = 2;
			n = bubble_fn_index - 2 * max_bubble_fns / 3;
		}
	}

	int n1, n2, n3;

	//For each group, the decomposition of the local index 'n' into the indices of
	//Lobatto functions and Legendre polynomials is different. Always, n3 is the
	//fastest-running index, n2 the medium fast, and n1 the slowest.
	switch (group) {
	case 0: //In this group, n1 goes from 0 to MAX_ELEMENT_ORDER, n2 and n3 from 2 to MAX_ELEMENT_ORDER + 1
		n1 = n / (MAX_ELEMENT_ORDER * MAX_ELEMENT_ORDER);
		n2 = (n / MAX_ELEMENT_ORDER) % MAX_ELEMENT_ORDER;
		n3 = n % MAX_ELEMENT_ORDER;
		n2 += 2; //index shift for the Lobatto functions
		n3 += 2; //index shift for the Lobatto functions
		//E[0] = legendre_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[n2](y) * lobatto_fn_tab_1d[n3](z);
		//E[1] = 0;
		//E[2] = 0;
		break;

	case 1: //In this group, n2 goes from 0 to MAX_ELEMENT_ORDER, n1 and n3 from 2 to MAX_ELEMENT_ORDER + 1
		n1 = n / ((MAX_ELEMENT_ORDER + 1) * MAX_ELEMENT_ORDER);
		n2 = (n / MAX_ELEMENT_ORDER ) % (MAX_ELEMENT_ORDER + 1);
		n3 = n % MAX_ELEMENT_ORDER;
		n1 += 2; //index shift for the Lobatto functions
		n3 += 2; //index shift for the Lobatto functions
		//E[0] = 0;
		//E[1] = lobatto_fn_tab_1d[n1](x) * legendre_fn_tab_1d[n2](y) * lobatto_fn_tab_1d[n3](z);
		//E[2] = 0;
		break;
	case 2: //In this group, n3 goes from 0 to MAX_ELEMENT_ORDER, n1 and n2 from 2 to MAX_ELEMENT_ORDER + 1
		n1 = n / ((MAX_ELEMENT_ORDER + 1) * MAX_ELEMENT_ORDER);
		n2 = (n / (MAX_ELEMENT_ORDER + 1)) % MAX_ELEMENT_ORDER;
		n3 = n % (MAX_ELEMENT_ORDER +1);
		n1 += 2; //index shift for the Lobatto functions
		n2 += 2; //index shift for the Lobatto functions
		//E[0] = 0;
		//E[1] = 0;
		//E[2] = lobatto_fn_tab_1d[n1](x) * lobatto_fn_tab_1d[n2](y) * legendre_fn_tab_1d[n3](z);
		break;
	default: ERROR("Wrong bubble function group in check_bubble_fn_index().");
	}
	indices[0] = n1;
	indices[1] = n2;
	indices[2] = n3;
	v_direction = group;
	which_legendre = group;

}

// fn_type... 0 for edge functions, 1 for face functions, 2 for bubble functions
// unit_index ... edge index if fn_type == 0, face index if fn_type == 1, -1 if fn_type == 2
void HCurlShapesetLobattoHex::check_fn_index(int index, int indices[3], int &fn_type, int &unit_index,
									int &which_legendre, int &v_direction)
{
	assert(index >= 0);

	int t_dir_1, t_dir_2;
	if(index < face_offset) {
		fn_type = SHAPE_FN_EDGE;
		check_edge_fn_index(index, indices, unit_index, which_legendre, t_dir_1, v_direction);
	}
	else if(index < bubble_offset) {
		fn_type = SHAPE_FN_FACE;
		check_face_fn_index(index - face_offset, indices, unit_index, which_legendre, t_dir_1, t_dir_2, v_direction);
	}
	else if (index < bubble_offset + 8 * max_bubble_fns) {
		fn_type = SHAPE_FN_BUBBLE;
		unit_index = -1;
		check_bubble_fn_index(index - bubble_offset, indices, which_legendre, v_direction);
	}
	else {
		ERROR("Type of shape function wrong in check_fn_index() : %d", index);
		assert(0);
	}
	assert(which_legendre < 3);
}

void HCurlShapesetLobattoHex::dump_fn_index(int index)
{
	int ind[3], fn_type, unit_index, which_leg, v_dir;
	check_fn_index(index, ind, fn_type, unit_index, which_leg, v_dir);
	printf("index %d -> (%d, %d, %d), type %d, unit index %d, which_legendre %d, vector direction %d\n", index, ind[0], ind[1], ind[2], fn_type, unit_index, which_leg, v_dir);
}


double HCurlShapesetLobattoHex::calc_any_value(int index, double x, double y, double z, int component, int which_der)
{
	int ori = GET_ORI_FROM_INDEX(index);
	int idx = GET_IDX_FROM_INDEX(index);

	int indices[3], oris[3];
	int fn_type, unit_index, which_legendre, v_direction;
	double point[3] = {x, y, z};

	check_fn_index(index, indices, fn_type, unit_index, which_legendre, v_direction);

//	printf("calc_any_value -> index %d, idx %d, ori %d, leg %d, v_dir %d\n", index, idx, ori, which_legendre, v_direction);
//	getchar();

	oris[0] = oris[1] = oris[2] = 0;
	if(fn_type == SHAPE_FN_EDGE){
		//edge fn
		assert(ori < 2);
		oris[v_direction] = ori;
	}
	else if(fn_type == SHAPE_FN_FACE){
		//face fn
		assert(ori < 8);
		int dir_1 = RefHex::get_face_tangent_direction(unit_index, 0);
		int dir_2 = RefHex::get_face_tangent_direction(unit_index, 1);

		if (ori % 2 == 1)
			oris[dir_1] = 1;
		if (ori % 4 >= 2)
			oris[dir_2] = 1;
		if (ori >= 4) {
			swapint(indices[dir_1], indices[dir_2]);
			swapint(oris[dir_1], oris[dir_2]);
			which_legendre = (which_legendre == dir_1) ? dir_2 : dir_1;
			//in our shapeset this is allways true
			v_direction = which_legendre;
		}
	}
	else{
		if(ori != 0)
			printf("index %d, fn_type %d, unit_index %d\n", index, fn_type, unit_index);
		assert(ori == 0);
	}

	if(v_direction != component)
		return 0.;

	if (oris[0] == 1) point[0] = -point[0];
	if (oris[1] == 1) point[1] = -point[1];
	if (oris[2] == 1) point[2] = -point[2];

	double value = 1.;

	for(int i = 0; i < 3; i++){
		if(which_legendre == i){
			if(which_der == i)
				value *= legendre_der_tab_1d[indices[i]](point[i]);
			else
				value *= legendre_fn_tab_1d[indices[i]](point[i]);
		}
		else{
			if(which_der == i)
				value *= lobatto_der_tab_1d[indices[i]](point[i]);
			else
				value *= lobatto_fn_tab_1d[indices[i]](point[i]);
		}
	}

	for(int i = 0; i < 3; i++)
		oris[i] = 1 - 2 * oris[i];

	if((which_der == -1) && (oris[component] == -1))  //POZOR, zmena 1->-1
		value = -value;
	if(which_der >= 0)
		value = oris[component] * oris[which_der] * value;

	return value;
}


double HCurlShapesetLobattoHex::calc_fn_value(int index, double x, double y, double z, int component)
{
	return calc_any_value(index, x, y, z, component, -1);
}

double HCurlShapesetLobattoHex::calc_dx_value(int index, double x, double y, double z, int component)
{
	return calc_any_value(index, x, y, z, component, 0);
}

double HCurlShapesetLobattoHex::calc_dy_value(int index, double x, double y, double z, int component)
{
	return calc_any_value(index, x, y, z, component, 1);
}

double HCurlShapesetLobattoHex::calc_dz_value(int index, double x, double y, double z, int component)
{
	return calc_any_value(index, x, y, z, component, 2);
}


/******************************************************************************
 *  DETERMINING INDICES OF SHAPE FUNCTIONS WHOSE POLYNOMIAL DEGREES ARE LESS  *
 *      THAN OR EQUAL TO A GIVEN CONSTANT (NEEDED IN ASSEMBLING PROCEDURE)    *
 ******************************************************************************/

//NOTE: If the user changes shape functions, this does not need to be changed,
//      as long as there is NUM_OF_EDGES*max_fns_per_edge edge functions, then
//      NUM_OF_FACES*max_fns_per_face face functions, and last max_bubble_fns bubble
//      functions. Functions compute_edge_indices(), compute_face_indices()
//      and compute_bubble_indices() use the functions check_edge_fn_index(),
//      check_face_fn_index(), check_bubble_fn_index().

//For a given "edge" and a given "order", this returns the indices of
//all edge functions associated with "edge" whose tangential polynomial
//degree is less or equal to "order"
void HCurlShapesetLobattoHex::compute_edge_indices(int edge, int wanted_ori, int order)
{
	int *indices = new int[get_num_edge_fns(order)];
	int position = 0;

	int fn_degree[3];
	int unit_index = -1;
	int which_legendre, t_direction, v_direction;
	//loop over all edge functions
	for(int edge_fn_index=0; edge_fn_index< Hex::NUM_EDGES * max_fns_per_edge * 8; edge_fn_index++) {

		int ori = GET_ORI_FROM_INDEX(edge_fn_index);
		int idx = GET_IDX_FROM_INDEX(edge_fn_index);
		if(ori != wanted_ori)
			continue;

		//t_direction: 0... edge is in x-direction, 1... edge is in y-direction, 2... edge is in z-direction
		int t_direction, v_direction, whcih_legendre;
		check_edge_fn_index(edge_fn_index, fn_degree, unit_index, which_legendre, t_direction, v_direction);
		//tangential degree along "edge"
		int t_degree = fn_degree[t_direction];
		if((unit_index == edge) && (t_degree <= order)) {
			if(position >= get_num_edge_fns(order))
				ERROR("Internal in compute_edge_indices(). position = %d, order = %d, num_edg_fns = %d", position, order, get_num_edge_fns(order));
			indices[position++] = edge_fn_index;
		}
	}

	assert(position == get_num_edge_fns(order));
	edge_indices[edge][wanted_ori][order] = indices;
}

void HCurlShapesetLobattoHex::compute_face_indices(int face, int wanted_ori, order2_t order)
{
	int order_1 = order.x;
	int order_2 = order.y;
	int *indices = new int[get_num_face_fns(order)];
	int position = 0;

	int fn_degree[3];
	int unit_index = -1;
	int which_legendre, t_direction_1, t_direction_2, v_direction;
	//loop over all face functions
	for(int face_fn_index=0; face_fn_index < Hex::NUM_FACES * max_fns_per_face * 8; face_fn_index++) {

		int ori = GET_ORI_FROM_INDEX(face_fn_index);
		int idx = GET_IDX_FROM_INDEX(face_fn_index);
		if(ori != wanted_ori)
			continue;

		//t_directions: 0,1... face is parallel to the x_1x_2 plane, etc.
		int t_direction_1, t_direction_2;
		check_face_fn_index(face_fn_index, fn_degree, unit_index, which_legendre, t_direction_1, t_direction_2, v_direction);

		//this is according to the definition in the book
		int which_lobatto = (t_direction_1 == which_legendre) ? t_direction_2 : t_direction_1;
		fn_degree[which_lobatto]--;

		if((unit_index == face) && (fn_degree[t_direction_1] <= order_1) && (fn_degree[t_direction_2] <= order_2))  {
			if(position >= get_num_face_fns(order))
				ERROR("Internal in compute_face_indices().");

			indices[position++] = face_offset + face_fn_index;
		}
	}

	assert(position == get_num_face_fns(order));

	// indices has to be sorted for all faces in the same way
	// fastest changing : lobatto order, then legendre order and slowest v_direction
	int fac_ind_1, fn_deg_1[3], which_legendre_1, t_dir_1_1, t_dir_2_1, v_dir_1, x;
	int fac_ind_2, fn_deg_2[3], which_legendre_2, t_dir_1_2, t_dir_2_2, v_dir_2;

	for (int i = 0; i < get_num_face_fns(order); i++){
		for (int j = 0; j < get_num_face_fns(order) - 1; j++){
			fac_ind_1 = indices[j] - face_offset;
			fac_ind_2 = indices[j+1] - face_offset;
			check_face_fn_index(fac_ind_1, fn_deg_1, x, which_legendre_1, t_dir_1_1, t_dir_2_1, v_dir_1);
			check_face_fn_index(fac_ind_2, fn_deg_2, x, which_legendre_2, t_dir_1_2, t_dir_2_2, v_dir_2);
			int which_lobatto_1 = (t_dir_1_1 == which_legendre_1) ? t_dir_2_1 : t_dir_1_1;
			int which_lobatto_2 = (t_dir_1_2 == which_legendre_2) ? t_dir_2_2 : t_dir_1_2;
			if((v_dir_1 == v_dir_2) && (fn_deg_1[which_legendre_1] > fn_deg_2[which_legendre_2])){
				int hlp = indices[j];
				indices[j] = indices[j+1];
				indices[j+1] = hlp;
			}
		}
	}

	face_indices[face][wanted_ori][order.get_idx()] = indices;
}

//For a given "order" (i.e., three directional orders order_i), this returns the
//indices of all bubble functions whose directional polynomial degree in
//axial direction x_i is less or equal to "order_i"
void HCurlShapesetLobattoHex::compute_bubble_indices(order3_t order)
{
	int order_1 = order.x;
	int order_2 = order.y;
	int order_3 = order.z;
	int *indices = new int[get_num_bubble_fns(order)];
	int position = 0;

	int fn_degree[3];
	int unit_index = -1;
	int which_legendre, v_direction;
	//loop over all bubble functions
	for (int bubble_fn_index = 0; bubble_fn_index < max_bubble_fns * 8; bubble_fn_index++) {

		int ori = GET_ORI_FROM_INDEX(bubble_fn_index);
		int idx = GET_IDX_FROM_INDEX(bubble_fn_index);
		if(ori != 0)
			continue;

		check_bubble_fn_index(bubble_fn_index, fn_degree, which_legendre, v_direction);

		//this is according to the definition in the book
		for(int i = 0; i < 3; i++)
			if(i != which_legendre)
				fn_degree[i]--;

		if(fn_degree[0] <= order_1 && fn_degree[1] <= order_2 && fn_degree[2] <= order_3) {
			if(position >= max_bubble_fns)
				ERROR("Internal in compute_bubble_indices().");

			indices[position++] = bubble_offset + bubble_fn_index;
		}
	}
	assert(position == get_num_bubble_fns(order));
	bubble_indices[order.get_idx()] = indices;
}


order3_t HCurlShapesetLobattoHex::get_order(int index) const
{
	int ori = GET_ORI_FROM_INDEX(index);
	int idx = GET_IDX_FROM_INDEX(index);

	CHECK_INDEX(idx);

	int indices[3];
	int fn_type, unit_index, which_legendre, v_direction;
	int ord;

	if (idx >= 0) {
		// check_fn_index(idx, indices, fn_type, unit_index, which_legendre, v_direction);
//		ord = index_to_order[index];
	}
	else {
		assert(ced_key.exists(-1 - index));
		CEDKey key = ced_key[-1 - index];
		ord = key.order;
	}
		// face function is turned due to orientation
/*	if (ori >= 4) {
		ord = turn_hex_face_order(ord);
	}
*/	return ord;

}

/***************************************************
 *                  CONSTRUCTOR                    *
 ***************************************************/

HCurlShapesetLobattoHex::HCurlShapesetLobattoHex()
{
	mode = MODE_HEXAHEDRON;

	num_components = 3;

	//TODO think of those constants. What do they mean?
//	max_edge_order = MAX_ELEMENT_ORDER + 1;
//	max_face_order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1);
//	max_order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1);

//	max_fns_per_edge = get_num_edge_fns(MAX_ELEMENT_ORDER);
//	max_fns_per_face = get_num_face_fns(MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER));
//	max_bubble_fns = get_num_bubble_fns(MAKE_HEX_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER));

//	printf("max_fns_per_edge %d, face %d, bubble %d\n", max_fns_per_edge, max_fns_per_face, max_bubble_fns);

	face_offset = max_fns_per_edge * Hex::NUM_EDGES * 8; //orientations
	bubble_offset = max_fns_per_face * Hex::NUM_FACES * 8 + face_offset;
	max_index = max_bubble_fns * 8 + bubble_offset;

	// fn, dx, dy, dz will be calculated on-the-fly
	shape_table_deleg[FN]  = HCurlShapesetLobattoHex::calc_fn_value;
	shape_table_deleg[DX]  = HCurlShapesetLobattoHex::calc_dx_value;
	shape_table_deleg[DY]  = HCurlShapesetLobattoHex::calc_dy_value;
	shape_table_deleg[DZ]  = HCurlShapesetLobattoHex::calc_dz_value;
	shape_table_deleg[DXY] = NULL;
	shape_table_deleg[DXZ] = NULL;
	shape_table_deleg[DYZ] = NULL;

	// index to order mapping (determining orders for numerical quadrature)
	// original H1 code is below
	// THESE ARRAYS SHOULD BE DECLARED GLOBALLY (AND COMPRESSED?)
//	index_to_order = new int[max_index];

	for (int shape_fn_index = 0; shape_fn_index < max_index; shape_fn_index++) {
		int indices[3];
		int fn_type, unit_index, which_legendre, v_direction;
		int ord[3];

		check_fn_index(shape_fn_index, indices, fn_type, unit_index, which_legendre, v_direction);
		for(int i = 0; i < 3; i++)
			ord[i] = (i == which_legendre) ? legendre_order_1d[indices[i]] : lobatto_order_1d[indices[i]];


		int ori = GET_ORI_FROM_INDEX(shape_fn_index);
		order3_t order(ord[0], ord[1], ord[2]);

		// abych to nedelal dvakrat...
//		if(ori > 3)
//			turn_hex_face_order(order);

		//attempt to avoid underintegration
		//TODO this should not be necessary...
/*		int mo = ord[0];
		if(ord[1] > mo) mo = ord[1];
		if(ord[2] > mo) mo = ord[2];
		order = MAKE_HEX_ORDER(mo, mo, mo);
*/
//		index_to_order[shape_fn_index] = order;
	}

}

/***************************************************
*                   DESTRUCTOR                    *
***************************************************/

HCurlShapesetLobattoHex::~HCurlShapesetLobattoHex()
{
	for (int edge = 0; edge < Hex::NUM_EDGES; edge++)
		for (int ori = 0; ori < NUM_EDGE_ORIS; ori++)
			for (Word_t idx = edge_indices[edge][ori].first(); idx != INVALID_IDX; idx = edge_indices[edge][ori].next(idx))
				delete [] edge_indices[edge][ori][idx];

	for (int face = 0; face < Hex::NUM_FACES; face++)
		for (int ori = 0; ori < NUM_FACE_ORIS; ori++)
			for (Word_t idx = face_indices[face][ori].first(); idx != INVALID_IDX; idx = face_indices[face][ori].next(idx))
				delete [] face_indices[face][ori][idx];

	for (Word_t idx = bubble_indices.first(); idx != INVALID_IDX; idx = bubble_indices.next(idx))
		delete [] bubble_indices[idx];
}

///
///
///
/// --- CED specific stuff ---
///
///
///

// implementation has to be redefined from standart
// for face functions we use only part of facefunctions defined for given order
double HCurlShapesetLobattoHex::get_constrained_value(int n, int index, double x, double y, double z, int component) {
	assert(ced_key.exists(-1 - index));
	CEDKey key = ced_key[-1 - index];

	//special threatment of face functions
	if(key.type == CED_KEY_TYPE_FACE)
		return get_constrained_face_value(n, index, x, y, z, component);

	CEDComb *comb = get_ced_comb(key);
	int *idx = get_ced_indices(key);

	assert(comb != NULL);
	assert(idx != NULL);

	double sum = 0.0;
	for (int i = 0; i < comb->n; i++)
		sum += comb->coef[i] * get_val(n, idx[i], x, y, z, component);

	return sum;
}

double HCurlShapesetLobattoHex::get_constrained_face_value(int n, int index, double x, double y, double z, int component) {
	assert(ced_key.exists(-1 - index));
	CEDKey key = ced_key[-1 - index];

	CEDComb *comb = get_ced_comb(key);
	int *idx = get_ced_indices(key);

	assert(comb != NULL);
	assert(idx != NULL);

//	int horder = GET_QUAD_ORDER_1(key.order);
//	int vorder = GET_QUAD_ORDER_2(key.order);

//	int n_fns = get_num_face_fns(key.order);
	int ind[3], fn_type, unit_index, which_legendre, v_direction;
	int i, ind_fn;

	double sum = 0.0;
/*	for (i = 0, ind_fn = 0; ind_fn < n_fns; ind_fn++){
		check_fn_index(idx[ind_fn], ind, fn_type, unit_index, which_legendre, v_direction);
		int dir0 = RefHex::get_face_tangent_direction(unit_index, 0);
		int dir1 = RefHex::get_face_tangent_direction(unit_index, 1);
		if((ind[dir0] <= horder) && (ind[dir1] <= vorder) && (get_facefn_variant(idx[ind_fn]) == key.facefn_variant)){
			sum += comb->coef[i++] * get_val(n, idx[ind_fn], x, y, z, component);
		}
	}
	assert(i == comb->n);
*/

	return sum;
}


//
// constraints are calculated on egde 0
//
CEDComb *HCurlShapesetLobattoHex::calc_constrained_edge_combination(int ori, order1_t order, Part part) {
	// determine the interval of the edge
	double hi, lo;
	get_interval_part(part.part, lo, hi);

	int n = get_num_edge_fns(order);							// total number of functions on the edge
	int *fn_idx = get_edge_indices(0, ori, order);				// indices of all functions on the edge
	int cing_idx = fn_idx[n - 1];

//	int ind[3], fn_type, unit_index, which_legendre, v_direction;
//	check_fn_index(cing_idx, ind, fn_type, unit_index, which_legendre, v_direction);

//	double f_lo = get_value(FN, cing_idx, lo, -1.0, -1.0, v_direction);		// fn. values at endpoints of the part
//	double f_hi = get_value(FN, cing_idx, hi, -1.0, -1.0, v_direction);		// depends on the ref. domain

	// in Hcurl, vectors has to be transformed, when mapping to different domain.
	double trans = 2 / (hi - lo);

	double **a = new_matrix<double>(n, n);
	double *b = new double[n];
	for (int i = 0; i < n; i++) {
		// chebyshev point
		double p = cos((i + 1) * M_PI / (n+1));
		double r = (p + 1.0) * 0.5;
		double s = 1.0 - r;

		// matrix row
		for (int j = 0; j < n; j++)
			a[i][j] = trans * get_value(FN, fn_idx[j], p, -1.0, -1.0, 0);	// depends on the ref. domain

		// rhs
		b[i] = get_value(FN, cing_idx, lo*s + hi*r, -1.0, -1.0, 0);// - f_lo*s - f_hi*r;	// depends on the ref. domain

//		printf(" p = % lf, x = % lf\n", p, lo*s + hi*r);
	}

/*	printf("a = \n");
	for (int i = 0; i < n; i++) {
	printf("  ");
	for (int j = 0; j < n; j++) {
	printf(" % lf", a[i][j]);
}
	printf("\n");
}
	printf("\n");

	printf("b = ");
	for (int i = 0; i < n; i++)
	printf(" % lf ", b[i]);
	printf("\n");
*/
	// solve the system
	double d;
	int *iperm = new int[n];
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

/*	printf(" -- x =");
	for (int i = 0; i < n; i++)
	printf(" % lf", b[i]);
	printf("\n");
*/
	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);
}

//
// constraints are calculated on face 5
//
CEDComb *HCurlShapesetLobattoHex::calc_constrained_edge_face_combination(int ori, order2_t order, Part part, int facefn_variant) {
	// determine the interval of the face
	double hi, lo;
/*	get_interval_part(part.fpart, lo, hi);
	// determine the position
	double x0;
	get_edge_part(part.epart, x0);

	int horder, vorder;
	horder = GET_QUAD_ORDER_1(order);
	vorder = GET_QUAD_ORDER_2(order);

	const int *const_edge_ori = RefHex::get_face_edge_orientation(part.ori, ori);
	int edge_ori[] = {const_edge_ori[0], const_edge_ori[1]};

	int n, m;								// total number of functions on the edge
	int *edge_fn_idx[2];
	int cheb_pts;    //number of chebychev points

	if (part.ori == PART_ORI_VERT) {
//		printf("  PART_ORI_VERT\n");
		n = get_num_edge_fns(vorder);
		m = get_num_edge_fns(horder);											// total number of functions on the edge
		cheb_pts = vorder + 1;
		edge_fn_idx[0] = get_edge_indices(0, edge_ori[0], horder);			// indices of all functions on the edge
		edge_fn_idx[1] = get_edge_indices(0, edge_ori[1], vorder);			// indices of all functions on the edge
	}
	else{
//		printf("  PART_ORI_HORZ\n");
		n = get_num_edge_fns(horder);
		m = get_num_edge_fns(vorder);											// total number of functions on the edge
		cheb_pts = horder + 1;
		edge_fn_idx[0] = get_edge_indices(0, edge_ori[1], vorder);			// indices of all functions on the edge
		edge_fn_idx[1] = get_edge_indices(0, edge_ori[0], horder);			// indices of all functions on the edge
	}



//	printf("horder %d, vorder %d, part.ori %d, facefn_var %d, cheb bodu %d\n", horder, vorder, part.ori, facefn_variant, cheb_pts);

	//facefn_variant .. vectors parallel with first or second direction, taken for the face with ori <=3
	int oriented_facefn_variant = ori <=3 ? facefn_variant : 1-facefn_variant; // orientation taken into account

	if (((part.ori == PART_ORI_HORZ) && (facefn_variant == 1)) ||
			  ((part.ori == PART_ORI_VERT)  && (facefn_variant == 0))){
		double *b = new double[n];
		memset(b, 0, n * sizeof(double));
		return new CEDComb(n, b);
	 }


	int ind[3], fn_type, unit_index, which_legendre, v_direction;
	int face_n = get_num_face_fns(order);//MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER));	// total number of functions on the face
	int *face_fn_idx = get_face_indices(5, ori, order);//MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER));	// indices of all functions onthe face
	int cing_idx, i;


	for(i = 0; i < face_n; i++){
		cing_idx = face_fn_idx[i];
		check_fn_index(cing_idx, ind, fn_type, unit_index, which_legendre, v_direction);
		if((ind[0] == horder) && (ind[1] == vorder) && (get_facefn_variant(cing_idx) == facefn_variant))
			break;
	}
	assert(i < face_n);

//	printf("cing index : ");
//	dump_fn_index(cing_idx);

	check_fn_index(cing_idx, ind, fn_type, unit_index, which_legendre, v_direction);
	int which_lobatto = ind[(which_legendre + 1) % 3] > 1 ? (which_legendre + 1) % 3 : (which_legendre + 2) % 3;

//	int n = get_num_edge_fns(ind[which_legendre]);								// total number of functions on the edge
//	int *edge_fn_idx = get_edge_indices(0, edge_ori[0], ind[which_legendre]);		// indices of all functions on the edge

//	printf("calc comb which leg %d, ind (%d, %d, %d), n %d\n", which_legendre, ind[0], ind[1], ind[2], n);
//	printf("   orientations %d, %d \n", edge_ori[0], edge_ori[1]);

	assert(v_direction != 2);

	// in Hcurl, vectors has to be transformed, when mapping to different domain.
	double trans = 2 / (hi - lo);

	double **a = new_matrix<double>(n, n);
	double *b = new double[n];
	for (int i = 0; i < n; i++) {
		// chebyshev point
		double p = cos((i+1) * M_PI / (cheb_pts + 1));
		double r = (p + 1.0) * 0.5;
		double s = 1.0 - r;

		// matrix row
		for (int j = 0; j < n; j++){
			a[i][j] = trans * get_value(FN, edge_fn_idx[1][j], p, -1.0, -1.0, 0);		// depends on the ref. domain
//			printf("   m[%d][%d] = %lf\n", i, j, a[i][j]);
		}

		// rhs
		b[i] = get_value(FN, edge_fn_idx[1][n - 1], lo*s + hi*r, -1.0, -1.0, 0);	// depends on the ref. domain
//		printf("   rhs[%d] = %lf, v bode %lf, chebychevuv bod %lf, pocet ch. b. %d\n", i, b[i], lo*s + hi*r, p, cheb_pts);
	}

	// solve the system
	double d;
	int *iperm = new int[n];
//	printf("solving system of size %d\n", n);
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

//	for(int i = 0; i < n; i++)
//		printf("   sol[%d] = %lf\n", i, b[i]);

	//double c = get_value(FN, edge_fn_idx[0][m - 1], x0, -1.0, -1.0, 0);			// depends on the ref. domain

	if(edge_ori[1] == 1)
		x0 = -x0;
	double c = lobatto_fn_tab_1d[ind[which_lobatto]](x0);

	//	if(edge_ori[1] == 1)
	//		c = -c;

//	printf("nasobim koeficientem %lf, vzatym z bodu %lf\n", c, x0);

	for (int i = 0; i < n; i++)
		b[i] *= c;

	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);
*/
	return NULL;
}




//
// constraints are calculated on face 5
//
//
//            edge2
//  v_hi +-----------+
//       |           |
// edge3 |           | edge1
//       |           |
//  v_lo +-----------+
//     h_lo  edge0  h_hi
//

#define calc(xp, yp) \
	xr = (xp + 1.0) * 0.5; \
	xs = 1.0 - xr; \
	yr = (yp + 1.0) * 0.5;\
	ys = 1.0 - yr; \
	cing_val = 	get_value(FN, fn_idx[n - 1], h_lo * xs + h_hi * xr, v_lo * ys + v_hi * yr, 1.0, trans_v_direction);\
	sub_val = get_constrained_value(FN, ced_edge_idx[0], xp, yp, 1.0, trans_v_direction) \
				+get_constrained_value(FN, ced_edge_idx[1], xp, yp, 1.0, trans_v_direction) \
				+get_constrained_value(FN, ced_edge_idx[2], xp, yp, 1.0, trans_v_direction) \
				+get_constrained_value(FN, ced_edge_idx[3], xp, yp, 1.0, trans_v_direction) ;

CEDComb *HCurlShapesetLobattoHex::calc_constrained_face_combination(int ori, order2_t order, Part part, int facefn_variant) {
/*
	// determine the horizontal interval of the face
	double h_hi, h_lo;
	get_interval_part(part.horz, h_lo, h_hi);
	// determine the vertical interval of the face
	double v_hi, v_lo;
	get_interval_part(part.vert, v_lo, v_hi);

	int i_h_hi, i_h_lo;
	get_limit_part(part.horz, i_h_lo, i_h_hi);

	int i_v_hi, i_v_lo;
	get_limit_part(part.vert, i_v_lo, i_v_hi);

//	const int *edge_ori = RefHex::get_face_edge_orientation(ori, 1235); //TODO nove rozhrani
	int *fn_idx = get_face_indices(5, ori, order);    // indices of all functions on the face
	int face_n = get_num_face_fns(order);    // total number of functions on the face
	int n;       //number of functions with vectors in given direction. Will be used as size of the system

//      ori = 0;
	int horder = GET_QUAD_ORDER_1(order);
	int vorder = GET_QUAD_ORDER_2(order);

//	printf("\n\n\n********START***********\n order (%d, %d), varianta %d\n\n", horder, vorder, facefn_variant);

	n = facefn_variant ? horder * (vorder + 1) : (horder + 1) * vorder;

	int ind[3], fn_type, unit_index, which_legendre, v_direction;
	int i, cing_idx = -1, idx, used_n;
	int used_fn_idx[n];

	for(used_n = 0, i = 0; i < face_n; i++){
		idx = fn_idx[i];
//		printf("checking %d/%d ... %d\n", i+1, face_n, idx);
//		dump_fn_index(idx);
		check_fn_index(idx, ind, fn_type, unit_index, which_legendre, v_direction);
		if((ind[0] <= horder) && (ind[1] <= vorder) && (get_facefn_variant(idx) == facefn_variant)){
			used_fn_idx[used_n++] = idx;
//			printf("adding\n");
		}
		if((ind[0] == horder) && (ind[1] == vorder) && (get_facefn_variant(idx) == facefn_variant)){
			cing_idx = idx;
//			printf("cing\n");
		}
	}
	assert(cing_idx >= 0);
//	printf("order (%d, %d), used_n %d, n %d\n", horder, vorder, used_n, n);
	//assert(i < face_n);
	//assert(used_n == n);
	n = used_n;  //nektere funkce maji vetsi rad nez ta omezujici, pro zatim takhle...

//	printf("cing index : ");
//	dump_fn_index(cing_idx);

	check_fn_index(cing_idx, ind, fn_type, unit_index, which_legendre, v_direction);
	int which_lobatto = ind[(which_legendre + 1) % 3] > 1 ? (which_legendre + 1) % 3 : (which_legendre + 2) % 3;

	//prepare to substract edge parts
	Part part0, part1, part2, part3;

	int TRANS_HORZ = ori <= 3 ? PART_ORI_HORZ : PART_ORI_VERT;
	int TRANS_VERT = ori <= 3 ? PART_ORI_VERT : PART_ORI_HORZ;

	part0.fpart = part.horz; part0.epart = i_v_lo; part0.ori = TRANS_HORZ;
	part1.fpart = part.vert; part1.epart = i_h_hi; part1.ori = TRANS_VERT;
	part2.fpart = part.horz; part2.epart = i_v_hi; part2.ori = TRANS_HORZ;
	part3.fpart = part.vert; part3.epart = i_h_lo; part3.ori = TRANS_VERT;

	int ced_edge_idx[4];
	ced_edge_idx[0] = get_constrained_edge_face_index( 8, ori, order, part0, facefn_variant);
	ced_edge_idx[1] = get_constrained_edge_face_index( 9, ori, order, part1, facefn_variant);
	ced_edge_idx[2] = get_constrained_edge_face_index(10, ori, order, part2, facefn_variant);
	ced_edge_idx[3] = get_constrained_edge_face_index(11, ori, order, part3, facefn_variant);

	double **a = new_matrix<double>(n, n);
	double *b = new double[n];

	int h_num_cheb, v_num_cheb; //numbers of chebychev points in each direction
	if(facefn_variant == 0){
		h_num_cheb = horder + 1;
		v_num_cheb = vorder - 1;
	}
	else{
		h_num_cheb = horder - 1;
		v_num_cheb = vorder + 1;
	}

//	printf("h_cheb %d, v_cheb %d\n", h_num_cheb, v_num_cheb);

	// v_direction when orientation is taken ito account
	int trans_v_direction = v_direction;
	if(ori > 3){
		swapint(h_num_cheb, v_num_cheb);
		trans_v_direction = 1 - trans_v_direction;
	}
	assert(h_num_cheb * v_num_cheb == n);

	// in Hcurl, vectors has to be transformed, when mapping to different domain.
	double trans;
	if(trans_v_direction == 0)
		trans = 2 / (h_hi - h_lo);
	else
		trans = 2 / (v_hi - v_lo);


//	printf("hpart %d -> (%d, %d), vpart %d -> (%d, %d)\n", part.horz, i_h_lo, i_h_hi, part.vert, i_v_lo, i_v_hi);
//	const double test[5][2] = {{0.3, -1.}, {1., -0.8}, {1., 0.8}, {-0.43, 1.}, {-1., 0.12}};
//	double xr, xs, yr, ys, cing_val, sub_val;
//	for(int i = 0; i < 5; i++){
//		calc(test[i][0], test[i][1]);
//		if(fabs(cing_val-sub_val) > 0.0000000001)
//			printf("XXXXX order %d (%d, %d), ori %d,  pt (%lf, %lf) -> (%lf, %lf) XXXXXXXXXXXXX\n", order, horder, vorder, ori, test[i][0], test[i][1], cing_val, sub_val);
//	}
//	printf("\n");


	for (int row = 0; row < n; row++) {
		int i = row % h_num_cheb;
		int j = row / h_num_cheb;

		double hp = cos((i+1) * M_PI / (h_num_cheb + 1));
		double vp = cos((j+1) * M_PI / (v_num_cheb + 1));

		double hr = (hp + 1.0) * 0.5;
		double hs = 1.0 - hr;

		double vr = (vp + 1.0) * 0.5;
		double vs = 1.0 - vr;

		for (int k = 0; k < n; k++){
			a[row][k] = trans * get_value(FN, used_fn_idx[k], hp, vp, 1.0, trans_v_direction);	// depends on the ref. domain
//			dump_fn_index(used_fn_idx[k]);
//			printf("   a[%d][%d] = %lf, (i, j) == (%d, %d), cheb body (%lf, %lf) ... funkce %d\n", row, k, a[row][k], i, j, hp, vp, used_fn_idx[k]);
		}

		// rhs
		b[row] = get_value(FN, used_fn_idx[n - 1], h_lo * hs + h_hi * hr, v_lo * vs + v_hi * vr, 1.0, trans_v_direction)
			// subtract residual of edge functions
				- get_constrained_value(FN, ced_edge_idx[0], hp, vp, 1.0, trans_v_direction) //* vs
				- get_constrained_value(FN, ced_edge_idx[1], hp, vp, 1.0, trans_v_direction) //* hr
				- get_constrained_value(FN, ced_edge_idx[2], hp, vp, 1.0, trans_v_direction) //* vr
				- get_constrained_value(FN, ced_edge_idx[3], hp, vp, 1.0, trans_v_direction) ;//* hs;

//		printf("   b[%d] = %lf  @ (%lf, %lf), val pred odectenim %lf, vzato z dir %d ... funkce %d\n", row, b[row], h_lo * hs + h_hi * hr, v_lo * vs + v_hi * vr, get_value(FN, used_fn_idx[n - 1], h_lo * hs + h_hi * hr, v_lo * vs + v_hi * vr, 1.0, trans_v_direction), trans_v_direction, used_fn_idx[n - 1]);
	}

	// solve the system
	double d;
	int *iperm = new int[n];
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

	// cleanup
	delete [] iperm;
	delete [] a;

//	for(int i=0; i < n; i++)
//		printf("   sol[%d] = %lf\n", i, b[i]);

	return new CEDComb(n, b);
*/
	return NULL;
}


