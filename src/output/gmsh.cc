//
// gmshoutputengine.cc
//
//

#include "../config.h"
#include "gmsh.h"
#include "../refdomain.h"
#include "../quadstd.h"
#include "../common.h"
#include <stdio.h>
#include <errno.h>

// size of the buffer that is used for copying files
#define BUFLEN							8192
#define FORMAT							"%.17g"

///
#define AVGTV(a, b) { 0.5 * (tv[(a)].x + tv[(b)].x), 0.5 * (tv[(a)].y + tv[(b)].y), 0.5 * (tv[(a)].z + tv[(b)].z) }

namespace Gmsh {

//// OutputQuad //////////////////////////////////////////////////////////////////////////////

/// Common ancestor for output quadratures. Extends the interface of Quad3D
///
/// @ingroup visualization
class OutputQuad : public Quad3D {
public:
	virtual QuadPt3D *get_points(int order) {
		if (tables[order] == NULL) calculate_view_points(order);
		return tables[order];
	}

	virtual int get_num_points(int order) {
		if (np[order] == 0) calculate_view_points(order);
		return np[order];
	}

	virtual int *get_subdiv_modes(int order) {
		if (subdiv_modes[order] == NULL) calculate_view_points(order);
		return subdiv_modes[order];
	}

	virtual int get_subdiv_num(int order) {
		if (subdiv_num[order] == 0) calculate_view_points(order);
		return subdiv_num[order];
	}

	virtual QuadPt3D *get_face_points(int face, int order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	virtual void set_output_precision(int p) { output_precision = p; }

protected:
	int  *subdiv_num;
	int **subdiv_modes;
	int output_precision;

	virtual void calculate_view_points(int order) = 0;
	virtual void recursive_division(const Point3D *ref_vtcs, QuadPt3D *table, int levels, int &idx) = 0;
};


//// OutputQuadTetra //////////////////////////////////////////////////////////////////////////////

/// Quadrature for visualizing the solution on tetrahedron
///
/// @ingroup visualization
class OutputQuadTetra : public OutputQuad {
public:
	OutputQuadTetra();
	virtual ~OutputQuadTetra();

protected:
	virtual void calculate_view_points(int order);
	virtual void recursive_division(const Point3D *ref_vtcs, QuadPt3D *table, int levels, int &idx);
};

OutputQuadTetra::OutputQuadTetra() {
#ifdef WITH_TETRA
//	ref_domain = &ref_domain_std_tetra;

	max_order = MAX_QUAD_ORDER;

	tables = new QuadPt3D *[max_order + 1];
	MEM_CHECK(tables);
	memset(tables, 0, (max_order + 1) * sizeof(QuadPt3D *));

	np = new int[max_order + 1];
	MEM_CHECK(np);
	memset(np, 0, (max_order + 1) * sizeof(int));

	subdiv_num = new int[max_order + 1];
	MEM_CHECK(subdiv_num);
	memset(subdiv_num, 0, (max_order + 1) * sizeof(int));

	subdiv_modes = new int *[max_order + 1];
	MEM_CHECK(subdiv_modes);
	memset(subdiv_modes, 0, (max_order + 1) * sizeof(int *));

	output_precision = 1;
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

OutputQuadTetra::~OutputQuadTetra() {
#ifdef WITH_TETRA
	for (int i = 0; i < max_order + 1; i++)
		delete [] tables[i];
	delete [] tables;

	delete [] np;

	for (int i = 0; i < max_order + 1; i++)
		delete [] subdiv_modes[i];
	delete [] subdiv_modes;
	delete [] subdiv_num;
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}


void OutputQuadTetra::calculate_view_points(int order) {
#ifdef WITH_TETRA
	if (tables[order] == NULL) {
		// check if the order is greater than 0, because we are taking log(o)
		if (order == 0) order++;

		// there should be refinement levels enough to catch the properties of order 'order' functions on an element
		// choose a different formula if this does not behave well
		int levels = int(log(order) / log(2)) + 1;

		// each refinement level means that a tetrahedron is divided into 8 subtetrahedra
		// i.e., there are 8^levels resulting tetrahedra => (8^levels)*10 points
		np[order] = (1 << (3 * levels)) * 10;
		subdiv_num[order] = (1 << (3 * levels));

		subdiv_modes[order] = new int[subdiv_num[order]];
		MEM_CHECK(subdiv_modes[order]);
		// the new subelements are tetrahedra only
		for (int i = 0; i < subdiv_num[order]; i++)
			subdiv_modes[order][i] = MODE_TETRAHEDRON;

		// compute the table of points recursively
		tables[order] = new QuadPt3D[np[order]];
		MEM_CHECK(tables[order]);
		int idx = 0;
		const Point3D *ref_vtcs = RefTetra::get_vertices();
		recursive_division(ref_vtcs, tables[order], levels, idx);
	}
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

void OutputQuadTetra::recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx) {
#ifdef WITH_TETRA
	if (levels == 0) {
		// vertices
		for (int i = 0; i < Tetra::NUM_VERTICES; i++) {
			table[idx].x = tv[i].x;
			table[idx].y = tv[i].y;
			table[idx].z = tv[i].z;
			table[idx].w = 1.0;
			idx++;
		}
		// edges
		Point3D edge_pt[] = {
			AVGTV(0, 1), AVGTV(1, 2), AVGTV(0, 2), AVGTV(0, 3), AVGTV(2, 3), AVGTV(1, 3)
		};
		for (int i = 0; i < Tetra::NUM_EDGES; i++) {
			table[idx].x = edge_pt[i].x;
			table[idx].y = edge_pt[i].y;
			table[idx].z = edge_pt[i].z;
			table[idx].w = 1.0;
			idx++;
		}
	}
	else {
		Point3D div_vtcs[8][4] = {
			{ Point3D(tv[0].x, tv[0].y, tv[0].z), AVGTV(0,1), AVGTV(0,2), AVGTV(0,3) },
			{ AVGTV(0,1), Point3D(tv[1].x, tv[1].y, tv[1].z), AVGTV(1,2), AVGTV(1,3) },
			{ AVGTV(0,2), AVGTV(1,2), Point3D(tv[2].x, tv[2].y, tv[2].z), AVGTV(2,3) },
			{ AVGTV(0,3), AVGTV(1,3), AVGTV(2,3), Point3D(tv[3].x, tv[3].y, tv[3].z) },
			{ AVGTV(0,2), AVGTV(0,1), AVGTV(1,2), AVGTV(1,3) },
			{ AVGTV(0,3), AVGTV(0,1), AVGTV(0,2), AVGTV(1,3) },
			{ AVGTV(0,3), AVGTV(0,2), AVGTV(2,3), AVGTV(1,3) },
			{ AVGTV(2,3), AVGTV(0,2), AVGTV(1,2), AVGTV(1,3) }
		};

		for (int i = 0; i < 8; i++)
			recursive_division(div_vtcs[i], table, levels - 1, idx);
	}
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

//// OutputQuadHex ////////////////////////////////////////////////////////////////////////////////

int get_principal_order(Order3 order) {
	int o1 = GET_HEX_ORDER_1(order);
	int o2 = GET_HEX_ORDER_2(order);
	int o3 = GET_HEX_ORDER_3(order);
	return (o1 > o2 ? (o1 > o3 ? o1 : o3) : (o2 > o3 ? o2 : o3));
}

/// Quadrature for visualizing the solution on hexahedron
///
/// @ingroup visualization
class OutputQuadHex : public OutputQuad {
public:
	OutputQuadHex();
	virtual ~OutputQuadHex();

protected:
	virtual void calculate_view_points(int order);
	virtual void recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx);
};

OutputQuadHex::OutputQuadHex() {
#ifdef WITH_HEX
	max_order = MAX_QUAD_ORDER;

	tables = new QuadPt3D *[max_order + 1];
	MEM_CHECK(tables);
	memset(tables, 0, (max_order + 1) * sizeof(QuadPt3D *));

	np = new int[max_order + 1];
	MEM_CHECK(np);
	memset(np, 0, (max_order + 1) * sizeof(int));

	subdiv_num = new int[max_order + 1];
	MEM_CHECK(subdiv_num);
	memset(subdiv_num, 0, (max_order + 1) * sizeof(int));

	subdiv_modes = new int *[max_order + 1];
	MEM_CHECK(subdiv_modes);
	memset(subdiv_modes, 0, (max_order + 1) * sizeof(int *));

	output_precision = 1;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
#ifdef WITH_HEX
	for (int i = 0; i < max_order + 1; i++)
		delete [] tables[i];
	delete [] tables;

	delete [] np;

	for (int i = 0; i < max_order + 1; i++)
		delete [] subdiv_modes[i];
	delete [] subdiv_modes;
	delete [] subdiv_num;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void OutputQuadHex::calculate_view_points(int order) {
#ifdef WITH_HEX
	assert(order != 0);

	if (tables[order] == NULL) {
//		int o = get_principal_order(order);
//		int levels = int(log(o) / log(2)) + output_precision;
		int o = 1;
		int levels = 1;

		np[o] = (1 << (3 * levels)) * 27;
		subdiv_num[o] = (1 << (3 * levels));

		subdiv_modes[o] = new int[subdiv_num[o]];
		MEM_CHECK(subdiv_modes[o]);
		// the new subelements are hexahedra only
		for (int i = 0; i < subdiv_num[o]; i++)
			subdiv_modes[o][i] = MODE_HEXAHEDRON;

		// compute the table of points recursively
		tables[o] = new QuadPt3D[np[o]];
		MEM_CHECK(tables[o]);
		int idx = 0;
		const Point3D *ref_vtcs = RefHex::get_vertices();
		recursive_division(ref_vtcs, tables[o], levels, idx);
	}
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void OutputQuadHex::recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx) {
#ifdef WITH_HEX
	if (levels == 0) {
		// vertices
		for (int i = 0; i < Hex::NUM_VERTICES; i++) {
			table[idx].x = tv[i].x;
			table[idx].y = tv[i].y;
			table[idx].z = tv[i].z;
			table[idx].w = 1.0;
			idx++;
		}
		// edges
		Point3D edge_pts[] = {
			AVGTV(0, 1), AVGTV(0, 3), AVGTV(0, 4), AVGTV(1, 2), AVGTV(1, 5), AVGTV(2, 3),
			AVGTV(2, 6), AVGTV(3, 7), AVGTV(4, 5), AVGTV(4, 7), AVGTV(5, 6), AVGTV(6, 7)
		};
		for (int i = 0; i < Hex::NUM_EDGES; i++) {
			table[idx].x = edge_pts[i].x;
			table[idx].y = edge_pts[i].y;
			table[idx].z = edge_pts[i].z;
			table[idx].w = 1.0;
			idx++;
		}
		// faces
		Point3D face_pts[] = {
			AVGTV(0, 2), AVGTV(0, 5), AVGTV(0, 7), AVGTV(1, 6), AVGTV(3, 6), AVGTV(4, 6)
		};
		for (int i = 0; i < Hex::NUM_FACES; i++) {
			table[idx].x = face_pts[i].x;
			table[idx].y = face_pts[i].y;
			table[idx].z = face_pts[i].z;
			table[idx].w = 1.0;
			idx++;
		}
		// volume
		Point3D pt_vol = AVGTV(0, 6);
		table[idx].x = pt_vol.x;
		table[idx].y = pt_vol.y;
		table[idx].z = pt_vol.z;
		table[idx].w = 1.0;
		idx++;
	}
	else {
		Point3D div_vtcs[8][8] = {
			{ { tv[0].x, tv[0].y, tv[0].z }, AVGTV(0,1), AVGTV(0,2), AVGTV(0,3), AVGTV(0,4), AVGTV(0,5), AVGTV(0,6), AVGTV(0,7) },
			{ AVGTV(0,1), { tv[1].x, tv[1].y, tv[1].z }, AVGTV(1,2), AVGTV(1,3), AVGTV(1,4), AVGTV(1,5), AVGTV(1,6), AVGTV(1,7) },
			{ AVGTV(0,2), AVGTV(1,2), { tv[2].x, tv[2].y, tv[2].z }, AVGTV(2,3), AVGTV(2,4), AVGTV(2,5), AVGTV(2,6), AVGTV(2,7) },
			{ AVGTV(0,3), AVGTV(1,3), AVGTV(2,3), { tv[3].x, tv[3].y, tv[3].z }, AVGTV(3,4), AVGTV(3,5), AVGTV(3,6), AVGTV(3,7) },
			{ AVGTV(0,4), AVGTV(1,4), AVGTV(2,4), AVGTV(3,4), { tv[4].x, tv[4].y, tv[4].z }, AVGTV(4,5), AVGTV(4,6), AVGTV(4,7) },
			{ AVGTV(0,5), AVGTV(1,5), AVGTV(2,5), AVGTV(3,5), AVGTV(4,5), { tv[5].x, tv[5].y, tv[5].z }, AVGTV(5,6), AVGTV(5,7) },
			{ AVGTV(0,6), AVGTV(1,6), AVGTV(2,6), AVGTV(3,6), AVGTV(4,6), AVGTV(5,6), { tv[6].x, tv[6].y, tv[6].z }, AVGTV(6,7) },
			{ AVGTV(0,7), AVGTV(1,7), AVGTV(2,7), AVGTV(3,7), AVGTV(4,7), AVGTV(5,7), AVGTV(6,7), { tv[7].x, tv[7].y, tv[7].z } },
		};

		for (int i = 0; i < 8; i++)
			recursive_division(div_vtcs[i], table, levels - 1, idx);
	}
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

//// OutputQuadPrism /////////////////////////////////////////////////////////////////////////////

/// TODO: output quad for prisms

} // namespace

//
#ifdef WITH_TETRA
	static Gmsh::OutputQuadTetra output_quad_tetra;
	#define OUTPUT_QUAD_TETRA		&output_quad_tetra
#else
	#define OUTPUT_QUAD_TETRA		NULL
#endif

#ifdef WITH_HEX
	static Gmsh::OutputQuadHex output_quad_hex;
	#define OUTPUT_QUAD_HEX			&output_quad_hex
#else
	#define OUTPUT_QUAD_HEX			NULL
#endif

static Gmsh::OutputQuad *output_quad[] = { OUTPUT_QUAD_TETRA, OUTPUT_QUAD_HEX, NULL };

///

GmshOutputEngine::GmshOutputEngine(FILE *file) {
	this->out_file = file;
}

GmshOutputEngine::~GmshOutputEngine() {
}

void GmshOutputEngine::dump_scalars(int mode, Point3D *pts, double *values, int num_pts) {
	const char *id;
	switch (mode) {
		case MODE_TETRAHEDRON: id = "SS2"; break;
		case MODE_HEXAHEDRON:  id = "SH2"; break;
		case MODE_PRISM:       ERROR("Unsupported mode."); break;
		default:               ERROR("Invalid mode."); break;
	}

	// write id
	fprintf(this->out_file, "\t%s(", id);

	// write points
	for (int j = 0; j < num_pts; j++)
		fprintf(this->out_file, FORMAT ", " FORMAT ", " FORMAT "%s", pts[j].x, pts[j].y, pts[j].z, j == num_pts - 1 ? "" : ", ");
	fprintf(this->out_file, ") { ");
	// write values
	for (int j = 0; j < num_pts; j++)
		fprintf(this->out_file, FORMAT "%s", values[j], j == num_pts - 1 ? "" : ", ");
	// end the row
	fprintf(this->out_file, " };\n");
}

void GmshOutputEngine::out(MeshFunction *fn, const char *name, int item/* = FN_VAL*/) {
	// prepare
	fprintf(this->out_file, "View \"%s\" {\n", name);

//	Space *space = fn->get_space();
	Mesh *mesh = fn->get_mesh();

	Word_t idx;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();
		// FIXME: get order from the space
		int order;
		switch (mode) {
			// FIXME
			case MODE_TETRAHEDRON: order = 1; break;
			case MODE_HEXAHEDRON: order = 1; /*MAKE_HEX_ORDER(1, 1, 1);*/ break;
				break;
/*			case MODE_TETRAHEDRON: order = space->get_element_order(idx); break;
			case MODE_HEXAHEDRON: order = get_principal_order(space->get_element_order(idx)); break;
*/			case MODE_PRISM: EXIT(ERR_NOT_IMPLEMENTED); break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}

		Gmsh::OutputQuad *quad = output_quad[mode];
		int subdiv_num = quad->get_subdiv_num(order);
		int *modes = quad->get_subdiv_modes(order);

		fn->set_active_element(element);
		fn->set_quad(quad);

		RefMap *refmap = fn->get_refmap();
		double *phys_x = refmap->get_phys_x(order);
		double *phys_y = refmap->get_phys_y(order);
		double *phys_z = refmap->get_phys_z(order);

		fn->set_quad_order(order, item);

		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);
		scalar *val;
		val = fn->get_values(a, b);
		int pt_idx = 0;
		// iterate through subelements and output them
		for (int i = 0; i < subdiv_num; i++) {
			int num_pts;
			switch (mode) {
				case MODE_TETRAHEDRON: num_pts = 10; break; // quadratic TETRA (see gmsh documentation)
				case MODE_HEXAHEDRON:  num_pts = 27; break; // quadratic HEX (see gmsh documentation)
				case MODE_PRISM: EXIT(ERR_NOT_IMPLEMENTED); break;
				default: EXIT(ERR_UNKNOWN_MODE); break;
			}

			// small buffers to hold values for one subelement
			Point3D phys_pt[num_pts];
			scalar *p_val = val + pt_idx;

			for (int j = 0; j < num_pts; j++, pt_idx++) {
				// physical coordinates of subelement
				phys_pt[j].x = phys_x[pt_idx];
				phys_pt[j].y = phys_y[pt_idx];
				phys_pt[j].z = phys_z[pt_idx];
			}

			dump_scalars(mode, phys_pt, p_val, num_pts);
		}
	}

	// finalize
	fprintf(this->out_file, "};\n");
}

void GmshOutputEngine::out(Mesh *mesh) {
	// see Gmsh documentation on details (http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html)
	Word_t idx;

	// header
	fprintf(this->out_file, "$MeshFormat\n");
	fprintf(this->out_file, "%.1lf %d %d\n", 2.0, 0, sizeof(double));
	fprintf(this->out_file, "$EndMeshFormat\n");

	// vertices
	fprintf(this->out_file, "$Nodes\n");
	fprintf(this->out_file, "%d\n", mesh->vertices.count());
	FOR_ALL_VERTICES(idx, mesh) {
		Vertex *v = mesh->vertices[idx];
		fprintf(this->out_file, "%d %lf %lf %lf\n", idx, v->x, v->y, v->z);
	}
	fprintf(this->out_file, "$EndNodes\n");

	// elements
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", mesh->get_num_active_elements());
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[element->get_num_of_vertices()];
		element->get_vertices(vtcs);

		switch (element->get_mode()) {
			case MODE_TETRAHEDRON:
				fprintf(this->out_file, "%ld 4 0 %d %d %d %d\n",
					element->id, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
				break;

			case MODE_HEXAHEDRON:
				fprintf(this->out_file, "%ld 5 0 %d %d %d %d %d %d %d %d\n",
					element->id, vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
				break;

			case MODE_PRISM:
				EXIT(ERR_NOT_IMPLEMENTED);
				break;

			default:
				EXIT(ERR_UNKNOWN_MODE);
				break;
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// edges
	// TODO: do not include edges twice or more
	// FIXME: Hex-specific
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);
	FOR_ALL_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[Edge::NUM_VERTICES];
		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			element->get_edge_vertices(iedge, vtcs);
			fprintf(this->out_file, "%ld 1 0 %d %d\n", mesh->get_edge_id(vtcs[0], vtcs[1]), vtcs[0], vtcs[1]);
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// faces
	// TODO: do not include faces twice
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", mesh->get_num_active_elements() * Hex::NUM_FACES);
	FOR_ALL_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[Quad::NUM_VERTICES];				// FIXME: HEX-specific
		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			element->get_face_vertices(iface, vtcs);
			fprintf(this->out_file, "%ld 3 0 %d %d %d %d\n", mesh->get_facet_id(element, iface), vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
		}
	}
	fprintf(this->out_file, "$EndElements\n");
}

void GmshOutputEngine::out_orders(Space *space, const char *name) {
	Mesh *mesh = space->get_mesh();

	// prepare
	fprintf(this->out_file, "$MeshFormat\n");
	fprintf(this->out_file, "%.1lf %d %d\n", 2.0, 0, sizeof(double));
	fprintf(this->out_file, "$EndMeshFormat\n");

	// HEX specific

	// nodes
	fprintf(this->out_file, "$Nodes\n");
	fprintf(this->out_file, "%d\n", mesh->get_num_active_elements() * 15);
	int id = 1;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int nv = Hex::NUM_VERTICES;
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			fprintf(this->out_file, "%d %lf %lf %lf\n", id++, v->x, v->y, v->z);
		}

		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			Word_t fvtcs[Quad::NUM_VERTICES];
			element->get_face_vertices(iface, fvtcs);
			Vertex *v[4] = { mesh->vertices[fvtcs[0]], mesh->vertices[fvtcs[1]], mesh->vertices[fvtcs[2]], mesh->vertices[fvtcs[3]] };
			Vertex fcenter((v[0]->x + v[2]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[2]->z) / 2.0);
			fprintf(this->out_file, "%d %lf %lf %lf\n", id++, fcenter.x, fcenter.y, fcenter.z);
		}

		Vertex *v[4] = { mesh->vertices[vtcs[0]], mesh->vertices[vtcs[1]], mesh->vertices[vtcs[3]], mesh->vertices[vtcs[4]] };
		Vertex center((v[0]->x + v[1]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[3]->z) / 2.0);
		fprintf(this->out_file, "%d %lf %lf %lf\n", id++, center.x, center.y, center.z);
	}
	fprintf(this->out_file, "$EndNodes\n");

	// elements
	fprintf(this->out_file, "$Elements\n");
	// FIXME: hex specific
	fprintf(this->out_file, "%d\n", mesh->get_num_active_elements() + (mesh->get_num_active_elements() * Hex::NUM_EDGES));
	id = 1;
	// trick: put the elements first so that they will be visible in gmsh
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[element->get_num_of_vertices()];
		element->get_vertices(vtcs);

		fprintf(this->out_file, "%ld 5 0 %d %d %d %d %d %d %d %d\n",
			id++, vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
	}
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[element->get_num_of_vertices()];
		element->get_vertices(vtcs);

		// orders
		int pyrel[][4] = {
			{  1,  9,  5, 10 },
			{  2, 11,  6,  9 },
			{  7,  8,  3, 11 },
			{  0, 10,  4,  8 },
			{  0, 10,  1, 12 },
			{  1,  9,  2, 12 },
			{  2, 11,  3, 12 },
			{  3, 12,  0,  8 },
			{  4, 10,  5, 13 },
			{  5,  9,  6, 13 },
			{  6, 11,  7, 13 },
			{  7, 13,  4,  8 }
		};

		for (int i = 0; i < Hex::NUM_EDGES; i++) {
			int base = (idx - 1) * 15;
			int v[4] = { base + pyrel[i][0] + 1, base + pyrel[i][1] + 1, base + pyrel[i][2] + 1, base + pyrel[i][3] + 1};
			fprintf(this->out_file, "%ld 3 0 %d %d %d %d\n", id++, v[0], v[1], v[2], v[3]);
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// node data
	fprintf(this->out_file, "$ElementNodeData\n");
	fprintf(this->out_file, "%d\n", 1);
	fprintf(this->out_file, "\"%s\"\n", name);
	fprintf(this->out_file, "%d\n", 0);
	fprintf(this->out_file, "%d\n", 3);
	fprintf(this->out_file, "0\n");		// time step (not used, but has to be there)
	fprintf(this->out_file, "1\n");		// 1 value per node
	fprintf(this->out_file, "%d\n", (mesh->get_num_active_elements() * Hex::NUM_EDGES));

	id = mesh->get_num_active_elements() + 1;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();
		// get order from the space
		int order = space->get_element_order(idx);

		for (int i = 0; i < Hex::NUM_EDGES; i++) {
			int dord;
			if (i == 0 || i == 1 || i == 2 || i == 3) dord = GET_HEX_ORDER_3(order);
			else if (i == 4 || i == 6 || i == 8 || i == 10) dord = GET_HEX_ORDER_1(order);
			else if (i == 5 || i == 7 || i == 9 || i == 11) dord = GET_HEX_ORDER_2(order);
			else assert(false);
			fprintf(this->out_file, "%ld 4 %d %d %d %d\n", id++, dord, dord, dord, dord);
		}

	}

	fprintf(this->out_file, "$EndElementNodeData\n");
}
