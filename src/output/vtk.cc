// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2007 - 2009 Pavel Kus <pavel.kus@gmail.com>
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
// vtkoutputengine.cc
//
//

#include "../h3dconfig.h"
#include "vtk.h"
#include "../refdomain.h"
#include "../quadstd.h"
#include "../common.h"

#include <stdio.h>
#include <errno.h>
#include <common/utils.h>

// size of the buffer that is used for copying files
#define BUFLEN							8192
#define FORMAT							"%.17g"

#define VTK_TETRA						10
#define VTK_HEXAHEDRON					12
#define VTK_WEDGE						13
#define VTK_PYRAMID                    14

namespace Vtk {

//// OutputQuad //////////////////////////////////////////////////////////////////////////////

/// Common ancestor for output quadratures. Extends the interface of Quad3D
///
/// @ingroup visualization
class OutputQuad : public Quad3D {
public:
	virtual QuadPt3D *get_points(order3_t order) {
		if (!tables.exists(order.get_idx())) calculate_view_points(order);
		return tables[order.get_idx()];
	}

	virtual int get_num_points(order3_t order) {
		if (!np.exists(order.get_idx())) calculate_view_points(order);
		return np[order.get_idx()];
	}

protected:
	virtual void calculate_view_points(order3_t order) = 0;
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
	virtual void calculate_view_points(order3_t order);
};

OutputQuadTetra::OutputQuadTetra() {
}

OutputQuadTetra::~OutputQuadTetra() {
}

void OutputQuadTetra::calculate_view_points(order3_t order) {
}

/// Quadrature for visualizing the solution on hexahedron
///
/// @ingroup visualization
class OutputQuadHex : public OutputQuad {
public:
	OutputQuadHex();
	virtual ~OutputQuadHex();

protected:
	virtual void calculate_view_points(order3_t order);
};

OutputQuadHex::OutputQuadHex() {
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
#ifdef WITH_HEX
	for (Word_t i = tables.first(); i != INVALID_IDX; i = tables.next(i))
		delete[] tables[i];
#endif
}

void OutputQuadHex::calculate_view_points(order3_t order) {
#ifdef WITH_HEX
	int o = order.get_idx();
	np[o] = (order.x + 1) * (order.y + 1) * (order.z + 1);
	tables[o] = new QuadPt3D[np[o]];
	double step_x, step_y, step_z;
	step_x = 2.0 / order.x;
	step_y = 2.0 / order.y;
	step_z = 2.0 / order.z;

	int n = 0;
	for (int k = 0; k < (int) order.x + 1; k++) {
		for (int l = 0; l < (int) order.y + 1; l++) {
			for (int m = 0; m < (int) order.z + 1; m++, n++) {
				assert(n < np[o]);
				tables[o][n].x = (step_x * k) - 1;
				tables[o][n].y = (step_y * l) - 1;
				tables[o][n].z = (step_z * m) - 1;
				tables[o][n].w = 1.0;	// not used
			}
		}
	}
#endif
}

//// OutputQuadPrism /////////////////////////////////////////////////////////////////////////////

/// TODO: output quad for prisms

} // namespace

//
#ifdef WITH_TETRA
static Vtk::OutputQuadTetra output_quad_tetra;
#define OUTPUT_QUAD_TETRA		&output_quad_tetra
#else
#define OUTPUT_QUAD_TETRA		NULL
#endif

#ifdef WITH_HEX
static Vtk::OutputQuadHex output_quad_hex;
#define OUTPUT_QUAD_HEX			&output_quad_hex
#else
#define OUTPUT_QUAD_HEX			NULL
#endif

static Vtk::OutputQuad *output_quad[] = { OUTPUT_QUAD_TETRA, OUTPUT_QUAD_HEX, NULL };

VtkOutputEngine::VtkOutputEngine(FILE *file) {
	this->out_file = file;
	this->has_points = false;
}

VtkOutputEngine::~VtkOutputEngine() {
}

order3_t VtkOutputEngine::get_order(int mode) {
	// FIXME: get order from the space and set sufficient division
	order3_t order(3, 3, 3);
	switch (mode) {
		case MODE_HEXAHEDRON:
			break;
		case MODE_TETRAHEDRON:
		case MODE_PRISM:
			EXIT(ERR_NOT_IMPLEMENTED); break;
		default:
			EXIT(ERR_UNKNOWN_MODE); break;
	}

	return order;
}

void VtkOutputEngine::dump_points(MeshFunction *fn) {
	Array<Point3D *> vertices;
	Array<int *> cells[3]; // 3 types of elements

	int type[] = { MODE_TETRAHEDRON, MODE_HEXAHEDRON, MODE_PRISM };
	int type_id[] = { VTK_TETRA, VTK_HEXAHEDRON, VTK_WEDGE };
	Mesh *mesh = fn->get_mesh();

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();

		Vtk::OutputQuad *quad = output_quad[mode];
		fn->set_active_element(element);
		fn->set_quad(quad);

		order3_t order = get_order(mode);

		// get coordinates of all points
		RefMap *refmap = fn->get_refmap();
		double *phys_x = refmap->get_phys_x(order);
		double *phys_y = refmap->get_phys_y(order);
		double *phys_z = refmap->get_phys_z(order);

		// insert points int the vertex array
		int np = quad->get_num_points(order);
		for (int i = 0; i < np; i++) {
			Point3D *pt = new Point3D;
			pt->x = phys_x[i];
			pt->y = phys_y[i];
			pt->z = phys_z[i];
			vertices.add(pt);
		}

		// insert cells
		switch (mode) {
			case MODE_HEXAHEDRON:
				for (unsigned int i = 0; i < order.x; i++) {
					for (unsigned int j = 0; j < order.y; j++) {
						for (unsigned int o = 0; o < order.z; o++) {
							int *cell = new int [Hex::NUM_VERTICES];
							cell[0] = (order.z + 1) * (i * (order.y + 1) + j) + o;
							cell[1] = cell[0] + ((order.y + 1) * (order.z + 1));
							cell[2] = cell[1] + (order.z + 1);
							cell[3] = cell[0] + (order.z + 1);
							cell[4] = cell[0] + 1;
							cell[5] = cell[1] + 1;
							cell[6] = cell[2] + 1;
							cell[7] = cell[3] + 1;

							cells[1].add(cell);
						}
					}
				}
				break;

			case MODE_PRISM:
			case MODE_TETRAHEDRON:
				EXIT(ERR_NOT_IMPLEMENTED); break;

			default: EXIT(ERR_UNKNOWN_MODE); break;
		} // switch
	}

	// DUMP ////////////////////////////////////////////////////////////////////////////////////////

	// start the output
	fprintf(this->out_file, "# vtk DataFile Version 2.0\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "ASCII\n");

	// dataset
	fprintf(this->out_file, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "POINTS %ld %s\n", vertices.count(), "float");
	for (Word_t i = vertices.first(); i != INVALID_IDX; i = vertices.next(i)) {
		Point3D *pt = vertices[i];
		fprintf(this->out_file, "%e %e %e\n", pt->x, pt->y, pt->z);
	}
	fprintf(this->out_file, "\n");

	fprintf(this->out_file, "CELLS %ld %ld\n", cells[0].count() + cells[1].count() + cells[2].count(),
		(Tetra::NUM_VERTICES + 1) * cells[0].count() +
		(Hex::NUM_VERTICES + 1) * cells[1].count() +
		(Prism::NUM_VERTICES + 1) * cells[2].count());
	for (unsigned int i = 0; i < countof(type); i++) { // 3 types of elements
		int pt_cnt[] = { Tetra::NUM_VERTICES, Hex::NUM_VERTICES, Prism::NUM_VERTICES };
		for (Word_t j = cells[type[i]].first(); j != INVALID_IDX; j = cells[type[i]].next(j)) {
			fprintf(this->out_file, "%d", pt_cnt[i]);
			for (int k = 0; k < pt_cnt[i]; k++)
				fprintf(this->out_file, " %d", cells[type[i]][j][k]);
			fprintf(this->out_file, "\n");
		}
	}
	fprintf(this->out_file, "\n");

	fprintf(this->out_file, "CELL_TYPES %ld\n", cells[0].count() + cells[1].count() + cells[2].count());
	for (unsigned int i = 0; i < countof(type); i++) { // 3 types of elements
		for (Word_t j = 0; j < cells[type[i]].count(); j++)
			fprintf(this->out_file, "%d\n", type_id[i]);
	}

	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "POINT_DATA %ld\n", vertices.count());

	// free allocated memory
	for (Word_t i = vertices.first(); i != INVALID_IDX; i = vertices.next(i))
		delete vertices[i];

	for (int i = 0; i < 3; i++) // 3 types of elements
		for (Word_t j = cells[i].first(); j != INVALID_IDX; j = cells[i].next(j))
			delete cells[i][j];
}


void VtkOutputEngine::out(MeshFunction *fn, const char *name, int item/* = FN_VAL_0*/) {
	if (!has_points) {
		// calculate points where we will evaluate the function 'fn'
		dump_points(fn);
		this->has_points = true;
	}

	assert(fn->get_num_components() == 1 || fn->get_num_components() == 3);

	int comp[COMPONENTS];		// components to output
	int nc;						// number of components to output
	int b = 0;
	if (fn->get_num_components() == COMPONENTS) {
		int a = 0;
		if ((item & FN_COMPONENT_0) && (item & FN_COMPONENT_1) && (item & FN_COMPONENT_2)) {
			mask_to_comp_val(item, a, b);
			for (int i = 0; i < COMPONENTS; i++) comp[i] = i;
			nc = 3;
		}
		else if ((item & FN_COMPONENT_0) > 0) {
			mask_to_comp_val(item & FN_COMPONENT_0, a, b);
			comp[0] = 0;
			nc = 1;
		}
		else if ((item & FN_COMPONENT_1) > 0) {
			mask_to_comp_val(item & FN_COMPONENT_1, a, b);
			comp[0] = 1;
			nc = 1;
		}
		else if ((item & FN_COMPONENT_2) > 0) {
			mask_to_comp_val(item & FN_COMPONENT_2, a, b);
			comp[0] = 2;
			nc = 1;
		}
		else {
			fprintf(this->out_file, "Unable to satisfy your request\n");
			return;					// Do not know what user wants
		}
	}
	else if (fn->get_num_components() == 1) {
		mask_to_comp_val(item & FN_COMPONENT_0, comp[0], b);
		nc = 1;
	}
	else {
		fprintf(this->out_file, "Unable to satisfy your request\n");
		return;					// Do not know what user wants
	}

	Mesh *mesh = fn->get_mesh();

	// header
	fprintf(this->out_file, "\n");
	if (nc == 1) {
		fprintf(this->out_file, "SCALARS %s %s %d\n", name, "float", 1);
		fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");
	}
	else if (nc == 3) {
		fprintf(this->out_file, "VECTORS %s %s\n", name, "float");
	}
	else assert(false);

	// values
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();

		order3_t order = get_order(mode);
		Vtk::OutputQuad *quad = output_quad[mode];
		fn->set_active_element(element);
		fn->set_quad(quad);
		fn->set_quad_order(ELEM_QORDER(order), item);
		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);
		scalar *val[3];
		for (int ic = 0; ic < nc; ic++)
			val[ic] = fn->get_values(ic, b);

		int np = quad->get_num_points(order);
		for (int i = 0; i < np; i++) {
			if (nc == 1) {				// scalar
#ifndef COMPLEX
				fprintf(this->out_file, "%e\n", val[0][i]);
#else
				assert(fabs(IMAG(val[0][i])) < 1e-12);
				fprintf(this->out_file, "%e\n", REAL(val[0][i]));
#endif
			}
			else if (nc == 3) {			// vector
#ifndef COMPLEX
				fprintf(this->out_file, "%e %e %e\n", val[0][i], val[1][i], val[2][i]);
#else
				assert(fabs(IMAG(val[0][i])) < 1e-12);
				assert(fabs(IMAG(val[1][i])) < 1e-12);
				assert(fabs(IMAG(val[2][i])) < 1e-12);
				fprintf(this->out_file, "%e %e %e\n", REAL(val[0][i]), REAL(val[1][i]), REAL(val[2][i]));
#endif
			}
		}
	}

}

void VtkOutputEngine::out(Mesh *mesh) {
	// Not implemented
	ERROR(ERR_NOT_IMPLEMENTED);
}

void VtkOutputEngine::out_orders(Space *space, const char *name) {
	Mesh *mesh = space->get_mesh();

	// FIXME: this function is hex-specific

	// local indices of subelements that build up the element
	int pyrel[][4] = {
		{  0, 10,  1, 12 },
		{  1,  9,  2, 12 },
		{  2, 11,  3, 12 },
		{  3, 12,  0,  8 },

		{  0, 10,  4,  8 },
		{  1,  9,  5, 10 },
		{  2, 11,  6,  9 },
		{  3, 11,  7,  8 },

		{  4, 10,  5, 13 },
		{  5,  9,  6, 13 },
		{  6, 11,  7, 13 },
		{  7, 13,  4,  8 }
	};

	// start the output
	fprintf(this->out_file, "# vtk DataFile Version 2.0\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "ASCII\n");

	// dataset
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(this->out_file, "POINTS %ld %s\n", mesh->get_num_active_elements() * 14, "float");
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int nv = Hex::NUM_VERTICES;
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			fprintf(this->out_file, "%e %e %e\n", v->x, v->y, v->z);
		}
		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			Word_t fvtcs[Quad::NUM_VERTICES];
			element->get_face_vertices(iface, fvtcs);
			Vertex *fv[4] = { mesh->vertices[fvtcs[0]], mesh->vertices[fvtcs[1]], mesh->vertices[fvtcs[2]], mesh->vertices[fvtcs[3]] };
			Vertex fctr((fv[0]->x + fv[2]->x) / 2.0, (fv[0]->y + fv[2]->y) / 2.0, (fv[0]->z + fv[2]->z) / 2.0);
			fprintf(this->out_file, "%e %e %e\n", fctr.x, fctr.y, fctr.z);
		}
	}
	fprintf(this->out_file, "\n");

	fprintf(this->out_file, "CELLS %ld %ld\n",
	        mesh->get_num_active_elements() * Hex::NUM_EDGES,
	        mesh->get_num_active_elements() * Hex::NUM_EDGES * 5
	);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[element->get_num_of_vertices()];
		element->get_vertices(vtcs);

		// orders
		for (int i = 0; i < Hex::NUM_EDGES; i++) {
			Word_t base = (idx - 1) * 14;
			Word_t v[4] = { base + pyrel[i][0], base + pyrel[i][1], base + pyrel[i][2], base + pyrel[i][3] };
			fprintf(this->out_file, "%d %ld %ld %ld %ld\n", 4, v[0], v[1], v[2], v[3]);
		}
	}

	fprintf(this->out_file, "\n");

	fprintf(this->out_file, "CELL_TYPES %ld\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);
	for (unsigned int i = 0; i < mesh->get_num_active_elements() * Hex::NUM_EDGES; i++) {
		fprintf(this->out_file, "%d\n", VTK_TETRA);
	}

	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "CELL_DATA %ld\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);

	// header
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "SCALARS %s %s %d\n", name, "float", 1);
	fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");

	// values
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();
		assert(mode == MODE_HEXAHEDRON);			// HEX-specific
		// get order from the space
		order3_t order = space->get_element_order(idx);

		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			int dord;
			if (iedge == 0 || iedge == 2 || iedge == 8 || iedge == 10) dord = order.x;
			else if (iedge == 1 || iedge == 3 || iedge == 9 || iedge == 11) dord = order.y;
			else if (iedge == 4 || iedge == 5 || iedge == 6 || iedge == 7) dord = order.z;
			else assert(false);
			fprintf(this->out_file, "%d\n", dord);
		}
	}
}
