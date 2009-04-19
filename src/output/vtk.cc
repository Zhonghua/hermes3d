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
#include <common/callstack.h>

// size of the buffer that is used for copying files
#define BUFLEN							8192
#define FORMAT							"%.17g"

#define VTK_TRIANGLE					5
#define VTK_QUAD						9

#define VTK_TETRA						10
#define VTK_HEXAHEDRON					12
#define VTK_WEDGE						13
#define VTK_PYRAMID						14

namespace Vtk {

//// OutputQuad //////////////////////////////////////////////////////////////////////////////

/// Common ancestor for output quadratures. Extends the interface of Quad3D
///
/// @ingroup visualization
class OutputQuad : public Quad3D {
public:
	virtual QuadPt3D *get_points(order3_t order) {
		_F_
		if (!tables.exists(order.get_idx())) calculate_view_points(order);
		return tables[order.get_idx()];
	}

	virtual int get_num_points(order3_t order) {
		_F_
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
	_F_
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	for (Word_t i = tables.first(); i != INVALID_IDX; i = tables.next(i))
		delete[] tables[i];
#endif
}

void OutputQuadHex::calculate_view_points(order3_t order) {
	_F_
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
	_F_
	this->out_file = file;
	this->has_points = false;
}

VtkOutputEngine::~VtkOutputEngine() {
	_F_
}

order3_t VtkOutputEngine::get_order(int mode) {
	_F_
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
	_F_
	Array<Point3D *> vertices;
	Array<int *> cells[3]; // 3 types of elements

	int type[] = { MODE_TETRAHEDRON, MODE_HEXAHEDRON, MODE_PRISM };
	int type_id[] = { VTK_TETRA, VTK_HEXAHEDRON, VTK_WEDGE };
	Mesh *mesh = fn->get_mesh();

	Word_t base = 0;
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

		// insert points in the vertex array
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
							cell[0] = base + (order.z + 1) * (i * (order.y + 1) + j) + o;
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
				base += (order.x + 1) * (order.y + 1) * (order.z + 1);
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
	_F_
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
		scalar *val[COMPONENTS];
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
	_F_

	// count the number of elements
	int ne = 0;
	int ne_sz = 0;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		switch (element->get_mode()) {
			case MODE_HEXAHEDRON: ne_sz += Hex::NUM_VERTICES + 1; break;
			case MODE_TETRAHEDRON: ne_sz += Tetra::NUM_VERTICES + 1; break;
			default: ERROR(ERR_NOT_IMPLEMENTED); break;
		}
		ne++;
	}

	fprintf(this->out_file, "# vtk DataFile Version 2.0\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "ASCII\n");

	// dataset
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "DATASET UNSTRUCTURED_GRID\n");

	// points
	fprintf(this->out_file, "POINTS %ld %s\n", mesh->vertices.count(), "float");
	for (Word_t i = mesh->vertices.first(); i != INVALID_IDX; i = mesh->vertices.next(i)) {
		Vertex *v = mesh->vertices[i];
		fprintf(this->out_file, "%e %e %e\n", v->x, v->y, v->z);
	}
	fprintf(this->out_file, "\n");

	// cells
	fprintf(this->out_file, "CELLS %d %d\n", ne, ne_sz);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		Word_t vtcs[Hex::NUM_VERTICES]; // FIXME: HEX-specific
		element->get_vertices(vtcs);

		switch (element->get_mode()) {
			case MODE_HEXAHEDRON:
				fprintf(this->out_file, "%d %ld %ld %ld %ld %ld %ld %ld %ld\n",
						Hex::NUM_VERTICES,
						vtcs[0] - 1, vtcs[1] - 1, vtcs[2] - 1, vtcs[3] - 1,
						vtcs[4] - 1, vtcs[5] - 1, vtcs[6] - 1, vtcs[7] - 1);
				break;

			case MODE_TETRAHEDRON:
				fprintf(this->out_file, "%d %ld %ld %ld %ld\n",
						Tetra::NUM_VERTICES,
						vtcs[0] - 1, vtcs[1] - 1, vtcs[2] - 1, vtcs[3] - 1);
				break;

			default:
				ERROR(ERR_NOT_IMPLEMENTED);
				break;
		}
	}
	fprintf(this->out_file, "\n");

	// cell types
	fprintf(this->out_file, "CELL_TYPES %d\n", ne);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		switch (element->get_mode()) {
			case MODE_HEXAHEDRON: fprintf(this->out_file, "%d\n", VTK_HEXAHEDRON); break;
			case MODE_TETRAHEDRON: fprintf(this->out_file, "%d\n", VTK_TETRA); break;
			default: ERROR(ERR_NOT_IMPLEMENTED); break;
		}
	}
	fprintf(this->out_file, "\n");

	// cell data
	fprintf(this->out_file, "CELL_DATA %d\n", ne);
	fprintf(this->out_file, "SCALARS %s %s %d\n", "mesh", "float", 1);
	fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");

	// values
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		fprintf(this->out_file, "%lf\n", 0.0);
	}
}


void VtkOutputEngine::out_bc(Mesh *mesh, const char *name) {
	_F_

	// count outer facets
	int fc = 0;
	int fc_sz = 0;
	FOR_ALL_FACETS(idx, mesh) {
		Facet *facet = mesh->facets[idx];
		if (facet->type == Facet::OUTER) {
			fc++;
			switch (facet->mode) {
				case MODE_TRIANGLE: fc_sz += Tri::NUM_VERTICES + 1; break;
				case MODE_QUAD: fc_sz += Quad::NUM_VERTICES + 1; break;
				default: EXIT(ERR_NOT_IMPLEMENTED); break;
			}
		}
	}

	fprintf(this->out_file, "# vtk DataFile Version 2.0\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "ASCII\n");

	// dataset
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "DATASET UNSTRUCTURED_GRID\n");

	// points
	fprintf(this->out_file, "POINTS %ld %s\n", mesh->vertices.count(), "float");
	for (Word_t i = mesh->vertices.first(); i != INVALID_IDX; i = mesh->vertices.next(i)) {
		Vertex *v = mesh->vertices[i];
		fprintf(this->out_file, "%e %e %e\n", v->x, v->y, v->z);
	}
	fprintf(this->out_file, "\n");

	// cells
	fprintf(this->out_file, "CELLS %ld %ld\n", fc, fc_sz);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		for (int iface = 0; iface < element->get_num_of_faces(); iface++) {
			Word_t vtcs[element->get_face_num_of_vertices(iface)];
			element->get_face_vertices(iface, vtcs);
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::INNER) continue;

			switch (facet->mode) {
				case MODE_TRIANGLE: fprintf(this->out_file, "3 %ld %ld %ld\n", vtcs[0] - 1, vtcs[1] - 1, vtcs[2] - 1); break;
				case MODE_QUAD: fprintf(this->out_file, "4 %ld %ld %ld %ld\n", vtcs[0] - 1, vtcs[1] - 1, vtcs[2] - 1, vtcs[3] - 1); break;
				default: EXIT(ERR_NOT_IMPLEMENTED); break;
			}
		}
	}
	fprintf(this->out_file, "\n");

	// cell types
	fprintf(this->out_file, "CELL_TYPES %ld\n", fc);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		for (int iface = 0; iface < element->get_num_of_faces(); iface++) {
			Word_t vtcs[element->get_face_num_of_vertices(iface)];
			element->get_face_vertices(iface, vtcs);
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::INNER) continue;

			switch (facet->mode) {
				case MODE_TRIANGLE: fprintf(this->out_file, "%d\n", VTK_TRIANGLE); break;
				case MODE_QUAD: fprintf(this->out_file, "%d\n", VTK_QUAD); break;
				default: EXIT(ERR_NOT_IMPLEMENTED); break;
			}
		}
	}
	fprintf(this->out_file, "\n");

	// cell data
	fprintf(this->out_file, "CELL_DATA %ld\n", fc);
	fprintf(this->out_file, "SCALARS %s %s %d\n", name, "float", 1);
	fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");

	// values
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		for (int iface = 0; iface < element->get_num_of_faces(); iface++) {
			Word_t vtcs[element->get_face_num_of_vertices(iface)];
			element->get_face_vertices(iface, vtcs);
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::INNER) continue;

			Boundary *bnd = mesh->boundaries[facet->right];
			fprintf(this->out_file, "%d\n", bnd->marker);
		}
	}
}

void VtkOutputEngine::out_orders(Space *space, const char *name) {
	_F_
	Mesh *mesh = space->get_mesh();

	fprintf(this->out_file, "# vtk DataFile Version 2.0\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "ASCII\n");

	// FIXME: this function is hex-specific

	// HEX specific
	Array<Vertex *> out_vtcs;	// vertices
	Array<int> vtx_pt;			// mapping from mesh vertex id to output vertex id
	MapHSOrd face_pts;			// id of points on faces
	MapHSOrd ctr_pts;			// id of points in the center

	// nodes
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int nv = Hex::NUM_VERTICES;
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			Word_t idx = out_vtcs.add(new Vertex(*v));
			vtx_pt[vtcs[i]] = idx;
		}

		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			Word_t fvtcs[Quad::NUM_VERTICES];
			element->get_face_vertices(iface, fvtcs);

			Word_t k[] = { fvtcs[0], fvtcs[1], fvtcs[2], fvtcs[3] };
			Word_t idx = INVALID_IDX;
			if (!face_pts.lookup(k, Quad::NUM_VERTICES, idx)) {
				// create new vertex
				Vertex *v[4] = { mesh->vertices[fvtcs[0]], mesh->vertices[fvtcs[1]], mesh->vertices[fvtcs[2]], mesh->vertices[fvtcs[3]] };
				Vertex *fcenter = new Vertex((v[0]->x + v[2]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[2]->z) / 2.0);
				Word_t idx = out_vtcs.add(fcenter);
				face_pts.set(k, Quad::NUM_VERTICES, idx);
			}
		}

		Word_t c[] = { vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7] };
		Word_t idx = INVALID_IDX;
		if (!ctr_pts.lookup(c, Hex::NUM_VERTICES, idx)) {
			// create new vertex
			Vertex *v[4] = { mesh->vertices[vtcs[0]], mesh->vertices[vtcs[1]], mesh->vertices[vtcs[3]], mesh->vertices[vtcs[4]] };
			Vertex *center = new Vertex((v[0]->x + v[1]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[3]->z) / 2.0);
			Word_t idx = out_vtcs.add(center);
			ctr_pts.set(c, Hex::NUM_VERTICES, idx);
		}
	}

	// dataset
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(this->out_file, "POINTS %ld %s\n", out_vtcs.count(), "float");
	for (Word_t i = out_vtcs.first(); i != INVALID_IDX; i = out_vtcs.next(i)) {
		Vertex *v = out_vtcs[i];
		fprintf(this->out_file, "%e %e %e\n", v->x, v->y, v->z);
		delete v;																	// we no longer need the vertex data
	}

	// faces that edges coincide with
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "CELLS %ld %ld\n",
	        mesh->get_num_active_elements() * Hex::NUM_EDGES,
	        mesh->get_num_active_elements() * Hex::NUM_EDGES * 5
	);

	int eface[][2] = {
		{ 2, 4 }, { 1, 4 }, { 3, 4 }, { 0, 4 },
		{ 2, 0 }, { 1, 2 }, { 3, 1 }, { 0, 3 },
		{ 2, 5 }, { 1, 5 }, { 3, 5 }, { 0, 5 }
	};
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[element->get_num_of_vertices()];
		element->get_vertices(vtcs);

		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			Word_t fvtcs[2][Quad::NUM_VERTICES];
			element->get_face_vertices(eface[iedge][0], fvtcs[0]);
			element->get_face_vertices(eface[iedge][1], fvtcs[1]);

			Word_t fidx[2] = { INVALID_IDX, INVALID_IDX };
			face_pts.lookup(fvtcs[0], Quad::NUM_VERTICES, fidx[0]);
			face_pts.lookup(fvtcs[1], Quad::NUM_VERTICES, fidx[1]);

			Word_t evtcs[2];
			element->get_edge_vertices(iedge, evtcs);
			Word_t v[4] = { vtx_pt[evtcs[0]], fidx[0], vtx_pt[evtcs[1]], fidx[1] };
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
