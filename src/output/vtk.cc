//
// vtkoutputengine.cc
//
//

#include "../config.h"
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

namespace Vtk {

//// OutputQuad //////////////////////////////////////////////////////////////////////////////

/// Common ancestor for output quadratures. Extends the interface of Quad3D
///
/// @ingroup visualization
class OutputQuad : public Quad3D {
public:
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
//	virtual void calculate_view_points(int order);
};

OutputQuadTetra::OutputQuadTetra() {
}

OutputQuadTetra::~OutputQuadTetra() {

}



/// Quadrature for visualizing the solution on hexahedron
///
/// @ingroup visualization
class OutputQuadHex : public OutputQuad {
public:
	OutputQuadHex();
	virtual ~OutputQuadHex();

	virtual QuadPt3D *get_points(int order) {
		if (tables[order] == NULL) calc_table(order);
		return tables[order];
	}

protected:
	void calc_table(int order);
};


OutputQuadHex::OutputQuadHex() {
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;
	max_order = MAKE_HEX_ORDER(MAX_QUAD_ORDER, MAX_QUAD_ORDER, MAX_QUAD_ORDER);
	max_edge_order = 0;				// not used
	max_face_order = 0;				// not used

	int num = max_order + 1;

	// allocate memory for tables (tables are calculated on demand)
	tables = new QuadPt3D *[num];
	memset(tables, 0, num * sizeof(QuadPt3D *));

	// number of integration points
	np = new int [num];
	memset(np, 0, num * sizeof(int));

	for (int i = 0; i <= MAX_QUAD_ORDER; i++) {
		for (int j = 0; j <= MAX_QUAD_ORDER; j++) {
			for (int o = 0; o <= MAX_QUAD_ORDER; o++) {
				int m = MAKE_HEX_ORDER(i, j, o);
				np[m] = (i + 1) * (j + 1) * (o + 1);
			}
		}
	}
#endif
}

OutputQuadHex::~OutputQuadHex() {
#ifdef WITH_HEX
	// free hexa points
	for (int i = 0; i <= MAX_QUAD_ORDER; i++) {
		for (int j = 0; j <= MAX_QUAD_ORDER; j++) {
			for (int o = 0; o <= MAX_QUAD_ORDER; o++) {
				int m = MAKE_HEX_ORDER(i, j, o);
				delete [] tables[m];
			}
		}
	}
	delete [] tables;
	delete [] np;
#endif
}

void OutputQuadHex::calc_table(int order) {
#ifdef WITH_HEX
	tables[order] = new QuadPt3D[np[order]];
	int i, j, o;
	i = GET_HEX_ORDER_1(order);
	j = GET_HEX_ORDER_2(order);
	o = GET_HEX_ORDER_3(order);
	double step_i, step_j, step_o;
	step_i = 2.0 / i;
	step_j = 2.0 / j;
	step_o = 2.0 / o;

	for (int k = 0, n = 0; k < i + 1; k++) {
		for (int l = 0; l < j + 1; l++) {
			for (int p = 0; p < o + 1; p++, n++) {
				assert(n < np[order]);
				tables[order][n].x = (step_i * k) - 1;
				tables[order][n].y = (step_j * l) - 1;
				tables[order][n].z = (step_o * p) - 1;
				tables[order][n].w = 0.0;	// not used
			}
		}
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

void VtkOutputEngine::dump_points(MeshFunction *fn) {
	Array<Point3D *> vertices;
	Array<int *> cells[3];				// 3 types of elements

	int type[] = { MODE_TETRAHEDRON, MODE_HEXAHEDRON, MODE_PRISM };
	int type_id[] = { VTK_TETRA, VTK_HEXAHEDRON, VTK_WEDGE };
	Mesh *mesh = fn->get_mesh();

	Word_t idx;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();

		// FIXME: get order from the space and set sufficient division
		int order;
		switch (mode) {
			case MODE_HEXAHEDRON: order = MAKE_HEX_ORDER(8, 8, 8); break;
		}

		Vtk::OutputQuad *quad = output_quad[mode];
		fn->set_active_element(element);
		fn->set_quad(quad);

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
			case MODE_HEXAHEDRON: {
				int size_i, size_j, size_o;
				size_i = GET_HEX_ORDER_1(order);
				size_j = GET_HEX_ORDER_2(order);
				size_o = GET_HEX_ORDER_3(order);

				for (int i = 0; i < size_i; i++) {
					for (int j = 0; j < size_j; j++) {
						for (int o = 0; o < size_o; o++) {
							int *cell = new int [Hex::NUM_VERTICES];
							cell[0] = (size_o + 1) * (i * (size_j + 1) + j) + o;
							cell[1] = cell[0] + ((size_j + 1) * (size_o + 1));
							cell[2] = cell[1] + (size_o + 1);
							cell[3] = cell[0] + (size_o + 1);
							cell[4] = cell[0] + 1;
							cell[5] = cell[1] + 1;
							cell[6] = cell[2] + 1;
							cell[7] = cell[3] + 1;

							cells[1].add(cell);
						}
					}
				}
			} break;

			case MODE_PRISM:
			case MODE_TETRAHEDRON: EXIT(ERR_NOT_IMPLEMENTED); break;

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
	fprintf(this->out_file, "POINTS %d %s\n", vertices.count(), "float");
	for (Word_t i = vertices.first(); i != INVALID_IDX; i = vertices.next(i)) {
		Point3D *pt = vertices[i];
		fprintf(this->out_file, "%e %e %e\n", pt->x, pt->y, pt->z);
	}
	fprintf(this->out_file, "\n");

	fprintf(this->out_file, "CELLS %d %d\n", cells[0].count() + cells[1].count() + cells[2].count(),
		(Tetra::NUM_VERTICES + 1) * cells[0].count() +
		(Hex::NUM_VERTICES + 1) * cells[1].count() +
		(Prism::NUM_VERTICES + 1) * cells[2].count());
	for (int i = 0; i < countof(type); i++) {			// 3 types of elements
		int pt_cnt[] = { Tetra::NUM_VERTICES, Hex::NUM_VERTICES, Prism::NUM_VERTICES };
		for (int j = cells[type[i]].first(); j != INVALID_IDX; j = cells[type[i]].next(j)) {
			fprintf(this->out_file, "%d", pt_cnt[i]);
			for (int k = 0; k < pt_cnt[i]; k++)
				fprintf(this->out_file, " %d", cells[type[i]][j][k]);
			fprintf(this->out_file, "\n");
		}
	}
	fprintf(this->out_file, "\n");

	fprintf(this->out_file, "CELL_TYPES %d\n", cells[0].count() + cells[1].count() + cells[2].count());
	for (int i = 0; i < countof(type); i++) {			// 3 types of elements
		for (int j = 0; j < cells[type[i]].count(); j++)
			fprintf(this->out_file, "%d\n", type_id[i]);
	}

	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "POINT_DATA %d\n", vertices.count());

	// free allocated memory
	for (Word_t i = vertices.first(); i != INVALID_IDX; i = vertices.next(i))
		delete vertices[i];

	for (int i = 0; i < 3; i++)			// 3 types of elements
		for (Word_t j = cells[i].first(); j != INVALID_IDX; j = cells[i].next(j))
			delete cells[i][j];
}

void VtkOutputEngine::dump_scalars(const char *name, MeshFunction *fn, int item) {
	Mesh *mesh = fn->get_mesh();

	// header
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "SCALARS %s %s %d\n", name, "float", 1);
	fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");

	// values
	Word_t idx;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();

		// FIXME: get order from the space and set sufficient division
		int order = MAKE_HEX_ORDER(8, 8, 8);

		switch (mode) {
			// FIXME
			case MODE_TETRAHEDRON:
			case MODE_HEXAHEDRON:
				break;
/*			case MODE_TETRAHEDRON: order = space->get_element_order(idx); break;
			case MODE_HEXAHEDRON: order = get_principal_order(space->get_element_order(idx)); break;
*/			case MODE_PRISM: EXIT(ERR_NOT_IMPLEMENTED); break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}

		Vtk::OutputQuad *quad = output_quad[mode];
		fn->set_active_element(element);
		fn->set_quad(quad);
		fn->set_quad_order(ELEM_QORDER(order), item);
		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);
		scalar *val = fn->get_values(a, b);

		int np = quad->get_num_points(order);
		for (int i = 0; i < np; i++){
			assert(fabs(IMAG(val[i])) < 1e-12);
			fprintf(this->out_file, "%e\n", REAL(val[i]));
		}
	}
}

void VtkOutputEngine::dump_vectors(const char *name, MeshFunction *fn, int item) {
	Mesh *mesh = fn->get_mesh();

	// header
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "VECTORS %s %s\n", name, "float");

	// values
	Word_t idx;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();

		// FIXME: get order from the space and set sufficient division
		int order = MAKE_HEX_ORDER(8, 8, 8);

		switch (mode) {
			// FIXME
			case MODE_TETRAHEDRON:
			case MODE_HEXAHEDRON:
				break;
/*			case MODE_TETRAHEDRON: order = space->get_element_order(idx); break;
			case MODE_HEXAHEDRON: order = get_principal_order(space->get_element_order(idx)); break;
*/			case MODE_PRISM: EXIT(ERR_NOT_IMPLEMENTED); break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}

		Vtk::OutputQuad *quad = output_quad[mode];
		fn->set_active_element(element);
		fn->set_quad(quad);
		fn->set_quad_order(ELEM_QORDER(order), item);
		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);
		scalar *val[3];
		for (int comp = 0; comp < 3; comp++)	// FIXME: magic number
			val[comp] = fn->get_values(comp, b);

		int np = quad->get_num_points(order);
		for (int i = 0; i < np; i++){
			assert(fabs(IMAG(val[0][i])) < 1e-12);
			assert(fabs(IMAG(val[1][i])) < 1e-12);
			assert(fabs(IMAG(val[2][i])) < 1e-12);
			fprintf(this->out_file, "%e %e %e\n", REAL(val[0][i]), REAL(val[1][i]), REAL(val[2][i]));
		}
	}
}

void VtkOutputEngine::out(MeshFunction *fn, const char *name, int item/* = FN_VAL_0*/) {
	if (!has_points) {
		// calculate points where we will evaluate the function 'fn'
		dump_points(fn);
		this->has_points = true;
	}

	assert(fn->get_num_components() == 1 || fn->get_num_components() == 3);
	// dump data
	// FIXME: better check for 'item' argument
	// if (item == FN_VAL_0 | FN_DX_1 | FN_DY_2) then the following condition will be TRUE
	// but who knows what will happen. For now, let's suppose the user is sane.
	if ((fn->get_num_components() == 3) && ((item & FN_COMPONENT_0) && (item & FN_COMPONENT_1) && (item & FN_COMPONENT_2)))
		dump_vectors(name, fn, item);
	else
		dump_scalars(name, fn, item);
}

void VtkOutputEngine::out(Mesh *mesh) {
	// Not implemented
}
