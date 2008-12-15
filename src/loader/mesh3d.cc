//
// mesh3dloader.cc
//
// Loading mesh from mesh3d format
//

#include "../config.h"

#include <string.h>
#include <common/trace.h>
#include <common/error.h>
#include "mesh3d.h"
#include "../mesh.h"
#include "../refdomain.h"

// maximal row length in bytes (used for reading the mesh3d-file)
#define MAX_ROW_LEN							1024

// exception error codes
#define E_CANT_OPEN_FILE					-1
#define E_READ_ERROR						-2

// number of markers on mesd3d file
#define MARKERS								1

Mesh3DReader::Mesh3DReader() {
}

Mesh3DReader::~Mesh3DReader() {
}

static bool read_num(char *row, int &vertex_count) {
	int n;
	if (sscanf(row, "%d", &n) == 1) {
		vertex_count = n;
		return true;
	} else
		return false;
}

static bool read_n_nums(char *row, int n, double values[]) {
	int i = 0;
	char delims[] = " \t\n\r";
	char *token = strtok(row, delims);
	while (token != NULL && i < n) {
		double d;
		sscanf(token, "%lf", &d);
		values[i++] = d;

		token = strtok(NULL, delims);
	}

	return (i == n);
}

static bool read_n_nums(char *row, int n, Word_t values[]) {
	int i = 0;
	char delims[] = " \t\n\r";
	char *token = strtok(row, delims);
	while (token != NULL && i < n) {
		int n;
		sscanf(token, "%d", &n);
		values[i++] = n;

		token = strtok(NULL, delims);
	}

	return (i == n);
}

static bool range_check(int max_index, Word_t *vs, int num_vs) {
	for (int i = 0; i < num_vs; i++) {
		if (vs[i] <= 0 || vs[i] > max_index)
			return false;
	}
	return true;
}

bool Mesh3DReader::load(const char *file_name, Mesh *mesh) {
	assert(mesh != NULL);

	FILE *file = fopen(file_name, "r");
	if (file == NULL)
		return false;

	try {
		enum EState {
			STATE_VERTICES_NUM,
			STATE_VERTICES,
			STATE_TETRAS_NUM,
			STATE_TETRAS,
			STATE_HEXES_NUM,
			STATE_HEXES,
			STATE_PRISMS_NUM,
			STATE_PRISMS,
			STATE_TRIS_NUM,
			STATE_TRIS,
			STATE_QUADS_NUM,
			STATE_QUADS,
			STATE_OK
		} state = STATE_VERTICES_NUM;

		int vertex_count = 0;
		int tetra_count = 0;
		int hex_count = 0;
		int prism_count = 0;
		int tri_count = 0;
		int quad_count = 0;
		double buffer[10];
		Word_t vs[10];

		int max_vertex_index = 0;

		char row[MAX_ROW_LEN] = {0};
		while (fgets(row, MAX_ROW_LEN, file) != NULL && state != STATE_OK) {
			if (row[0] == '#') continue; // comment

			// trim right
			int i = strlen(row) - 1;
			for (; i >= 0; i--) {
				if (row[i] != ' ' && row[i] != '\r' && row[i] != '\n' && row[i] != '\t')
					break;
			}
			row[i + 1] = '\0';
			if (strlen(row) <= 0) continue; // skip empty lines

			switch (state) {
				case STATE_VERTICES_NUM:
					if (read_num(row, vertex_count)) {
						state = STATE_VERTICES;
						if (vertex_count <= 0)
							throw E_READ_ERROR;
						max_vertex_index = vertex_count; //vertices are counted from 1 in mesh3d format
					}
					break;

				case STATE_VERTICES:
					if (read_n_nums(row, Vertex::NUM_COORDS, buffer)) {
						mesh->add_vertex(buffer[0], buffer[1], buffer[2]);

						vertex_count--;
						if (vertex_count == 0)
							state = STATE_TETRAS_NUM;
					}
					break;

				case STATE_TETRAS_NUM:
					if (read_num(row, tetra_count)) {
						state = STATE_TETRAS;
						if (tetra_count <= 0)
							state = STATE_HEXES_NUM;
					}
					break;

				case STATE_TETRAS:
				if (read_n_nums(row, Tetra::NUM_VERTICES, vs)) {
					if (!range_check(max_vertex_index, vs, Tetra::NUM_VERTICES)) {
						ERROR("Invalid vertex index was found in the section defining tetras.");
						throw E_READ_ERROR;
					}

					for (int i = 0; i < Tetra::NUM_VERTICES; i++) vs[i]--;
					Tetra *tetra = mesh->add_tetra(vs);

					tetra_count--;
					if (tetra_count == 0)
						state = STATE_HEXES_NUM;
				}
				break;

				case STATE_HEXES_NUM:
					if (read_num(row, hex_count)) {
						state = STATE_HEXES;
						if (hex_count <= 0)
							state = STATE_PRISMS_NUM;
					}
					break;

				case STATE_HEXES:
					if (read_n_nums(row, Hex::NUM_VERTICES, vs)) {
						if (!range_check(max_vertex_index, vs, Hex::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining hexes.");
							throw E_READ_ERROR;
						}

						for (int i = 0; i < Hex::NUM_VERTICES; i++) vs[i]--;
						Hex *hex = mesh->add_hex(vs);

						hex_count--;
						if (hex_count == 0)
							state = STATE_PRISMS_NUM;
					}
					break;

				case STATE_PRISMS_NUM:
					if (read_num(row, prism_count)) {
						state = STATE_PRISMS;
						if (prism_count <= 0)
							state = STATE_TRIS_NUM;
					}
					break;

				case STATE_PRISMS:
					if (read_n_nums(row, Prism::NUM_VERTICES, vs)) {
						if (!range_check(max_vertex_index, vs, Prism::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining prisms.");
							throw E_READ_ERROR;
						}
						// add prism
						for (int i = 0; i < Prism::NUM_VERTICES; i++) vs[i]--;
						Prism *prism = mesh->add_prism(vs);
						prism_count--;
						if (prism_count == 0)
							state = STATE_TRIS_NUM;
					}
					break;

				case STATE_TRIS_NUM:
					if (read_num(row, tri_count)) {
						state = STATE_TRIS;
						if (tri_count <= 0)
							state = STATE_QUADS_NUM;
					}
					break;

				case STATE_TRIS:
					if (read_n_nums(row, Tri::NUM_VERTICES + MARKERS, vs)) {
						if (!range_check(max_vertex_index, vs, Tri::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining tris.");
							throw E_READ_ERROR;
						}

						Word_t facet_idxs[Tri::NUM_VERTICES] = {vs[0] - 1, vs[1] - 1, vs[2] - 1};
						mesh->add_tri_boundary(facet_idxs, vs[Tri::NUM_VERTICES]);
						tri_count--;
						if (tri_count == 0)
							state = STATE_QUADS_NUM;
					}
					else {
						ERROR("Not enough information for tris. You probably forgot to define boundary condition.");
						throw E_READ_ERROR;
					}
					break;

				case STATE_QUADS_NUM:
					if (read_num(row, quad_count)) {
						state = STATE_QUADS;
						if (quad_count <= 0)
							state = STATE_QUADS;
					}
					break;

				case STATE_QUADS:
					if (read_n_nums(row, Quad::NUM_VERTICES + MARKERS, vs)) {
						if (!range_check(max_vertex_index, vs, Quad::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining quads.");
							throw E_READ_ERROR;
						}

						Word_t facet_idxs[Quad::NUM_VERTICES] = { vs[0] - 1, vs[1] - 1, vs[2] - 1, vs[3] - 1 };
						mesh->add_quad_boundary(facet_idxs, vs[Quad::NUM_VERTICES]);
						quad_count--;
						if (quad_count == 0)
							state = STATE_OK;
					}
					else {
						ERROR("Not enough information for quads. You probably forgot to define boundary condition.");
						throw E_READ_ERROR;
					}
					break;

				case STATE_OK:
				break;
			}
		}

		// check if all "outer" faces have defined boundary condition
		for (Word_t i = mesh->facets.first(); i != INVALID_IDX; i = mesh->facets.next(i)) {
			Facet *facet = mesh->facets.get(i);

			if ((facet->left == INVALID_IDX) || (facet->right == INVALID_IDX)) {
				ERROR("Not all outer faces have defined boundary condition.");
				throw E_READ_ERROR;
			}
		}

		mesh->ugh();
	}
	catch (int e) {
		fclose(file);
		return false;
	}

	fclose(file);

	return true;
}

bool Mesh3DReader::save(const char *file_name, Mesh *mesh) {
	assert(mesh != NULL);

	FILE *file = fopen(file_name, "w");
	if (file == NULL)
		return false;

	int count;

	// save vertices
	fprintf(file, "# vertices\n");
	fprintf(file, "%d\n", mesh->vertices.count());
	for (Word_t i = mesh->vertices.first(); i != INVALID_IDX; i = mesh->vertices.next(i)) {
		Vertex *v = mesh->vertices[i];
		fprintf(file, "%lf %lf %lf\n", v->x, v->y, v->z);
	}
	fprintf(file, "\n");

	// elements
	Array<Element *> tet, hex, pri;
    for (Word_t i = mesh->elements.first(); i != INVALID_IDX; i = mesh->elements.next(i)) {
    	Element *elem = mesh->elements[i];
    	if (elem->active) {
	    	switch (elem->get_mode()) {
	    		case MODE_TETRAHEDRON: tet.add(elem); break;
	    		case MODE_HEXAHEDRON: hex.add(elem); break;
	    		case MODE_PRISM: pri.add(elem); break;
	    	}
	    }
    }

	// save tetras
	fprintf(file, "# tetras\n");
	fprintf(file, "%d\n", tet.count());
    for (Word_t i = tet.first(); i != INVALID_IDX; i = tet.next(i)) {
    	Word_t vtcs[Tetra::NUM_VERTICES];
    	tet[i]->get_vertices(vtcs);
		fprintf(file, "%d %d %d %d\n",
			vtcs[0] + 1, vtcs[1] + 1, vtcs[2] + 1, vtcs[3] + 1);
    }
	fprintf(file, "\n");

	// save hexes
	fprintf(file, "# hexes\n");
	fprintf(file, "%d\n", hex.count());
    for (Word_t i = hex.first(); i != INVALID_IDX; i = hex.next(i)) {
    	Word_t vtcs[Hex::NUM_VERTICES];
    	hex[i]->get_vertices(vtcs);
		fprintf(file, "%d %d %d %d %d %d %d %d\n",
			vtcs[0] + 1, vtcs[1] + 1, vtcs[2] + 1, vtcs[3] + 1, vtcs[4] + 1, vtcs[5] + 1, vtcs[6] + 1, vtcs[7] + 1);
    }
	fprintf(file, "\n");

	// save prisms
	fprintf(file, "# prisms\n");
	fprintf(file, "%d\n", pri.count());
    for (Word_t i = pri.first(); i != INVALID_IDX; i = pri.next(i)) {
    	Word_t vtcs[Prism::NUM_VERTICES];
    	pri[i]->get_vertices(vtcs);
		fprintf(file, "%d %d %d %d %d %d\n",
			vtcs[0] + 1, vtcs[1] + 1, vtcs[2] + 1, vtcs[3] + 1, vtcs[4] + 1, vtcs[5] + 1);
    }
	fprintf(file, "\n");

	// boundaries
	Array<Facet *> tri_facets, quad_facets;
    for (Word_t i = mesh->facets.first(); i != INVALID_IDX; i = mesh->facets.next(i)) {
    	Facet *facet = mesh->facets.get(i);
    	if (facet->type == Facet::OUTER && mesh->elements[facet->left]->active) {
	    	switch (facet->type) {
    			case MODE_TRIANGLE: tri_facets.add(facet); break;
    			case MODE_QUAD: quad_facets.add(facet); break;
    		}
    	}
    }

	// tris
	fprintf(file, "# tris\n");
	fprintf(file, "%d\n", tri_facets.count());
	for (Word_t i = tri_facets.first(); i != INVALID_IDX; i = tri_facets.next(i)) {
		Facet *facet = tri_facets[i];
    	Boundary *bnd = mesh->boundaries[facet->right];
    	Element *elem = mesh->elements[facet->left];

    	Word_t vtcs[Tri::NUM_VERTICES];
    	elem->get_face_vertices(facet->left_face_num, vtcs);

		fprintf(file, "%d %d %d     %d\n",
			vtcs[0] + 1, vtcs[1] + 1, vtcs[2] + 1, bnd->marker);
	}
	fprintf(file, "\n");

	// quads
	fprintf(file, "# quads\n");
	fprintf(file, "%d\n", quad_facets.count());
	for (Word_t i = quad_facets.first(); i != INVALID_IDX; i = quad_facets.next(i)) {
		Facet *facet = quad_facets[i];
    	Boundary *bnd = mesh->boundaries[facet->right];
    	Element *elem = mesh->elements[facet->left];

    	Word_t vtcs[Quad::NUM_VERTICES];
    	elem->get_face_vertices(facet->left_face_num, vtcs);

		fprintf(file, "%d %d %d %d     %d\n",
			vtcs[0] + 1, vtcs[1] + 1, vtcs[2] + 1, vtcs[3] + 1, bnd->marker);
	}
	fprintf(file, "\n");

	// the end
	fclose(file);

	return true;
}
