// This file is part of Hermes3D
//
// Copyright (c) 20010 hp-FEM group at the University of Nevada, Reno (UNR).
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

#ifndef CURVED_H_
#define CURVED_H_

#include "common.h"
#include "mesh.h"
#include "order.h"

/// Represents a curved face.
///
/// The structure contains data for representing a curved face (bicubic patch)
struct Curve {
	static const int PTS = 16;			// 16 control points

	Point3D kij[PTS];					// control point: K_i,j
};

/// CurvMap is a structure storing complete information on the curved face of
/// an element. There are two variants of this structure. The first if for
/// top-level (master mesh) elements.
///
struct CurvMap {
	CurvMap() { coefs = NULL; nc = 0; }
	CurvMap(CurvMap *cm);
	~CurvMap();

	// this structure defines a curved mapping of an element; it has two
	// modes, depending on the value of 'toplevel'
	bool toplevel;
	union {
		// if toplevel=true, this structure belongs to a base mesh element
		// and the array 'nurbs' points to (up to six) NURBS curved faces
		Curve *curves[6];
		struct {
			// if toplevel=false, this structure belongs to a refined element
			// and 'parent' points to the base mesh element CurvMap structure;
			Element *parent;
			uint64_t part;
		};
	};

	// current polynomial degree of the refmap approximation
	order3_t order;

	// finally here are the coefficients of the higher-order basis functions
	// that constitute the projected reference mapping:
	int nc;			// number of coefficients
	Vertex *coefs;	// array of the coefficients

	// this is called for every curvilinear element when it is created
	// or when it is necessary to re-calculate coefficients for another
	// order: 'e' is a pointer to the element to which this CurvMap
	// belongs to. First, old "coefs" are removed if they are not NULL,
	// then new coefficients are projected.
	void update_refmap_coefs(Mesh *mesh, Element *e);
};

#endif /* CURVED_H_ */
