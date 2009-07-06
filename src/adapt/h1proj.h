// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 David Andrs <dandrs@unr.edu>
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

#ifndef _ADAPT_H1_PROJECTION_H_
#define _ADAPT_H1_PROJECTION_H_

#include "proj.h"

/// H1 projection
///
/// FIXME: hex specific
///
/// @ingroup hp-adaptivity
class H1Projection : public Projection {
public:
	H1Projection(Solution *afn, Element *e, Shapeset *ss) : Projection(afn, e, ss) {}

	virtual double get_error(int split, int son, order3_t order);

protected:
	virtual void calc_vertex_proj(int split, int son);
	virtual void calc_edge_proj(int edge, int split, int son, order3_t order);
	virtual void calc_face_proj(int face, int split, int son, order3_t order);
	virtual void calc_bubble_proj(int split, int son, order3_t order);

	static double mdx[8];
	static double mdy[8];
	static double mdz[8];
};


#endif
