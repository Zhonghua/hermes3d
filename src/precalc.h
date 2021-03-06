// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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

#ifndef _PRECALC_H_
#define _PRECALC_H_

#include "function.h"
#include "shapeset.h"

//
//
//
class PrecalcShapeset : public RealFunction {
public:
	/// Constructs a standard precalculated shapeset class.
	/// @param shapeset [in] Pointer to the shapeset to be precalculated.
	PrecalcShapeset(Shapeset *shapeset);

	/// Constructs a slave precalculated shapeset class. The slave instance
	/// does not hold any precalculated tables. Instead it refers to those
	/// contained in the master instance. This is used for test functions
	/// when calling bilinear forms.
	/// @param master_pss [in] Master precalculated shapeset pointer.
	PrecalcShapeset(PrecalcShapeset *master_pss);

	/// Destructor.
	virtual ~PrecalcShapeset();

	virtual void set_quad(Quad3D *quad_3d);

	/// Frees all precalculated tables.
	virtual void free();

	/// Ensures subsequent calls to get_active_element() will be returning 'e'.
	/// Switches the class to the appropriate mode (triangle, quad).
	virtual void set_active_element(Element *e);

	/// Activates a shape function given by its index. The values of the shape function
	/// can then be obtained by setting the required integration rule order by calling
	/// set_quad_order() and after that calling get_values(), get_dx_values(), etc.
	/// @param index [in] Shape index.
	void set_active_shape(int index);

	/// @return Index of the active shape (can be negative if the shape is constrained).
	int get_active_shape() const { return index; };

	/// @return Pointer to the shapeset which is being precalculated.
	Shapeset *get_shapeset() const { return shapeset; }

	///
	void set_master_transform();

	void dump_info(int quad, FILE *file);

protected:
	Shapeset *shapeset;

	void *tables;		/// primary Judy array of shapes

	int index;			/// index of active shape

	PrecalcShapeset *master_pss;

	bool is_slave() const { return master_pss != NULL; }

	virtual void precalculate(qorder_t order, int mask);

	/// Forces a transform without using push_transform() etc.
	/// Used by the Solution class. <b>For internal use only</b>.
	void force_transform(uint64 sub_idx, Trf *ctm) {
		this->sub_idx = sub_idx;
		this->ctm = ctm;
	}

	friend class Solution;
	friend class RefMap;
};


#endif
