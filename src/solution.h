// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2009 Pavel Kus <pavel.kus@gmail.com>
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

#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "function.h"
#include "precalc.h"
#include "space.h"
#include "asmlist.h"
#include "refmap.h"

/// @defgroup solutions Solutions
///
/// TODO: description

typedef
	scalar (*exact_fn_t)(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz);

/// Represents a function defined on a mesh.
///
/// MeshFunction is a base class for all classes representing an arbitrary function
/// superimposed on a mesh (ie., domain). These include the Solution, ExactSolution
/// and Filter classes, which define the concrete behavior and the way the function
/// is (pre)calculated. Any such function can later be visualized.
///
/// (This is an abstract class and cannot be instantiated.)
///
/// @ingroup solutions
class MeshFunction : public ScalarFunction {
public:
	MeshFunction(Mesh *mesh);
	virtual ~MeshFunction();

	virtual void set_quad(Quad3D *quad);
	virtual void set_active_element(Element *e);

	Mesh *get_mesh() const { return mesh; }
	RefMap *get_refmap() { update_refmap(); return refmap; }

protected:
	Mesh *mesh;
	RefMap *refmap;

	int mode, seq;
	bool noinc;

public:
	/// For internal use only.
	void force_transform(MeshFunction *mf)
		{ ScalarFunction::force_transform(mf->get_transform(), mf->get_ctm()); }
	void update_refmap()
		{ refmap->force_transform(sub_idx, ctm); }
};


/// Represents the solution of a PDE.
///
/// Solution represents the solution of a PDE. Given a space and a solution vector,
/// it calculates the appripriate linear combination of basis function at the specified
/// element and integration points.
///
/// @ingroup solutions
class Solution : public MeshFunction {
public:
	Solution(Mesh *mesh);
	virtual ~Solution();
	virtual void free();

	void set_space_and_pss(Space *space, PrecalcShapeset *pss);
	void set_solution_vector(scalar *vec, bool owner);
	void set_zero_vector();

	virtual void set_quad(Quad3D *quad);
	virtual void set_active_element(Element *e);

	void save_solution_vector(char *filename, int ndofs);
	void load_solution_vector(char *filename, int ndofs);

	/// \return The solution value a the given reference domain point (x, y, z)
	/// This is a special-purpose function, use the cached versions get_fn_values,
	/// get_dx_values, etc., instead for normal pusposes.
	/// \param [in] x,y,z - x, y, z coordinate of the point
	/// \param [in] which (0..value, 1..dx, 2..dy, 3..dxx, 4..dyy, 5..dxy) -- 0 by default
	/// \param [in] component - [0..1] solution component, 0 by default
	scalar get_sln_value(double x, double y, double z, EValueType which = FN, int component = 0);

	/// Enables or disables transformation of the solution derivatives (H1 case)
	/// or values (vector (Hcurl) case). This means FN_DX_0 and FN_DY_0 or
	/// FN_VAL_0 and FN_VAL_1 will or will not be returned premultiplied by the reference
	/// mapping matrix. The default is enabled (true).
	void enable_transform(bool enable);

	Space *get_space() const { return space; }
	PrecalcShapeset *get_pss() const { return pss; }
	scalar *get_solution_vector() const { return vec; }

protected:
	static const int NUM_ELEMENTS = 4;

	Space *space;
	PrecalcShapeset *pss;
	PrecalcShapeset *slave_pss;

	scalar *vec;
	bool owner;
	bool transform;

	AsmList al[NUM_ELEMENTS];        	///< assembly lists for last used elements
	void *tables[8][NUM_ELEMENTS];   	///< precalculated tables for last used elements
	Element *elems[8][NUM_ELEMENTS];
	int cur_elem, oldest[8];

	virtual void precalculate(qorder_t qord, int mask);

	void free_tables();
};


/// Represents an exact solution to a PDE.
///
/// ExactSolution represents an arbitrary user-specified function defined on
/// a domain (mesh), typically an exact solution to a PDE. This can be used to
/// compare an approximate solution with an exact solution (see DiffFilter).
///
/// @ingroup solutions
class ExactSolution : public MeshFunction {
public:
	ExactSolution(Mesh *mesh, exact_fn_t fn0, exact_fn_t fn1 = NULL, exact_fn_t fn2 = NULL);
	virtual ~ExactSolution();
	virtual void free();

	virtual void set_active_element(Element *e);

protected:
	exact_fn_t fn[COMPONENTS];
	void *tables[8];

	virtual void precalculate(qorder_t qord, int mask);
};


/// Represents a constant function.
///
/// ConstantSolution represents a constant function on a domain (mesh).
///
/// @ingroup solutions
class ConstantSolution : public ExactSolution {
public:
	ConstantSolution(Mesh *mesh, scalar c0, scalar c1 = 0.0, scalar c2 = 0.0)
		: ExactSolution(mesh, NULL, NULL) { c[0] = c0; c[1] = c1; c[2] = c2; }

protected:
	scalar c[COMPONENTS];

	virtual void precalculate(qorder_t qord, int mask);
};


/// Represents the zero function.
///
/// ZeroSolution represents the zero function on a domain (mesh).
///
/// @ingroup solutions
class ZeroSolution : public ConstantSolution {
public:
	ZeroSolution(Mesh *mesh) : ConstantSolution(mesh, 0.0) { }
};


#endif
