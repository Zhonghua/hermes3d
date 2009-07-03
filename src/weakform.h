// This file is part of Hermes3D
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2009 David Andrs <andrsd@gmail.com>
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
// This file was taken from hermes2d and adjusted for hermes3d
//

#ifndef _WEAKFORM_H_
#define _WEAKFORM_H_

#include "function.h"
#include "forms.h"


// Bilinear form symmetry flag, see WeakForm::add_biform
enum SymFlag {
	ANTISYM = -1,
	UNSYM = 0,
	SYM = 1
};


/// Represents the weak formulation of a problem.
///
/// The WeakForm class represents the weak formulation of a system of linear PDEs.
/// The number of equations ("neq") in the system is fixed and is passed to the constructor.
/// The weak formulation of the system A(U,V) = L(V) has a block structure. A(U,V) is
/// a (neq x neq) matrix of bilinear forms a_mn(u,v) and L(V) is a neq-component vector
/// of linear forms l(v). U and V are the vectors of basis and test functions.
///
///
/// @ingroup assembling
class WeakForm {
	typedef scalar (*biform_val_t) (int n, double *wt, fn_t<double> *u, fn_t<double> *v, geom_t<double> *e, user_data_t<scalar> *);
	typedef ord_t (*biform_ord_t) (int n, double *wt, fn_t<ord_t> *u, fn_t<ord_t> *v, geom_t<ord_t> *e, user_data_t<ord_t> *);

	typedef scalar (*liform_val_t)(int n, double *wt, fn_t<double> *v, geom_t<double> *e, user_data_t<scalar> *);
	typedef ord_t (*liform_ord_t)(int n, double *wt, fn_t<ord_t> *v, geom_t<ord_t> *e, user_data_t<ord_t> *);

public:
	WeakForm(int neq);
	virtual ~WeakForm();

	int def_area(int n, ...);

	void add_biform(int i, int j, biform_val_t fn, biform_ord_t ord, SymFlag sym = UNSYM, int area = ANY, int nx = 0, ...);
	void add_biform_surf(int i, int j, biform_val_t fn, biform_ord_t ord, int area = ANY, int nx = 0, ...);

	void add_liform(int i, liform_val_t fn, liform_ord_t ord, int area = ANY, int nx = 0, ...);
	void add_liform_surf(int i, liform_val_t fn, liform_ord_t ord, int area = ANY, int nx = 0, ...);

	void set_ext_fns(void *fn, int nx, ...);

	order3_t get_int_order();

protected:
	int neq;

	struct Area  {
		std::vector<int> markers;
	};

	std::vector<Area> areas;

	struct BiFormVol {
		int i, j, sym, area;
		biform_val_t fn;						// callback for evaluating the form
		biform_ord_t ord;						// callback to determine the integration order
		std::vector<MeshFunction *> ext;		// external functions
	};
	struct BiFormSurf {
		int i, j, area;
		biform_val_t fn;
		biform_ord_t ord;
		std::vector<MeshFunction *> ext;
	};
	struct LiFormVol {
		int i, area;
		liform_val_t fn;
		liform_ord_t ord;
		std::vector<MeshFunction *> ext;
	};
	struct LiFormSurf {
		int i, area;
		liform_val_t fn;
		liform_ord_t ord;
		std::vector<MeshFunction *> ext;
	};

	std::vector<BiFormVol>  bfvol;
	std::vector<BiFormSurf> bfsurf;
	std::vector<LiFormVol>  lfvol;
	std::vector<LiFormSurf> lfsurf;


	struct Stage {
		std::vector<int> idx;
		std::vector<Mesh *> meshes;
		std::vector<Transformable *> fns;
		std::vector<MeshFunction *> ext;

		std::vector<BiFormVol *>  bfvol;
		std::vector<BiFormSurf *> bfsurf;
		std::vector<LiFormVol *>  lfvol;
		std::vector<LiFormSurf *> lfsurf;

		std::set<int> idx_set;
		std::set<unsigned> seq_set;
		std::set<MeshFunction *> ext_set;
	};

	void get_stages(Space **spaces, std::vector<Stage> &stages, bool rhsonly);
	bool **get_blocks();

	bool is_in_area(int marker, int area) const
		{ return area >= 0 ? area == marker : is_in_area_2(marker, area); }

	bool is_sym() const { return false; /* not impl. yet */ }

private:
	Stage *find_stage(std::vector<Stage> &stages, int ii, int jj,
	                  Mesh *m1, Mesh *m2, std::vector<MeshFunction *> &ext);

	bool is_in_area_2(int marker, int area) const;


	friend class LinProblem;
};

#endif /* WEAKFORM_H_ */
