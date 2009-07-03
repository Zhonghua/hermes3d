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

#include "common.h"
#include "weakform.h"
#include "matrix.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

WeakForm::WeakForm(int neq) {
	_F_
	this->neq = neq;
}

WeakForm::~WeakForm() {
	_F_
}

#define init_ext_fns \
	va_list ap; va_start(ap, nx); \
	for (int i = 0; i < nx; i++) \
		form.ext.push_back(va_arg(ap, MeshFunction*)); \
	va_end(ap)


void WeakForm::add_biform(int i, int j, biform_val_t fn, biform_ord_t ord, SymFlag sym, int area, int nx, ...) {
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) ERROR("Invalid equation number.");
	if (sym < -1 || sym > 1) ERROR("\"sym\" must be ANTISYM, UNSYM or SYM.");
	if (sym < 0 && i == j) ERROR("Only off-diagonal forms can be antisymmetric.");
	if (area != ANY && area < 0 && -area > areas.size()) ERROR("Invalid area number.");
	if (bfvol.size() > 100) WARNING("Large number of forms (> 100). Is this the intent?");

	BiFormVol form = { i, j, sym, area, fn, ord };
	init_ext_fns;
	bfvol.push_back(form);
}

void WeakForm::add_biform_surf(int i, int j, biform_val_t fn, biform_ord_t ord, int area, int nx, ...) {
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) ERROR("Invalid equation number.");
	if (area != ANY && area < 0 && -area > areas.size()) ERROR("Invalid area number.");

	BiFormSurf form = { i, j, area, fn, ord };
	init_ext_fns;
	bfsurf.push_back(form);
}

void WeakForm::add_liform(int i, liform_val_t fn, liform_ord_t ord, int area, int nx, ...) {
	_F_
	if (i < 0 || i >= neq) ERROR("Invalid equation number.");
	if (area != ANY && area < 0 && -area > areas.size()) ERROR("Invalid area number.");

	LiFormVol form = { i, area, fn, ord };
	init_ext_fns;
	lfvol.push_back(form);
}

void WeakForm::add_liform_surf(int i, liform_val_t fn, liform_ord_t ord, int area, int nx, ...) {
	_F_
	if (i < 0 || i >= neq) ERROR("Invalid equation number.");
	if (area != ANY && area < 0 && -area > areas.size()) ERROR("Invalid area number.");

	LiFormSurf form = { i, area, fn, ord };
	init_ext_fns;
	lfsurf.push_back(form);
}


void WeakForm::set_ext_fns(void *fn, int nx, ...) {
	EXIT(ERR_NOT_IMPLEMENTED);
}


//// stages ////////////////////////////////////////////////////////////////////////////////////////

/// Constructs a list of assembling stages. Each stage contains a list of forms
/// that share the same meshes. Each stage is then assembled separately. This
/// improves the performance of multi-mesh assembling.
///
void WeakForm::get_stages(Space **spaces, std::vector<WeakForm::Stage> &stages, bool rhsonly) {
	_F_
	int i, j;
	stages.clear();

	if (!rhsonly) {
		// process volume biforms
		for (i = 0; i < bfvol.size(); i++) {
			int ii = bfvol[i].i, jj = bfvol[i].j;
			Mesh *m1 = spaces[ii]->get_mesh();
			Mesh *m2 = spaces[jj]->get_mesh();
			Stage *s = find_stage(stages, ii, jj, m1, m2, bfvol[i].ext);
			s->bfvol.push_back(&bfvol[i]);
		}

		// process surface biforms
		for (i = 0; i < bfsurf.size(); i++) {
			int ii = bfsurf[i].i, jj = bfsurf[i].j;
			Mesh *m1 = spaces[ii]->get_mesh();
			Mesh *m2 = spaces[jj]->get_mesh();
			Stage *s = find_stage(stages, ii, jj, m1, m2, bfsurf[i].ext);
			s->bfsurf.push_back(&bfsurf[i]);
		}
	}

	// process volume liforms
	for (i = 0; i < lfvol.size(); i++) {
		int ii = lfvol[i].i;
		Mesh *m = spaces[ii]->get_mesh();
		Stage *s = find_stage(stages, ii, ii, m, m, lfvol[i].ext);
		s->lfvol.push_back(&lfvol[i]);
	}

	// process surface liforms
	for (i = 0; i < lfsurf.size(); i++) {
		int ii = lfsurf[i].i;
		Mesh *m = spaces[ii]->get_mesh();
		Stage *s = find_stage(stages, ii, ii, m, m, lfsurf[i].ext);
		s->lfsurf.push_back(&lfsurf[i]);
	}

	// helper macro for iterating in a set
	#define SET_FOR_EACH(myset, type) \
		for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

	// initialize the arrays meshes and fns needed by Traverse for each stage
	for (i = 0; i < stages.size(); i++) {
		Stage *s = &stages[i];
		SET_FOR_EACH(s->idx_set, int) {
			s->idx.push_back(*it);
			s->meshes.push_back(spaces[*it]->get_mesh());
			s->fns.push_back(NULL);
		}
		SET_FOR_EACH(s->ext_set, MeshFunction *) {
			s->ext.push_back(*it);
			s->meshes.push_back((*it)->get_mesh());
			s->fns.push_back(*it);
		}
		s->idx_set.clear();
		s->seq_set.clear();
		s->ext_set.clear();
	}
}


/// Finds an assembling stage with the same set of meshes as [m1, m2, ext]. If no such
/// stage can be found, a new one is created and returned.
///
WeakForm::Stage *WeakForm::find_stage(std::vector<WeakForm::Stage> &stages, int ii, int jj,
                                      Mesh *m1, Mesh *m2, std::vector<MeshFunction *> &ext)
{
	_F_
	// first create a list of meshes the form uses
	std::set<unsigned> seq;
	seq.insert(m1->get_seq());
	seq.insert(m2->get_seq());
	for (int i = 0; i < ext.size(); i++)
		seq.insert(ext[i]->get_mesh()->get_seq());

	// find a suitable existing stage for the form
	Stage *s = NULL;
	for (int i = 0; i < stages.size(); i++)
		if (seq.size() == stages[i].seq_set.size() && equal(seq.begin(), seq.end(), stages[i].seq_set.begin())) {
			s = &stages[i];
			break;
		}

	// create a new stage if not found
	if (s == NULL) {
		Stage newstage;
		stages.push_back(newstage);
		s = &stages.back();
		s->seq_set = seq;
	}

	// update and return the stage
	for (int i = 0; i < ext.size(); i++)
		s->ext_set.insert(ext[i]);
	s->idx_set.insert(ii);
	s->idx_set.insert(jj);

	return s;
}


/// Returns a (neq x neq) array containing true in each element, if the corresponding
/// block of weak forms is used, and false otherwise.
///
bool **WeakForm::get_blocks() {
	_F_
	bool **blocks = new_matrix<bool>(neq, neq);
	int i, j;
	for (i = 0; i < neq; i++)
		for (j = 0; j < neq; j++)
			blocks[i][j] = false;

	for (i = 0; i < bfvol.size(); i++) {
		blocks[bfvol[i].i][bfvol[i].j] = true;
		if (bfvol[i].sym)
			blocks[bfvol[i].j][bfvol[i].i] = true;
	}

	for (i = 0; i < bfsurf.size(); i++)
		blocks[bfsurf[i].i][bfsurf[i].j] = true;

	return blocks;
}


//// areas /////////////////////////////////////////////////////////////////////////////////////////

int WeakForm::def_area(int n, ...) {
	_F_
	Area newarea;
	va_list ap; va_start(ap, n);
	for (int i = 0; i < n; i++)
		newarea.markers.push_back(va_arg(ap, int));
	va_end(ap);

	areas.push_back(newarea);
	return -areas.size();
}


bool WeakForm::is_in_area_2(int marker, int area) const {
	_F_
	if (-area > areas.size()) ERROR("Invalid area number.");
	const Area *a = &areas[-area - 1];

	for (int i = 0; i < a->markers.size(); i++)
		if (a->markers[i] == marker)
			return true;

	return false;
}
