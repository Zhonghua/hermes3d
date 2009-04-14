// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2008 Miroslav Simko <msimko@miners.utep.edu>
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

#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

#include "../h3dconfig.h"
#include "pardiso.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int pardisoinit_(void *, int *, int *);

extern int
    pardiso_(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *);

#define PARDISOINIT pardisoinit_
#define PARDISO pardiso_

#ifdef __cplusplus
}
#endif

PardisoLinearSolver::PardisoLinearSolver() {
	_F_
#ifdef WITH_PARDISO
#else
	EXIT(ERR_PARDISO_NOT_COMPILED);
#endif
}

PardisoLinearSolver::~PardisoLinearSolver() {
	_F_
#ifdef WITH_PARDISO
#endif
}

bool PardisoLinearSolver::solve() {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return false;
}

