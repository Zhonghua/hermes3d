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

#ifndef _LINSOLVER_H_
#define _LINSOLVER_H_

#include "common.h"

/// @defgroup linearsolvers Linear solvers
///
/// TODO: description


/// Abstract class for encapsulation of a linear solver
///
///
/// TODO: Adjust interface to support faster update of matrix and rhs
///
/// @ingroup linearsolvers
class LinearSolver {
public:
	LinearSolver() { sln = NULL; }
	virtual ~LinearSolver() { delete [] sln; }

	virtual bool solve() = 0;
	scalar *get_solution() { return sln; }

protected:
	scalar *sln;
};

#endif
