// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
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

#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include "common.h"
#include "matrix.h"

/// @defgroup eigensolvers Eigenvalue solvers
///
/// TODO: description


/// Abstract class for encapsulation of a eigenvalue solver
///
///
/// @ingroup eigensolvers
class EigenSolver {
public:
	virtual ~EigenSolver() { }

	/// solve the system
	///
	/// @return true if the system was solved successfully
	virtual bool solve() = 0;
};

#endif
