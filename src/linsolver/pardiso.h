// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandr@unr.edu>
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


#ifndef PARDISOSOLVER_H_
#define PARDISOSOLVER_H_

#include "../linsolver.h"

/// Encapsulation of pardiso linear solver
///
/// @ingroup linearsolvers
class PardisoLinearSolver : public LinearSolver {
public:
	PardisoLinearSolver();
	virtual ~PardisoLinearSolver();

	virtual bool solve();
protected:
};

#endif /*PARDISOSOLVER_H_*/
