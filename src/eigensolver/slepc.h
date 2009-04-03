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


#ifndef _SLEPC_H_
#define _SLEPC_H_

#include "../h3dconfig.h"
#include "../eigensolver.h"
#include "../linsolver/petsc.h"

// define an alias for Slepc matrix
#define SlepcMatrix		PetscMatrix


class SlepcEigenSolver : public EigenSolver {
public:
	SlepcEigenSolver(SlepcMatrix &a, SlepcMatrix &b);
	virtual ~SlepcEigenSolver();

	virtual bool solve();

protected:
	SlepcMatrix &a;
	SlepcMatrix &b;
};

#endif /* _SLEPC_H_ */
