// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2007 - 2009 Pavel Kus <pavel.kus@gmail.com>
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

#ifndef _HERMES_3D_
#define _HERMES_3D_

#include "common/array.h"
#include "common/bitarray.h"
#include "common/map.h"
#include "common/maphs.h"
#include "common/mapord.h"
#include "common/trace.h"
#include "common/utils.h"
#include "common/timer.h"

#include "common.h"

// mesh
#include "mesh.h"

// mesh loaders
#include "meshloader.h"
#include "loader/mesh3d.h"
#include "loader/hdf5.h"

// spaces
#include "space.h"
#include "space/h1.h"
#include "space/hcurl.h"

#include "order.h"
// quadrature
#include "quad.h"
#include "quadstd.h"

#include "refmap.h"
#include "integrals/h1.h"
#include "integrals/hcurl.h"

#include "refdomain.h"

#include "precalc.h"

// shapesets
#include "shapeset.h"
#include "shapeset/common.h"
#include "shapeset/h1lobattotetra.h"
#include "shapeset/h1lobattohex.h"
#include "shapeset/hcurllobattohex.h"

#include "norm.h"

// output
#include "output.h"
#include "output/gmsh.h"
#include "output/vtk.h"
#include "output/graph.h"

#include "asmlist.h"
#include "discretization.h"
#include "solution.h"
#include "filter.h"

// linear solvers
#include "linsolver.h"
#include "linsolver/umfpack.h"
#include "linsolver/pardiso.h"
#include "linsolver/petsc.h"

#endif
