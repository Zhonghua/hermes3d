//
// General include file for HERMES3D
//

#ifndef _HERMES_3D_
#define _HERMES_3D_

#include "common/array.h"
#include "common/bitarray.h"
#include "common/map.h"
#include "common/maphs.h"
#include "common/mapord.h"
#include "common/trace.h"
#include "common/utils.h"

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
