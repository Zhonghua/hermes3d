#ifndef _NORM_H_
#define _NORM_H_

#include "solution.h"

/// @defgroup norms Norms
///
/// @{

double calc_error(double (*fn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*), MeshFunction *sln1, MeshFunction *sln2);
double calc_norm(double (*fn)(MeshFunction*, RefMap*), MeshFunction *sln);

double h1_error(MeshFunction *sln1, MeshFunction *sln2);
double h1_norm(MeshFunction *sln);

double l2_error(MeshFunction *sln1, MeshFunction *sln2);
double l2_norm(MeshFunction *sln);

double hcurl_error(MeshFunction *sln1, MeshFunction *sln2);
double hcurl_norm(MeshFunction *sln);

double l2_error_hcurl(MeshFunction *sln1, MeshFunction *sln2);
double l2_norm_hcurl(MeshFunction *sln);

#if 0
double l2_error_norm(Solution *sln1, Solution *sln2);
double l2_norm(Solution *sln);
double l2_error_norm_exact(Solution *sln, scalar (*exact)(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz));

double h1_error_norm(Solution *sln1, Solution *sln2);
double h1_norm(Solution *sln);
double h1_error_norm_exact(Solution *sln, scalar (*exact)(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz));
#endif

/// @}

#endif
