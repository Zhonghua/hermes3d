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

/// @}

#endif
