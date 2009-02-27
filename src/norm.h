// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2009 Pavel Kus <pavel.kus@gmail.com>
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
