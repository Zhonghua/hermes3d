// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 David Andrs <dandrs@unr.edu>
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

#include "../h3dconfig.h"
#include "../function.h"
#include "../solution.h"
//#include "../shapeset/h1sinhex.h"
#include "proj.h"
#include <common/callstack.h>

int Projection::vtx_son[27][8] = {
	{ 0, 1, 2, 3, 4, 5, 6, 7 }, // NONE
	{ 0, 0, 0, 0, 0, 0, 0, 0 }, // XYZ
	{ 1, 1, 1, 1, 1, 1, 1, 1 },
	{ 2, 2, 2, 2, 2, 2, 2, 2 },
	{ 3, 3, 3, 3, 3, 3, 3, 3 },
	{ 4, 4, 4, 4, 4, 4, 4, 4 },
	{ 5, 5, 5, 5, 5, 5, 5, 5 },
	{ 6, 6, 6, 6, 6, 6, 6, 6 },
	{ 7, 7, 7, 7, 7, 7, 7, 7 },
	{ 0, 0, 0, 0, 4, 4, 4, 4 },	// XY
	{ 1, 1, 1, 1, 5, 5, 5, 5 },
	{ 2, 2, 2, 2, 6, 6, 6, 6 },
	{ 3, 3, 3, 3, 7, 7, 7, 7 },
	{ 0, 0, 3, 3, 0, 0, 3, 3 }, // XZ
	{ 1, 1, 2, 2, 1, 1, 2, 2 },
	{ 5, 5, 6, 6, 5, 5, 6, 6 },
	{ 4, 4, 7, 7, 4, 4, 7, 7 },
	{ 0, 1, 1, 0, 0, 1, 1, 0 }, // YZ
	{ 3, 2, 2, 3, 3, 2, 2, 3 },
	{ 7, 6, 6, 7, 7, 6, 6, 7 },
	{ 4, 5, 5, 4, 4, 5, 5, 4 },
	{ 0, 0, 3, 3, 4, 4, 7, 7 }, // X
	{ 1, 1, 2, 2, 5, 5, 6, 6 },
	{ 0, 1, 1, 0, 4, 5, 5, 4 }, // Y
	{ 3, 2, 2, 3, 7, 6, 6, 7 },
	{ 0, 1, 2, 3, 0, 1, 2, 3 }, // Z
	{ 4, 5, 6, 7, 4, 5, 6, 7 }
};

int Projection::edge_son[27][12][2] = {
	{ { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } }, // NONE
	{ { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 } }, // XYZ
	{ { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 } },
	{ { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 } },
	{ { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 } },
	{ { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 } },
	{ { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 } },
	{ { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 } },
	{ { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 } },
	{ { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0,-1 }, { 0, 4 }, { 0, 4 }, { 0, 4 }, { 0, 4 }, { 4,-1 }, { 4,-1 }, { 4,-1 }, { 4,-1 } }, // XY
	{ { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1,-1 }, { 1, 5 }, { 1, 5 }, { 1, 5 }, { 1, 5 }, { 5,-1 }, { 5,-1 }, { 5,-1 }, { 5,-1 } },
	{ { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2,-1 }, { 2, 6 }, { 2, 6 }, { 2, 6 }, { 2, 6 }, { 6,-1 }, { 6,-1 }, { 6,-1 }, { 6,-1 } },
	{ { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3,-1 }, { 3, 7 }, { 3, 7 }, { 3, 7 }, { 3, 7 }, { 7,-1 }, { 7,-1 }, { 7,-1 }, { 7,-1 } },
	{ { 0,-1 }, { 0, 3 }, { 3,-1 }, { 3, 0 }, { 0,-1 }, { 0,-1 }, { 3,-1 }, { 3,-1 }, { 0,-1 }, { 0, 3 }, { 3,-1 }, { 3, 0 } }, // XZ
	{ { 1,-1 }, { 1, 2 }, { 2,-1 }, { 2, 1 }, { 1,-1 }, { 1,-1 }, { 2,-1 }, { 2,-1 }, { 1,-1 }, { 1, 2 }, { 2,-1 }, { 2, 1 } },
	{ { 5,-1 }, { 5, 6 }, { 6,-1 }, { 6, 5 }, { 5,-1 }, { 5,-1 }, { 6,-1 }, { 6,-1 }, { 5,-1 }, { 5, 6 }, { 6,-1 }, { 6, 5 } },
	{ { 4,-1 }, { 4, 7 }, { 7,-1 }, { 7, 4 }, { 4,-1 }, { 4,-1 }, { 7,-1 }, { 7,-1 }, { 4,-1 }, { 4, 7 }, { 7,-1 }, { 7, 4 } },
	{ { 0, 1 }, { 1,-1 }, { 1, 0 }, { 0,-1 }, { 0,-1 }, { 1,-1 }, { 1,-1 }, { 0,-1 }, { 0, 1 }, { 1,-1 }, { 1, 0 }, { 0,-1 } }, // YZ
	{ { 3, 2 }, { 2,-1 }, { 2, 3 }, { 3,-1 }, { 3,-1 }, { 2,-1 }, { 2,-1 }, { 3,-1 }, { 3, 2 }, { 2,-1 }, { 2, 3 }, { 3,-1 } },
	{ { 7, 6 }, { 6,-1 }, { 6, 7 }, { 7,-1 }, { 7,-1 }, { 6,-1 }, { 6,-1 }, { 7,-1 }, { 7, 6 }, { 6,-1 }, { 6, 7 }, { 7,-1 } },
	{ { 4, 5 }, { 5,-1 }, { 5, 4 }, { 4,-1 }, { 4,-1 }, { 5,-1 }, { 5,-1 }, { 4,-1 }, { 4, 5 }, { 5,-1 }, { 5, 4 }, { 4,-1 } },
	{ { 0,-1 }, { 0, 3 }, { 3,-1 }, { 3, 0 }, { 0, 4 }, { 0, 4 }, { 3, 7 }, { 3, 7 }, { 4,-1 }, { 4, 7 }, { 7,-1 }, { 7, 4 } }, // X
	{ { 1,-1 }, { 1, 2 }, { 2,-1 }, { 2, 1 }, { 1, 5 }, { 1, 5 }, { 2, 6 }, { 2, 6 }, { 5,-1 }, { 5, 6 }, { 6,-1 }, { 6, 5 } },
	{ { 0, 1 }, { 1,-1 }, { 1, 0 }, { 0,-1 }, { 0, 4 }, { 1, 5 }, { 1, 5 }, { 0, 4 }, { 4, 5 }, { 5,-1 }, { 5, 4 }, { 4,-1 } }, // Y
	{ { 3, 2 }, { 2,-1 }, { 2, 3 }, { 3,-1 }, { 3, 7 }, { 2, 6 }, { 2, 6 }, { 3, 7 }, { 7, 6 }, { 6,-1 }, { 6, 7 }, { 7,-1 } },
	{ { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0,-1 }, { 1,-1 }, { 2,-1 }, { 3,-1 }, { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } }, // Z
	{ { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 4,-1 }, { 5,-1 }, { 6,-1 }, { 7,-1 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } },
};

int Projection::edge_ns[8][12] = {
	{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }, // NONE
	{ 1, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1, 2 }, // X
	{ 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1 }, // Y
	{ 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2 }, // Z
	{ 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1 }, // XY
	{ 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2 }, // XZ
	{ 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 1 }, // YZ
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, // XYZ
};

int Projection::edge_trf[8][12][2] = {
	{ {  0,  1 }, {  1,  2 }, { 2,   3 }, {  3,  0 }, {  0,  4 }, {  1,  5 }, {  2,  6 }, {  3,  7 }, {  4,  5 }, {  5,  6 }, {  6,  7 }, {  7,  4 } }, // NONE
	{ { 16, -1 }, { 16, 17 }, { 17, -1 }, { 17, 16 }, { 16, 19 }, { 16, 19 }, { 17, 18 }, { 17, 18 }, { 19, -1 }, { 19, 18 }, { 18, -1 }, { 18, 19 } }, // X
	{ { 12, 13 }, { 13, -1 }, { 13, 12 }, { 12, -1 }, { 12, 15 }, { 13, 14 }, { 13, 14 }, { 12, 15 }, { 15, 14 }, { 14, -1 }, { 14, 15 }, { 15, -1 } }, // Y
	{ {  8,  9 }, {  9, 10 }, { 10, 11 }, { 11,  8 }, {  8, -1 }, {  9, -1 }, { 10, -1 }, { 11, -1 }, {  8,  9 }, {  9, 10 }, { 10, 11 }, { 11,  8 } }, // Z
	{ { 24, -1 }, { 24, -1 }, { 24, -1 }, { 24, -1 }, { 24, 25 }, { 24, 25 }, { 24, 25 }, { 24, 25 }, { 25, -1 }, { 25, -1 }, { 25, -1 }, { 25, -1 } }, // XY
	{ { 22, -1 }, { 22, 23 }, { 23, -1 }, { 23, 22 }, { 22, -1 }, { 22, -1 }, { 23, -1 }, { 23, -1 }, { 22, -1 }, { 22, 23 }, { 23, -1 }, { 23, 22 } }, // XZ
	{ { 20, 21 }, { 21, -1 }, { 21, 20 }, { 20, -1 }, { 20, -1 }, { 21, -1 }, { 21, -1 }, { 20, -1 }, { 20, 21 }, { 21, -1 }, { 21, 20 }, { 20, -1 } }, // YZ
	{ { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 } }, // XYZ
};

int Projection::face_son[27][6][4] = {
	{ { 0, 3, 7, 4 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 }, { 3, 2, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } }, // NONE
	{ { 0,-1,-1,-1 }, { 0,-1,-1,-1 }, { 0,-1,-1,-1 }, { 0,-1,-1,-1 }, { 0,-1,-1,-1 }, { 0,-1,-1,-1 } }, // XYZ
	{ { 1,-1,-1,-1 }, { 1,-1,-1,-1 }, { 1,-1,-1,-1 }, { 1,-1,-1,-1 }, { 1,-1,-1,-1 }, { 1,-1,-1,-1 } },
	{ { 2,-1,-1,-1 }, { 2,-1,-1,-1 }, { 2,-1,-1,-1 }, { 2,-1,-1,-1 }, { 2,-1,-1,-1 }, { 2,-1,-1,-1 } },
	{ { 3,-1,-1,-1 }, { 3,-1,-1,-1 }, { 3,-1,-1,-1 }, { 3,-1,-1,-1 }, { 3,-1,-1,-1 }, { 3,-1,-1,-1 } },
	{ { 4,-1,-1,-1 }, { 4,-1,-1,-1 }, { 4,-1,-1,-1 }, { 4,-1,-1,-1 }, { 4,-1,-1,-1 }, { 4,-1,-1,-1 } },
	{ { 5,-1,-1,-1 }, { 5,-1,-1,-1 }, { 5,-1,-1,-1 }, { 5,-1,-1,-1 }, { 5,-1,-1,-1 }, { 5,-1,-1,-1 } },
	{ { 6,-1,-1,-1 }, { 6,-1,-1,-1 }, { 6,-1,-1,-1 }, { 6,-1,-1,-1 }, { 6,-1,-1,-1 }, { 6,-1,-1,-1 } },
	{ { 7,-1,-1,-1 }, { 7,-1,-1,-1 }, { 7,-1,-1,-1 }, { 7,-1,-1,-1 }, { 7,-1,-1,-1 }, { 7,-1,-1,-1 } },
	{ { 0, 4,-1,-1 }, { 0, 4,-1,-1 }, { 0, 4,-1,-1 }, { 0, 4,-1,-1 }, { 0,-1,-1,-1 }, { 4,-1,-1,-1 } }, // XY
	{ { 1, 5,-1,-1 }, { 1, 5,-1,-1 }, { 1, 5,-1,-1 }, { 1, 5,-1,-1 }, { 1,-1,-1,-1 }, { 5,-1,-1,-1 } },
	{ { 2, 6,-1,-1 }, { 2, 6,-1,-1 }, { 2, 6,-1,-1 }, { 2, 6,-1,-1 }, { 2,-1,-1,-1 }, { 6,-1,-1,-1 } },
	{ { 3, 7,-1,-1 }, { 3, 7,-1,-1 }, { 3, 7,-1,-1 }, { 3, 7,-1,-1 }, { 3,-1,-1,-1 }, { 7,-1,-1,-1 } },
	{ { 0, 3,-1,-1 }, { 0, 3,-1,-1 }, { 0,-1,-1,-1 }, { 3,-1,-1,-1 }, { 0, 3,-1,-1 }, { 0, 3,-1,-1 } }, // XZ
	{ { 1, 2,-1,-1 }, { 1, 2,-1,-1 }, { 1,-1,-1,-1 }, { 2,-1,-1,-1 }, { 1, 2,-1,-1 }, { 1, 2,-1,-1 } },
	{ { 5, 6,-1,-1 }, { 5, 6,-1,-1 }, { 5,-1,-1,-1 }, { 6,-1,-1,-1 }, { 5, 6,-1,-1 }, { 5, 6,-1,-1 } },
	{ { 4, 7,-1,-1 }, { 4, 7,-1,-1 }, { 4,-1,-1,-1 }, { 7,-1,-1,-1 }, { 4, 7,-1,-1 }, { 4, 7,-1,-1 } },
	{ { 0,-1,-1,-1 }, { 1,-1,-1,-1 }, { 0, 1,-1,-1 }, { 0, 1,-1,-1 }, { 0, 1,-1,-1 }, { 0, 1,-1,-1 } }, // YZ
	{ { 3,-1,-1,-1 }, { 2,-1,-1,-1 }, { 3, 2,-1,-1 }, { 3, 2,-1,-1 }, { 3, 2,-1,-1 }, { 3, 2,-1,-1 } },
	{ { 7,-1,-1,-1 }, { 6,-1,-1,-1 }, { 7, 6,-1,-1 }, { 7, 6,-1,-1 }, { 7, 6,-1,-1 }, { 7, 6,-1,-1 } },
	{ { 4,-1,-1,-1 }, { 5,-1,-1,-1 }, { 4, 5,-1,-1 }, { 4, 5,-1,-1 }, { 4, 5,-1,-1 }, { 4, 5,-1,-1 } },
	{ { 0, 3, 7, 4 }, { 0, 3, 7, 4 }, { 0, 4,-1,-1 }, { 3, 7,-1,-1 }, { 0, 3,-1,-1 }, { 4, 7,-1,-1 } }, // X
	{ { 1, 2, 6, 5 }, { 1, 2, 6, 5 }, { 1, 5,-1,-1 }, { 2, 6,-1,-1 }, { 1, 2,-1,-1 }, { 5, 6,-1,-1 } },
	{ { 0, 4,-1,-1 }, { 1, 5,-1,-1 }, { 0, 1, 5, 4 }, { 0, 1, 5, 4 }, { 0, 1,-1,-1 }, { 4, 5,-1,-1 } }, // Y
	{ { 3, 7,-1,-1 }, { 2, 6,-1,-1 }, { 3, 2, 6, 7 }, { 3, 2, 6, 7 }, { 3, 2,-1,-1 }, { 7, 6,-1,-1 } },
	{ { 0, 3,-1,-1 }, { 1, 2,-1,-1 }, { 0, 1,-1,-1 }, { 3, 2,-1,-1 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3 } }, // Z
	{ { 4, 7,-1,-1 }, { 5, 6,-1,-1 }, { 4, 5,-1,-1 }, { 7, 6,-1,-1 }, { 4, 5, 6, 7 }, { 4, 5, 6, 7 } },
};

int Projection::face_trf[8][6][4] = {
	{ {  0,  3,  7,  4 }, {  1,  2,  6,  5 }, {  0,  1,  5,  4 }, {  3,  2,  6,  7 }, {  0,  1,  2,  3 }, {  4,  5,  6,  7 } }, // NONE
	{ { 16, 17, 18, 19 }, { 16, 17, 18, 19 }, { 16, 19, -1, -1 }, { 17, 18, -1, -1 }, { 16, 17, -1, -1 }, { 19, 18, -1, -1 } }, // X
	{ { 12, 15, -1, -1 }, { 13, 14, -1, -1 }, { 12, 13, 14, 15 }, { 12, 13, 14, 15 }, { 12, 13, -1, -1 }, { 15, 14, -1, -1 } }, // Y
	{ {  8, 11, -1, -1 }, {  9, 10, -1, -1 }, {  8,  9, -1, -1 }, { 11, 10, -1, -1 }, {  8,  9, 10, 11 }, {  8,  9, 10, 11 } }, // Z
	{ { 24, 25, -1, -1 }, { 24, 25, -1, -1 }, { 24, 25, -1, -1 }, { 24, 25, -1, -1 }, { 24, -1, -1, -1 }, { 25, -1, -1, -1 } }, // XY
	{ { 22, 23, -1, -1 }, { 22, 23, -1, -1 }, { 22, -1, -1, -1 }, { 23, -1, -1, -1 }, { 22, 23, -1, -1 }, { 22, 23, -1, -1 } }, // XZ
	{ { 20, -1, -1, -1 }, { 21, -1, -1, -1 }, { 20, 21, -1, -1 }, { 20, 21, -1, -1 }, { 20, 21, -1, -1 }, { 20, 21, -1, -1 } }, // YZ
	{ { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 } }, // XYZ
};

int Projection::face_ns[8][6] = {
	{ 4, 4, 4, 4, 4, 4 },	// NONE
	{ 4, 4, 2, 2, 2, 2 },	// X
	{ 2, 2, 4, 4, 2, 2 },	// Y
	{ 2, 2, 2, 2, 4, 4 },	// Z
	{ 2, 2, 2, 2, 1, 1 },	// XY
	{ 2, 2, 1, 1, 2, 2 },	// XZ
	{ 1, 1, 2, 2, 2, 2 },	// YZ
	{ 1, 1, 1, 1, 1, 1 },	// XYZ
};


int Projection::int_son[27][8] = {
	{ 0,  1,  2,  3,  4,  5,  6,  7 }, // NONE
	{ 0, -1, -1, -1, -1, -1, -1, -1 }, // XYZ
	{ 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 0,  4, -1, -1, -1, -1, -1, -1 }, // XY
	{ 1,  5, -1, -1, -1, -1, -1, -1 },
	{ 2,  6, -1, -1, -1, -1, -1, -1 },
	{ 3,  7, -1, -1, -1, -1, -1, -1 },
	{ 0,  3, -1, -1, -1, -1, -1, -1 }, // XZ
	{ 1,  2, -1, -1, -1, -1, -1, -1 },
	{ 5,  6, -1, -1, -1, -1, -1, -1 },
	{ 4,  7, -1, -1, -1, -1, -1, -1 },
	{ 0,  1, -1, -1, -1, -1, -1, -1 }, // YZ
	{ 3,  2, -1, -1, -1, -1, -1, -1 },
	{ 7,  6, -1, -1, -1, -1, -1, -1 },
	{ 4,  5, -1, -1, -1, -1, -1, -1 },
	{ 0,  3,  7,  4, -1, -1, -1, -1 }, // X
	{ 1,  2,  6,  5, -1, -1, -1, -1 },
	{ 0,  1,  5,  4, -1, -1, -1, -1 }, // Y
	{ 3,  2,  6,  7, -1, -1, -1, -1 },
	{ 0,  1,  2,  3, -1, -1, -1, -1 }, // Z
	{ 4,  5,  6,  7, -1, -1, -1, -1 },
};

int Projection::int_ns[8] = { 8, 4, 4, 4, 2, 2, 2, 1 };

int Projection::int_trf[8][8] = {
	{  0,  1,  2,  3,  4,  5,  6,  7 }, // NONE
	{ 16, 17, 18, 19, -1, -1, -1, -1 }, // X
	{ 12, 13, 14, 15, -1, -1, -1, -1 }, // Y
	{  8,  9, 10, 11, -1, -1, -1, -1 }, // Z
	{ 24, 25, -1, -1, -1, -1, -1, -1 }, // XY
	{ 22, 23, -1, -1, -1, -1, -1, -1 }, // XZ
	{ 20, 21, -1, -1, -1, -1, -1, -1 }, // YZ
	{ -1, -1, -1, -1, -1, -1, -1, -1 }, // XYZ
};


Trf *Projection::get_trf(int trf) {
	_F_
	static Trf none = { { 1, 1, 1 }, { 0, 0, 0 } };
	if (trf == -1) return &none;
	else return Transformable::hex_trf + trf;			// FIXME: HEX-specific
}

Projection::Projection(Solution *afn, Element *e, Shapeset *ss) {
	_F_
	this->sln = afn;
	this->ss = ss;

	mesh = sln->get_mesh();
	// TODO: check that the element 'e' is not active and has 8 sons
	base_elem = mesh->elements[e->id];

	quad = get_quadrature(e->get_mode());

	fu = new ShapeFunction(ss);
	fv = new ShapeFunction(ss);

	fu->set_active_element(base_elem);
	fv->set_active_element(base_elem);

	// null
	vertex_proj = NULL;
	for (int i = 0; i < Hex::NUM_EDGES; i++) edge_proj[i] = NULL;
	for (int i = 0; i < Hex::NUM_FACES; i++) face_proj[i] = NULL;
	bubble_proj = NULL;
	proj = NULL;
	proj_fns = 0;
}

Projection::~Projection() {
	_F_
	delete fu;
	delete fv;

	free_proj();
}

void Projection::free_proj() {
	_F_
	delete [] vertex_proj;
	for (int i = 0; i < Hex::NUM_EDGES; i++) delete [] edge_proj[i];
	for (int i = 0; i < Hex::NUM_FACES; i++) delete [] face_proj[i];
	delete [] bubble_proj;

	delete [] proj;
}

void Projection::calc_projection(int split, int son, order3_t order) {
	_F_
	free_proj();
	calc_vertex_proj(split, son + 1);
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
		calc_edge_proj(iedge, split, son + 1, order);
	for (int iface = 0; iface < Hex::NUM_FACES; iface++)
		calc_face_proj(iface, split, son + 1, order);
	calc_bubble_proj(split, son + 1, order);

	proj_fns = (order.x + 1) * (order.y + 1) * (order.z + 1);
	proj = new ProjItem *[proj_fns];

	int m = 0;
	for (int i = 0; i < Hex::NUM_VERTICES; i++, m++)
		proj[m] = vertex_proj + i;
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		order1_t edge_order = order.get_edge_order(iedge);
		int edge_fns = edge_order - 1;
		for (int i = 0; i < edge_fns; i++, m++)
			proj[m] = edge_proj[iedge] + i;
	}
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		order2_t face_order = order.get_face_order(iface);
		int face_fns = (face_order.x - 1) * (face_order.y - 1);
		for (int i = 0; i < face_fns; i++, m++)
			proj[m] = face_proj[iface] + i;
	}
	int bubble_fns = (order.x - 1) * (order.y - 1) * (order.z - 1);
	for (int i = 0; i < bubble_fns; i++, m++)
		proj[m] = bubble_proj + i;
}
