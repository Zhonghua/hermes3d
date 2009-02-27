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

#ifndef _COMMON_H_
#define _COMMON_H_

// common headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <cstdarg>
#include <algorithm>			// std::min, std::max

#ifdef COMPLEX

#include <complex>

typedef std::complex<double> complex;
typedef complex scalar;
typedef complex complex2[2];
#define CONJ(a)				(std::conj(a))
#define REAL(a)				(std::real(a))
#define IMAG(a)				(std::imag(a))
#define ABS(a)				(std::abs(a))
#define SCALAR_FMT			"(%lf, %lf)"
#define SCALAR(a)			std::real(a), std::imag(a)

#else

typedef double scalar;
#define CONJ(a)				(a)
#define REAL(a)				(a)
#define IMAG(a)				(0)
#define ABS(a)				(fabs(a))
#define SCALAR_FMT			"%lf"
#define SCALAR(a)			(a)

#endif


// 1D element modes
enum EMode1D {
	MODE_LINE = 0
};

// 2D element modes
enum EMode2D {
	MODE_TRIANGLE = 0,
	MODE_QUAD = 1
};

// 3D element modes
enum EMode3D {
	MODE_TETRAHEDRON = 0,
	MODE_HEXAHEDRON = 1,
	MODE_PRISM = 2
};


// points
struct Point1D {
	double x;			// coordinates of a point
};

struct Point2D {
	double x, y;		// coordinates of a point
};

struct Point3D {
	double x, y, z;		// coordinates of a point
};

inline double dot_product(const Point3D &a, const Point3D &b) { return a.x * b.x + a.y * b.y + a.z * b.z;}

inline Point3D cross_product(Point3D a, Point3D b) {
	Point3D r = {
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	};
	return r;
}

inline Point3D lin_comb(Point3D a, double coef_a, Point3D b, double coef_b) {
	Point3D r = {
		coef_a * a.x + coef_b * b.x,
		coef_a * a.y + coef_b * b.y,
		coef_a * a.z + coef_b * b.z,
	};
	return r;
}

inline double norm(const Point3D &pt) { return sqrt(dot_product(pt, pt)); }

inline Point3D normalize(const Point3D &pt) {
	double n = norm(pt);
	Point3D res = { pt.x / n, pt.y / n, pt.z / n };
	return res;
}


struct Vector3D {
	scalar x, y, z;		// coordinates of a point

	Vector3D() {
		x = y = z = 0;
	}

	Vector3D(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void set(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	scalar dot_product(Vector3D vec2) {return x * vec2.x + y * vec2.y + z * vec2.z;};
	scalar dot_product(Point3D vec2) {return x * vec2.x + y * vec2.y + z * vec2.z;};
	double norm() { return sqrt(REAL(dot_product(*this)));};
	void cross_product(Vector3D a, Vector3D b) {
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
	}
	void cross_product(Point3D a, Vector3D b) {
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
	}
	void cross_product(Vector3D a, Point3D b) {
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
	}

	void normalize(){
		double n = norm();
		x /= n;
		y /= n;
		z /= n;
	}

	void subtract(Vector3D b){
		x -= b.x;
		y -= b.y;
		z -= b.z;
	}
};



typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double2x2[2][2];
typedef double double3x3[3][3];
typedef int int2[2];
typedef scalar scalar3[3];
typedef unsigned long long int uint64;

// maximal polynomial order of elements
#define MAX_ELEMENT_ORDER							10

// Dirichlet lift is a special DOF with nubmer -1
#define DIRICHLET_DOF								-1

//

inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }

#ifdef COMPLEX
inline double sqr(complex x) { return std::norm(x); }
#endif


#endif
