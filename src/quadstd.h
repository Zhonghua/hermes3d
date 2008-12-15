#ifndef _QUAD_STD_H_
#define _QUAD_STD_H_

#include "quad.h"

//
// 1D numeric quadrature //////////////////////////////////////////////////////
//

/// Numerical quadrature for 1D
///
/// @ingroup quadrature
class Quad1DStd : public Quad1D {
public:
	Quad1DStd();
};


//
// 2D numeric quadrature //////////////////////////////////////////////////////
//

/// Numeric quadrature for 2D triangle
///
/// @ingroup quadrature
class QuadStdTri : public Quad2D {
public:
	QuadStdTri();
	~QuadStdTri();
};

/// Numerical quadrature for 2D quadrilateral
///
/// @ingroup quadrature
class QuadStdQuad : public Quad2D {
public:
	QuadStdQuad();
	~QuadStdQuad();
};


//
// 3D numeric quadrature //////////////////////////////////////////////////////
//

/// Numerical quadrature for 3D hexahedron 
///
/// @ingroup quadrature
class QuadStdHex : public Quad3D {
public:
	QuadStdHex();
	~QuadStdHex();

	virtual QuadPt3D *get_points(int order) {
		if (tables[order] == NULL) calc_table(order);
		return tables[order]; 
	}

	virtual QuadPt3D *get_face_points(int face, int order) {
		if (face_tables[face][order] == NULL) calc_face_table(face, order);
		return face_tables[face][order]; 
	}

protected:
	void calc_table(int order);
	void calc_face_table(int face, int order);
	///
	Order3 lower_order_same_accuracy(Order3 ord);
};


/// Numerical quadrature for 3D tetrahedron
///
/// @ingroup quadrature
class QuadStdTetra : public Quad3D {
public:
	QuadStdTetra();
	~QuadStdTetra();
};


/// Numerical quadrature for 3D prisms
///
/// @ingroup quadrature
class QuadStdPrism : public Quad3D {
public:
	QuadStdPrism();
	~QuadStdPrism();
};

#endif
