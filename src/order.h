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

#ifndef _ORDER_H_
#define _ORDER_H_

#include "common.h"
#include <common/error.h>


// maximal order of quadratures for 1D
#define MAX_QUAD_ORDER								24
// maximal order of quadratures for triangle
#define MAX_QUAD_ORDER_TRI							20
// maximal order of quadratures for tetra
#define MAX_QUAD_ORDER_TETRA						20


// orders

// 1D polynomial order
typedef
	int order1_t;


// 2D polynomial order
struct order2_t {
	order2_t() { type = 3; order = 31; }
	order2_t(int order) { type = MODE_TRIANGLE; this->order = order; }
	order2_t(int x, int y) { type = MODE_QUAD; this->x = x; this->y = y; }

	unsigned type:2;		// EMode2D
	union {
		// tri
		struct {
			unsigned order:5;
		};
		// quad
		struct {
			// directional orders
			unsigned x:5;
			unsigned y:5;
		};
	};

	bool invalid() { return (type == 3); }

	// Operators

	order2_t operator+(const order2_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	return order2_t(this->order + o.order);
			case MODE_QUAD: return order2_t(this->x + o.x, this->y + o.y);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order2_t(-1);
	}

	order2_t operator+=(const order2_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	this->order += o.order; break;
			case MODE_QUAD: this->x += o.x; this->y += o.y; break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	order2_t operator*(const int c) {
		switch (type) {
			case MODE_TRIANGLE: return order2_t(c * this->order);
			case MODE_QUAD: return order2_t(c * this->x, c * this->y);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order2_t(-1);
	}

	order2_t operator*(const order2_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	return order2_t(this->order * o.order);
			case MODE_QUAD: return order2_t(this->x * o.x, this->y * o.y);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order2_t(-1);
	}

	order2_t operator*=(const int c) {
		switch (type) {
			case MODE_TRIANGLE:	this->order *= c; break;
			case MODE_QUAD: this->x *= c; this->y *= c; break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	order2_t operator*=(const order2_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	this->order *= o.order; break;
			case MODE_QUAD: this->x *= o.x; this->y *= o.y; break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	// relation operators
	bool operator==(const order2_t &o) {
		if (this->type != o.type) return false;
		switch (this->type) {
			case MODE_TRIANGLE: return this->order == o.order;
			case MODE_QUAD: return (this->x == o.x) && (this->y == o.y);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	bool operator!=(const order2_t &o) {
		if (this->type != o.type) return true;
		switch (this->type) {
			case MODE_TRIANGLE: return this->order != o.order;
			case MODE_QUAD: return (this->x != o.x) || (this->y != o.y);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	const char *str() {
		static char s[64];
		switch (type) {
			case MODE_TRIANGLE: sprintf(s, "(%d)", this->order); break;
			case MODE_QUAD: sprintf(s, "(%d, %d)", this->x, this->y); break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return s;
	}

	int get_idx() {
		switch (type) {
			case MODE_TRIANGLE: return (this->type << 10) | this->order;
			case MODE_QUAD: return (((this->type << 5) | this->y) << 5) | this->x;
			default: assert(false); EXIT(ERR_UNKNOWN_MODE); break;
		}
		return -1;
	}

	static order2_t from_int(int o) {
		int type = (o >> 10) & 0x03;
		switch (type) {
			case MODE_TRIANGLE: return order2_t(o & 0x1F); break;
			case MODE_QUAD: return order2_t(o & 0x1F, (o >> 5) & 0x1F); break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order2_t(-1);
	}
};

inline order2_t max(order2_t a, order2_t b) {
	assert(a.type == b.type);
	switch (a.type) {
		case MODE_TRIANGLE: return order2_t(std::max(a.order, b.order));
		case MODE_QUAD: return order2_t(std::max(a.x, b.x), std::max(a.y, b.y));
		default: EXIT(ERR_UNKNOWN_MODE); break;
	}
	return order2_t(-1);
}


// 3D polynomial order
//
// all 1s mean invalid (not set) - see default ctor
//
struct order3_t {
	order3_t() { type = 7; }
	order3_t(int order) { type = MODE_TETRAHEDRON; this->order = order; }
	order3_t(int x, int y, int z) { type = MODE_HEXAHEDRON; this->x = x; this->y = y; this->z = z; }

	unsigned type:3;		// EMode3D
	union {
		// tetra
		struct {
			unsigned order: 15;
		};
		// hex
		struct {
			// directional orders
			unsigned x:5;
			unsigned y:5;
			unsigned z:5;
		};
	};

	bool invalid() { return (type == 7); }

	// Operators

	order3_t operator+(const order3_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON:	return order3_t(this->order + o.order);
			case MODE_HEXAHEDRON: return order3_t(this->x + o.x, this->y + o.y, this->z + o.z);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order3_t(-1);
	}

	order3_t operator+=(const order3_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON:	this->order += o.order; break;
			case MODE_HEXAHEDRON: this->x += o.x; this->y += o.y; this->z += o.z; break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	order3_t operator*(const int c) {
		switch (type) {
			case MODE_TETRAHEDRON:	return order3_t(c * this->order);
			case MODE_HEXAHEDRON: return order3_t(c * this->x, c * this->y, c * this->z);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order3_t(-1);
	}

	order3_t operator*(const order3_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON:	return order3_t(this->order * o.order);
			case MODE_HEXAHEDRON: return order3_t(this->x * o.x, this->y * o.y, this->z * o.z);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order3_t(-1);
	}

	order3_t operator*=(const int c) {
		switch (type) {
			case MODE_TETRAHEDRON:	this->order *= c; break;
			case MODE_HEXAHEDRON: this->x *= c; this->y *= c; this->z *= c; break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	order3_t operator*=(const order3_t &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON: this->order *= o.order; break;
			case MODE_HEXAHEDRON: this->x *= o.x; this->y *= o.y; this->z *= o.z; break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	// relation operators
	bool operator==(const order3_t &o) {
		if (this->type != o.type) return false;
		switch (this->type) {
			case MODE_TETRAHEDRON: return this->order == o.order;
			case MODE_HEXAHEDRON: return (this->x == o.x) && (this->y == o.y) && (this->z == o.z);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	bool operator!=(const order3_t &o) {
		if (this->type != o.type) return true;
		switch (this->type) {
			case MODE_TETRAHEDRON: return this->order != o.order;
			case MODE_HEXAHEDRON: return (this->x != o.x) || (this->y != o.y) || (this->z != o.z);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	const char *str() {
		static char s[64];
		switch (type) {
			case MODE_TETRAHEDRON: sprintf(s, "(%d)", this->order); break;
			case MODE_HEXAHEDRON: sprintf(s, "(%d, %d, %d)", this->x, this->y, this->z); break;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return s;
	}

	int get_idx() {
		assert(!invalid());
		switch (type) {
			case MODE_TETRAHEDRON: return ((this->type) << 15) | this->order;
			case MODE_HEXAHEDRON: return (((((this->type << 5) | this->z) << 5) | this->y) << 5) | this->x;
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return -1;
	}

	static order3_t from_int(int o) {
		int type = (o >> 15) & 0x07;
		switch (type) {
			case MODE_TETRAHEDRON: return order3_t(o & 0x7FFF);
			case MODE_HEXAHEDRON: return order3_t(o & 0x1F, (o >> 5) & 0x1F, (o >> 10) & 0x1F);
			default: EXIT(ERR_UNKNOWN_MODE); break;
		}
		return order3_t(-1);
	}

	order1_t get_edge_order(int edge) {
		switch (type) {
			case MODE_TETRAHEDRON: return this->order;
			case MODE_HEXAHEDRON:
				if ((edge == 0) || (edge == 2) || (edge == 10) || (edge == 8)) return this->x;
				else if((edge == 1) || (edge == 3) || (edge == 11) || (edge == 9)) return this->y;
				else if((edge == 4) || (edge == 5) || (edge == 6) || (edge == 7)) return this->z;
				else EXIT(ERR_EDGE_INDEX_OUT_OF_RANGE);
			default:
				EXIT(ERR_UNKNOWN_MODE);
				break;
		}
		return -1;
	}

	order2_t get_face_order(int face) {
		switch (type) {
			case MODE_TETRAHEDRON: return this->order;
			case MODE_HEXAHEDRON:
				if ((face == 0) || (face == 1)) return order2_t(this->y, this->z);
				else if ((face == 2) || (face == 3)) return order2_t(this->x, this->z);
				else if ((face == 4) || (face == 5)) return order2_t(this->x, this->y);
				else EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
			default:
				EXIT(ERR_UNKNOWN_MODE);
				break;
		}
		return order2_t(-1);
	}

	void limit() {
#ifndef DEBUG_ORDER
		switch (type) {
			case MODE_TETRAHEDRON:
				if (this->order > MAX_QUAD_ORDER_TETRA) this->order = MAX_QUAD_ORDER_TETRA;
				break;

			case MODE_HEXAHEDRON:
				if (this->x > MAX_QUAD_ORDER) this->x = MAX_QUAD_ORDER;
				if (this->y > MAX_QUAD_ORDER) this->y = MAX_QUAD_ORDER;
				if (this->z > MAX_QUAD_ORDER) this->z = MAX_QUAD_ORDER;
				break;

			default:
				EXIT(ERR_UNKNOWN_MODE);
				break;
		}
#else
		set_maximal();
#endif
	}

	void set_maximal() {
		switch (type) {
			case MODE_TETRAHEDRON: this->order = MAX_QUAD_ORDER_TETRA; break;

			case MODE_HEXAHEDRON:
				this->x = MAX_QUAD_ORDER;
				this->y = MAX_QUAD_ORDER;
				this->z = MAX_QUAD_ORDER;
				break;

			default:
				EXIT(ERR_UNKNOWN_MODE);
				break;
		}
	}
};

inline order3_t max(order3_t a, order3_t b) {
	assert(a.type == b.type);
	switch (a.type) {
		case MODE_TETRAHEDRON: return order3_t(std::max(a.order, b.order));
		case MODE_HEXAHEDRON: return order3_t(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
		default: EXIT(ERR_UNKNOWN_MODE); break;
	}
	return order3_t(-1);
}

inline order3_t turn_hex_face_order(order3_t ord) {
	int o1 = ord.x;
	int o2 = ord.y;
	int o3 = ord.z;
	if (o1 <= 1) std::swap(o2, o3);
	if (o2 <= 1) std::swap(o1, o3);
	if (o3 <= 1) std::swap(o1, o2);
	return order3_t(o1, o2, o3);
}

#endif
