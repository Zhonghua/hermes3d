// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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

#ifndef _MAPHS_H_
#define _MAPHS_H_

#include <assert.h>

/// \file maphs.h

#include <Judy.h>

#ifndef INVALID_IDX
#define INVALID_IDX					((Word_t) -1)
#endif

/// Implementation of a dynamic array, using an array-of-bytes of Length as an Index and a word as a Value.
/// C++ encapsulation of JudyHS functions.
class MapHS {
protected:
	void *judy;

public:
	MapHS() {
		judy = NULL;
	}

	virtual ~MapHS() {
		remove_all();
	}

	/// Lookup for item with key \c key.
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Size of \c key.
	/// \param[in] item Item to insert
	/// \return
	/// 	\li true if the key exists, item then contains the value,
	/// 	\li false otherwise
	bool lookup(char *key, int length, Word_t &item) const {
		void *pval;
		// check if the key exists
		JHSG(pval, judy, key, length);
		if (pval == NULL) {
			return false;
		}
		else {
			item = *(Word_t *) pval;
			return true;
		}
	}

	/// Add a new (key, item) pair into the map
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Size of \c key.
	/// \param[in] item Item to insert
	bool set(char *key, int length, Word_t item) {
		void *pval;
		JHSG(pval, judy, key, item);
		if (pval == NULL) {
			// insert new item
			JHSI(pval, judy, key, length);
			if (pval == NULL) return false;
		}
		*(Word_t *) pval = item;
		return true;
	}

	/// Delete an item with key \c key from the map.
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Size of \c key.
	bool remove(char *key, int length) {
		void *pval;
		JHSG(pval, judy, key, length);
		if (pval == NULL) {
			return false;
		}
		else {
			int rc;
			JHSD(rc, judy, key, length);
			return (rc == 1);
		}
	}

	/// Remove all items from the array
	void remove_all() {
		int val;
		JHSFA(val, judy);
	}
};


#endif
