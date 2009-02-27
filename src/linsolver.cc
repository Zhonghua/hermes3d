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

#include "config.h"
#include "linsolver.h"
#include <common/error.h>

using namespace std;

int LinearSolver::sort_and_store_indices(Page *page, int *buffer, int *max) {
	// gather all pages in the buffer, deleting them along the way
	int *end = buffer;
	while (page != NULL) {
		memcpy(end, page->idx, sizeof(int) * page->count);
		end += page->count;
		Page *tmp = page;
		page = page->next;
		delete tmp;
	}

	// sort the indices and remove duplicities
	qsort_int(buffer, end - buffer);
	int *q = buffer;
	for (int *p = buffer, last = -1; p < end; p++)
		if (*p != last)
			*q++ = last = *p;

	return q - buffer;
}

int LinearSolver::get_num_indices(Page **pages, int ndofs) {
	int total = 0;
	for (int i = 0; i < ndofs; i++)
		for (Page *page = pages[i]; page != NULL; page = page->next)
			total += page->count;

	return total;
}

void LinearSolver::insert_value(int *Ai, scalar *Ax, int Alen, int idx, scalar value) {
	if (idx >= 0) {
		register int lo = 0, hi = Alen - 1, mid;

		while (1) {
			mid = (lo + hi) >> 1;

			if (idx < Ai[mid]) hi = mid - 1;
			else if (idx > Ai[mid]) lo = mid + 1;
			else break;

			if (lo > hi) EXIT(ERR_FAILURE, "Sparse matrix entry not found.");
		}

		Ax[mid] += value;
	}
}
