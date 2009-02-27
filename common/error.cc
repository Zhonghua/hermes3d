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

//
// error.cc
//

#include "error.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

static const char *h_str_error[] = {
	"Out of memory.",
	"Not enough parameters.",
	"Can not open file '%s'.",
	NULL,							// mesh problem
	"Not yet implemened.",
	"Unknown mode (mode = %d).",
	"Hermes3D was not built with PETSc support.",
	"Hermes3D was not built with UMFPACK support.",
	"Hermes3D was not built with HDF5 support.",
	"Hermes3D was not built with MPI support.",
	"Face index out of range.",
	"Edge index out of range.",
	"Hermes3D was not built with tetra elements.",
	"Hermes3D was not built with hex elements.",
	"Hermes3D was not built with prism elements.",
	"Hermes3D was not built with Pardiso support.",
	"Unknown refinement type (refinement = %d)."
};

// errors

void h_exit(int err_code, int line, const char *func, const char *file, char const *fmt, ...) {
	fprintf(stderr, "ERROR: %s:%d: %s: ", file, line, func);
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, "\n");

	if (err_code == -1)
		exit(err_code);
	else
		exit(ERR_BASE - err_code);
}

void h_exit(int err_code, int line, const char *func, const char *file, ...) {
	fprintf(stderr, "ERROR: %s:%d: %s: ", file, line, func);
	va_list ap;
	va_start(ap, file);
	if (h_str_error[err_code] != NULL)
		vfprintf(stderr, h_str_error[err_code], ap);
	va_end(ap);
	fprintf(stderr, "\n");

	if (err_code == -1)
		exit(err_code);
	else
		exit(ERR_BASE - err_code);
}

void h_error(int line, const char *func, const char *file, char const *fmt, ...) {
	fprintf(stderr, "ERROR: %s:%d: %s: ", file, line, func);
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, "\n");
}

void h_error(int line, const char *func, const char *file, int err_code, ...) {
	fprintf(stderr, "ERROR: %s:%d: %s: ", file, line, func);
	va_list ap;
	va_start(ap, err_code);
	if (h_str_error[err_code] != NULL)
		vfprintf(stderr, h_str_error[err_code], ap);
	va_end(ap);
	fprintf(stderr, "\n");
}

void h_syserror(int line, const char *func, const char *file, char const *fmt, ...) {
	fprintf(stderr, "ERROR: %s:%d: %s: ", file, line, func);
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, " %s\n", strerror(errno));
}

void h_warning(const char *func, const char *fmt, ...) {
	fprintf(stderr, "WARNING: %s: ", func);
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, "\n");
}

void h_mem_check(int line, const char *func, const char *file, void *var, ...) {
	if (var == NULL) {
		fprintf(stderr, "ERROR: %s:%d: %s: ", file, line, func);
		fprintf(stderr, " %s\n", h_str_error[ERR_OUT_OF_MEMORY]);
	}
}
