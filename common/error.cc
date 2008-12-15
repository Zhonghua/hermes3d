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
