#include <stdarg.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "error.h"

/// Find the largest value
///
/// @return The largest value among numbers passed as arguments
/// @param[in] count - the number of values passed into this function (all of them are searched)
///
int maxn(int count, ...) {
	va_list ap;
	va_start(ap, count);
	int mx = INT_MIN;

	for (int i = 0; i < count; i++) {
		int num = va_arg(ap, int);
		if (num > mx)
			mx = num;
	}

	return mx;
}

void hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream) {
	if (fwrite(ptr, size, nitems, stream) != nitems || ferror(stream))
		EXIT(ERR_FAILURE, "Error writing to file: %s", strerror(ferror(stream)));
}

void hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream) {
	if (fread(ptr, size, nitems, stream) != nitems || ferror(stream))
		EXIT(ERR_FAILURE, "Error reading file: %s", strerror(ferror(stream)));
}
