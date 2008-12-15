#ifndef _UTILS_H_
#define _UTILS_H_

//
// Miscelaneous utils
//

int maxn(int count, ...);

inline int max(int a, int b) {
	return a > b ? a : b;
}

void hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream);
void hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream);

#define countof(a) 								(sizeof(a)/sizeof(a[0]))

#endif
