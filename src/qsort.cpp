// The QuickSort routine from glibc-2.5 modified for sorting int arrays.

#include <limits.h>
#include <stdlib.h>

#define SWAP(a, b) { int c = *(a); *(a) = *(b); *(b) = c; }

#define H3D_MAX_THRESH 4

typedef struct {
	int *lo;
	int *hi;
} stack_node;

#define H3D_STACK_SIZE      (CHAR_BIT * sizeof(size_t))
#define PUSH(low, high) ((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define POP(low, high)  ((void) (--top, (low = top->lo), (high = top->hi)))
#define H3D_STACK_NOT_EMPTY (stack < top)

#define min(x, y) ((x) < (y) ? (x) : (y))

void qsort_int(int *pbase, size_t total_elems) {
	register int *base_ptr = pbase;

	if (total_elems == 0)
		return;

	if (total_elems > H3D_MAX_THRESH) {
		int *lo = base_ptr;
		int *hi = lo + total_elems - 1;
		stack_node stack[H3D_STACK_SIZE];
		stack_node *top = stack;

		PUSH(NULL, NULL);

		while (H3D_STACK_NOT_EMPTY) {
			int *left_ptr;
			int *right_ptr;

			// Select median value from among LO, MID, and HI. Rearrange
			// LO and HI so the three values are sorted. This lowers the
			// probability of picking a pathological pivot value and
			// skips a comparison for both the LEFT_PTR and RIGHT_PTR in
			// the while loops.

			int *mid = lo + ((hi - lo) >> 1);

			if (*mid < *lo) SWAP (mid, lo);
			if (*hi < *mid) SWAP (mid, hi) else
				goto jump_over;
			if (*mid < *lo) SWAP(mid, lo);
			jump_over:

			left_ptr = lo + 1;
			right_ptr = hi - 1;

			// Here's the famous ``collapse the walls'' section of quicksort.
			// Gotta like those tight inner loops!  They are the main reason
			// that this algorithm runs much faster than others.
			do {
				while (*left_ptr < *mid)
					left_ptr++;

				while (*mid < *right_ptr)
					right_ptr--;

				if (left_ptr < right_ptr) {
					SWAP(left_ptr, right_ptr);
					if (mid == left_ptr)
						mid = right_ptr;
					else if (mid == right_ptr)
						mid = left_ptr;
					left_ptr++;
					right_ptr--;
				}
				else if (left_ptr == right_ptr) {
					left_ptr++;
					right_ptr--;
					break;
				}
			} while (left_ptr <= right_ptr);

			// Set up pointers for next iteration.  First determine whether
			// left and right partitions are below the threshold size.  If so,
			// ignore one or both.  Otherwise, push the larger partition's
			// bounds on the stack and continue sorting the smaller one.

			if ((size_t) (right_ptr - lo) <= H3D_MAX_THRESH) {
				if ((size_t) (hi - left_ptr) <= H3D_MAX_THRESH)
					// Ignore both small partitions.
					POP(lo, hi);
				else
					// Ignore small left partition.
					lo = left_ptr;
			} else if ((size_t) (hi - left_ptr) <= H3D_MAX_THRESH)
				// Ignore small right partition.
				hi = right_ptr;
			else if ((right_ptr - lo) > (hi - left_ptr)) {
				// Push larger left partition indices.
				PUSH(lo, right_ptr);
				lo = left_ptr;
			} else {
				// Push larger right partition indices.
				PUSH(left_ptr, hi);
				hi = right_ptr;
			}
		}
	}

	// Once the BASE_PTR array is partially sorted by quicksort the rest
	// is completely sorted using insertion sort, since this is efficient
	// for partitions below H3D_MAX_THRESH size. BASE_PTR points to the beginning
	// of the array to sort, and END_PTR points at the very last element in
	// the array (*not* one beyond it!).

	{
		int *const end_ptr = base_ptr + total_elems - 1;
		int *tmp_ptr = base_ptr;
		int *thresh= min(end_ptr, base_ptr + H3D_MAX_THRESH);
		register int *run_ptr;

		// Find smallest element in first threshold and place it at the
		// array's beginning.  This is the smallest array element,
		// and the operation speeds up insertion sort's inner loop.

		for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
			if (*run_ptr < *tmp_ptr)
				tmp_ptr = run_ptr;

		if (tmp_ptr != base_ptr)SWAP(tmp_ptr, base_ptr);

		// Insertion sort, running from left-hand-side up to right-hand-side.
		run_ptr = base_ptr + 1;
		while (++run_ptr <= end_ptr) {
			tmp_ptr = run_ptr - 1;
			while (*run_ptr < *tmp_ptr)
				tmp_ptr--;

			tmp_ptr++;
			if (tmp_ptr != run_ptr) {
				int c = *run_ptr;
				register int *hi, *lo;

				for (hi = lo = run_ptr; --lo >= tmp_ptr; hi = lo)
					*hi = *lo;
				*hi = c;
			}
		}
	}
}

/* test code

 const int max = 10000;
 int array[max];
 int i, j, n;
 for (i = 0; i < 1000; i++)
 {
 n = rand() % max;
 for (j = 0; j < n; j++)
 array[j] = j;
 for (j = 0; j < n; j++)
 swap(array[j], array[rand() % n]);
 qsort_int(array, n);
 for (j = 0; j < n; j++)
 if (array[j] != j)
 warning("failed.");
 }
 info("passed");

 */
