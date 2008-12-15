#ifndef _ERROR_H_
#define _ERROR_H_

//
// Error handling
//
// It is important to handle as much errors as possible (it will help to debug the code).
// When something bad happened, you will know were. Use at least EXIT and ERROR macros.
//
// In order to add a new error:
// 1. Assign a new identifier (ERR_xxx) and add it to the section error codes
// 2. If the error has static description
//    - add the description to the h_str_error array in error.cc
//    If the error state requires a description
//    - add the NULL to the h_str_error array in error.cc
//
//    The index of the description in h_str_error array (in error.cc) must match the
//    numerical value of error ID.
//

// error codes
#define ERR_FAILURE							-1

#define ERR_SUCCESS							0

#define ERR_BASE							-1000
#define ERR_OUT_OF_MEMORY					0
#define ERR_NOT_ENOUGH_PARAMS				1
#define ERR_CAN_NOT_OPEN_FILE				2
#define ERR_MESH_PROBLEM                    3
#define ERR_NOT_IMPLEMENTED                 4
#define ERR_UNKNOWN_MODE                    5
#define ERR_PETSC_NOT_COMPILED              6
#define ERR_UMFPACK_NOT_COMPILED            7
#define ERR_HDF5_NOT_COMPILED               8
#define ERR_MPI_NOT_COMPILED                9
#define ERR_FACE_INDEX_OUT_OF_RANGE         10
#define ERR_EDGE_INDEX_OUT_OF_RANGE         11
#define ERR_TETRA_NOT_COMPILED              12
#define ERR_HEX_NOT_COMPILED                13
#define ERR_PRISM_NOT_COMPILED              14
#define ERR_PARDISO_NOT_COMPILED            15
#define ERR_UNKNOWN_REFINEMENT_TYPE         16

// error handling functions

// use EXIT for unrecoverable error occurs (out of memory conditions)
#define EXIT(errcode, ...) h_exit(errcode, __LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)
void h_exit(int err_code, int line, const char *func, const char *file, char const *fmt, ...);
void h_exit(int err_code, int line, const char *func, const char *file, ...);

// notify the user about errors by using ERROR macro
#define ERROR(...) h_error(__LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)
void h_error(int line, const char *func, const char *file, char const *fmt, ...);
void h_error(int line, const char *func, const char *file, int err_code, ...);

// use SYSERROR for system errors like error opening/closing file
#define SYSERROR(...) h_error(__LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)
void h_syserror(int line, const char *func, const char *file, char const *fmt, ...);

// notify the user about warning (execution continues)
#define WARNING(...) h_warning(__PRETTY_FUNCTION__,  ## __VA_ARGS__)
void h_warning(const char *func, const char *fmt, ...);

#define MEM_CHECK(var, ...) h_mem_check(__LINE__, __PRETTY_FUNCTION__, __FILE__, var, ## __VA_ARGS__)
void h_mem_check(int line, const char *func, const char *file, void *var, ...);

#endif
