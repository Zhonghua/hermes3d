// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

#ifndef _ERROR_H_
#define _ERROR_H_

#include "compat-util.h"

//
// Error handling
//
// It is important to handle as much errors as possible (it will help to debug the code).
// When something bad happened, you will know were. Use at least EXIT and ERROR macros.
//

// error codes
#define ERR_FAILURE							-1
#define ERR_SUCCESS							0

// error handling functions

/// Report unrecoverable errors where you need to report the location of the error
/// It also reports call stack
#define EXIT(...) h_exit(__LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)
void h_exit(int line, const char *func, const char *file, char const *fmt, ...) NORETURN;

/// Report unrecoverable error (no call stack or location dumped)
void die(char const *fmt, ...) NORETURN;

/// Notify the user about an error (the execution continues), neither location or call stack
/// is dumped
void error(char const *fmt, ...);

/// Notify the user about warning (the execution continues), neither location or call stack
/// is dumped
void warning(const char *warn, ...);

/// Check that memory allocation was ok, it not, report an error (also dump call stack) and
/// terminate
#define MEM_CHECK(var) h_mem_check(__LINE__, __PRETTY_FUNCTION__, __FILE__, var)
void h_mem_check(int line, const char *func, const char *file, void *var);

#endif
