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

#ifndef _H_
#define _H_

#include <sys/time.h>
#include <sys/resource.h>

/// \class Timer
///
/// TODO: Measure time that CPU spent on the task
///
class Timer {
public:
	Timer();						// default constructor
	Timer(const char *name);
	virtual ~Timer();				// destructor

	/// start the timer
	/// @param[in] rst - Reset timer of true
	void start(bool rst = true);
	/// stop the timer
	void stop();
	/// reset the timer
	void reset();

	const char *get_name() { return name; }
	double get_seconds();
	const char *get_human_time();

protected:
	/// name of the timer (can be NULL)
	char *name;
	/// time when the timer was started/resumed
	struct timeval st;
	/// accumulator
	struct timeval accum;
};


#endif
