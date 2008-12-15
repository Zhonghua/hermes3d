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
	void start();
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
