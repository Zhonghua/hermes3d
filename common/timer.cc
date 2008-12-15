//
// timer.cc
//

/*
#include <config.h>
#include <common/memcheck.h>
#ifdef DEBUG
#define new DEBUG_NEW
#define malloc DEBUG_MALLOC
#define free DEBUG_FREE
#endif
*/
#include "timer.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// Timer //////////////////////////////////////////////////////////////////////

Timer::Timer() {
	reset();
	this->name = NULL;
}

Timer::Timer(const char *name) {
	reset();
	this->name = new char [strlen(name) + 1];
	strcpy(this->name, name);
}

Timer::~Timer() {
	delete [] name;
}

void Timer::start() {
	gettimeofday(&this->st, NULL);
}

void Timer::stop() {
	struct timeval tv;
	gettimeofday(&tv, NULL);

	accum.tv_sec += tv.tv_sec - st.tv_sec;
	accum.tv_usec += tv.tv_usec - st.tv_usec;
}

void Timer::reset() {
	timerclear(&this->accum);
}

double Timer::get_seconds() {
	return accum.tv_sec + (accum.tv_usec / 1000000.0);
}

const char *Timer::get_human_time() {
	#define MAX_STR					200
	// FIXME: not thread safe!
	static char str[MAX_STR] = "";
	double secs = get_seconds();
	int hours = (int) secs / (3600);
	int mins = (int) fmod(secs, 3600) / 60;
	sprintf(str, "%dh %dm %02.2lfs", hours, mins, fmod(fmod(secs, 3600), 60));

	return str;
}

