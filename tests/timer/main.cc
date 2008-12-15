/*
 * main.cc
 *
 * Test of timer
 *
 */

#include "config.h"
#include <hermes3d.h>

#include <common/trace.h>
#include <common/error.h>
#include <common/timer.h>
#include <unistd.h>

// helpers ////////////////////////////////////////////////////////////////////

bool test_print(bool value, const char *msg, bool correct) {
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

void print_result(bool value) {
	if (value) {
		printf("OK\n");
	}
	else {
		printf("failed\n");
	}
}

//
// tests themselves
//

bool test_timer_rough() {
	#define SECS		3

	// test #1 running timer fo 3 secs
	printf("* Running timer for %d secs...", SECS);
	fflush(stdout);

	Timer t("LowPrecision timer");
	t.start();
	int r = sleep(SECS);
	t.stop();

	bool result = (int) t.get_seconds() == SECS;
	print_result(result);
	if (!result) return false;

	return true;
}

bool test_timer_hiprecision() {
	// test #1 running timer fo 123 ms
	printf("* Running timer for 123 ms...");
	fflush(stdout);

	Timer t("HiPrecision timer");
	t.start();
	int r = usleep(123000);
	t.stop();

	// take only 2 decimal digits (the rest is beyond all recognition until we will measure time that CPU really spent on the
	// task)
	int v = (int) (t.get_seconds() * 100);
	bool result = v == 12;
	print_result(result);
	if (!result) return false;

	return true;
}

bool test_timer_output() {
	Timer t;
	t.start();
	sleep(2);
	t.stop();

	printf("* Time: %lf secs (%s)\n", t.get_seconds(), t.get_human_time());

	return true;
}


int main(int argc, char *argv[]) {
//	TRACE_START("trace.txt");
	DEBUG_OUTPUT_OFF;
	SET_VERBOSE_LEVEL(0);

	bool ret = ERR_SUCCESS;
	try {
		if (!test_timer_rough()) throw ERR_FAILURE;
		if (!test_timer_hiprecision()) throw ERR_FAILURE;
		if (!test_timer_output()) throw ERR_FAILURE;
	}
	catch (int e) {
		ret = false;
	}

	return ret;
}

