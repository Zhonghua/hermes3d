//
// output.cc
//
//

#include "config.h"
#include "output.h"
#include "refdomain.h"
#include "quadstd.h"
#include "common.h"



//// Output ///////////////////////////////////////////////////////////////////////////////////////

Output::Output(OutputEngine *engine) {
	this->engine = engine;
}

Output::~Output() {
}

