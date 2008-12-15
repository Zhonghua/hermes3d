#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "solution.h"

/// @defgroup visualization Visualization
///
/// 

/// Abstract class for deriving classes which format the data for visualization 
///
/// @ingroup visualization
class OutputEngine {
public:
	/// Run the output with specified output engine
	///
	/// @return true if ok
	/// @param[in] fn A function that will be visualized
	virtual void out(MeshFunction *fn, const char *name, int item = FN_VAL_0) = 0;
	virtual void out(Mesh *mesh) = 0;
};


#endif
