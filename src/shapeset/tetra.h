#ifndef _SHAPESET_TETRA_H_
#define _SHAPESET_TETRA_H_

// tetra specific macros

// validation macros
#define CHECK_VERTEX(vertex)    	assert(vertex >= 0 && vertex < 4)

#define CHECK_EDGE(edge)	      	assert(edge >= 0 && edge < 6)
#define CHECK_EDGE_ORI(ori)	      	assert(ori == 0 || ori == 1)

#define CHECK_FACE(face)			assert(face >= 0 && face < 4)
#define CHECK_FACE_MODE(mode)		assert(mode == MODE_TRIANGLE)
#define CHECK_FACE_ORI(ori)	      	assert(ori >= 0 || ori <= 6)

//#define CHECK_PART(p)				assert(p >= 0)

#endif /* SHAPESET_TETRA_H_ */
