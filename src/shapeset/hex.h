#ifndef _SHAPESET_HEX_H_
//#define _SHAPESET_HEX_H_

// hex specific macros

// validation macros
#define CHECK_VERTEX(vertex)    	assert(vertex >= 0 && vertex < 8)

#define CHECK_EDGE(edge)	      	assert(edge >= 0 && edge < 12)
#define CHECK_EDGE_ORI(ori)	      	assert(ori == 0 || ori == 1)
#define CHECK_EDGE_ORDER(o)  		assert((o) >= 0 && (o) <= max_edge_order)

#define CHECK_FACE(face)			assert(face >= 0 && face < 6)
#define CHECK_FACE_MODE(mode)		assert(mode == MODE_QUAD)
#define CHECK_FACE_ORI(ori)	      	assert(ori >= 0 || ori <= 8)
#define CHECK_FACE_ORDER(o)  		assert((o) >= 0 && (o) <= max_face_order)

#define CHECK_PART(p)				assert(p >= 0)


#endif /* _SHAPESET_HEX_H_ */
