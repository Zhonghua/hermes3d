#ifndef _TRAVERSE_H_
#define _TRAVERSE_H_

struct Mesh;
class  Transformable;
struct State;
struct Box;


struct UniData {
	Element *e;
	uint64 idx;
};

#define SPLIT_NONE							0x0000
#define SPLIT_HEX_X							0x0001
#define SPLIT_HEX_Y							0x0002
#define SPLIT_HEX_Z							0x0004
#define SPLIT_HEX_XY						SPLIT_HEX_X | SPLIT_HEX_Y
#define SPLIT_HEX_XZ						SPLIT_HEX_X | SPLIT_HEX_Z
#define SPLIT_HEX_YZ						SPLIT_HEX_Y | SPLIT_HEX_Z
#define SPLIT_HEX_XYZ						SPLIT_HEX_X | SPLIT_HEX_Y | SPLIT_HEX_Z


/// Traverse is a multi-mesh traversal utility class. Given N meshes sharing the
/// same base mesh it walks through all (pseudo-)elements of the union of all
/// the N meshes.
///
class Traverse {
public:
	void begin(int n, Mesh **meshes, Transformable **fn = NULL);
	void finish();

	Element **get_next_state(bool *bnd, FacePos *ep);
	Element  *get_base() const { return base; }

	UniData **construct_union_mesh(Mesh *unimesh);

private:
	int num;
	Mesh **meshes;
	Transformable **fn;

	State *stack;
	int top, size;

	int id;
	Element *base;
	int (*sons)[8];
	uint64 *subs;

	UniData **unidata;
	int udsize;

	State *push_state();
	void set_boundary_info(State *s, bool *bnd, FacePos *ep);
	void union_recurrent(Box *cr, Element **e, Box *er, uint64 *idx, Element *uni);
//	uint64 hex_init_idx(Box *cr, Box *er);

	void hex_union_rec(Box *cr, Element **e, Box *er, uint64 *idx, Element *uni);
	void hex_push_son_states(State *s);

	Mesh *unimesh;
};



#endif
