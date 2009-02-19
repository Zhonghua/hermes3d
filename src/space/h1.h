#ifndef _SPACE_H1_H_
#define _SPACE_H1_H_

#include "../space.h"

/// H1 space
///
/// @ingroup spaces
class H1Space : public Space {
public:
	H1Space(Mesh *mesh, Shapeset *ss);
	virtual ~H1Space();

	virtual Space *dup(Mesh *mesh) const;

   	virtual void get_element_assembly_list(Element *e, AsmList *al);
	virtual void get_boundary_assembly_list(Element *e, int face, AsmList *al);

protected:
	virtual int get_vertex_ndofs();
	virtual int get_edge_ndofs(order1_t order);
	virtual int get_face_ndofs(order2_t order);
	virtual int get_element_ndofs(order3_t order);

	virtual void assign_dofs_internal();

	virtual void calc_vertex_boundary_projection(Element *elem, int ivertex);
	virtual void calc_edge_boundary_projection(Element *elem, int iedge);
	virtual void calc_face_boundary_projection(Element *elem, int iface);

	virtual void update_constraints();
	void update_constrained_nodes(Word_t fid);
};

#endif
