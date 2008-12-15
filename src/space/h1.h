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
	virtual int get_vertex_ndofs(Order0 order);
	virtual int get_edge_ndofs(Order1 order);
	virtual int get_face_ndofs(Facet *face, Order2 order);
	virtual int get_element_ndofs(Element *elem, Order3 order);

	virtual int assign_dofs_internal(int first_dof = 0, int stride = 1);

	virtual void calc_vertex_boundary_projection(Element *elem, int ivertex);
	virtual void calc_edge_boundary_projection(Element *elem, int iedge);
	virtual void calc_face_boundary_projection(Element *elem, int iface);

	virtual void update_constraints();
	void update_constrained_nodes(Word_t fid);
};

#endif
