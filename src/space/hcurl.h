#ifndef _SPACE_HCURL_H_
#define _SPACE_HCURL_H_

#include "../space.h"

///
/// \class HCurlSpace
/// \brief
///
class HCurlSpace : public Space {
public:
	HCurlSpace(Mesh *mesh, Shapeset *ss);
	virtual ~HCurlSpace();

	virtual Space *dup(Mesh *mesh) const;

	virtual void get_element_assembly_list(Element *e, AsmList *al);
	virtual void get_boundary_assembly_list(Element *e, int face, AsmList *al);

protected:
	virtual int get_vertex_ndofs(int order);
	virtual int get_edge_ndofs(int order);
	virtual int get_face_ndofs(Facet *face, order2_t order);
	virtual int get_element_ndofs(Element *elem, order3_t order);

	virtual int assign_dofs_internal(int first_dof = 0, int strid = 1);

	// for now we do not implement boundary projections of nonzero functions
	// it does not have much physical sense (even though some artifical situations
	// using nonzero dirichlet bc could be considered) and its implementation is
	// not as straightforward as in 2D. It is given by the fact, that there are
	// two tangential components and simple condition E * t = g does not suffice.
	// So the only essential bc, that we consider so far is "perfect conductor"
	// condition E * t = 0.
	virtual void calc_vertex_boundary_projection(Element *elem, int ivertex);
	virtual void calc_edge_boundary_projection(Element *elem, int iedge);
	virtual void calc_face_boundary_projection(Element *elem, int iface);

	virtual void update_constraints();
};

#endif
