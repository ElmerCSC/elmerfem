#if !defined( MESH_VERTEX_H )
#define MESH_VERTEX_H

#include "TriangleElement.h"

#include <set>

class Node;

class Vertex : public TriangleElement
{
public:
	Vertex();
	Vertex( Node* n1, Node* n2, Node* n3 );
	void reset( Node* n1, Node* n2, Node* n3 );
	
	void setVertices(Vertex *v1, Vertex *v2, Vertex *v3);
	void replaceVertexWith(const int ind, Vertex* vtx) { vertices[ind] = vtx; }
	
	Vertex* firstAcceptableFound(double tx, double ty);
	Vertex* vertexWith( Node* a, Node* b);
	Vertex* vertexOwning( Node* a, Node* b);
	int isDestroyedBy(double tx, double ty, bool extOK);
	
	int id() const { return tag; }
	
	bool isDeleted() { return deleted; }
	void makeDeleted() { deleted = true; }
	void unDelete() { deleted = false; }
	
	void labelBodies(  const int bd, std::set< std::pair< int, int > >& bounds );
	
	bool isExternal() { return body == 0; }
	
	bool rightSize( double ref ) 
	{
		const double DELTA = 1.25;
		
		ratio = radius / ref;
		if( ratio < DELTA ) return true;
		return false;
	}
	
	virtual double orderingValue() { return ratio; }
	
	bool isAtHeap() { return heap != -1; }
	void setHeap( const int p ) { heap = p; }
	int atHeap() { return heap; }
	
protected:
	
	int tag;
	
	double radius;
	double ratio;
	
	bool deleted;
	
	int heap;
public:	
	
	Vertex *vertices[3];
	double vx, vy;
};

#endif /* MESH_VERTEX_H */
