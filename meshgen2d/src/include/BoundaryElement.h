#if !defined( MESH_BOUNDARYELEMENT_H )
#define MESH_BOUNDARYELEMENT_H

#include "Node.h"
#include "Vertex.h"
#include <iostream>
#include <fstream>

class BoundaryElement
{
public:
	friend std::ostream& operator<< (std::ostream& o, const BoundaryElement& A);
	BoundaryElement( int edgeTag, Node* n1, Node* n2 )
	{
		edge = edgeTag;
		a = n1;
		b = n2;
		c = NULL;
		left = right = 0;
		flipped = false;
		newTag();
	}
	
	void newTag();
	
	void flip( bool status ) { flipped = status; }
	
	Vertex* setHolder( Vertex *v ) 
	{
		Vertex *holder;
		
		if( !flipped ) 
		{
			holder = v->vertexWith(a, b);
			left = holder->elementId();
		}
		else 
		{
			holder = v->vertexWith(b, a);
			right = holder->elementId();
		}
		
		return holder;
	}
	
	void setLeft( int id ) { (flipped?right:left) = id; }
	void setRight( int id ) { (flipped?left:right) = id; }

	Node *from() { return a; }
	Node *to() { return b; }
	Node *middle() { return c; }

	void addMiddleNode(Node *n) { c = n; }
	
private:
	int edge;
	Node *a, *b, *c;
	int left, right;
	bool flipped;
	int tag;
};

std::ostream& operator<< (std::ostream& o, const BoundaryElement& A);

#endif /* MESH_BOUNDARYELEMENT_H */
