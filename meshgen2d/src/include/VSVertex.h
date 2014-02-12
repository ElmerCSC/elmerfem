#if !defined( MESH_VSVERTEX_H )
#define MESH_VSVERTEX_H

#include "PQ.h"
#include "Vertex.h"
#include "BGMesh.h"

enum	vertexType { ACCEPTED=0, ACTIVE, WAITING, CRYSTAL };
class VSVertex : public Vertex
{
public:
	VSVertex() : Vertex() { type = WAITING; }
	VSVertex( Node* n1, Node* n2, Node* n3 ) : Vertex(n1, n2, n3) { type = WAITING; }
	void reset( Node* n1, Node* n2, Node* n3 )
	{
		Vertex::reset(n1,n2,n3);
		type = WAITING;
	}
	
	virtual double orderingValue() { return radius; }
	double ratioValue() { return ratio; }
	bool isWaiting() { return type == WAITING; }
	bool isAccepted() { return type == ACCEPTED; }
	bool isActive() { return type == ACTIVE; }
	
	void makeWaiting() { type = WAITING; }
	void makeAccepted() { type = ACCEPTED; }
	void makeActive() { type = ACTIVE; }
	
	bool testIfActive( pq& actives ) 
	{
		VSVertex *v1 = static_cast<VSVertex *>( vertices[0] );
		VSVertex *v2 = static_cast<VSVertex *>( vertices[1] );
		VSVertex *v3 = static_cast<VSVertex *>( vertices[2] );
		
		if( v1->isExternal() || v2->isExternal() ||
		    v3->isExternal() ||
		    v1->isAccepted() || v2->isAccepted() ||
		    v3->isAccepted() )
		{
			type = ACTIVE;
			actives.insert( this );
			return true;
		}
		return false;
	}
	
	void radiate( pq& actives ) 
	{
		VSVertex *v1 = static_cast<VSVertex *>( vertices[0] );
		VSVertex *v2 = static_cast<VSVertex *>( vertices[1] );
		VSVertex *v3 = static_cast<VSVertex *>( vertices[2] );
		
		if( !v1->isExternal() && v1->isWaiting() )
		{
			v1->makeActive();
			actives.insert( v1 );
		}
		if( !v2->isExternal() && v2->isWaiting() )
		{
			v2->makeActive();
			actives.insert( v2 );
		}
		if( !v3->isExternal() && v3->isWaiting() )
		{
			v3->makeActive();
			actives.insert( v3 );
		}
	}
	
	bool borderTest()
	{
		VSVertex *v1 = static_cast<VSVertex *>( vertices[0] );
		VSVertex *v2 = static_cast<VSVertex *>( vertices[1] );
		VSVertex *v3 = static_cast<VSVertex *>( vertices[2] );
		
		if( (v1->isAccepted() || v1->isExternal()) &&
		    (v2->isAccepted() || v2->isExternal()) &&
		    (v3->isAccepted() || v3->isExternal()))
		{
			makeAccepted();
			return true;
		}
		return false;
	}
	
	int computeNewCoordinates( BGMesh& bg, double& nx, double& ny );
	
protected:
	vertexType type;
};
#endif /* MESH_VSVERTEX_H */
