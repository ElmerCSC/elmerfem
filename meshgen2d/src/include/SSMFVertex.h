#if !defined( MESH_SSMFVERTEX_H )
#define MESH_SSMFVERTEX_H

#include "VSVertex.h"

class SSMFVertex : public VSVertex
{
public:
	SSMFVertex() : VSVertex() { }
	SSMFVertex( Node* n1, Node* n2, Node* n3 ) : VSVertex(n1, n2, n3) { }
	
	bool testIfActive( pq& actives ) 
	{
		VSVertex *v1 = static_cast<VSVertex *>( vertices[0] );
		VSVertex *v2 = static_cast<VSVertex *>( vertices[1] );
		VSVertex *v3 = static_cast<VSVertex *>( vertices[2] );
		
		if( v1->isAccepted() || v2->isAccepted() ||	v3->isAccepted() )
		{
			type = ACTIVE;
			actives.insert( this );
			return true;
		}
		return false;
	}
};

#endif /* MESH_SSMFVERTEX_H */
