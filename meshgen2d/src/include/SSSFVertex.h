#if !defined( MESH_SSSFVERTEX_H )
#define MESH_SSSFVERTEX_H

#include "VSVertex.h"

class SSSFVertex : public VSVertex
{
public:
	SSSFVertex() : VSVertex() { }
	SSSFVertex( Node* n1, Node* n2, Node* n3 ) : VSVertex(n1, n2, n3) { }
	
	bool testIfActive( pq& actives, pq& crystals ) 
	{
		SSSFVertex *v1 = static_cast<SSSFVertex *>( vertices[0] );
		SSSFVertex *v2 = static_cast<SSSFVertex *>( vertices[1] );
		SSSFVertex *v3 = static_cast<SSSFVertex *>( vertices[2] );
		
		if( v1->isAccepted() || v2->isAccepted() || v3->isAccepted() )
		{
			makeActive();
			if( isActive() )
				actives.insert( this );
			else if( isCrystal() )
				crystals.insert( this );
			return true;
		}
		return false;
	}
	
	void radiate( pq& actives, pq& crystals ) 
	{
		SSSFVertex *v1 = static_cast<SSSFVertex *>( vertices[0] );
		SSSFVertex *v2 = static_cast<SSSFVertex *>( vertices[1] );
		SSSFVertex *v3 = static_cast<SSSFVertex *>( vertices[2] );
		
		if( !v1->isExternal() && v1->isWaiting() )
		{
			v1->makeActive();
			if( v1->isActive() )
				actives.insert( v1 );
			else if( v1->isCrystal() )
				crystals.insert( v1 );
		}
		if( !v2->isExternal() && v2->isWaiting() )
		{
			v2->makeActive();
			if( v2->isActive() )
				actives.insert( v2 );
			else if( v2->isCrystal() )
				crystals.insert( v2 );
		}
		if( !v3->isExternal() && v3->isWaiting() )
		{
			v3->makeActive();
			if( v3->isActive() )
				actives.insert( v3 );
			else if( v3->isCrystal() )
				crystals.insert( v3 );
		}
	}
	
	void makeAccepted()
	{
		MeshNode *n1 = static_cast<MeshNode *>( nodes[0] );
		MeshNode *n2 = static_cast<MeshNode *>( nodes[1] );
		MeshNode *n3 = static_cast<MeshNode *>( nodes[2] );
		
		type = ACCEPTED;
		n1->putOnCrysralIfNotFixed();
		n2->putOnCrysralIfNotFixed();
		n3->putOnCrysralIfNotFixed();
	}
	
	void makeActive()
	{
		MeshNode *n1 = static_cast<MeshNode *>( nodes[0] );
		MeshNode *n2 = static_cast<MeshNode *>( nodes[1] );
		MeshNode *n3 = static_cast<MeshNode *>( nodes[2] );
		
		if( n1->isCrystal() && n2->isCrystal() && n3->isCrystal() )
		{
			// 		type = CRYSTAL;
			// 	return;
			
			SSSFVertex *v1 = static_cast<SSSFVertex *>( vertices[0] );
			SSSFVertex *v2 = static_cast<SSSFVertex *>( vertices[1] );
			SSSFVertex *v3 = static_cast<SSSFVertex *>( vertices[2] );
			
			int counter = 0;
			if(v1 != (SSSFVertex *)0) if(v1->isAccepted()) ++counter;
			if(v2 != (SSSFVertex *)0) if(v2->isAccepted()) ++counter;
			if(v3 != (SSSFVertex *)0) if(v3->isAccepted()) ++counter;
			if(counter >= 2)
			{
				type = CRYSTAL;
				return;
			}
		}
		type = ACTIVE;
	}
	
	bool isCrystal() { return type == CRYSTAL; }
	bool wasCrystal()
	{
		MeshNode *n1 = static_cast<MeshNode *>( nodes[0] );
		MeshNode *n2 = static_cast<MeshNode *>( nodes[1] );
		MeshNode *n3 = static_cast<MeshNode *>( nodes[2] );
		
		if( n1->isCrystal() && n2->isCrystal() && n3->isCrystal() )
			return true;
		return false;
	}
};

#endif /* MESH_SSSFVERTEX_H */


