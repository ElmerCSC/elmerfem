#include "VoronoiVertex.h"
#include "Border.h"
#include <math.h>
#include <iostream>

void VoronoiVertex::
discretize(NodeMap& fixedNodes,	NodeMap& allNodes, std::list< Element* >& allElements )
{
	triangulate();
	
	removeBoundaryConnectors();
	generate();
	exportNodes( allNodes, allElements );
}

void VoronoiVertex::
generate()
{
	std::cout << "Generating" << std::endl;
	
	const double sqrt3 = 1.7320508075688772935274463415059;
	
	double ref;
	std::list< Vertex* >::iterator vxIt;
	for( vxIt = allVertices.begin(); vxIt != allVertices.end(); ++vxIt )
	{
		Vertex *vtx = (*vxIt);
		if( !(vtx->isDeleted() || vtx->isExternal() ))
		{
			ref = bg->interpolate( vtx->vx, vtx->vy );
			if( !vtx->rightSize( ref / sqrt3 ) )
			{
				actives.insert( vtx );
			}
		}
	}
	
	size_t steps = 1;
	
	while( actives.size() > 0 )
	{
		Vertex *vtx = actives.first();

		MeshNode *nd = new MeshNode( vtx->vx, vtx->vy );
		
		int i = 0;
		
		if( addSite( vtx, nd, false, false ) == true )
		{
			if( deleted.size() == 1 )
			{
				// This is a terrible hack and we get a lot of it wrong!
				vtx->unDelete();
				deleted.clear();
				delete nd;
			}
			else
			{
				actives.removeRelevants( deleted );
				
				int len = newVertices.size();
				for( int i = 0; i < len; ++i )
				{
					newVertices[i]->setBody( 1 );
					ref = bg->interpolate( newVertices[i]->vx, newVertices[i]->vy );
					if( !newVertices[i]->rightSize( ref / sqrt3 ) )
					{
						actives.insert( newVertices[i] );
					}
				}
			}
		}
		
		recycle();
	}
}

