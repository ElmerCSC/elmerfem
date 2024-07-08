#include "VoronoiSegment.h"
#include "Border.h"
#include "VSVertex.h"

#include <iostream>
#include <math.h>
#include <fstream>
#include <stdio.h>
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

void VoronoiSegment::
discretize( NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements )
{
	triangulate();
	removeBoundaryConnectors();
	
	generate();
	exportNodes( allNodes, allElements );
}

void VoronoiSegment::
generate()
{
	const double sqrt3 = 1.7320508075688772935274463415059;
	
	std::cout << "Generating" << std::endl;
	
	std::list< Vertex* >::iterator vxIt;
	for( vxIt = allVertices.begin(); vxIt != allVertices.end(); ++vxIt )
	{
		VSVertex *vtx = static_cast<VSVertex *>(*vxIt);
		if( vtx->isDeleted() ) continue;
		
		if( vtx->isExternal() )
		{
			for( int i = 0; i < 3; ++i )
			{
				VSVertex *vs = static_cast<VSVertex *>( vtx->vertices[i] );
				if( vs == (VSVertex *) 0 ) continue;
				
				if( !vs->isExternal() && vs->isWaiting() )
				{
					vs->makeActive();
					actives.insert( vs );
				}
			}
		}
		else
		{
			double ref = bg->interpolate( vtx->vx, vtx->vy ) / sqrt3;
			if( vtx->rightSize( ref ) )
			{
				if( vtx->isActive() ) actives.remove( vtx );
				vtx->makeAccepted();
				for( int i = 0; i < 3; ++i )
				{
					VSVertex *vs = static_cast<VSVertex *>( vtx->vertices[i] );
					if( !vs->isExternal() && vs->isWaiting() )
					{
						vs->makeActive();
						actives.insert( vs );
					}
				}
			}
		}
	}

	while( actives.size() > 0 )
	{
		VSVertex *vtx = static_cast<VSVertex *>( actives.first() );
		
		if( !vtx->borderTest( ) )
		{
			double nx, ny;
			if( vtx->computeNewCoordinates( *bg, nx, ny ) > 0 )
			{
				MeshNode *nd = new MeshNode( nx, ny );
				
				if( addSite( vtx, nd, false, true ) == true )
				{
					if( deleted.size() == 1 )
					{
						// This is a terrible hack and we get a lot of it wrong!
						vtx->unDelete();
						vtx->makeAccepted();
						vtx->radiate( actives );
						deleted.erase( deleted.begin(), deleted.end() );
					}
					else 
					{
						actives.removeRelevants( deleted );
						
						int len = newVertices.size(), i;
						for( i = 0; i < len; ++i )
						{
							VSVertex *vs = static_cast<VSVertex *>( newVertices[i] );
							vs->setBody( tag );
						}
						
						for( i = 0; i < len; ++i )
						{
							VSVertex *vs = static_cast<VSVertex *>( newVertices[i] );
							double ref = bg->interpolate( vs->vx, vs->vy ) / sqrt3;
							if( vs->rightSize( ref ) )
							{
								vs->makeAccepted();
							}
							else
							{
								if( !vs->testIfActive( actives ) )
								{
									vs->makeWaiting();
								}
							}
						}
						
						for( i = 0; i < len; ++i )
						{
							VSVertex *vs = static_cast<VSVertex *>( newVertices[i] );
							if( vs->isWaiting() )
							{
								vs->testIfActive( actives );
							}
							else if( vs->isAccepted() )
							{
								vs->radiate( actives );
							}
						}
					}
				}
				// 		else std::cout << "Bad node" << std::endl;
				
				recycle();
			}
		}
	}
}

void VoronoiSegment::
getVertices( std::vector< Vertex * >& v, const int count )
{
	int i;
	for( i = 0; i < count; ++i )
	{
		if( vertexStore.empty() ) break;
		v.push_back( vertexStore.back() );
		vertexStore.pop_back();
	}
	for( ; i < count; ++i )
	{
		VSVertex *vtx = new VSVertex;
		v.push_back( vtx );
		allVertices.push_back( vtx );
	}	
}
