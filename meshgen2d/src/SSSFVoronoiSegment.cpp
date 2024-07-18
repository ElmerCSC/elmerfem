#include "SSSFVoronoiSegment.h"
#include "Border.h"
#include "SSSFVertex.h"
#include "MGError.h"

#include <iostream>

#include <math.h>
#include <fstream>
#include <stdio.h>
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include "coreGeometry.h"

void SSSFVoronoiSegment::
discretize( NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements )
{
	triangulate();
	removeBoundaryConnectors();
	
	makeSeed(fixedNodes);

	generate();
	exportNodes( allNodes, allElements );
}

void SSSFVoronoiSegment::
makeSeed(NodeMap& fixedNodes)
{
	SSSFVertex *vtx;
	double nx, ny;
	
	if( explicitSeed )
	{
		int j;
		for (j = 0; j < 3; j++)
		{
			NodeMapIt it = fixedNodes.find(seedTags[j]);
			if (it == fixedNodes.end())
			{
				std::cerr << "Seed node " << seedTags[j] << "does not exist!" << std::endl;
				exit(1);
			}
			seed[j] = static_cast<MeshNode *>(it->second);
		}
		
		for( j = 0; j < 3; ++j )
		{		
			Vertex *v = root->firstAcceptableFound( seed[j]->x, seed[j]->y);
			
			if( addSite( v, seed[j], false ) == true )
			{
				int len = newVertices.size();
				for( int i = 0; i < len; ++i )
				{
					vtx = static_cast<SSSFVertex *>( newVertices[i] );
					vtx->setBody( tag );
				}
				recycle();
			}
			else
			{
				blm_error("Could not insert the seed node", seedTags[j]);
			}
		}
	}
	else
	{
		Node *a, *b;
		
		baseEdge->midNodes( a, b, baseDirection );
		
		seed[0] = static_cast<MeshNode *>(a);
		seed[1] = static_cast<MeshNode *>(b);
		
		vtx = static_cast<SSSFVertex *>( root->vertexWith( seed[0], seed[1] ) );
		vtx->computeNewCoordinates( *bg, nx, ny );
		seed[2] = new MeshNode( nx, ny );
		
		if( addSite( vtx, seed[2], false ) == true )
		{
			int len = newVertices.size();
			for( int i = 0; i < len; ++i )
			{
				vtx = static_cast<SSSFVertex *>( newVertices[i] );
				vtx->setBody( tag );
			}
			recycle();
		}
		else
		{
			std::cerr << "Could not insert the seed!" << std::endl;
			exit(1);
		}
	}	
	
	vtx = static_cast<SSSFVertex *>( root->vertexWith( seed[0], seed[1] ) );
	if (vtx == NULL)
	{
		std::cerr << "Could not insert the seed!" << std::endl;
		exit(1);
	}
	
	vtx->makeAccepted();
	vtx->radiate( actives, crystals );
}

void SSSFVoronoiSegment::
generate()
{
	SSSFVertex *vs;
	int i;
	double nx, ny;
	double ref;
	
	const double sqrt3 = 1.7320508075688772935274463415059;
	
	std::cout << "Generating" << std::endl;
	
	nx = ny = 0.0;

	int steps = 0;
	while( actives.size() > 0 || crystals.size() > 0)
	{
		bool active = true;
		SSSFVertex *vtx;
		
		if( crystals.size() > 0) 
		{
			active = false;
			vtx = static_cast<SSSFVertex *>( crystals.first() );
		}
		else
		{
			vtx = static_cast<SSSFVertex *>( actives.first() );
		}
		if( !vtx->borderTest( ) )
		{
			if( vtx->computeNewCoordinates( *bg, nx, ny ) > 0 )
			{
				MeshNode *nd = new MeshNode( nx, ny );
				
				if( addSite( vtx, nd, false, false ) == true )
				{
					if( deleted.size() == 1 )
					{
						// This is a terrible hack and we get a lot of it wrong!
						vtx->unDelete();
						vtx->makeAccepted();
						vtx->radiate( actives, crystals );
						deleted.erase( deleted.begin(), deleted.end() );
					}
					else 
					{
						std::list< Vertex* >::iterator vxIt;
						for( vxIt = deleted.begin(); vxIt != deleted.end(); ++vxIt )
						{
							vs = static_cast<SSSFVertex *>( *vxIt );
							if( vs->isCrystal() && vs->isAtHeap() )
								crystals.remove( vs );
							else if( vs->isActive() && vs->isAtHeap() )
								actives.remove( vs );
						}
						
						int len = newVertices.size();
						for( i = 0; i < len; ++i )
						{
							vs = static_cast<SSSFVertex *>( newVertices[i] );
							vs->setBody( tag );
						}
						
						for( i = 0; i < len; ++i )
						{
							vs = static_cast<SSSFVertex *>( newVertices[i] );
							ref = bg->interpolate( vs->vx, vs->vy );
							ref /= sqrt3;
							
							if( vs->rightSize( ref ) )
							{
								vs->makeAccepted();
							}
							else
							{
								vs->makeWaiting();
							}
						}
						
						for( i = 0; i < len; ++i )
						{
							vs = static_cast<SSSFVertex *>( newVertices[i] );
							if( vs->isWaiting() )
							{
								vs->testIfActive( actives, crystals );
							}
							else if( vs->isAccepted() )
							{
								vs->radiate( actives, crystals );
							}
						}
					}				
				}
				recycle();
			}
		}
		++steps;
	}
}

void SSSFVoronoiSegment::
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
		SSSFVertex *vtx = new SSSFVertex;
		v.push_back( vtx );
		allVertices.push_back( vtx );
	}	
}
