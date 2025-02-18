#include "BGMesh.h"
#include "BGVertex.h"
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <stdio.h>

void BGTriangleMesh::
initialize(std::vector< GeometryNode* >& nds, std::vector< GeometryEdge* >& eds)
{
	initialized = true;
	
	nodes = nds;
	
	int i;
	int len = nodes.size();
	
	for( i = 0; i < len; ++i )
	{
		border.push_back( nodes[i] );
	}
	
	makeWorld();
	
	len = border.size();
	int *indirect = new int[len];
	for( i = 0; i < len; ++i)
	{
		indirect[i] = i;
	}
	std::random_shuffle(indirect, indirect + len);
	
	for( i = 0; i < len; ++i)
	{
		int ind = indirect[i];
		Vertex *vtx = root->firstAcceptableFound( border[ind]->x, border[ind]->y);
		addSite( vtx, border[ind], true );
		GeometryNode *nd = static_cast<GeometryNode *>( border[ind] );
		recycle();
	}
	
	delete [] indirect;
	
	// We have to set deltas for the BGVertices and compute the interpolation parameters
	std::list< Vertex* >::iterator vxIt;
	for( vxIt = allVertices.begin(); vxIt != allVertices.end(); ++vxIt )
	{
		BGVertex *vtx = static_cast< BGVertex * >(*vxIt);
		if( !(vtx->isDeleted()) )
		{
			vtx->initInterpolation( );
		}
	}
	dump();		
}

void BGTriangleMesh::
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
		BGVertex *vtx = new BGVertex;
		v.push_back( vtx );
		allVertices.push_back( vtx );
	}	
}

double BGTriangleMesh::
interpolate( const double ix, const double iy )
{
	BGVertex *vtx = static_cast< BGVertex * >( root );
	return vtx->interpolate( ix, iy );
}

void blm_error(const char* msg, const char *opt = "");

void BGTriangleMesh::
dump()
{
}
