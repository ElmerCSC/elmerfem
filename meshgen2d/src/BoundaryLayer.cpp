#include "BoundaryLayer.h"
#include <algorithm>
#include "Border.h"
#include "BGMesh.h"
#include "MGError.h"

#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#ifdef _NO_STD_MINMAX
	#include "minmaxpatch.h"
#endif

void BoundaryLayer::initialize()
{
	bounds->collectGeometryNodes( nodes );
	bounds->collectGeometryEdges( edges, directions );
	if(!bg->isInitialized()) bg->initialize( nodes, edges );
}

/*
  Basically a do nothing. The boundaries have been already discretized.
*/
void BoundaryLayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	exportNodes( allNodes, allElements );
}

void BoundaryLayer::
exportNodes( NodeMap& allNodes, std::list< Element* >& allElements )
{	
	int i, len;
	std::vector< Node * > nodes;
	bounds->collectNodes( nodes );

	len = nodes.size();

	for ( i = 0; i < len; i++)
	{
		allNodes[nodes[i]->tag] = nodes[i];
	}

	std::vector< BoundaryElement * > bels;
	bounds->collectBoundaryElements( bels );
	
	len = bels.size();
	
	for( i = 0; i < len; ++i )
	{
		bels[i]->setRight(-1);
	}
}
