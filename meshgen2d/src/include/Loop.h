#if !defined( MESH_LOOP_H )
#define MESH_LOOP_H

#include "GeometryEdge.h"
#include <set>
#include <vector>
#include <algorithm>

#include "../../config.h"
#if defined(_NO_STD_MINMAX) || defined(WIN32)
	#include "minmaxpatch.h"
#endif

class Loop
{
public:
	Loop() { }
	
	void addEdge( const int dir, GeometryEdge *ed )
	{
		direction.push_back( dir );
		edges.push_back( ed );
		if( dir == 1 )
			nodes.push_back( ed->base() );
		else
			nodes.push_back( ed->end() );
	}
	
	void collectNodes( std::vector< Node * >& chain )
	{
		int i;
		std::vector<Node*> strip;
		int len = edges.size();
		for( i = 0; i < len; ++i )
		{
			edges[i]->exportNodes(strip, direction[i]);
			strip.pop_back();
		}
		
		std::copy(strip.begin(), strip.end(), std::back_inserter( chain ));
	}
	
	void collectNodesAndPairs( std::vector< Node * >& chain, std::set< std::pair< int, int > >& links )
	{
		int i;
		std::vector<Node*> strip;
		int len = edges.size();
		for( i = 0; i < len; ++i )
		{
			edges[i]->exportNodes(strip, direction[i]);
			strip.pop_back();
		}
		
		len = strip.size();
		for( i = 0; i < len; ++i )
		{
			int t1 = strip[i]->tag;
			int t2 = strip[(i+1)%len]->tag;
			const std::pair<int, int> ip = std::make_pair(std::min(t1,t2),std::max(t1,t2));
			links.insert( ip );
		}
		
		std::copy(strip.begin(), strip.end(), std::back_inserter( chain ));
	}
	
	void collectGeometryNodes( std::vector< GeometryNode * >& chain )
	{
		int i;
		int len = edges.size();
		for( i = 0; i < len; ++i )
		{
			edges[i]->exportGeometryNodes( chain, direction[i] );
			chain.pop_back();
		}
	}
	
	void collectBoundaryElements( std::vector< BoundaryElement * >& chain )
	{
		int i;
		int len = edges.size();
		for( i = 0; i < len; ++i )
		{
			if( !edges[i]->isVirtual() )
				edges[i]->elements(chain, direction[i]);
		}
	}
	
	std::vector< int > direction;
	std::vector< GeometryNode * > nodes;
	std::vector< GeometryEdge * > edges;
};

#endif /* MESH_LOOP_H */
