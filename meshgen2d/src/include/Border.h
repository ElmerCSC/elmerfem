#if !defined( MESH_BORDER_H )
#define MESH_BORDER_H

#include <vector>
#include "Loop.h"

class Border
{
public:
	Border( const int loopCount ) 
	{ 
		loops.reserve( loopCount );
		for( int i = 0; i < loopCount; ++i )
		{
			loops.push_back( new Loop );
		} 
	}
	void addLoopEdge( const int loopTag, const int direction, GeometryEdge *ed )
	{
		loops[loopTag]->addEdge( direction, ed );
	}
	
	void collectNodes( std::vector< Node * >& chain)
	{
		int loopCount = loops.size();
		for( int i = 0; i < loopCount; ++i )
		{
			loops[i]->collectNodes( chain );
		} 
	}
	
	void collectNodesAndPairs( std::vector< Node * >& chain, std::set< std::pair< int, int > >& links )
	{
		int loopCount = loops.size();
		for( int i = 0; i < loopCount; ++i )
		{
			loops[i]->collectNodesAndPairs( chain, links );
		} 
	}
	
	void collectGeometryNodes( std::vector< GeometryNode * >& chain )
	{
		int loopCount = loops.size();
		for( int i = 0; i < loopCount; ++i )
		{
			loops[i]->collectGeometryNodes( chain );
		} 
	}
	
	void collectGeometryEdges( std::vector< GeometryEdge* >& eds, std::vector<int> &dirs )
	{
		for (int i = 0; i < loops.size(); i++)
		{
			eds.insert(eds.end(), loops[i]->edges.begin(), loops[i]->edges.end());
			dirs.insert(dirs.end(), loops[i]->direction.begin(), loops[i]->direction.end());
		}
	}

	void collectBoundaryElements( std::vector< BoundaryElement * >& chain )
	{
		int loopCount = loops.size();
		for( int i = 0; i < loopCount; ++i )
		{
			loops[i]->collectBoundaryElements( chain );
		} 
	}
	
	void copyLoop( std::vector< GeometryNode * >& nds, std::vector< GeometryEdge* >& eds, std::vector<int> &dirs)
	{
		nds = loops[0]->nodes;
		eds = loops[0]->edges;
		dirs = loops[0]->direction;
	}
	
	void copyLoops( std::vector< GeometryNode * >& nds, std::vector< GeometryEdge* >& eds, std::vector<int> &dirs)
	{
		for (int i = 0; i < loops.size(); i++)
		{
			nds.insert(nds.end(), loops[i]->nodes.begin(), loops[i]->nodes.end());
			eds.insert(eds.end(), loops[i]->edges.begin(), loops[i]->edges.end());
			dirs.insert(dirs.end(), loops[i]->direction.begin(), loops[i]->direction.end());
		}
	}
	
protected:
	std::vector< Loop * > loops;
};
#endif /* MESH_BORDER_H */
