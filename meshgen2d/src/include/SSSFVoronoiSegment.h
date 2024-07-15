#if !defined( MESH_SSSFVORONOISEGMENT_H )
#define MESH_SSSFVORONOISEGMENT_H 

#include "SSVoronoiSegment.h"
#include "PQ.h"

#include "BGMesh.h"
#include "SSSFVertex.h"

class SSSFVoronoiSegment : public SSVoronoiSegment
{
public:
	SSSFVoronoiSegment( const int t ) : SSVoronoiSegment(t) { }
	SSSFVoronoiSegment( const int t, BGMesh *bgMesh ) : 
	SSVoronoiSegment(t, bgMesh) { }
	
	void discretize( NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements );
	
protected:
	void makeSeed(NodeMap& fixedNodes);
	virtual void getVertices( std::vector< Vertex * >& v, const int count );
	void generate();
	
	pq crystals;
};
#endif /* MESH_SSSFVORONOISEGMENT_H  */
