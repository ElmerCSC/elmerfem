#if !defined( MESH_VORONOISEGMENT_H )
#define MESH_VORONOISEGMENT_H 

#include "Connect.h"
#include "PQ.h"

#include "BGMesh.h"
#include "VSVertex.h"

class VoronoiSegment : public Connect
{
public:
	VoronoiSegment( const int t ) : Connect(t) { }
	VoronoiSegment( const int t, BGMesh *bgMesh ) : 
	Connect(t, bgMesh) { }
	
	void discretize( NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements );
	
protected:
	virtual void getVertices( std::vector< Vertex * >& v, const int count );
	void generate();
	
	pq actives;
};
#endif /* MESH_VORONOISEGMENT_H	*/
