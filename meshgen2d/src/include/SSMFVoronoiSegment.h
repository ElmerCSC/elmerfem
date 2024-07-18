#if !defined( MESH_SSMFVORONOISEGMENT_H )
#define MESH_SSMFVORONOISEGMENT_H 

#include "SSVoronoiSegment.h"
#include "PQ.h"

#include "BGMesh.h"
#include "SSMFVertex.h"

class SSMFVoronoiSegment : public SSVoronoiSegment
{
public:
	SSMFVoronoiSegment( const int t ) : SSVoronoiSegment(t) { }
	SSMFVoronoiSegment( const int t, BGMesh *bgMesh ) : 
	SSVoronoiSegment(t, bgMesh) { }
	
	void discretize( NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements );
	
protected:
	void makeSeed(NodeMap& fixedNodes);
	virtual void getVertices( std::vector< Vertex * >& v, const int count );
	void generate();
};
#endif /* MESH_SSMFVORONOISEGMENT_H  */
