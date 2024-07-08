#if !defined( MESH_VORONOIVERTEX_H )
#define MESH_VORONOIVERTEX_H 

#include "Connect.h"
#include "PQ.h"

#include "BGMesh.h"

class VoronoiVertex : public Connect
{
public:
	VoronoiVertex( const int t ) : Connect(t) { }
	VoronoiVertex( const int t, BGMesh *bgMesh ) : 
	Connect(t, bgMesh) { }
	
	void discretize(NodeMap& fixedNodes,  NodeMap& allNodes, std::list< Element* >& allElements );
	
protected:
	void generate();
	
	pq actives;
};
#endif /* MESH_VORONOIVERTEX_H  */
