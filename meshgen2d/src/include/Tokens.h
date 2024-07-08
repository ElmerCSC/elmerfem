#if !defined( MESH_MESHTOKENS_H )
#define MESH_MESHTOKENS_H

#include <list>
#include <vector>
#include <string>

class NodeToken
{
public:
	NodeToken() { }
	int tag;
	int boundaryTag;
	double x,y,delta;
};


class EdgeToken
{
public:
	EdgeToken() { segCount = 0; }
	
	int from() const { return nodes.front(); }
	int to() const { return nodes.back(); }
	
	int tag;
	int boundaryTag;
	
	std::vector<int> nodes;
	
	int segCount;
	double delta;
};

class LoopToken
{
public:
	LoopToken() { }
	int tag;
	int direction;
	std::vector<int> edges;
};

class BGMeshToken
{
public:
	BGMeshToken() { }
	std::string type;
	std::vector< NodeToken > nodes;
};

class SeedToken
{
public:
	SeedToken() { }
	std::string type;
	// union{ 
	std::vector<int> nodes;
	int edge;	
	// };
};

enum { CONNECT, BOUNDARY_MESH, VORONOI_VERTEX, VORONOI_SEGMENT, SSSF_VORONOI_SEGMENT, SSMF_VORONOI_SEGMENT,
QUAD_GRID, TRIANGLE_NE_GRID, TRIANGLE_NW_GRID, TRIANGLE_UJ_NE_GRID, TRIANGLE_UJ_NW_GRID,
TRIANGLE_FB_NE_GRID, TRIANGLE_FB_NW_GRID };

class LayerToken
{
public:
	LayerToken() { bg = NULL; seed = NULL; gridh = 0; gridv = 0; }
	int tag;
	int type;
	SeedToken *seed;
	std::vector< LoopToken * > loops;
	std::vector< NodeToken * > fixedNodes;
	BGMeshToken *bg;
	double delta;
	int gridh, gridv;
};

class BodyToken
{
public:
	BodyToken() { }
	int tag;
	std::vector< LayerToken * > layers;
	double delta;
	bool parabolic;
};

class BoundaryToken
{
public:
	BoundaryToken() { }
	std::vector<int> outerBoundaries;
	std::vector<int> innerBoundaries;
};

#endif /* MESH_MESHTOKENS_H */
