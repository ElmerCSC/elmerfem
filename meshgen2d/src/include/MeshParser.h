#if !defined( MESH_MESHPARSER_H )
#define MESH_MESHPARSER_H 

#include <fstream>
#include <vector>
#include <map>

#include "Tokens.h"

class MeshParser
{
public:
	friend class Mesh;
	MeshParser( const char *path );
	MeshParser( const char *path, std::map<int, BGMeshToken*> *ebgs);
	~MeshParser();

	void readInputFile();
	int readExplicitBGMesh(BGMeshToken *bg);
	
protected:
	int ignoreComment();
	int readHeader();
	int readNodes();
	int readEdges();
	int readBodies();
	int readLayers( const int layerCount );
	int readSeed();
	int readLoops( const int loopCount );
	int readFixedNodes( const int fixedNodeCount );
	int readBGMesh();
	int readBoundary();
	
	int readWord( const char *ref );
	int readInteger( int& ref );
	int readDouble( double& ref );
	int readString( std::string& ref );
	
	std::ifstream s;
	
	int nodeCount;
	int edgeCount;
	int bodyCount;
	
	double delta;
	double scale;
	
	std::vector< NodeToken * > nodes;
	std::vector< EdgeToken * > edges;
	std::vector< BodyToken * > bodies;
	BoundaryToken *boundary;

	std::map<int, BGMeshToken *> *bgmeshes;
	
	BodyToken *currentBody;
	LayerToken *currentLayer;
	SeedToken *currentSeed;
};

#endif /* MESH_MESHPARSER_H */
