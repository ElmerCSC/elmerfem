#if !defined( MESH_BGMESH_H )
#define MESH_BGMESH_H

#include "Connect.h"

class BGMesh
{
public:
	BGMesh() { initialized = false; }
	
	virtual double interpolate( const double ix, const double iy ) = 0;
	virtual void initialize(std::vector< GeometryNode* >& nodes, std::vector< GeometryEdge* >& edges) = 0;
	
	bool isInitialized() { return initialized; }
	
protected:
	bool initialized;
};

class BGTriangleMesh : public Connect, public BGMesh
{
public:
	BGTriangleMesh() : Connect(-1) { }
	
	virtual void discretize(NodeMap& fixedNodes, NodeMap& allNodes, 
		std::list< Element* >& allElements) { }
	
	virtual void initialize(std::vector< GeometryNode* >& nodes, std::vector< GeometryEdge* >& edges);
	
	void getVertices( std::vector< Vertex * >& v, const int count );
	
	virtual double interpolate( const double ix, const double iy );
	void dump();
};

class BGGridMesh : public BGMesh
{
public:
	BGGridMesh(double scale) { cellsize = scale; }
	
	virtual double interpolate( const double ix, const double iy );
	virtual void initialize(std::vector< GeometryNode* >& nodes, std::vector< GeometryEdge* >& edges);
	
protected:
	int nh, nv;
	double cellsize;
	double ox, oy, width, height;
	double *grid;
};

#endif /* MESH_BGMESH_H */
