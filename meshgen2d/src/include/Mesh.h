#if !defined( MESH_MESH_H )
#define MESH_MESH_H

#include "MeshParser.h"
#include <map>
#include <fstream>

class GeometryNode;
class GeometryEdge;
class Node;
class Element;
class BoundaryElement;
class Body;
class BGMesh;

typedef std::map< int, Body * > BodyMap;
typedef std::map< int, Body * >::iterator BodyMapIt;
typedef std::map< int, GeometryNode * > GeometryNodeMap;
typedef std::map< int, GeometryNode * >::iterator GeometryNodeMapIt;
typedef std::map< int, GeometryEdge * > GeometryEdgeMap;
typedef std::map< int, GeometryEdge * >::iterator GeometryEdgeMapIt;
typedef std::map< int, Node * > NodeMap;
typedef std::map< int, Node * >::iterator NodeMapIt;

class Mesh
{
public:
	Mesh();
	void convertEdgeFormat( MeshParser& parser );
	void discretize();
	void createMiddleNodes();
	
	void outputFlat(std::ofstream& o);
	
	void output( std::ofstream& o ); 
	
	void outputHeader( std::ofstream& o ); 
	void outputNodes( std::ofstream& o ); 
	void outputElements( std::ofstream& o );
	void outputBoundary( std::ofstream& o );
	
protected:
	BodyMap bodies;
	GeometryNodeMap geometryNodes;
	GeometryEdgeMap geometryEdges;
	NodeMap fixedNodes;
	NodeMap meshNodes;
	std::list< Element * > elements;
	std::vector< BoundaryElement * > boundaryElements;
};

#endif /* MESH_MESH_H */
