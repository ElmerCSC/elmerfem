#if !defined( MESH_CONNECT_H )
#define MESH_CONNECT_H

#include "Layer.h"

class Node;

class Vertex;
class PQ;
class Border;
class Triple;

#include <set>
#include <vector>
#include <list>

class Connect : public Layer
{
public:
	Connect(const int t) : Layer(t) { bg = NULL; }
	Connect(const int t, BGMesh *bgMesh) : Layer(t) { bg = bgMesh; }
	~Connect() { }

	virtual void initialize();
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
	void setBounds( Border *bd ) { bounds = bd; }
	
protected:
	
	void makeWorld();
	bool addSite( Vertex *vtx, Node *nd, bool externalOK, bool SplitOne = true, bool newPos = false );
	void recycle();
	void exportNodes( NodeMap& allNodes, std::list< Element* >& allElements );
	void triangulate();
	virtual void getVertices( std::vector< Vertex * >& v, const int count );
	Triple *getTriple();
	
	void removeBoundaryConnectors();
	
	Vertex *root;
	Node *world[4];
	
	std::vector< Node* > border;
	std::set< std::pair< int, int > > links;
	
	std::list< Vertex* > deleted;
	
	std::vector< Vertex* > newVertices;
	std::list< Vertex* > vertexStore;
	
	std::list< Triple * > hull;
	std::list< Triple * > tripleStore;
	
	std::list< Vertex* > allVertices;
	
	std::set< int, std::less<int> > existingTags;
	int firstNew;

	BGMesh *bg;
};

#endif /* MESH_CONNECT_H */
