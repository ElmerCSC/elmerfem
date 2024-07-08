#if !defined( BL_LAYER_H )
#define BL_LAYER_H

#include <map>
#include <vector>
#include <list>

#include "GeometryNode.h"
#include "GeometryEdge.h"
#include "Element.h"
#include "Border.h"

typedef std::map< int, Node * > NodeMap;
typedef std::map< int, Node * >::iterator NodeMapIt;

class Border;

class Layer
{
public:
	Layer(const int t) { tag = t; }
	~Layer() {}
	
	void copyNodesAndEdgesFrom(std::vector< GeometryNode* >& nds, std::vector< GeometryEdge* >& eds)
	{
		nodes = nds; edges = eds;
	}
	
	void collectBoundaryElements(std::vector< BoundaryElement* >& elements)
	{
		bounds->collectBoundaryElements(elements);
	}
	
	void addFixedNode(GeometryNode *node)
	{
		nodes.push_back(node);
		MeshNode *nd = new MeshNode(*node);
		fixedGeometryNodes.push_back(nd);
	}
	
	virtual void sanityCheck() { }

	virtual void initialize() = 0;

	virtual void discretize(NodeMap& fixedNodes, NodeMap& allNodes, 
		std::list< Element* >& allElements) = 0;
	virtual void setBounds( Border *bd ) { bounds = bd; }
	
protected:
	int tag;
	std::vector< GeometryNode* > nodes; 
	std::vector< Node* > fixedGeometryNodes;
	std::vector< GeometryEdge* > edges;
	std::vector< int > directions;
	Border *bounds;
};

#endif /* BL_LAYER_H */
