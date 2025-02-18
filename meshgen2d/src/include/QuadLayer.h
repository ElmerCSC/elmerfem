#if !defined( BL_QUADLAYER_H )
#define BL_QUADLAYER_H
#include <iostream>

#include "Layer.h"
#include "Border.h"
#include "MGError.h"

class Node;
class MeshNode;

typedef std::map< int, Node * > NodeMap;
typedef std::map< int, Node * >::iterator NodeMapIt;

class QuadLayer : public Layer
{
public:
	QuadLayer(const int t) : Layer(t) { grid = (MeshNode **) 0; }
	
	void sanityCheck()
	{
		if( edges.size() == 4 )
		{
			if(!(edges[0]->size() == edges[2]->size() && edges[1]->size() == edges[3]->size()))
			{
				blm_error("Grid mismatch");
			}
		}
		else
		{
			blm_error("Not a grid", tag);
		}
	}
	
	virtual void setBounds( Border *bd ) 
	{ 
		bounds = bd;
		// bounds is simply a single loop
		bounds->copyLoop( nodes, edges, directions );
	}
	
	virtual void initialize() { }

	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
protected:
	void makeGrid( NodeMap& allNodes, const int m, const int n );
	MeshNode **grid;
};
#endif /* BL_QUADLAYER_H */
