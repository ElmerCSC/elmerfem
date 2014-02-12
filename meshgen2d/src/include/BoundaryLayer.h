#if !defined( BOUNDARY_LAYER_H )
#define BOUNDARY_LAYER_H

#include "Layer.h"

class Border;

#include <vector>

class BoundaryLayer : public Layer
{
public:
	BoundaryLayer(const int t) : Layer(t) { bg = NULL; }
	BoundaryLayer(const int t, BGMesh *bgMesh) : Layer(t) { bg = bgMesh; }
	~BoundaryLayer() { }

	virtual void initialize();
	void discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements);
	void setBounds( Border *bd ) { bounds = bd; }
	
protected:
	
	void exportNodes( NodeMap& allNodes, std::list< Element* >& allElements );
	std::vector< Node* > border;
		
	BGMesh *bg;
};

#endif /* BOUNDARY_LAYER_H */
