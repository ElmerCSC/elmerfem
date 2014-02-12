#if !defined( MESH_BODY_H )
#define MESH_BODY_H

#include <vector>
#include "Layer.h"

class Body
{
public:
	Body(const int t, bool p) { tag = t; parabolic = p; }
	void addLayer( Layer* l) { layers.push_back( l ); }
	void discretize(NodeMap& fixedNodes,  NodeMap& allNodes, std::list< Element* >& allElements)
	{
		std::list< Element* > elements;
		int len = layers.size();
		int i;
		for( i = 0; i < len; ++i )
		{
			layers[i]->sanityCheck();
			layers[i]->discretize( fixedNodes, allNodes, elements );
		}
		
		std::list< Element* >::iterator e;
		for (e = elements.begin(); e != elements.end(); e++)
			(*e)->setBody(tag);
		
		allElements.insert(allElements.end(), elements.begin(), elements.end());
	}

	bool isParabolic() const { return parabolic; }
	
	int tag;
protected:
	std::vector< Layer * > layers;
	bool parabolic;
};

#endif /* MESH_BODY_H */
