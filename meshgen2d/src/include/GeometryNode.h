#if !defined( BL_GEOMETRYNODE_H )
#define BL_GEOMETRYNODE_H

#include "Node.h"
#include <algorithm>

class GeometryNode : public Node
{
public:
	GeometryNode() { tag = -1; x = y = 0; delta = 0; }
	GeometryNode(const int t, const double xc, const double yc) 
	{ tag = t; x = xc; y = yc; delta = 0; }
	~GeometryNode() { }
	
	GeometryNode& operator= (const GeometryNode& nd)
	{
		tag = nd.tag;
		x = nd.x;
		y = nd.y;
		
		return *this;
	}
	
	void setDelta( const double d )
	{
		delta = d;
	}
	
	int boundaryTag;
	double delta;
};

#endif /* BL_GEOMETRYNODE_H */
