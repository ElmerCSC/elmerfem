#include "MeshNode.h"

static int nextTag = 1;

MeshNode::MeshNode()
{
	tag = nextTag;
	x = 0;
	y = 0;
	fixed = NEUTRAL;
	++nextTag;
}


MeshNode::MeshNode(const double xc, const double yc)
{
	tag = nextTag;
	x = xc;
	y = yc;
	fixed = NEUTRAL;
	++nextTag;
}

MeshNode::MeshNode(const int t, const double xc, const double yc)
{
	tag = t;
	x = xc;
	y = yc;
	fixed = NEUTRAL;
}

MeshNode::MeshNode(const GeometryNode& nd)
{
	tag = nd.tag;
	x = nd.x;
	y = nd.y;
	fixed = NEUTRAL;
	if (tag >= nextTag)
		nextTag = tag + 1;
}
