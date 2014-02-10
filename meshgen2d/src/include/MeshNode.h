#if !defined( BL_MESHNODE_H )
#define BL_MESHNODE_H

#include "Node.h"
#include "GeometryNode.h"

enum node_type { NEUTRAL = 0, FIXED, CRYSTALNODE };
class MeshNode : public Node
{
public:
	MeshNode();
	MeshNode(const int t, const double xc, const double yc);
	MeshNode(const double xc, const double yc);
	MeshNode(const GeometryNode& nd);
	~MeshNode() { }
	
	void fix() { fixed = FIXED; }
	bool isFixed() { return fixed == FIXED; }
	void putOnCrysralIfNotFixed() { if( fixed != FIXED ) fixed = CRYSTALNODE; }
	bool isCrystal() { return fixed == CRYSTALNODE; }
	
private:
	node_type fixed;
};

#endif /* BL_MESHNODE_H */
