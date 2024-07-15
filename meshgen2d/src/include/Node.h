#if !defined( BL_NODE_H )
#define BL_NODE_H

#include <fstream>

class Node
{
public:
	Node() { tag = 0; x = y = 0; boundarynode = false; }
	
	int tag;
	bool boundarynode;
	double x,y;
};

std::ostream& operator<< (std::ostream& o, const Node& A);

#endif /* BL_NODE_H */
