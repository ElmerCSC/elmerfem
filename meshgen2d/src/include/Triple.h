#if !defined( MESH_TRIPLE_H )
#define MESH_TRIPLE_H

class Node;
class Vertex;

class Triple
{
public:
	Triple() { }
	Node *from;
	Node *to;
	Vertex *vtx;
};


#endif /* MESH_TRIPLE_H */
