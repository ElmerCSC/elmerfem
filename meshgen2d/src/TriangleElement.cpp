#include "TriangleElement.h"

TriangleElement::TriangleElement(Node *n1, Node *n2, Node *n3)  : Element() 
{
	sz = 3;
	nodes = new Node*[3];
	nodes[0] = n1;
	nodes[1] = n2;
	nodes[2] = n3;
}
