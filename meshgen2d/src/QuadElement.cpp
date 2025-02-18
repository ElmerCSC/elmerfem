#include "QuadElement.h"

QuadElement::QuadElement(Node *n1, Node *n2, Node *n3, Node *n4)	: Element() 
{
	sz = 4;
	nodes = new Node*[4];
	nodes[0] = n1;
	nodes[1] = n2;
	nodes[2] = n3;
	nodes[3] = n4;
}
