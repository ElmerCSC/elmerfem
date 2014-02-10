#if !defined( BL_TRIANGLEELEMENT )
#define BL_TRIANGLEELEMENT

#include "Element.h"

class TriangleElement : public Element
{
public:
	TriangleElement( const int t ) : Element( t )  
	{ sz = 3; nodes = new Node*[3]; }
	TriangleElement(Node *n1, Node *n2, Node *n);

	virtual int elementType() const
	{
		return 300 + sz;
	}
	
private:
};

#endif /* BL_TRIANGLEELEMENT */
