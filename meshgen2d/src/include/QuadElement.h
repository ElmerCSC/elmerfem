#if !defined( BL_QUADELEMENT )
#define BL_QUADELEMENT

#include "Element.h"

class QuadElement : public Element
{
public:
	QuadElement(Node *n1, Node *n2, Node *n3, Node *n4);

	virtual int elementType() const
	{
		return 400 + sz;
	}

private:
};

#endif /* BL_QUADELEMENT */
