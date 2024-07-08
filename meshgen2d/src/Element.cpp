#include "Element.h"
#include <iostream>
#include <algorithm>

#ifdef _NO_STD_MINMAX
	#include "minmaxpatch.h"
#endif

static int nextTag = 1;

Element::Element()
{
	newTag();
}

void Element::
newTag()
{
	elementTag = nextTag;
	++nextTag;
}

bool Element::isBoundaryConnector( const std::set< std::pair< int, int > > &links ) const
{
	int n = 0;
	
	for (int i = 0; i < sz - 1; i++)
		for (int j = i + 1; j < sz; j++)
		{
			int t1 = nodes[i]->tag, t2 = nodes[j]->tag;
			if ( nodes[i]->boundarynode && nodes[j]->boundarynode &&
			     links.find(std::make_pair(std::min(t1, t2), std::max(t1, t2))) != links.end() )
				n++;
		}
		
	return n > 1;
}

std::ostream& operator<< (std::ostream& o, const Element& A)
{
	int body = A.body;
	
	o << A.elementTag << ' ' << body << ' ' << A.elementType();
	
	for(int i = 0; i < A.sz; ++i)
	{
		o << ' ' << A.nodes[i]->tag;
	}
	o << '\n';

	return o;
}
