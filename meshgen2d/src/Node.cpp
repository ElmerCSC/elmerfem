#include "Node.h"

std::ostream& operator<< (std::ostream& o, const Node& A)
{
	o << A.tag << ' ' << -1 << ' ' << A.x << ' ' << A.y << " 0\n";
	return o;
}
