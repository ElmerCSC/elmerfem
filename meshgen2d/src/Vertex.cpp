#include <math.h>

#include "Vertex.h"
#include "MeshNode.h"
#include "coreGeometry.h"

#include <algorithm>
#ifdef _NO_STD_MINMAX
	#include "minmaxpatch.h"
#endif

static int nextTag = 1;

Vertex::Vertex()	: TriangleElement( -1 )
{
	tag = nextTag;  
	heap = -1; 
	++nextTag;
}

Vertex::Vertex(Node *n1, Node *n2, Node *n3)  : TriangleElement( -1 )
{
	tag = nextTag; 	
	++nextTag;
	
	reset(n1, n2, n3);
}

void Vertex::
reset(Node *n1, Node *n2, Node *n3)
{
	nodes[0] = n1; nodes[1]= n2; nodes[2] = n3;
	
	vertexLocation(n1->x ,n1->y, n2->x ,n2->y, n3->x, n3->y, &vx, &vy);
	radius = distance(vx, vy, n2->x, n2->y);
	
	deleted = false;
	heap = -1;
	body = 0;
}

void Vertex::
setVertices(Vertex *v1, Vertex *v2, Vertex *v3)
{
	vertices[0] = v1; vertices[1] = v2; vertices[2] = v3;
}

Vertex* Vertex::
firstAcceptableFound(double tx, double ty)
{
	int i;
	Vertex *curr = this;
	Vertex *prev = (Vertex *)0;
	
	int s[3] = {2,3,1};
	int e[3] = {1,2,3};
	
	double x[4],y[4];
	
	x[0] = tx; y[0] = ty;
	x[1] = nodes[0]->x; y[1] = nodes[0]->y;
	x[2] = nodes[1]->x; y[2] = nodes[1]->y;
	x[3] = nodes[2]->x; y[3] = nodes[2]->y;
	
	do //while(!inCircle(x,y))
	{
		for(i = 0; i < 3; ++i)
		{
			if((curr->vertices[i]) && (curr->vertices[i] != prev) &&
				orientation(tx,ty,x[s[i]],y[s[i]],x[e[i]],y[e[i]]))
			{
				prev = curr;
				curr = curr->vertices[i];
				
				x[0] = tx; y[0] = ty;
				x[1] = curr->nodes[0]->x; y[1] = curr->nodes[0]->y;
				x[2] = curr->nodes[1]->x; y[2] = curr->nodes[1]->y;
				x[3] = curr->nodes[2]->x; y[3] = curr->nodes[2]->y;
				break;
			}
		}	
	} while (i < 3);
	
	return curr;
}

int Vertex::
isDestroyedBy(double tx, double ty, bool extOK)
{
	double x[4],y[4];
	
	x[0] = tx; y[0] = ty;
	x[1] = nodes[0]->x; y[1] = nodes[0]->y;
	x[2] = nodes[1]->x; y[2] = nodes[1]->y;
	x[3] = nodes[2]->x; y[3] = nodes[2]->y;
	if(inCircle(x,y))
	{
		if(isExternal())
		{
			if(!extOK)
				return -1;
			else
				return 1;
		}
		else return 1;
	}
	return 0;
}

Vertex* Vertex::
vertexOwning( Node* a, Node* b)
{
	bool found = false;
	int i, holder = -1;
	Vertex *curr = this;
	Vertex *prev = (Vertex *)0;
	
	double tx = (a->x + b->x)/2.0;
	double ty = (a->y + b->y)/2.0;
	
	int s[3] = {2,3,1};
	int e[3] = {1,2,3};
	
	double x[4],y[4];
	
	x[0] = tx; y[0] = ty;
	x[1] = nodes[0]->x; y[1] = nodes[0]->y;
	x[2] = nodes[1]->x; y[2] = nodes[1]->y;
	x[3] = nodes[2]->x; y[3] = nodes[2]->y;
	
	// What if we already have Node* a?
	for(i = 0; i < 3; ++i)
	{
		if(nodes[i] == a)
		{
			found = true;
			holder = i;
			if(nodes[(i+1)%3] == b) return curr;
			if(nodes[(i - 1 + 3)%3] == b) return curr;
			break;
		}
	}	
	
	while(!found)
	{
		for(i = 0; i < 3; ++i)
		{
			if((curr->vertices[i]) && (curr->vertices[i] != prev) &&
				orientation(tx,ty,x[s[i]],y[s[i]],x[e[i]],y[e[i]]))
			{
				prev = curr;
				curr = curr->vertices[i];
				
				x[0] = tx; y[0] = ty;
				x[1] = curr->nodes[0]->x; y[1] = curr->nodes[0]->y;
				x[2] = curr->nodes[1]->x; y[2] = curr->nodes[1]->y;
				x[3] = curr->nodes[2]->x; y[3] = curr->nodes[2]->y;
				break;
			}
		}
		
		for(i = 0; i < 3; ++i)
		{
			if(curr->nodes[i] == a)
			{
				found = true;
				holder = i;
				if(curr->nodes[(i+1)%3] == b) return curr;
				if(curr->nodes[(i - 1 + 3)%3] == b) return curr;
				break;
			}
		}	
	}
	// So, curr owns Node* a.
	// Where is Node* b?
	// Let's go around Node* a using a linked list.
	
	Vertex* owner = (Vertex*)0;
	Vertex* start = curr;
	do {
		curr = curr->vertices[holder];
		for(i = 0; i < 3; ++i)
		{
			if(curr->nodes[i] == a)
			{
				holder = i;
				break;
			}
		}
		if(curr->nodes[(i+1)%3] == b) return curr;
		if(curr->nodes[(i - 1 + 3)%3] == b) return curr;
	} while(curr != start);
	
	if((curr == start) && (owner == (Vertex*)0)) return (Vertex*)0;
	return owner;
}

Vertex* Vertex::
vertexWith( Node* a, Node* b)
{
	bool found = false;
	int i;
	Vertex *curr = this;
	Vertex *prev = (Vertex *)0;
	
	double tx = (a->x + b->x)/2.0;
	double ty = (a->y + b->y)/2.0;
	
	int s[3] = {2,3,1};
	int e[3] = {1,2,3};
	
	double x[4],y[4];
	
	x[0] = tx; y[0] = ty;
	x[1] = nodes[0]->x; y[1] = nodes[0]->y;
	x[2] = nodes[1]->x; y[2] = nodes[1]->y;
	x[3] = nodes[2]->x; y[3] = nodes[2]->y;
	
	// What if we already have Node* a?
	for(i = 0; i < 3; ++i)
	{
		if(nodes[i] == a)
		{
			if(nodes[(i+1)%3] == b) return curr;
			if(nodes[(i - 1 + 3)%3] == b) return curr->vertices[(i - 1 + 3)%3];
		}
	}	
	
	while(!found)
	{
		for(i = 0; i < 3; ++i)
		{
			if((curr->vertices[i]) && (curr->vertices[i] != prev) &&
				orientation(tx,ty,x[s[i]],y[s[i]],x[e[i]],y[e[i]]))
			{
				prev = curr;
				curr = curr->vertices[i];
				
				x[0] = tx; y[0] = ty;
				x[1] = curr->nodes[0]->x; y[1] = curr->nodes[0]->y;
				x[2] = curr->nodes[1]->x; y[2] = curr->nodes[1]->y;
				x[3] = curr->nodes[2]->x; y[3] = curr->nodes[2]->y;
				break;
			}
		}
		
		for(int j = 0; j < 3; ++j)
		{
			if(curr->nodes[j] == a && curr->nodes[(j+1)%3] == b)
			{
				found = true;
				break;
			}
			else if(curr->nodes[j] == a && curr->nodes[((j-1)+3)%3] == b)
			{
				curr = curr->vertices[(j - 1 + 3)%3];
				found = true;
				break;
			}
		}
		
		if (!found && i == 3) // This is wrong!!!!
			return NULL;
	}
	return curr;
}

void Vertex::
labelBodies( const int bd, std::set< std::pair< int, int > >& bounds )
{
	int i;
	if(body == bd) return;
	
	body = bd;

	for(i = 0; i < 3; ++i)
	{
		int t1 = nodes[i]->tag;
		int t2 = nodes[(i+1)%3]->tag;
		
		if(bounds.find( std::make_pair(std::min(t1,t2),std::max(t1,t2)) ) == bounds.end())
		{
			vertices[i]->labelBodies( bd, bounds );
		}
	}
}
