#include "Connect.h"
#include <algorithm>
#include "Vertex.h"
#include "TriangleElement.h"
#include "Border.h"
#include "BGMesh.h"
#include "Triple.h"
#include "PQ.h"
#include "MGError.h"

#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#ifdef _NO_STD_MINMAX
	#include "minmaxpatch.h"
#endif

void Connect::
makeWorld()
{
	double maxx, maxy, minx, miny;
	int len = border.size();
	
	maxx = minx = border[0]->x;
	maxy = miny = border[0]->y;
	for( int i = 1; i < len; ++i)
	{
		double tx = border[i]->x;
		double ty = border[i]->y;
		
		if( minx > tx) minx = tx;
		else if( maxx < tx) maxx = tx;
		
		if( miny > ty) miny = ty;
		else if( maxy < ty) maxy = ty;
	}
	
	double dx = maxx - minx;
	double dy = maxy - miny;
	
	double diam = std::max( dx, dy );
	
	double cx = (maxx + minx) / 2.0;
	double cy = (maxy + miny) / 2.0;
	
	world[0] = new MeshNode( -1, cx - dx, cy - dy );
	world[1] = new MeshNode( -2, cx + dx, cy - dy );
	world[2] = new MeshNode( -3, cx + dx, cy + dy );
	world[3] = new MeshNode( -4, cx - dx, cy + dy );
	
	std::vector< Vertex * > box;
	getVertices( box, 2 );
	
	box[0]->reset( world[0], world[1], world[2] );
	box[1]->reset( world[2], world[3], world[0] );
	
	box[0]->setVertices((Vertex *)0, (Vertex *)0, box[1]);
	box[1]->setVertices((Vertex *)0, (Vertex *)0, box[0]);
	
	root = box[0];
}

void Connect::initialize()
{
	bounds->collectGeometryNodes( nodes );
	bounds->collectGeometryEdges( edges, directions );
	if(!bg->isInitialized()) bg->initialize( nodes, edges );
}

void Connect::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	triangulate();

	removeBoundaryConnectors();
	
	exportNodes( allNodes, allElements );
}


void Connect::
triangulate()
{
	int i;
	int len;
	
	bounds->collectNodesAndPairs( border, links );
	
	makeWorld();
	
	for (i = 0; i < fixedGeometryNodes.size(); i++)
	{
		MeshNode *nd = static_cast< MeshNode* >( fixedGeometryNodes[i] );
		Vertex *vtx = root->firstAcceptableFound( nd->x, nd->y);
		addSite( vtx, nd, true );
		nd->fix();
		recycle();
	}
	
	len = border.size();
	int *indirect = new int[len];
	for( i = 0; i < len; ++i)
	{
		indirect[i] = i;
	}
	std::random_shuffle(indirect, indirect + len);
	
//	std::cout << "Inserting border" << std::endl;
	for( i = 0; i < len; ++i)
	{
		int ind = indirect[i];
		Vertex *vtx = root->firstAcceptableFound( border[ind]->x, border[ind]->y);
		addSite( vtx, border[ind], true );
		MeshNode *nd = static_cast< MeshNode* >( border[ind] );
		nd->fix();
		recycle();
	}
	
	delete [] indirect;
	
	// Verify borders (no recovery!)
//	std::cout << "Verifying borders" << std::endl;
	
	len = border.size();
	
	Vertex *verify = root;
	for(i = 0; i < len; ++i)
	{
		int t1 = border[i]->tag, t2 = border[(i + 1) % len]->tag;
		if( links.find( std::make_pair( std::min(t1, t2), std::max(t1, t2) ) ) != 
			links.end() )
		{ 
			verify = verify->vertexOwning( border[i], border[(i + 1) % len] );
			if(!verify)
			{
				blm_error("Initial triangulation is incomplete!",
					"Please readjust your mesh parameters.");
			}
		}
	} 
//	std::cout << "Verifying borders completed" << std::endl;
	
	//
	Vertex *vtx = root->vertexWith( border[0], border[1] );
	
	vtx->labelBodies( tag, links );
}

void Connect::
exportNodes( NodeMap& allNodes, std::list< Element* >& allElements )
{	
	std::list< Vertex* >::iterator vxIt;
	for( vxIt = allVertices.begin(); vxIt != allVertices.end(); ++vxIt )
	{
		Vertex *vtx = (*vxIt);
		if( !(vtx->isDeleted() || vtx->isExternal()))
		{
			vtx->newTag();
			allElements.push_back( vtx );
			
			for(int i = 0; i < 3; ++i)
			{
				allNodes[ vtx->nodeAt(i)->tag ] = vtx->nodeAt(i);
			}
		}
	}
	
//	std::cout << "Nodes OK!" << std::endl;
	
	std::vector< BoundaryElement * > bels;
	bounds->collectBoundaryElements( bels );
	
//	std::cout << "BoundaryElement I OK!" << std::endl;
	
	int i, len;
	Vertex *v = root;
	len = bels.size();
	
	for( i = 0; i < len; ++i )
	{
		v = bels[i]->setHolder( v );
	}
//	std::cout << "BoundaryElement II OK!" << std::endl;
}

bool Connect::
addSite( Vertex *vx, Node *nd, bool externalOK, bool splitOne, bool force)
{
	int i, count;
	Vertex *vtx;
	Triple *trip, *tmp;
	std::list< Triple * > tripleQueue;
	
	hull.clear();
	
	for( i = 2; i >= 0; --i )
	{
		trip = getTriple();
		trip->from = vx->nodeAt(i);
		trip->to = vx->nodeAt((i+1)%3);
		trip->vtx = vx->vertices[i];
		
		tripleQueue.push_back( trip );
	}
	
	vx->makeDeleted();
	deleted.push_back( vx );
	
	while( tripleQueue.size() > 0 )
	{
		trip = tripleQueue.back();
		tripleQueue.pop_back();
		
		vtx = trip->vtx;
		if(vtx == (Vertex *)0)
		{
			hull.push_back( trip );
			continue;
		}
		if( vtx->isDeleted() )
		{
		/*
		If a vertex is deleted, we must have visited it during this loop.
		Therefore, there must be a triple with reversed direction in the stack list.
		We must remove both triples.
			*/
			tmp = tripleQueue.back();
			if(tmp->from == trip->to && tmp->to == trip->from)
			{
				tripleQueue.pop_back();
			}
			else
			{
				tripleQueue.pop_front();
			}
		}
		else
		{
		/*
		Return status:
		=  1 OK, destroyed
		=  0 OK, not destroyed
		= -1 NOT OK, external element destroyed
		
		 If retval = -1, abort the operation, clear buffers and indicate the same thing 
		 for the upper level.
			*/
			int retval = vtx->isDestroyedBy(nd->x, nd->y, externalOK);
			if(retval == 1)
			{
				vtx->makeDeleted();
				deleted.push_back( vtx );
				
				for(i = 0; i < 3; ++i) if(vtx->nodeAt(i) == trip->from) break;
				
				tmp = getTriple();
				tmp->from = vtx->nodeAt((i+1)%3);
				tmp->to = vtx->nodeAt((i+2)%3);
				tmp->vtx = vtx->vertices[(i+1)%3];
				
				tripleQueue.push_back( tmp );
				
				tmp = getTriple();
				tmp->from = vtx->nodeAt(i);
				tmp->to = vtx->nodeAt((i+1)%3);
				tmp->vtx = vtx->vertices[i];
				
				tripleQueue.push_back( tmp );
			}
			else if(retval == 0 || force)
			{
				hull.push_back( trip );
			}
			else
			{
				std::list< Vertex* >::iterator deletedIt;
				for( deletedIt = deleted.begin(); deletedIt != deleted.end();  ++deletedIt )
				{
					(*deletedIt)->unDelete();
				}
				deleted.clear();
				return false;
			}
		}
	}
	if (!externalOK && !splitOne && deleted.size() == 1)
		return true;
	
	// These parts should be separated!	
	count = hull.size();
	
	getVertices( newVertices, count );
	
	std::list<Triple *>::const_iterator hit;
	for(hit = hull.begin(), i = 0; i < count; ++i, ++hit)
	{
		tmp = *hit;
		newVertices[i]->reset( nd, tmp->from, tmp->to );
	}
	
	for(hit = hull.begin(), i = 0; i < count; ++i, ++hit)
	{
		tmp = *hit;
		vtx = newVertices[i];
		vtx->setVertices(newVertices[(i+count-1)%count], 
			tmp->vtx,
			newVertices[(i+1)%count]);
		
		if(tmp->vtx != (Vertex *)0)
		{
			if(tmp->vtx->nodeAt(0) == tmp->to)
				tmp->vtx->replaceVertexWith(0, vtx);
			else if(tmp->vtx->nodeAt(1) == tmp->to)
				tmp->vtx->replaceVertexWith(1, vtx);
			else if(tmp->vtx->nodeAt(2) == tmp->to)
				tmp->vtx->replaceVertexWith(2, vtx);
		}
	}
	
	root = newVertices.back();
	
	return true;		 
}

void Connect::
recycle()
{
	std::copy( deleted.begin(), deleted.end(), std::back_inserter( vertexStore ) );
	deleted.erase( deleted.begin(), deleted.end() );
	newVertices.erase( newVertices.begin(), newVertices.end() );
	std::copy( hull.begin(), hull.end(), std::back_inserter( tripleStore ) );
	hull.erase( hull.begin(), hull.end() );
}

void Connect::
getVertices( std::vector< Vertex * >& v, const int count )
{
	int i;
	
	for( i = 0; i < count; ++i )
	{
		if( vertexStore.empty() ) break; 
		v.push_back( vertexStore.back() );
		vertexStore.pop_back();
	}
	
	for( ; i < count; ++i )
	{
		Vertex *vtx = new Vertex;
		v.push_back( vtx );
		allVertices.push_back( vtx );
	}	
}

Triple * Connect::
getTriple()
{
	Triple *tmp;
	if( tripleStore.empty() ) return new Triple;
	else
	{
		tmp = tripleStore.back();
		tripleStore.pop_back();
	}
	return tmp;
}

#include <queue>
#include <set>
#include "coreGeometry.h"

class Segment
{
public:
	Node *a, *b;
	double length;
};

class SegCompare
{
public:
	bool operator() (const Segment *s0, const Segment *s1)
	{
		return s0->length < s1->length;
	}
};

void Connect::
removeBoundaryConnectors()
{
	const double sqrt3 = 1.7320508075688772935274463415059;
	
	pq actives;

	double ref;
	std::list< Vertex* >::iterator vxIt;
	for( vxIt = allVertices.begin(); vxIt != allVertices.end(); ++vxIt )
	{
		Vertex *vtx = (*vxIt);
		if( !(vtx->isDeleted() || vtx->isExternal() ))
		{
			if( vtx->isBoundaryConnector(links) )
			{
				vtx->rightSize( 1.0 );
				actives.insert( vtx );
			}
		}
	}
	
	while( actives.size() > 0 )
	{
		Vertex *vtx = actives.first();

		if (!vtx->isBoundaryConnector(links))
			continue;
		
		int i;
		double x, y;

		for (i = 0; i < 3; i++)
		{
			int t0 = vtx->nodeAt(i)->tag, t1 = vtx->nodeAt((i+1)%3)->tag;
			if (links.find(std::make_pair(std::min(t0, t1), std::max(t0, t1))) == links.end())
			{
				double dx0 = vtx->nodeAt(i)->x - vtx->nodeAt((i+2)%3)->x;
				double dy0 = vtx->nodeAt(i)->y - vtx->nodeAt((i+2)%3)->y;
				double dx1 = vtx->nodeAt((i+1)%3)->x - vtx->nodeAt((i+2)%3)->x;
				double dy1 = vtx->nodeAt((i+1)%3)->y - vtx->nodeAt((i+2)%3)->y;

				// r1 = sqrt(dx0*dx0+dy0*dy0)
				// r2 = sqrt(dx1*dx1+dy1*dy1)
				// theta1 = atan(dy0/dx0)
				// theta2 = atan(dy1/dx1)
				// (x + iy) = sqrt(r1*r2)*e^(i*(theta1+theta2)/2)

				double x2y2 = dx0 * dx1 - dy0 * dy1; // x*x-y*y
				double xy = dx0 * dy1 + dx1 * dy0; // 2*x*y

				y = sqrt(0.5 * (sqrt(x2y2 * x2y2 + xy * xy) - x2y2));
				x = 0.5 * xy / y;

				if (dx0 * x + dy0 * y < 0.0)
				{
					x = -x;
					y = -y;
				}

				x = vtx->nodeAt((i+2)%3)->x + x;
				y = vtx->nodeAt((i+2)%3)->y + y;

				break;
			}
		}

		MeshNode *nd = new MeshNode( x, y ); //vtx->vx, vtx->vy );
		
		if( addSite( vtx, nd, false, true, false) == true )
		{
			actives.removeRelevants( deleted );
			
			int len = newVertices.size();
			for( int i = 0; i < len; ++i )
			{
				newVertices[i]->setBody( 1 );
				if( newVertices[i]->isBoundaryConnector(links) )
				{
					newVertices[i]->rightSize( 1.0 );
					actives.insert( newVertices[i] );
				}
			}
		}
		
		recycle();
	}

	return;
}
