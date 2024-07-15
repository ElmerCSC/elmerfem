#include "BGVertex.h"
#include "GeometryNode.h"
#include "coreGeometry.h"
#include <iostream>

void BGVertex::
initInterpolation( )
{
	int i;
	
	double a[9];
	double b[3];
	
	for ( i = 0; i < 3; ++i )
	{
		int i1 = (i + 1) % 3, i2 = (i + 2) % 3;
		a[i] = nodes[i1]->x * nodes[i2]->y - nodes[i2]->x * nodes[i1]->y;
		a[3 + i] = nodes[i1]->y - nodes[i2]->y;
		a[6 + i] = nodes[i2]->x - nodes[i1]->x;
	}
	
	double det = a[0] + a[1] + a[2];
	
	int n = 0;
	double avg = 0.0;
	
	for( i = 0; i < 3; ++i )
	{
		int t = nodes[i]->tag;
		if( t < 0 ) continue;
		GeometryNode *nd = static_cast< GeometryNode * >( nodes[i] );
		b[i] = MAP(nd->delta);
		avg += b[i];
		n++;
	}
	
	avg /= n;
	for (i = 0; i < 3; ++i)
	{
		if (nodes[i]->tag < 0)
			b[i] = avg;
	}
	
	for ( i = 0; i < 3; ++i )
	{
		coeff[i] = (a[3*i] * b[0] + a[3*i+1] * b[1] + a[3*i+2] * b[2]) / det;
	}
}

double BGVertex::
interpolate( const double tx, const double ty )
{
	bool found = false;
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
	
	while( !found )
	{
		Vertex *test = curr;
		for(i = 0; i < 3; ++i)
		{
			if((curr->vertices[i]) && (curr->vertices[i] != prev) &&
				orientation(tx,ty,x[s[i]],y[s[i]],x[e[i]],y[e[i]]))
			{
				prev = curr;
				curr = curr->vertices[i];
				
				x[0] = tx; y[0] = ty;
				x[1] = curr->nodeAt(0)->x; y[1] = curr->nodeAt(0)->y;
				x[2] = curr->nodeAt(1)->x; y[2] = curr->nodeAt(1)->y;
				x[3] = curr->nodeAt(2)->x; y[3] = curr->nodeAt(2)->y;
				break;
			}
		}
		if( test == curr ) found = true;
	}
	
	BGVertex *vtx = static_cast< BGVertex * >( curr );
	double val =  UNMAP(vtx->coeff[0] + vtx->coeff[1] * tx + vtx->coeff[2] * ty);
	
	return val;
}
