#include "VSVertex.h"
#include <iostream>
#include <math.h>
#include "coreGeometry.h"
#include <algorithm>

#ifdef _NO_STD_MINMAX
	#include "minmaxpatch.h"
#endif

//#define ORIGINAL

#if defined(ORIGINAL)

int VSVertex::
computeNewCoordinates( BGMesh& bg, double& nx, double& ny )
{
	double rho;
	double mx,my;
	double P;
	double Q;
	double ex,ey;
	double len;
	double dist;
	double rhohat;
	int i,pick;
	Node *a,*b;
	bool touched = false;
	
	len = 2 * radius;
	pick = 0;
	
	for(i = 0; i < 3; ++i)
	{
		VSVertex *vtx = static_cast<VSVertex *>( vertices[i] );
		if(vtx->isExternal() || vtx->isAccepted())
		{
			dist = distance(nodes[i]->x, nodes[i]->y,
				nodes[(i + 1) % 3]->x, 
				nodes[(i + 1) % 3]->y);
			if(dist < len)
			{
				len = dist;
				pick = i;
				touched = true;
			}
		}
	}
	
	if(!touched)
	{
		makeWaiting();
		return 0;
	}
	a = nodes[pick];
	b = nodes[(pick + 1)%3];
	
	dist = distance(vx,vy,vertices[pick]->vx,vertices[pick]->vy);
	ex = (vx - vertices[pick]->vx) / dist;
	ey = (vy - vertices[pick]->vy) / dist;
	
	mx = (a->x + b->x) / 2.0;
	my = (a->y + b->y) / 2.0;
	
	rho = bg.interpolate( mx, my );
	rho /= sqrt( 3.0 );
	
	Q = distance(vx,vy,mx,my);
	P = len / 2.0;
	rhohat = std::min(std::max(rho,P),((P*P + Q*Q) / (2.0 * Q)));
	double diff = rhohat*rhohat - P*P;
	if(diff < 0) dist = rhohat;
	else dist = rhohat + sqrt(diff);
	
	nx = mx + dist * ex;
	ny = my + dist * ey;
	return 1;
}
#else

const double htor = 1.0 / sqrt( 3.0 );
const double rtoh = 0.5 * sqrt( 3.0 );

#if 1
int VSVertex::
computeNewCoordinates( BGMesh& bg, double& lx, double& ly )
{
	double h = 0;
	double mx,my;
	double len;
	double dist;
	double nx, ny;
	int i,pick;
	Node *a,*b,*c;
	bool touched = false;
	
	nx = 0; ny = 0;
	len = 2 * radius;
	pick = 0;
	
	for(i = 0; i < 3; ++i)
	{
		VSVertex *vtx = static_cast<VSVertex *>( vertices[i] );
		if(vtx->isExternal() || vtx->isAccepted())
		{
			dist = distance(nodes[i]->x, nodes[i]->y,
				nodes[(i + 1) % 3]->x, 
				nodes[(i + 1) % 3]->y);
			if(dist < len)
			{
				len = dist;
				pick = i;
				touched = true;
			}
		}
	}
	
	if(!touched)
	{
		makeWaiting();
		return 0;
	}
	
	a = nodes[pick];
	b = nodes[(pick + 1)%3];
	c = nodes[(pick + 2)%3];
	
	mx = (a->x + b->x) / 2.0;
	my = (a->y + b->y) / 2.0;
	
	h = bg.interpolate( mx, my );
	
	ny = b->x - a->x;
	nx = b->y - a->y;
	nx = -nx;
	nx /= len;
	ny /= len;
	
	double error, A, B, C = h, D = h;
	
	do
	{
		A = (len * len + C * C - D * D) / (2.0 * len);
		
		if (A > C || (len - A) > D)
		{
			lx = vx;
			ly = vy;
			return 5;
		}
		
		B = sqrt(C * C - A * A);
		
		lx = a->x + A * ny + B * nx;
		ly = a->y - A * nx + B * ny;
		
		double vx1 = (lx + a->x) / 2.0, vy1 = (ly + a->y) / 2.0,
			vx2 = (lx + b->x) / 2.0, vy2 = (ly + b->y) / 2.0;
		
		double h1 = bg.interpolate(vx1, vy1);
		double h2 = bg.interpolate(vx2, vy2);
		
		double e1 = fabs(distance(lx, ly, a->x, a->y) - h1) / h1, e2 = fabs(distance(lx, ly, b->x, b->y) - h2) / h2;
		
		C = 0.75 * C + 0.25 * h1;
		D = 0.75 * D + 0.25 * h2;
		
		error = e1 + e2;
	} while (error > 1e-2);
	
	if (C > radius && D > radius)
	{
		lx = vx;
		ly = vy;
	}
	
	return 4;
}

#else
int VSVertex::
computeNewCoordinates( BGMesh& bg, double& lx, double& ly )
{
	double h = 0;
	double rho;
	double mx,my;
	double P;
	double Q;
	double ex,ey;
	double len;
	double dist;
	double rhohat;
	double nx, ny;
	int i,pick;
	Node *a,*b,*c;
	bool touched = false;
	
	nx = 0; ny = 0;
	len = 2 * radius;
	pick = 0;
	
	for(i = 0; i < 3; ++i)
	{
		VSVertex *vtx = static_cast<VSVertex *>( vertices[i] );
		if(vtx->isExternal() || vtx->isAccepted())
		{
			dist = distance(nodes[i]->x, nodes[i]->y,
				nodes[(i + 1) % 3]->x, 
				nodes[(i + 1) % 3]->y);
			if(dist < len)
			{
				len = dist;
				pick = i;
				touched = true;
			}
		}
	}
	
	if(!touched)
	{
		makeWaiting();
		return 0;
	}
	
	a = nodes[pick];
	b = nodes[(pick + 1)%3];
	c = nodes[(pick + 2)%3];
	
	mx = (a->x + b->x) / 2.0;
	my = (a->y + b->y) / 2.0;
	
	h = bg.interpolate( mx, my );
	//  rho = h / std::sqrt( 3.0 );
	rho = h * htor;
	
	Q = distance(vx,vy,mx,my);
	P = len / 2.0;
	
	double minrad = (P*P + Q*Q) / (2.0 * Q);
	
	if(minrad < rho || minrad < P)
	{
		lx = vx;
		ly = vy;
		return 2;
	}
	
	ny = b->x - a->x;
	nx = b->y - a->y;
	nx = -nx;
	nx /= len;
	ny /= len;
	
	if(P > rho)
	{	
		lx = mx + P * nx;
		ly = my + P * ny;
		return 3;
	}
	
	//  dist = 0.5 * std::sqrt( 3.0 ) * h; 
	dist = h * rtoh;
	lx = mx + dist * nx;
	ly = my + dist * ny;
	
	return 1;
}
#endif

#endif
