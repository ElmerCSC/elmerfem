#include "QuadLayer.h"
#include "QuadElement.h"
#include "coreGeometry.h"
#include <iostream>

void getGridline(NodeMap& allNodes, std::vector<Node*> &strip, MeshNode** grid, int offset, double *places)
{
	double prevx = strip[0]->x, prevy = strip[0]->y, totLength = 0.0;
	int len = strip.size(), i;
	for (i = 0; i < len; i++)
	{
		totLength += distance(prevx, prevy, strip[i]->x, strip[i]->y);
		places[i] = totLength;
		prevx = strip[i]->x;
		prevy = strip[i]->y;
		grid[i * offset] = static_cast< MeshNode * >(strip[i]);
		allNodes[strip[i]->tag] = strip[i];
	}
	
	for (i = 0; i < len; i++)
		places[i] /= totLength;
}

// grid generation changed for curved bodies by Reino Ruusu 24.1.2000
void QuadLayer::
makeGrid( NodeMap& allNodes, const int m, const int n)
{
	int i,j;
	
	grid = new MeshNode*[m*n];
	
	// relative place of nodes on edges (0..1)
	double *places[4];
	places[0] = new double[n];
	places[1] = new double[m];
	places[2] = new double[n];
	places[3] = new double[m];
	
	std::vector<Node*> strip;
	
	// edge 0 is the column m
	edges[0]->exportNodes(strip, directions[0]);
	getGridline(allNodes, strip, &grid[m-1], m, places[0]);
	strip.clear();
	
	// edge 1 is the row n (reversed)
	edges[1]->exportNodes(strip, -directions[1]);
	getGridline(allNodes, strip, &grid[(n-1)*m], 1, places[1]);
	strip.clear();
	
	// edge 2 is the column 1 (reversed)
	edges[2]->exportNodes(strip, -directions[2]);
	getGridline(allNodes, strip, &grid[0], m, places[2]);
	strip.clear();
	
	// edge 3 is the row 1
	edges[3]->exportNodes(strip, directions[3]);
	getGridline(allNodes, strip, &grid[0], 1, places[3]);
	strip.clear();
	
	double *_offh = new double[8*m],
		*_offv = new double[8*n];
	double *xhll = _offh+0*m, *xhlr = _offh+1*m, *xhul = _offh+2*m, *xhur = _offh+3*m;
	double *yhll = _offh+4*m, *yhlr = _offh+5*m, *yhul = _offh+6*m, *yhur = _offh+7*m;
	double *xvll = _offv+0*n, *xvlr = _offv+1*n, *xvul = _offv+2*n, *xvur = _offv+3*n;
	double *yvll = _offv+4*n, *yvlr = _offv+5*n, *yvul = _offv+6*n, *yvur = _offv+7*n;
	
	// Precalculate offsets from corner points
	// rows
	for (i = 0; i < m; i++)
	{
		// LL
		xhll[i] = grid[i]->x - grid[0]->x;
		yhll[i] = grid[i]->y - grid[0]->y;
		// LR
		xhlr[i] = grid[i]->x - grid[m-1]->x;
		yhlr[i] = grid[i]->y - grid[m-1]->y;
		// UL
		xhul[i] = grid[(n-1)*m+i]->x - grid[(n-1)*m]->x;
		yhul[i] = grid[(n-1)*m+i]->y - grid[(n-1)*m]->y;
		// UR
		xhur[i] = grid[(n-1)*m+i]->x - grid[(n-1)*m+m-1]->x;
		yhur[i] = grid[(n-1)*m+i]->y - grid[(n-1)*m+m-1]->y;
	}
	
	// columns
	for (j = 0; j < n; j++)
	{
		// LL
		xvll[j] = grid[j*m]->x - grid[0]->x;
		yvll[j] = grid[j*m]->y - grid[0]->y;
		// UL
		xvul[j] = grid[j*m]->x - grid[(n-1)*m]->x;
		yvul[j] = grid[j*m]->y - grid[(n-1)*m]->y;
		// LR
		xvlr[j] = grid[j*m+m-1]->x - grid[m-1]->x;
		yvlr[j] = grid[j*m+m-1]->y - grid[m-1]->y;
		// UR
		xvur[j] = grid[j*m+m-1]->x - grid[(n-1)*m+m-1]->x;
		yvur[j] = grid[j*m+m-1]->y - grid[(n-1)*m+m-1]->y;
	}
	
	for( j = 1; j < (n - 1); ++j)
	{
		for( i = 1; i < (m - 1); ++i)
		{
			double pl = places[2][j],
				pr = places[0][j],
				pb = places[3][i],
				pt = places[1][i];
			
				/*   Calculate weights for interpolation
				wt = wl * pl + wr * pr
				wb = 1 - wt
				wr = wb * pb + wt * pt
				wl = 1 - wr
			*/
			double wt = ((pr - pl) * pb + pl) / (1.0 - (pr - pl) * (pt - pb));
			double wb = 1.0 - wt;
			double wr = wb * pb + wt * pt;
			double wl = 1.0 - wr;
			
			// Interpolate coordinates from the four gridline endpoints
			Node *left = grid[j*m], *right = grid[j*m+m-1],
				*bottom = grid[i], *top = grid[(n-1)*m+i];
			
			double x, y;
			x = ((wl * (left->x + wb * xhll[i] + wt * xhul[i]) +
				wr * (right->x + wb * xhlr[i] + wt * xhur[i])) +
				(wb * (bottom->x + wl * xvll[j] + wr * xvlr[j]) +
				wt * (top->x + wl * xvul[j] + wr * xvur[j]))) / 2.0;
			
			y = ((wl * (left->y + wb * yhll[i] + wt * yhul[i]) +
				wr * (right->y + wb * yhlr[i] + wt * yhur[i])) +
				(wb * (bottom->y + wl * yvll[j] + wr * yvlr[j]) +
				wt * (top->y + wl * yvul[j] + wr * yvur[j]))) / 2.0;
			
			grid[ j*m+i ] = new MeshNode( x, y );
			allNodes[grid[ j*m+i ]->tag] = grid[ j*m+i ];
		}
	}

	delete[] _offh;
	delete[] _offv;
	for (i = 0; i < 4; i++)
		delete[] places[i];
}

void QuadLayer::
discretize(NodeMap& fixedNodes, NodeMap& allNodes, std::list< Element* >& allElements)
{
	int i, j;
	
	int m = edges[1]->size();
	int n = edges[0]->size();
	
	QuadElement **elements = new QuadElement*[(n-1)*(m-1)];
	
	makeGrid( allNodes, m, n ); 
	for( j = 0; j < (n - 1); ++j)
		for( i = 0; i < (m - 1); ++i)
		{
			QuadElement *el = new QuadElement( grid[j*m+i+1], grid[(j+1)*m+i+1], grid[(j+1)*m+i], grid[j*m+i]);
			elements[j*(m-1)+i] = el;
			allElements.push_back(el);
		}
		
	// set boundary element relations
	std::vector<BoundaryElement*> bels;
	// edge 0 is the column m
	edges[0]->elements(bels, directions[0]);
	for (j = 0; j < n - 1; j++)
	{
		bels[j]->setLeft(elements[j*(m-1)+m-2]->elementId());
	}
	
	bels.clear();
	// edge 1 is the row n (reversed)
	edges[1]->elements(bels, -directions[1]);
	for (i = 0; i < m - 1; i++)
	{
		bels[i]->setRight(elements[(n-2)*(m-1)+i]->elementId());
	}
	
	bels.clear();
	// edge 2 is the column 1 (reversed)
	edges[2]->elements(bels, -directions[2]);
	for (j = 0; j < n - 1; j++)
	{
		bels[j]->setRight(elements[j*(m-1)]->elementId());
	}
	
	bels.clear();
	// edge 3 is the row 1
	edges[3]->elements(bels, directions[3]);
	for (i = 0; i < m - 1; i++)
	{
		bels[i]->setLeft(elements[i]->elementId());
	}
	
	delete [] grid;
	delete [] elements;
}
