#include "BGMesh.h"
#include "BGVertex.h"
#include <math.h>
#include <algorithm>
#include "coreGeometry.h"

void BGGridMesh::
initialize(std::vector< GeometryNode* >& nodes, std::vector< GeometryEdge* >& edges)
{
	initialized = true;
	
	std::cout << "Calculating background mesh" << std::endl;
	
	int i;
	int len = nodes.size();
	
	double mean = MAP(nodes[0]->delta);
	double mindelta = nodes[0]->delta;
	double minx = nodes[0]->x, maxx = minx, miny = nodes[0]->y, maxy = miny;
	
	for( i = 1; i < len; ++i )
	{
		mean += MAP(nodes[i]->delta);
		
		if (nodes[i]->delta < mindelta)
			mindelta = nodes[i]->delta;
		
		if (nodes[i]->x < minx)
			minx = nodes[i]->x;
		if (nodes[i]->x > maxx)
			maxx = nodes[i]->x;
		if (nodes[i]->y < miny)
			miny = nodes[i]->y;
		if (nodes[i]->y > maxy)
			maxy = nodes[i]->y;
	}
	
	mean /= len;
	
	width = maxx - minx;
	height = maxy - miny;
	ox = minx;
	oy = miny;
	
	cellsize = mindelta * cellsize;
	nh = (int)ceil(width / cellsize) + 1;
	width = (nh - 1) * cellsize;
	nv = (int)ceil(height / cellsize) + 1;
	height = (nv - 1) * cellsize;
	
	nh += 2;
	nv += 2;
	ox -= cellsize;
	oy -= cellsize;
	width += 2.0 * cellsize;
	height += 2.0 * cellsize;
	
	grid = new double[nh * nv];
	bool *mask = new bool[nh * nv];
	
	for (i = 0; i < nh * nv; i++)
	{
		grid[i] = mean;
		mask[i] = false;
	}
	
	for (i = 0; i < len; i++)
	{
		int h = (int)((nodes[i]->x - ox) / cellsize + 0.5);
		int v = (int)((nodes[i]->y - oy) / cellsize + 0.5);
		grid[v*nh+h] = MAP(nodes[i]->delta);
		mask[v*nh+h] = true;
	}
	
	len = edges.size();
	for (i = 0; i < len; i++)
	{
		std::vector<int> x, y;
		std::vector<double> delta;
		edges[i]->getGridPoints(ox, oy, cellsize, x, y, delta);
		for (int j = 0; j < x.size(); j++)
		{
			grid[y[j]*nh+x[j]] = delta[j];
			mask[y[j]*nh+x[j]] = true;
		}
	}
	
	std::vector<std::pair<int, int> > stack;
	stack.push_back(std::make_pair(0, 0));
	while (!stack.empty())
	{
		std::pair<int, int> p = stack.back();
		stack.pop_back();
		
		int h = p.first, v = p.second;
		mask[v*nh+h] = true;
		
		if (h > 0 && !mask[v*nh+h-1])
			stack.push_back(std::make_pair(h-1, v));
		if (h < nh - 1 && !mask[v*nh+h+1])
			stack.push_back(std::make_pair(h+1, v));
		if (v > 0 && !mask[(v-1)*nh+h])
			stack.push_back(std::make_pair(h, v-1));
		if (v < nv - 1 && !mask[(v+1)*nh+h])
			stack.push_back(std::make_pair(h, v+1));
	}
	
	int niters = 0;
	double error = 1e100, last;
	do
	{
		niters++;
		last = error;
		error = 0.0;
		
		int i, j;
		
		for (int c = 0; c < 200; c++)
		{
			for (j = 0; j < nv; j++)
			{
				for (i = 0; i < nh; i++)
				{
					if (mask[j*nh+i]) continue;
					
					double v = 0.0;
					int n = 0;
					
					if (i > 0)
					{
						v += grid[j*nh+i-1];
						n++;
					}
					if (i < nh - 1)
					{
						v += grid[j*nh+i+1];
						n++;
					}
					if (j > 0)
					{
						v += grid[(j-1)*nh+i];
						n++;
					}
					if (j < nv - 1)
					{
						v += grid[(j+1)*nh+i];
						n++;
					}
					
					grid[j*nh+i] = v / n;
				}
			}
		}
		
		for (j = 0; j < nv; j++)
		{
			for (i = 0; i < nh; i++)
			{
				if (mask[j*nh+i]) continue;
				
				double residual = 0.0;
				
				double v = grid[j*nh+i];
				if (i > 0)
					residual += grid[j*nh+i-1] - v;
				if (i < nh - 1)
					residual += grid[j*nh+i+1] - v;
				if (j > 0)
					residual += grid[(j-1)*nh+i] - v;
				if (j < nv - 1)
					residual += grid[(j+1)*nh+i] - v;
				
				residual /= cellsize;
				
				double e = fabs(residual);
				if (e > error)
					error = e;
			}
		}
	} while (error > 0.00001);
}

double BGGridMesh::
interpolate( const double ix, const double iy )
{
	int h = (int)((ix - ox) / cellsize);
	int v = (int)((iy - oy) / cellsize);
	
	if (h < 0) h = 0;
	if (h > nh - 1) h = nh - 1;
	if (v < 0) v = 0;
	if (v > nv - 1) v = nv - 1;
	
	double xoff = (ix - (ox + h * cellsize)) / cellsize,
		yoff = (iy - (oy + v * cellsize)) / cellsize;
	
	return UNMAP((1.0 - xoff) * (1.0 - yoff) * grid[v*nh+h] +
	             xoff * (1.0 - yoff) * grid[v*nh+h+1] +
	             (1.0 - xoff) * yoff * grid[(v+1)*nh+h] +
	             xoff * yoff * grid[(v+1)*nh+h+1]);
}
