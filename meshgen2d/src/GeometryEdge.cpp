#include "GeometryEdge.h"
#include "BGMesh.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include "coreGeometry.h"

static int nextTag = 1;

GeometryEdge::GeometryEdge( int nSegments )
{
  tag = nextTag;

  segments = nSegments;

  boundaryTag = 0;
  type = OUTER;
  ++nextTag;
}

GeometryEdge::GeometryEdge( const int t, int nSegments )
{
  tag = t;

  segments = nSegments;

  boundaryTag = 0;
  type = OUTER;
  ++nextTag;
}

void GeometryEdge::
exportNodes(std::vector<Node*>& strip, int direction)
{
  if( direction < 0 )
  {
    std::reverse_copy( nodes.begin(), nodes.end(), std::back_inserter( strip ) );
  }
  else
  {
    std::copy( nodes.begin(), nodes.end(), std::back_inserter( strip ) );
  }
}

void GeometryEdge::
elements(std::vector< BoundaryElement * >& strip, int direction)
{
  int i, len = bels.size();

  if( len == 0 ) // Not constructed yet
  {
    int nlen = nodes.size();
    for( i = 0; i < (nlen - 1); ++i )
    {
      bels.push_back( new BoundaryElement( boundaryTag, nodes[i], nodes[i+1] ) );
    }

    len = bels.size();
  }

  if( direction < 0 )
  {
    for( i = 0; i < len; ++i ) bels[i]->flip( true );
    std::reverse_copy( bels.begin(), bels.end(), std::back_inserter( strip ) );
  }
  else
  {
    for( i = 0; i < len; ++i ) bels[i]->flip( false );
    std::copy( bels.begin(), bels.end(), std::back_inserter( strip ) );
  }
}

void GeometryEdge::
midNodes( Node*& a, Node*& b, int dir )
{
  int len = nodes.size();
  int ta, tb;

  ta = len / 2 - 1;
  tb = ta + 1;

  if( dir == 1 )
  {
    a = nodes[ta];
    b = nodes[tb];
  }
  else
  {
    a = nodes[tb];
    b = nodes[ta];
  }
}

double GeometryEdge::interpolate(double x, double y)
{
  double value = 0.0;
  for (std::vector<BGMesh*>::iterator it = bgMeshes.begin(); it != bgMeshes.end(); it++)
  {
    double h = (*it)->interpolate(x, y);
    if (value == 0.0 || h < value)
      value = h;
  }

  return value;
}

void GeometryEdge::
discretize( NodeMap &allNodes )
{
  if (segments == 0)
  {
    int len = dots.size() - 1;
    for (int i = 0; i < len; i++)
    {
      discretizeGradedSegment(allNodes, dots[i], dots[i+1]);
      if (i < len - 1)
        nodes.pop_back();
    }
  }
  else
  {
    int len = dots.size() - 1, i;
    if (segments < len) segments = len;

    double totLength = 0.0;
    double *lengths = new double[len];

    for (i = 0; i < len; i++)
    {
      lengths[i] = distance(dots[i]->x, dots[i]->y, dots[i+1]->x, dots[i+1]->y);
      totLength += lengths[i];
    }

    int *segs = new int[len];
    for (i = 0; i < len; i++)
      segs[i] = 1;

    int remaining = segments - len;
    while (remaining > 0)
    {
      int longest = 0;
      double length = lengths[0] / segs[0];
      for (i = 1; i < len; i++)
      {
        double l = lengths[i] / segs[i];
        if (l > length)
        {
          length = l;
          longest = i;
        }
      }

      segs[longest]++;
      remaining--;
    }

    delete [] lengths;

    for (i = 0; i < len; i++)
    {
      discretizeConstantSegment(segs[i], allNodes, dots[i], dots[i+1]);
      if (i < len - 1)
        nodes.pop_back();
    }

    delete [] segs;
  }
}

void GeometryEdge::
discretizeConstantSegment( int nSeg, NodeMap& allNodes, GeometryNode *from, GeometryNode *to )
{
  Node *first = allNodes[from->tag];
  Node *last = allNodes[to->tag];

  double dx = last->x - first->x;
  double dy = last->y - first->y;
  double length = sqrt(dx*dx + dy*dy);
  double sx = dx / length;
  double sy = dy / length;
  double step = length / nSeg;
  double bx = first->x;
  double by = first->y;

  from->setDelta( step );
  to->setDelta( step );

  nodes.push_back( first );
  for( int i = 1; i < nSeg; ++i )
  {
    MeshNode *nd = new MeshNode(bx+i*step*sx, by+i*step*sy);
    nd->boundarynode = true;
    nodes.push_back( nd );
  }
  nodes.push_back( last );
}

void GeometryEdge::
discretizeGradedSegment( NodeMap& allNodes, GeometryNode *from, GeometryNode *to )
{
  int i;
  double lineLength = distance( from->x, from->y, to->x, to->y );
  double grading[2];

  double sx = (to->x - from->x) / lineLength;
  double sy = (to->y - from->y) / lineLength;

  double limit = lineLength;
  double at = 0.0;
  double predicted = interpolate(from->x, from->y);

  std::vector< double > lengths;
  int step = 0;

  while( 1 )
  {
    double prev;
    do
    {
      prev = predicted;
      double center = at + predicted / 2.0;
      // Avoid interpolating outside the domain
      if (center > lineLength)
        center = lineLength;
      predicted = 0.5 * (predicted + interpolate(from->x + center * sx, from->y + center * sy));
    }
    while( fabs((prev - predicted)/lineLength) > 0.0000001 );

    at += predicted;

    if (at > limit - 0.01 * predicted)
      break;

    lengths.push_back( at );
    ++step;
  }

  double mult = lineLength / at;
  for( i = 0; i < step; ++i)
  {
    lengths[i] *= mult;
  }

  Node *first = allNodes[from->tag];
  Node *last = allNodes[to->tag];
  nodes.push_back( first );
  for( i = 0; i < step; ++i )
  {
    MeshNode *nd = new MeshNode( lengths[i]*sx + from->x, lengths[i]*sy + from->y );
    nd->boundarynode = true;
    nodes.push_back( nd );
  }
  nodes.push_back( last );
}

std::ostream& operator<< (std::ostream& o, const GeometryEdge& A)
{
  const int len = A.nodes.size();
  o << A.tag << ' ' << len << ' ';
  for(int i = 0; i < len; ++i)
  {
    o << A.nodes[i]->tag << ' ';
  }
  o << std::endl;
  return o;
}

void GeometryEdge::
getGridPoints(double ox, double oy, double cellsize, std::vector<int> &x, std::vector<int> &y, std::vector<double> &delta)
{
  int i, len = dots.size();
  for (i = 0; i < len - 1; i++)
  {
    double ax = dots[i]->x, ay = dots[i]->y;
    double bx = dots[i+1]->x, by = dots[i+1]->y;

    double ad = MAP(dots[i]->delta), bd = MAP(dots[i+1]->delta);

    int h = (int)((ax - ox) / cellsize + 0.5);
    int v = (int)((ay - oy) / cellsize + 0.5);

    x.push_back(h);
    y.push_back(v);
    delta.push_back(ad);

    double dx = bx - ax;
    double dy = by - ay;
    double d = sqrt(dx * dx + dy * dy);
    dx /= d;
    dy /= d;
    double dd = bd - ad;

    int endh = (int)((bx - ox) / cellsize + 0.5);
    int endv = (int)((by - oy) / cellsize + 0.5);

    if (endh == h && endv == v)
      continue;

    int hstep = endh > h ? 1 : -1;
    int vstep = endv > v ? 1 : -1;

    if (abs(endh - h) > abs(endv - v))
    {
      int pos = 0, num = abs(endv - v), div = abs(endh - h);

      while ((h += hstep) != endh)
      {
        pos += num;
        if (pos > div / 2)
        {
          v += vstep;
          pos -= div;
        }

        double nx = ox + h * cellsize;
        double ny = oy + v * cellsize;
        double nd = ad + ((nx - ax) * dx + (ny - ay) * dy) / d * dd;

        x.push_back(h);
        y.push_back(v);
        delta.push_back(nd);
      }
    }
    else
    {
      int pos = 0, num = abs(endh - h), div = abs(endv - v);

      while ((v += vstep) != endv)
      {
        pos += num;
        if (pos > div / 2)
        {
          h += hstep;
          pos -= div;
        }

        double nx = ox + h * cellsize;
        double ny = oy + v * cellsize;
        double nd = ad + ((nx - ax) * dx + (ny - ay) * dy) / d * dd;

        x.push_back(h);
        y.push_back(v);
        delta.push_back(nd);
      }
    }

    if (i == len - 2)
    {
      x.push_back(endh);
      y.push_back(endv);
      delta.push_back(bd);
    }
  }
}
