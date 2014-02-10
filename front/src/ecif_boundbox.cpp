/*  
   ElmerFront - A graphical user interface of Elmer software
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_boundbox.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include <algorithm>

#include "ecif_boundbox.h"
#include "ecif_geometry.h"

BoundBox::BoundBox()
{
  x1 = y1 = z1 = MAX_RANGE;
  x2 = y2 = z2 = MIN_RANGE;
}
 

BoundBox::BoundBox(double min_value, double max_value)
{
  x1 = y1 = z1 = min_value;
  x2 = y2 = z2 = max_value;
}

BoundBox::BoundBox(RangeVector rv)
{
  //x-range
  if (rv[0] <= rv[1]) { 
    x1=rv[0];
    x2=rv[1];
  }
  else {
    x1=rv[1];
    x2=rv[0];
  }
  //y-range
  if (rv[2] <= rv[3]) {
    y1=rv[2];
    y2=rv[3];
  }
  else {
    y1=rv[3];
    y2=rv[2];
  }
  //z-range
  if (rv[4] <= rv[5]) {
    z1=rv[4];
    z2=rv[5];
  }
  else {
    z1=rv[5];
    z2=rv[4];
  }
}


bool 
BoundBox::contains(BoundBox* other_box)
{
  //If this-box contains *other_box*
  return ( x1 <= (other_box->x1 + GAP_TOLERANCE) &&
           x2 >= (other_box->x2 - GAP_TOLERANCE) &&
           y1 <= (other_box->y1 + GAP_TOLERANCE) && 
           y2 >= (other_box->y2 - GAP_TOLERANCE) &&
           z1 <= (other_box->z1 + GAP_TOLERANCE) &&
           z2 >= (other_box->z2 - GAP_TOLERANCE) );
}


//Extend range by point
void
BoundBox::extendByPoint(GcPoint* point)
{
  Point3 p;
  
  point->getPoint(p);
  
  extendByPoint(p);
}


//Extend range by point
void
BoundBox::extendByPoint(Point3 p)
{
  if ( x1 > p[0] )
    x1 = p[0];

  if ( x2 < p[0] )
    x2 = p[0];

  if ( y1 > p[1] )
    y1 = p[1];

  if ( y2 < p[1] )
    y2 = p[1];

  if ( z1 > p[2] )
    z1 = p[2];

  if ( z2 < p[2] )
    z2 = p[2];
}


//Extend range by new coordinates
//Note: argument *rv* must be ordered properly!
void
BoundBox::extendByRange(RangeVector rv)
{
  // Update maximum-box values
  setMaximumPair(x1, x2, rv[0], rv[1]);
  setMaximumPair(y1, y2, rv[2], rv[3]);
  setMaximumPair(z1, z2, rv[4], rv[5]);
}


// NOTE! Works currently only for 2D!!
// Method gets a bounding box axis as corner values.
// Corners are numbered ccw starting from the near sw-corner.
// z=0: 1,2,3,4; z=1: 5,6,7,8
//
void
BoundBox::getBoxAxis(int crn1, int crn2, double* start, double* end1, double* end2)
{
  // Deafaults for 2D
  start[2] = end1[2] = 0.0;
  end2[0] = end2[1] = end2[2] = 0.0;

  // Start corner
  switch (crn1) {
  case 1:
    start[0] = x1; start[1] = y1; break;
  case 2:
    start[0] = x1; start[1] = y2; break;
  case 3:
    start[0] = x2; start[1] = y2; break;
  case 4:
    start[0] = x2; start[1] = y1; break;
  }
  // End corner
  switch (crn2) {
  case 1:
    end1[0] = x1; end1[1] = y1; break;
  case 2:
    end1[0] = x1; end1[1] = y2; break;
  case 3:
    end1[0] = x2; end1[1] = y2; break;
  case 4:
    end1[0] = x2; end1[1] = y1; break;
  }
}


// Get box xyz-dimensions
void
BoundBox::getSize(double& dim_x, double& dim_y, double& dim_z)
{
  dim_x = x2 - x1;
  dim_y = y2 - y1;
  dim_z = z2 - z1;
}


void
BoundBox::getRangeVector(RangeVector rv)
{
  rv[0] = x1; rv[1] = x2;
  rv[2] = y1; rv[3] = y2;
  rv[4] = z1; rv[5] = z2;
}


ostream& 
BoundBox::output(ostream& out) 
{
  out << "(" << x1 << ", " << x2 << ") ";
  out << "(" << y1 << ", " << y2 << ") ";
  out << "(" << z1 << ", " << z2 << ") ";
  
  return out;
}


bool 
BoundBox::overlap(BoundBox* other_box)
{
  using namespace std;

  double xd = min(x2, other_box->x2) - max(x1, other_box->x1);
  double yd = min(y2, other_box->y2) - max(y1, other_box->y1);
  double zd = min(z2, other_box->z2) - max(z1, other_box->z1);

  double mindif = min(min(xd,yd),zd);
  double maxdif = max(max(xd,yd),zd);

  return ( (mindif > -1 * GAP_TOLERANCE)  && (maxdif > GAP_TOLERANCE) );

}


// Restrict minimum-pair box by new coordinates
// This type of box is maintaining two smallest values
// for vertex-values for each coordinate
// First of values in each pair is always the smaller.
void
BoundBox::restrictByRange(RangeVector rv)
{ 
  // Update box-values.
  setMinimumPair(x1, x2, rv[0], rv[1]);
  setMinimumPair(y1, y2, rv[2], rv[3]);
  setMinimumPair(z1, z2, rv[4], rv[5]);
}


// NOTE: We suppouse that c1,c2 is ordered pair!!!
void
BoundBox::setMaximumPair( double& c1, double& c2, double r1, double r2)
{
  if ( isGreater(r1, r2) ) {
    double tmp = r1;
    r1 = r2;
    r2 = tmp;
  }

  // Lower value
  if ( r1 != -NSVD && r1 != NSVD && r1 < c1) 
    c1 = r1;

  // Upper value
  if ( r2 != -NSVD && r2 != NSVD && r2 > c2) 
    c2 = r2;
}


// NOTE: We suppouse that c1,c2 is ordered pair!!!
void
BoundBox::setMinimumPair(double& c1, double& c2, double r1, double r2)
{
  if ( isGreater(r1, r2) ) {
    double tmp = r1;
    r1 = r2;
    r2 = tmp;
  }

  // At first we select the lowest value.
  double tmp = r1;
  if ( isEqual(c1,NSVD) || isLess(r1,c1) )  {
    tmp = c1;
    c1 = r1;
  }

  // If c1 was changed and it originally was unitialized
  // we have to change the value of *tmp*
  if ( isEqual(tmp,NSVD) ) {
    tmp = r2;
  }

  // Now we have to select the lowest of the remaining three, such
  // that it is different from *c1*.
  if ( !isEqual(tmp,c1) && isLess(tmp,c2) ) {
    c2 = tmp;
  }

  if ( !isEqual(r2,c1) && isLess(r2,c2) ) {
    c2 = r2;
  }
}
