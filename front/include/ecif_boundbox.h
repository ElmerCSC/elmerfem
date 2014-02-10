/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_boundbox.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   A Base class for geometry object's bounding box. 

************************************************************************/

#ifndef _ECIF_BOUNDBOX_
#define _ECIF_BOUNDBOX_

#include "ecif_def.h"


class BoundBox
{     
public:
  BoundBox();
  BoundBox(double min_value, double max_value);
  BoundBox(RangeVector range_values);
  bool contains(BoundBox* other_box);
  void extendByPoint(GcPoint* point);
  void extendByPoint(Point3 point);
  void extendByRange(RangeVector);
  void getBoxAxis(int crn1, int crn2, double* start, double* end1, double* end2);
  void getSize(double& dim_x, double& dim_y, double& dim_z);
  void getRangeVector(RangeVector rv);
  ostream& output(ostream& out);
  bool overlap(BoundBox* other_box);
  void restrictByRange(RangeVector);
  
private:
  double x1, x2, y1, y2, z1, z2;
  void setMaximumPair(double& c1, double& c2, double r1, double r2);
  void setMinimumPair(double& c1, double& c2, double r1, double r2);
} ;

#endif
