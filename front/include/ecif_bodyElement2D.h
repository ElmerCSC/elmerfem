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
Module:     ecif_bodyelement2D.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Derived from Bodyelement-class. 
            Handles 2D-bodyelements.

************************************************************************/

#ifndef _ECIF_BODYELEMENT2D_
#define _ECIF_BODYELEMENT2D_

#include "ecif_bodyElement.h"

                                 
class BodyElement2D : public BodyElement
{     
public:
  BodyElement2D();
  BodyElement2D(int tag);
  BodyElement2D(int v1_id, int v2_id);
  BodyElement2D(int nof_vertices, int* vertex_ids);
  BodyElement2D(int v1_id, int v2_id, Geometry* pGmtr, int code = 0, char* name = 0);
  BodyElement2D(int v1_id, int v2_id, ecif_EdgeGeometry_X* params);
  BodyElement2D(ecif_Element_X& trx_element);
  BodyElement2D(int tag, int parent1_tag, int parent2_tag,int nof_fem_elements);
  BodyElement2D(int parent1_tag, int parent2_tag, int nof_fem_elements);
  ~BodyElement2D();

  void addCoveringElement(BodyElement* se, beStatus se_stat);
  void checkBoundaryDiscretization(int mesh_index);
  bool checkOuterBoundaries();
  int compareOrientation(BodyElement* oe);
  BodyElement* createElement(int v1_id, int v2_id, ecif_geometryType gt);
  //void draw(Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop);
  int findMeshBorderNodes(int buf_size, int* ids_buffer);
  int getMifGeometryTag(int index) const;
  int getNofMifGeometries() const;
  double getParamArea(Geometry* gp);
  ParamValues* getParamValues(Geometry* gp);
  BodyElement* getSubElement(int index);
  void init(char* be_name = NULL);
  static void initClass(Model* model);
  //bool isBemBoundary();
  //bool isInnerBoundary();
  bool isOk();
  bool isOnSameAxis(GcPoint& p1, GcPoint& p2);
  matchType matchToLinear(BodyElement* be2, BodyElement*& common);
  matchType matchToNurbs(BodyElement* be2, BodyElement*& common);
  ostream& output_mif(ostream& out);
  void setMifTag(int& next_tag);

protected:
  int calcDirection(BodyElement* sub_element);
  BodyElementList* findOuterBoundary();
  int getLastTag() { return last_tag; }
  bool isBefore(BodyElement* se1, int dir1,
                BodyElement* se2, int dir2,
                Geometry* gp);
  int newTag() { return ++last_tag;}
  void setLastTag(int ltag) { last_tag = ltag;}

  int nofMifTags;
  int* mifTags;

  static int last_tag;

} ;

#endif
