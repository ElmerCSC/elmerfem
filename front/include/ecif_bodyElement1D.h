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
Module:     ecif_bodyelement1D.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   One dimensional bodyelemet ie. vertex. 

************************************************************************/

#ifndef _ECIF_BODYELEMENT1D_
#define _ECIF_BODYELEMENT1D_

#include "ecif_bodyElement.h"


class BodyElement1D : public BodyElement
{   
public:
  BodyElement1D();
  BodyElement1D(int tag);
  BodyElement1D(int tag, GcPoint* point);
  BodyElement1D(GcPoint* point);
  BodyElement1D(ecif_Element_X& trx_element);
  BodyElement1D(ecif_Vertex_X& trx_vertex);
  ~BodyElement1D();

  int addMeshElement(int elem_id, short direction) { nodeId = elem_id; return 1;}
  bool checkOuterBoundaries() { return false;}
  void convertVertexIds2Vertices() {}
  BodyElement* createElement(int nof_vertices, int* vertex_ids, ecif_geometryType gt);
  int compareOrientation(BodyElement* oe) { return 1;}
  void draw(Renderer* renderer, flagName geometry_type, int body_id);
  virtual void draw(Renderer* renderer, flagName geometry_type, int body_id, int direction) {draw(renderer,geometry_type,body_id);}
  virtual void draw(Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop) {draw(renderer,geometry_type,body_id);}
  matchType findCommonBoundary(BodyElement* other_element,
                               BodyElement*& common );
  int findMeshBorderNodes(int buf_size, int* ids_buffer);
  int getMeshElementId(int index);
  int getNofVertices() { return 1; }
  BodyElementList* getOuterBoundary() { return NULL;}
  double getParamArea(Geometry* gp) { return 0.0;}
  ParamValues* getParamValues(Geometry* gp) { return NULL;}
  bool hasInside(BodyElement* other_element) { return false;}
  void init(char* be_name = NULL);
  static void initClass(Model* model);
  bool isBemBoundary() { return false; }
  bool isInnerBoundary() {return false;}
  bool isOk() {return true;}
  void markActiveMeshObjects(bool*& active_object_flags) {};
  //ostream& output_emf(ostream& out, short indent_size, short indent_level, bool isOnSymmAxis);
  ostream& output_mif(ostream& out);
  GcPoint* param2Point(double u_p, double v_p = 0) { return NULL;}  
  ParamPair* point2Param(GcPoint* p) { return NULL;}

protected:
  int getLastTag() {return last_tag;}
  void initLabelData();
  int newTag() { return ++last_tag;}
  int nofCurrentMeshHSources;
  void setLastTag(int ltag) { last_tag = ltag; }

  static int last_tag;
  int nodeId;
};
 
  
#endif

