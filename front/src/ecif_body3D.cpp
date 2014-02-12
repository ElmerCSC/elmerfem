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
Module:     ecif_body3D.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_body3D.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement3D.h"
#include "ecif_bodyElementLoop.h"
#include "ecif_model.h"
#include "ecif_renderer.h"
#include "ecif_userinterface.h"

 
// Constructors.

Body3D::Body3D()
  : Body()
{
}

//--Body-elements etc are not known
Body3D::Body3D(bodyGmtrType body_type, int ext_tag, char* name, colorIndices color)
  : Body(body_type, ext_tag, name, color)
{
}

//--Component ids are known, but not more.
Body3D::Body3D(ecif_Body_X& trx_body, bool add_default_layer)
  : Body(trx_body, add_default_layer)
{
}


Body3D::Body3D(bodyGmtrType body_type, int int_tag, int ext_tag,
               int nof_mesh_elements, int* mesh_element_ids)
               :Body(body_type, int_tag, ext_tag, NULL)
{
  maxNofMeshElements = nof_mesh_elements;
  nofMeshElements = nof_mesh_elements;
  meshElementIds = mesh_element_ids;
}


// Returns nof elements added
// NOTE: Not implmented yet!
int
Body3D::addAllPendingVertices(int layer)
{
  return 0;
}



// Method checks how face-loop is oriented.
int
Body3D::calcDirection()
{
  // Default is ccw.
  int result = 1;

  return result;
}


// NOTE: Not yet proper version for the 3D geometry!!!
// Return true if body is ok, otherwise false.
bool
Body3D::check()
{
  if (nofLayers == 0) return false;

  initName();

  // For a mesh body
  // ===============
  if (gmtrType == MESH_BODY) {
    return Body::check();

  // For a cad body
  // ==============
  } else {
    // Same as mesh body!
    // NOTE: We should do something real here!!!
    return Body::check();
  }
}

