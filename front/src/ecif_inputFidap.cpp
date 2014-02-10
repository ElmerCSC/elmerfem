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
Module:     ecif_inputFidap.cpp
Language:   C++
Date:       27.04.99
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_body2D.h"
#include "ecif_body3D.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputFidap.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_timer.h"
#include "ecif_userinterface.h"

extern char read_buffer[];

// Indices to the special reorder table for element types
const int fidapReorderTableIndex_203 = 0;
const int fidapReorderTableIndex_510 = 1;
const int fidapReorderTableIndex_718 = 2;
const int fidapReorderTableIndex_808 = 3;
const int fidapReorderTableIndex_809 = 4;
const int fidapReorderTableIndex_827 = 5;

// We have to reorder Fidap nodes to match Elmer numbering
const int fidapReorderTable[][MAX_NOF_NODES] = {

  {0,2,1}, // Line 203

  {0,2,5,9,1,4,3,6,7,8}, // Tetra 510

  // NOTE: 718 is delivered to Elmer as 706, but here is the complete reorder! (ok?)
  {0,2,5,12,14,17,1,4,3,6,7,8,10,11,9,13,16,15},     // Wedge 718

  {0,1,3,2,4,5,7,6},     // Brick 808

  {0,1,3,2,4,5,7,6, 8},  // Brick 809 (1 in center)

  // NOTE: 827 is delivered to Elmer as 820, but here is the complete reorder!
  {0,2,8,6,18,20,26,24,1,5,7,3,9,11,17,15,19,23,25,21,10,14,16,12,4,22,13} // Brick 827
};


InputFidap::InputFidap(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
  for (short i = 0; i < 1 + MAX_NOF_ELEM_CODES; i++)
    elementCodeCounters[i] = 0;

  nofElementGroups = 0;
  dataDimension = ECIF_ND;
}


InputFidap::~InputFidap()
{
}


Body*
InputFidap::findExistingBody(FidapElementGroupInfo* gi)
{
  int index = 0;
  while (true) {
    Body* body = model->getBody(index++);
    if (body==NULL) break;
    if ( LibFront::ncEqual((char*)body->getName(), gi->entityName) ) {
      return body;
    }
  }

  return NULL;
}


BodyElement*
InputFidap::findExistingBoundary(FidapElementGroupInfo* gi)
{
  int index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    if ( LibFront::ncEqual((char*)be->getName(), gi->entityName) ) {
      return be;
    }
  }

  return NULL;
}


Rc
InputFidap::findMaxExternalIds(int& max_ext_nd_id, int& max_ext_el_id)
{
  enum fidapKeywordId keyword = FIDAP_NO_ID;

  bool count_only = true;

	while ( !infile.eof() ) {

    // Read new line
    readFileLine(infile, read_buffer);

    // Check if it was an interesting keyword line
    keyword = next_is_keyword_line(read_buffer);

    if ( keyword == FIDAP_NO_ID)
      continue;

    // Read element group
    // ==================
    if ( keyword == FIDAP_ELEMENT_GROUPS_ID ) {

      int element_group_counter = 0;

      while ( element_group_counter < nofElementGroups &&
              !infile.eof()
              ) {

        element_group_counter++;

        // Read one group
        FidapElementGroupInfo gi;

        // Check that group is ok
        if ( !readElementGroupInfo(gi) )
          return ECIF_ERROR;

        if ( !setElementType(gi) )
          return ECIF_ERROR;

        // Read group elements
        // ===================
        int element_counter = 0;

        readElements(gi, element_counter, max_ext_el_id, count_only);
      }
    }

    // Read nodes
    // ==========
    else if ( keyword == FIDAP_NODAL_COORDINATES_ID ) {
      int node_counter = 0;
      readNodes(nofNodes, node_counter, max_ext_nd_id, count_only);
    }

  } // while !eof

  return ECIF_OK;

}


// NOTE: Fidap mesh dimension is read from the header and
// stored in the class attribute 'dimension'
//
enum ecif_modelDimension
InputFidap::findMeshModelDimension()
{
  ecif_modelDimension box_dim = ECIF_ND;
  ecif_modelDimension dim = ECIF_ND;

#if 0
  RangeVector rv;
  bool x_is_const = false;
  bool y_is_const = false;
  bool z_is_const = false;

  meshBoundingBox.getRangeVector(rv);

  // Find model dimension by boundary box
  if ( isEqual(rv[0], rv[1]) ) {
    x_is_const = true;
  }

  if ( isEqual(rv[2], rv[3]) ) {
    y_is_const = true;
  }

  if ( isEqual(rv[4], rv[5]) ) {
    z_is_const = true;
  }

  if ( !(x_is_const || y_is_const || z_is_const) ) {
    box_dim = ECIF_3D;
  } else if ( !( x_is_const || y_is_const) && z_is_const ) {
    box_dim = ECIF_2D;
  } else if ( !(x_is_const && y_is_const && z_is_const) ) {
    box_dim = ECIF_1D;
  } else {
    box_dim = ECIF_ND;
  }
#endif

  box_dim = modelDimension;

  // Find actual dimension
  dim = findModelDimension(box_dim);

  return dim;
}


enum fidapKeywordId
InputFidap::next_is_keyword_line(char* line_buffer)
{
  if ( 0 == strncmp("NODAL COORDINATES", line_buffer, 17) )
    return FIDAP_NODAL_COORDINATES_ID;

  else if ( 0 == strncmp("ELEMENT GROUPS", line_buffer, 14) )
    return FIDAP_ELEMENT_GROUPS_ID;

  else if ( 0 == strncmp("VERSION", line_buffer, 7) )
    return FIDAP_VERSION_ID;

  else
    return FIDAP_NO_ID;
}


Rc
InputFidap::readElementGroup(int& body_counter, int& bndr_counter,
                             int& element_counter, int& max_ext_elem_id,
                             bool count_only)
{
  int element_group_counter = 0;

  while ( element_group_counter < nofElementGroups &&
          !infile.eof()
          ) {

    element_group_counter++;

    // Read one group
    FidapElementGroupInfo gi;

    // Check that group is ok
    if ( !readElementGroupInfo(gi) )
      return ECIF_ERROR;

    if ( !setElementType(gi) )
      return ECIF_ERROR;

    gi.parentTag = NO_INDEX;

    // Check elment type bulk/boundary
    model->checkMeshElementType(gi.elementType,
                                gi.isBulk,
                                gi.isBoundary,
                                gi.isEdge);

    // Bulk elements --> add elements to a body
    // =============
    if (gi.isBulk) {

      nofBulkElements += gi.nofElements;

      // Try to find an existing body
      // ----------------------------
      Body* body = findExistingBody(&gi);

      // Create a new body
      // -----------------
      if (body == NULL) {

        int int_id = ++body_counter; // Body internal id!
        int ext_id = gi.id;          // Body external id!

        if ( modelDimension == ECIF_2D ) {
          body = new Body2D(MESH_BODY, int_id, ext_id, 0, NULL);

        } else if ( modelDimension == ECIF_3D ) {
          body = new Body3D(MESH_BODY, int_id, ext_id, 0, NULL);
        }

        if (body != NULL) {
          model->addBody(body);
          body->setName(gi.entityName);
        }
      }

      // Ok, home found/created for the bulk elements!
      if ( body != NULL ) {
        gi.parentTag = body->Tag();
      }

    // Boundary elements --> add elements to a boundary
    // =================
    } else if (gi.isBoundary) {

      nofBoundaryElements += gi.nofElements;

      // Try to find an existing boundary
      // --------------------------------
      BodyElement* be = findExistingBoundary(&gi);

      // Create a new boundary
      // ---------------------
      if ( be == NULL ) {

        int int_id = ++bndr_counter; // Internal id!

        // Create a proper boundary (body element)
        if ( modelDimension == ECIF_2D ) {
          be = new BodyElement2D(int_id, NO_INDEX, NO_INDEX, gi.nofElements);

        } else if ( modelDimension == ECIF_3D ) {
          be = new BodyElement3D(int_id, NO_INDEX, NO_INDEX, gi.nofElements);
        }

        if (be != NULL) {
          model->addBodyElement(be);
          be->setName(gi.entityName);
        }
      }

      // Ok, home found/created for the boundary elements!
      if ( be != NULL ) {
        gi.parentTag = be->Tag();
      }

    }

    // Read group elements
    // ===================
    if ( ECIF_OK != readElements(gi, element_counter, max_ext_elem_id, count_only) ) {
      return ECIF_ERROR;
    }
  }

  return ECIF_OK;
}


// Read element group info
// Returns true if read was succesful
//
bool
InputFidap::readElementGroupInfo(FidapElementGroupInfo& gi)
{
  static char key_buffer[80];
  static char value_buffer[80];

  eat_white(infile);

  // First group info line
  readFileLine(infile, read_buffer);
  strstream strm1;
  strm1 << read_buffer;
  while( !strm1.eof() ) {

    ws(strm1);

    strm1.getline(key_buffer, 80, ':');

    if ( 0 == strncmp(key_buffer, "GROUP", 5) )
      strm1 >> gi.id;
    if ( 0 == strncmp(key_buffer, "ELEMENTS", 8) )
      strm1 >> gi.nofElements;
    if ( 0 == strncmp(key_buffer, "NODES", 5) )
      strm1 >> gi.nofNodes;
    if ( 0 == strncmp(key_buffer, "GEOMETRY", 8) )
      strm1 >> gi.fidap_geometry;
    if ( 0 == strncmp(key_buffer, "TYPE", 4) )
      strm1 >> gi.fidap_type;
  }

  // Second group info line
  readFileLine(infile, read_buffer);
  strstream strm2;
  strm2 << read_buffer;
  while( !strm2.eof() ) {
    ws(strm2);
    strm2.getline(key_buffer, 80, ':');
    if ( 0 == strncmp(key_buffer, "ENTITY NAME", 11) )
      strm2.getline(gi.entityName, 80 , '\n');
  }

  LibFront::trim(gi.entityName);

  return true;
}


Rc
InputFidap::readElements(FidapElementGroupInfo& gi,
                         int& element_counter, int& max_ext_id,
                         bool count_only)
{
  Rc rc;
  int ext_el_id, ext_nd_id;

  static int ext_node_ids[MAX_NOF_NODES];
  static int ext_node_ids2[MAX_NOF_NODES]; // For reordering!

  meshElementCode elem_code = model->convertElementType(gi.elementType);

  if ( elem_code == MEC_000 ) {
    return ECIF_ERROR;
  }

  short nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

  short i, j;

  eat_white(infile);

  int group_elem_counter = 0;

	while ( group_elem_counter < gi.nofElements && !infile.eof() ) {

    if ( !count_only ) {
      theControlCenter->update(element_counter, GUI_UPDATE_INTERVAL);
    }

    int node_count = 0;

    while ( node_count < nof_nodes ) {

      // Element id
      if ( node_count == 0 ) {
        infile >> ext_el_id;
      }

      infile >> ext_nd_id;

      ext_node_ids[node_count] = ext_node_ids2[node_count] = ext_nd_id;
      node_count++;

    }

    if ( max_ext_id < ext_el_id ) {
      max_ext_id = ext_el_id;
    }

    group_elem_counter++;
    element_counter++;

    if ( count_only ) continue;

    // Apply possible preset reordering vector
    if ( gi.fidap_reorder != NULL ) {

      for (i = 0; i < nof_nodes; i++) {
        ext_node_ids[i] = ext_node_ids2[ gi.fidap_reorder[i] ];
      }

    // Or apply possible generic reorder for non-ok element types
    } else {
      switch (gi.elementType) {

      // These are ok
      case 303: case 304:
      case 404: case 405:
      case 504: case 505:
      case 706: case 707:
        break;

      // Others have a generic reordering
      default:
        for (i = 0, j = 0; j < nof_nodes - 1; i++, j += 2) {
          ext_node_ids[i] = ext_node_ids2[j];
        }

        for (j = 1; j < nof_nodes; i++, j += 2) {
          ext_node_ids[i] = ext_node_ids2[j];
        }
      break;
      }
    }

    elementCodeCounters[elem_code]++;

    // Add element to the model
    rc = model->addMeshInputElement(gi.elementType, ext_el_id, gi.parentTag, gi.id, ext_node_ids);

    if (gi.isBulk) {
      nofInputBulkElements++;
    } else if (gi.isBoundary) {
      nofInputBoundaryElements++;
    } else if (gi.isEdge) {
      nofInputEdgeElements++;
    } else {
      nofInputVertexElements++;
    }

    if (rc != ECIF_OK) {
      return rc;
    }

  } // while not eof()

	return ECIF_OK;
}


bool
InputFidap::readNumericKeywordValue(char* line_buffer, char* keyword, int& value)
{
  const int buffer_len = 80;
  static char buffer[buffer_len];

  if ( !readKeywordValue(line_buffer, keyword, buffer, buffer_len) )
    return false;

  value = atol(buffer);

  return true;
}


bool
InputFidap::readStringKeywordValue(char* line_buffer, char* keyword, char* value_buffer)
{
  const int buffer_len = 80;
  static char buffer[buffer_len];

  if ( !readKeywordValue(line_buffer, keyword, buffer, buffer_len) )
    return false;

  strcpy(buffer, value_buffer);

  return true;
}


bool
InputFidap::readKeyword(char* line_buffer, char* buffer, int buffer_len)
{
  strstream strm;
  strm << line_buffer << ends;

  if ( strm.eof() )
    return false;

  strm.getline(buffer, buffer_len, ':');

  // Trim trailing blanks
  int len = strlen(buffer);
  int pos = len - 1;

  while ( pos >= 0 && buffer[pos] == ' ' )
    pos--;

  buffer[1 + pos] = '\0';

  return true;
}


bool
InputFidap::readKeywordValue(char* line_buffer, char* keyword, char* buffer, int buffer_len)
{
  int strm_buf_len = strlen(line_buffer);
  char* strm_buf = new char[1 + strm_buf_len];

  strstream strm(strm_buf, 1 + strm_buf_len, ios::app);
  strm << line_buffer << ends;

  // Work buffer
  char* work_buf = new char[1 + strm_buf_len];
  char* tmp = work_buf;

  bool found = false;

  // Try to read keyword value
  while ( !strm.eof() ) {

    // Read next keyword
    strm.getline(tmp, strm_buf_len, ':');

    // Skip leading blanks
    LibFront::trimLeft(tmp);

     // If this was searched keyword
    if ( LibFront::ncEqualPartial(tmp, keyword) ) {
      found = true;

      // Jump over the searched keyword
      tmp += strlen(keyword);

      // Make a new stream from what is left
      reset(strm);
      strm << tmp << ends;
      break;
    }
  }

  if (!found)
    return false;

  // Pick until next keyword to get the keyword value
  // for the previous keyword!
  strm.getline(tmp, strm_buf_len, ':');

  // Skip backword the next keyword (until first blank)
  int len = strlen(tmp);
  while (len > 0 && tmp[len - 1] != ' ' )
    len--;
  tmp[len] = '\0';

  // Trim blanks
  LibFront::trim(tmp);

  // Check that value fits to the final buffer
  if ( strlen(tmp) <= buffer_len ) {
    strcpy(buffer, tmp);
    found = true;
  }
  else {
    found = false;
  }

  delete[] strm_buf;
  delete[] work_buf;

  return found;
}


// NOTE: Model dimension is stored in this stage in the class
// attribute 'dimension' which is read from the mesh file header
//
Rc
InputFidap::readMeshData(int& element_counter, int& node_counter)
{
  enum fidapKeywordId keyword = FIDAP_NO_ID;

  bool count_only = false;

  UserInterface* gui = theControlCenter->getGui();

  int body_counter = 0;
  int bndr_counter = 0;

  int max_ext_el_id = -1;
  int max_ext_nd_id = -1;

  while ( !infile.eof() ) {

    // Read new line
    readFileLine(infile, read_buffer);

    // Check if it was an interesting keyword line
    keyword = next_is_keyword_line(read_buffer);

    if ( keyword == FIDAP_NO_ID)
      continue;

    // Read elements (group)
    // =====================
    if ( keyword == FIDAP_ELEMENT_GROUPS_ID ) {
      if ( ECIF_OK != readElementGroup(body_counter, bndr_counter, element_counter, max_ext_el_id, count_only) ) {
        return ECIF_ERROR;
      }
    }

    // Read nodes
    // ==========
    else if ( keyword == FIDAP_NODAL_COORDINATES_ID ) {
      if ( ECIF_OK != readNodes(nofNodes, node_counter, max_ext_nd_id, count_only) ) {
        return ECIF_ERROR;
      }
      model->setMeshNodes();
    }

  } // while !eof

	return ECIF_OK;
}



// Create mesh geometry from Fidap neutral input file
//
// --Data should be in the following order in the input file:
// Nodes, Elements, Element groups
//
// --Those bulk and boundary elements which are not listed in groups, are processed later
// and the bodies/boundaries based on them are then created
//
// --Model dimension is known when nodes are read
//
bool
InputFidap::readMeshGeometry()
{
  int dimension;

  //---Header: Read dimension, nof nodes and elements
  //
  if (!readMeshHeader(dimension)) {
    modelDimension = ECIF_ND;
    infile.close();
    return false;
  }

  if ( nofElements == 0 || nofNodes == 0 ) {
    modelDimension = ECIF_ND;
    return false;
  }

  if ( dimension == 1 ) {
    dataDimension = ECIF_1D;
  } else if ( dimension == 2 ) {
    dataDimension = ECIF_2D;
  } else if (dimension == 3 ) {
    dataDimension = ECIF_3D;
  } else {
    dataDimension = ECIF_ND;
    return false;
  }

  // Check input dimension, it cannot be smaller than
  // data dimension!

  // 1D not actually supported!
  if ( inputDimension == ECIF_1D && dataDimension == ECIF_2D ) {
    inputDimension = ECIF_2D;
  }

  if ( inputDimension == ECIF_2D && dataDimension == ECIF_3D ) {
    inputDimension = ECIF_3D;
  }

  modelDimension = inputDimension;

  model->setModelDimension(modelDimension);

  findMaxExternalIds(maxExternalNodeId, maxExternalElementId);

  //---Allocate mesh tables

  // NOTE: This is total of bulk+bndr elements given explicitely in the file
  model->allocateMeshInputElements(nofElements, maxExternalElementId);

  model->allocateMeshNodes(nofNodes, maxExternalNodeId);

  //---Read nodes and element groups
  Rc rc;
  infile.clear();
  infile.seekg(0);
  int element_counter = 0;
  int node_counter = 0;
  rc = readMeshData(element_counter, node_counter);

  if (rc != ECIF_OK) {
    modelDimension = ECIF_ND;
    return false;
  }


  model->setModelDimension(modelDimension);

  if (modelDimension == ECIF_ND) {
    return false;
  }

  return true;
}


// All Fidap neutral files have similar header section
//
// NOTE: In header Nof-elements is sum of listed bulk and bndr elements!
//
bool
InputFidap::readMeshHeader(int& dimension)
{
  readFileLine(infile, read_buffer);

  if ( 0 != strncmp("** FIDAP NEUTRAL FILE", read_buffer, 21) ) {
    return false;
  }

  readFileLine(infile, read_buffer, 5);

  strstream strm;
  strm << read_buffer;

  if ( !(strm >> nofNodes >> nofElements >> nofElementGroups >> dimension) ) {
    return false;
  }

  return true;
}


Rc
InputFidap::readNodes(int nof_nodes,
                      int& node_counter, int& max_ext_id,
                      bool count_only)
{
  Rc rc;

  int ext_nd_id;
  static Point3 point;

  // Init point buffer
  point[0] = point[1] = point[2] = 0.0;

  bool is_3D = false;

  if ( dataDimension == ECIF_3D ) {
    is_3D = true;
  }

  // Read nof-nodes entries
  while ( !infile.eof() && node_counter < nof_nodes ) {

    theControlCenter->update(node_counter, GUI_UPDATE_INTERVAL);

    infile >> ext_nd_id;

    infile >> point[0] >> point[1];

    if (is_3D) {
      infile >> point[2];
    }

    node_counter++;

    if ( max_ext_id < ext_nd_id ) {
      max_ext_id = ext_nd_id;
    }

    // Update only counter etc.
    //
    if ( count_only ) {
      meshBoundingBox.extendByPoint(point);

    // Add node to the model
    //
    } else {
      rc = model->addMeshNode(node_counter - 1, ext_nd_id, point);
    }


  } // while not eof()

  // Model dimension can now be set!
  //
  if ( !count_only ) {
    modelDimension = findMeshModelDimension();
  }

	return ECIF_OK;
}


// Find ELMER element code for an Fidad element
bool
InputFidap::setElementType(FidapElementGroupInfo& gi)
{
  int base = 200;
  gi.fidap_reorder = NULL;

  switch (gi.fidap_geometry) {
  case 0:
    if ( dataDimension == ECIF_3D ) base = 400; break;
  case 1: base = 400; break;
  case 2: base = 300; break;
  case 3: base = 800; break;
  case 4: base = 700; break;
  case 5: base = 500; break;
  default:
    return false;
  }

  int type = base + gi.nofNodes;

  // Types needing specific reordering or mapping to
  // "lower" level types

  if (type == 203) {
    gi.fidap_reorder = fidapReorderTable[fidapReorderTableIndex_203];

  } else if (type == 510) {
    gi.fidap_reorder = fidapReorderTable[fidapReorderTableIndex_510];

  } else if (type >= 718 && type < 800) {
    if (type == 718)
      gi.fidap_reorder = fidapReorderTable[fidapReorderTableIndex_718];
    type = 706;

  } else if (type == 808) {
    gi.fidap_reorder = fidapReorderTable[fidapReorderTableIndex_808];

  } else if (type == 809) {
    gi.fidap_reorder = fidapReorderTable[fidapReorderTableIndex_809];

  } else if (type >= 820) {
    if (type == 827)
      gi.fidap_reorder = fidapReorderTable[fidapReorderTableIndex_827];
    type = 820;
  }

  gi.elementType = type;

  return true;

}
