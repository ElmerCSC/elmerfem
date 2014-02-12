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
Module:     ecif_inputAbaqus.cpp
Language:   C++
Date:       01.10.98
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
#include "ecif_inputAbaqus.h" 
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_timer.h"
#include "ecif_userinterface.h"

extern char read_buffer[];

InputAbaqus::InputAbaqus(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
  for (short i = 0; i < 1 + MAX_NOF_ELEM_CODES; i++) {
    elementCodeCounters[i] = 0;
  }

  nofElementSets = 0;
  nofNodeSets = 0;
}


InputAbaqus::~InputAbaqus()
{
}


AbaqusElementSet*
InputAbaqus::addElementSet(char* elem_set_name)
{
  if (nofElementSets == MAX_NOF_ABAQUS_ELEMENT_SETS) {
    return NULL;
  }

  AbaqusElementSet* es = &(elementSets[nofElementSets]);

  update_dyna_string(es->name, elem_set_name);
  
  nofElementSets++;

  // Set ids (used for bodies etc.)
  es->extParentTag = nofElementSets;

  return es;
}


AbaqusElementSet*
InputAbaqus::findElementSet(char* elem_set_name)
{
  AbaqusElementSet* eset = NULL;

  // Check if this is a stored body name
  for (int i = 0; i < nofElementSets; i++) {
    eset = &(elementSets[i]);

    if ( LibFront::ncEqual(eset->name, elem_set_name) ) {
      return eset;
    }
  }

  return NULL;
}


enum ecif_modelDimension 
InputAbaqus::findMeshModelDimension()
{
  ecif_modelDimension dim = ECIF_ND;
  ecif_modelDimension box_dim = ECIF_ND;
  ecif_modelDimension elem_dim = ECIF_ND;
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

  // Find actual dimension
  //
  dim = findModelDimension(box_dim);

  return dim;
}


void 
InputAbaqus::countNofBoundaryAndBulkElements()
{
  for (short i = 0; i < MAX_NOF_ELEM_CODES; i++) {

    int nof_elements = elementCodeCounters[i];

    if (nof_elements == 0) continue;

    //--These are 2D bulk elements of 3D boundary elements
    if (i = MEC_000 && i < MEC_504) {

      //-2D bulk elements or 3D boundary model 'bulk' elements
      if ( modelDimension == ECIF_2D ) {
        nofBulkElements += nof_elements;

      //-3D boundary elements
      } else {
        nofBoundaryElements += nof_elements;
      }

    //--These are 3D bulk elements
    } else if ( i >= MEC_504 ) {
      nofBulkElements += nof_elements;
    }

  }
}



// Find ELMER element code for an Abaqus element name
//
// NOTE: Currently only solid Abaqus elements and not even all
// of them are supported, they really have elements!!!
//
int
InputAbaqus::findElementCode(char* abq_name)
{
  // NOTE: Keep the falses at the end of each if statements
  //       This way it is safer to insert new types

  //---BEAMS
  // Linear (2 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "DC1D2")   ||
       false
     ) {
    return 202;
  }

  // Parabolic (3 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "DC1D3")   ||
       false
     ) {
    return 203;
  }

  //---TRIANGLES
  // Linear (3 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "CPE3")    ||
       LibFront::ncEqualPartial(abq_name, "CPS3")    ||
       LibFront::ncEqualPartial(abq_name, "DC2D3")   ||
       LibFront::ncEqualPartial(abq_name, "CAX3")    ||
       LibFront::ncEqualPartial(abq_name, "CGAX3")   ||
       LibFront::ncEqualPartial(abq_name, "DCAX3")   ||
       LibFront::ncEqualPartial(abq_name, "S3")      ||
       false
     ) {
    return 303;
  }

  // Parabolic (6 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "CPE6")    ||
       LibFront::ncEqualPartial(abq_name, "CPS6")    ||
       LibFront::ncEqualPartial(abq_name, "DC2D6")   ||
       LibFront::ncEqualPartial(abq_name, "CAX6")    ||
       LibFront::ncEqualPartial(abq_name, "CGAX6")   ||
       LibFront::ncEqualPartial(abq_name, "DCAX6")   ||
       LibFront::ncEqualPartial(abq_name, "S6")      ||
       false
     ) {
    return 306;
  }

  //---QUADRILATERALS
  // Linear (4 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "CPE4")    ||
       LibFront::ncEqualPartial(abq_name, "CPS4")    ||
       LibFront::ncEqualPartial(abq_name, "DC2D4")   ||
       LibFront::ncEqualPartial(abq_name, "AC2D4")   ||
       LibFront::ncEqualPartial(abq_name, "CAX4")    ||
       LibFront::ncEqualPartial(abq_name, "CGAX4")   ||
       LibFront::ncEqualPartial(abq_name, "DCAX4")   ||
       LibFront::ncEqualPartial(abq_name, "S4")      ||
       false
     ) {
    return 404;
  }

  // Parabolic (8 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "CPE8")    ||
       LibFront::ncEqualPartial(abq_name, "CPS8")    ||
       LibFront::ncEqualPartial(abq_name, "C2D8")    ||
       LibFront::ncEqualPartial(abq_name, "DC2D8")   ||
       LibFront::ncEqualPartial(abq_name, "AC2D8")   ||
       LibFront::ncEqualPartial(abq_name, "CAX8")    ||
       LibFront::ncEqualPartial(abq_name, "CGAX8")   ||
       LibFront::ncEqualPartial(abq_name, "DCAX8")   ||
       LibFront::ncEqualPartial(abq_name, "S8")      ||
       false
    ) {
    return 408;
  }

  //---TETRAS
  // Linear (4 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "C3D4")    ||
       LibFront::ncEqualPartial(abq_name, "DC3D4")   ||
       false
     ) {
    return 504;
  }

  // Parabolic (10 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "C3D10")    ||
       LibFront::ncEqualPartial(abq_name, "DC3D10")   ||
       false
     ) {
    return 510;
  }

  //---WEDGE
  // Linear (6 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "C3D6")    ||
       LibFront::ncEqualPartial(abq_name, "DC3D6")   ||
       false
     ) {
    return 706;
  }

  //---BRICKS
  // Linear (8 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "C3D8")    ||
       LibFront::ncEqualPartial(abq_name, "DC3D8")   ||
       false
     ) {
    return 808;
  }

  // Parabolic (20 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "C3D20")   ||
       LibFront::ncEqualPartial(abq_name, "DC3D20")  ||
       false
     ) {
    return 820;
  }

  // Parabolic (27 nodes)
  if ( LibFront::ncEqualPartial(abq_name, "C3D27")   ||
       LibFront::ncEqualPartial(abq_name, "DC3D27")  ||
       false
     ) {
    return 827;
  }

  //---Unknown or unsupported!!!
  return (int)NO_INDEX;
}


// Check if element (set) name is for an axisymmetric
// solid Abaqus element type
bool
InputAbaqus::isAxisymmetricElement(char* abq_name)
{
  if ( LibFront::ncEqualPartial(abq_name, "CAX")     ||
       LibFront::ncEqualPartial(abq_name, "CGAX")    ||
       LibFront::ncEqualPartial(abq_name, "DCAX")    
     )
    return true;

  return false;
}


bool
InputAbaqus::next_is_keyword_line()
{
  eat_white(infile);

  bool is_keyword_line = false;

  char c = infile.get();

  // If first is '*', check if next is not
  // '*' and we really have a keyword line
  if ( c == '*' && infile.peek() != '*' ) {
    is_keyword_line = true;
  }

  infile.putback(c);

  return is_keyword_line;

}



Rc
InputAbaqus::readElements(AbaqusElementSet* es, int elem_type,
                          bool count_only, int& element_counter)
{
  Rc rc;
  int ext_elem_id, ext_parent_tag;

  static int ext_node_ids[MAX_NOF_NODES];
  char number_buffer[81];

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  // Initially eof flag set true
  strm.clear(ios::eofbit);

  meshElementCode elem_code = model->convertElementType(elem_type);

  if ( elem_code == MEC_000 ) {
    return ECIF_ERROR;
  }

  short nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
  short nof_missing_nodes = 0;

  if ( es != NULL ) {
    ext_parent_tag = es->extParentTag;
  } else {
    ext_parent_tag = NO_INDEX;
  }

	while ( !infile.eof() ) {
 
    theControlCenter->update(element_counter, GUI_UPDATE_INTERVAL);

    // Read new line if needed
    if ( !infile.eof() && strm.eof() ) {

      skip_comment_lines();

      // If data is at end!!!
      if ( next_is_keyword_line() ) {
        break;
      }

      readFileLine(infile, read_buffer);

      // Change comma-separator to blanko
      for (int i = 0; i < strlen(read_buffer); i++) {

        if (read_buffer[i] == ',') {
          read_buffer[i] = ' ';
        }
      }

      reset(strm);
      strm << read_buffer << ends;
    }

    // Read possible missing node ids for the element
    while ( !strm.eof() && nof_missing_nodes > 0 ) {

      strm >> number_buffer;

      if (strlen(number_buffer) == 0) {
        strm.clear(ios::eofbit);
        break;
      }

      ext_node_ids[nof_nodes - nof_missing_nodes] = atol(number_buffer);

      nof_missing_nodes--;
    }

    // Read new element id if no node ids are missing for the
    // previous
    if ( !strm.eof() && nof_missing_nodes == 0 ) {

      strm >> number_buffer;

      // No more numbers in the stream
      if (strlen(number_buffer) == 0) {
        strm.clear(ios::eofbit);

      // Start new element
      } else {

        ext_elem_id = atol(number_buffer);
        element_counter++;
        elementCodeCounters[elem_code]++;

        if (es != NULL) {
          es->nofElements++;
        }

        nof_missing_nodes = nof_nodes;
      }
    }

    // Update max-external-id
    if ( maxExternalElementId < ext_elem_id ) {
      maxExternalElementId = ext_elem_id;
    }

    // Check element type
    bool is_bulk;
    bool is_bndr;
    bool is_edge;

    model->checkMeshElementType(elem_type, is_bulk, is_bndr, is_edge);
    
    // Add element to the model's input-element table if not in count-only mode
    //
    if ( !count_only && nof_missing_nodes == 0) {

      // NOTE: Parabolic beam is the only Abaqus element where
      // we have to reorder nodes!!!
      // (Why on the earth they did this exception!?!)
      if ( elem_type == 203 ) {
        int tmp = ext_node_ids[1]; // Pick middle node
        ext_node_ids[1] = ext_node_ids[2];
        ext_node_ids[2] = tmp; // Put middle node at the end
      }
     
      rc = model->addMeshInputElement(elem_type,
                                      ext_elem_id, NO_INDEX, ext_parent_tag,
                                      ext_node_ids);

      if (rc != ECIF_OK) return rc;

      if ( is_bulk ) {
        nofInputBulkElements++;
      } else if (is_bndr ) {
        nofInputBoundaryElements++;
      } else if (is_edge ) {
        nofInputEdgeElements++;
      } else {
        nofInputVertexElements++;
      }
    }
 
  } // while not eof()

  delete[] strm_buffer;

	return ECIF_OK;
}


bool
InputAbaqus::readElementSet(AbaqusElementSet* es, int& body_counter, int& bndr_counter, bool count_only)
{
  if (es == NULL ) return false;

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  // Initially eof flag set true
  strm.clear(ios::eofbit);

  char number_buffer[80];
  int parent_ids[2] = {NO_INDEX, NO_INDEX};
  
  // Only counting elements
  // ----------------------
  if ( count_only ) {
    es->nofElements = 0;

  // Create body/boundary
  // --------------------
  } else {
    
    int ext_elem_id = es->firstElementId;
    int int_elem_id = model->getMeshInputElementIdExt2Int(ext_elem_id);
    int elem_type = model->getMeshInputElementType(int_elem_id);

    bool is_bulk;
    bool is_bndr;
    bool is_edge;

    model->checkMeshElementType(elem_type, is_bulk, is_bndr, is_edge);
    
    // Create a new body
    // -----------------
    if ( is_bulk ) {

      Body* body = NULL;
 
      int int_id = ++body_counter; // Body internal id!
      int ext_id = es->extParentTag;   // Body external id!

      if ( modelDimension == ECIF_2D ) {
        body = new Body2D(MESH_BODY, int_id, ext_id, 0, NULL);

      } else if ( modelDimension == ECIF_3D ) {
        body = new Body3D(MESH_BODY, int_id, ext_id, 0, NULL);
      }

      if (body != NULL) {
        model->addBody(body);
        body->setName(es->name);
      }
    
      // Ok, home found/created for the bulk elements!
      if ( body != NULL ) {
        es->parentTag = body->Tag();
      }
      
    // Create a new boundary
    // ---------------------
    } else if ( is_bndr ) {

      BodyElement* be = NULL;
      
      int int_id = ++bndr_counter; // Internal id!

      // Create a proper boundary (body element)
      if ( modelDimension == ECIF_2D ) {
        be = new BodyElement2D(int_id, NO_INDEX, NO_INDEX, es->nofElements);

      } else if ( modelDimension == ECIF_3D ) {
        be = new BodyElement3D(int_id, NO_INDEX, NO_INDEX, es->nofElements);
      }

      if (be != NULL) {
        model->addBodyElement(be);
        be->setName(es->name);
      }

      // Ok, home found/created for the boundary elements!
      if ( be != NULL ) {
        es->parentTag = be->Tag();
      }
    }

  } // Create body/boundary


  // Read element ids
  // ----------------
  while ( !infile.eof() ) {

    skip_comment_lines();

    // If set data is at end, stop!
    //
    if ( next_is_keyword_line() ) break;

    theControlCenter->update(es->nofElements, GUI_UPDATE_INTERVAL);

    // Read new line if needed
    if ( !infile.eof() && strm.eof() ) {
      readFileLine(infile, read_buffer);

      // Change comma-separators to blankos
      for (int i = 0; i < strlen(read_buffer); i++) {
        if (read_buffer[i] == ',')
          read_buffer[i] = ' ';
      }

      reset(strm);
      strm << read_buffer << ends;
    }

    // Read element ids in the set
    //
    while ( !strm.eof() ) {
      strm >> number_buffer;

      if (strlen(number_buffer) == 0) {
        strm.clear(ios::eofbit);
        break;
      }

      int ext_elem_id = atol(number_buffer);

      if ( count_only ) {
        es->nofElements++;
        
        // Store first element id, so that we can find the
        // type of the boundary!
        //
        if ( es->nofElements == 1 ) {
          es->firstElementId = ext_elem_id;
        }

        continue;
      }

      // Ok, input element tables have been created, set parent tags
      // to the element
      //
      int int_elem_id = model->getMeshInputElementIdExt2Int(ext_elem_id);
      model->setMeshInputElementParentTag(int_elem_id, es->parentTag);
      model->setMeshInputElementExtParentTag(int_elem_id, es->extParentTag);

    }

  } // while not eof()

  delete[] strm_buffer;


  return true;
}


// Returns true if a new element set was encountered
// buffer: element set name
bool
InputAbaqus::readElementSetName(istrstream* strm, char* buffer, int buffer_len)
{
  if ( !readKeywordValue(strm, "ELSET=", buffer, buffer_len) ) {
    return false;
  }

  // Check if element set (name) is already stored
  AbaqusElementSet* es = findElementSet(buffer);

  // Add new set
  if ( es == NULL ) {
    es = addElementSet(buffer);
  }

  if (es == NULL)
    return false;
  else
    return true;
}
 


bool
InputAbaqus::readElementType(istrstream* strm, char* buffer, int buffer_len)
{
  return readKeywordValue(strm, "TYPE=", buffer, buffer_len);;
}


bool
InputAbaqus::readKeyword(istrstream* strm, char* buffer, int buffer_len)
{
  if ( strm->eof() )
    return false;

  strm->getline(buffer, buffer_len, ',');

  // Trim trailing blanks
  int len = strlen(buffer);
  int pos = len - 1;

  while ( pos >= 0 && buffer[pos] == ' ' )
    pos--;

  buffer[1 + pos] = '\0';

  return true;
}


bool
InputAbaqus::readKeywordValue(istrstream* istrm, char* keyword, char* buffer, int buffer_len)
{
  // Work stream
  char* tmp = ((strstream*)istrm)->str(); // A stupid cast for Sgi!
  int tmp_buf_len = strlen(tmp);
  char* tmp_buf = new char[1 + tmp_buf_len];
  strstream strm(tmp_buf, 1 + tmp_buf_len, ios::app);
  strm << tmp << ends;

  // Work buffer
  char* tmp_buf2 = new char[1 + tmp_buf_len];

  bool found = false;

  // Try to read keyword value
  while ( !strm.eof() ) {

    // Read next keyword
    strm.getline(tmp_buf2, tmp_buf_len, ',');

    // Skip leading blanks
    char* tmp = tmp_buf2;
    int len2 = strlen(tmp_buf2);
    int pos = 0;
    while ( pos < len2 && tmp[pos] == ' ')
      pos++;
    tmp += pos;

    // If this was searched eyword
    if ( LibFront::ncEqualPartial(tmp, keyword) ) {
      found = true;
      // Jump over the keyword
      tmp += strlen(keyword);

      // Make a new stream from what is left
      reset(strm);
      strm << tmp << ends;
      break;
    }
  }

  if (!found)
    return false;

  // Pick keyword value to tmp_buf
  strm.getline(tmp_buf, tmp_buf_len, ',');

  // Trim trailing blanks
  int len = strlen(tmp_buf);
  while (len > 0 && tmp_buf[len - 1] == ' ' )
    len--;
  tmp_buf[len] = '\0';

  // Check that value fits to the final buffer
  if ( strlen(tmp_buf) <= buffer_len ) {
    strcpy(buffer, tmp_buf);
    found = true;
  }
  else {
    found = false;
  }

  delete[] tmp_buf;
  delete[] tmp_buf2;

  return found;
}


// Reads nodes, elements etc.
Rc
InputAbaqus::readMeshData(int& element_counter, int& node_counter, bool count_only)
{
  Rc rc;

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  int keyword_buffer_len = 80;
  char* keyword_buffer = new char[keyword_buffer_len];

  int name_buffer_len = 80;
  char* name_buffer = new char[name_buffer_len];

  UserInterface* gui = theControlCenter->getGui();

  int body_counter = 0;
  int bndr_counter = 0;

  // Read file
  //
  while ( !infile.eof() ) {

    skip_comment_lines();

    bool keyword_found = false;

    if ( next_is_keyword_line() ) 
      keyword_found = true;

    // Read new line
    readFileLine(infile, read_buffer);

    if (!keyword_found)
      continue;

    char* tmp = read_buffer;
    tmp++;

    reset(strm);
    strm << tmp << ends;


    istrstream* istrm = (istrstream*) &strm;

    if ( !readKeyword(istrm, keyword_buffer, 80) )
      keyword_found = false;

    // Uups, something went wrong!!!
    if (!keyword_found)
      continue;

    // Elements
    // ========
    if ( LibFront::ncEqual(keyword_buffer, (char*)ABQ_ELEMENT) ) {

      // Read element type
      if ( !readElementType(istrm, name_buffer, name_buffer_len) ) {
        return ECIF_ERROR;
      }

      int elem_code = findElementCode(name_buffer);

      if (elem_code == (int)NO_INDEX) {

        strstream strm;
        strm << "***ERROR: Unsupported Abaqus element type: ";
        strm << name_buffer << ends;
        gui->showMsg(strm.str());

        return ECIF_ERROR;
      }

      AbaqusElementSet* es = NULL;

      // Read possible element set (=body) name

      // NOTE: Element sets are also given in separete sections where only
      // element ids are given
      //
      if ( readElementSetName(istrm, name_buffer, name_buffer_len) ) {
        es = findElementSet(name_buffer);
      }

      rc = readElements(es, elem_code, count_only, element_counter);

      if ( rc != ECIF_OK ) return rc;
    }

    // Nodes
    // =====
    else if ( LibFront::ncEqual(keyword_buffer, (char*)ABQ_NODE) ) {

      readNodeSetName(istrm, name_buffer, name_buffer_len);
      readNodes(count_only, node_counter);

      if (!count_only) {
        model->setMeshNodes();

      // We can now set the model dimension!
      } else {
        modelDimension = findMeshModelDimension();
      }
    }

    // Element set
    // ===========
    else if ( LibFront::ncEqual(keyword_buffer, (char*)ABQ_ELSET) ) {
      AbaqusElementSet* es = NULL;

      if ( readElementSetName(istrm, name_buffer, name_buffer_len) ) {
        es = findElementSet(name_buffer);
      }

      readElementSet(es, body_counter, bndr_counter, count_only);
    }

    // Node set (node set section is not actually used!)
    // ========
    else if ( LibFront::ncEqual(keyword_buffer, (char*)ABQ_NSET) ) {
      readNodeSetName(istrm, name_buffer, name_buffer_len);
      readNodeSet(name_buffer, count_only);
    }

	}

  delete[] keyword_buffer;
  delete[] name_buffer;
  delete[] strm_buffer;

	return ECIF_OK;
}


// Create mesh geometry from Abaqus input file
//
// --Data should be in the following order in the input file:
// Nodes, Elements, Element sets
//
// --File is read in two phases:
// 1. phase: calculate nof nodes, nof elements and nof elements in sets
// 2. phase: install nodes, elements and create bodies/boundaries by element set
//
// --Those bulk and boundary elements which are not listed in sets, are processed later
// and the bodies/boundaries based on them are then created
//
// --Model dimension is known when nodes are read
//
bool
InputAbaqus::readMeshGeometry()
{
  // For the very first, check that file makes sense!
  if (!readMeshHeader()) {
    modelDimension = ECIF_ND;
    infile.close();
    return false;
  }

  Rc rc;

  int element_counter;
  int node_counter;
  bool count_only;

  //--- 1. pass. Calc nof nodes and elements, create element sets,
  //
  element_counter = 0;
  node_counter = 0;
  count_only = true;

  rc = readMeshData(element_counter, node_counter, count_only);

  nofElements = element_counter;
  nofNodes = node_counter;

  if ( rc != ECIF_OK || nofElements == 0 || nofNodes == 0 ) {
    modelDimension = ECIF_ND;
    return false;
  }

  // Model dimension is known when nodes have been read!
  model->setModelDimension(modelDimension);

  //---Allocate temporary mesh elements table (for bulk+bndr explicitely in file)
  model->allocateMeshInputElements(nofElements, maxExternalElementId);

  model->allocateMeshNodes(nofNodes, maxExternalNodeId);

  //--- 2. pass. Install nodes and elements, create bodies/boundaries by sets
  //
  infile.clear();
  infile.seekg(0);
  element_counter = 0;
  node_counter = 0;
  count_only = false;

  rc = readMeshData(element_counter, node_counter, count_only);

  if (rc != ECIF_OK) {
    modelDimension = ECIF_ND;
    return false;
  }

  return true;
}


// Check Abaqus-type header
bool
InputAbaqus::readMeshHeader()
{
  return true;
}


Rc
InputAbaqus::readNodes(bool count_only, int& node_counter)
{
  Rc rc;

  int external_nd_id;
  static Point3 point;
  char number_buffer[81];

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);
  // Initially eof flag set true
  strm.clear(ios::eofbit);

  short nof_coord = 3;
  short nof_missing_coord = 0;
  
  // Init point buffer
  point[0] = point[1] = point[2] = 0.0;

  // Read nodes until a non-data line is encountered
	while ( !infile.eof() ) {

    theControlCenter->update(node_counter, GUI_UPDATE_INTERVAL);

    // Read new line if needed
    if ( !infile.eof() && strm.eof() ) {

      skip_comment_lines();
      // Check that data is not at end
      if ( next_is_keyword_line() )
        break;

      readFileLine(infile, read_buffer);

      // Change comma-separator to blankos
      for (int i = 0; i < strlen(read_buffer); i++) {
        if (read_buffer[i] == ',')
          read_buffer[i] = ' ';
      }

      reset(strm);
      strm << read_buffer << ends;
    }

    // Read possible missing coordinates for the node
    while ( !strm.eof() && nof_missing_coord > 0 ) {
      strm >> number_buffer;

      if (strlen(number_buffer) == 0) {
        strm.clear(ios::eofbit);
        break;
      }

      point[nof_coord - nof_missing_coord] = atof(number_buffer);
      nof_missing_coord--;
    }

    // Read possible node id if no coordinates are missing
    if ( !strm.eof() && nof_missing_coord == 0 ) {
      strm >> number_buffer;

      // No more numbers in the stream
      if (strlen(number_buffer) == 0) {
        strm.clear(ios::eofbit);

      // Start new node
      } else {
        external_nd_id = atol(number_buffer);
        node_counter++;

        if (external_nd_id > maxExternalNodeId)
          maxExternalNodeId = external_nd_id;

        nof_missing_coord = nof_coord;
      }

    }

    // Add node to the model if everything is ready
    if ( !count_only && nof_missing_coord == 0 )
      rc = model->addMeshNode(node_counter - 1, external_nd_id, point);
    
    if ( count_only ) {
      meshBoundingBox.extendByPoint(point);
    }

  } // while not eof()

  delete[] strm_buffer;

	return ECIF_OK;
}


bool
InputAbaqus::readNodeSet(char* node_set_name, bool count_only)
{
  return true;
}


bool
InputAbaqus::readNodeSetName(istrstream* strm, char* buffer, int buffer_len)
{
  return readKeywordValue(strm, "NSET=", buffer, buffer_len);;
}


void
InputAbaqus::skip_comment_lines()
{
  while ( !infile.eof() ) {

    eat_white(infile);

    char c = infile.get();

    // If first is '*', check if next is also
    // a '*' and we really have a comment line
    // and we skip it
    if ( c == '*' && infile.peek() == '*' ) {
      readFileLine(infile, read_buffer);
    }
    // Otherwise we put character back and stop
    else {
      infile.putback(c);
      break;
    }
  }

}

