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
Module:     ecif_inputIdeas.cpp
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
#include "ecif_bodyElement1D.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputIdeas.h"
#include "ecif_inputIdeasWF.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_timer.h"
#include "ecif_userinterface.h"

extern char read_buffer[];



// Indices to the special reorder table for element types
const int ideasReorderTableIndex_203 = 0;
const int ideasReorderTableIndex_306 = 1;
const int ideasReorderTableIndex_408 = 2;
const int ideasReorderTableIndex_510 = 3;
const int ideasReorderTableIndex_706 = 4;
const int ideasReorderTableIndex_808 = 5;
const int ideasReorderTableIndex_820 = 6;

// We have to reorder Ideas nodes to match Elmer numbering
const int ideasReorderTable[][MAX_NOF_NODES] = {

  {0,2,1},                                            // 0: Line 203
  {0,2,4,1,3,5},                                      // 1: Triangle 306
  {0,2,4,6,1,3,5,7},                                  // 2: Quadri 408
  {0,2,4,9,1,3,5,6,7,8},                              // 3: Tetra 510
  {0,2,1,3,5,4},                                      // 4: Wedge 706 2. version
  //{0,3,2,1,4,7,6,5},                                  // 5: Brick 808 1. version
  {0,1,2,3,4,5,6,7},                                  // 5: Brick 808 2. version
  {0,2,4,6,12,14,16,18,1,3,5,7,8,9,10,11,13,15,17,19} // 6: Brick 820
};



InputIdeas::InputIdeas(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
  elementCodeCounters = new int[1 + MAX_IDEAS_ELEMENT_CODE];

  for (short i = 0; i < 1 + MAX_IDEAS_ELEMENT_CODE; i++) {
    elementCodeCounters[i] = 0;
  }

  elementInfos = NULL;

  maxExternalBulkElemId = -1;
  maxExternalBndrElemId = -1;
  maxMaterialId = -1;
}


InputIdeas::~InputIdeas()
{
  delete[] elementCodeCounters;
  delete[] elementInfos;
}


// Method checks if end-of-dataset symbol is in the input-line.
// Ideas: '-1' in columns 5-6.
// Input argument is single line from a file
// Returns true if eod-symbol found.
bool
InputIdeas::endofDataset(char* s)
{
  if ( strlen(s) >= 6 && s[4] == '-' && s[5] == '1')
    return true;
  else
    return false;
}

// Method checks if end-of-dataset symbol is in the infile-file.
// Ideas: '-1' in columns 5-6.
// Input argument is the whole file.
// Returns true if eod-symbol found.
bool
InputIdeas::endofDataset()
{
  while (! infile.eof()) {
    readFileLine(infile, read_buffer);
    if ( endofDataset(read_buffer) )
      return true;
  }
  return false;
}


//***Find the next Object Header.
//  Returns the number for the dataset.
// Returns 0 if end-of-file reached.
int
InputIdeas::findDataset(int new_flag)
{
  //If next -1 in the inputfile means start of a new dataset (newflg=1)
  // or if it means end of a dataset (newflg = 0).
  int newflg = new_flag;
  int result;
  while (1) {
    if (infile.eof())
      return 0;
    // All interesting dataset must be here!
    switch (result = readDatasetNbr(newflg)) {
    case NEW_DATASET:
      newflg = 1;
      break;
    case DS_UNITS:
    case DS_WIRE_FRAME_CURVES:
    case DS_NODES_781:
    case DS_NODES_2411:
    case DS_ELEMENTS_780:
    case DS_ELEMENTS_2412:
    case DS_BODIES_AND_BOUNDARIES_2430:
    case DS_BODIES_AND_BOUNDARIES_2432:
    case DS_BODIES_AND_BOUNDARIES_2435:
      return result;
    default:
      newflg = 0;
      break;
    }
  }

}


enum ecif_modelDimension
InputIdeas::findCadModelDimension()
{
  // NOTE: Currently only 2D-cad models !!!
  return ECIF_2D;
}


enum ecif_modelDimension
InputIdeas::findMeshModelDimension() {

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
  dim = findModelDimension(box_dim);

  return dim;
}


void
InputIdeas::countNofBoundaryAndBulkElements() {


  for (short i = 0; i < 1 + MAX_IDEAS_ELEMENT_CODE; i++) {

    int nof_elements = elementCodeCounters[i];

    if (nof_elements == 0) continue;

    //--These are 2D bulk elements or 3D boundary elements
    if (i >= 1 && i <= 99) {

      //-2D bulk elements or 3D boundary model 'bulk' elements
      if ( modelDimension == ECIF_2D ) {
        nofBulkElements += nof_elements;

      //-3D boundary elements
      } else {
        nofBoundaryElements += nof_elements;
      }

    //--These are 3D bulk elements
    } else if ( i >= 100) {
      nofBulkElements += nof_elements;
    }
  }

  nofElements = nofBulkElements + nofBoundaryElements;
}



// Method reads all bodies in the input file.
bool
InputIdeas::readCadGeometry()
{
  InputIdeas* input = NULL;

  // Find out correct dimension for the geometry if its is unknown
  if (modelDimension == ECIF_ND) {
    modelDimension = findCadModelDimension();
  }

  // Read-method depends on the type of the input file and dimension.
  // If reading of input is unsuccesfule (modelType = ECIF_ND), we stop
  while (!infile.eof() && modelDimension != ECIF_ND) {
    int dataset;
    switch(dataset = findDataset()) {
    case DS_WIRE_FRAME_CURVES:
      input = new InputIdeasWF(modelDimension, infile, infileName);
      if (!input->readCadGeometry()) {
        modelDimension = ECIF_ND;
      }
      break;
    default:
      break;
    }
  }

  // After reading all body-information, data can be checked!
  //if (!model->checkBodies())
  //  modelDimension = ECIF_ND;

  delete input;

  return (modelDimension != ECIF_ND);
}


/*
// Method reads color-code information from the input file.
int
InputIdeas::readColor(ifstream& infile)
{
  int color_id;
  RGBfloat* color = new RGBfloat[4];

  // Color number (0,1,2,...)
  readFileLine(infile, read_buffer);
  istrstream  strline1(read_buffer);
  strline1 >> color_id;

  // Color's RGB value,
  readFileLine(infile, read_buffer);
  istrstream  strline2(read_buffer);
  strline2 >> color[0] >> color[1] >> color[2];
  color[3] = 1.0f;

  // Color is added to color-table.
  (*colorTable)[color_id] = color;

  // End-of-dataset code must be found!
  if (endofDataset(infile, read_buffer))
    return 1;
  else
      return 0;
}
*/


// Method reads dataset number from Ideas univ. file.
// Returns the number of the dataset or 0 if nothing was found.
// Datasets are marked by characters '-1' in the columns 5 and 6.
// Parameter *newflag* tells if next -1 means start of a new dataset
// (newflag = 1) or end of the previous (newflag = 0).
int
InputIdeas::readDatasetNbr(int newflag)
{

  bool newset = endofDataset();
  // If we must pass over to the beginning of the next dataset
  if (newset && !newflag)
    newset = endofDataset();

  // Start of dataset was found.
  if (newset) {
    readFileLine(infile, read_buffer);
    istrstream strline(read_buffer);
      int nbr;
      strline >> nbr;
    return nbr;
  }
  // There was no new dataset
  else
    return 0;
}


// Methods checks the type of body-element in the input file.
// To be defined in derived classes
ecif_geometryType
InputIdeas::readGeomType(char* s)
{
  return ECIF_NODIM;
}


// Method reads header information from the input file.
// (this is model level information (units etc.))
bool
InputIdeas::readCadHeader()
{
  int i;

  // Read units dataset
  if (findDataset() != DS_UNITS) {
    cerr << ": cannot find units" << endl;
    return false;
  }

	readFileLine(infile, read_buffer);
	istrstream strline1(read_buffer);

  if (!(strline1 >> model->unit.code))
    return false;

  // Read unit string
  char *cptr;
  char c;
  cptr = model->unit.str; i = 0;

  while (((c = strline1.get())) && (i < MAXBODYNAME)) {
    *cptr = c;
    cptr++; i++;
  }
  *cptr = '\0';

  if (!(strline1 >> model->unit.descr))
    model->unit.descr = 0;

  // Read coordinate scaling factor
  // Note: we want to use factors like x = factor_x * x
	readFileLine(infile, read_buffer);
	istrstream strline2(read_buffer);
  eat_white(strline2);
  strline2.getline(read_buffer, 80, ' ');

  double scale = (double) atof(read_buffer);

  // Set same scaling factor for x,y,z
  for (i = 0; i < 3; i++) {
    if ( scale != 0 )
      model->unit.conv[i] = 1 / scale;
  }

  // Scaling factor for weight component
  model->unit.conv[3] = 1.0;

  // End-of-dataset code must be found!
  if (endofDataset())
    return true;
  else
    return false;
}


// Method reads line body-element (line,plane).
// To be defined in derived classes
bool
InputIdeas::readLine(Body* body, char* buffer)
{
  return false;
}

 // Method reads a nurbs-spline formed body-element.
// To be defined in derived classes
bool
InputIdeas::readNurbs(Body* body, char* buffer)
{
  return false;
}


// Method reads one point from input-string *s*.
GcPoint*
InputIdeas::readPoint(char* s)
{
  double p[3];
  istrstream strline(s);

  for (int i = 0; i < 3; i++) {

    eat_white(strline);
    strline.getline(read_buffer, 80, ' ');

    p[i] = (double) atof(read_buffer);

    p[i] *= model->unit.conv[i];
  }

  GcPoint* point = new GcPoint(p);

  return point;
}


// Method reads one point from input-string *s*.
bool
InputIdeas::readPoint(char* s, Point3& p)
{
  p[0] = p[1] = p[2] = 0.0;

  for (int i = 0; i < 3; i++) {

    while ( *s != '\0' && *s == ' ' )
      s++;

    p[i] = (double) atof(s);

    p[i] *= model->unit.conv[i];

    while ( *s != '\0' && *s != ' ' )
      s++;
  }

  return true;
}


// Method reads one vertex from input-string *s*.
BodyElement*
InputIdeas::readVertex(char* s)
{
  static Point3 point;
  static GcPoint gpoint;

  // Create vertex from the point in the input file.
  if ( !readPoint(s, point) )
    return NULL;

  gpoint.setPosit(point[0], point[1], point[2]);

  BodyElement* vertex = new BodyElement1D(&gpoint);

  return vertex;
}


//*********************************************

// NOTE: All interesting dataset must be added also into
// method: InputIdeas::findDataset !!!***!!!
bool
InputIdeas::readMeshGeometry()
{
  // For the very first, check that file makes sense!
  if (!readMeshHeader()) {
    modelDimension = ECIF_ND;
    infile.close();
    return false;
  }

  Rc rc;
  int data_set, new_flag;

  // -Calc nof elements and nodes
  // -Read body and boundary sets
  // ----------------------------
  int nof_nodes = 0;
  int nof_elements = 0;
  bool nof_elements_read = false;
  bool nof_nodes_read = false;
  bool has_element_sets = false;

  new_flag = 1;

  while ( !doBreak() && !infile.eof() ) {

    data_set = findDataset(new_flag);
	  switch (data_set) {

	  case DS_NODES_781:
	  case DS_NODES_2411:
	  	rc =   readNofMeshNodes(data_set, nof_nodes, maxExternalNodeId);
      nof_nodes_read = true;
      new_flag = 1;
		  break;

	  case DS_ELEMENTS_780:
	  case DS_ELEMENTS_2412:
		  rc = readNofMeshElements(data_set, nof_elements, maxExternalElementId, maxMaterialId);
      nof_elements_read = true;
      new_flag = 1;
      modelDimension = findMeshModelDimension();
      countNofBoundaryAndBulkElements();
		  break;

	  case DS_BODIES_AND_BOUNDARIES_2430:
	  case DS_BODIES_AND_BOUNDARIES_2432:
	  case DS_BODIES_AND_BOUNDARIES_2435:
		  has_element_sets = true;
      new_flag = 1;
		  break;

    default:
      new_flag = 0;
		  break;
    }
  }

  if ( doBreak() || rc != ECIF_OK || !nof_nodes_read || !nof_elements_read ) {
    modelDimension = ECIF_ND;
    return false;
  }

  model->setModelDimension(modelDimension);

  // Allocate input tables
  // ---------------------
  model->allocateMeshNodes(nof_nodes, maxExternalNodeId);
  model->allocateMeshInputElements(nofElements, maxExternalElementId);

  //model->allocateMeshBulkElements(nofBulkElements, maxExternalBulkElemId);
  //model->allocateMeshBulkElements(nofBulkElements, maxExternalElementId);

  // NOTE: We store boundary elements first in a temporary table, because we want to
  // store also the material id for these elements, this is for the case they would
  // be BEM boundaries and need finally some body id!
  //
  //model->allocateMeshInputElements(nofBoundaryElements, maxExternalBndrElemId);
  //model->allocateMeshInputElements(nofBoundaryElements, maxExternalElementId);


  // Rewind to the beg of the file
  infile.clear();
  infile.seekg(0);

  bool elements_read = false;
  bool nodes_read = false;
  new_flag = 1;

  // Read mesh nodes and elements
  // ----------------------------
  while ( !doBreak() &&!(elements_read && nodes_read) && !infile.eof() ) {

    int data_set = findDataset(new_flag);

	  switch (data_set) {

	  case DS_NODES_781:
	  case DS_NODES_2411:
	  	rc = readMeshNodes(data_set);
      model->setMeshNodes();
      nodes_read = true;
      new_flag = 1;
		  break;

	  case DS_ELEMENTS_780:
	  case DS_ELEMENTS_2412:
		  rc = readMeshElements(data_set);
      elements_read = true;
      new_flag = 1;
		  break;

    default:
      new_flag = 0;
		  break;
    }
  }

  if (doBreak() || rc != ECIF_OK) {
    modelDimension = ECIF_ND;
	  model->setModelDimension(ECIF_ND);
    return false;
  }

  if ( !has_element_sets ) return true;

  // Rewind to the beg of the file
  infile.clear();
  infile.seekg(0);
  new_flag = 1;

  // Read element sets
  // -----------------
  while ( !doBreak() && !infile.eof() ) {

    data_set = findDataset(new_flag);
	  switch (data_set) {

	  case DS_BODIES_AND_BOUNDARIES_2430:
	  case DS_BODIES_AND_BOUNDARIES_2432:
	  case DS_BODIES_AND_BOUNDARIES_2435:
		  rc = readMeshBodiesAndBoundaries(data_set);
      new_flag = 1;
		  break;

    default:
      new_flag = 0;
		  break;
    }
  }

  if ( rc != ECIF_OK) {
    modelDimension = ECIF_ND;
	  model->setModelDimension(ECIF_ND);
    return false;
  }

  return true;
}


// All Universal files have similar header section
bool
InputIdeas::readMeshHeader()
{
  return readCadHeader();
}


// Read bodies and boundaries given separetely
// like in dataset 2430, 2435
//
// Create corresponding bodies and boundaries!
//
// NOTE: Element set must be after elements in the file, otherwise
// we do not know the max element ids when reading the set and this
// is necessay in order to separte bulk/boundary elements!!!
//
Rc
InputIdeas::readMeshBodiesAndBoundaries(int data_set)
{
	Rc rc;

  char name_buffer[256];

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  int elems_per_row;

  // How many elements per row given in file
  if ( data_set == DS_BODIES_AND_BOUNDARIES_2435 ) {
    elems_per_row = 2;
  } else {
    elems_per_row = 4;
  }

  int set_id;
  int element_id, type_id;
  int nof_elements;
  int dummy;

  int material_id;

  int body_counter = 0;  // Internal body tag
  int bndr_counter = 0;  // Internal boundary tag

  int bulk_elements_counter = 0;
  int bndr_elements_counter = 0;

  // Read first line of the first element set
  readFileLine(infile, read_buffer);

	// Read all (element-id,element-value) pairs.
	while ( !endofDataset(read_buffer) ) {

    reset(strm);
		strm << read_buffer << ends;

    //-Read set id
    //
		strm >> set_id;

    //-Skip six dummy numbers (0)
    for (short i = 0; i < 6; i++)
		  strm >> dummy;
    //-Read nof elements in the set
    strm >> nof_elements;

    //-Read set name
	  readFileLine(infile, read_buffer);
    strcpy(name_buffer, read_buffer);

    bool bulk_elements = false;
    bool bndr_elements = false;

    //-Read element ids in the set
    //
    int counter = 0;
    while (counter < nof_elements) {

      //-Read one line of the set
	    readFileLine(infile, read_buffer);
      reset(strm);
      strm << read_buffer << ends;

      // Read 4 pairs of (8, elementid) for DS_BODIES_AND_BOUNDARIES_2430
      // Read 4 pairs of (8, elementid) for DS_BODIES_AND_BOUNDARIES_2432
      // Read 2 pairs of (8, elementId, dummy, dummy) for DS_BODIES_AND_BOUNDARIES_2435
      for (int i = 0; i < elems_per_row; i++) {
        strm >> type_id; // Nbr 8
        strm >> element_id;

        if ( data_set == DS_BODIES_AND_BOUNDARIES_2435 ) {
          strm >> dummy >> dummy;
        }

        int inp_elem_id = model->getMeshInputElementIdExt2Int(element_id);

        if (counter == 0) {
          create_dyna_string(elementSetInfos[set_id].name, name_buffer);
          elementSetInfos[set_id].nofElements = nof_elements;
        }

        // Bulk element set
        // ----------------
        if ( bulk_elements || model->isMeshInputBulkElement(element_id) ) {

          // If set starts
          if (bulk_elements == false) {
            bulk_elements_counter += nof_elements;
            bulk_elements = true;
            body_counter++;
            material_id = body_counter + maxMaterialId;
          }

          // If this is the first element in the set
          if (counter == 0) {

            Body* body = NULL;
            if (modelDimension == ECIF_2D)
              body = new Body2D(MESH_BODY, body_counter, material_id, nof_elements, NULL);
            else if (modelDimension == ECIF_3D)
              body = new Body3D(MESH_BODY, body_counter, material_id, nof_elements, NULL);

            if (body != NULL) {
              model->addBody(body);
              body->setName(name_buffer);
              //model->setMeshBodyExt2IntFlag(material_id);
            }
          }

          model->setMeshInputElementParentTag(inp_elem_id, body_counter);
          model->setMeshInputElementExtParentTag(inp_elem_id, material_id);
        }

        // Boundary element set
        // --------------------
        else {

          // If set starts
          if (bndr_elements == false) {
            bndr_elements_counter += nof_elements;
            bndr_elements = true;
            bndr_counter++;
          }

          // If this is the first element in the set
          if (counter == 0) {

            BodyElement* be = NULL;
            if (modelDimension == ECIF_2D)
              be = new BodyElement2D(bndr_counter, NO_INDEX, NO_INDEX, nof_elements);
            else if (modelDimension == ECIF_3D)
              be = new BodyElement3D(bndr_counter, NO_INDEX, NO_INDEX, nof_elements);

            if (be != NULL) {
              model->addBodyElement(be);
              be->setName(name_buffer);
            }
          }

          model->setMeshInputElementParentTag(inp_elem_id, bndr_counter);
        }

        counter++;

      } // end for

    } // end while set elements

    // Read first line of the next element set
    readFileLine(infile, read_buffer);

  }

  delete[] strm_buffer;

	return ECIF_OK;
}


// NOTE: Ideas element definition consist of two lines
Rc
InputIdeas::readMeshElements(int data_set)
{
	Rc rc;
	int ext_elem_id, ext_node_id;
	int elem_type, material_id, nof_nodes;
  meshElementCode elem_code;
	int dummy;

  static int ext_node_ids[MAX_NOF_NODES];
  static int ext_node_ids2[MAX_NOF_NODES]; // For reordering!

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  // Read first line of the first element
	readFileLine(infile, read_buffer);

	// Read all (element-id,element-value) pairs.
	while ( !endofDataset(read_buffer) ) {

		// First line per element
    reset(strm);
		strm << read_buffer << ends;

    //--Element id, element type
		strm >> ext_elem_id >> elem_type;

	  //--Material id
		if (data_set == DS_ELEMENTS_780) {
		  strm >> dummy >> dummy >> dummy >> material_id >> dummy;
    }
		else if (data_set == DS_ELEMENTS_2412) {
		  strm >> dummy >> material_id >> dummy;
    }

    //--Nof nodes per element
    strm >> nof_nodes;

    bool reorder_nodes = false;
    int reorder_table_index = NO_INDEX;
    meshElementCode reverse_table_index = MEC_000;

	  //---Recode element type
    switch (elem_type) {

    case 21:  elem_type = 202;   // Linear beam
      break;

    case 24:  elem_type = 203;   // Parabolic beam
      reorder_table_index = ideasReorderTableIndex_203;
      break;

    case 41:  elem_type = 303;   // Linear triangle (plane stress)
    case 51:  elem_type = 303;   // Linear triangle (plane strain)
    case 81:  elem_type = 303;   // Linear triangle (axi solid)
    case 91:  elem_type = 303;   // Linear triangle (thin sell)
      reverse_table_index = MEC_303;
      break;

    case 42:  elem_type = 306;   // Parabolic triangle (plane stress)
    case 52:  elem_type = 306;   // Parabolic triangle (plane strain)
    case 82:  elem_type = 306;   // Parabolic triangle (axi solid)
    case 92:  elem_type = 306;   // Parabolic triangle (thin sell)
      reorder_table_index = ideasReorderTableIndex_306;
      reverse_table_index = MEC_306;
      break;

    case 44:  elem_type = 404;   // Linear quadrilateral (plane stress)
    case 54:  elem_type = 404;   // Linear quadrilateral (plane strain)
    case 84:  elem_type = 404;   // Linear quadrilateral (axi)
    case 94:  elem_type = 404;   // Linear quadrilateral (thin sell)
      reverse_table_index = MEC_404;
      break;

    case 45:  elem_type = 408;   // Parabolic quadrilateral (plane stress)
    case 55:  elem_type = 408;   // Parabolic quadrilateral (plane strain)
    case 85:  elem_type = 408;   // Parabolic quadrilateral (axi)
    case 95:  elem_type = 408;   // Parabolic quadrilateral (thin sell)
      reorder_table_index = ideasReorderTableIndex_408;
      reverse_table_index = MEC_408;
      break;

    case 111: elem_type = 504;   // Linear tetrahedron
      break;

    case 118: elem_type = 510;   // Parabolic tetrahedron
      reorder_table_index = ideasReorderTableIndex_510;
      break;

    case 112: elem_type = 706;   // Linear wedge
      reorder_table_index = ideasReorderTableIndex_706;
      break;

    case 115: elem_type = 808;   // Linear brick (8 nodes)
      reorder_table_index = ideasReorderTableIndex_808;
      break;

    case 116: elem_type = 820;   // Parabolic brick (8 nodes)
      reorder_table_index = ideasReorderTableIndex_820;
      break;

    // Unsupported element type!
    default:
      return ECIF_ERROR;
    }

    elem_code = model->convertElementType(elem_type);

    // Check element type
    //
    bool is_bulk = true;
    if ( modelDimension == ECIF_2D ) {
      if ( elem_type < 300 ) is_bulk = false;
    } else if ( modelDimension == ECIF_3D ) {
      if ( elem_type < 500 ) is_bulk = false;
    }

    // Update external element counters
    // NOTE: These are needed when reading element sets!
    //
    if ( is_bulk ) {
      if ( ext_elem_id > maxExternalBulkElemId ) {
        maxExternalBulkElemId = ext_elem_id;
      }
    } else {
      if ( ext_elem_id > maxExternalBndrElemId ) {
        maxExternalBndrElemId = ext_elem_id;
      }
    }

		//--Element nodes
    int count = 0;
    while ( count < nof_nodes ) {

		  readFileLine(infile, read_buffer);

      reset(strm);
		  strm << read_buffer << ends;

      while ( count < nof_nodes && (strm >> ext_node_id) ) {
        //strm >> ext_node_id;
		    ext_node_ids[count] = ext_node_ids2[count] = ext_node_id;
        count++;
      }
    }

    //--Reorder nodes, if necessary (middle nodes to the end!)
    if ( reorder_table_index != NO_INDEX ) {

      for (int i = 0; i < nof_nodes; i++) {
        ext_node_ids[i] = ext_node_ids2[ ideasReorderTable[reorder_table_index][i] ];
      }
    }

    //-Store element data to the model

    //--Bulk element
    if ( is_bulk ) {
      rc =	model->addMeshInputElement(elem_type, ext_elem_id, NO_INDEX, material_id, ext_node_ids);
      nofInputBulkElements++;

    //--Boundary element
    } else {

      // Reverse ccw element into cw order as boundary element!
      //
      if ( false && reverse_table_index != MEC_000 ) {
        int* rindices = (int*)MeshElementReversedNodeIndices[reverse_table_index];

        int i;
        for (i = 0;  i < nof_nodes; i++) {
          ext_node_ids2[i] = ext_node_ids[rindices[i]];
        }

        for (i = 0;  i < nof_nodes; i++) {
          ext_node_ids[i] = ext_node_ids2[i];
        }
      }

      rc = model->addMeshInputElement(elem_type, ext_elem_id, NO_INDEX, material_id, ext_node_ids);
      nofInputBoundaryElements++;
    }

		if (rc != ECIF_OK) {
  		return rc;
	  	break;
		}

    // Read first line of the next element
		readFileLine(infile, read_buffer);
	}

  delete[] strm_buffer;

	return ECIF_OK;
}


Rc
InputIdeas::readMeshNodes(int data_set)
{
	Rc rc;
	int nd_id;
  static Point3 point;

	readFileLine(infile, read_buffer);

	// Read (node-id,node point) line-pairs.
  int internal_id = 0;
	while (! endofDataset(read_buffer)) {

		// First line per node
    nd_id = atol(read_buffer);

		// Second line per node
		readFileLine(infile, read_buffer);
    readPoint(read_buffer, point);

 		rc = model->addMeshNode(internal_id, nd_id, point);
    internal_id++;

    theControlCenter->update(internal_id, GUI_UPDATE_INTERVAL);

		if (rc != ECIF_OK)
			return rc;
		readFileLine(infile, read_buffer);
	}

	return ECIF_OK;
}


Rc
InputIdeas::readNofMeshElements(int data_set, int& nof_elements,
                                int& max_element_id, int& max_material_id)
{
  int el_id, material_id, dummy, nof_nodes;
  short el_code;

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  //NOTE: this works only if there are two lines per element!

  //---First line of the first element
  readFileLine(infile, read_buffer);
	while (! endofDataset(read_buffer)) {

    reset(strm);
    strm << read_buffer << ends;

    //-Read element id
	  strm >> el_id;
    if (el_id > max_element_id)
      max_element_id = el_id;

    //-Read element code
	  strm >> el_code;
    elementCodeCounters[el_code]++;

	  //--Material id and the others not currently interesting
		if (data_set == DS_ELEMENTS_780) {
		  strm >> dummy >> dummy >> dummy >> material_id >> dummy;
    }
		else if (data_set == DS_ELEMENTS_2412) {
		  strm >> dummy >> material_id >> dummy;
    }

    if ( material_id > max_material_id ) {
      max_material_id = material_id;
    }

    //--Nof nodes per element
    strm >> nof_nodes;

    //--Read element nodes
    int count = 0;
    while ( count < nof_nodes ) {

		  readFileLine(infile, read_buffer);

      reset(strm);
		  strm << read_buffer << ends;

      while ( count < nof_nodes && (strm >> dummy) ) {
        count++;
      }
    }

    nof_elements++;

    theControlCenter->update(nof_elements, GUI_UPDATE_INTERVAL);

    //---Read first line for the next element
		readFileLine(infile, read_buffer);
	}

  delete[] strm_buffer;

	return ECIF_OK;
}


Rc
InputIdeas::readNofMeshNodes(int data_set, int& nof_nodes, int& max_node_id)
{
  Point3 point;

  //NOTE: this works when there are two lines per node!
  // First line of the first node
  readFileLine(infile, read_buffer);
  int nd_id;

	while (! endofDataset(read_buffer)) {

	  nd_id = atol(read_buffer);
    if (nd_id > max_node_id)
      max_node_id = nd_id;

    // Second line of the node
		readFileLine(infile, read_buffer);
    readPoint(read_buffer, point);

    meshBoundingBox.extendByPoint(point);

    nof_nodes++;

    theControlCenter->update(nof_nodes, GUI_UPDATE_INTERVAL);

    // Read first line for the next node
		readFileLine(infile, read_buffer);
	}

  return ECIF_OK;
}
