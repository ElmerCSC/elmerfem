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
Module:     ecif_inputIges.cpp
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
#include "ecif_geometry.h"
#include "ecif_inputIges.h"
#include "ecif_model.h"

const int IGES_LINE_LEN = 80;
const int SEC_ID_POS = 72;
const int DIR_FLD_LEN = 8;
const int DATA_PART_LEN = 64;
char DATA_END = ';';
char DATA_SEP = ',';

// Global buffers
char fieldBuffer[1+DIR_FLD_LEN];
char lineBuffer[1+IGES_LINE_LEN];
char dataBuffer[1+DATA_PART_LEN];
char dataFieldBuffer[1+DATA_PART_LEN];

int  statusNbrBuffer[4];
char statusFldBuffer[3];


// Constructor
IgesStatusField::IgesStatusField()
{
  blankFlag = 0;
  subordFlag  = 0;
  useFlag   = 0;
  hrchFlag    = 0;
}


// Split character format data-field into status number values
void
IgesStatusField::setValues(char* data_str)
{
  for(int i=0; i<4; i++) {
    statusFldBuffer[0] = data_str[2*i];
    statusFldBuffer[1] = data_str[2*i + 1];
    statusNbrBuffer[i] = atoi(statusFldBuffer);
  }
  blankFlag = statusNbrBuffer[0];
  subordFlag  = statusNbrBuffer[1];
  useFlag   = statusNbrBuffer[2];
  hrchFlag    = statusNbrBuffer[3];
}


// Constructor
IgesDirectoryEntry::IgesDirectoryEntry()
{
  id = 0;
  entNbr = 0;
  startLine = 0;
  nofLines = 0;
  colorNbr = DEFAULT_COLOR_INDEX;
  formNbr = 0;
  transfId = 0;
  canBeBody = false;
  referingId = 0;
  body = NULL;
}


// Constructor
InputIges::InputIges(enum ecif_modelDimension m_dim,
                     ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
  modelDimension = inputDimension;

  directory = new IgesDirectory;
  paramSecLine  = 0;
  paramSecStart = 0;
}


// Method adds new directory entry to the container "directory"
int
InputIges::addToDirectory(IgesDirectoryEntry* dir_entry)
{
  int dir_id = dir_entry->id;
  (*directory)[dir_id] = dir_entry;

  return dir_id;
}


// Method checks if argument strin stream buffer is 'empty', it
// means if there are no data values in the line.
// NOTE: we suppouse that there are no blanks among data values !!!***!!!
bool
InputIges::dataLineStrmIsEmpty(istrstream& data_line)
{
  char c = data_line.peek();
  if (c == EOF || c == ' ')
    return true;
  return false;
}


bool
InputIges::checkEntryReferences()
{
  //--Go to the start position of the parameter (= data) section.
  infile.seekg(paramSecStart - infile.tellg(), ios::cur);
  paramSecLine = 1;

  bool check_only_status = true;

  //--Loop through stored dircetory entities and update references to
  //  other directory entries
  //  This is needed when we decide which entries defines bodies!
  IgesDirectory::iterator pos = directory->begin();

  bool ok;

  while (pos != directory->end()) {

    IgesDirectoryEntry* de = (*pos++).second;

    switch (de->entNbr) {

    case 100: ok = read_100(de, check_only_status); break;
    case 102: ok = read_102(de, check_only_status); break;
    case 106: ok = read_106(de, check_only_status); break;
    case 141: ok = read_141(de, check_only_status); break;
    case 142: ok = read_142(de, check_only_status); break;
    case 143: ok = read_143(de, check_only_status); break;
    case 144: ok = read_144(de, check_only_status); break;
    default: ok = true; break;
    }

    if (!ok) {
      return false;
    }
  }

  return true;
}


bool
InputIges::createBodies()
{
  IgesDirectory::iterator pos;

  // First create bodies
  pos = directory->begin();
  while (pos != directory->end()) {

    IgesDirectoryEntry* de = (*pos++).second;

    if ( !de->canBeBody || de->referingId > 0 ) {
      continue;
    }

    Body* body;

    if ( modelDimension == ECIF_2D ) {
      body = new Body2D();
    } else {
      body = new Body3D();
    }

    model->addBody(body);
    de->body = body;
  }

  // Next set parent bodies for the component elements
  pos = directory->begin();
  while (pos != directory->end()) {

    IgesDirectoryEntry* de = (*pos++).second;

    if ( de->referingId <= 0 ) {
      continue;
    }

    IgesDirectoryEntry* ref_de = getDirectoryEntry(de->referingId);

    de->body = ref_de->body;
  }

  return true;
}


// Method finds out the dimension of a Iges format Cad-model.
// One way to do this is to check if one of the (x,y,z) values
// is constant for all vertices. In this case model would be 2D.
// Otherwise it is 3D
//
enum ecif_modelDimension
InputIges::findCadModelDimension() {

  // Currently only 2D
  return ECIF_2D;

  //return ECIF_3D;
}


// Method reads one data field from string stream argument "data-line"
// into buffer "buffer".
void
InputIges::getDataField(istrstream*& data_line, char* buffer)
{
  // If there is no more data in the data-buffer, read next data-line
  if (dataLineStrmIsEmpty(*data_line)) {
    delete data_line;
    readDataLine(lineBuffer, dataBuffer);
    data_line = new istrstream(dataBuffer);
  }

  data_line->getline(buffer, DATA_PART_LEN, DATA_SEP);
}


// Methods retrieves directory entry from the container "directory"
// using the argument "entry-id" as teh key.
IgesDirectoryEntry*
InputIges::getDirectoryEntry(int entry_id)
{
  IgesDirectoryEntry* de =  (*directory)[entry_id];
  return de;
}


// Method reads directory field number "fld_nbr" into buffer "fld_buffer"
// Data in in the buffer "line_buffer"
void
InputIges::getDirectoryField(int fld_nbr, char* line_buffer, char* fld_buffer)
{
  strncpy(fld_buffer, line_buffer + (fld_nbr-1) * DIR_FLD_LEN, DIR_FLD_LEN);
}


// Method locates the starting position of the paramter lines
// (data lines) where the argument directory entry is pointing to.
// First line is read int global buffer "lineBuffer"
void
InputIges::locateParamEntry(IgesDirectoryEntry* de)
{
  int delta = de->startLine - paramSecLine;
  if (delta  == 0)
    return;
  else if (delta > 0) {
    paramSecLine += delta;
    readFileLine(infile, lineBuffer, delta);
    return;
  }
  else {
    infile.seekg(paramSecStart, ios::beg);
    paramSecLine = de->startLine;
    readFileLine(infile, lineBuffer, paramSecLine - 1);
  }
}


// Read those entries which define body elements
bool
InputIges::readBodyElements()
{
  //-Go to the start position of the parameter (= data) section.
  infile.seekg(paramSecStart - infile.tellg(), ios::cur);
  paramSecLine = 1;

  // Loop through directory entries
  IgesDirectory::iterator pos = directory->begin();

  while (pos != directory->end()) {

    IgesDirectoryEntry* de = (*pos++).second;

    // If this entry is not in the model at all!
    if ( de->body == NULL ) {
      continue;
    }

    bool ok;

    switch(de->entNbr) {
    case 100: ok = read_100(de); break;
    case 106: ok = read_106(de); break;
    case 110: ok = read_110(de); break;
    case 126: ok = read_126(de); break;
    case 128: ok = read_128(de); break;
    default: ok = true; break;
    }

    if (!ok) {
      return false;
    }
  }

  return true;
}




// Method reads all bodies in the Iges input file.
// Returns the number of the bodies read.
bool
InputIges::readCadGeometry()
{
  bool result;

  //--Read directory section and find the start position of the
  //  of parameters-section in the file.
  //  Paramter-section line-counter is set to 0.
  //  NOTE: This counter is updated by: readDataLine and
  //  locateParamEntry functions.
  paramSecStart = readDirectory();

  //--If model type is unknown
  if (modelDimension == ECIF_ND) {
    modelDimension = findCadModelDimension();
  }

  //--Update entry references ans states
  checkEntryReferences();

  //--Create bodies in each entry which defines a body
  createBodies();

  //--Read body elements
  readBodyElements();


  // After reading all body-information, data can be checked!
  //result = this->model->checkBodies();

  //  if (!modelIsOk)
  //    exit(1);
  return true;
}


// Method reads header information from the input file.
// (this is model level information (units etc.))
bool
InputIges::readCadHeader()
{
  return true;
}


// Method reads one data-line into buffer "data_buf" from Iges Cad-file.
// NOTE: end character (like ;) for the group of data lines
// is not read from the file.
void
InputIges::readDataLine(char* line_buf, char* data_buf)
{
  readFileLine(infile, line_buf);

  char* tmp = line_buf;
  int i = 0;
  // Copy data-part into data_buf. Don't copy data-end
  // separator character (normally ';')
  while (i < DATA_PART_LEN && *tmp != DATA_END)
    data_buf[i++] = *tmp++;

  data_buf[i] = '\0';

  // If we working om parameters-section, increase the
  // line counter.
  if (paramSecLine > 0)
    paramSecLine++;
}


// Method reads all (interesting)directory entries from Iges Cad-file
// and stores them into attribute "directory" (container)
int
InputIges::readDirectory()
{
  fieldBuffer[1+DIR_FLD_LEN] = '\0';
  lineBuffer[SEC_ID_POS] = ' ';

  while (!infile.eof() && lineBuffer[SEC_ID_POS] != 'D')
    readFileLine(infile, lineBuffer);

  if (infile.eof())
    return 0;

  IgesDirectoryEntry* de;
  int cur_file_pos = infile.tellg();

  while (!infile.eof() && lineBuffer[SEC_ID_POS] == 'D') {
    de = readDirectoryEntry(lineBuffer);

    //-Add geometry-use types to directory
    //if (de->status.useFlag == 0)
    //  addToDirectory(de);

    // Currently all entries are added to directory.
    // The problem is: what entries could be excluded? !!!
    addToDirectory(de);

    //-Read the first line of next entry
    cur_file_pos = infile.tellg();
    readFileLine(infile, lineBuffer);
  }
  return cur_file_pos;
}


IgesDirectoryEntry*
InputIges::readDirectoryEntry(char* dir_line)
{
  char* tmp;
  IgesDirectoryEntry* dir_entry = new IgesDirectoryEntry;

  // First directory-entry line
  // ==========================

  //--Entity type number
  getDirectoryField(1, dir_line, fieldBuffer);
  dir_entry->entNbr = atoi(fieldBuffer);

  //--Id
  getDirectoryField(10, dir_line, fieldBuffer);
  tmp = fieldBuffer;
  dir_entry->id = atoi(++tmp);

  //--Start-line in par-section
  getDirectoryField(2, dir_line, fieldBuffer);
  dir_entry->startLine = atol(fieldBuffer);

  //--Tranformation matrix id
  getDirectoryField(7, dir_line, fieldBuffer);
  dir_entry->transfId = atol(fieldBuffer);

  //--Status field
  getDirectoryField(9, dir_line, fieldBuffer);
  dir_entry->status.setValues(fieldBuffer);


  // Second directory-entry line
  // ===========================
  readFileLine(infile, dir_line);

  //--Color number
  getDirectoryField(3, dir_line, fieldBuffer);
  dir_entry->colorNbr = (colorIndices)atoi(fieldBuffer);

  //--Number of lines in par-section
  getDirectoryField(4, dir_line, fieldBuffer);
  dir_entry->nofLines = atoi(fieldBuffer);

  //--Form number
  getDirectoryField(5, dir_line, fieldBuffer);
  dir_entry->formNbr = atoi(fieldBuffer);

  //--Check if this entry could be a body
  switch (dir_entry->entNbr) {
  case 100:
  case 102:
  case 106:
  case 141:
  case 142:
  case 143:
  case 144:
    dir_entry->canBeBody = true;
    break;
  }

  return dir_entry;
}


int
InputIges::readDoubleFields(istrstream*& data_line, int nof_fields, double* buffer)
{
  //NOTE: we read data-lines until necessary number of
  //data fields is read from the file
  int counter;
  for(counter = 0; counter < nof_fields; counter++) {
    getDataField(data_line, dataFieldBuffer);
    buffer[counter] = atof(dataFieldBuffer);
  }

  return counter;
}


int
InputIges::readIntFields(istrstream*& data_line, int nof_fields, int* buffer)
{
  //NOTE: we read data-lines until necessary number of
  //data fields is read from the file
  int counter;
  for(counter = 0; counter < nof_fields; counter++) {
    getDataField(data_line, dataFieldBuffer);
    buffer[counter] = atol(dataFieldBuffer);
  }

  return counter;
}


// Method reads line body-element.
bool
InputIges::readLine(Body* body)
{
  return false;
}


// Method reads a nurbs-spline formed body-element.
bool
InputIges::readNurbs(Body* body)
{
  return false;
}


// Method reads one point from input-string *s*.
bool
InputIges::readPoint(istrstream*& data_line, Point3& p)
{
  for(int i = 0; i < 3; i++) {

    getDataField(data_line, dataFieldBuffer);

    p[i] = atof(dataFieldBuffer);

    p[i] *= model->unit.conv[i];
  }

  return true;
}


// Method reads one vertex from input-string *s*.
BodyElement*
InputIges::readVertex(istrstream*& data_line)
{
  static Point3 point;
  static GcPoint gpoint;

  if ( !readPoint(data_line, point) )
    return NULL;

  gpoint.setPosit(point[0], point[1], point[2]);

  // Create vertex from the point in the input file.
  BodyElement* vertex = new BodyElement1D(&gpoint);

  return vertex;
}


//******************************************
//****** IGES ENTITY READING ROUTINES ******
//******************************************

// Circular arc entity (Type 100)
bool
InputIges::read_100(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  static ecif_EdgeGeometry_X edge;
  init_trx_data(edge);

  edge.type = ECIF_CIRCLE;

  int i,j;
  int int_fld;
  double dbl_fldb[1];
  char chr_fld;

  bool is_closed;

  static Point3 points[3];

  locateParamEntry(de);

  //-Read one data-line
  readDataLine(lineBuffer, dataBuffer);
  istrstream* data_line = new istrstream(dataBuffer);

  //-Read center and two points
  *data_line >> int_fld >> chr_fld; // Entity-id
  readDoubleFields(data_line, 1, dbl_fldb); // z-inclination

  for(i = 0; i < 3; i++) {

    points[i][2] = 0.0;

    for (j = 0; j < 2; j++) {
      readDoubleFields(data_line, 1, dbl_fldb);
      points[i][j] = dbl_fldb[0];
    }
  }

  if ( isZero(dist3(points[1], points[2])) ) {
    is_closed = true;
  } else {
    is_closed = false;
  }

  if ( modelDimension == ECIF_2D && is_closed ) {
    de->canBeBody = true;
  }

  if ( check_only_status ) {
    return result;
  }

  // Create body element
  // ===================

  edge.location = new Point3[1];
  edge.start = new Point3[1];
  edge.end = new Point3[1];

  for (i = 0; i < 3; i++) {
    edge.location[0][i] = points[0][i];
    edge.start[0][i] = points[1][i];
    edge.end[0][i] = points[2][i];
  }

  int body_layer = 0;
  createBodyElement2D(de->body, body_layer, edge);

  reset_trx_data(edge);

  return result;
}


// Composite curve entity (Type 102)
bool
InputIges::read_102(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  int i;
  int int_fld, int_fldb[1];
  char chr_fld;

  istrstream* data_line;

  readDataLine(lineBuffer, dataBuffer);
  data_line = new istrstream(dataBuffer);


  //-Read number of curves in the composite
  int nof_crvs;
  *data_line >> int_fld >> chr_fld; // Entry-id
  *data_line >> nof_crvs >> chr_fld;

  int* crv_id_table = new int[nof_crvs];

  IgesDirectoryEntry* comp_de;

  //-Read component dir-ids and update the reference
  for(i = 0; i < nof_crvs; i++) {
    readIntFields(data_line, 1, int_fldb);
    crv_id_table[i] = int_fldb[0];

    comp_de = getDirectoryEntry(int_fldb[0]);
    comp_de->referingId = de->id;
  }

  delete [] crv_id_table;

  return result;
}


// Copious data entity (Type 106)
// Path of data points
bool
InputIges::read_106(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  bool is_closed = ( de->formNbr == 63 );

  if ( modelDimension == ECIF_2D && is_closed ) {
    de->canBeBody = true;
  }

  if ( check_only_status ) {
    return result;
  }

  // Read geometry data
  // ==================
  int i,j;
  int int_fld;
  char chr_fld;
  double dbl_fld[1];

  istrstream* data_line;

  readDataLine(lineBuffer, dataBuffer);
  data_line = new istrstream(dataBuffer);

  //-Read format of points
  int read_size;
  int point_size;
  int point_type;
  int nof_points;
  double common_z;

  *data_line >> int_fld >> chr_fld;         // Entry-id
  *data_line >> point_type >> chr_fld;      // Point type
  *data_line >> nof_points >> chr_fld;      // Nof tuples
  readDoubleFields(data_line, 1, dbl_fld);  // Common z coordinate (for type 1 points)

  if ( point_type == 1 ) {
    point_size = 2;
    read_size = 2;
  } else if ( point_type == 2 ) {
    point_size = 3;
    read_size = 3;
  } else if ( point_type == 3 ) {
    point_size = 6;
    read_size = 3;
  }

  Point3* points = new Point3[nof_points];

  //-Read points
  int read_count = 0;
  for(i = 0; i < nof_points; i++) {

    points[i][0] = points[i][1] = points[i][2] = 0.0;

    for (j = 0; j < read_size; j++) {

      readDoubleFields(data_line, 1, dbl_fld);
      read_count++;

      // If we are reading relevant point data
      if ( read_count <= read_size ) {
        points[i][read_count - 1] = dbl_fld[0];
      }

      // Reset counter if all point data read
      read_count = read_count % point_size;
    }
  }

  IgesDirectoryEntry* comp_de;

  //-We are reading a surface into a 3D model
  if ( modelDimension == ECIF_3D ) {
    // body->addElement(from this face-element etc....);
    result = false;
  }

  //-We are reading a body into a 2D model!
  else if ( modelDimension == ECIF_2D ) {
    int body_layer = 0;
    createBodyElement2D(de->body, body_layer, nof_points, points);
  }

  delete[] points;

  return result;
}


// Line entity (Type 110)
bool
InputIges::read_110(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  int i,j;
  int int_fld;
  char chr_fld;

  static Point3 points[2];

  locateParamEntry(de);

  //-Read one data-line
  readDataLine(lineBuffer, dataBuffer);
  istrstream* data_line = new istrstream(dataBuffer);

  //-Read two points
  *data_line >> int_fld >> chr_fld; // Entity-id

  for(i = 0; i < 2; i++) {
    readPoint(data_line, points[i]);
  }

  int body_layer = 0;
  createBodyElement2D(de->body, body_layer, points[0], points[1]);

  return result;
}


// Rational B-spline curve entity (Type 126)
bool
InputIges::read_126(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  static ecif_EdgeGeometry_X edge;
  init_trx_data(edge);

  locateParamEntry(de);
  int i,j, int_fld;
  char chr_fld;

  //-Read one data-line
  readDataLine(lineBuffer, dataBuffer);
  istrstream* data_line = new istrstream(dataBuffer);

  int nof_cp, nof_kp, degree;
   int planar_flg, open_flg, polyn_flg, periodic_flg;

  *data_line >> int_fld >> chr_fld;   //-Entry-id
  *data_line >> nof_cp >> chr_fld;      //-Nof control-points - 1
  nof_cp += 1;
  *data_line >> degree >> chr_fld;      //-Degree of basis functions
  nof_kp = nof_cp - degree + 1;       //-Number of knot-points
  nof_kp += 2 * degree;           // multiplicity at ends added!!
  *data_line >> planar_flg >> chr_fld;  //-0= non-planar
  *data_line >> open_flg >> chr_fld;    //-0= open curve
  *data_line >> polyn_flg >> chr_fld;   //-0=rational, 1=polynomial
  *data_line >> periodic_flg >> chr_fld;  //-0=non-periodic

  //-Read data into a temporary data-vector
  //Data length: (control-points (x,y,z,w)+ knots
  int data_len = 4 * nof_cp + nof_kp;
  double* crv_data = new double[data_len];
   int nof_values = readDoubleFields(data_line, data_len, crv_data);

  //-Copy data from crv_data into proper components
  double* knots = new double[nof_kp];
  Point4* ct_points = new Point4[nof_cp];

  int pos = 0;

  //-Knot points are first
  for(i=0; i < nof_kp; i++)
    knots[i] = crv_data[pos++];

  //-Then come weights
  for(i=0; i < nof_cp; i++) {
    ct_points[i][3] = crv_data[pos++];
  }

  //-Then come control points

  edge.start = new Point3[1];
  edge.end = new Point3[1];

  for(i=0; i < nof_cp; i++) {

    for(j=0; j < 3; j++) {

      //-First control point is the start-point!
      if ( i == 0 ) {
        edge.start[0][j] = crv_data[pos];
      }

      //-Last control point is the end-point!
      if ( i == nof_cp - 1 ) {
        edge.end[0][j] = crv_data[pos];
      }

      //NOTE: In Iges points are not in homogenous form
      // We must transform them into (xw,yw,zw,w) format
      ct_points[i][j] = ct_points[i][3] * crv_data[pos++];
    }
  }

  edge.isRational = (polyn_flg == 0);
  edge.degree = degree;
  edge.nofKnots = nof_kp;
  edge.knots = knots;
  edge.nofCpoints = nof_cp;
  edge.cpoints = ct_points;

  int body_layer = 0;
  createBodyElement2D(de->body, body_layer, edge);

  delete[] crv_data;

  reset_trx_data(edge);

  return result;
}


// Rational B-spline surface entity (Type 128)
bool
InputIges::read_128(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  int i,j;
  int int_fld;
  char chr_fld;

  //-Read one data-line
  readDataLine(lineBuffer, dataBuffer);
  istrstream* data_line = new istrstream(dataBuffer);

  int nof_cp, polyn_flg;
  int nof_cp_u, nof_kp_u, degree_u;
  int open_flg_u, periodic_flg_u;
  int nof_cp_v, nof_kp_v, degree_v;
  int open_flg_v, periodic_flg_v;

  *data_line >> int_fld >> chr_fld;   //-Entry-id

   //---u-parameter
  *data_line >> nof_cp_u >> chr_fld;    //-Nof control-points - 1
  nof_cp_u += 1;
  *data_line >> degree_u >> chr_fld;    //-Degree of basis functions
  nof_kp_u = nof_cp_u - degree_u + 1;   //-Number of knot-points
  nof_kp_u += 2 * degree_u;         // multiplicity at ends added!!

   //---v-parameter
  *data_line >> nof_cp_v >> chr_fld;    //-Nof control-points - 1
  nof_cp_v += 1;
  *data_line >> degree_v >> chr_fld;    //-Degree of basis functions
  nof_kp_v = nof_cp_v - degree_v + 1;   //-Number of knot-points
  nof_kp_v += 2 * degree_v;         // multiplicity at ends added!!

   //---other parameters
  *data_line >> open_flg_u >> chr_fld;    //-0= open curve
  *data_line >> open_flg_v >> chr_fld;    //-0= open curve
  *data_line >> polyn_flg >> chr_fld;      //-0=rational, 1=polynomial
  *data_line >> periodic_flg_u >> chr_fld;  //-0=non-periodic
  *data_line >> periodic_flg_v >> chr_fld;  //-0=non-periodic

  //-Read data into a temporary data-vector
  //Data length: (control-points (x,y,z,w)+ knots
   // total number of control-point (cp_u * cp_v array)
  nof_cp = nof_cp_u * nof_cp_v;
  int data_len =  nof_kp_u + nof_kp_v + 4 * nof_cp ;
  double* crv_data = new double[data_len];
  int nof_values = readDoubleFields(data_line, data_len, crv_data);

  //-Copy data from crv_data into proper components
  double* knots_u = new double[nof_kp_u];
  double* knots_v = new double[nof_kp_v];
  Point4* cpoints = new Point4[nof_cp];

  int pos = 0;

  //-Knot-u points are first
  for(i=0; i < nof_kp_u; i++)
    knots_u[i] = crv_data[pos++];

  //-Knot-u points are next
  for(i=0; i < nof_kp_v; i++)
    knots_v[i] = crv_data[pos++];

  //-Weights are next
  for(i=0; i < nof_cp; i++) {
    cpoints[i][3] = crv_data[pos++];
  }

  //-Control points are next
  //NOTE: In Iges points are not in homogenous form
  // We must transform them into (xw,yw,zw,w) format
  for(i=0; i < nof_cp; i++) {
    for(j=0; j < 3; j++)
      cpoints[i][j] = cpoints[i][3] * crv_data[pos++];
  }

  ecif_FaceGeometry_X sP;
  sP.isRational = (polyn_flg == 0);
  sP.degree_u = degree_u;
  sP.degree_v = degree_v;
  sP.nofKnots_u = nof_kp_u;
  sP.nofKnots_v = nof_kp_v;
  sP.knots_u = knots_v;
  sP.knots_v = knots_v;
  sP.nofCpoints_u = nof_cp_u;
  sP.nofCpoints_v = nof_cp_v;
  sP.nofCpoints = nof_cp;
  sP.cpoints = cpoints;

  delete [] crv_data;

  //-Create nurbs-curve element
  //--Corner points
  BodyElement* vertices[4];
  int  vertex_ids[4];

  static GcPoint point;

  int corner_ids[4];
  corner_ids[0] = 0;
  corner_ids[1] = nof_cp_u - 1;
  corner_ids[2] = (nof_cp_u - 1) * nof_cp_v;
  corner_ids[3] = nof_cp_u * nof_cp_v - 1;

  // Create model vertex elements from nurbs corner points
  for (i = 0; i < 4; i++) {
    int c_id = corner_ids[i];
    double w = sP.cpoints[c_id][3];

    point.setPos(X, sP.cpoints[c_id][0]/w);
    point.setPos(Y, sP.cpoints[c_id][1]/w);
    point.setPos(Z, sP.cpoints[c_id][2]/w);

    vertices[i] = new BodyElement1D(&point);
    vertex_ids[i] = vertices[i]->Id();

    model->addBodyElement(vertices[i]);
  }

  // Create a new element into the 3D-body
  BodyElement* body_elm = new BodyElement3D(4, vertex_ids, &sP);
  de->body->addElement(0, body_elm);
  model->addBodyElement(body_elm);

  return result;
}


// Boundary entity (Type 141)
// This entry point to a untrimmed surface and to all curve elements
// which define this surface as boundaries.
bool
InputIges::read_141(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  int i, j;
  int int_fld, int_fldb[1];
  char chr_fld;

  istrstream* data_line;

  readDataLine(lineBuffer, dataBuffer);
  data_line = new istrstream(dataBuffer);

   //-Skip 4 (nbr sep) pairs:
   // Entry id
   // Type of srf representation (0=in trimmimg space, 1=in model space)
   // Preferred representation (0=unspec., 1=model space, 2=param. spcace, 3=equal)
   // Pointer to the DE of the untrimmed surface
  for(i = 0; i < 4; i++) {
    *data_line >> int_fld >> chr_fld;
  }

  //-Read number of curves in the boundary
  int nof_crvs;
  *data_line >> nof_crvs >> chr_fld;
  int* crv_id_table = new int[nof_crvs];

  IgesDirectoryEntry* comp_de;

  //-Read component ids and update references
  for(i = 0; i < nof_crvs; i++) {
    for(j = 0; j < 4; j++) {

      readIntFields(data_line, 1, int_fldb);

      if (j == 0) {
        crv_id_table[i] = int_fldb[0];
        comp_de = getDirectoryEntry(int_fldb[0]);
        comp_de->referingId = de->id;
      }
    }
  }

  delete [] crv_id_table;
  return result;
}


// Curve on a parametric surface entity (Type 142)
bool
InputIges::read_142(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  int int_fld;
  char chr_fld;
  istrstream* data_line;

  readDataLine(lineBuffer, dataBuffer);
  data_line = new istrstream(dataBuffer);

  int prmc_entry;
  int srf_entry;
  int srfc_entry;
  *data_line >> int_fld >> chr_fld;     //-Entry-number
  *data_line >> int_fld >> chr_fld;     //-Type of the curve
  *data_line >> srf_entry >> chr_fld;   //-Pointer to surface
  *data_line >> prmc_entry >> chr_fld;  //-Pointer to crv in param space
  *data_line >> srfc_entry >> chr_fld;  //-Pointer to crv on the surface
  *data_line >> int_fld >> chr_fld;     //-Preferred format

  // Update reference to the component
  IgesDirectoryEntry* comp_de = getDirectoryEntry(srfc_entry);
  comp_de->referingId = de->id;

  return result;
}


// This entry defines a bounded surface entity (Type 143)
// It points to the surface and to the entry which defines the
// bounding curves (fex. to an entry 141)
bool
InputIges::read_143(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  int i;
  int int_fld;
  char chr_fld;
  istrstream* data_line;

  readDataLine(lineBuffer, dataBuffer);
  data_line = new istrstream(dataBuffer);

  int nof_srfs, nof_bnds, nof_values;

  *data_line >> int_fld >> chr_fld;   // entry-number

  *data_line >> nof_srfs >> chr_fld;    // Nof surface entities
  int* srf_id_table = new int[nof_srfs];
  nof_values = readIntFields(data_line, nof_srfs, srf_id_table);

  *data_line >> nof_bnds >> chr_fld;     // Nof boundary entities
  int* bnd_id_table = new int[nof_bnds];
  nof_values = readIntFields(data_line, nof_bnds, bnd_id_table);

  IgesDirectoryEntry* comp_de;
  int nof_refs;
  int* refs_table;

  // In 2D boundaries
  if ( modelDimension == ECIF_2D ) {
    nof_refs = nof_bnds;
    refs_table = bnd_id_table;

  // In 3D surfaces
  } else {
    nof_refs = nof_srfs;
    refs_table = srf_id_table;
  }

  for(i = 0; i < nof_refs; i++) {
    comp_de = getDirectoryEntry(refs_table[i]);
    comp_de->referingId = de->id;
  }

  delete [] bnd_id_table;
  delete [] srf_id_table;

  return result;
}


// Trimmed (parametric) surface entity (Type 144)
bool
InputIges::read_144(IgesDirectoryEntry* de, bool check_only_status)
{
  bool result = true;

  locateParamEntry(de);

  int i;
  int int_fld, int_fldb[1];
  char chr_fld;

  istrstream* data_line;

  readDataLine(lineBuffer, dataBuffer);
  data_line = new istrstream(dataBuffer);

  int srf_entry;
  int obnd_entry;
  int nof_ibnds;

  *data_line >> int_fld >> chr_fld;     //-Entry-number
  *data_line >> srf_entry >> chr_fld;   //-Pointer to surface
  *data_line >> int_fld >> chr_fld;     //-Type of bounded surface
                                          // 0=the outer bndr is the bndr of D
                                          // 1=otherwise
  *data_line >> nof_ibnds >> chr_fld;   //-Nof inner bndrs defined by simple closed curves
  *data_line >> obnd_entry >> chr_fld;  //-Pointer to outer bndr of the surface
                                        // defining the inner bndr of the srf.

  IgesDirectoryEntry* comp_de;

  //-Read outer boundary on the surface (if defined!)
  if (obnd_entry > 0) {
    comp_de = getDirectoryEntry(obnd_entry);
    comp_de->referingId = de->id;
  }

  //-Read inner boundary defining curves (if defined!)
  if (nof_ibnds > 0) {

    int* bnd_id_table = new int[nof_ibnds];
    int ibnd_entry;

    //-Read component row-ids and corresponding data-rows
    for(i = 0; i < nof_ibnds; i++) {
      readIntFields(data_line, 1, int_fldb);
      bnd_id_table[i] = int_fldb[0];

      comp_de = getDirectoryEntry(int_fldb[0]);
      comp_de->referingId = de->id;
    }

    delete [] bnd_id_table;
  }

  return result;
}
