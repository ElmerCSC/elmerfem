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
Module:     ecif_inputIdeasWF.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation for Ideas WireFrame bodies

************************************************************************/

#include "ecif_body2D.h"
#include "ecif_body3D.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_geometry.h"
#include "ecif_inputIdeasWF.h" 
#include "ecif_model.h"

extern char read_buffer[];

  
//*** Constants for reading Ideas data in Nurbs-Format.
// Number of data-items (points) per line in the input file.
const int NF_COLS   = 3;
// Length of nurbs-description vector.
const int NF_SIZE   = 12;
// Indecies in the nurbs-description vector (0 - NF_SIZE-1)
const int NFI_N     = 0;  // Number of control-points.
const int NFI_NK    = 1;  // Number of knot-points (N+K).
const int NFI_K     = 3;  // Order of basis.
const int NFI_LEN   = 10; // Length of curve's control data.
const int NFI_DIM   = 11; // Dimension (of control points) (2-4).


InputIdeasWF::InputIdeasWF(enum ecif_modelDimension m_dim,
                           ifstream& in_file, char* in_filename):
        InputIdeas(m_dim, in_file, in_filename)
{
  modelDimension = inputDimension;
}


//***Read all wireframe bodies from a Ideas univ. file dataset (in fact a 2D case!!)
bool
InputIdeasWF::readCadGeometry()
{
  while (!infile.eof()) {
    // Read all bodies.
    int status = readWireFrameBody(); 
    // If end-of-dataset was encountered.
    if (status == -1)
      break;
  }   
  return true;
}


// Method reads body-name from Ideas univ. file
char*
InputIdeasWF::readBodyName(char* s)
{
  int dummy;
  char* name = new char[MAXBODYNAME+1];
  char* cptr;

  istrstream strline(s);  
  // Body 'name' is field-4 (simply color code)
  strline >> dummy >> dummy >> dummy;   
  cptr = name;
  int i = 0; char c;
  // Trim left all blanks
  while (strline.get() == ' ');
  // read all non-blanks into name
  while (((c = strline.get()) != ' ') && (i < MAXBODYNAME)) {
    *cptr = c;
    cptr++; i++;
  }
  *cptr = '\0';
  return name;
}


// Method reads body-nbr from Ideas univ. file
int
InputIdeasWF::readBodyNbr(char* s) 
{
  int dummy;
  int nbr;

  istrstream strline(s);
  // Body nbr is field-4 (color code)
    strline >> dummy >> dummy >> dummy >> nbr;
  
  return nbr;
}

                                               
// Method reads curve-type info from Ideas univ. file (801)
ecif_geometryType
InputIdeasWF::readGeomType(char* s)
{
  ecif_geometryType g_type = ECIF_LINE;
  int dummy;
  int cnstr_flag;
  int nurbs_flag;

  istrstream strline(s);
  // Points/constraints-exist flag is field-3 (0=no, 1=yes).
  strline >> dummy >> dummy >> cnstr_flag;  
  // Nurbs-format flag is field-5 (0=no, 1=yes)
  strline >> dummy >> nurbs_flag;   

  if (nurbs_flag == 1)
    g_type = ECIF_NURBS;

  return g_type;
}
 

// Method reads a line-segment-element from Ideas univ. file (801)
bool
InputIdeasWF::readLine(Body* body, char* buffer)
{
  static Point3 p1, p2;

  // First read Record-3 away.
  // (it is : 1 1 and doesn't contain useful info)
  readFileLine(infile, buffer);

  // Next read the two vertices from infile.

  // -first vertex-point
  readFileLine(infile,buffer);
  readPoint(buffer, p1);

  // -second vertex-point
  readFileLine(infile,buffer);
  readPoint(buffer, p2);

  int body_layer = 0;

  // Create a new 2D element into the body
  createBodyElement2D(body, body_layer, p1, p2);

  return true;
}


// Method reads a spline-segment-element in nurbs-form from Ideas univ. file (801)
bool
InputIdeasWF::readNurbs(Body* body, char* buffer)
{
  static ecif_EdgeGeometry_X edge;
  init_trx_data(edge);

  int i,j,pos;
  // *** Control-point section.
  int nf_vec[NF_SIZE];

  // Read curve-desrcribing vector.
  readFileLine(infile, buffer);
  istrstream strline(buffer);
  for (i = 0; i < NF_SIZE; i++) {
    strline >> nf_vec[i];
  }
  
  // Pick curve-defining factors from nf-vector.
  int nf_dim  = nf_vec[NFI_DIM];  // Dimension of control-points.
  int nf_len  = nf_vec[NFI_LEN];  // Length of data.
  int nf_n  = nf_vec[NFI_N];    // Nof control points.
  int nf_nk   = nf_vec[NFI_NK];   // Nof knot points.
  int nf_k  = nf_vec[NFI_K];    // Order of basis.

  // Read nurbs-curve defining data into a temporary data-vector.
  double* ct_data = new double[nf_len]; // a temporary data vector.
  
  int nf_rows; // nof data-rows in this section.
  if ( NF_COLS > 0 )
    nf_rows = nf_len / NF_COLS;
  else
    nf_rows = nf_len;

  pos = 0;

  for (i = 0; i < nf_rows; i++) {
    readFileLine(infile, buffer);
    istrstream strline(buffer);
    for (j = 0; j < NF_COLS; j++) {
      strline >> ct_data[pos++];
    }
  }

  // Now create final data-structures for nurbs-curve and
  // read them from temporary data-vector *ct_data*.

  // *** Control points (at the beginning of the data)
  Point4* ct_points = new Point4[nf_n];

  pos = 0;
  for (i = 0; i < nf_n; i++) {

    //We use always 4D-points
    int pos = 4 * i;
    ct_points[i][0] = 0.0;
    ct_points[i][1] = 0.0;
    ct_points[i][2] = 0.0;
    ct_points[i][3] = 1.0;

    for (j = 0; j < nf_dim; j++) {
      //Note: coordinates are scaled by unit
      ct_points[i][j] = model->unit.conv[j] * ct_data[pos++];
    }
  }

  // *** Next come knot-points
  double* knots = new double[nf_nk]; 

  pos = nf_n * nf_dim;
  for (i = 0; i < nf_nk; i++) {
    knots[i] = ct_data[pos + i];
  }

  //---Create a nurbs-curve body-element

  //-Two vertices
  // first vertex from the first control-point
  // second vertex from the last control-point
  edge.start = new Point3[1]; 
  edge.end = new Point3[1];

  for (i = 0; i < 3; i++) {
    edge.start[0][i] = ct_points[0][i] / ct_points[0][3];
    edge.end[0][i] = ct_points[nf_n - 1][i] / ct_points[nf_n - 1][3];
  }

  //-Other parameters
  edge.type = ECIF_NURBS;
  edge.isRational = (nf_dim == 4); // Is this ok !!!###!!!
  edge.degree = nf_k - 1;
  edge.nofKnots = nf_nk;
  edge.nofCpoints = nf_n;
  edge.knots = knots;
  edge.cpoints = ct_points;
  
  //-Create a new 2D element into the body
  int body_layer = 0;
  createBodyElement2D(body, body_layer, edge);

  // *** Defining points section.
  // Read nof defining point
  // NOT IN USE in practice. Used only to read 'off' the data from file !!!***!!!
  int df_n;
  readFileLine(infile, buffer);
  strline >> df_n;

  // Read defining-points desrcibing vector. All data on a single line.
  // This is a table where there are three items per defining-point:
  // 1. curve passes through point flag (0=no, 1=yes).
  // 2. tangent/derivative vector specified (0=no, 1=tangent, 2=derivative).
  // 3. curvature specified (0=no, 1=yes).
  readFileLine(infile, buffer);

  // Jump to the end of defining points section.
  readFileLine(infile, buffer, df_n);

  reset_trx_data(edge);

  return true;
}


// Method reads one wireframe body-element from Ideas univ. file dataset
// and adds body-element into body/model.   
// Returns: 1=new body was created, 0=old body was updateds, -1 = end-of-odataset.
int
InputIdeasWF::readWireFrameBody()
{
  Body* body;

  // We want to start from this position, when the body-read loop is started!
  //streampos cur_pos = infile.tellg();
 
  // =====Read body ID, CODE and NAME (from color code!) ( Record-1).
  readFileLine(infile, read_buffer);

  // However, check first that ***end-of-datset*** (-1) is not encountered
  // In this case dataset is empty!
  if (endofDataset(read_buffer))
    return -1;

  // Now read id, code and name.
  //char* name = readBodyName(read_buffer);
  int id_nbr = readBodyNbr(read_buffer);
  char* name = NULL;

  // =====Read CURVE-TYPE information ( Record-2).
  readFileLine(infile, read_buffer);
  ecif_geometryType g_type = readGeomType(read_buffer);

  // Back to body's starting position!
  //infile.seekg(cur_pos);

  // We check if we are reading new pieces for an old body!
  // Checking is based on a *local table* where we have id-numbers
  // from the file (id_nbr) and system (automagically) generated
  // technical id-numbers (body->Tag()) as pairs.
  bool isNewBody = false;
  IdNumberTable::iterator itr = bodyNumbers.find(id_nbr);

  // Old body
  if (itr != bodyNumbers.end())
    body = model->getBodyByTag((*itr).second);

  // Create a new body  
  else {
    colorIndices color = (colorIndices)id_nbr; // id_nbr is read from color code!

    //2D geometry
    if (modelDimension == ECIF_2D)
      body = new Body2D(GEOM_BODY, id_nbr, name, color);

    //3D geometry
    else if (modelDimension == ECIF_3D)
        body = new Body3D(GEOM_BODY, id_nbr, name, color);

    model->addBody(body);

    isNewBody = true;
    bodyNumbers[id_nbr] = body->Tag();
  }

  // How to read depends on the geometry-type:
  switch (g_type) {
  case ECIF_LINE:
    readLine(body, read_buffer);
    break;
  case ECIF_NURBS:
    readNurbs(body, read_buffer);
    break;
  default:
    break;
  }

  if (isNewBody)
    return 1;
  else
    return 0;
}
