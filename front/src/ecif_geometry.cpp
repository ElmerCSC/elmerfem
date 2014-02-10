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
Module:     ecif_geometry.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_def_stl.h"
#include "ecif_boundbox.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement1D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_model.h"
#include "ecif_nurbs.h"
#include "ecif_renderer.h"
#include "ecif_userinterface.h"

// Typedefs for safer calling
typedef int (Callbackp dllCircleFunc)(enum egf_FuncType ft, int argc, double* argv,
                             egf_Point3* start, egf_Point3* end,
                             int& nof_circles, egf_Circle*& circles);

typedef int (Callbackp dllPolyLineFunc)(enum egf_FuncType ft, int argc, double* argv,
                             egf_Point3* start, egf_Point3* end,
                             int& nof_polylines, egf_PolyLine*& polylines);


// Static class variables
//
Model* Geometry::model = NULL;

ecif_geometryType GcCircle::type = ECIF_CIRCLE;
ecif_geometryType GcLine::type = ECIF_LINE;
ecif_geometryType GcPlane::type = ECIF_LINE;
ecif_geometryType GcPoint::type = ECIF_NODIM;
ecif_geometryType GcPolyLine::type = ECIF_POLYLINE;
ecif_geometryType GcNurbsCurve::type = ECIF_NURBS;
ecif_geometryType GcNurbsSurface::type = ECIF_NURBS;

//extern const float colorValues[][4];


RayHit::RayHit() {
  count = 0;
  min_value = 1.0e+100;
}


RayHit::~RayHit() {
  Point3List::iterator itr;

  for (itr = points.begin(); itr != points.end(); itr++) {
    delete *itr;
  }
  points.clear();
}


//***************************
//*****Geometry methods *****
//***************************

Geometry::Geometry() {
 boundbox = NULL;
 geometryOk = true;
}


Geometry::~Geometry()
{
  delete boundbox;
}


void
Geometry::calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u, double delta_v)
{
  nof_points_u = 0;
  nof_points_v = 0;
  points = NULL;
}


// Just cover function to call a more specific form
void
Geometry::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction)
{
  draw(renderer, dmode, dstate);
}


// Just cover function to call a more specific form
void
Geometry::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  draw(renderer, dmode, dstate, direction);
}


// Just cover function to call a more specific form
void
Geometry::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id, bool is_first_loop)
{
  draw(renderer, dmode, dstate, direction, elem_id);
}

// Just cover function to call a more specific form
void
Geometry::draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction)
{
  if ( gmtr_index < 0 ) return;

  if ( isMultiGeometry() ) {
    draw(gmtr_index, renderer, dmode, dstate);
  } else {
    draw(renderer, dmode, dstate, direction);
  }
}


// Just cover function to call a more specific form
void
Geometry::draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  if ( gmtr_index < 0 ) return;

  if ( isMultiGeometry() ) {
    draw(gmtr_index, renderer, dmode, dstate, direction);
  } else {
    draw(renderer, dmode, dstate, direction, elem_id);
  }
}


// Just cover function to call a more specific form
void
Geometry::draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id, bool is_first_loop)
{
  if ( gmtr_index < 0 ) return;

  if ( isMultiGeometry() ) {
    draw(gmtr_index, renderer, dmode, dstate, direction, elem_id);
  } else {
    draw(renderer, dmode, dstate, direction, elem_id, is_first_loop);
  }
}


void
Geometry::getBoundaryPoints(int& count, BoundaryPoint**& points)
{
  points = NULL;
  count = 0;
}

 
void
Geometry::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 0;
  types = NULL;
  valuesU = NULL;
  valuesV = NULL;
  useFixedN = NULL;
}


// Finds 'nice' positions for the label point
void
Geometry::getLabelPoint(Point3& point)
{
  GcPoint* p = param2Point(0.5, 0.5);

  if (p != NULL) {
    point[0] = p->Pos(X);
    point[1] = p->Pos(Y);
    point[2] = p->Pos(Z);
  }

  else {
    point[0]  = point[1] = point[2] = NSVD;
  }
}


void
Geometry::getMifTags(int& nof_tags, int*& tags)
{
  nof_tags = 0;
  tags = NULL;
}




//**************************************************************************
//                        Geometry1D methods
//**************************************************************************

Geometry1D::Geometry1D()
{
}


Geometry1D::~Geometry1D()
{
}

 

//****************************
//***** GcPoint methods *****
//****************************

GcPoint::GcPoint()
{
    posit[X] = 0;
    posit[Y] = 0;
    posit[Z] = 0;
}


GcPoint
GcPoint::average(GcPoint &p2, GcPoint&p3)
{
  GcPoint res;

  for (char i = 0; i < MAX_DIMENSION; i++)
    res.posit[i] = (posit[i] + p2.posit[i] + p3.posit[i]) / 3;
  return res;
}


void
GcPoint::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate)
{
  renderer->drawPoint(dmode, dstate, &posit);
}


// Searching for space: if GcPoint is outside space extend the space
void
GcPoint::find_space(GcPoint &pmin, GcPoint &pmax)
{
  for (char i = 0; i < MAX_DIMENSION; i++) {
    if (pmin.posit[i] > posit[i])
      pmin.posit[i] = posit[i];
    if (pmax.posit[i] < posit[i])
      pmax.posit[i] = posit[i];
  }
}


// Finds 'nice' positions for the label point
void
GcPoint::getLabelPoint(Point3& point)
{
  point[0] = posit[0];
  point[1] = posit[1];
  point[2] = posit[2];
}


int
GcPoint::hashKey()
{
  return int(posit[X]);
}


bool
GcPoint::isBetween(GcPoint* p1, GcPoint* p2)
{
  GcPoint tmp1 = *this - *p1;
  GcPoint tmp2 = *p2 - *this;
  GcPoint tmp  = *p2 - *p1;
  return ((tmp1 < tmp) && (tmp2 < tmp));
}

bool
GcPoint::isParallel(GcPoint* other_vector)
{
  //We check if cross product is zero-vector.
  GcPoint* ov = other_vector;
  double tmp1 = fabs(posit[Y] * ov->posit[Z] - posit[Z] * ov->posit[Y]);
  double tmp2 = fabs(posit[X] * ov->posit[Z] - posit[Z] * ov->posit[X]);
  double tmp3 = fabs(posit[X] * ov->posit[Y] - posit[Y] * ov->posit[X]);
  return ((tmp1 < POINT_EPSILON) && (tmp2 < POINT_EPSILON) && (tmp3 < POINT_EPSILON));
}


void
GcPoint::normalize(double norm[3], double shift[3])
{
  for (int i = 0; i < 3 ; i++) {
    if ( norm[i] > 0 )
      posit[i] = (posit[i] + shift[i]) / norm[i];
  }
}


void
GcPoint::de_normalize(double norm[3], double shift[3])
{
  posit[0] = (posit[0] * norm[0]) - shift[0];
  posit[1] = (posit[1] * norm[1]) - shift[1];
  posit[2] = (posit[2] * norm[2]) - shift[2];
}



ostream&
GcPoint::output_emf(ostream& out, short indent_size, short indent_level)
{
  LibFront::output_vector(out, indent_size, indent_level, "Point", NULL, 3, posit, false);

  return out;
}


// *** System operators
bool
GcPoint::operator==(const GcPoint& p2)
{
  GcPoint hlp;
  double tmp;

  hlp = *this - p2;
  if ((hlp.posit[X] < POINT_EPSILON) && (hlp.posit[X] > -POINT_EPSILON) &&
      (hlp.posit[Y] < POINT_EPSILON) && (hlp.posit[Y] > -POINT_EPSILON) &&
      (hlp.posit[Z] < POINT_EPSILON) && (hlp.posit[Z] > -POINT_EPSILON))
    return true;
  else
    return false;

  tmp = sqrt(hlp.posit[X] * hlp.posit[X] + hlp.posit[Y] * hlp.posit[Y]
        + hlp.posit[Z] * hlp.posit[Z]);

  if (tmp < POINT_EPSILON && tmp > -POINT_EPSILON)
    return true;
  else
    return false;

}


ostream &
operator<<(ostream& stream, const GcPoint& p)
{
  stream << p.posit[X] << " ";
  stream << p.posit[Y] << " ";
  stream << p.posit[Z];
  return stream;
}

// *** Arithmetic operators
GcPoint
operator+(const GcPoint& p1, const GcPoint& p2)
{
  double p[MAX_DIMENSION];

  for (short i = 0; i < MAX_DIMENSION; i++)
    p[i] = p1.posit[i] + p2.posit[i];
  return GcPoint(p[X], p[Y], p[Z]);
}

GcPoint
operator-(const GcPoint& p1, const GcPoint& p2)
{
  double p[MAX_DIMENSION];

  for (short i = 0; i < MAX_DIMENSION; i++)
    p[i] = p1.posit[i] - p2.posit[i];
  return GcPoint(p[X], p[Y], p[Z]);
}

GcPoint
operator*(int num, const GcPoint& p)
{
  double tmp[MAX_DIMENSION];

  for (short i = 0; i < MAX_DIMENSION; i++)
    tmp[i] = num * p.posit[i];
  return GcPoint(tmp[X], tmp[Y], tmp[Z]);
}

GcPoint
operator/(const GcPoint& p, int num)
{
  double tmp[MAX_DIMENSION];

  for (short i = 0; i < MAX_DIMENSION; i++)
    if ( num != 0 )
      tmp[i] = p.posit[i] / num;

  return GcPoint(tmp[X], tmp[Y], tmp[Z]);
}


// *** Relational operators
bool
operator==(const GcPoint &p1, const GcPoint &p2)
{
  GcPoint hlp;

  hlp = p1 - p2;
  if ((hlp.posit[X] < POINT_EPSILON) && (hlp.posit[X] > -POINT_EPSILON) &&
      (hlp.posit[Y] < POINT_EPSILON) && (hlp.posit[Y] > -POINT_EPSILON) &&
      (hlp.posit[Z] < POINT_EPSILON) && (hlp.posit[Z] > -POINT_EPSILON))
    return true;
  else
    return false;
}


bool
operator!=(const GcPoint& p1, const GcPoint& p2)
{
  if (p1 == p2)
    return true;
  else
    return false;
}


bool
operator<(const GcPoint& p1, const GcPoint& p2)
{
  if ((p1.Pos(X) < p2.Pos(X)) && (p1.Pos(Y) < p2.Pos(Y)) &&
      (p1.Pos(Z) < p2.Pos(Z)))
    return true;
  else
    return false;
}


bool
operator>(const GcPoint& p1, const GcPoint& p2)
{
  if ((p1.posit[X] > p2.posit[X]) && (p1.posit[Y] > p2.posit[Y]) &&
      (p1.posit[Z] > p2.posit[Z]))
    return true;
  else
    return false;
}


bool
operator<=(const GcPoint& p1, const GcPoint& p2)
{
  if ((p1.posit[X] <= p2.posit[X]) && (p1.posit[Y] <= p2.posit[Y]) &&
      (p1.posit[Z] <= p2.posit[Z]))
    return 1;
  else
    return 0;
}


bool
operator>=(const GcPoint& p1, const GcPoint& p2)
{
  if ((p1.posit[X] >= p2.posit[X]) && (p1.posit[Y] >= p2.posit[Y]) &&
      (p1.posit[Z] >= p2.posit[Z]))
    return 1;
  else
    return 0;
}
//*end GcPoint methods



//******************************
// Methods for TransMatrix class

TransfMatrix::TransfMatrix()
{
  AA = 1; EE = 1; II = 1;
  BB = 0; CC = 0; DD = 0; FF = 0; GG = 0; HH = 0; JJ = 0; KK = 0; LL = 0;
}

TransfMatrix::TransfMatrix(double a, double b, double c, double d, double e, double f,
                           double g, double h, double i, double j, double k, double l)
{
  AA = a; BB = b; CC = c; DD = d; EE = e; FF = f;
  GG = g; HH = h; II = i; JJ = j; KK = k; LL = l;
}

TransfMatrix::TransfMatrix(TransfMatrix &m)
{
  *this = m;
}

void
TransfMatrix::transform(GcPoint& oldp, GcPoint& newp)
{
  double tmpx = AA * oldp.Pos(X) + DD * oldp.Pos(Y) + GG * oldp.Pos(Z) + JJ;
  double tmpy = BB * oldp.Pos(X) + EE * oldp.Pos(Y) + HH * oldp.Pos(Z) + KK;
  double tmpz = CC * oldp.Pos(X) + FF * oldp.Pos(Y) + II * oldp.Pos(Z) + LL;
  newp = GcPoint(tmpx, tmpy, tmpz);
}

TransfMatrix &
TransfMatrix::operator=(const TransfMatrix& m)
{
  AA = m.AA; BB = m.BB; CC = m.CC;
  DD = m.DD; EE = m.EE; FF = m.FF;
  GG = m.GG; HH = m.HH; II = m.II;
  JJ = m.JJ; KK = m.KK; LL = m.LL;
  return *this;
}

 

//**************************************************************************
//                        Geometry2D methods
//**************************************************************************

Geometry2D::Geometry2D()
{
  nofBoundaryPoints = 0;
  boundaryPoints = NULL;
}


Geometry2D::~Geometry2D()
{
}


// Linearize geometry and store points in boundaryPoints array
//
void
Geometry2D::calcBoundaryPoints()
{
  int i;

  for (i = 0; i < nofBoundaryPoints; i++) {
    delete boundaryPoints[i];
  }
  delete[] boundaryPoints;

  nofBoundaryPoints = 0;
  boundaryPoints = NULL;

  int nof_points;
  GcPoint** points;

  double delta = getParamDeltaU();

  calcLinearizingPoints(nof_points, points, delta);

  model->updateMinimumEdgeSize(nof_points, points);

  if (nof_points == 0) return;

  boundaryPoints = new BoundaryPoint*[nof_points];

  for (i = 0; i < nof_points; i++) {

    BoundaryPoint* bp = new BoundaryPoint;

    boundaryPoints[i] = bp;

    // Check if bp is an existing vertex!
    //BodyElement* v = model->getVertex(points[i]);
    BodyElement* v = model->findVertex(points[i]);

    if ( v != NULL ) {
      bp->tag = v->Tag();
      bp->vertexTag = v->Tag();
    }

    //bp->boundaryTag = boundary_tag;
    bp->point = points[i];
  }

  // Handle possibly closed geometry
  // -------------------------------

  //-Copy start points to end point if U-closed
  //
  if ( isClosedU() ) {
    boundaryPoints[nof_points - 1]->copy(*boundaryPoints[0]);
    boundaryPoints[nof_points - 1]->activeInMeshing = false;
  }

  nofBoundaryPoints = nof_points;

  // NOTE: Delete only array, not the points!!!
  delete[] points;
}


void
Geometry2D::calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u)
{
  nof_points = 0;
  points = NULL;
}


void
Geometry2D::calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u, double delta_v)
{
  nof_points_u = 0;
  nof_points_v = 0;
  points = NULL;

  calcLinearizingPoints(nof_points_u, points, delta_u);
}


void
Geometry2D::getBoundaryPoints(int& count, BoundaryPoint**& points)
{
  points = NULL;
  count = nofBoundaryPoints;

  if ( count > 0 ) {
    points = new BoundaryPoint*[count];

    for (int i = 0; i < count; i++) {
      points[i] = boundaryPoints[i];
    }
  }
}


void
Geometry2D::getDiscretizationData(int nof_components, char*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 0;
  types = NULL;
  valuesU = NULL;
  valuesV = NULL;
  useFixedN = NULL;
}


void
Geometry2D::getMifTags(int& nof_tags, int*& tags)
{
  nof_tags = 1;
  tags = new int[1];
  tags[0] = mifTag;
}


// Output boundary point tags (for Elmer Mesh input file)
//
ostream&
Geometry2D::outputBoundaryPointTags(ostream& out, int indent_size)
{
  int max_per_line = 20;

  //--Normal vertices
  if ( nofBoundaryPoints == 0 ) {
    for (int i = 0; i < nofVertices;  i++) {
      out << vertexTags[i] << " ";
    }

  //--Boundary points
  } else {

    int counter = 0;
    BoundaryPoint* bp;

    for (int i = 0; i < nofBoundaryPoints;  i++) {

      if ( isClosedU() && i == nofBoundaryPoints - 1 ) {
        bp = boundaryPoints[0];
      } else {
        bp = boundaryPoints[i];
      }

      if ( bp != NULL ) {
        out << bp->tag << " ";
        counter++;
      }

      // Continue from next line
      if ( counter >= max_per_line ) {
        out << endl;
        indent(out, indent_size);
        counter = 0;
      }
    }
  }

  return out;
}

// Output boundary points for Elmer Mesh input file (mif-file)
//
ostream&
Geometry2D::output_mif(ostream& out, const char* hdr, bool useFxdMeshN)
{
  int nof_points = 0;

  //--Find nof nodes (verices or linearizing points)
  //
  if ( nofBoundaryPoints > 0 ) {
    nof_points = nofBoundaryPoints;
  } else {
    nof_points = nofVertices;
  }

  // Print positions before and after printing fixed data
  // (for getting the indent size for ids)
  int pos1, pos2;

  pos1 = out.tellp();

  //---Mif-tag and boundary tag
  //
  out << "EdgeId: " << mifTag << " " << hdr;
  
  // If fixed discretation
  if ( useFxdMeshN && this->useFixedMeshN() ) {
    out << " N: " << nof_points - 1;
  }

  //---Nof points and point (vertex or bpoint) tags
  //
  out << "  " << nof_points << "  ";
  pos2 = out.tellp();

  outputBoundaryPointTags(out, pos2 - pos1);

  out << endl;

  return out;
}


//****************************
//***** GcCircle methods *****
//****************************

// Constructor (for egf-files)
//
GcCircle::GcCircle(ecif_ElementComponent_X& tx, IdList& vertex_tags)
{
  Point3** points = NULL;

  nofVertices = tx.nof_vertices;
  vertexTags = NULL;

  if ( nofVertices > 0 ) {

    points = new Point3*[nofVertices];

    vertexTags = new int[nofVertices];

    for (int i = 0; i < nofVertices; i++) {

      BodyElement* v = model->getVertexByTag(tx.vertex_tags[i]);
      vertex_tags.push_back(v->Tag());
      vertexTags[i] = v->Tag();

      points[i] = ((GcPoint*)v->getGeometry())->getPoint();
    }
  }

  if ( !create(nofVertices, points, tx.geometry.edge, tx.lin_delta_type, tx.lin_delta[0]) ){
    geometryOk = false;
    return;
  }

  delete[] points;

  data.useFixedMeshN = (bool)tx.use_fixed_mesh_n;

  init();

  copyMatcValueTable(tx.matcTable, data.matcTable);
  copyMatcValueTable(tx.geometry.edge->matcTable, data.matcTable);
}


// Constructor
// NOTE: no vertices
//
GcCircle::GcCircle(Point3& center ,Point3& start, Point3& stop, bool do_linearize)
{
  nofVertices = 0;
  vertexTags = NULL;

  Point3** points = new Point3*[2];
  points[0] = &start;
  points[1] = &stop;

  ecif_EdgeGeometry_X* param = new ecif_EdgeGeometry_X;
  init_trx_data(*param);
  param->location = new Point3[1];

  param->location[0][0] = center[0];
  param->location[0][1] = center[1];
  param->location[0][2] = center[2];

  if ( !create(2, points, param, LIN_DELTA_NONE, -1.0, do_linearize) ) {
    geometryOk = false;
  }

  delete[] points;

  if (!geometryOk) return;

  init(do_linearize);

  reset_trx_data(*param);
  delete param;
}


// Constructor
// NOTE: no vertices
//
GcCircle::GcCircle(int nof_points, Point3** points, ecif_EdgeGeometry_X* params, bool do_linearize)
{
  int i;

  nofVertices = 0;
  vertexTags = NULL;

  if ( !create(nof_points, points, params, LIN_DELTA_NONE, -1.0, do_linearize) ) {
    geometryOk = false;
  }

  if (!geometryOk) return;

  init(do_linearize);
}


// Constructor
//
GcCircle::GcCircle(BodyElement* vertex1, BodyElement* vertex2, ecif_EdgeGeometry_X* edge, bool do_linearize)
{
  if ( vertex1 == vertex2 ) {
    nofVertices = 1;
  } else {
    nofVertices = 2;
  }

  vertexTags = new int[2];
  vertexTags[0] = vertex1->Tag();
  vertexTags[1] = vertex2->Tag();

  Point3** points = new Point3*[2];
  points[0] = ((GcPoint*)vertex1->getGeometry())->getPoint();
  points[1] = ((GcPoint*)vertex2->getGeometry())->getPoint();
  
  if ( !create(nofVertices, points, edge, LIN_DELTA_NONE, -1.0, do_linearize) ) {
    geometryOk = false;
  }

  delete[] points;

  if (!geometryOk) return;

  init(do_linearize);
}


// Create circle
bool
GcCircle::create(int nof_points, Point3** points, ecif_EdgeGeometry_X* params,
                 enum linDeltaType lin_delta_type, double lin_delta, bool do_linearize)
{
  bool center_given = false;
  bool radius_given = false;

  Point3 p1, p2;

  // Pick parameters
  // ---------------

  //--Circle arc
  if ( nof_points == 2 ) {
    data.start = new GcPoint(points[0]);
    data.end = new GcPoint(points[1]);

    data.start->getPoint(p1);
    data.end->getPoint(p2);

    if ( isZero(dist3(p1, p2)) ) {
      data.isClosed = true;
    } else {
      data.isClosed = false;
    }

  //--Full circle
  } else {
    data.isClosed = true;
    if ( params->start != NULL ) {
      data.start = new GcPoint(params->start);
    }
    if ( params->end != NULL ) {
      data.end = new GcPoint(params->end);
    }
    data.start_u = 0.0;
    data.end_u = TWO_PI;
  }

  // If center given
  if (params->location != NULL) {
    data.centerGiven = true;
    copy3(*params->location, data.center);
  }

  // If radius given
  if (params->radius1 > 0 ) {
    data.radiusGiven = true;
    data.radius = params->radius1;
  }

  // Copy possible defining points
  data.nofDefiningPoints = params->nofDefiningPoints;
  if ( data.nofDefiningPoints  > 0 ) {
    data.definingPoints = new Point3[data.nofDefiningPoints];
    for (int i = 0; i < data.nofDefiningPoints; i++) {
      for (int j = 0; j < 3; j++) {
        data.definingPoints[i][j] = params->definingPoints[i][j];
      }
    }
  }

  data.deltaType = lin_delta_type;
  data.delta = lin_delta;


  // Set final data parameters
  // -------------------------
  if ( !updateData() ) return false;
  
  setDeltaU(lin_delta_type, lin_delta);

  nofBoundaryPoints = 0;
  boundaryPoints = NULL;

  return true;
}


void
GcCircle::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 1;
  types = new enum linDeltaType[1];
  types[0] = data.deltaType;

  valuesU = new double[1];
  valuesU[0] = data.delta;

  valuesV = NULL;

  useFixedN = new bool[1];
  useFixedN[0] = data.useFixedMeshN;
}


bool
GcCircle::getLine(int index, GcLine*& line)
{   
  if ( index < 0 || index >= nofBoundaryPoints - 1) {
    return false;
  }

  Point3* p1 = boundaryPoints[index]->point->getPoint();
  Point3* p2 = boundaryPoints[index + 1]->point->getPoint();
   
  line = new GcLine(p1, p2);

  return true;
}


void
GcCircle::init(bool do_linearize)
{
  if (do_linearize) {
    calcBoundaryPoints();
  }

  boundbox = calcBoundBox();
}


// Destructor.
GcCircle::~GcCircle()
{
  delete[] vertexTags;
}


// Calculates bounding box for the circle.
BoundBox*
GcCircle::calcBoundBox()
{
  RangeVector range;
  BoundBox* bbox = new BoundBox;

  Point3* bp1;
  Point3* bp2;

  if ( nofBoundaryPoints == 0 ) {
    calcBoundaryPoints();
  }

  // Arc
  if ( !data.isClosed ) {

    for (int i = 0; i < nofBoundaryPoints; i++) {
      bbox->extendByPoint(boundaryPoints[i]->point);
    }

  // Full circle
  } else {
    RangeVector rv;
    rv[0] = data.center[0] - data.radius;
    rv[1] = data.center[0] + data.radius;
    rv[2] = data.center[1] - data.radius;
    rv[3] = data.center[1] + data.radius;
    rv[4] = 0;
    rv[5] = 0;
    bbox->extendByRange(rv);
  }

  return bbox;
}


// Calculate linearized boundary
// Store result in the points-buffer
// NOTE: Buffer is allocated here, client should take care of its deletion!!!
//
void
GcCircle::calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u)
{
  nof_points = 0;
  points = NULL;

  PointList all_points;

  // NOTE: Start point is added already here
  //
  all_points.push_back(data.start);

  calcLinearizingPoints(data.start_u, data.end_u, delta_u, all_points);

  // NOTE: End point is added not until here
  //
  if ( data.isClosed ) {
    all_points.push_back(data.start);

  } else {
    all_points.push_back(data.end);
  }

  //---Collect all points to the result
  nof_points = all_points.size();
  points = new GcPoint*[nof_points];

  PointList::iterator itr = all_points.begin();

  for (int i = 0; i < nof_points; i++, itr++) {
    points[i] = *itr;
  }

}


// Helper function 
// 
void
GcCircle::calcLinearizingPoints(double start_u, double end_u, double delta_u, PointList& all_points)
{
  //-If no positive delta_u given, use the default
  if (delta_u <= 0) {
    delta_u = deltaU;
  }

  //-Create points
  double pos_u = start_u + delta_u;

  while (true) {
    
    // At end!
    if ( isEqual(pos_u, end_u) || isGreater(pos_u, end_u)) {
      break;
    }

    double x = data.radius * cos(pos_u) + data.center[0];
    double y = data.radius * sin(pos_u) + data.center[1];
    double z = 0.0;

    all_points.push_back(new GcPoint(x, y, z));

    pos_u += delta_u;
  }
}


void
GcCircle::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  const Point3* p1;
  const Point3* p2;
  int n;

  // A pair of consequtive points as a line
  if ( direction == 1 ) {
    for (n = 0; n < nofBoundaryPoints-1; n++) {
      p1 = boundaryPoints[n]->point->getPoint();
      p2 = boundaryPoints[n+1]->point->getPoint();
      renderer->drawLine(dmode, dstate, 1, p1, p2, elem_id);
    }
  } else {
    for (n = nofBoundaryPoints-1; n > 0;  n--) {
      p1 = boundaryPoints[n]->point->getPoint();
      p2 = boundaryPoints[n-1]->point->getPoint();
      renderer->drawLine(dmode, dstate, 1, p1, p2, elem_id);
    }
  }
}


#if 0
void
GcCircle::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  const Point3* p1;
  const Point3* p2;
  int n;

  // A pair of consequtive points as a line
  if ( direction == 1 ) {
    for (n = 0; n < nofLinearizingPoints-1; n++) {
      p1 = linearizingPoints[n]->getPoint();
      p2 = linearizingPoints[n+1]->getPoint();
      renderer->drawLine(dmode, dstate, 1, p1, p2, elem_id);
    }
  } else {
    for (n = nofLinearizingPoints-1; n > 0;  n--) {
      p1 = linearizingPoints[n]->getPoint();
      p2 = linearizingPoints[n-1]->getPoint();
      renderer->drawLine(dmode, dstate, 1, p1, p2, elem_id);
    }
  }
}
#endif


#if 0
// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns a structure where there is the nof intersection and
// the minimum x-coordinate for the intersections.
RayHit*
GcCircle::isectionsWithXdir(GcPoint* xdir_startp, bool& negative_on_left)
{
  RayHit* result = new RayHit;

  negative_on_left = false;

  result->count = 0;
  result->min_value = 1.0e+100;

  int nof_values;
  double* values;

  isectionsWithXdir(xdir_startp, nof_values, values, negative_on_left);

  result->count = nof_values;

  for (int i = 0; i < nof_values; i++) {

    if ( isLess(values[i], result->min_value) ) {
      result->min_value = values[i];
    }
  }

  if (result->count == 0) {
    delete result;
    result = NULL;
  }

  delete[] values;

  return result;
}
#endif


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns a structure where there is the nof intersection and
// the minimum x-coordinate for the intersections.
RayHit*
GcCircle::isectionsWithXdir(GcPoint* xdir_startp, bool& negative_on_left)
{
  RayHit* result = new RayHit;

  int index = 0;
  GcLine* line;

  while (true) {
    // No more lines!
    if ( !getLine(index++, line) ) {
      break;
    }
    
    bool neg_on_left;
    RayHit* hits = line->isectionsWithXdir(xdir_startp, neg_on_left);

    if ( hits != NULL ) {

      result->count++;
      result->points.push_back(hits->points.front());

      if ( hits->min_value < result->min_value ) {
        result->min_value = hits->min_value;

        // NOTE: Negative on left indicator is picked
        // from the "leftmost" hit, this way we can conclude the
        // possible loop orientation correctly
        negative_on_left = neg_on_left;
      }

    }
    
    delete line;
    delete hits;
  }

  if (result->count == 0) {
    delete result;
    result = NULL;
  }

  return result;
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns nof-intersection and the corresponding x-values in the arguments
void
GcCircle::isectionsWithXdir(GcPoint* xdir_startp, int& nof_values , double*& values)
{
  nof_values = 0;

  // Intersection point(s):
  //  line_x = center[0] + radius * sin(u)
  //  line_y = center[1] + radius * cos(u)

  double px = xdir_startp->Pos(X);
  double py = xdir_startp->Pos(Y);

  // If no chances!
  if ( isGreater(py, data.center[1] + data.radius) ||
       isLess(py, data.center[1] - data.radius)
     ) {
    return;
  }

  values = new double[2];

  // Line-y is constant:
  // cos(u) = (line_y - center[1]) / radius
  double u1 = asin( (py - data.center[1]) / data.radius );
  double u2 = PI - u1; // Other possible solution

  // X values
  // ========
  if ( !isLess(u1, data.start_u) || !isGreater(u1, data.end_u) ) {
    values[nof_values++] = data.center[0] + data.radius * cos(u1);
  }

  if ( !isLess(u2, data.start_u) || !isGreater(u2, data.end_u) ) {
    values[nof_values++] = data.center[0] + data.radius * cos(u2);
  }

  // Make the nearest intersection point the first
  if ( nof_values == 2 ) {

    double d1 = fabs(values[0] - px);
    double d2 = fabs(values[1] - px);

    if ( isLess(d2, d1) ) {
      double tmp = values[0];
      values[0] = values[1];
      values[1] = tmp;
    }
  }
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and y-axis and which starts from point *startp* towards right.
// Returns nof-intersection and the corresponding y-values in the arguments
void
GcCircle::isectionsWithYdir(GcPoint* ydir_startp, int& nof_values , double*& values)
{
  nof_values = 0;

  // Intersection point(s):
  //  line_x = center[0] + radius * sin(u)
  //  line_y = center[1] + radius * cos(u)

  double px = ydir_startp->Pos(X);
  double py = ydir_startp->Pos(Y);

  // If no chances!
  if ( isGreater(px, data.center[0] + data.radius) ||
       isLess(px, data.center[0] - data.radius)
     ) {
    return;
  }

  values = new double[2];

  // Line-x is constant:
  // sin(u) = (line_x - center[0]) / radius
  double u1 = acos( (px - data.center[0]) / data.radius );
  double u2 = -u1; // Other possible solution

  // Y-values
  // ========
  if ( !isLess(u1, data.start_u) || !isGreater(u1, data.end_u) ) {
    values[nof_values++] = data.center[1] + data.radius * sin(u1);
  }

  if ( !isLess(u2, data.start_u) || !isGreater(u2, data.end_u) ) {
    values[nof_values++] = data.center[1] + data.radius * sin(u2);
  }
}


ostream&
GcCircle::output_emf(ostream& out, short indent_size, short indent_level)
{
  const char* def;

  LibFront::output_scalar(out, indent_size, indent_level, EMF_GEOMETRY, NULL, "Circle", true);

  if ( nofVertices > 0 ) {
    LibFront::output_vector(out, indent_size, indent_level, EMF_VERTICES, NULL, nofVertices, vertexTags, false);
  }
 
  def = getMatcString(data.matcTable, EMF_CENTER);

  if ( model->keepMatcDefinitions() && def != NULL ) {
    LibFront::output_matcDef(out, indent_size, indent_level, EMF_CENTER, NULL, def);
  } else {
    LibFront::output_vector(out, indent_size, indent_level, EMF_CENTER, NULL, 3, data.center, false);
  }

  def = getMatcString(data.matcTable, EMF_RADIUS);

  if ( model->keepMatcDefinitions() && def != NULL ) {
    LibFront::output_matcDef(out, indent_size, indent_level, EMF_RADIUS, NULL, def);
  } else {
    LibFront::output_scalar(out, indent_size, indent_level, EMF_RADIUS, NULL, data.radius);
  }
 
  if ( data.deltaType != LIN_DELTA_NONE && data.delta > 0.0 ) {
    switch ( data.deltaType ) {
    case LIN_DELTA_H:
      LibFront::output_scalar(out, indent_size, indent_level, EMF_DELTA_H, NULL, data.delta);
      break;
    case LIN_DELTA_N:
      LibFront::output_scalar(out, indent_size, indent_level, EMF_DELTA_N, NULL, int(data.delta));
      break;
    case LIN_DELTA_U:
      LibFront::output_scalar(out, indent_size, indent_level, EMF_DELTA_U, NULL, data.delta);
      break;
    }
  }

  if ( data.useFixedMeshN ) {
    LibFront::output_scalar(out, indent_size, indent_level, EMF_USE_MESH_N, NULL, 1);
  }

  return out;
}


// Calculate the point corresponding parametric values given as arguments.
// NOTE: Argument parameters are in [0,1], but circle's own parametric
// representation is between [0, 2PI]
void
GcCircle::param2Point(double u_p, double v_p, Point3& point)
{
  double u = data.start_u + u_p * (data.end_u - data.start_u);
  param2Point_own(u, point);
}


// Calculate the point corresponding paramtric values given as arguments.
// NOTE: Argument parameters are in [0,1], but circle's own parametric
// representation is between [0, 2PI]
GcPoint*
GcCircle::param2Point(double u_p, double v_p)
{
  Point3 pp;

  double u = data.start_u + u_p * (data.end_u - data.start_u);

  param2Point_own(u, pp);

  GcPoint* p = new GcPoint(pp);

  return p;
}


// Calculate the point corresponding circle's own paramtric values given as argument.
// NOTE: Argument parameters are in [0,2PI]
bool
GcCircle::param2Point_own(double u, Point3& point)
{
  if (u < data.start_u || u > data.end_u)
  return false;

  point[0] = data.radius * cos(u) + data.center[0];
  point[1] = data.radius * sin(u) + data.center[1];
  point[2] = 0.0;

  return true;
}


// Parameter value for the point (0...1)
void
GcCircle::point2Param(Point3*point, double& u, double& v)
{
  Point3 p;

  if ( isClosedU() ) {
    Point3 sp;
    param2Point(0, 0, sp);
    diff3(*point, sp, p);

  } else {
    diff3(*point, *data.start->getPoint(), p);
  }

  double len = length3(p);

  if ( isZero(len) ) {
    u = 0.0;
  } else {
    // "Center angle" is twice the "arc angle"
    u = 2 * asin(len / (2 * data.radius));
  }

  if ( p[1] < 0 ) {
    u = TWO_PI - u;
  }

  // Normed 0...1
  u = u / (data.end_u - data.start_u);

  v = 0.0;
}


// Calculates the paramtetric values corresponding the argument-point.
ParamPair*
GcCircle::point2Param(GcPoint* point)
{
  double u, v;

  point2Param(point->getPoint(), u, v);

  // Point not in the circle
  if ( u < 0 ) {
    return NULL;

  } else {
    ParamPair* ppair = new ParamPair;
    ppair->u = u;
    ppair->v = v;

    return ppair;
  }
}


// NOTE: Do not update data-structure values in this function
// This should be called after data-structure values are set!!!
//
void
GcCircle::setDeltaU(enum linDeltaType lin_delta_type, double lin_delta)
{
  double arc_angle = data.end_u - data.start_u;
  
  // Max delta-h arc is limited to a half circle
  double max_angle = min(PI, arc_angle);
  double max_h = 2 * data.radius * sin(max_angle / 2);

  double delta_h;
  int delta_n;
  double delta_u;

  // Linearization
  switch (lin_delta_type) {

  case LIN_DELTA_H:
    // Check that given delta-h is not too large
    delta_h = min(max_h, lin_delta);
    delta_u = 2 * asin( delta_h / (data.radius * 2) );
    break;

  case LIN_DELTA_N:
    delta_u = arc_angle / lin_delta;
    break;
  
  case LIN_DELTA_U:
    delta_u = lin_delta * arc_angle;
    break;
  
  case LIN_DELTA_NONE:
  default:
    int n = int(arc_angle / (TWO_PI / 72) + 0.5);
    delta_u = arc_angle / n;
    break;
  }

  // Make an equidistance division
  delta_n = int(0.5 + arc_angle / delta_u);

  // If full circle, we need at least three edges to 
  // linearize it
  if ( data.isClosed ) {
    delta_n = (delta_n < 3)?3:delta_n;
  }

  deltaU = arc_angle / delta_n;

  deltaH = 2 * data.radius * sin(0.5 * deltaU);
}


void
GcCircle::setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN)
{
  if ( nof_components = 0 ) return;

  data.deltaType = types[0];
  data.delta = valuesU[0];
  data.useFixedMeshN = useFixedN[0];
  
  // Update also deltaU which is always in radians
  //
  setDeltaU(data.deltaType, data.delta);
}


bool
GcCircle::updateData()
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm;

  // Calc center if not given
  // ------------------------
  // NOTE: We need 2 or 3 points for this, depending if the radius is
  // given or not!
  //
  if (!data.centerGiven) {

    Point3 *p1, *p2, *p3;
    p1 = p2 = p3 = NULL;

    //--Pick the two points anyway needed
    //
    // For an arc we can use the two vertices
    if (!data.isClosed) {
      p1 = data.start->getPoint();
      p2 = data.end->getPoint();

    // For a full circle we must have at least two defining points
    } else if ( data.nofDefiningPoints >= 2 ) {
      p1 = &(data.definingPoints[0]);
      p2 = &(data.definingPoints[1]);

    // ERROR: No way to calculate the Center!
    } else {
      strm << "Cannot calculate the circle center!" << ends;
      gui->showMsg(strm.str());
      return false;
    }

    //--Next use the radius or a third point
    //
    // If radius given, direction and two point is needed
    if (data.radiusGiven) {
      circle_center(data.radius, *p1, *p2, data.center);

    // Otherwise we need the third point (we pick the last defining point for p3,
    // first two must anyway be given and they were already used!)
    } else if (data.nofDefiningPoints >= 3 ) {
      p3 = &(data.definingPoints[data.nofDefiningPoints - 1]);
      circle_center(*p1, *p2, *p3, data.center);

      // ERROR: No way to calculate the Center!
    } else {
      strm << "Cannot calculate the circle center!" << ends;
      gui->showMsg(strm.str());
      return false;
    }
  }

  // Calc radius if not given
  // ------------------------
  if ( !data.radiusGiven ) {
    
    Point3* p = NULL;

    if ( data.nofDefiningPoints > 0 ) {
      p = &data.definingPoints[0];

    } else if ( data.start != NULL ) {
      p = data.start->getPoint();

    } else if ( data.end != NULL ) {
      p = data.end->getPoint();
    }

    // ERROR: No way to calculate the Radius!
    if ( p == NULL ) {
      strm << "Cannot calculate the circle radius!" << ends;
      gui->showMsg(strm.str());
      return false;
    }

    Point3 len;
    diff3(*p, data.center, len);
    data.radius = length3(len);
  }

  // Set parametrics
  // ---------------
  if ( data.radius > 0 ) {

    //--Circle arc
    if ( !data.isClosed ) {

      // Normalize start
      double start_x = (data.start->Pos(X) - data.center[X]) / data.radius;
      double start_y = (data.start->Pos(Y) - data.center[Y]) / data.radius;
      data.start_u = acos(start_x) * (start_y < 0?-1:1);

      // Normalize end
      double end_x = (data.end->Pos(X) - data.center[X]) / data.radius;
      double end_y = (data.end->Pos(Y) - data.center[Y]) / data.radius;
      data.end_u = acos(end_x) * (end_y < 0?-1:1);

      if (data.end_u < data.start_u) {
        data.end_u += TWO_PI;
      }

    //--Full circle
    } else {

      Point3 p;
      param2Point(0.0, 0.0, p);

      data.start = new GcPoint(p);
      data.end = new GcPoint(p);
    }

  // ERROR: Radius is zero!
  } else {
    strm << "Circle radius is zero!" << ends;
    gui->showMsg(strm.str());
    return false;
  }
  
  return true;
}


// Update geometry. Relevant when Matc-parameters have been changed!
//
bool
GcCircle::updateGeometry(int parent_tag, IdList& vertex_tags)
{
  static char buffer[1025];
  const char* def;
  char* val;

  def = getMatcString(data.matcTable, EMF_CENTER);
  if ( def != NULL ) {
    val = mtc_domath((char*)def);
    if ( val != NULL && val[0] != '\0' && !LibFront::isMatcError(val) ) {
      strstream strm;
      strm << val << ends;
      strm >> data.center[0] >> data.center[1] >> data.center[2];
    }
  }

  def = getMatcString(data.matcTable, EMF_RADIUS);
  if ( def != NULL ) {
    val = mtc_domath((char*)def);
    if ( val != NULL && val[0] != '\0' && !LibFront::isMatcError(val) ) {
      strstream strm;
      strm << val << ends;
      strm >> data.radius;
    }
  }

  def = getMatcString(data.matcTable, EMF_START_POINT);
  if ( def != NULL ) {
    val = mtc_domath((char*)def);
    if ( val != NULL && val[0] != '\0' && !LibFront::isMatcError(val) ) {
      data.start->setPos(val);
    }
  }

  def = getMatcString(data.matcTable, EMF_END_POINT);
  if ( def != NULL ) {
    val = mtc_domath((char*)def);
    if ( val != NULL && val[0] != '\0' && !LibFront::isMatcError(val) ) {
      data.end->setPos(val);
    }
  }

  def = getMatcString(data.matcTable, EMF_MESH_H);
  if ( def != NULL ) {
    val = mtc_domath((char*)def);
    if ( val != NULL && val[0] != '\0' && !LibFront::isMatcError(val) ) {
      strstream strm;
      strm << val << ends;
      strm >> data.delta;
    }
  }

  def = getMatcString(data.matcTable, EMF_MESH_N);
  if ( def != NULL ) {
    val = mtc_domath((char*)def);
    if ( val != NULL && val[0] != '\0' && !LibFront::isMatcError(val) ) {
      strstream strm;
      strm << val << ends;
      strm >> data.delta;
    }
  }
 
  if ( !updateData() ) return false;

  setDeltaU(data.deltaType, data.delta);

  init(true);

  calcBoundaryPoints();

  return true;
}



//**************************
//*****GcLine methods *****
//**************************

GcLine::GcLine(ecif_ElementComponent_X& tx, IdList& vertex_tags)
{
  nofVertices = 2;
  vertexTags = new int[2];

  GcPoint** points = new GcPoint*[2];

  for (int i = 0; i < 2; i++) {

    int vtag = tx.vertex_tags[i];
    vertexTags[i] = vtag;
    vertex_tags.push_back(vtag);

    BodyElement* v = model->getVertexByTag(vtag);
    points[i] = (GcPoint*)v->getGeometry();
  }
  
  data.hasVertexTies = true;
  
  data.useFixedMeshN = (bool)tx.use_fixed_mesh_n;

  create(points[0], points[1]);

  delete[] points;

  copyMatcValueTable(tx.matcTable, data.matcTable);
}


// Constructor
GcLine::GcLine(BodyElement* vertex1, BodyElement* vertex2)
{
  nofVertices = 2;
  vertexTags = new int[2];

  vertexTags[0] = vertex1->Tag();
  vertexTags[1] = vertex2->Tag();

  GcPoint* gp1 = (GcPoint*)vertex1->getGeometry();
  GcPoint* gp2 = (GcPoint*)vertex2->getGeometry();

  data.hasVertexTies = true;

  create(gp1, gp2);
}


// Constructor
// NOTE: no vertices
GcLine::GcLine(Point3* p1, Point3* p2)
{
  GcPoint* gp1 = new GcPoint(p1);
  GcPoint* gp2 = new GcPoint(p2);
  create(gp1, gp2);
}



// Create line
void
GcLine::create(GcPoint* gp1, GcPoint* gp2)
{
  data.start = gp1;
  data.end   = gp2;
  
  updateData();
}


#if 0
// Create line
void
GcLine::create(Point3* p1, Point3* p2)
{
  data.start = new GcPoint(p1);
  data.end   = new GcPoint(p2);

  Point3 len;
  diff3(*(data.end->getPoint()), *(data.start->getPoint()), len);

  data.length = length3(len);
  data.dir = new GcPoint(len);

  data.normal[0] = len[0];
  data.normal[1] = len[1];
  data.normal[2] = len[2];

  normalize(data.normal);
  double tmp = data.normal[0];
  data.normal[0] = -data.normal[1];
  data.normal[1] = tmp;

  boundbox = calcBoundBox();
}
#endif


// Destructor.
GcLine::~GcLine()
{
  //delete[] vertexTags;
}


// Calculates bounding box for the line.
BoundBox*
GcLine::calcBoundBox()
{
  RangeVector range;

  range[0] = data.start->Pos(X);
  range[1] = data.end->Pos(X);
  range[2] = data.start->Pos(Y);
  range[3] = data.end->Pos(Y);
  range[4] = data.start->Pos(Z);
  range[5] = data.end->Pos(Z);

  BoundBox* bbox = new BoundBox(range);

  return bbox;
}


#if 0
void
GcLine::calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u, double delta_v)
{
  PointList all_points;

  calcEdgePoints(data.start_u, data.end_u, all_points);

  //---Add line end-point
  all_points.push_back(data.end);

  //---Collect all points to the result
  nof_points_u = all_points.size();
  nof_points_v = 0;
  points = new GcPoint*[nof_points_u];

  PointList::iterator itr = all_points.begin();

  for (int i = 0; i < nof_points_u; i++, itr++) {
    points[i] = *itr;
  }

}
#endif

#if 0
void
GcLine::calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u, double delta_v)
{
  nof_points_v = 0;

  if ( nofPatterns == 0 ) {
    nof_points_u = 2;
    points = new GcPoint*[2];
    points[0] = data.start;
    points[1] = data.end;
    return;
  }

  // Points with patterns
  // ====================
  PointList all_points;

  all_points.push_back(data.start);

  GcPoint* ebase = data.start;

  calcEdgePointsWithPatterns(this, ebase, nofPatterns, patterns, all_points);

  //---Remove last point from the list if it at the line end point
  GcPoint* lastp = *(--all_points.end());

  if ( ((*data.end) == *lastp) ) {
    all_points.pop_back();
  }

  //---Add line end-point
  all_points.push_back(data.end);

  //---Collect all points to the result
  nof_points_u = all_points.size();
  points = new GcPoint*[nof_points_u];

  PointList::iterator itr = all_points.begin();

  for (int i = 0; i < nof_points_u; i++, itr++) {
    points[i] = *itr;
  }

}
#endif

void
GcLine::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  const Point3* p1;
  const Point3* p2;

  p1 = data.start->getPoint();
  p2 = data.end->getPoint();
  renderer->drawLine(dmode, dstate, direction, p1, p2, elem_id);
}


void
GcLine::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 1;

  types = NULL;
  valuesU = NULL;
  valuesV = NULL;

  useFixedN = new bool[1];
  useFixedN[0] = data.useFixedMeshN;
}



// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns a structure where there is the nof intersection and
// the minimum x-coordinate for the intersections.
// (For a line-segment (as this one) nof intersections is of course 0 or 1)
//
RayHit*
GcLine::isectionsWithXdir(GcPoint* xdir_startp, bool& negative_on_left)
{
  // Parameters t,w are solutions for the following equations:
  // csx + t = clx + klx * w;
  // csy = cly + kly * w;
  // csz = clz + kly * w;

  negative_on_left = false;

  // Constant parameters for *xdir_startp*-line:
  double csx = xdir_startp->Pos(X);
  double csy = xdir_startp->Pos(Y);
  double csz = xdir_startp->Pos(Z);

  // Constant parameters for *this*-line:
  double clx = data.start->Pos(X);
  double cly = data.start->Pos(Y);
  double clz = data.start->Pos(Z);

  // Deflection parameters for *this*-line:
  double klx = data.dir->Pos(X);
  double kly = data.dir->Pos(Y);
  double klz = data.dir->Pos(Z);

  // Lines are parallel (klx !=0 and kly = 0)
  if ( ( klx < -EPSILON || klx > EPSILON) &&
       ( kly > -EPSILON && kly < EPSILON)
     ) {
    return NULL;
  }

  // Intersection-point's domain-parametervalue on *this*-line:
  double w;
  if ( kly != 0 ) {
    w = (csy - cly) / kly;
  } else {
    w = 0.0;
  }

  // If w not in (0,1), intersection-point can't be on *this*-line.
  if (w < -EPSILON || w > 1+EPSILON) return NULL;

  // Intersection-point's domain-parametetrvalue on *xdir_startp*-line.
  // It should be positive because we are looking only to right:
  double t = (clx - csx) + klx * w;

  if (t < -EPSILON ) return NULL;

  RayHit* result = new RayHit;

  Point3* p = new Point3[1];
  (*p)[0] = clx + klx * w;
  (*p)[1] = cly + kly * w;
  (*p)[2] = clz + klz * w;

  result->count = 1;
  result->points.push_back(p);
  result->min_value = csx + t;

  negative_on_left = kly > 0;

  return result;
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp*.
// Returns the intersection x-coordinate value in the argument
// (For a line-segment (as this one) nof intersections is of course 0 or 1)
void
GcLine::isectionsWithXdir(GcPoint* xdir_startp, int& nof_values, double*& values)
{
  // Parameters t,w are solution for the following equations:
  // csx + t = clx + klx * w;
  // csy = cly + kly * w;
  // csz = clz + kly * w;

  nof_values = 0;

  // Constant parameters for *xdir_startp*-line:
  double csx = xdir_startp->Pos(X);
  double csy = xdir_startp->Pos(Y);
  double csz = xdir_startp->Pos(Z);

  // Constant parameters for *this*-line:
  double clx = data.start->Pos(X);
  double cly = data.start->Pos(Y);
  double clz = data.start->Pos(Z);

  // Deflection parameters for *this*-line:
  double klx = data.dir->Pos(X);
  double kly = data.dir->Pos(Y);
  double klz = data.dir->Pos(Z);

  // Lines are parallel (klx !=0 and kly = 0)
  if ( !isZero(klx) && isZero(kly) ) {
    return;
  }

  // Intersection-point's domain-parametervalue on *this*-line:
  double w;

  if ( !isZero(kly) ) {
    w = (csy - cly) / kly;
  } else {
    w = 0.0;
  }

  // If w not in (0,1), intersection-point can't be on *this*-line.
  if ( isLess(w, 0.0) || isGreater(w, 1.0) ) {
    return;
  }

  // Intersection-point's domain-parametetrvalue on *xdir_startp*-line.
  double t = (clx - csx) + klx * w;

  // Set result values
  nof_values = 1;
  values = new double[1];
  values[0] = csx + t;
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and y-axis and which starts from point *startp*.
// Returns the intersection x-coordinate value in the argument
// (For a line-segment (as this one) nof intersections is of course 0 or 1)
void
GcLine::isectionsWithYdir(GcPoint* ydir_startp, int& nof_values, double*& values)
{
  // Parameters t,w are soluition for the following equations:
  // csx + t = clx + klx * w;
  // csy = cly + kly * w;
  // csz = clz + kly * w;

  nof_values = 0;

  // Constant parameters for *xdir_startp*-line:
  double csx = ydir_startp->Pos(X);
  double csy = ydir_startp->Pos(Y);
  double csz = ydir_startp->Pos(Z);

  // Constant parameters for *this*-line:
  double clx = data.start->Pos(X);
  double cly = data.start->Pos(Y);
  double clz = data.start->Pos(Z);

  // Deflection parameters for *this*-line:
  double klx = data.dir->Pos(X);
  double kly = data.dir->Pos(Y);
  double klz = data.dir->Pos(Z);

  // Lines are parallel (klx =0 and kly != 0)
  if ( isZero(klx) && !isZero(kly) ) {
    return;
  }

  // Intersection-point's domain-parametervalue on *this*-line:
  double w;

  if ( !isZero(klx) ) {
    w = (csx - clx) / klx;
  } else {
    w = 0.0;
  }

  // If w not in (0,1), intersection-point can't be on *this*-line.
  if ( isLess(w, 0.0) || isGreater(w, 1.0) ) {
    return;
  }

  // Intersection-point's domain-parametetrvalue on *xdir_startp*-line.
  double t = (cly - csy) + kly * w;

  // Set result values
  nof_values = 1;
  values = new double[1];
  values[0] = csy + t;
}


bool
GcLine::isOnSameAxis(GcPoint& p1, GcPoint& p2)
{
  GcPoint dir_12 = p1 - p2;
  GcPoint dir_1s = p1 - *data.start;

  return ( dir_12.isParallel(data.dir) &&
           dir_1s.isParallel(data.dir) );
}


ostream&
GcLine::output_emf(ostream& out, short indent_size, short indent_level)
{
  LibFront::output_vector(out, indent_size, indent_level, "Vertices", NULL, nofVertices, vertexTags, false);

  if ( data.useFixedMeshN ) {
    LibFront::output_scalar(out, indent_size, indent_level, EMF_USE_MESH_N, NULL, 1);
  }

  return out;
}


// Calculate the point corresponding the paramtric values given as arguments.
void
GcLine::param2Point(double u_p, double v_p, Point3& point)
{
  // For line only parametr *u_p* is relevant.
  point[0] = data.start->Pos(X) + ( u_p * data.dir->Pos(X));
  point[1] = data.start->Pos(Y) + ( u_p * data.dir->Pos(Y));
  point[2] = data.start->Pos(Z) + ( u_p * data.dir->Pos(Z));
}


// Calculate the point corresponding the paramtric values given as arguments.
GcPoint*
GcLine::param2Point(double u_p, double v_p)
{
  // For line only parametr *u_p* is relevant.
  double x = data.start->Pos(X) + ( u_p * data.dir->Pos(X));
  double y = data.start->Pos(Y) + ( u_p * data.dir->Pos(Y));
  double z = data.start->Pos(Z) + ( u_p * data.dir->Pos(Z));

  GcPoint* p = new GcPoint(x,y,z);

  return p;
}


void
GcLine::point2Param(Point3* point, double& u, double& v)
{
  Point3* start = data.start->getPoint();
  Point3* end = data.end->getPoint();

  // Start point
  if ( samepoint(*point, *start) ) {
    u = 0;

  // End point
  } else if ( samepoint(*point, *end) ) {
    u = 1;

  // Point is possibly within line
  } else {

    Point3 len;
    Point3 part;

    diff3(*end, *start, len);
    diff3(*point, *start, part);

    // We return negative value if the point is not
    // geometrically on the line
    Point3 p;
    cross3(len, part, p);

    if ( !isZero(length3(p)) ) {
      u = -1.0;

    // Now we can pick the parameter value from the first nonzero component.
    } else if ( !isZero(len[0]) ) {
      u = part[0] / len[0];

    } else if ( !isZero(len[1]) ) {
      u = part[1] / len[1];

    } else {
      u = 0.0;
    }
  }

  v = 0.0;
}


// Calculates the paramtetric values corresponding the argument-point.
ParamPair*
GcLine::point2Param(GcPoint* p)
{
  double u, v;

  point2Param(p->getPoint(), u, v);

  // Point not in the line!
  if ( u < 0 ) {
    return NULL;

  } else {
    // Create the result-object.
    ParamPair* ppair = new ParamPair;
    ppair->u = u;
    ppair->v = v;

    return ppair;
  }
}


void
GcLine::setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN)
{
  if ( nof_components = 0 ) return;

  data.useFixedMeshN = useFixedN[0];
}

bool
GcLine::updateData()
{
  Point3 len;
  diff3(*(data.end->getPoint()), *(data.start->getPoint()), len);

  data.length = length3(len);
  data.dir = new GcPoint(len);

  data.normal[0] = len[0];
  data.normal[1] = len[1];
  data.normal[2] = len[2];

  normalize(data.normal);
  double tmp = data.normal[0];
  data.normal[0] = -data.normal[1];
  data.normal[1] = tmp;

  boundbox = calcBoundBox();

  return true;
}


// Update geometry. Relevant when Matc-parameters have been changed!
//
bool
GcLine::updateGeometry(int parent_tag, IdList& vertex_tags)
{
  return updateData();
}


//****************************
//***** GcMulti2D methods ****
//****************************

GcMulti2D::GcMulti2D(ecif_Element_X& tx, IdList& vertex_tags)
{
  int i,j;

  data.nofComponents = tx.nof_components;
  data.components = new Geometry2D*[data.nofComponents];

  int final_nof_cmpnts = 0;

  for (i = 0; i < data.nofComponents; i++) {

    data.components[i] = NULL;

    Geometry2D* ptrGmtr = NULL;

    switch ( tx.components[i]->gmtr_type ) {

    case ECIF_CIRCLE:

      //-Function circle
      //
      if ( tx.components[i]->isFunction ) {

        if ( tx.components[i]->isCpp )
          ptrGmtr = createFunctionCircleC(tx.tag, *tx.components[i], vertex_tags);

        else if (tx.components[i]->isF95 )
          ;//ptrGmtr = createFunctionCircleF(tx.tag, *tx.components[i], vertex_tags);

        else if (tx.components[i]->isMatc )
          ptrGmtr = createFunctionCircleM(tx.tag, *tx.components[i], vertex_tags);

      //-Normal circle
      //
      } else {
        ptrGmtr = new GcCircle(*tx.components[i], vertex_tags);
      }
      break;

    case ECIF_LINE:
      ptrGmtr = new GcLine(*tx.components[i], vertex_tags);
      break;

    case ECIF_POLYLINE:

      //-Function polyline
      //
      if ( tx.components[i]->isFunction ) {

        if ( tx.components[i]->isCpp )
          ptrGmtr = createFunctionPolyLineC(tx.tag, *tx.components[i], vertex_tags);

        else if (tx.components[i]->isF95 )
          ;//ptrGmtr = createFunctionPolyLineF(tx.tag, *tx.components[i], vertex_tags);

        else if (tx.components[i]->isMatc )
          ptrGmtr = createFunctionPolyLineM(tx.tag, *tx.components[i], vertex_tags);

      //-Normal polyline
      //
      } else {
        ptrGmtr = new GcPolyLine(*tx.components[i], vertex_tags);
      }
      break;

    case ECIF_NURBS:
      ptrGmtr = new GcNurbsCurve(*tx.components[i], vertex_tags);
      break;

    default:
      ptrGmtr = NULL;
      break;
    }

    // Error in geometry!
    //
    if ( ptrGmtr == NULL || !ptrGmtr->geometryIsOk() ) {
      delete ptrGmtr;
      geometryOk = false;
      return;
    }

    data.components[i] = ptrGmtr;
    
    int ccount = ptrGmtr->getNofComponents();
    final_nof_cmpnts += ccount;
  }

  // Init other data
  init();

  // Collect all (simple) components into one array
  //
  Geometry2D** tmp = new Geometry2D*[final_nof_cmpnts];

  data.isDeletable = new bool[final_nof_cmpnts];
  
  int idx = 0;
  for (i = 0; i < data.nofComponents; i++) {

    Geometry2D* pg = data.components[i];

    // If a component is simple, take as it is
    //
    if ( !pg->isMultiGeometry() ) {
      tmp[idx] = pg;
      data.isDeletable[idx] = false;
      idx++;

    // If a component is a Multi2D, copy its subcomponents!
    // NOTE: Sub components should have been already unpurged
    // recursively!
    //
    } else {
      for (j = 0; j < pg->getNofComponents(); j++) {
        tmp[idx] = (Geometry2D*)pg->getComponent(j);
        data.isDeletable[idx] = true;
        idx++;
      }
    }
  }

  // Store original components as components to be output to
  // emf file
  //
  data.nofEmfComponents = data.nofComponents;
  data.emfComponents = data.components;

  // These are components which all are geometry primitives (polylines,circles etc)
  data.nofComponents = final_nof_cmpnts;
  data.components = tmp;

  // Finally linearize and calc boundbox
  //
  updateData();
}


// Destructor.
GcMulti2D::~GcMulti2D()
{
  int i;

  for (i = 0; i < data.nofComponents; i++) {
    delete data.components[i];
    data.components[i] = NULL;
  }

  for (i = 0; i < data.nofEmfComponents; i++) {
    //delete data.emfComponents[i];
    //data.emfComponents[i] = NULL;
  }
}


// Calculates bounding box for the line.
BoundBox*
GcMulti2D::calcBoundBox()
{
  return data.components[0]->calcBoundBox();
}


// Calculate linearized boundary
void
GcMulti2D::calcBoundaryPoints()
{
  for (int i = 0; i < data.nofComponents; i++) {
    data.components[i]->calcBoundaryPoints();
  }
}


// Calculate linearized boundary
// Store result in the points-buffer
// NOTE: Buffer is allocated here, client should take care of its
// deletion!!!
void
GcMulti2D::calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u)
{
  for (int i = 0; i < data.nofComponents; i++) {
    data.components[i]->calcLinearizingPoints(nof_points, points, delta_u);
  }
}



// ================================================================
//                  C++-function calls
// ================================================================

//===========================
// Create egf function-CIRCLE
// ==========================
// C++-version
//
GcMulti2D*
GcMulti2D::createFunctionCircleC(int edge_tag, ecif_ElementComponent_X& txc,
                                IdList& vertex_tags)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm, strm1, strm2;
  int i,j;

  bool is_ok = true;


  Hdll hDLL = NULL;   // Handle to DLL library
  Hfunc hFunc = NULL; // Handle to dll function
  char* err_msg = NULL;
    
  is_ok = loadDllFunction(txc.libraryName, txc.functionName, hDLL, hFunc, err_msg);

  int nof_circles = 0;
  egf_Circle* circles = NULL;
  dllCircleFunc userFunc;

  if ( is_ok ) {
    userFunc = (dllCircleFunc)hFunc;
  } else {
    strm2 << err_msg << endl;
  }

  // Ok, try to call the circle fucntion
  if (is_ok ) {
    int rc = userFunc(EGF_CIRCLE, txc.argc, txc.argv, txc.startPoint, txc.endPoint, nof_circles, circles);
  
    // Error: function cal failed
    if ( rc != 0 ) {
      strm2 << "Function call failed with return code " << rc;
      is_ok = false;
    }
  }

  if ( is_ok && nof_circles == 0 ) {
    is_ok = false;
    strm2 << "Nof circles is zero!";
  }

  // Error
  //
  if ( !is_ok ) {
    strm1 << "***ERROR in C++ circle function for edge " << edge_tag << ends;
    strm2 << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str());

    if ( hDLL != NULL ) {
      closeDllLibrary(hDLL);
    }

    return NULL;
  }

  GcMulti2D* ptrG = createFunctionCircle(nof_circles, circles, txc, vertex_tags);
  
  ptrG->data.funcType = ECIF_CPP;

  if ( hDLL != NULL ) {
    closeDllLibrary(hDLL);
  }
  
  if ( ptrG == NULL ) {
    strm1 << "***ERROR in Front when processing circle function for edge " << edge_tag << ends;
    gui->showMsg(strm1.str());
    return NULL;
  }

  return ptrG;

}


// ============================
// Create egf function-POLYLINE
// ============================
// C++-version
//
GcMulti2D*
GcMulti2D::createFunctionPolyLineC(int edge_tag, ecif_ElementComponent_X& txc,
                                   IdList& vertex_tags)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm, strm1, strm2;
  int i,j;

  Hdll hDLL = NULL;   // Handle to DLL library
  Hfunc hFunc = NULL; // Handle to dll function
  char* err_msg = NULL;
    
  bool is_ok = loadDllFunction(txc.libraryName, txc.functionName, hDLL, hFunc, err_msg);

  int nof_polylines = 0;
  egf_PolyLine* polylines = NULL;
  dllPolyLineFunc userFunc;

  if ( is_ok ) {
    userFunc = (dllPolyLineFunc)hFunc;
  } else {
    strm2 << err_msg << endl;
  }

  if ( is_ok ) {
    int rc = userFunc(EGF_POLYLINE, txc.argc, txc.argv, txc.startPoint, txc.endPoint, nof_polylines, polylines);
    
    // Error: function call failed
    if ( rc != 0 ) {
      strm2 << "Function call failed with return code " << rc;
      is_ok = false;
    }
  }

  if ( is_ok && nof_polylines == 0 ) {
    is_ok = false;
    strm2 << "Nof polylines is zero!";
  }

  // Error in library or function call, stop
  // =======================================
  if ( !is_ok ) {
    strm1 << "***ERROR in C++ polyline function for edge " << edge_tag << ends;
    strm2 << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str());

    if ( hDLL != NULL ) {
      closeDllLibrary(hDLL);
    }

    return NULL;
  }


  GcMulti2D* ptrG = createFunctionPolyLine(nof_polylines, polylines, txc, vertex_tags);

  ptrG->data.funcType = ECIF_CPP;

  if ( hDLL != NULL ) {
    closeDllLibrary(hDLL);
  }

  if ( ptrG == NULL ) {
    strm1 << "***ERROR in Front when processing polyline function for edge " << edge_tag << ends;
    gui->showMsg(strm1.str());
    return NULL;
  }

  return ptrG;
}


#if 0
// ================================================================
//                  Fortan95-function calls
// ================================================================

// ==========================
// Create egf function-CIRCLE
// ==========================
// F95-version
//
GcMulti2D*
GcMulti2D::createFunctionCircleF(int edge_tag, ecif_ElementComponent_X& txc,
                                IdList& vertex_tags)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm, strm1, strm2;
  int i,j;

  bool is_ok = true;

  Hdll hDLL; // Handle to DLL

  int nof_circles = 0;
  egf_Circle* circles = NULL;

  // Load library
  hDLL = LoadLibrary(txc.libraryName);
  
  // Error, library not found
  //
  if (hDLL == NULL) {
    is_ok = false;
    strm2 << "Cannot open dll-library: " << txc.libraryName;
  }

  // Typedef for a safer function call
  typedef int (Callbackp dllFunc)(int argc, double* argv,
                                  egf_Point3* start, egf_Point3* end,
                                  int& nof_circles, egf_Circle*& circles);

  dllFunc userFunc;   // Function pointer to the user's function

  // Try first decorated name (_name@n)
  strm << "_" << txc.functionName << "@24" << ends;
  userFunc = (dllFunc)GetProcAddress(hDLL, strm.str());

  // Next try undecorated name
  if ( userFunc == NULL ) {
    userFunc = (dllFunc)GetProcAddress(hDLL, txc.functionName);
  }

  // 
  if (is_ok ) {
    
    // Error: function not found
    if ( userFunc == NULL ) {
      is_ok = false;
      strm2 << "Cannot open function " << txc.functionName;

    // Ok, try to call the function
    } else {
      int rc = userFunc(txc.argc, txc.argv, txc.startPoint, txc.endPoint, nof_circles, circles);
    
      // Error: function cal failed
      if ( rc != 0 ) {
        strm2 << "Function call failed with return code " << rc;
        is_ok = false;
      }
    }
  }

  if ( is_ok && nof_circles == 0 ) {
    is_ok = false;
    strm2 << "Nof circles is zero!";
  }

  // Error
  //
  if ( !is_ok ) {
    strm1 << "***ERROR in F95 circle function for edge " << edge_tag << ends;
    strm2 << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str());

    if ( hDLL != NULL ) {
      closeDllLibrary(hDLL);
    }

    return NULL;
  }


  GcMulti2D* ptrG = createFunctionCircle(nof_circles, circles, txc, vertex_tags);

  ptrG->data.funcType = ECIF_F95;

  if ( hDLL != NULL ) {
    closeDllLibrary(hDLL);
  }

  if ( ptrG == NULL ) {
    strm1 << "***ERROR in Front when processing polyline function for edge " << edge_tag << ends;
    gui->showMsg(strm1.str());
    return NULL;
  }

  return ptrG;
}


// ============================
// Create egf function-POLYLINE
// ============================
// F95-version
//
GcMulti2D*
GcMulti2D::createFunctionPolyLineF(int edge_tag, ecif_ElementComponent_X& txc,
                                IdList& vertex_tags)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm, strm1, strm2;
  int i,j;

  Hdll hDLL = NULL; // Handle to DLL

  // Typedef for a safer function call
  typedef int (Callbackp dllFunc)(int argc, double* argv,
                                  egf_Point3* start, egf_Point3* end,
                                  int& nof_polylines, egf_PolyLine*& polylines);
  dllFunc userFunc = NULL;  // Function pointer to the user's function

  int nof_polylines = 0;
  egf_PolyLine* polylines = NULL;

  bool is_ok = true;

  // Load the library
  hDLL = LoadLibrary(txc.libraryName);

  // Error, library not found
  if (hDLL == NULL) {
    is_ok = false;
    strm2 << "Cannot open dll-library: " << txc.libraryName;
  }
  
  if ( is_ok ) {
    // Try first  decorated name (_name@n)
    strm << "_" << txc.functionName << "@24" << ends;
    userFunc = (dllFunc)GetProcAddress(hDLL, strm.str());

    // Next try undecorated name
    if ( userFunc == NULL ) {
      userFunc = (dllFunc)GetProcAddress(hDLL, txc.functionName);
    }
  }
 
  if ( is_ok ) {

    // Error: function not found
    if ( userFunc == NULL ) {
    is_ok = false;
    strm2 << "Cannot open function " << txc.functionName;

    // Ok, try to call the function
    } else {
      
      int rc = userFunc(txc.argc, txc.argv, txc.startPoint, txc.endPoint, nof_polylines, polylines);
      
      // Error: function call failed
      if ( rc != 0 ) {
        strm2 << "Function call failed with return code " << rc;
        is_ok = false;
      }
    }
  }

  if ( is_ok && nof_polylines == 0 ) {
    is_ok = false;
    strm2 << "Nof polylines is zero!";
  }

  // Error in library or function call, stop
  // =======================================
  if ( !is_ok ) {
    strm1 << "***ERROR in F95 polyline function for edge " << edge_tag << ends;
    strm2 << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str());

    if ( hDLL != NULL ) {
      closeDllLibrary(hDLL);
    }

    return NULL;
  }


  GcMulti2D* ptrG = createFunctionPolyLine(nof_polylines, polylines, txc, vertex_tags);

  ptrG->data.funcType = ECIF_F95;

  if ( hDLL != NULL ) {
    closeDllLibrary(hDLL);
  }


  if ( ptrG == NULL ) {
    strm1 << "***ERROR in Front when processing polyline function for edge " << edge_tag << ends;
    gui->showMsg(strm1.str());
    return NULL;
  }

  return ptrG;
}
#endif


// ================================================================
//                    Matc-function calls
// ================================================================

// ==========================
// Create egf function-CIRCLE
// ==========================
// Matc-version
//
GcMulti2D*
GcMulti2D::createFunctionCircleM(int edge_tag, ecif_ElementComponent_X& txc,
                                IdList& vertex_tags)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm, strm1, strm2;

  char* matc_result = NULL;
 
  int i,j;

  LibFront::initMatcFrmt();  

  bool is_ok = true;

  int nof_circles = 0;
  egf_Circle* circles = NULL;

  // Check that matc function exists (call help("name")
  //
  strm << "help(\"";
  // Copy function name until first left parnethesis
  for (i = 0; i < strlen(txc.functionName) ; i++) {
    char c = txc.functionName[i];
    if ( c != '(' ) {
      strm << c;
    } else {
      break;
    }
  }
  strm << "\")" << ends;

  char* fnm = strm.str();

  if ( LibFront::isMatcError(mtc_domath(fnm)) ) {
    is_ok = false;
    strm2 << "MATC function: " << txc.functionName << " not defined!";
  }

  // Ok, try to call the function
  if (is_ok ) {
    matc_result = mtc_domath(txc.functionName);
    //xc.argc, txc.argv, txc.startPoint, txc.endPoint, nof_circles, circles);
    
    // Error: function call failed
    if ( matc_result == NULL || matc_result[0] == '\0' ) {
      strm2 << "***MATC ERROR: no value from circle function call: " << txc.functionName;
      is_ok = false;
    } else if ( LibFront::isMatcError(matc_result) ) {
      LibFront::formatMatcError(matc_result);
      strm2 << matc_result;
      is_ok = false;
    } else {
      is_ok = getMatcCircleData(matc_result, nof_circles, circles, strm2);
    }
  }

  if ( is_ok && nof_circles == 0 ) {
    is_ok = false;
    strm2 << "Nof circles is zero!";
  }

  // Error
  //
  if ( !is_ok ) {
    strm1 << "***ERROR in Matc circle function for edge " << edge_tag << ends;
    strm2 << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str());
    return NULL;
  }
  
  GcMulti2D* ptrG = createFunctionCircle(nof_circles, circles, txc, vertex_tags);

  ptrG->data.funcType = ECIF_MATC;

  if ( ptrG == NULL ) {
    strm1 << "***ERROR in Front when processing circle function for edge " << edge_tag << ends;
    gui->showMsg(strm1.str());
    return NULL;
  }

  return ptrG;
}


// ============================
// Create egf function-polyline
// ============================
// Matc-version
//
GcMulti2D*
GcMulti2D::createFunctionPolyLineM(int edge_tag, ecif_ElementComponent_X& txc,
                                IdList& vertex_tags)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm, strm1, strm2;
  int i,j;

  char* matc_result = NULL;

  LibFront::initMatcFrmt();  

  int nof_polylines = 0;
  egf_PolyLine* polylines = NULL;

  bool is_ok = true;

  // Check that matc function exists (call help("name")
  //
  strm << "help(\"";
  // Copy function name until first left parnethesis
  for (i = 0; i < strlen(txc.functionName); i++) {
    char c = txc.functionName[i];
    if ( c != '(' ) {
      strm << c;
    } else {
      break;
    }
  }
  strm << "\")" << ends;

  char* fnm = strm.str();

  if ( LibFront::isMatcError(mtc_domath(fnm)) ) {
    is_ok = false;
    strm2 << "MATC function: " << txc.functionName << " not defined!";
  }

 
  // Ok, try to call the polyline function
  if (is_ok ) {
    matc_result = mtc_domath(txc.functionName);

    LibFront::trim(matc_result);

    //int rc = userFunc(txc.argc, txc.argv, txc.startPoint, txc.endPoint, nof_polylines, polylines);
    
    // Error: function call failed
    if ( matc_result == NULL || matc_result[0] == '\0' ) {
      strm2 << "***MATC ERROR: no value from polyline function call: " << txc.functionName;
      is_ok = false;
    } else if ( LibFront::isMatcError(matc_result) ) {
      strm2 << matc_result;
      is_ok = false;
    } else {
      is_ok = getMatcPolyLineData(matc_result, nof_polylines, polylines, strm2);
    }
  }

  if ( is_ok && nof_polylines == 0 ) {
    is_ok = false;
    strm2 << "Nof polylines is zero!";
  }

  // Error in library or function call, stop
  // =======================================
  if ( !is_ok ) {
    strm1 << "***ERROR in polyline function for edge " << edge_tag << ends;
    strm2 << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str());

    return NULL;
  }

  GcMulti2D* ptrG = createFunctionPolyLine(nof_polylines, polylines, txc, vertex_tags);

  ptrG->data.funcType = ECIF_MATC;

  if ( ptrG == NULL ) {
    strm1 << "***ERROR in Front when processing polyline function for edge " << edge_tag << ends;
    gui->showMsg(strm1.str());
    return NULL;
  }

  return ptrG;
}


// Matc circle-function helper
//
bool
GcMulti2D::getMatcCircleData(const char* matc_result,
                             int& nof_circles, egf_Circle*& circles,
                             strstream& msg_strm)
{
  int i,k;

  nof_circles = 0;
  circles = NULL;
  
  int dim;
  int radius_given;
  int center_given;
  int start_given;
  int end_given;
  int vflags_given;
  
  strstream strm;
  strm << matc_result << ends;
 
  // These must always be in the matc-result
  strm >> nof_circles >> dim;
  strm >> radius_given >> center_given >> start_given >> end_given >> vflags_given;

  // Not enough common parameter data!
  //
  int st = strm.fail();

  if ( 0 != strm.fail() ) {
    msg_strm << "Not enough parameter values for a Matc function circle";
    return false;
  }

  if ( nof_circles == 0 ) return false;

  if ( dim < 1 || dim > 3 ) {
    msg_strm << "Illegal point dimension value for a Matc function circle (" << dim << ")";
    return false;
  }

  circles = new egf_Circle[nof_circles];

  for (i = 0; i < nof_circles; i++) {
 
    // Not enough data!
    //
    if ( 0 != strm.eof() ) {
      delete[] circles;
      circles = NULL;
      nof_circles = 0;
      msg_strm << "Not enough data values for a Matc function circle";
      return false;
    }
    
    egf_Circle* c = &(circles[i]);

    if ( 0 != radius_given ) {
      c->radiusGiven = true;
      strm >> c->radius;
    }

    if ( 0 != center_given ) {
      c->centerGiven = true;
      for (k = 0; k < dim; k++) {
        strm >> c->center[k];
      }
    }

    if ( 0 != start_given ) {
      c->startGiven = true;
      for (k = 0; k < dim; k++) {
        strm >> c->start[k];
      }
    }

    if ( 0 != end_given ) {
      c->endGiven = true;
      for (k = 0; k < dim; k++) {
        strm >> c->end[k];
      }
    }

    if ( 0 != vflags_given ) {
      strm >> c->vflags[0] >> c->vflags[1];
    }
 
    // Not enough data!
    //
    if ( 0 != strm.fail() ) {
      delete[] circles;
      circles = NULL;
      nof_circles = 0;
      msg_strm << "Not enough data values for a Matc function circle";
      return false;
    }
 
  }

  return true;
}
 

// Matc polyline-function helper
//
bool
GcMulti2D::getMatcPolyLineData(const char* matc_result,
                            int& nof_polylines, egf_PolyLine*& polylines,
                            strstream& msg_strm)
{
  int i,k, nof_points;

  int dim;
  int vflags_given;

  nof_polylines = 0;
  polylines = NULL;

  strstream strm;
  strm << matc_result << ends;

  // These must always be in the matc-result
  strm >> nof_polylines >> dim >> vflags_given;

  // Not enough common parameter data!
  //
  if ( 0 != strm.fail() ) {
    msg_strm << "Not enough parameter values for Matc a function polyline";
    return false;
  }

  if ( nof_polylines == 0 ) return false;

  if ( dim < 1 || dim > 3 ) {
    msg_strm << "Illegal point dimension value for a Matc function polyline (" << dim << ")";
    return false;
  }

  polylines = new egf_PolyLine[nof_polylines];
  
  // Read all polylines
  //
  for (i = 0; i < nof_polylines; i++) {

    // Not enough data!
    //
    if ( strm.eof() ) {
      delete[] polylines;
      polylines = NULL;
      nof_polylines = 0;
      msg_strm << "Not enough data values for a Matc function polyline";
      return false;
    }
    
    egf_PolyLine* c = &(polylines[i]);

    strm >> nof_points;
    
    c->points = new Point3[nof_points];
    
    if ( 0 != vflags_given ) {
      c->vflags = new int[nof_points];
    }
      
    // Read all polyline points and possible vertex flags
    //
    for (int n = 0; n < nof_points; n++) {
      for (k = 0; k < dim; k++) {
        strm >> c->points[n][k];
      }
      if ( 0 != vflags_given ) {
        strm >> c->vflags[n];
      }
    }

    if ( 0 != strm.fail() ) {
      delete[] polylines;
      polylines = NULL;
      nof_polylines = 0;
      msg_strm << "Not enough data values for a Matc function polyline";
      return false;
    }
  }

  return true;
}



// ================================================================
//                       Common calls
// ================================================================


// =============================
// Common create function-circle 
// =============================
//
GcMulti2D*
GcMulti2D::createFunctionCircle(int nof_circles, egf_Circle* circles,
                                ecif_ElementComponent_X& txc,
                                IdList& vertex_tags)
{

// Try to create a new circle element
// ----------------------------------
try
{
  ecif_Element_X ntx;
  init_trx_data(ntx);

  ntx.tplg_type = ECIF_EDGE;
  ntx.nof_components = nof_circles;
  ntx.components = new ecif_ElementComponent_X*[nof_circles];

  GcPoint point;
  int i,k;

  for (i = 0; i < nof_circles; i++) {
    
    // Create new sub-circle
    ecif_ElementComponent_X* ntxc = new ecif_ElementComponent_X;
    init_trx_data(*ntxc);
    ntx.components[i] = ntxc;

    ntxc->gmtr_type = ECIF_CIRCLE;

    ntxc->lin_delta_type = txc.lin_delta_type;
    ntxc->lin_delta[0] = txc.lin_delta[0];
    ntxc->lin_delta[1] = txc.lin_delta[1];

    ntxc->geometry.edge = new ecif_EdgeGeometry_X;
    init_trx_data(*(ntxc->geometry.edge));

    // Arguments
    egf_Circle* c = &(circles[i]);

    ntxc->geometry.edge->radius1 = c->radius;
    
    if ( c->centerGiven ) {
      ntxc->geometry.edge->location = new Point3[1];
      for (k = 0; k < 3; k++) {
        (*ntxc->geometry.edge->location)[k] = c->center[k];
      }
    }

    if ( c->startGiven ) {
      ntxc->geometry.edge->start = new Point3[1];
      for (k = 0; k < 3; k++) {
        (*ntxc->geometry.edge->start)[k] = c->start[k];
      }
    }

    if ( c->endGiven ) {
      ntxc->geometry.edge->end = new Point3[1];
      for (k = 0; k < 3; k++) {
        (*ntxc->geometry.edge->end)[k] = c->end[k];
      }
    }

    if ( !c->startGiven ||
         !c->endGiven   ||
         isZero(dist3(c->start, c->end))
       ) {
      ntxc->geometry.edge->isClosed = true;
    }

    // Store vertex tags if start/end points flagged as vertex
    ntxc->vertex_tags = new int[2];
    ntxc->vertex_tags[0] = NO_INDEX;
    ntxc->vertex_tags[1] = NO_INDEX;
    //
    //-Start point
    if ( c->startGiven && 0 != c->vflags[0] ) {
      point.setPos(c->start[0], c->start[1], c->start[2]);
      BodyElement* v = model->findVertex(&point);
      //-Create new vertex if needed
      if ( v == NULL ) {
        v = new BodyElement1D(&point);
        model->addBodyElement(v);
      }
      ntxc->vertex_tags[0] = v->Tag();
    }
    //-End point
    if ( c->endGiven && 0 != c->vflags[1] ) {
      point.setPos(c->end[0], c->end[1], c->end[2]);
      BodyElement* v = model->findVertex(&point);
      //-Create new vertex if needed
      if ( v == NULL ) {
        v = new BodyElement1D(&point);
        model->addBodyElement(v);
      }
      ntxc->vertex_tags[1] = v->Tag();
    }

  } // all circles

  // Ok, create new Multi2D geometry from circles
  GcMulti2D* ptrG = new GcMulti2D(ntx, vertex_tags);
  
  if ( !ptrG->geometryIsOk() ) {
    delete ptrG;
    return NULL;
  }
  
  ptrG->data.delta = txc.lin_delta[0];
  ptrG->data.deltaType = txc.lin_delta_type;
  ptrG->data.useFixedMeshN = txc.use_fixed_mesh_n;

  ptrG->data.gmtrType = ECIF_CIRCLE;
  ptrG->data.functionName = NULL;
  ptrG->data.libraryName = NULL;
  update_dyna_string(ptrG->data.functionName, txc.functionName);
  update_dyna_string(ptrG->data.libraryName, txc.libraryName);

  ptrG->data.funcArgC = txc.argc;
  ptrG->data.funcArgV = new double[txc.argc];
  for (i = 0; i < txc.argc; i++) {
    ptrG->data.funcArgV[i] = txc.argv[i];
  }
  
  copyMatcValueTable(txc.matcTable, ptrG->data.matcTable);

  reset_trx_data(ntx);

  return ptrG;
}


// Error in circle element
// -----------------------
catch (...)
{
  return NULL;
}

} // End createFunctionCircle()



// ===============================
// Common create function-polyline 
// ===============================
//
GcMulti2D*
GcMulti2D::createFunctionPolyLine(int nof_polylines, egf_PolyLine* polylines,
                                  ecif_ElementComponent_X& txc,
                                  IdList& vertex_tags)
{
  int i,j;

// Try to create a new polyline element
// ---------------------------------
try
{
  ecif_Element_X ntx;
  init_trx_data(ntx);

  ntx.tplg_type = ECIF_EDGE;
  ntx.nof_components = nof_polylines;
  ntx.components = new ecif_ElementComponent_X*[nof_polylines];

  GcPoint point;

  for (i = 0; i < nof_polylines; i++) {
    
    // Create new sub-polyline
    ecif_ElementComponent_X* ntxc = new ecif_ElementComponent_X;
    init_trx_data(*ntxc);
    ntx.components[i] = ntxc;

    egf_PolyLine* c = &(polylines[i]);
    int nof_points = c->nof_points;

    ntxc->gmtr_type = ECIF_POLYLINE;

    ntxc->lin_delta_type = txc.lin_delta_type;
    ntxc->lin_delta[0] = txc.lin_delta[0];
    ntxc->lin_delta[1] = txc.lin_delta[1];

    ntxc->geometry.edge = new ecif_EdgeGeometry_X;
    init_trx_data(*(ntxc->geometry.edge));
    
    if ( isZero(dist3(c->points[0], c->points[nof_points-1])) ) {
      ntxc->geometry.edge->isClosed = true;
    }
    ntxc->geometry.edge->nofDefiningPoints = nof_points;
    ntxc->geometry.edge->definingPoints = new Point3[nof_points];
    ntxc->geometry.edge->pointVertexFlags = new bool[nof_points];

    // Count nof vertices
    int nof_vertices = 0;

    for (j = 0; j < nof_points; j++) {
      if ( c->vflags[i] ) nof_vertices++;
    }

    ntxc->nof_vertices = nof_vertices;
    ntxc->vertex_tags = new int[nof_vertices];
    
    // Loop all polyline points
    //
    for (j = 0; j < nof_points; j++) {
      
      // If point flagged as vertex, store coordinates
      // also in the point-object (for finding the vertex)
      for (short k = 0; k < 3; k++) {

        double p = c->points[j][k];

        if ( 1 == c->vflags[i] ) {
          point.setPosit(k, p);
        }

        ntxc->geometry.edge->definingPoints[j][k] = p;
      }
      
      // Store vertex tag if point is flagged as vertex
      //
      if ( 1 == c->vflags[j] ) {
        ntxc->geometry.edge->pointVertexFlags[j] = true;
        BodyElement* v = model->findVertex(&point);

        //-Create new vertex if needed
        if ( v == NULL ) {
          v = new BodyElement1D(&point);
          model->addBodyElement(v);
        }
        ntxc->vertex_tags[j] = v->Tag();

      } else {
        ntxc->geometry.edge->pointVertexFlags[j] = false;
      }

    } // all polyline points


  } // all polylines

  // Create new Multi2D elements from polylines

  GcMulti2D* ptrG = new GcMulti2D(ntx, vertex_tags);
  
  ptrG->data.gmtrType = ECIF_POLYLINE;
  ptrG->data.functionName = NULL;
  ptrG->data.libraryName = NULL;
  update_dyna_string(ptrG->data.functionName, txc.functionName);
  update_dyna_string(ptrG->data.libraryName, txc.libraryName);
  
  ptrG->data.funcArgC = txc.argc;
  ptrG->data.funcArgV = new double[txc.argc];
  for (i = 0; i < txc.argc; i++) {
    ptrG->data.funcArgV[i] = txc.argv[i];
  }

  copyMatcValueTable(txc.matcTable, ptrG->data.matcTable);

  reset_trx_data(ntx);

  return ptrG;
}

// Error in polyline element
// ----------------------
catch (...)
{
  return NULL;
}

} // End createFunctionPolyLine()



void
GcMulti2D::draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  if ( gmtr_index < 0 || gmtr_index > data.nofComponents - 1 ) return;

  Geometry* pg = data.components[gmtr_index];
  
  if ( 1 == pg->getNofComponents() ) {
    pg->draw(renderer, dmode, dstate, direction, elem_id);
  } else {
    pg->draw(gmtr_index, renderer, dmode, dstate, direction, elem_id);
  }
}


void
GcMulti2D::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id, bool is_first_loop)
{
  for (int i = 0; i < data.nofComponents; i++) {
    data.components[i]->draw(renderer, dmode, dstate, direction, elem_id);

    if ( i < data.nofComponents - 1 ) {
      renderer->stopDrawingCadBodyElementLoop();
      renderer->startDrawingCadBodyElementLoop(is_first_loop);
    }
  }
}


void
GcMulti2D::getBoundaryPoints(int& count, BoundaryPoint**& points)
{
  count = getNofBoundaryPoints();

  points = new BoundaryPoint*[count];

  int index = 0;

  for (int i = 0; i < data.nofComponents; i++) {
    int countc;
    BoundaryPoint** pointsc;
    data.components[i]->getBoundaryPoints(countc, pointsc);

    for (int j = 0; j < countc; j++) {
      points[index++] = pointsc[j];
    }

    delete[] pointsc;
    pointsc = NULL;
  }
}


// NOTE: Argument flags tells if only one of the emf components (ie. components defined
// egf- and emf-files) or if the one of the all primitive components should be given
//
// NOTE: The component can be different is some the emf components are functions which create multiple
// primitive components (circles,polylines)
//
Geometry*
GcMulti2D::getComponent(int index, bool only_emf_components)
{
  if ( only_emf_components ) {
    if ( index < 0 || index > data.nofEmfComponents - 1 ) return NULL;
    return data.emfComponents[index];

  } else {
    if ( index < 0 || index > data.nofComponents - 1 ) return NULL;
    return data.components[index];
  }
}


enum linDeltaType
GcMulti2D::getDeltaType()
{
  return data.emfComponents[0]->getDeltaType();
}


double
GcMulti2D::getDeltaU()
{
  return data.emfComponents[0]->getDeltaU(); 
}


void
GcMulti2D::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = data.nofEmfComponents;

  types = new linDeltaType[nof_components];
  valuesU = new double[nof_components];
  useFixedN = new bool[nof_components];

  for (int i = 0; i < nof_components; i++) {

    Geometry* sc = data.emfComponents[i];

    types[i] = sc->getDeltaType();
    valuesU[i] = sc->getDeltaU();
    useFixedN[i] = sc->useFixedMeshN();
  }
}


// NOTE: First simple component defines label's position!
//
void
GcMulti2D::getLabelPoint(Point3& point)
{
  data.components[0]->getLabelPoint(point);
}


void
GcMulti2D::getMifTags(int& nof_tags, int*& tags)
{
  nof_tags = getNofComponents();

  tags = new int[nof_tags];

  int index = 0;

  for (int i = 0; i < data.nofComponents; i++) {
    int nof_tagsc;
    int* tagsc;
    data.components[i]->getMifTags(nof_tagsc, tagsc);

    for (int j = 0; j < nof_tagsc; j++) {
      tags[index++] = tagsc[j];
    }

    delete[] tagsc;
    tagsc = NULL;
  }

}



int
GcMulti2D::getNofBoundaryPoints()
{
  int count = 0;

  for (int i = 0; i < data.nofComponents; i++) {
    count += data.components[i]->getNofBoundaryPoints();
  }

  return count;
}


// NOTE: Argument flags tells if only the nof emf components (ie. components defined
// egf- and emf-files) or if the nof all primitive components shouyld be given
//
// NOTE: This number can be different if some of the emf componets are functions which create multiple
// primitive components (circles,polylines)
//
int
GcMulti2D::getNofComponents(bool only_emf_components)
{

  if ( only_emf_components ) {
    return data.nofEmfComponents;
  } else {
    return data.nofComponents;
  }
}


ecif_geometryType
GcMulti2D::getType()
{
  if ( data.nofComponents < 1 ) {
    return ECIF_NODIM;
  } else if ( data.nofComponents > 1 ) {
    return ECIF_MULTI2D; 
  } else {
    return data.components[0]->getType();
  }
}


void
GcMulti2D::init()
{
  data.nofEmfComponents = 0;
  data.emfComponents = NULL;

  data.gmtrType = ECIF_NODIM;
  data.funcType = ECIF_NOFUNCTION;
  data.funcArgC = 0;
  data.funcArgV = NULL;
  data.funcStartPoint = NULL;
  data.funcEndPoint = NULL;
  data.funcStartVertex = NO_INDEX;
  data.funcEndVertex = NO_INDEX;
  data.functionName = NULL;
  data.libraryName = NULL;
}


bool
GcMulti2D::isClosedU()
{
  return data.components[0]->isClosedU();
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns a structure where there is the nof intersection and
// the minimum x-coordinate for the intersections.
RayHit*
GcMulti2D::isectionsWithXdir(GcPoint* xdir_startp, bool& negative_on_left)
{
  return data.components[0]->isectionsWithXdir(xdir_startp, negative_on_left);
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns nof-intersection and the corresponding x-values in the arguments
void
GcMulti2D::isectionsWithXdir(GcPoint* xdir_startp, int& nof_values , double*& values)
{
  data.components[0]->isectionsWithXdir(xdir_startp, nof_values , values);
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and y-axis and which starts from point *startp* towards right.
// Returns nof-intersection and the corresponding y-values in the arguments
void
GcMulti2D::isectionsWithYdir(GcPoint* ydir_startp, int& nof_values , double*& values)
{
  data.components[0]->isectionsWithYdir(ydir_startp, nof_values , values);
}


bool
GcMulti2D::isOnSameAxis(GcPoint& p1, GcPoint& p2)
{
  if ( data.nofComponents < 1 || data.nofComponents > 1 ) {
    return false;
  } else {
    return data.components[0]->isOnSameAxis(p1, p2);
  }
}


ostream&
GcMulti2D::output_emf(ostream& out, short indent_size, short indent_level)
{
  for (int i = 0; i < data.nofEmfComponents; i++) {

    Geometry2D* pg = data.emfComponents[i];

    // Simple geometry
    //
    if ( !pg->isMultiGeometry() ) {
      pg->output_emf(out, indent_size, indent_level);

    // Multi geometry
    //
    } else {
      GcMulti2D* pg = (GcMulti2D*)data.emfComponents[i];

      if ( pg->data.functionName != NULL && pg->data.libraryName != NULL ) {
        pg->output_1_emf(out, indent_size, indent_level);
      } else {
        pg->output_n_emf(out, indent_size, indent_level);
      }
    }
  }

  return out;
}


ostream&
GcMulti2D::output_1_emf(ostream& out, short indent_size, short indent_level)
{
  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;

  const char* def;

  // Geometry type
  def = getMatcString(data.matcTable, EMF_GEOMETRY);
  if ( model->keepMatcDefinitions() && def != NULL ) 
    LibFront::output_matcDef(out, indent_size, indent_level, EMF_GEOMETRY, NULL, def);
  else if ( data.gmtrType == ECIF_CIRCLE ) 
    LibFront::output_scalar(out, is, il, EMF_GEOMETRY, NULL, "Circle");
  else if ( data.gmtrType == ECIF_POLYLINE ) 
    LibFront::output_scalar(out, is, il, EMF_GEOMETRY, NULL, "PolyLine");

  // Argument vector
  if ( data.funcArgC > 0 ) {
    def = getMatcString(data.matcTable, EMF_ARGUMENTS);
    if ( model->keepMatcDefinitions() && def != NULL ) 
      LibFront::output_matcDef(out, indent_size, indent_level, EMF_ARGUMENTS, NULL, def);
    else
      LibFront::output_vector(out, is, il, EMF_ARGUMENTS, NULL, data.funcArgC, data.funcArgV, false);
  }

  // Function definition (Func-name Library-name pair)
  def = getMatcString(data.matcTable, EMF_FUNCTION);
    if ( model->keepMatcDefinitions() && def != NULL ) 
      LibFront::output_matcDef(out, indent_size, indent_level, EMF_FUNCTION, NULL, def);
    else if ( data.functionName != NULL && data.libraryName != NULL ) {
        LibFront::indent(out, is, il);
        out << "Function " 
            << QM << data.functionName << QM
            << " "
            << QM << data.libraryName << QM
            << endl;
    }

  // Start,end points/vertices
  if ( data.funcStartPoint != NULL ) {
    def = getMatcString(data.matcTable, EMF_START_POINT);
    if ( model->keepMatcDefinitions() && def != NULL ) 
      LibFront::output_matcDef(out, indent_size, indent_level, EMF_START_POINT, NULL, def);
    else
      LibFront::output_vector(out, is, il, EMF_START_POINT, NULL, 3, *data.funcStartPoint, false);
  }

  if ( data.funcEndPoint != NULL ) {
    def = getMatcString(data.matcTable, EMF_END_POINT);
    if ( model->keepMatcDefinitions() && def != NULL ) 
      LibFront::output_matcDef(out, indent_size, indent_level, EMF_END_POINT, NULL, def);
    else
      LibFront::output_vector(out, is, il, EMF_END_POINT, NULL, 3, *data.funcEndPoint, false);
  }

  if ( data.funcStartVertex != NO_INDEX ) {
    def = getMatcString(data.matcTable, EMF_START_VERTEX);
    if ( model->keepMatcDefinitions() && def != NULL ) 
      LibFront::output_matcDef(out, indent_size, indent_level, EMF_START_VERTEX, NULL, def);
    else
      LibFront::output_scalar(out, is, il, EMF_START_VERTEX, NULL, data.funcStartVertex);
  }

  if ( data.funcEndVertex != NO_INDEX ) {
    def = getMatcString(data.matcTable, EMF_END_VERTEX);
    if ( model->keepMatcDefinitions() && def != NULL ) 
      LibFront::output_matcDef(out, indent_size, indent_level, EMF_END_VERTEX, NULL, def);
    else
      LibFront::output_scalar(out, is, il, EMF_END_VERTEX, NULL, data.funcEndVertex);
  }

  return out;
}


ostream&
GcMulti2D::output_n_emf(ostream& out, short indent_size, short indent_level)
{
  for (int i = 0; i < data.nofComponents; i++) {
    data.components[i]->output_emf(out, indent_size, indent_level);
  }

  return out;
}


ostream&
GcMulti2D::output_mif(ostream& out, const char* header, bool useGmtrMeshN)
{
  for (int i = 0; i < data.nofComponents; i++) {
      data.components[i]->output_mif(out, header, useGmtrMeshN);
  }

  return out;
}


// Calculate the point corresponding parametric values given as arguments.
// NOTE: Argument parameters are in [0,1], but circle's own parametric
// representation is between [0, 2PI]
void
GcMulti2D::param2Point(double u_p, double v_p, Point3& point)
{
  data.components[0]->param2Point(u_p, v_p, point);
}


// Calculate the point corresponding paramtric values given as arguments.
// NOTE: Argument parameters are in [0,1], but circle's own parametric
// representation is between [0, 2PI]
GcPoint*
GcMulti2D::param2Point(double u_p, double v_p)
{
  return data.components[0]->param2Point(u_p, v_p);
}


// Parameter value for the point (0...1)
void
GcMulti2D::point2Param(Point3* point, double& u, double& v)
{
  data.components[0]->point2Param(point, u, v);
}


// Calculates the paramtetric values corresponding the argument-point.
ParamPair*
GcMulti2D::point2Param(GcPoint* point)
{
  return data.components[0]->point2Param(point);
}


void
GcMulti2D::setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN)
{
  if ( nof_components != data.nofEmfComponents ) return;

  for (int i = 0; i < nof_components; i++) {

    int sub_nof_components = data.emfComponents[i]->getNofComponents();

    if ( sub_nof_components > 1 ) {
      for (int j = 0; j < sub_nof_components; j++) {
        Geometry* sc = data.emfComponents[i]->getComponent(j);
        sc->setDiscretizationData(1, &types[i], &valuesU[i], NULL, &useFixedN[i]);
      }

    } else {
      Geometry* sc = data.emfComponents[i];
      sc->setDiscretizationData(1, &types[i], &valuesU[i], NULL, &useFixedN[i]);
    }
  }
}


// Copy component's data to a transfer-structure for geoemtry update
//
// NOTE: By copying the data into a transfer-structure we can use the
// same create-circle/polyline functions as for original egf/emf data!
//
void
GcMulti2D::setGeometryUpdateData(ecif_ElementComponent_X& txc, ecif_Multi2D data)
{
  static char func[1024];

  const char* def;
  const char* res;

  //txc.lin_delta;
  //txc.lin_delta_type;
  //txc.nof_vertices;
  //txc.vertex_tags;
  //txc.gmtr_type;
  //txc.geometry;

  // Function related stuff
  if ( data.funcType != ECIF_NOFUNCTION ) {
    txc.isFunction = true;

    switch (data.funcType) {
    case ECIF_CPP: txc.isCpp = true; break;
    case ECIF_F95: txc.isF95 = true; break;
    case ECIF_MATC: txc.isMatc = true; break;
    }
  }

  def = getMatcString(data.matcTable, EMF_START_POINT);
  if ( def == NULL ) {
    txc.startPoint = data.funcStartPoint;
  } else {
    res = LibFront::evalMatcString(def);
    if (!LibFront::isMatcError(res)) {
      txc.startPoint = create3(res);
    }
  }

  def = getMatcString(data.matcTable, EMF_END_POINT);
  if ( def == NULL ) {
    txc.endPoint = data.funcEndPoint;
  } else {
    res = LibFront::evalMatcString(def);
    if (!LibFront::isMatcError(res)) {
      txc.endPoint = create3(res);
    }
  }

  def = getMatcString(data.matcTable, EMF_START_VERTEX);
  if ( def == NULL ) {
    txc.startVertex = data.funcStartVertex;
  } else {
    res = LibFront::evalMatcString(def);
    if (!LibFront::isMatcError(res)) {
      txc.startVertex = atoi(res);
    }
  }

  def = getMatcString(data.matcTable, EMF_END_VERTEX);
  if ( def == NULL ) {
    txc.endVertex = data.funcEndVertex;
  } else {
    res = LibFront::evalMatcString(def);
    if (!LibFront::isMatcError(res)) {
      txc.endVertex = atoi(res);
    }
  }

  def = getMatcString(data.matcTable, EMF_FUNCTION);
  if ( def == NULL ) {
    update_dyna_string(txc.functionName, data.functionName);
    update_dyna_string(txc.libraryName, data.libraryName);
  } else {
    res = LibFront::evalMatcString(def);
    if (!LibFront::isMatcError(res)) {
      strstream strm;
      strm << res << ends;
      strm >> func;
      update_dyna_string(txc.functionName, func);
      if (!strm.eof()) {
        update_dyna_string(txc.libraryName, strm.str());
      }
    }
  }

  def = getMatcString(data.matcTable, EMF_ARGUMENTS);
  if ( def == NULL ) {
    txc.argc = data.funcArgC;
    txc.argv = new double[txc.argc];
    for (int i = 0; i < txc.argc; i++) {
      txc.argv = data.funcArgV;
    }
  } else {
    res = LibFront::evalMatcString(def);
    if (!LibFront::isMatcError(res)) {
      strstream strm;
      strm << LibFront::trim((char*)res) << ends;
      int count = 0;
      double value;
      while (true) {
        strm >> value;
        if (0 != strm.fail()) break;
        count++;
      }

      if (count > 0) {
        txc.argc = count;
        txc.argv = new double[count];
        strstream strm;
        strm << res << ends;
        int count = 0;
        while (true) {
          strm >> value;
          if (0 != strm.fail()) break;
          txc.argv[count++] = value;
        }
        txc.argc = count;
      }
    }
  }

  // Finally copy Matc-expressions
  copyMatcValueTable(data.matcTable, txc.matcTable);
}


void
GcMulti2D::setMifTag(int& next_tag)
{
  for (int i = 0; i < data.nofComponents; i++) {
    data.components[i]->setMifTag(next_tag);
  }
}


bool
GcMulti2D::updateData()
{
  calcBoundaryPoints();
  boundbox = calcBoundBox();

  return true;
}


// Update geometry. Relevant when Matc-parameters have been changed!
//
bool
GcMulti2D::updateGeometry(int parent_tag, IdList& vertex_tags)
{
  int i, j;

  ecif_ElementComponent_X txc;
  init_trx_data(txc);

  int final_nof_cmpnts = 0;

  Geometry2D* newGmtr = NULL;

  for (i = 0; i < data.nofEmfComponents; i++) {

    Geometry2D* cmp = data.emfComponents[i];
    
    // Recreate a multi-geometry component (=functions)
    // ----------------------------------
    //
    if ( cmp->isMultiGeometry() ) {

      GcMulti2D* oldGmtr = (GcMulti2D*)cmp;
    
      reset_trx_data(txc, ECIF_EDGE);
      //init_trx_data(txc);
    
      setGeometryUpdateData(txc, oldGmtr->data);

      switch ( oldGmtr->data.gmtrType ) {

      //---Function circle
      // 
      case ECIF_CIRCLE:
      
        if ( oldGmtr->data.funcType == ECIF_CPP )
          newGmtr = createFunctionCircleC(parent_tag, txc, vertex_tags);

        else if ( oldGmtr->data.funcType == ECIF_F95 )
          ;//newGmtr = createFunctionCircleF(parent_tag, txc, vertex_tags);
        
        else if ( oldGmtr->data.funcType == ECIF_MATC )
          newGmtr = createFunctionCircleM(parent_tag, txc, vertex_tags);
        
        break;

      //---Function polyline
      //
      case ECIF_POLYLINE:

        if ( oldGmtr->data.funcType == ECIF_CPP )
          newGmtr = createFunctionPolyLineC(parent_tag, txc, vertex_tags);
        
        else if ( oldGmtr->data.funcType == ECIF_F95 )
          ;//newGmtr = createFunctionPolyLineF(parent_tag, txc, vertex_tags);
        
        else if ( oldGmtr->data.funcType == ECIF_MATC )
          newGmtr = createFunctionPolyLineM(parent_tag, txc, vertex_tags);
        
        break;

      default:
        newGmtr = NULL;
        break;
      }
    
    // Update simple components (and keep it, no copy is made!)
    // ------------------------
    //
    } else {
      newGmtr = cmp;
      newGmtr->updateGeometry(parent_tag, vertex_tags);
    }

    // Error in geometry!
    //
    if ( newGmtr == NULL || !newGmtr->geometryIsOk() ) {
      delete newGmtr;
      geometryOk = false;
      return false;
    }

    data.emfComponents[i] = newGmtr;
    
    int ccount = newGmtr->getNofComponents();
    final_nof_cmpnts += ccount;

  } // Each emf-component
  
  // Delete old simple geometry
  // NOTE: old simple emfComponents should not be deleted here, because
  // they were not copied above!!!
  //
  for (i = 0; i < data.nofComponents; i++) {
    if ( data.isDeletable[i] ) {
      delete data.components[i];
    }
  }
  delete[] data.components;
  delete[] data.isDeletable;
  data.nofComponents = 0;

  // Collect all new (simple) components into one array
  //
  Geometry2D** tmp = new Geometry2D*[final_nof_cmpnts];

  data.isDeletable = new bool[final_nof_cmpnts];

  int idx = 0;
  for (i = 0; i < data.nofEmfComponents; i++) {

    Geometry2D* pg = data.emfComponents[i];

    // If a component is simple, take as it is
    //
    if ( !pg->isMultiGeometry() ) {
      tmp[idx] = pg;
      data.isDeletable[idx] = false;
      idx++;

    // If a component is a Multi2D, copy its subcomponents!
    // NOTE: Sub components should have been already unpurged
    // recursively!
    //
    } else {
      for (j = 0; j < pg->getNofComponents(); j++) {
        tmp[idx] = (Geometry2D*)pg->getComponent(j);
        data.isDeletable[idx] = true;
        idx++;
      }
    }
  }

  // These are components which all are geometry primitves (polylines,circles etc)
  //
  data.nofComponents = final_nof_cmpnts;
  data.components = tmp;

  return updateData();
}


bool
GcMulti2D::useFixedMeshN()
{
  // Hack hack hack!!!
  try {
    return data.emfComponents[0]->useFixedMeshN();
  }

  catch (...) {
    return false;
  }
}



//*******************************
//***** GcNurbsCurve methods *****
//*******************************

// Constructor
// Emf-style call (In practice, nurbs are not acutally input from emf!)
GcNurbsCurve::GcNurbsCurve(ecif_ElementComponent_X& tx, IdList& vertex_tags)
{
  nofVertices = tx.nof_vertices;
  vertexTags = NULL;

  if ( nofVertices > 0 ) {

    vertexTags = new int[nofVertices];

    for (int i = 0; i < nofVertices; i++) {

      BodyElement* v = model->getVertexByTag(tx.vertex_tags[i]);
      vertex_tags.push_back(v->Tag());
      vertexTags[i] = v->Tag();
    }
  }

  data.deltaType = tx.lin_delta_type;
  data.delta = tx.lin_delta[0];
  data.useFixedMeshN = (bool)tx.use_fixed_mesh_n;

  create(tx.geometry.edge);

  copyMatcValueTable(tx.matcTable, data.matcTable);
}


// Constructor
// NOTE: No vertices
GcNurbsCurve::GcNurbsCurve(ecif_EdgeGeometry_X* params, bool do_linearize)
{
  nofVertices = 0;
  vertexTags = NULL;
  create(params, do_linearize);
}


// Constructor
// Iges/Ideas-style call
GcNurbsCurve::GcNurbsCurve(BodyElement* vertex1, BodyElement* vertex2,
                           ecif_EdgeGeometry_X* params, bool do_linearize)
{
  nofVertices = 2;
  vertexTags = new int[2];

  vertexTags[0] = vertex1->Tag();
  vertexTags[1] = vertex2->Tag();

  create(params, do_linearize);
}


// Create nurbs curve
void
GcNurbsCurve::create(ecif_EdgeGeometry_X* params, bool do_linearize)
{
  int i,j;
  data.degree = params->degree;
  data.nofKnots = params->nofKnots;
  data.nofCpoints = params->nofCpoints;

  data.knots = new double[data.nofKnots];

  for (i=0; i < data.nofKnots; i++)
    data.knots[i] = params->knots[i];

  data.cpoints = new Point4[data.nofCpoints];

  for (i=0; i < data.nofCpoints; i++) {
    for(j=0; j < 4; j++)
      data.cpoints[i][j]  = params->cpoints[i][j];
  }

  // Linearization factor
  deltaU = 1.0 / 36;

  nofLinearizingPoints = 0;
  linearizingPoints = NULL;

  nofBoundaryPoints = 0;
  boundaryPoints = NULL;

  if (do_linearize) {
    linearize();
  }

  boundbox = calcBoundBox();
}


// Destructor.
GcNurbsCurve::~GcNurbsCurve()
{
  delete[] vertexTags;
}


// Calculates the bounding box for a nurbs curve.
BoundBox*
GcNurbsCurve::calcBoundBox()
{
  int i;

  // We calculate the simple min-max box from control points.
  // (based on the convex hull proprety of the control polygon for nurbs)
  RangeVector values;

  for (i = 0; i < MAX_DIMENSION; i++) {
    values[2*i] = MAX_RANGE;
    values[2*i + 1] = MIN_RANGE;
  }

  const int DIM2 = 2; // We are in 2D !!!***!!!

  for (i = 0; i < data.nofCpoints; i++) {

    for (int k = 0; k < DIM2; k++) {
      if (data.cpoints[i][k] < values[2*k])
        values[2*k] = data.cpoints[i][k];
      if (data.cpoints[i][k] > values[2*k + 1])
        values[2*k + 1] = data.cpoints[i][k];
    }
  }

  BoundBox* bbox = new BoundBox(values);

  return bbox;
}


// Calculate linearized boundary
// Store result in the points-buffer
// NOTE: Buffer is allocated here, client should take care of its
// deletion!!!
//
// NOTE: Pure carbage!!!
//
void
GcNurbsCurve::calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u)
{
  nof_points = 2;

  points = new GcPoint*[nof_points];

  points[0] = data.start;
  points[1] = data.end;
}


void
GcNurbsCurve::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
   renderer->drawNurbsCrv(dmode, dstate, direction, data, elem_id);
}


void
GcNurbsCurve::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 1;
  types = new enum linDeltaType[1];
  types[0] = data.deltaType;

  valuesU = new double[1];
  valuesU[0] = data.delta;

  valuesV = NULL;

  useFixedN = new bool[1];
  useFixedN[0] = data.useFixedMeshN;
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// (For a line-segment (*this*) nof intersections is of course 0 or 1)
RayHit*
GcNurbsCurve::isectionsWithXdir(GcPoint* xdir_startp, bool& negative_on_left)
{
  // Currently nonsense!  !!!***!!!
  // We create a line from end-points and test with that!
  GcLine line(data.start->getPoint(), data.end->getPoint());
  return line.isectionsWithXdir(xdir_startp, negative_on_left);
}


void
GcNurbsCurve::linearize()
{
  nofLinearizingPoints = 0;
  delete[] linearizingPoints;
  linearizingPoints = NULL;

  calcLinearizingPoints(nofLinearizingPoints, linearizingPoints, deltaU);

  model->updateMinimumEdgeSize(nofLinearizingPoints, linearizingPoints);
}


ostream&
GcNurbsCurve::output_emf(ostream& out, short indent_size, short indent_level)
{
  return out;
}


// Calculate the point corresponding paramtric values given as arguments.
void
GcNurbsCurve::param2Point(double u_p, double v_p, Point3& point)
{
  nurbsCurveValue(&data, u_p, point);
}


// Calculate the point corresponding paramtric values given as arguments.
GcPoint*
GcNurbsCurve::param2Point(double u_p, double v_p)
{
  Point3 point;
  nurbsCurveValue(&data, u_p, point);

  GcPoint* p = new GcPoint(point);

  return p;
}


// Calculates the paramtetric values corresponding the argument-point.
void
GcNurbsCurve::point2Param(Point3* p, double& u, double& v)
{
  // Currently nearly nonsense!  !!!***!!! :
  // Curve is suppoused to go through first and last control-points!
  // Only these "end-points" are checked!
  // Other cases are suppoused to be OFF the curve!
  Point4& c1 = data.cpoints[0];
  Point4& c2 = data.cpoints[data.nofCpoints - 1];

  Point3 p1;
  Point3 p2;

  for (int i = 0; i < 3; i++) {
    p1[i] = c1[i];
    p2[i] = c2[i];
  }

  //If point doesn't p1 match one of the "end-points"
  //we return NULL
  if ( samepoint(*p, p1) ) {
    u = 0;
  } else if ( samepoint(*p, p2) ) {
    u = 1;
  } else {
    u = -1.0;
  }

  v = 0.0;
}


// Calculates the paramtetric values corresponding the argument-point.
ParamPair*
GcNurbsCurve::point2Param(GcPoint* p)
{
  // Currently nearly nonsense!  !!!***!!! :
  // Curve is suppoused to go through first and last control-points!
  // Only these "end-points" are checked!
  // Other cases are suppoused to be OFF the curve!
  double u = 0;
  Point4& c1 = data.cpoints[0];
  Point4& c2 = data.cpoints[data.nofCpoints - 1];
  GcPoint p1(c1[0], c1[1], c1[2]);
  GcPoint p2(c2[0], c2[1], c2[2]);

  //If point doesn't p1 match one of the "end-points"
  //we return NULL
  if (*p == p1)
    u = 0;
  else if (*p == p2)
    u = 1;
  else
    return NULL;

  // Now we can create non-null result-object.
  ParamPair* par = new ParamPair;
  par->u = u;
  par->v = 0;
  return par;
}


// Checks if point is on the curve.
// Currently nonsense!  !!!***!!!
ParamPair*
GcNurbsCurve::point2Param(GcPoint* p, int direction,
                 double start_u, double start_v)
{
  // Currently nonsense!  !!!***!!!
  ParamPair* par = new ParamPair;
  par->u = start_u;
  par->v = start_v;
  return par;
}


void
GcNurbsCurve::setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN)
{
  if ( nof_components = 0 ) return;

  data.deltaType = types[0];
  data.delta = valuesU[0];
  data.useFixedMeshN = useFixedN[0];
  
  // Update also deltaU which is always in radians
  //
  //setDeltaU(data.deltaType, data.delta);
}


//*******************************
//*****GcPolyLine methods *****
//*******************************

GcPolyLine::GcPolyLine(ecif_ElementComponent_X& tx, IdList& vertex_tags)
{
  nofVertices = tx.nof_vertices;
  vertexTags = new int[nofVertices];
  
  int i;

  for (i = 0; i < nofVertices; i++) {
    vertexTags[i] = tx.vertex_tags[i];
    vertex_tags.push_back(tx.vertex_tags[i]);
  }

  GcPoint** points = NULL;
  int nof_points = tx.geometry.edge->nofDefiningPoints;

  // Defining points given
  // ---------------------
  if ( nof_points > 0 ) {
    points = new GcPoint*[nof_points];

    for (i = 0; i < nof_points; i++) {
      points[i] = new GcPoint(&tx.geometry.edge->definingPoints[i]);
    }

  // Defining vertices given
  // -----------------------
  } else {
    nof_points = nofVertices;
    points = new GcPoint*[nof_points];

    for (i = 0; i < nofVertices; i++) {
      BodyElement* v = model->getVertexByTag(tx.vertex_tags[i]);
      vertexTags[i] = v->Tag();
      vertex_tags.push_back(v->Tag());
      points[i] = (GcPoint*)v->getGeometry();
    }
    
    data.hasVertexTies = true;
  }

  data.useFixedMeshN = (bool)tx.use_fixed_mesh_n;

  create(nof_points, points);

  delete[] points;

  copyMatcValueTable(tx.matcTable, data.matcTable);
}


// Constructor.
GcPolyLine::GcPolyLine(int nof_points, GcPoint** gc_points)
{
  nofVertices = 0;
  vertexTags = NULL;

  create(nof_points, gc_points);
}


GcPolyLine::GcPolyLine(int nof_vertices, BodyElement** vertices)
{
  nofVertices = nof_vertices;
  vertexTags = new int[nofVertices];

  GcPoint** points = new GcPoint*[nofVertices];

  for (int i = 0; i < nofVertices; i++) {
    BodyElement* v = vertices[i];
    vertexTags[i] = v->Tag();

    points[i] = (GcPoint*)v->getGeometry();
  }

  data.hasVertexTies = true;

  create(nofVertices, points);

  delete[] points;
}


// Destructor.
GcPolyLine::~GcPolyLine()
{
  delete[] vertexTags;
}


// Create polyline
void
GcPolyLine::create(int nof_points, GcPoint** points)
{
  data.nofPoints = nof_points;
  data.points = new GcPoint*[data.nofPoints];

  for (int i = 0; i < data.nofPoints; i++)
    data.points[i] = points[i];

  nofBoundaryPoints = 0;
  boundaryPoints = NULL;

  updateData();
}



// Calculates bounding box for the polyline
BoundBox*
GcPolyLine::calcBoundBox()
{
  RangeVector range;
  BoundBox* bbox = new BoundBox();

  for (int n = 1; n < data.nofPoints; n++) {
    // A pair of consequtive points as a line
    GcPoint* start = data.points[n - 1];
    GcPoint* end = data.points[n];

    range[0] = start->Pos(X);
    range[1] = end->Pos(X);
    range[2] = start->Pos(Y);
    range[3] = end->Pos(Y);
    range[4] = start->Pos(Z);
    range[5] = end->Pos(Z);

    bbox->extendByRange(range);
  }

  return bbox;
}


void
GcPolyLine::calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u)
{
  nof_points = data.nofPoints;
  points = new GcPoint*[data.nofPoints];

  for (int i = 0; i < data.nofPoints; i++) {
    points[i] = data.points[i];
  }
}


void
GcPolyLine::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
  const Point3* p1;
  const Point3* p2;

  // Positive direction: from start to end
  if ( direction == 1 ) {
    for (int n = 0; n < data.nofPoints -1; n++) {
      // A pair of consequtive points as a line
      p1 = data.points[n]->getPoint();
      p2 = data.points[n + 1]->getPoint();
      renderer->drawLine(dmode, dstate, 1, p1, p2, elem_id);
    }

  // Negative direction: from end to start
  } else {
    for (int n = data.nofPoints - 1; n > 0; n--) {
      // A pair of consequtive points as a line
      p1 = data.points[n]->getPoint();
      p2 = data.points[n - 1]->getPoint();
      renderer->drawLine(dmode, dstate, 1, p1, p2, elem_id);
    }
  }
}


void
GcPolyLine::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 1;

  types = NULL;
  valuesU = NULL;
  valuesV = NULL;

  useFixedN = new bool[1];
  useFixedN[0] = data.useFixedMeshN;
}


void
GcPolyLine::getLabelPoint(Point3& point)
{
  // Label in the middle of the center segment

  Point3 p1;
  Point3 p2;

  data.points[data.nofPoints/2 - 1]->getPoint(p1);
  data.points[data.nofPoints/2]->getPoint(p2);

  point[0] = 0.5 * (p1[0] + p2[0]);
  point[1] = 0.5 * (p1[1] + p2[1]);
  point[2] = 0.0;
}


bool
GcPolyLine::getLine(int index, GcLine*& line)
{   
  if ( index < 0 || index > getNofLines() - 1) {
    return false;
  }

  Point3* p1 = data.points[index]->getPoint();
  Point3* p2 = data.points[index + 1]->getPoint();
   
  line = new GcLine(p1, p2);

  return true;
}


bool
GcPolyLine::isClosedU()
{
  if ( data.nofPoints > 3 &&
       *data.points[0] == *data.points[data.nofPoints - 1]
     ) {
    return true;

  } else {
    return false;
  }
}


// Counts the number of intersections with a line which is parallel to
// xy-plane and x-axis and which starts from point *startp* towards right.
// Returns a structure where there is the nof intersection and
// the minimum x-coordinate for the intersections.
RayHit*
GcPolyLine::isectionsWithXdir(GcPoint* xdir_startp, bool& negative_on_left)
{
  RayHit* result = new RayHit;

  int index = 0;
  GcLine* line;

  while (true) {
    // No more lines!
    if ( !getLine(index++, line) ) {
      break;
    }
    
    bool neg_on_left;
    RayHit* hits = line->isectionsWithXdir(xdir_startp, neg_on_left);

    if ( hits != NULL ) {

      result->count++;
      result->points.push_back(hits->points.front());

      if ( hits->min_value < result->min_value ) {
        result->min_value = hits->min_value;

        // NOTE: Negative on left indicator is picked
        // from the "leftmost" hit, this way we can conclude the
        // possible loop orientation correctly
        negative_on_left = neg_on_left;
      }

    }
    
    delete line;
    delete hits;
  }

  if (result->count == 0) {
    delete result;
    result = NULL;
  }

  return result;
}


ostream&
GcPolyLine::output_emf(ostream& out, short indent_size, short indent_level)
{
  LibFront::output_vector(out, indent_size, indent_level, EMF_VERTICES, NULL, nofVertices, vertexTags, false);

  if ( data.useFixedMeshN ) {
    LibFront::output_scalar(out, indent_size, indent_level, EMF_USE_MESH_N, NULL, 1);
  }

  return out;
}


// Calculate the point corresponding the paramtric values given as arguments.
// Note param value is n + u_p, where n is line index and u_p is the paramter value
// within the line
//
void
GcPolyLine::param2Point(double u_p, double v_p, Point3& point)
{
  point[0] = point[1] = point[2] = 0.0;

  int index = int(u_p);

  GcLine* line = NULL;
  
  getLine(index, line);

  if ( line == NULL ) {
    return;
  }

  line->param2Point(u_p - index, v_p, point);

  delete line;
}


ParamPair*
GcPolyLine::point2Param(GcPoint* p)
{
  GcLine* line = NULL;
  int index = 0;

  while (1) {

    getLine(index++, line);

    if ( line == NULL ) {
      break;
    }

    ParamPair* pair = line->point2Param(p);

    if ( pair != NULL ) {
      return pair;
    }

    delete line;
  }

  return NULL;
}


void
GcPolyLine::setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN)
{
  if ( nof_components = 0 ) return;

  data.useFixedMeshN = useFixedN[0];
}

bool
GcPolyLine::updateData()
{
  boundbox = calcBoundBox();

  return true;
}


// Update geometry. Relevant when Matc-parameters have been changed!
//
bool
GcPolyLine::updateGeometry(int parent_tag, IdList& vertex_tags)
{
  return updateData();
}



//**************************************************************************
//                        Geometry3D methods
//**************************************************************************
//
// NOTE: Badly incomplete (3D Cad geometry not supported!)

Geometry3D::Geometry3D()
{
}

Geometry3D::~Geometry3D()
{
}


// Linearize geometry and store points in boundaryPoints array

void
Geometry3D::calcBoundaryPoints(int boundary_tag)
{
  int i;

  for (i = 0; i < nofBoundaryPointsU * nofBoundaryPointsV ; i++) {
    delete boundaryPoints[i];
  }
  delete[] boundaryPoints;

  nofBoundaryPointsU = 0;
  nofBoundaryPointsV = 0;
  boundaryPoints = NULL;

  int nof_points_u;
  int nof_points_v;
  GcPoint** points;

  double delta_u = getParamDeltaU();
  double delta_v = getParamDeltaV();

  calcLinearizingPoints(nof_points_u, nof_points_v, points, delta_u, delta_v);

  int nof_points = 0;

  nof_points = nof_points_u * nof_points_v;

  if (nof_points == 0) return;

  boundaryPoints = new BoundaryPoint*[nof_points];

  for (i = 0; i < nof_points; i++) {

    BoundaryPoint* bp = new BoundaryPoint;

    boundaryPoints[i] = bp;

    // Check if bp is an existing vertex!
    BodyElement* v = model->getVertex(points[i]);

    if ( v != NULL ) {
      bp->tag = v->Tag();
      bp->vertexTag = v->Tag();

    }

    //bp->boundaryTag = boundary_tag;
    bp->point = points[i];

  }

  // Handle possibly closed geometry
  // -------------------------------
  int nof_u = nof_points_u;
  int nof_v = nof_points_v;

  //-Copy start points to end points if U-closed
  //
  if ( isClosedU() ) {

    for (int i = 1; i <= nof_v; i++) {
      int pos0 = (i - 1) * nof_u;
      int pos1 = pos0 + nof_u - 1;

      BoundaryPoint* tmp = boundaryPoints[pos1];
      boundaryPoints[pos1] = boundaryPoints[pos0];

      delete tmp;
    }
  }

  //-Copy start points to end points if V-closed
  //
  if ( isClosedV() ) {

    for (int i = 1; i <= nof_u; i++) {
      int pos0 = (i - 1) * nof_v;
      int pos1 = pos0 + nof_v - 1;

      BoundaryPoint* tmp = boundaryPoints[pos1];
      boundaryPoints[pos1] = boundaryPoints[pos0];

      delete tmp;
    }
  }

  nofBoundaryPointsU = nof_points_u;
  nofBoundaryPointsV = nof_points_v;

  // NOTE: Delete only the array,not the points!
  //
  delete[] points;
}



//**********************************
//***** GcNurbsSurface methods *****
//**********************************
//
// NOTE: Badly incomplete (3D Cad geometry not supported!)

GcNurbsSurface::GcNurbsSurface(ecif_ElementComponent_X& tx, IdList& vertex_tags)
{
  nofVertices = 0;
  vertexTags = NULL;

  copyMatcValueTable(tx.matcTable, data.matcTable);
}


// Constructor.
GcNurbsSurface::GcNurbsSurface(ecif_FaceGeometry_X* params, bool do_linearize)
{
  nofVertices = 0;
  vertexTags = NULL;
  create(params, do_linearize);
}


// Create surface
void
GcNurbsSurface::create(ecif_FaceGeometry_X* params, bool do_linearize)
{
  int i,j;

  data.isRational = params->isRational;

  data.degree_u = params->degree_u;
  data.nofKnots_u = params->nofKnots_u;
  data.knots_u = new double[data.nofKnots_u];
  for (i=0; i < data.nofKnots_u; i++)
    data.knots_u[i] = params->knots_u[i];

  data.degree_v = params->degree_v;
  data.nofKnots_v = params->nofKnots_v;
  data.knots_v = new double[data.nofKnots_v];
  for (i=0; i < data.nofKnots_v; i++)
    data.knots_v[i] = params->knots_v[i];

  data.nofCpoints_u = params->nofCpoints_u;
  data.nofCpoints_v = params->nofCpoints_v;
  data.nofCpoints = params->nofCpoints;

  data.cpoints = new Point4[data.nofCpoints];

  for (i=0; i < data.nofCpoints; i++) {
    for(j=0; j < 4; j++)
      data.cpoints[i][j] = params->cpoints[i][j];
  }

  // Linearization factors
  deltaU = 1.0 / 36;
  deltaV = 1.0 / 36;

  nofLinearizingPointsU = 0;
  nofLinearizingPointsV = 0;
  linearizingPoints = NULL;

  nofBoundaryPointsU = 0;
  nofBoundaryPointsV = 0;
  boundaryPoints = NULL;

  if (do_linearize) {
    linearize();
  }

  boundbox = calcBoundBox();
}


// Destructor.
GcNurbsSurface::~GcNurbsSurface()
{
  delete[] vertexTags;
}


BoundBox*
GcNurbsSurface::calcBoundBox()
{
  int i;

  // We calculate the simple min-max box from control points.
  // (based on the convex hull proprety of the control polygon for nurbs)
  RangeVector range;

  for (i = 0; i < MAX_DIMENSION; i++) {
    range[2*i] = MAX_RANGE;
    range[2*i + 1] = MIN_RANGE;
  }

  const int DIM3 = 3; // We are in 3D !!!***!!!

  for (i = 0; i < data.nofCpoints; i++) {

    for (int k = 0; k < DIM3; k++) {
      if (data.cpoints[i][k] < range[2*k])
        range[2*k] = data.cpoints[i][k];
      if (data.cpoints[i][k] > range[2*k + 1])
        range[2*k + 1] = data.cpoints[i][k];
    }
  }

  BoundBox* bbox = new BoundBox(range);

  return bbox;
}


// Calculate linearized boundary
// Store result in the points-buffer
// NOTE: Buffer is allocated here, client should take care of its
// deletion!!!
// NOTE: Pure carbage!!!
void
GcNurbsSurface::calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u, double delta_v)
{
  nof_points_u = 2;
  nof_points_v = 2;

  points = new GcPoint*[nof_points_u * nof_points_v];

  for (int i = 0; i < nof_points_u * nof_points_v ; i++) {
    points[i] = data.corners[i];
  }
}


void
GcNurbsSurface::draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id)
{
   renderer->drawNurbsSrf(dmode, dstate, direction, data, elem_id);
}


void
GcNurbsSurface::linearize()
{
  nofLinearizingPointsU = 0;
  nofLinearizingPointsV = 0;
  delete[] linearizingPoints;
  linearizingPoints = NULL;

  calcLinearizingPoints(nofLinearizingPointsU, nofLinearizingPointsV, linearizingPoints, deltaU, deltaV);

  model->updateMinimumEdgeSize(nofLinearizingPointsU * nofLinearizingPointsV, linearizingPoints);

}



ostream&
GcNurbsSurface::output_emf(ostream& out, short indent_size, short indent_level)
{
  return out;
}


//**************************
//*****GcPlane methods *****
//**************************
//
// NOTE: Badly incomplete (3D Cad geometry not supported!)

GcPlane::GcPlane(ecif_ElementComponent_X& tx, IdList& vertex_tags)
{
  nofVertices = 0;
  vertexTags = NULL;
  copyMatcValueTable(tx.matcTable, data.matcTable);
}

// Constructor.
GcPlane::GcPlane(int nof_points, GcPoint** points)
{
  nofVertices = 0;
  vertexTags = NULL;
  create(nof_points, points);
}


// Create plane
void
GcPlane::create(int nof_points, GcPoint** points)
{
  if (nof_points != 4 )
    return;

  for (int i = 0; i < 4; i++ ) {
    data.corners[i] = points[i];
  }

  boundbox = calcBoundBox();
}


// Destructor.
GcPlane::~GcPlane()
{
  delete[] vertexTags;
}


BoundBox*
GcPlane::calcBoundBox()
{
   return 0;
}


ostream&
GcPlane::output_emf(ostream& out, short indent_size, short indent_level)
{

  return out;
}




// ==========================================================================
// ==========================================================================
//                      Geometry data structs
// ==========================================================================
// ==========================================================================


// ==========================================================================
//                            2D classes
// ==========================================================================


// Struct for circle data
//
ecif_Circle::ecif_Circle() { init();}
ecif_Circle::~ecif_Circle() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  purgeMatcValueTable(matcTable);
}

void
ecif_Circle::init() {
  hasVertexTies = false;
  center[0] = center[1] = center[2] = 0.0;
  radius = 0.0;
  isClosed = true;
  start = NULL;
  end = NULL;
  nofDefiningPoints = 0;
  definingPoints = NULL;
  start_u = 0.0;
  end_u = 0.0;
  centerGiven = false;
  radiusGiven = false;
  delta = 0.0;
  deltaType = LIN_DELTA_NONE;
  useFixedMeshN = false;

}


// Struct for curve data
//
ecif_Curve::ecif_Curve() { init();}

ecif_Curve::~ecif_Curve() {
  if (!hasVertexTies) {
    for (int i = 0; i < nofPoints; i++) {
      delete points[i];
    }
  }
  delete[] points;
  delete[] vertex_flags;
  purgeMatcValueTable(matcTable);
}
void
ecif_Curve::init() {
  hasVertexTies = false;
  isClosed = true;
  nofPoints = 0;
  points = NULL;
  vertex_flags = NULL;
}


// Struct for ellipse data
//
ecif_Ellipse::ecif_Ellipse() { init();}
ecif_Ellipse::~ecif_Ellipse() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  purgeMatcValueTable(matcTable);
}
void
ecif_Ellipse::init() {
  hasVertexTies = false;
  center[0] = center[1] = center[2] = 0.0;
  radius1 = 0.0;
  radius2 = 0.0;
  isClosed = true;
  start = NULL;
  end = NULL;
}


// Struct for hyperbola data
//
ecif_Hyperbola::ecif_Hyperbola() { init();}
ecif_Hyperbola::~ecif_Hyperbola() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  purgeMatcValueTable(matcTable);
}
void
ecif_Hyperbola::init() {
  hasVertexTies = false;
  center[0] = center[1] = center[2] = 0.0;
  direction[0] = direction[1] = direction[2] = 0.0;
  radius1 = 0.0;
  radius2 = 0.0;
  start = NULL;
  end = NULL;
}


// Struct for line data
//
ecif_Line::ecif_Line() { init();}
ecif_Line::~ecif_Line() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  purgeMatcValueTable(matcTable);
}
void
ecif_Line::init() {
  hasVertexTies = false;
  onSymmAxis = 0;
  normal[0] = normal[1] = normal[2] = 0.0;
  length = 0.0;
  start = NULL;
  end = NULL;
  dir = NULL;
  useFixedMeshN = false;
}


// Struct for multi-geometry data
//
ecif_Multi2D::ecif_Multi2D() { init();}
ecif_Multi2D::~ecif_Multi2D() {
  delete[] emfComponents;
  delete[] components;
  delete[] functionName;
  delete[] libraryName;
  purgeMatcValueTable(matcTable);
}
void
ecif_Multi2D::init() {
  hasVertexTies = false;

  nofEmfComponents = 0;
  emfComponents = NULL;

  nofComponents = 0;
  isDeletable = NULL;
  components = NULL;

  gmtrType = ECIF_NODIM;

  funcType = ECIF_NOFUNCTION;
  funcArgC = 0;
  funcArgV = NULL;
  funcStartPoint = NULL;
  funcEndPoint = NULL;
  funcStartVertex = NO_INDEX;
  funcEndVertex = NO_INDEX;
  functionName = NULL;
  libraryName = NULL;

  // NOTE: These are valid when there is only one emf component!
  deltaType = LIN_DELTA_NONE;
  delta = 0.0;
  useFixedMeshN = false;
}


// Struct for nurbs curve data

ecif_NurbsCurve::ecif_NurbsCurve() { init();}

ecif_NurbsCurve::~ecif_NurbsCurve() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  delete[] knots;
  delete[] cpoints;
  purgeMatcValueTable(matcTable);
}
void
ecif_NurbsCurve::init() {
  hasVertexTies = false;
  isRational = 0;
  degree = 0;
  nofKnots = 0;
  nofCpoints = 0;
  knots = NULL;
  cpoints = NULL;
  start = NULL;
  end = NULL;
  useFixedMeshN = false;
}


// Struct for parabola data
//
ecif_Parabola::ecif_Parabola() { init();}
ecif_Parabola::~ecif_Parabola() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  purgeMatcValueTable(matcTable);
}
void
ecif_Parabola::init() {
  hasVertexTies = false;
  apex[0] = apex[1] = apex[2] = 0.0;
  direction[0] = direction[1] = direction[2] = 0.0;
  focalLength = 0.0;
  start = NULL;
  end = NULL;
}
  

// Struct for polyline data
//
ecif_PolyLine::ecif_PolyLine() { init();}
ecif_PolyLine::~ecif_PolyLine() {
  if (!hasVertexTies) {
    for (int i = 0; i < nofPoints; i++) {
      delete points[i];
    }
  }
  delete[] points;
  purgeMatcValueTable(matcTable);
}
void
ecif_PolyLine::init() {
  hasVertexTies = false;
  nofPoints = 0;
  points = NULL;
  useFixedMeshN = false;
}

// Struct for spline curve data
//
ecif_SplineCurve::ecif_SplineCurve() { init();}
ecif_SplineCurve::~ecif_SplineCurve() {
  if (!hasVertexTies) {
    delete start;
    delete end;
  }
  delete[] cpoints;
  purgeMatcValueTable(matcTable);
}
void
ecif_SplineCurve::init() {
  hasVertexTies = false;
  degree = 0;
  nofCpoints = 0;
  cpoints = NULL;
  start = NULL;
  end = NULL;
}
