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
Module:     ecif_geometry.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Classe for geometrical forms (surface, edge, point).

************************************************************************/
#ifndef _ECIF_GEOMETRY_
#define _ECIF_GEOMETRY_

#include "front_egfdefs.h"

#include "ecif_def.h"
#include "ecif_def_stl.h"
#include "ecif_def_trx.h"
#include "ecif_boundbox.h"

class Geometry2D;
class BodyElement;
class Model;
class Renderer;


// Used when returning results from ray-tracing tests
struct RayHit {
  RayHit();
  ~RayHit();
  int count;  // nof hits with the ray
  Point3List points;
  double min_value; // min. coordinate value (X,Y or Z) of the hits
};


class Geometry
{
friend class BodyElement;
public:
  Geometry();
  virtual ~Geometry();
  virtual BoundBox* calcBoundBox() {return NULL; }
  virtual void calcBoundaryPoints() {}
  virtual void calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u = 0.0, double delta_v = 0.0);
  virtual void draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate) {}
  virtual void draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction);
  virtual void draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id);
  virtual void draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id, bool is_first_loop);
  virtual void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate) {}
  virtual void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction);
  virtual void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id);
  virtual void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id, bool is_first_loop);
  virtual bool geometryIsOk() { return geometryOk; }
  BoundBox* getBoundBox() {return boundbox;}
  virtual void getBoundaryPoints(int& count, BoundaryPoint**& points);
  virtual Geometry* getComponent(int index, bool only_emf_components = false) { return NULL; }
  virtual double getDeltaU() { return -1.0; }
  virtual double getDeltaV() { return -1.0; }
  virtual linDeltaType getDeltaType() { return LIN_DELTA_NONE; }
  virtual void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  virtual void getLabelPoint(Point3& point);
  virtual double getMeshH() { return -1.0; }
  virtual void getMifTags(int& nof_tags, int*& tags);
  virtual int getNofBoundaryPoints() { return 0; }
  virtual int getNofComponents(bool only_emf_components = false) { return 1; }
  virtual int getNofVertices() { return nofVertices; }
  virtual double getParamDeltaU() { return -1.0; }
  virtual double getParamDeltaV() { return -1.0; }
  void getRangeVector(RangeVector rv) {boundbox->getRangeVector(rv);}
  virtual ecif_geometryType getType() = 0;  // ==> Abstract class!
  static void initClass(Model* model) {Geometry::model = model;}
  virtual bool isClosedU() { return false; }
  virtual bool isClosedV() { return false; }
  virtual bool isDiscretized() { return false; }
  virtual RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left) {negative_on_left = false; return NULL;}
  virtual void isectionsWithXdir(GcPoint* startp, int& nof_values, double*& values) {nof_values = 0;}
  virtual void isectionsWithYdir(GcPoint* startp, int& nof_values, double*& values) {nof_values = 0;}
  virtual bool isMultiGeometry() {return false;}
  virtual bool isOnSameAxis(GcPoint& p1, GcPoint& p2) {return false;}
  virtual ostream& output_emf(ostream& out, short indent_size, short indent_level) { return out;}
  virtual ostream& output_mif(ostream& out, const char* header, bool useGmtrMeshN) { return out;}
  virtual void param2Point(double u_p, double v_p, Point3& point) {}
  virtual GcPoint* param2Point(double u_p, double v_p) {return NULL;}
  virtual void point2Param(Point3* point, double& u_p, double& v_p) {}
  virtual ParamPair* point2Param(GcPoint* p) {return NULL;}
  virtual void setID(int point_id) {}
  virtual void setDeltaU(double delta_u) {}
  virtual void setDeltaU(enum linDeltaType lin_delta_type, double lin_delta) {}
  virtual void setDeltaV(double delta_v) {}
  virtual void setDeltaV(enum linDeltaType lin_delta_type, double lin_delta) {}
  virtual void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN) {}
  virtual void setMifTag(int& next_tag) {}
  virtual bool updateGeometry(int parent_tag, IdList& vertex_tags) {return true;}
  virtual bool useFixedMeshN() { return false; }

protected:
  virtual ostream& outputBoundaryPointTags(ostream& out, int indent_size) {return out;}

  BoundBox* boundbox;
  bool geometryOk;
  static Model* model;

  int nofVertices;
  int* vertexTags;
};


// ==========================================================================
//                            1D classes
// ==========================================================================

class Geometry1D : public Geometry
{
friend class BodyElement;
public:
  Geometry1D();
  virtual ~Geometry1D();

protected:

};


class GcPoint: public Geometry1D
{
  friend class BodyElement;
  friend GcPoint operator+(const GcPoint &p1, const GcPoint &p2);
  friend GcPoint operator-(const GcPoint &p1, const GcPoint &p2);
  friend GcPoint operator*(int num, const GcPoint &p);
  friend GcPoint operator/(const GcPoint &p, int num);
  friend bool operator==(const GcPoint &p1, const GcPoint &p2);
  friend bool operator!=(const GcPoint &p1, const GcPoint &p2);
  friend bool operator<(const GcPoint &p1, const GcPoint &p2);
  friend bool operator>(const GcPoint &p1, const GcPoint &p2);
  friend bool operator<=(const GcPoint &p1, const GcPoint &p2);
  friend bool operator>=(const GcPoint &p1, const GcPoint &p2);
  friend ostream &operator<<(ostream &stream, const GcPoint &p);

public:
  GcPoint();
  GcPoint(double x, double y, double z);
  GcPoint(double p[3]);
  GcPoint(Point3* p);
  GcPoint(const GcPoint& p);
  GcPoint(const GcPoint* p);
  ~GcPoint(){};
  GcPoint average(GcPoint &p2, GcPoint&p3);
  void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate);
  void find_space(GcPoint &pmin, GcPoint &pmax);
  void getLabelPoint(Point3& point);
  void getPoint(Point3& p) {p[0] = posit[0], p[1] = posit[1]; p[2] = posit[2];}
  void getPoint(Point3* p) {(*p)[0] = posit[0], (*p)[1] = posit[1]; (*p)[2] = posit[2];}
  Point3* const getPoint() { return &posit;}
  ecif_geometryType getType() {return GcPoint::type;}
  int hashKey();
  bool isBetween(GcPoint* p1, GcPoint* p2);
  bool isParallel(GcPoint* other_vector);
  void normalize(double norm[3], double shift[3]);
  void de_normalize(double norm[], double shift[]);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  double Posit(short i) const;
  double Pos(short i) const { return posit[i]; }  // Unsafe !!
  void setPosit(short i, double value);
  void setPosit(double x, double y, double z);
  void setPosit(const char* values);
  void setPos(short i, double value);
  void setPos(double x, double y, double z);
  void setPos(const char* values);
  bool operator==(const GcPoint &p2);

protected:
  Point3 posit;
  static enum ecif_geometryType type;
};//***** end class GcPoint


inline
GcPoint::GcPoint(double x, double y, double z)
{
  type = ECIF_NODIM;
  posit[X] = x;
  posit[Y] = y;
  posit[Z] = z;
}


inline
GcPoint::GcPoint(double p[3])
{
  type = ECIF_NODIM;
  posit[X] = p[0];
  posit[Y] = p[1];
  posit[Z] = p[2];
}


inline
GcPoint::GcPoint(Point3* p)
{
  type = ECIF_NODIM;
  posit[X] = (*p)[0];
  posit[Y] = (*p)[1];
  posit[Z] = (*p)[2];
}


inline
GcPoint::GcPoint(const GcPoint& p)
{
  type = ECIF_NODIM;
  posit[X] = p.posit[X];
  posit[Y] = p.posit[Y];
  posit[Z] = p.posit[Z];
}

inline
GcPoint::GcPoint(const GcPoint* p)
{
  type = ECIF_NODIM;
  posit[X] = p->posit[X];
  posit[Y] = p->posit[Y];
  posit[Z] = p->posit[Z];
}

inline double
GcPoint::Posit(short i) const
{
  if (i >= 0 && i <= MAX_DIMENSION)
    return posit[i];
  else
    return 0;
}


inline void
GcPoint::setPosit(short i, double value)
{
  if (i >= 0 && i <= MAX_DIMENSION)
    posit[i] = value;
}

inline void
GcPoint::setPosit(double x, double y, double z)
{
  posit[0] = x;
  posit[1] = y;
  posit[2] = z;
}


inline void
GcPoint::setPosit(const char* values)
{
  strstream strm;
  strm << values << ends;
  strm >> posit[0];
  strm >> posit[1];
  strm >> posit[2];
}

inline void
GcPoint::setPos(short i, double value)
{
  posit[i] = value;
}

inline void
GcPoint::setPos(double x, double y, double z)
{
  posit[0] = x;
  posit[1] = y;
  posit[2] = z;
}

inline void
GcPoint::setPos(const char* values)
{
  strstream strm;
  strm << values << ends;
  strm >> posit[0];
  strm >> posit[1];
  strm >> posit[2];
}


GcPoint operator+(const GcPoint &p1, const GcPoint &p2);
GcPoint operator-(const GcPoint &p1, const GcPoint &p2);
GcPoint operator*(int num, const GcPoint &p);
GcPoint operator/(const GcPoint &p, int num);
bool operator==(const GcPoint &p1, const GcPoint &p2);
bool operator!=(const GcPoint &p1, const GcPoint &p2);
bool operator<(const GcPoint &p1, const GcPoint &p2);
bool operator>(const GcPoint &p1, const GcPoint &p2);
bool operator<=(const GcPoint &p1, const GcPoint &p2);
bool operator>=(const GcPoint &p1, const GcPoint &p2);
ostream &operator<<(ostream &stream, const GcPoint &p);


class TransfMatrix
{
public:
  TransfMatrix();
  TransfMatrix(double a, double b, double c, double d, double e, double f,
               double g, double h, double i, double j, double k, double l);
  TransfMatrix(TransfMatrix& m);
  void transform(GcPoint& oldp, GcPoint& newp);
  TransfMatrix &operator=(const TransfMatrix& m);

private:
  double AA, BB, CC;  // rotation     A D G J     xold     xnew
  double DD, EE, FF;  //              B E H K  x  yold  =  ynew
  double GG, HH, II;  //              C F I L     zold     znew
  double JJ, KK, LL;  // translation  0 0 0 1       1        1
};//***end class TransfMatrix



// ==========================================================================
//                            2D classes
// ==========================================================================

// =====================
// Geometry data structs
// =====================

// Struct for circle data
//
struct ecif_Circle
{
  ecif_Circle();
  ~ecif_Circle();
  void init();

  bool hasVertexTies;
  bool centerGiven;
  bool radiusGiven;
  Point3 center;
  double radius;
  GcPoint* start;
  GcPoint* end;
  int nofDefiningPoints;
  Point3* definingPoints;
  int isClosed;
  double start_u;
  double end_u;
  enum linDeltaType deltaType;
  double delta;  // NOTE: This is user given discr. paramter, type is defined by deltaTyped when linearizing, in radians!
  MatcValueTable matcTable;
  bool useFixedMeshN;
};


// Struct for curve data
//
struct ecif_Curve
{
  ecif_Curve();
  ~ecif_Curve();
  void init();

  bool hasVertexTies;
  int isClosed;
  int nofPoints;
  GcPoint** points;
  int* vertex_flags;
  MatcValueTable matcTable;
};


// Struct for ellipse data
//
struct ecif_Ellipse
{
  ecif_Ellipse();
  ~ecif_Ellipse();
  void init();

  bool hasVertexTies;
  Point3 center;
  Point3 direction;
  double radius1;
  double radius2;
  int isClosed;
  GcPoint* start;
  GcPoint* end;
  MatcValueTable matcTable;
};


// Struct for hyperbola data
//
struct ecif_Hyperbola
{
  ecif_Hyperbola();
  ~ecif_Hyperbola();
  void init();

  bool hasVertexTies;
  Point3 center;
  Point3 direction;
  double radius1;
  double radius2;
  GcPoint* start;
  GcPoint* end;
  MatcValueTable matcTable;
};

// Struct for line data
struct ecif_Line
{
  ecif_Line();
  ~ecif_Line();
  void init();

  bool hasVertexTies;
  int onSymmAxis;       // 1 <--> is on symmetry-axis
  GcPoint* start;
  GcPoint* end;
  GcPoint* dir;        // NOTE: end - start (not unit vector!)
  double length;
  Point3 normal;
  MatcValueTable matcTable;
  bool useFixedMeshN;
};


struct ecif_Multi2D
{
  ecif_Multi2D();
  ~ecif_Multi2D();
  void init();

  int nofEmfComponents;
  Geometry2D** emfComponents;

  linDeltaType deltaType;
  double delta;
  bool useFixedMeshN;

  int nofComponents;
  bool* isDeletable;  // These flags 'protect' simple componets which emf-components!!!
  Geometry2D** components;

  // Possible function definition related stuff
  //
  // NOTE: Function defines either a circle(s) or
  // a curve(s). We do not store any function
  // definition in a single circle or curve.
  // So, a function based circle(s) and curve(s) always
  // becomes a Multi2D geometry!
  //
  bool hasVertexTies;
  ecif_geometryType gmtrType;
  ecif_functionType funcType;
  int funcArgC;
  double* funcArgV;
  Point3* funcStartPoint;
  Point3* funcEndPoint;
  int funcStartVertex;
  int funcEndVertex;
  char* functionName;
  char* libraryName;

  MatcValueTable matcTable;
};


// Struct for nurbs curve data
//
struct ecif_NurbsCurve
{
  ecif_NurbsCurve();
  ~ecif_NurbsCurve();
  void init();

  bool hasVertexTies;
  int isRational;
  int degree;
  int nofKnots;
  int nofCpoints;
  double* knots;
  Point4* cpoints; /* NOTE: is (nofCcpoints,4)-array*/
  GcPoint* start;
  GcPoint* end;
  enum linDeltaType deltaType;
  double delta;  // NOTE: This is used when linearizing, in radians, NOT IN USE!
  MatcValueTable matcTable;
  bool useFixedMeshN;
};


// Struct for parabola data
//
struct ecif_Parabola
{
  ecif_Parabola();
  ~ecif_Parabola();
  void init();

  bool hasVertexTies;
  Point3 apex;
  Point3 direction;
  double focalLength;
  GcPoint* start;
  GcPoint* end;
  MatcValueTable matcTable;
};


// Struct for polyline data
//
  struct ecif_PolyLine
{
  ecif_PolyLine();
  ~ecif_PolyLine();
  void init();

  bool hasVertexTies;
  int nofPoints;
  GcPoint** points;
  MatcValueTable matcTable;
  bool useFixedMeshN;
};


// Struct for spline curve data
//
struct ecif_SplineCurve
{
  ecif_SplineCurve();
  ~ecif_SplineCurve();
  void init();

  bool hasVertexTies;
  int degree;
  int nofCpoints;
  Point4* cpoints; /* NOTE: is (nofCpoints,4)-array*/
  GcPoint* start;
  GcPoint* end;
  MatcValueTable matcTable;
};



// 2D geometry classes
// ===================

class Geometry2D : public Geometry
{
friend class BodyElement;
public:
  Geometry2D();
  virtual ~Geometry2D();
  virtual void calcBoundaryPoints();
  virtual void calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u = 0.0, double delta_v = 0.0);
  virtual void calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u = 0.0);
  virtual void getBoundaryPoints(int& count, BoundaryPoint**& points);
  virtual void getDiscretizationData(int nof_components, char*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  virtual void getMifTags(int& nof_tags, int*& tags);
  virtual int getNofBoundaryPoints() { return nofBoundaryPoints; }
  virtual ecif_geometryType getType() { return ECIF_NODIM; }
  virtual ostream& output_mif(ostream& out, const char* header, bool useGmtrMeshN);
  virtual void setDiscretizationData(int nof_components, char* types, double* valuesU, double* valuesV, bool* useFixedN) {}
  virtual void setMifTag(int& next_tag) { mifTag = next_tag++; }

protected:
  virtual GcPoint* getEndPoint() { return NULL;}
  virtual GcPoint* getStartPoint() { return NULL;}
  virtual ostream& outputBoundaryPointTags(ostream& out, int indent_size);

  int nofBoundaryPoints;
  BoundaryPoint** boundaryPoints;
  int mifTag;

  static enum ecif_geometryType type;

};


class GcCircle: public Geometry2D
{
friend class BodyElement;
public:
  GcCircle(ecif_ElementComponent_X& tx, IdList& vertex_tags);
  GcCircle(Point3& center, Point3& start, Point3& stop, bool do_linearize = true);
  GcCircle(int nof_points, Point3** points, ecif_EdgeGeometry_X* params, bool do_linearize = true);
  GcCircle(BodyElement* vertex1, BodyElement* vertex2, ecif_EdgeGeometry_X* params, bool do_linearize = true);
  ~GcCircle();
  BoundBox* calcBoundBox();
  void calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta = 0.0);
  void calcLinearizingPoints(double start, double end, double delta, PointList& all_points);
  void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0);
  double getDeltaU() { return data.delta; }
  linDeltaType getDeltaType() { return data.deltaType; }
  void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  double getMeshH() { return deltaH; }
  double getParamDeltaU() { return deltaU; }
  ecif_geometryType getType() {return GcCircle::type;}
  bool isClosedU() { return data.isClosed; }
  virtual bool isDiscretized() { return true; }
  RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left);
  void isectionsWithXdir(GcPoint* startp, int& nof_values, double*& values);
  void isectionsWithYdir(GcPoint* startp, int& nof_values, double*& values);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  void param2Point(double u_p, double v_p, Point3& point);
  GcPoint* param2Point(double u_p, double v_p);
  void point2Param(Point3* point, double& u_p, double& v_p);
  ParamPair* point2Param(GcPoint* p);
  void setDeltaU(double delta_u) { deltaU = delta_u;}
  void setDeltaU(enum linDeltaType lin_delta_type, double lin_delta);
  void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN);
  bool updateGeometry(int parent_tag, IdList& vertex_tags);
  bool useFixedMeshN() { return data.useFixedMeshN; }

protected:
  bool create(int nof_points, Point3** points, ecif_EdgeGeometry_X* params,
              enum linDeltaType lin_delta_type, double lin_delta, bool do_linearize = true);
  void init(bool do_linearize = true);
  GcPoint* getEndPoint() { return data.end;}
  GcPoint* getStartPoint() { return data.start;}
  bool getLine(int index, GcLine*& line);
  bool updateData();

  ecif_Circle data;
  double deltaU;
  double deltaH;
  static enum ecif_geometryType type;

  bool param2Point_own(double u, Point3& point);
}; //***** end class GcCircle


class GcLine: public Geometry2D
{
friend class BodyElement;
public:
  GcLine(ecif_ElementComponent_X& tx, IdList& vertex_tags);
  GcLine(Point3* p1, Point3* p2);
  GcLine(BodyElement* vertex1, BodyElement* vertex2);
  ~GcLine();
  BoundBox* calcBoundBox();
  //void calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u = 0.0, double delta_v = 0.0);
  void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0);
  const GcPoint* getDir() { return data.dir;}
  void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  ecif_geometryType getType() {return GcLine::type;}
  RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left);
  void isectionsWithXdir(GcPoint* startp, int& nof_values, double*& values);
  void isectionsWithYdir(GcPoint* startp, int& nof_values, double*& values);
  bool isOnSameAxis(GcPoint& p1, GcPoint& p2);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  void param2Point(double u_p, double v_p, Point3& point);
  GcPoint* param2Point(double u_p, double v_p);
  void point2Param(Point3* point, double& u_p, double& v_p);
  ParamPair* point2Param(GcPoint* p);
  void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN);
  bool updateGeometry(int parent_tag, IdList& vertex_tags);

protected:
  void create(GcPoint* gp1, GcPoint* gp2);
  GcPoint* getEndPoint() { return data.end;}
  GcPoint* getStartPoint() { return data.start;}
  bool updateData();

  ecif_Line data;
  static enum ecif_geometryType type;

}; //***** end class GcLine


class GcMulti2D : public Geometry2D
{
friend class BodyElement;
public:
  GcMulti2D();
  GcMulti2D(ecif_Element_X& tx, IdList& vertex_tags);
  virtual ~GcMulti2D();
  virtual BoundBox* calcBoundBox();
  virtual void calcBoundaryPoints();
  virtual void calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u = 0.0);
  virtual void draw(int gmtr_index, Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0);
  virtual void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0, bool is_first_loop = true);
  virtual void getBoundaryPoints(int& count, BoundaryPoint**& points);
  virtual Geometry* getComponent(int index, bool only_emf_components = false);
  virtual linDeltaType getDeltaType();
  virtual double getDeltaU();
  virtual void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  virtual void getLabelPoint(Point3& point);
  virtual void getMifTags(int& nof_tags, int*& tags);
  virtual int getNofBoundaryPoints();
  virtual int getNofComponents(bool only_emf_components = false);
  virtual ecif_geometryType getType();
  virtual bool isClosedU();
  virtual RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left);
  virtual void isectionsWithXdir(GcPoint* startp, int& nof_values, double*& values);
  virtual void isectionsWithYdir(GcPoint* startp, int& nof_values, double*& values);
  virtual bool isMultiGeometry() {return data.nofComponents > 1; }
  virtual bool isOnSameAxis(GcPoint& p1, GcPoint& p2);
  virtual ostream& output_emf(ostream& out, short indent_size, short indent_level);
  virtual ostream& output_mif(ostream& out, const char* header, bool useGmtrMeshN);
  virtual void param2Point(double u_p, double v_p, Point3& point);
  virtual GcPoint* param2Point(double u_p, double v_p);
  virtual void point2Param(Point3* point, double& u_p, double& v_p);
  virtual ParamPair* point2Param(GcPoint* p);
  virtual void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN);
  virtual void setMifTag(int& next_tag);
  bool updateGeometry(int parent_tag, IdList& vertex_tags);
  virtual bool useFixedMeshN();

protected:
  GcMulti2D* createFunctionCircleC(int edge_tag, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  //GcMulti2D* createFunctionCircleF(int edge_tag, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  GcMulti2D* createFunctionCircleM(int edge_tag, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  GcMulti2D* createFunctionCircle(int nof_circles, egf_Circle* circles, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  GcMulti2D* createFunctionPolyLineC(int edge_tag, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  //GcMulti2D* createFunctionPolyLineF(int edge_tag, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  GcMulti2D* createFunctionPolyLineM(int edge_tag, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  GcMulti2D* createFunctionPolyLine(int nof_polylines, egf_PolyLine* polylines, ecif_ElementComponent_X& txc, IdList& vertex_tags);
  bool getMatcCircleData(const char* matc_result, int& nof_circles, egf_Circle*& circles, strstream& strm);
  bool getMatcPolyLineData(const char* matc_result, int& nof_polylines, egf_PolyLine*& polylines, strstream& strm);
  void init();
  ostream& output_1_emf(ostream& out, short indent_size, short indent_level);
  ostream& output_n_emf(ostream& out, short indent_size, short indent_level);
  void setGeometryUpdateData(ecif_ElementComponent_X& txc, ecif_Multi2D data);
  bool updateData();

  ecif_Multi2D data;
  static enum ecif_geometryType type;
}; //***** end class GcMulti2D


class GcNurbsCurve: public Geometry2D
{
friend class BodyElement;
public:
  GcNurbsCurve(ecif_ElementComponent_X& tx, IdList& vertex_tags);
  GcNurbsCurve(ecif_EdgeGeometry_X* params, bool do_linearize = true);
  GcNurbsCurve(BodyElement* vertex1, BodyElement* vertex2,
               ecif_EdgeGeometry_X* params, bool do_linearize = true);
  ~GcNurbsCurve();
  BoundBox* calcBoundBox();
  void calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta = 0.0);
  void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0);
  linDeltaType getDeltaType() { return data.deltaType; }
  double getDeltaU() { return data.delta; }
  void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  double getParamDeltaU() { return deltaU; }
  ecif_geometryType getType() {return GcNurbsCurve::type;}
  RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left);
  bool isLinear() {return false;}
  void linearize();
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  void param2Point(double u_p, double v_p, Point3& point);
  GcPoint* param2Point(double u_p, double v_p);
  void point2Param(Point3* point, double& u_p, double& v_p);
  ParamPair* point2Param(GcPoint* p);
  ParamPair* point2Param(GcPoint* p, int direction = 1,
              double start_u = 0, double start_v = 0);
  void setDeltaU(double delta_u) { deltaU = delta_u;}
  void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN);
  bool useGeometryMeshN() { return true; }

protected:
  void create(ecif_EdgeGeometry_X* params, bool do_linearize = true);
  GcPoint* getEndPoint() { return data.end;}
  GcPoint* getStartPoint() { return data.start;}

  ecif_NurbsCurve data;
  double deltaU;

  int nofLinearizingPoints;
  GcPoint** linearizingPoints;

  int nofBoundaryPoints;
  BoundaryPoint** boundaryPoints;

  static enum ecif_geometryType type;
}; //***** end class GcNurbsCurve


class GcPolyLine: public Geometry2D
{
friend class BodyElement;
public:
  GcPolyLine(ecif_ElementComponent_X& tx, IdList& vertex_tags);
  GcPolyLine(int nof_points, GcPoint** points);
  GcPolyLine(int nof_vertices, BodyElement** vertices);
  ~GcPolyLine();
  BoundBox* calcBoundBox();
  void calcLinearizingPoints(int& nof_points, GcPoint**& points, double delta_u = 0.0);
  void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0);
  void getLabelPoint(Point3& point);
  void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  ecif_geometryType getType() {return GcPolyLine::type;}
  bool isClosedU();
  RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left);
  //void isectionsWithXdir(GcPoint* startp, int& nof_values, double*& values);
  //void isectionsWithYdir(GcPoint* startp, int& nof_values, double*& values);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  void param2Point(double u_p, double v_p, Point3& point);
  ParamPair* point2Param(GcPoint* p);
  void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN);
  bool updateGeometry(int parent_tag, IdList& vertex_tags);
  bool useFixedMeshN() { return data.useFixedMeshN; }

protected:
  void create(int nof_points, GcPoint** line_points);
  GcPoint* getEndPoint() { return data.points[data.nofPoints - 1];}
  bool getLine(int index, GcLine*& line);
  int getNofLines() { return data.nofPoints - 1;}
  GcPoint* getStartPoint() { return data.points[0];}
  bool updateData();

  ecif_PolyLine data;
  static enum ecif_geometryType type;

}; //***** end class GcPolyLine



// ==========================================================================
//                            3D classes
// ==========================================================================

// =====================
// Geometry data structs
// =====================

// Struct for plane data
//
struct ecif_Plane
{
  ecif_Plane() { init();}
  ~ecif_Plane() {
    purgeMatcValueTable(matcTable);
  }

  void init() {
    hasVertexTies = false;
    onSymmPlane = 0;
  }

  bool hasVertexTies;
  int onSymmPlane;       /* 1 <--> is on symmetry-plane */
  GcPoint* corners[4];
  MatcValueTable matcTable;
};


// Struct for nurbs surface data
//
struct ecif_NurbsSurface
{
  ecif_NurbsSurface() { init();}
  ~ecif_NurbsSurface() {
    delete[] knots_u;
    delete[] knots_v;
    delete[] cpoints;
    purgeMatcValueTable(matcTable);
  }

  void init() {
    hasVertexTies = false;
    isRational = 0;
    degree_u = 0;
    degree_v = 0;
    nofKnots_u = 0;
    nofKnots_v = 0;
    nofCpoints_u = 0;
    nofCpoints_v = 0;
    nofCpoints = 0;
    knots_u = NULL;
    knots_v = NULL;
    cpoints = NULL;
  }

  bool hasVertexTies;
  int isRational;
  int degree_u;
  int degree_v;
  int nofKnots_u;
  int nofKnots_v;
  int nofCpoints_u;
  int nofCpoints_v;
  int nofCpoints;
  double* knots_u;
  double* knots_v;
  Point4* cpoints; /* NOTE: is (nofCpoints,4)-array*/
  GcPoint* corners[4];
  MatcValueTable matcTable;
};


// 3D geometry classes
// ===================

class Geometry3D : public Geometry
{
friend class BodyElement;
public:
  Geometry3D();
  virtual ~Geometry3D();
  virtual void calcBoundaryPoints(int boundary_tag);

protected:

  int nofBoundaryPointsU;
  int nofBoundaryPointsV;
  BoundaryPoint** boundaryPoints;
};


class GcNurbsSurface: public Geometry3D
{
friend class BodyElement;
public:
  GcNurbsSurface(ecif_ElementComponent_X& tx, IdList& vertex_tags);
  GcNurbsSurface(ecif_FaceGeometry_X* params, bool do_linearize = true);
  ~GcNurbsSurface();
  BoundBox* calcBoundBox();
  void calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u = 0.0, double delta_v = 0.0);
  void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0);
  ecif_geometryType getType() {return GcNurbsSurface::type;}
  bool isLinear() {return false;}
  void linearize();
  //RayHit* isectionsWithXdir(GcPoint* startp);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  //GcPoint* param2Point(double u_p, double v_p = 0);
  //void point2Param(Point3* point, double& u_p, double& v_p);
  //ParamPair* point2Param(GcPoint* p);
  //ParamPair* point2Param(GcPoint* p, int direction = 1,
            //  double start_u = 0, double start_v = 0);
  //bool rightOnNegative();
  void setDeltaU(double delta_u) { deltaU = delta_u;}
  void setDeltaV(double delta_v) { deltaV = delta_v;}

protected:
  void create(ecif_FaceGeometry_X* params, bool do_linearize = true);

  void calcLinearizingPoints( double delta_u, double delta_v);
  ecif_NurbsSurface data;
  double deltaU;
  double deltaV;

  int nofLinearizingPointsU;
  int nofLinearizingPointsV;
  GcPoint** linearizingPoints;

  static enum ecif_geometryType type;
}; //***** end class VttNurbsSrf


class GcPlane: public Geometry3D
{
friend class BodyElement;
public:
  GcPlane(ecif_ElementComponent_X& tx, IdList& vertex_tags);
  GcPlane(int nof_points, GcPoint** points);
  ~GcPlane();
  BoundBox* calcBoundBox();
  //void calcLinearizingPoints(int& nof_points_u, int& nof_points_v, GcPoint**& points, double delta_u = 0.0, double delta_v = 0.0);
  //void draw(Renderer* renderer, objectDrawingMode dmode, objectDrawingState dstate, int direction, int elem_id = 0, bool is_first_loop = true) {};
  ecif_geometryType getType() {return GcPlane::type;}
  //RayHit* isectionsWithXdir(GcPoint* startp);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  //GcPoint* param2Point(double u_p, double v_p = 0);
  //ParamPair* point2Param(GcPoint* p);
  //void point2Param(Point3* point, double& u_p, double& v_p);
  //bool rightOnNegative();
protected:
  void create(int nof_points, GcPoint** points);

  ecif_Plane data;

  static enum ecif_geometryType type;
}; //***** end class GcLine


#endif
