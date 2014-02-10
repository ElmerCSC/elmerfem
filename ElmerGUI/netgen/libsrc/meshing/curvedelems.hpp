#ifndef CURVEDELEMS
#define CURVEDELEMS

/**************************************************************************/
/* File:   curvedelems.hpp                                                */
/* Author: Robert Gaisbauer                                               */
/* Date:   27. Sep. 02 (second version: 30. Jan. 03)                      */
/**************************************************************************/

#include "bisect.hpp"
#include <iostream>

#define EPSILON 1e-20



void ComputeGaussRule (int n, ARRAY<double> & xi, ARRAY<double> & wi);





// ----------------------------------------------------------------------------
//      CurvedElements
// ----------------------------------------------------------------------------

class CurvedElements
{
  const Mesh & mesh;
  const MeshTopology & top;

  bool isHighOrder;
  int nvisualsubsecs;
  int nIntegrationPoints;

  ARRAY<int> edgeorder;
  ARRAY<int> faceorder;

  /*

  ARRAY< Vec<3> > edgecoeffs;
  ARRAY< Vec<3> > facecoeffs;

  ARRAY<int> edgecoeffsindex;
  ARRAY<int> facecoeffsindex;

  */

  inline Vec<3> GetEdgeCoeff (int edgenr, int k);
  inline Vec<3> GetFaceCoeff (int facenr, int k);

  
  void CalcSegmentTransformation (double xi, int segnr,
				  Point<3> * x = NULL, Vec<3> * dxdxi = NULL);

  void CalcSurfaceTransformation (Point<2> xi, int elnr,
				  Point<3> * x = NULL, Mat<3,2> * dxdxi = NULL);

  void CalcElementTransformation (Point<3> xi, int elnr,
				  Point<3> * x = NULL, Mat<3,3> * dxdxi = NULL);



public:

  Refinement * refinement;

  ARRAY< Vec<3> > edgecoeffs;
  ARRAY< Vec<3> > facecoeffs;

  ARRAY<int> edgecoeffsindex;
  ARRAY<int> facecoeffsindex;





  CurvedElements (const Mesh & amesh);
  ~CurvedElements();

  bool IsHighOrder() const
  { return isHighOrder; };
  void SetHighOrder () { isHighOrder = 1; }


  int GetNVisualSubsecs() const
  { return nvisualsubsecs; };

  const class Mesh & GetMesh() const
  { return mesh; };

  void BuildCurvedElements(Refinement * ref, int polydeg, bool rational=false);

  int GetEdgeOrder (int edgenr) const
  { return edgeorder[edgenr]; };

  int GetFaceOrder (int facenr) const
  { return faceorder[facenr]; };

  int IsEdgeCurved (int edgenr) const;

  int IsFaceCurved (int facenr) const;

  int IsSurfaceElementCurved (int elnr) const;

  int IsElementCurved (int elnr) const;


  void CalcSegmentTransformation (double xi, int segnr,
				  Point<3> & x)
  { CalcSegmentTransformation (xi, segnr, &x, NULL); };

  void CalcSegmentTransformation (double xi, int segnr,
				  Vec<3> & dxdxi)
  { CalcSegmentTransformation (xi, segnr, NULL, &dxdxi); };

  void CalcSegmentTransformation (double xi, int segnr,
				  Point<3> & x, Vec<3> & dxdxi)
  { CalcSegmentTransformation (xi, segnr, &x, &dxdxi); };


  void CalcSurfaceTransformation (const Point<2> & xi, int elnr,
				  Point<3> & x)
  { CalcSurfaceTransformation (xi, elnr, &x, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, int elnr,
				  Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcSurfaceTransformation (const Point<2> & xi, int elnr,
				  Point<3> & x, Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi); };


  void CalcElementTransformation (const Point<3> & xi, int elnr,
				  Point<3> & x)
  { CalcElementTransformation (xi, elnr, &x, NULL); };

  void CalcElementTransformation (const Point<3> & xi, int elnr,
				  Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcElementTransformation (const Point<3> & xi, int elnr,
				  Point<3> & x, Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi); };




  void CalcMultiPointSegmentTransformation (ARRAY<double> * xi, int segnr,
					    ARRAY<Point<3> > * x,
					    ARRAY<Vec<3> > * dxdxi);

  void CalcMultiPointSurfaceTransformation (ARRAY< Point<2> > * xi, int elnr,
					    ARRAY< Point<3> > * x,
					    ARRAY< Mat<3,2> > * dxdxi);

  void CalcMultiPointElementTransformation (ARRAY< Point<3> > * xi, int elnr,
					    ARRAY< Point<3> > * x,
					    ARRAY< Mat<3,3> > * dxdxi);

};



// ----------------------------------------------------------------------------
//      PolynomialBasis
// ----------------------------------------------------------------------------

class PolynomialBasis
{
  int order;
  int maxorder;
  ArrayMem<double,20> f;
  ArrayMem<double,20> df;
  ArrayMem<double,20> ddf;

  ArrayMem<double,20> lp;
  ArrayMem<double,20> dlp;

  inline void CalcLegendrePolynomials (double x);
  // P_i(x/t) t^i
  inline void CalcScaledLegendrePolynomials (double x, double t);
  inline void CalcDLegendrePolynomials (double x);

public:

  PolynomialBasis ()
  { maxorder = -1; };

  ~PolynomialBasis ()
  {};

  void SetOrder (int aorder)
  {
    order = aorder;
    if (order > maxorder)
      {
	maxorder = order;
	f.SetSize(order-1);
	df.SetSize(order-1);
	ddf.SetSize(order-1);
	lp.SetSize(order+1);
	dlp.SetSize(order);
      };
  };

  inline void CalcF (double x);
  inline void CalcDf (double x);
  inline void CalcDDf (double x);

  inline void CalcFDf (double x);

  // compute F_i(x/t) t^i
  inline void CalcFScaled (double x, double t);
  static inline void CalcFScaled (int p, double x, double t, double * values);

  double GetF (int p) { return f[p-2]; };
  double GetDf (int p) { return df[p-2]; };
  double GetDDf (int p) { return ddf[p-2]; };
};



// ----------------------------------------------------------------------------
//      BaseFiniteElement
// ----------------------------------------------------------------------------

template <int DIM>
class BaseFiniteElement
{
protected:

  Point<DIM> xi;
  int elnr;
  const CurvedElements & curv;
  const Mesh & mesh;
  const MeshTopology & top;

public:

  BaseFiniteElement(const CurvedElements & acurv)
    : curv(acurv), mesh(curv.GetMesh()), top(mesh.GetTopology())
  {};

  virtual ~BaseFiniteElement()
  {};

  void SetElementNumber (int aelnr)
  { elnr = aelnr; }; // 1-based arrays in netgen

  virtual void SetReferencePoint (Point<DIM> axi)
  { xi = axi; };
};



// ----------------------------------------------------------------------------
//      BaseFiniteElement1D
// ----------------------------------------------------------------------------

class BaseFiniteElement1D : public BaseFiniteElement<1>
{
protected:
  PolynomialBasis b;

  int vertexnr[2];
  int edgenr;
  int edgeorient;
  int edgeorder;

  int maxedgeorder;

  double vshape[2];
  double vdshape[2];
  ArrayMem<double,20> eshape;
  ArrayMem<double,20> edshape;
  ArrayMem<double,20> eddshape;

public:

  BaseFiniteElement1D (const CurvedElements & acurv) : BaseFiniteElement<1>(acurv)
  { maxedgeorder = 1; };

  virtual ~BaseFiniteElement1D()
  {};

  int GetVertexNr (int v)
  { return vertexnr[v]; };

  int GetEdgeNr ()
  { return edgenr; };

  int GetEdgeOrder ()
  { return edgeorder; };

  int GetEdgeOrientation ()
  { return edgeorient; };

  void CalcVertexShapes();
  void CalcEdgeShapes();
  void CalcEdgeLaplaceShapes();

  double GetVertexShape (int v)
  { return vshape[v]; };

  double GetEdgeShape (int index)
  { return eshape[index]; };

  double GetVertexDShape (int v)
  { return vdshape[v]; };

  double GetEdgeDShape (int index)
  { return edshape[index]; };

  double GetEdgeLaplaceShape (int index)
  { return eddshape[index]; };

};




// ----------------------------------------------------------------------------
//      FESegm
// ----------------------------------------------------------------------------

class FESegm : public BaseFiniteElement1D
{

public:

  FESegm(const CurvedElements & acurv) : BaseFiniteElement1D(acurv)
  {};

  virtual ~FESegm()
  {};

  void SetElementNumber (int aelnr)
  { 
    BaseFiniteElement<1> :: SetElementNumber (aelnr);
    Segment s = mesh.LineSegment(elnr);
    vertexnr[0] = s.p1;
    vertexnr[1] = s.p2;
    edgenr = top.GetSegmentEdge(elnr);
    edgeorient = top.GetSegmentEdgeOrientation(elnr);
    edgeorder = curv.GetEdgeOrder(edgenr-1); // 1-based arrays in netgen

    if (edgeorder > maxedgeorder)
      {
	maxedgeorder = edgeorder;
	eshape.SetSize(maxedgeorder-1);
	edshape.SetSize(maxedgeorder-1);
	eddshape.SetSize(maxedgeorder-1);
      }
  };

};



// ----------------------------------------------------------------------------
//      FEEdge
// ----------------------------------------------------------------------------

class FEEdge : public BaseFiniteElement1D
{

public:

  FEEdge(const CurvedElements & acurv) : BaseFiniteElement1D(acurv)
  {};

  virtual ~FEEdge()
  {};

  void SetElementNumber (int aelnr)
  { 
    BaseFiniteElement<1> :: SetElementNumber (aelnr);
    top.GetEdgeVertices (elnr, vertexnr[0], vertexnr[1]);
    edgenr = elnr;
    edgeorient = 1;
    edgeorder = curv.GetEdgeOrder(edgenr-1); // 1-based arrays in netgen

    if (edgeorder > maxedgeorder)
      {
	maxedgeorder = edgeorder;
	eshape.SetSize(maxedgeorder-1);
	edshape.SetSize(maxedgeorder-1);
	eddshape.SetSize(maxedgeorder-1);
      }
  };
    
};



// ----------------------------------------------------------------------------
//      BaseFiniteElement2D
// ----------------------------------------------------------------------------

class BaseFiniteElement2D : public BaseFiniteElement<2>
{
protected:

  int nvertices;
  int nedges;

  int vertexnr[4];
  int edgenr[4];
  int edgeorient[4];
  int edgeorder[4];
  int facenr;
  int faceorient;
  int faceorder;
 
  int nfaceshapes;

  int maxedgeorder;
  int maxfaceorder;

  PolynomialBasis b1, b2;

  double vshape[4];
  Vec<2> vdshape[4];
  ArrayMem<double,80> eshape;
  ArrayMem< Vec<2>,80> edshape;
  ArrayMem<double,400> fshape;
  ArrayMem<Vec<2>,400> fdshape;
  ArrayMem<double,400> fddshape;

  virtual void CalcNFaceShapes () = 0;

public:

  BaseFiniteElement2D (const CurvedElements & acurv) : BaseFiniteElement<2>(acurv)
  { maxedgeorder = maxfaceorder = -1; };

    virtual ~BaseFiniteElement2D()
	{};

  void SetElementNumber (int aelnr);

  virtual void SetVertexSingularity (int v, int exponent) = 0;

  int GetVertexNr (int v)
  { return vertexnr[v]; };

  int GetEdgeNr (int e)
  { return edgenr[e]; };

  int GetFaceNr ()
  { return facenr; };

  int GetEdgeOrder (int e)
  { return edgeorder[e]; };

  int GetFaceOrder ()
  { return faceorder; }

  int GetNVertices ()
  { return nvertices; };

  int GetNEdges ()
  { return nedges; };

  int GetNFaceShapes ()
  { return nfaceshapes; };

  int IsCurved ()
  {
    bool iscurved = 0;
    int e;

    for (e = 0; e < GetNEdges(); e++)
      iscurved = iscurved || (GetEdgeOrder(e) > 1);

    return iscurved || (GetFaceOrder() > 1);
  }

  virtual void CalcVertexShapes() = 0;
  virtual void CalcEdgeShapes() = 0; 
  virtual void CalcFaceShapes() = 0;

  virtual void CalcFaceLaplaceShapes() = 0;

  double GetVertexShape (int v)
  { return vshape[v]; };

  double GetEdgeShape (int index)
  { return eshape[index]; };

  double GetFaceShape (int index)
  { return fshape[index]; };

  Vec<2> GetVertexDShape (int v)
  { return vdshape[v]; };

  Vec<2> GetEdgeDShape (int index)
  { return edshape[index]; };

  Vec<2> GetFaceDShape (int index)
  { return fdshape[index]; };

  double GetFaceLaplaceShape (int index)
  { return fddshape[index]; };
};



// ----------------------------------------------------------------------------
//      FETrig
// ----------------------------------------------------------------------------

class FETrig : public BaseFiniteElement2D
{
  Point<3> lambda;
  Mat<3,2> dlambda;

  const ELEMENT_EDGE * eledge;
  const ELEMENT_FACE * elface;

  virtual void CalcNFaceShapes ()
  { nfaceshapes = ((faceorder-1)*(faceorder-2))/2; };

public:

  FETrig (const CurvedElements & acurv) : BaseFiniteElement2D(acurv)
  {
    nvertices = 3;
    nedges = 3;
    eledge = MeshTopology :: GetEdges (TRIG);
    elface = MeshTopology :: GetFaces (TRIG);
  };

    virtual ~FETrig()
	{};

  virtual void SetReferencePoint (Point<2> axi);

  virtual void SetVertexSingularity (int v, int exponent);

  virtual void CalcVertexShapes();
  virtual void CalcEdgeShapes();
  virtual void CalcFaceShapes();

  virtual void CalcFaceLaplaceShapes();
};



// ----------------------------------------------------------------------------
//      FEQuad
// ----------------------------------------------------------------------------

class FEQuad : public BaseFiniteElement2D
{
  const ELEMENT_FACE * elface;

  virtual void CalcNFaceShapes ()
  { nfaceshapes = (faceorder-1)*(faceorder-1); };

public:

  FEQuad (const CurvedElements & acurv) : BaseFiniteElement2D(acurv)
  {
    nvertices = 4;
    nedges = 4;
    elface = MeshTopology :: GetFaces (QUAD);
  };

    virtual ~FEQuad()
	{};

  virtual void SetVertexSingularity (int /* v */, int /* exponent */)
	{};

  virtual void CalcVertexShapes();
  virtual void CalcEdgeShapes();
  virtual void CalcFaceShapes();

  virtual void CalcFaceLaplaceShapes();
};




// ----------------------------------------------------------------------------
//      BaseFiniteElement3D
// ----------------------------------------------------------------------------

class BaseFiniteElement3D : public BaseFiniteElement<3>
{
protected:

  int nvertices;
  int nedges;
  int nfaces;

  int vertexnr[8];
  int edgenr[12];
  int edgeorient[12];
  int edgeorder[12];
  int facenr[6];
  int faceorient[6];
  int faceorder[6];
  int surfacenr[6];
  // int surfaceorient[6];

  int nfaceshapes[6];

  int maxedgeorder;
  int maxfaceorder;

  PolynomialBasis b1, b2;

  double vshape[8];
  Vec<3> vdshape[8];
  ArrayMem<double,120> eshape;
  ArrayMem<Vec<3>,120> edshape;
  ArrayMem<double,300> fshape;
  ArrayMem<Vec<3>,300> fdshape;

  virtual void CalcNFaceShapes () = 0;

public:

  int locmaxedgeorder;
  int locmaxfaceorder;

  BaseFiniteElement3D (const CurvedElements & acurv) : BaseFiniteElement<3>(acurv)
  { maxedgeorder = maxfaceorder = -1; };

  void SetElementNumber (int aelnr);

  int GetVertexNr (int v)
  { return vertexnr[v]; };

  int GetEdgeNr (int e)
  { return edgenr[e]; };

  int GetFaceNr (int f)
  { return facenr[f]; };

  int GetNFaceShapes (int f)
  { return nfaceshapes[f]; };

  int GetEdgeOrder (int e)
  { return edgeorder[e]; };

  int GetFaceOrder (int f)
  { return faceorder[f]; };

  int GetNVertices ()
  { return nvertices; };

  int GetNEdges ()
  { return nedges; };

  int GetNFaces ()
  { return nfaces; };

  int IsCurved ()
  {
    bool iscurved = 0;
    int e, f;

    for (e = 0; e < GetNEdges(); e++)
      iscurved = iscurved || (GetEdgeOrder(e) > 1);

    for (f = 0; f < GetNFaces(); f++)
      iscurved = iscurved || (GetFaceOrder(f) > 1);

    return iscurved;
  }

  virtual void CalcVertexShapes() = 0;
  virtual void CalcVertexShapesOnly()
  { CalcVertexShapes(); }

  virtual void CalcEdgeShapes() = 0;
  virtual void CalcEdgeShapesOnly() 
    { CalcEdgeShapes(); }

  virtual void CalcFaceShapes() = 0;

  double GetVertexShape (int v)
  { return vshape[v]; };

  double GetEdgeShape (int index)
  { return eshape[index]; };

  double GetFaceShape (int index)
  { return fshape[index]; };

  Vec<3> GetVertexDShape (int v)
  { return vdshape[v]; };

  Vec<3> GetEdgeDShape (int index)
  { return edshape[index]; };

  Vec<3> GetFaceDShape (int index)
  { return fdshape[index]; };
};



// ----------------------------------------------------------------------------
//      FETet
// ----------------------------------------------------------------------------

class FETet : public BaseFiniteElement3D
{
  Point<4> lambda;
  Mat<4,3> dlambda;

  const ELEMENT_EDGE * eledge;
  const ELEMENT_FACE * elface;

  virtual void CalcNFaceShapes ()
  {
    for (int f = 0; f < nfaces; f++)
      nfaceshapes[f] = ((faceorder[f]-1)*(faceorder[f]-2))/2;
  };

public:

  FETet (const CurvedElements & acurv) : BaseFiniteElement3D(acurv)
  {
    nvertices = 4;
    nedges = 6;
    nfaces = 4;
    eledge = MeshTopology :: GetEdges (TET);
    elface = MeshTopology :: GetFaces (TET);
  };

  void SetReferencePoint (Point<3> axi);

  virtual void CalcVertexShapes();
  virtual void CalcVertexShapesOnly();
  virtual void CalcEdgeShapes();
  virtual void CalcEdgeShapesOnly();
  virtual void CalcFaceShapes();
};



// ----------------------------------------------------------------------------
//      FEPrism
// ----------------------------------------------------------------------------

class FEPrism : public BaseFiniteElement3D
{
  Point<4> lambda;   // mixed barycentric coordinates
  Mat<4,3> dlambda;

  const ELEMENT_EDGE * eledge;
  const ELEMENT_FACE * elface;

  virtual void CalcNFaceShapes ()
  {
    int f;
    for (f = 0; f < 2; f++)
      nfaceshapes[f] = ((faceorder[f]-1)*(faceorder[f]-2))/2;
    for (f = 2; f < nfaces; f++)
      nfaceshapes[f] = (faceorder[f]-1)*(faceorder[f]-1);
  };

public:

  FEPrism (const CurvedElements & acurv) : BaseFiniteElement3D(acurv)
  {
    nvertices = 6;
    nedges = 9;
    nfaces = 5;
    eledge = MeshTopology :: GetEdges (PRISM);
    elface = MeshTopology :: GetFaces (PRISM);
  };

  void SetReferencePoint (Point<3> axi);

  virtual void CalcVertexShapes();
  virtual void CalcEdgeShapes();
  virtual void CalcFaceShapes();
};




// ----------------------------------------------------------------------------
//      FEPyramid
// ----------------------------------------------------------------------------

class FEPyramid : public BaseFiniteElement3D
{

  const ELEMENT_EDGE * eledge;
  const ELEMENT_FACE * elface;

  virtual void CalcNFaceShapes ()
  {
    int f;
    for (f = 0; f < 4; f++)
      nfaceshapes[f] = ((faceorder[f]-1)*(faceorder[f]-2))/2;
    for (f = 4; f < nfaces; f++)
      nfaceshapes[f] = (faceorder[f]-1)*(faceorder[f]-1);
  };

public:

  FEPyramid (const CurvedElements & acurv) : BaseFiniteElement3D(acurv)
  {
    nvertices = 5;
    nedges = 8;
    nfaces = 5;
    eledge = MeshTopology :: GetEdges (PYRAMID);
    elface = MeshTopology :: GetFaces (PYRAMID);
  };

  void SetReferencePoint (Point<3> axi);

  virtual void CalcVertexShapes();
  virtual void CalcEdgeShapes();
  virtual void CalcFaceShapes();
};




// ----------------------------------------------------------------------------
//      FEHex
// ----------------------------------------------------------------------------

class FEHex : public BaseFiniteElement3D
{

  const ELEMENT_EDGE * eledge;
  const ELEMENT_FACE * elface;

  virtual void CalcNFaceShapes ()
  {
    int f;
    for (f = 0; f < 6; f++)
      nfaceshapes[f] = (faceorder[f]-1)*(faceorder[f]-1);
  };

public:

  FEHex (const CurvedElements & acurv) : BaseFiniteElement3D(acurv)
  {
    nvertices = 8;
    nedges = 12;
    nfaces = 6;
    eledge = MeshTopology :: GetEdges (HEX);
    elface = MeshTopology :: GetFaces (HEX);
  };

  void SetReferencePoint (Point<3> axi);

  virtual void CalcVertexShapes();
  virtual void CalcEdgeShapes();
  virtual void CalcFaceShapes();
};




#endif
