#ifndef CURVEDELEMS_NEWH
#define CURVEDELEMS_NEWH

/**************************************************************************/
/* File:   curvedelems.hpp                                                */
/* Author: Robert Gaisbauer (first version)                               */
/*         redesign by Joachim Schoeberl                                  */
/* Date:   27. Sep. 02, Feb 2006                                          */
/**************************************************************************/




class Refinement;


class CurvedElements
{
  const Mesh & mesh;

  ARRAY<int> edgeorder;
  ARRAY<int> faceorder;

  ARRAY<int> edgecoeffsindex;
  ARRAY<int> facecoeffsindex;

  ARRAY< Vec<3> > edgecoeffs;
  ARRAY< Vec<3> > facecoeffs;

  ARRAY< double > edgeweight;  // for rational 2nd order splines

  int order;
  bool rational;

  bool ishighorder;

public:
  CurvedElements (const Mesh & amesh);
  ~CurvedElements();

  // bool IsHighOrder() const { return order > 1; }
  bool IsHighOrder() const { return ishighorder; }

  void SetHighOrder (int aorder) { order=aorder; }
  
  void BuildCurvedElements(Refinement * ref, int aorder, bool arational = false);

  int GetOrder () { return order; }


  bool IsSegmentCurved (SegmentIndex segnr) const;
  bool IsSurfaceElementCurved (SurfaceElementIndex sei) const;
  bool IsElementCurved (ElementIndex ei) const;


  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x)
  { CalcSegmentTransformation (xi, segnr, &x, NULL); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Vec<3> & dxdxi)
  { CalcSegmentTransformation (xi, segnr, NULL, &dxdxi); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x, Vec<3> & dxdxi)
  { CalcSegmentTransformation (xi, segnr, &x, &dxdxi, NULL); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x, Vec<3> & dxdxi, bool & curved)
  { CalcSegmentTransformation (xi, segnr, &x, &dxdxi, &curved); };



  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x)
  { CalcSurfaceTransformation (xi, elnr, &x, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x, Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x, Mat<3,2> & dxdxi, bool & curved)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi, &curved); };





  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x)
  { CalcElementTransformation (xi, elnr, &x, NULL); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x, Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , NULL */ ); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x, Mat<3,3> & dxdxi,
                                  void * buffer, bool valid)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi, /* NULL, */ buffer, valid ); };

  // void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
  // 				  Point<3> & x, Mat<3,3> & dxdxi) // , bool & curved)
  //   { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , &curved * ); }



  void CalcMultiPointSegmentTransformation (ARRAY<double> * xi, SegmentIndex segnr,
					    ARRAY<Point<3> > * x,
					    ARRAY<Vec<3> > * dxdxi);

  void CalcMultiPointSurfaceTransformation (ARRAY< Point<2> > * xi, SurfaceElementIndex elnr,
					    ARRAY< Point<3> > * x,
					    ARRAY< Mat<3,2> > * dxdxi);

  void CalcMultiPointElementTransformation (ARRAY< Point<3> > * xi, ElementIndex elnr,
					    ARRAY< Point<3> > * x,
					    ARRAY< Mat<3,3> > * dxdxi);

  void CalcMultiPointElementTransformation (ElementIndex elnr, int n,
                                            const double * xi, int sxi,
                                            double * x, int sx,
                                            double * dxdxi, int sdxdxi);




private:
  
  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> * x = NULL, Vec<3> * dxdxi = NULL, bool * curved = NULL);

  void CalcSurfaceTransformation (Point<2> xi, SurfaceElementIndex elnr,
				  Point<3> * x = NULL, Mat<3,2> * dxdxi = NULL, bool * curved = NULL);

  void CalcElementTransformation (Point<3> xi, ElementIndex elnr,
				  Point<3> * x = NULL, Mat<3,3> * dxdxi = NULL, // bool * curved = NULL,
                                  void * buffer = NULL, bool valid = 0);






  class SegmentInfo
  {
  public:
    SegmentIndex elnr;
    int order;
    int nv;
    int ndof;
    int edgenr;
  };

  void CalcElementShapes (SegmentInfo &  elnr, double xi, Vector & shapes) const;
  void GetCoefficients (SegmentInfo & elnr, ARRAY<Vec<3> > & coefs) const;
  void CalcElementDShapes (SegmentInfo & elnr, double xi, Vector & dshapes) const;


  class ElementInfo
  {
  public:
    ElementIndex elnr;
    int order;
    int nv;
    int ndof;
    int nedges;
    int nfaces;
    int edgenrs[12];
    int facenrs[6];
    Mat<3> hdxdxi;
    Vec<3> hcoefs[10]; // enough for second order tets
  };


  void CalcElementShapes (ElementInfo & info, const Point<3> & xi, Vector & shapes) const;
  void GetCoefficients (ElementInfo & info, Vec<3> * coefs) const;
  void CalcElementDShapes (ElementInfo & info, const Point<3> & xi, MatrixFixWidth<3> & dshapes) const;

  
  class SurfaceElementInfo
  {
  public:
    SurfaceElementIndex elnr;
    int order;
    int nv;
    int ndof;
    ArrayMem<int,4> edgenrs;
    int facenr;
  };

  void CalcElementShapes (SurfaceElementInfo & elinfo, const Point<2> & xi, Vector & shapes) const;
  void GetCoefficients (SurfaceElementInfo & elinfo, ARRAY<Vec<3> > & coefs) const;
  void CalcElementDShapes (SurfaceElementInfo & elinfo, const Point<2> & xi, DenseMatrix & dshapes) const;
};



#endif
