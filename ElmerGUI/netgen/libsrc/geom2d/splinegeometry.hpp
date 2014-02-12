/*


JS, Nov 2007


The 2D/3D template-base classes should go into the libsrc/gprim directory

in geom2d only 2D - Geometry classes (with material properties etc.)


*/





#ifndef _FILE_SPLINEGEOMETRY
#define _FILE_SPLINEGEOMETRY
#include "../csg/csgparser.hpp"

/// 
extern void LoadBoundarySplines (const char * filename,
				 ARRAY < GeomPoint<2> > & geompoints,
				 ARRAY < SplineSeg<2>* > & splines, 
				 double & elto0);
///
extern void PartitionBoundary (const ARRAY < SplineSeg<2>* > & splines,
			       double h, double elto0,
			       Mesh & mesh2d);


// allow to turn off messages: cover all couts !!
extern int printmessage_importance;

template < int D >
class SplineGeometry
{
  ARRAY < GeomPoint<D> > geompoints;
  ARRAY < SplineSeg<D>* > splines;
  double elto0;
  ARRAY<char*> materials;
  ARRAY<string*> bcnames;
  ARRAY<double> maxh;
  ARRAY<bool> quadmeshing;
  ARRAY<bool> tensormeshing;

private:
  void AppendSegment(SplineSeg<D> * spline, const int leftdomain, const int rightdomain,
		     const int bc,
		     const double reffac, const bool hprefleft, const bool hprefright,
		     const int copyfrom);

public:
  ~SplineGeometry();

  int Load (const ARRAY<double> & raw_data, const int startpos = 0);
  void Load (const char * filename);
  void CSGLoad (CSGScanner & scan);

  void LoadData( ifstream & infile );
  void LoadDataNew ( ifstream & infile );
  void LoadDataV2 ( ifstream & infile );

  void PartitionBoundary (double h, Mesh & mesh2d);

  void GetRawData (ARRAY<double> & raw_data) const;

  void CopyEdgeMesh (int from, int to, Mesh & mesh2d, Point3dTree & searchtree);

  const ARRAY<SplineSeg<D>*> & GetSplines () const
  { return splines; }

  int GetNSplines (void) const { return splines.Size(); }
  string GetSplineType (const int i) const { return splines[i]->GetType(); }
  SplineSeg<D> & GetSpline (const int i) {return *splines[i];}
  const SplineSeg<D> & GetSpline (const int i) const {return *splines[i];}

  void GetBoundingBox (Box<D> & box) const;
  Box<D> GetBoundingBox () const 
  { Box<D> box; GetBoundingBox (box); return box; }

  int GetNP () const { return geompoints.Size(); }
  const GeomPoint<D> & GetPoint(int i) const { return geompoints[i]; }

  void SetGrading (const double grading);
  void AppendPoint (const double x, const double y, const double reffac = 1., const bool hpref = false);
  void AppendPoint (const Point<D> & p, const double reffac = 1., const bool hpref = false);
  
  void AppendLineSegment (const int n1, const int n2,
			  const int leftdomain, const int rightdomain, const int bc = -1,
			  const double reffac = 1.,
			  const bool hprefleft = false, const bool hprefright = false,
			  const int copyfrom = -1);
  void AppendSplineSegment (const int n1, const int n2, const int n3,
			    const int leftdomain, const int rightdomain, const int bc = -1,
			    const double reffac = 1.,
			    const bool hprefleft = false, const bool hprefright = false,
			    const int copyfrom = -1);
  void AppendCircleSegment (const int n1, const int n2, const int n3,
			    const int leftdomain, const int rightdomain, const int bc = -1,
			    const double reffac = 1.,
			    const bool hprefleft = false, const bool hprefright = false,
			    const int copyfrom = -1);
  void AppendDiscretePointsSegment (const ARRAY< Point<D> > & points, 
				    const int leftdomain, const int rightdomain, const int bc = -1,
				    const double reffac = 1.,
				    const bool hprefleft = false, const bool hprefright = false,
				    const int copyfrom = -1);
  void TestComment ( ifstream & infile ) ;
  void	GetMaterial( const int  domnr, char* & material );

  double GetDomainMaxh ( const int domnr );
  bool GetDomainQuadMeshing ( int domnr ) 
  { 
    if ( quadmeshing.Size() ) return quadmeshing[domnr-1]; 
    else return false;
  }
  bool GetDomainTensorMeshing ( int domnr ) 
  { 
    if ( tensormeshing.Size() ) return tensormeshing[domnr-1]; 
    else return false;
  }

  string GetBCName ( const int bcnr ) const;

  string * BCNamePtr ( const int bcnr );
};


void MeshFromSpline2D (SplineGeometry<2> & geometry,
		       Mesh *& mesh, 
		       MeshingParameters & mp);



typedef SplineGeometry<2> SplineGeometry2d;


#endif // _FILE_SPLINEGEOMETRY
