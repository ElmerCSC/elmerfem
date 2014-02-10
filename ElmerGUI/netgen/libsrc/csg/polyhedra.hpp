#ifndef FILE_POLYHEDRA
#define FILE_POLYHEDRA


/**************************************************************************/
/* File:   polyhedra.hh                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Mar. 2000                                                  */
/**************************************************************************/

/*

  Polyhedral primitive
  
*/

class Polyhedra : public Primitive
{
  class Face {
  public:
    int pnums[3];
    int planenr;

    int inputnr;

    Box<3> bbox;
    //    Point<3> center;
    Vec<3> v1, v2;   // edges
    Vec<3> w1, w2;   // pseudo-inverse
    Vec<3> n;        // normal to face
    Vec<3> nn;       // normed normal

    Face () { ; }
    Face (int pi1, int pi2, int pi3, 
	  const ARRAY<Point<3> > & points,
	  int ainputnr);
  };

  ARRAY<Point<3> > points;
  ARRAY<Face> faces;
  ARRAY<Plane*> planes;
  Box<3> poly_bbox;

  double eps_base1;

public:
  Polyhedra ();
  virtual ~Polyhedra ();
  static Primitive * CreateDefault ();

  virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
  virtual INSOLID_TYPE PointInSolid (const Point<3> & p,
				     double eps) const;
  virtual INSOLID_TYPE VecInSolid (const Point<3> & p,
				   const Vec<3> & v,
				   double eps) const;

  // checks if lim s->0 lim t->0  p + t(v1 + s v2) in solid
  virtual INSOLID_TYPE VecInSolid2 (const Point<3> & p,
				    const Vec<3> & v1,
				    const Vec<3> & v2,
				    double eps) const;

  virtual void GetTangentialSurfaceIndices (const Point<3> & p, 
					    ARRAY<int> & surfind, double eps) const;


  virtual void GetTangentialVecSurfaceIndices2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
						ARRAY<int> & surfind, double eps) const;

  virtual void CalcSpecialPoints (ARRAY<Point<3> > & pts) const;
  virtual void AnalyzeSpecialPoint (const Point<3> & pt, 
				    ARRAY<Point<3> > & specpts) const;
  virtual Vec<3> SpecialPointTangentialVector (const Point<3> & p, int s1, int s2) const;

  virtual int GetNSurfaces() const 
    { return planes.Size(); }
  virtual Surface & GetSurface (int i) 
    { return *planes[i]; }
  virtual const Surface & GetSurface (int i) const
    { return *planes[i]; }

  virtual void GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const;
  virtual void SetPrimitiveData (ARRAY<double> & coeffs);

  virtual void Reduce (const BoxSphere<3> & box);
  virtual void UnReduce ();

  int AddPoint (const Point<3> & p);
  int AddFace (int pi1, int pi2, int pi3, int inputnum);

  void GetPolySurfs(ARRAY < ARRAY<int> * > & polysurfs);
  
protected:
  int FaceBoxIntersection (int fnr, const BoxSphere<3> & box) const;
  //  void CalcData();
};

#endif
