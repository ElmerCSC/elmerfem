#ifndef FILE_BRICK
#define FILE_BRICK


/**************************************************************************/
/* File:   brick.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   11. Mar. 98                                                    */
/**************************************************************************/

/*

  brick geometry, has several surfaces
  
*/



class Parallelogram3d : public Surface
{
  Point<3> p1, p2, p3, p4;
  Vec<3> v12, v13;
  Vec<3> n;

public:
  Parallelogram3d (Point<3> ap1, Point<3> ap2, Point<3> ap3);
  virtual ~Parallelogram3d ();

  void SetPoints (Point<3> ap1, Point<3> ap2, Point<3> ap3);

  virtual int IsIdentic (const Surface & s2, int & inv, double eps) const; 

  virtual double CalcFunctionValue (const Point<3> & point) const;
  virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
  virtual void CalcHesse (const Point<3> & point, Mat<3> & hesse) const;
  virtual double HesseNorm () const;

  virtual Point<3> GetSurfacePoint () const;
  virtual void Print (ostream & str) const;
  
  virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					 const Box<3> & boundingbox, 
					 double facets) const;

protected:
  void CalcData();
};


class Brick : public Primitive
{
  Point<3> p1, p2, p3, p4;
  Vec<3> v12, v13, v14;
  // ARRAY<OneSurfacePrimitive*> faces;
  ARRAY<Plane*> faces;

public:
  Brick (Point<3> ap1, Point<3> ap2, Point<3> ap3, Point<3> ap4);
  virtual ~Brick ();
  static Primitive * CreateDefault ();

  virtual Primitive * Copy () const;
  virtual void Transform (Transformation<3> & trans);


  virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;

  virtual INSOLID_TYPE PointInSolid (const Point<3> & p,
				     double eps) const;
  virtual INSOLID_TYPE VecInSolid (const Point<3> & p,
				   const Vec<3> & v,
				   double eps) const;
  virtual INSOLID_TYPE VecInSolid2 (const Point<3> & p,
				    const Vec<3> & v1,
				    const Vec<3> & v2,
				    double eps) const;

  virtual INSOLID_TYPE VecInSolid3 (const Point<3> & p,
				    const Vec<3> & v1,
				    const Vec<3> & v2,
				    double eps) const;

  virtual INSOLID_TYPE VecInSolid4 (const Point<3> & p,
				    const Vec<3> & v,
				    const Vec<3> & v2,
				    const Vec<3> & m,
				    double eps) const;


  virtual int GetNSurfaces() const 
    { return 6; }
  virtual Surface & GetSurface (int i) 
    { return *faces[i]; }
  virtual const Surface & GetSurface (int i) const
    { return *faces[i]; }


  virtual void GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const;
  virtual void SetPrimitiveData (ARRAY<double> & coeffs);

  virtual void Reduce (const BoxSphere<3> & box);
  virtual void UnReduce ();

protected:
  void CalcData();
};


class OrthoBrick : public Brick 
{
protected:
  Point<3> pmin, pmax;
public:
  OrthoBrick (const Point<3> & ap1, const Point<3> & ap2);
  
  virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
  virtual void Reduce (const BoxSphere<3> & box);
};

#endif
