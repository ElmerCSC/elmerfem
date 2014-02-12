#ifndef _EXTRUSION_HPP
#define _EXTRUSION_HPP


class Extrusion;

class ExtrusionFace : public Surface
{
private:
  const SplineSeg<2> * profile;
  const SplineGeometry<3> * path;
  Vec<3> glob_z_direction;

  bool deletable;
  
  ARRAY< const SplineSeg3<3> * > spline3_path;
  ARRAY< const LineSeg<3> * > line_path;
  
  mutable ARRAY < Vec<3> > x_dir, y_dir, z_dir, loc_z_dir;
  mutable ARRAY < Point<3> > p0;

  mutable Vec<3> profile_tangent;
  mutable double profile_par;
  
  mutable Vector profile_spline_coeff;

  mutable int latest_seg;
  mutable double latest_t;
  mutable Point<2> latest_point2d;
  mutable Point<3> latest_point3d;


private:
  void Orthogonalize(const Vec<3> & v1, Vec<3> & v2) const;

  void Init(void);

public:
  double CalcProj(const Point<3> & point3d, Point<2> & point2d,
		  const int seg) const;
  void CalcProj(const Point<3> & point3d, Point<2> & point2d,
		int & seg, double & t) const;

public:
  ExtrusionFace(const SplineSeg<2> * profile_in,
		const SplineGeometry<3> * path_in,
		const Vec<3> & z_direction);

  ExtrusionFace(const ARRAY<double> & raw_data);
    

  ~ExtrusionFace();
  
  virtual int IsIdentic (const Surface & s2, int & inv, double eps) const; 
  
  virtual double CalcFunctionValue (const Point<3> & point) const;
  virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
  virtual void CalcHesse (const Point<3> & point, Mat<3> & hesse) const;
  virtual double HesseNorm () const;

  virtual double MaxCurvature () const;
  //virtual double MaxCurvatureLoc (const Point<3> & /* c */ , 
  //				  double /* rad */) const;

  virtual void Project (Point<3> & p) const;

  virtual Point<3> GetSurfacePoint () const;
  virtual void Print (ostream & str) const;
  
  virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					 const Box<3> & boundingbox, 
					 double facets) const;

  const SplineGeometry<3> & GetPath(void) const {return *path;}
  const SplineSeg<2> & GetProfile(void) const {return *profile;}

  bool BoxIntersectsFace(const Box<3> & box) const;

  void LineIntersections ( const Point<3> & p,
			   const Vec<3> & v,
			   const double eps,
			   int & before,
			   int & after,
			   bool & intersecting ) const;

  INSOLID_TYPE VecInFace ( const Point<3> & p,
			   const Vec<3> & v,
			   const double eps ) const;

  const Vec<3> & GetYDir ( void ) const {return y_dir[latest_seg];}
  const Vec<3> & GetProfileTangent (void) const {return profile_tangent;}
  double GetProfilePar(void) const {return profile_par;}

  void GetRawData(ARRAY<double> & data) const;
};



class Extrusion : public Primitive
{
private:
  const SplineGeometry<3> & path;
  const SplineGeometry<2> & profile;

  const Vec<3> & z_direction;

  ARRAY<ExtrusionFace*> faces;

  mutable int latestfacenum;

public:
  Extrusion(const SplineGeometry<3> & path_in,
	    const SplineGeometry<2> & profile_in,
	    const Vec<3> & z_dir);
  ~Extrusion();
  virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
  virtual INSOLID_TYPE PointInSolid (const Point<3> & p,
				     double eps) const;
  INSOLID_TYPE PointInSolid (const Point<3> & p,
			     double eps,
			     ARRAY<int> * const facenums) const;
  virtual INSOLID_TYPE VecInSolid (const Point<3> & p,
				   const Vec<3> & v,
				   double eps) const;

  // checks if lim s->0 lim t->0  p + t(v1 + s v2) in solid
  virtual INSOLID_TYPE VecInSolid2 (const Point<3> & p,
				    const Vec<3> & v1,
				    const Vec<3> & v2,
				    double eps) const;

  
  virtual int GetNSurfaces() const;
  virtual Surface & GetSurface (int i = 0);
  virtual const Surface & GetSurface (int i = 0) const;


  virtual void Reduce (const BoxSphere<3> & box);
  virtual void UnReduce ();

};


#endif //_EXTRUSION_HPP
