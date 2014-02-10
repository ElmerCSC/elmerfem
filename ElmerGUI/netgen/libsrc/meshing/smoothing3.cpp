#include <mystdlib.h>

#include "meshing.hpp"
#ifdef SOLIDGEOM
#include <csg.hpp>
#endif
#include <opti.hpp>


namespace netgen
{
  

  double MinFunctionSum :: Func (const Vector & x) const
  {
    double retval = 0;
    for(int i=0; i<functions.Size(); i++)
      retval += functions[i]->Func(x);
      
    return retval;
  }

  void MinFunctionSum :: Grad (const Vector & x, Vector & g) const
  {
    g = 0.;
    static Vector gi(3);
    for(int i=0; i<functions.Size(); i++)
      {
	functions[i]->Grad(x,gi);
	for(int j=0; j<g.Size(); j++)
	  g[j] += gi[j];
      }
  }
      

  double MinFunctionSum :: FuncGrad (const Vector & x, Vector & g) const
  {
    double retval = 0;
    g = 0.;
    static Vector gi(3);
    for(int i=0; i<functions.Size(); i++)
      {
	retval += functions[i]->FuncGrad(x,gi);
	for(int j=0; j<g.Size(); j++)
	  g[j] += gi[j];
      }
    return retval;
  }

  double MinFunctionSum :: FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    double retval = 0;
    deriv = 0.;
    double derivi;
    for(int i=0; i<functions.Size(); i++)
      {
	retval += functions[i]->FuncDeriv(x,dir,derivi);
	deriv += derivi;
      }
    return retval;
  }

  double MinFunctionSum :: GradStopping (const Vector & x) const
  {
    double minfs(0), mini;
    for(int i=0; i<functions.Size(); i++)
      {
	mini = functions[i]->GradStopping(x);
	if(i==0 || mini < minfs)
	  minfs = mini;
      }
    return minfs;
  }


  void MinFunctionSum :: AddFunction(MinFunction & fun)
  {
    functions.Append(&fun);
  }
  
  const MinFunction & MinFunctionSum :: Function(int i) const
  {
    return *functions[i];
  }
  MinFunction & MinFunctionSum :: Function(int i)
  {
    return *functions[i];
  }


  PointFunction1 :: PointFunction1 (Mesh::T_POINTS & apoints, 
				    const ARRAY<INDEX_3> & afaces,
				    double ah)
    : points(apoints), faces(afaces)
  {
    h = ah;
  }
  

  double PointFunction1 :: Func (const Vector & vp) const
  {
    double badness = 0;
    Point<3> pp(vp(0), vp(1), vp(2));

    for (int j = 0; j < faces.Size(); j++)
      {
	const INDEX_3 & el = faces[j];

	double bad = CalcTetBadness (points[el.I1()], 
				     points[el.I3()], 
				     points[el.I2()], 
				     pp, 0);
	badness += bad;
      }
 
    return badness;
  }


  double PointFunction1 :: 
  FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    static Vector hx(3);
    static double eps = 1e-6;

    double dirlen = dir.L2Norm();
    if (dirlen < 1e-14)
      {
	deriv = 0;
	return Func(x);
      }

    hx.Set(1, x);
    hx.Add(eps * h / dirlen, dir);
    double fr = Func (hx);
    hx.Set(1, x);
    hx.Add(-eps * h / dirlen, dir);
    double fl = Func (hx);

    deriv = (fr - fl) / (2 * eps * h) * dirlen;

    return Func(x);
  }


  double PointFunction1 :: FuncGrad (const Vector & x, Vector & g) const
  {
    static Vector hx(3);
    static double eps = 1e-6;

    hx = x;
    for (int i = 1; i <= 3; i++)
      {
	hx.Elem(i) = x.Get(i) + eps * h;
	double fr = Func (hx);
	hx.Elem(i) = x.Get(i) - eps * h;
	double fl = Func (hx);
	hx.Elem(i) = x.Get(i);

	g.Elem(i) = (fr - fl) / (2 * eps * h);
      }

    return Func(x);
  }

  double PointFunction1 :: GradStopping (const Vector & x) const
  {
    double f = Func(x);
    return 1e-8 * f * f;
  }




  /* Cheap Functional depending of inner point inside triangular surface */

  // is it used ????
  class CheapPointFunction1 : public MinFunction
  {
    Mesh::T_POINTS & points;
    const ARRAY<INDEX_3> & faces;
    DenseMatrix m;
    double h;
  public:
    CheapPointFunction1 (Mesh::T_POINTS & apoints, 
			 const ARRAY<INDEX_3> & afaces,
			 double ah);
  
    virtual double Func (const Vector & x) const;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
  };

  CheapPointFunction1 :: CheapPointFunction1 (Mesh::T_POINTS & apoints, 
					      const ARRAY<INDEX_3> & afaces,
					      double ah)
    : points(apoints), faces(afaces)
  {
    h = ah;
  

    int nf = faces.Size();

    m.SetSize (nf, 4);
  
    for (int i = 1; i <= nf; i++)
      {
	const Point3d & p1 = points[faces.Get(i).I1()];
	const Point3d & p2 = points[faces.Get(i).I2()];
	const Point3d & p3 = points[faces.Get(i).I3()];
	Vec3d v1 (p1, p2);
	Vec3d v2 (p1, p3);
	Vec3d n;
	Cross (v1, v2, n);
	n /= n.Length();

	m.Elem(i, 1) = n.X();
	m.Elem(i, 2) = n.Y();
	m.Elem(i, 3) = n.Z();
	m.Elem(i, 4) = - (n.X() * p1.X() + n.Y() * p1.Y() + n.Z() * p1.Z());
      } 
  }
  
  double CheapPointFunction1 :: Func (const Vector & vp) const
  {

    /*
      int j;
      double badness = 0;
      Point3d pp(vp.Get(1), vp.Get(2), vp.Get(3));

      for (j = 1; j <= faces.Size(); j++)
      {
      const INDEX_3 & el = faces.Get(j);

      double bad = CalcTetBadness (points.Get(el.I1()), 
      points.Get(el.I3()), 
      points.Get(el.I2()), 
      pp, 0);
      badness += bad;
      }
    */

    int i;
    double badness = 0;
    static Vector hv(4);
    static Vector res;
    res.SetSize (m.Height());

    for (i = 1;i <= 3; i++)
      hv.Elem(i) = vp.Get(i);
    hv.Elem(4) = 1;
    m.Mult (hv, res);

    for (i = 1; i <= res.Size(); i++)
      {
	if (res.Get(i) < 1e-10)
	  badness += 1e24;
	else
	  badness += 1 / res.Get(i);
      }
 
    return badness;
  }


  double CheapPointFunction1 :: FuncGrad (const Vector & x, Vector & g) const
  {
    static Vector hx(3);
    static double eps = 1e-6;

    hx = x;
    for (int i = 1; i <= 3; i++)
      {
	hx.Elem(i) = x.Get(i) + eps * h;
	double fr = Func (hx);
	hx.Elem(i) = x.Get(i) - eps * h;
	double fl = Func (hx);
	hx.Elem(i) = x.Get(i);

	g.Elem(i) = (fr - fl) / (2 * eps * h);
      }

    return Func(x);
  }














  /* ************* PointFunction **************************** */


  class PointFunction 
  {
  public:
    Mesh::T_POINTS & points;
    const Mesh::T_VOLELEMENTS & elements;
    TABLE<INDEX,PointIndex::BASE> elementsonpoint;
    PointIndex actpind;
    double h;
  
  public:
    PointFunction (Mesh::T_POINTS & apoints, 
		   const Mesh::T_VOLELEMENTS & aelements);
  
    virtual void SetPointIndex (PointIndex aactpind);
    void SetLocalH (double ah) { h = ah; }
    double GetLocalH () const { return h; }
    virtual double PointFunctionValue (const Point<3> & pp) const;
    virtual double PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const;
    virtual double PointFunctionValueDeriv (const Point<3> & pp, const Vec<3> & dir, double & deriv) const;

    int MovePointToInner ();
  };


  PointFunction :: PointFunction (Mesh::T_POINTS & apoints, 
				  const Mesh::T_VOLELEMENTS & aelements)
    : points(apoints), elements(aelements), elementsonpoint(apoints.Size())
  {
    INDEX i;
    int j;
  
    for (i = 1; i <= elements.Size(); i++)
      {
	if (elements.Get(i).NP() == 4)
	  for (j = 1; j <= elements.Get(i).NP(); j++)
	    elementsonpoint.Add (elements.Get(i).PNum(j), i);  
      }
  }

  void PointFunction :: SetPointIndex (PointIndex aactpind)
  {
    actpind = aactpind; 
  }  

  double PointFunction :: PointFunctionValue (const Point<3> & pp) const
  {
    int j;
    INDEX eli;
    const Element * el;
    double badness;
    //  ARRAY<const Point3d*> p(4);
    Point<3> hp;

    badness = 0;

    hp = points[actpind];
    points[actpind] = Point<3> (pp);

    for (j = 0; j < elementsonpoint[actpind].Size(); j++)
      {
	eli = elementsonpoint[actpind][j];
	el = &elements.Get(eli);
	badness += CalcTetBadness (points[el->PNum(1)], 
				   points[el->PNum(2)], 
				   points[el->PNum(3)], 
				   points[el->PNum(4)], -1);
      }
  
    points[actpind] = Point<3> (hp); 
    return badness;
  }


  double PointFunction :: PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const
  {
    double f;//, delta = h * 1e-6;
    // Point<3> hpp;

    f = PointFunctionValue (pp);

    /*
      hpp = pp;
      hpp.X() = pp.X() + delta;
      fr = PointFunctionValue (hpp);
      hpp.X() = pp.X() - delta;
      fl = PointFunctionValue (hpp);
      grad.Elem(1) = (fr - fl) / (2 * delta);

      hpp = pp;
      hpp.Y() = pp.Y() + delta;
      fr = PointFunctionValue (hpp);
      hpp.Y() = pp.Y() - delta;
      fl = PointFunctionValue (hpp);
      grad.Elem(2) = (fr - fl) / (2 * delta);

      hpp = pp;
      hpp.Z() = pp.Z() + delta;
      fr = PointFunctionValue (hpp);
      hpp.Z() = pp.Z() - delta;
      fl = PointFunctionValue (hpp);
      grad.Elem(3) = (fr - fl) / (2 * delta);
    */


    // new gradient calculation
    //  double badness;
    Point<3> hp;
    Vec<3> vgradi, vgrad(0,0,0);

    //  badness = 0;

    hp = points[actpind];
    points[actpind] = Point<3> (pp);

    for (int j = 0; j < elementsonpoint[actpind].Size(); j++)
      {
	INDEX eli = elementsonpoint[actpind][j];
	const Element & el = elements.Get(eli);

	for (int k = 1; k <= 4; k++)
	  if (el.PNum(k) == actpind)
	    {
	      CalcTetBadnessGrad (points[el.PNum(1)], 
				  points[el.PNum(2)], 
				  points[el.PNum(3)], 
				  points[el.PNum(4)], -1, k, vgradi);

	      vgrad += vgradi;
	    }
      }

    points[actpind] = Point<3> (hp); 

    grad = vgrad;
    return f;
  }


  double PointFunction :: PointFunctionValueDeriv (const Point<3> & pp, const Vec<3> & dir,
						   double & deriv) const
  {
    double f;
    // Point<3> hpp;

    Vec<3> dirn (dir);
    //double ldir = dir.Length();

    int j, k;
    INDEX eli;
    //  double badness;
    Point<3> hp;
    Vec<3> vgradi, vgrad(0,0,0);

    //  badness = 0;

    hp = points[actpind];
    points[actpind] = pp;
    f = 0;

    for (j = 0; j < elementsonpoint[actpind].Size(); j++)
      {
	eli = elementsonpoint[actpind][j];
	const Element & el = elements.Get(eli);

	for (k = 1; k <= 4; k++)
	  if (el.PNum(k) == actpind)
	    {
	      f += CalcTetBadnessGrad (points[el.PNum(1)], 
				       points[el.PNum(2)], 
				       points[el.PNum(3)], 
				       points[el.PNum(4)], -1, k, vgradi);

	      vgrad += vgradi;
	    }
      }

    points[actpind] = Point<3> (hp); 
    deriv = dir * vgrad;
    return f;
  }

  int PointFunction :: MovePointToInner ()
  {
    // try point movement 
    ARRAY<Element2d> faces;
  
    for (int j = 0; j < elementsonpoint[actpind].Size(); j++)
      {
	const Element & el = 
	  elements.Get(elementsonpoint[actpind][j]);
      
	for (int k = 1; k <= 4; k++)
	  if (el.PNum(k) == actpind)
	    {
	      Element2d face;
	      el.GetFace (k, face);
	      Swap (face.PNum(2), face.PNum(3));
	      faces.Append (face);
	    }
      }
  
    Point3d hp;
    int hi = FindInnerPoint (points, faces, hp);
    if (hi)
      {
	// cout << "inner point found" << endl;
	points[actpind] = Point<3> (hp);
      }
    else
      ;
    //      cout << "no inner point found" << endl;

    /*
    Point3d hp2;
    int hi2 = FindInnerPoint (points, faces, hp2);
    if (hi2)
      {
	cout << "new: inner point found" << endl;
      }
    else
      cout << "new: no inner point found" << endl;
  
    (*testout) << "hi(orig) = " << hi << ", hi(new) = " << hi2;
    if (hi != hi2) (*testout) << "hi different" << endl;
    */

    return hi;
  }






  class CheapPointFunction : public PointFunction
  {
    DenseMatrix m;
  public:
    CheapPointFunction (Mesh::T_POINTS & apoints, 
			const Mesh::T_VOLELEMENTS & aelements);
    virtual void SetPointIndex (PointIndex aactpind);
    virtual double PointFunctionValue (const Point<3> & pp) const;
    virtual double PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const;
  };


  CheapPointFunction :: CheapPointFunction (Mesh::T_POINTS & apoints, 
					    const Mesh::T_VOLELEMENTS & aelements)
    : PointFunction (apoints, aelements)
  {
    ;
  }


  void CheapPointFunction :: SetPointIndex (PointIndex aactpind)
  {
    actpind = aactpind; 

    int ne = elementsonpoint[actpind].Size();
    int i, j;
    int pi1, pi2, pi3;

    m.SetSize (ne, 4);

    for (i = 0; i < ne; i++)
      {
	pi1 = 0;
	pi2 = 0;
	pi3 = 0;

	const Element & el = elements.Get (elementsonpoint[actpind][i]);
	for (j = 1; j <= 4; j++)
	  if (el.PNum(j) != actpind)
	    {
	      pi3 = pi2;
	      pi2 = pi1;
	      pi1 = el.PNum(j);
	    }

	const Point3d & p1 = points[pi1];
	Vec3d v1 (p1, points[pi2]);
	Vec3d v2 (p1, points[pi3]);
	Vec3d n;
	Cross (v1, v2, n);
	n /= n.Length();

	Vec3d v (p1, points[actpind]);
	double c = v * n;
      
	if (c < 0)
	  n *= -1;    
      
	// n is inner normal

	m.Elem(i+1, 1) = n.X();
	m.Elem(i+1, 2) = n.Y();
	m.Elem(i+1, 3) = n.Z();
	m.Elem(i+1, 4) = - (n.X() * p1.X() + n.Y() * p1.Y() + n.Z() * p1.Z());
      }
  }

  double CheapPointFunction :: PointFunctionValue (const Point<3> & pp) const
  {
    static Vector p4(4);
    static Vector di;
    int n = m.Height();

    p4.Elem(1) = pp(0);
    p4.Elem(2) = pp(1);
    p4.Elem(3) = pp(2);
    p4.Elem(4) = 1;

    di.SetSize (n);
    m.Mult (p4, di);
  
    double sum = 0;
    for (int i = 1; i <= n; i++)
      {
	if (di.Get(i) > 0)
	  sum += 1 / di.Get(i);
	else
	  return 1e16;
      }
    return sum;
  }




  double CheapPointFunction :: PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const
  {
    static Vector p4(4);
    static Vector di;

    int n = m.Height();

    p4.Elem(1) = pp(0);
    p4.Elem(2) = pp(1);
    p4.Elem(3) = pp(2);
    p4.Elem(4) = 1;

    di.SetSize (n);
    m.Mult (p4, di);
  
    double sum = 0;
    grad = 0;
    for (int i = 1; i <= n; i++)
      {
	if (di.Get(i) > 0)
	  {
	    double idi = 1 / di.Get(i);
	    sum += idi;
	    grad(0) -= idi * idi * m.Get(i, 1);
	    grad(1) -= idi * idi * m.Get(i, 2);
	    grad(2) -= idi * idi * m.Get(i, 3);
	  }
	else
	  {
	    return 1e16;
	  }
      }
    return sum;
  }








  class Opti3FreeMinFunction : public MinFunction
  { 
    const PointFunction & pf;
    Point<3> sp1;
  
  public:
    Opti3FreeMinFunction (const PointFunction & apf);
    void SetPoint (const Point<3> & asp1) { sp1 = asp1; }
    virtual double Func (const Vector & x) const;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;  
    virtual double GradStopping (const Vector & x) const;
    virtual void ApproximateHesse (const Vector & x,
				   DenseMatrix & hesse) const;
  };

  Opti3FreeMinFunction :: Opti3FreeMinFunction (const PointFunction & apf)
    : pf(apf)
  {
    ;
  }

  double Opti3FreeMinFunction :: Func (const Vector & x) const
  {
    Point<3> pp;
    for (int j = 0; j < 3; j++)
      pp(j) = sp1(j) + x(j);
    return pf.PointFunctionValue (pp);
  }
  
  double Opti3FreeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
  {
    Vec<3> vgrad;
    Point<3> pp;

    for (int j = 0; j < 3; j++)
      pp(j) = sp1(j) + x(j);

    double val = pf.PointFunctionValueGrad (pp, vgrad);

    for (int j = 0; j < 3; j++)
      grad(j) = vgrad(j);

    return val;
  }

  double Opti3FreeMinFunction :: FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    Point<3> pp;

    for (int j = 0; j < 3; j++)
      pp(j) = sp1(j) + x(j);

    Vec<3> vdir;
    for (int j = 0; j < 3; j++)
      vdir(j) = dir(j);

    return pf.PointFunctionValueDeriv (pp, vdir, deriv);
  }
  
  double Opti3FreeMinFunction :: GradStopping (const Vector & x) const
  {
    double f = Func(x);
    return 1e-3 * f / pf.GetLocalH();
  }


  void Opti3FreeMinFunction :: ApproximateHesse (const Vector & x,
						 DenseMatrix & hesse) const
  {
    int n = x.Size();

    static Vector hx;
    hx.SetSize(n);

    double eps = 1e-8;
    double f, f11, f22; //, f12, f21

    f = Func(x);
  
    for (int i = 1; i <= n; i++)
      {
	for (int j = 1; j < i; j++)
	  {
	    /*
	      hx = x;
	      hx.Elem(i) = x.Get(i) + eps;
	      hx.Elem(j) = x.Get(j) + eps;
	      f11 = Func(hx);
	      hx.Elem(i) = x.Get(i) + eps;
	      hx.Elem(j) = x.Get(j) - eps;
	      f12 = Func(hx);
	      hx.Elem(i) = x.Get(i) - eps;
	      hx.Elem(j) = x.Get(j) + eps;
	      f21 = Func(hx);
	      hx.Elem(i) = x.Get(i) - eps;
	      hx.Elem(j) = x.Get(j) - eps;
	      f22 = Func(hx);
	    */
	    hesse.Elem(i, j) = hesse.Elem(j, i) = 0;
	    //	    (f11 + f22 - f12 - f21) / (2 * eps * eps);
	  }

	hx = x;
	hx.Elem(i) = x.Get(i) + eps;
	f11 = Func(hx);
	hx.Elem(i) = x.Get(i) - eps;
	f22 = Func(hx);

	hesse.Elem(i, i) = (f11 + f22 - 2 * f) / (eps * eps) + 1e-12;
      }
  }






#ifdef SOLIDGEOM
  class Opti3SurfaceMinFunction : public MinFunction
  {
    const PointFunction & pf;
    Point3d sp1;
    const Surface * surf;
    Vec3d t1, t2;
  
  public:
    Opti3SurfaceMinFunction (const PointFunction & apf);
  
    void SetPoint (const Surface * asurf, const Point3d & asp1);

    void CalcNewPoint (const Vector & x, Point3d & np) const; 
    virtual double Func (const Vector & x) const;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
  };


  Opti3SurfaceMinFunction :: Opti3SurfaceMinFunction (const PointFunction & apf)
    : MinFunction(), pf(apf)
  {
    ;
  }

  void Opti3SurfaceMinFunction :: SetPoint (const Surface * asurf, const Point3d & asp1)
  { 
    Vec3d n;
    sp1 = asp1; 
    surf = asurf;
  
    Vec<3> hn;
    surf -> GetNormalVector (sp1, hn);
    n = hn;

    n.GetNormal (t1);
    t1 /= t1.Length();
    t2 = Cross (n, t1);
  }

  
  void Opti3SurfaceMinFunction :: CalcNewPoint (const Vector & x, 
						Point3d & np) const
  {
    np.X() = sp1.X() + x.Get(1) * t1.X() + x.Get(2) * t2.X();
    np.Y() = sp1.Y() + x.Get(1) * t1.Y() + x.Get(2) * t2.Y();
    np.Z() = sp1.Z() + x.Get(1) * t1.Z() + x.Get(2) * t2.Z();

    Point<3> hnp = np;
    surf -> Project (hnp);
    np = hnp;
  }


  double Opti3SurfaceMinFunction :: Func (const Vector & x) const
  {
    Point3d pp1;

    CalcNewPoint (x, pp1);
    return pf.PointFunctionValue (pp1);
  }



  double Opti3SurfaceMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
  {
    Vec3d n, vgrad;
    Point3d pp1;
    double badness;
    static Vector freegrad(3);

    CalcNewPoint (x, pp1);

    badness = pf.PointFunctionValueGrad (pp1, freegrad);
    vgrad.X() = freegrad.Get(1);
    vgrad.Y() = freegrad.Get(2);
    vgrad.Z() = freegrad.Get(3);

    Vec<3> hn;
    surf -> GetNormalVector (pp1, hn);
    n = hn;

    vgrad -= (vgrad * n) * n;

    grad.Elem(1) = vgrad * t1;
    grad.Elem(2) = vgrad * t2;
    
    return badness;
  }
#endif
  
  
  


  
  
  
#ifdef SOLIDGEOM
  class Opti3EdgeMinFunction : public MinFunction
  {
    const PointFunction & pf;
    Point3d sp1;
    const Surface *surf1, *surf2;
    Vec3d t1;
  
  public:
    Opti3EdgeMinFunction (const PointFunction & apf);
  
    void SetPoint (const Surface * asurf1, const Surface * asurf2,
		   const Point3d & asp1);
    void CalcNewPoint (const Vector & x, Point3d & np) const; 
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double Func (const Vector & x) const;
  };

  Opti3EdgeMinFunction :: Opti3EdgeMinFunction (const PointFunction & apf)
    : MinFunction(), pf(apf)
  {
    ;
  }
  
  void Opti3EdgeMinFunction :: SetPoint (const Surface * asurf1, 
					 const Surface * asurf2, 
					 const Point3d & asp1) 
  { 
    Vec3d n1, n2;
    sp1 = asp1; 
    surf1 = asurf1;
    surf2 = asurf2;

    Vec<3> hn1, hn2;
    surf1 -> GetNormalVector (sp1, hn1);
    surf2 -> GetNormalVector (sp1, hn2);
    n1 = hn1;
    n2 = hn2;
    t1 = Cross (n1, n2);
  }

  void Opti3EdgeMinFunction :: CalcNewPoint (const Vector & x,
					     Point3d & np) const
{
  np.X() = sp1.X() + x.Get(1) * t1.X();
  np.Y() = sp1.Y() + x.Get(1) * t1.Y();
  np.Z() = sp1.Z() + x.Get(1) * t1.Z();
  Point<3> hnp = np;
  ProjectToEdge (surf1, surf2, hnp);
  np = hnp;
}   

double Opti3EdgeMinFunction :: Func (const Vector & x) const
{
  Vector g(x.Size());
  return FuncGrad (x, g);
}


double Opti3EdgeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
{
  Vec3d n1, n2, v1, vgrad;
  Point3d pp1;
  double badness;
  static Vector freegrad(3);

  CalcNewPoint (x, pp1);


  badness = pf.PointFunctionValueGrad (pp1, freegrad);

  vgrad.X() = freegrad.Get(1);
  vgrad.Y() = freegrad.Get(2);
  vgrad.Z() = freegrad.Get(3);

  Vec<3> hn1, hn2;
  surf1 -> GetNormalVector (pp1, hn1);
  surf2 -> GetNormalVector (pp1, hn2);
  n1 = hn1;
  n2 = hn2;

  v1 = Cross (n1, n2);
  v1 /= v1.Length();

  grad.Elem(1) = (vgrad * v1) * (t1 * v1);
  return badness;
}
#endif




double CalcBad (const Mesh::T_POINTS & points, const Element & elem,
		double h)
{
  if (elem.GetType() == TET)
    return CalcTetBadness (points[elem.PNum(1)], 
			   points[elem.PNum(2)],  
			   points[elem.PNum(3)],  
			   points[elem.PNum(4)], h);  
  return 0;
}


extern double teterrpow;
double CalcTotalBad (const Mesh::T_POINTS & points, 
		     const Mesh::T_VOLELEMENTS & elements)
{
  int i;
  double sum = 0;
  double elbad;
  
  tets_in_qualclass.SetSize(20);
  for (i = 1; i <= 20; i++)
    tets_in_qualclass.Elem(i) = 0;


  for (i = 1; i <= elements.Size(); i++)
    {
      elbad = pow (max2(CalcBad (points, elements.Get(i), 0),1e-10),
		   1/teterrpow);

      int qualclass = int (20 / elbad + 1);
      if (qualclass < 1) qualclass = 1;
      if (qualclass > 20) qualclass = 20;
      tets_in_qualclass.Elem(qualclass)++;

      sum += elbad;
    }
  return sum;
}

int WrongOrientation (const Mesh::T_POINTS & points, const Element & el)
{
  const Point3d & p1 = points[el.PNum(1)];
  const Point3d & p2 = points[el.PNum(2)];
  const Point3d & p3 = points[el.PNum(3)];
  const Point3d & p4 = points[el.PNum(4)];

  Vec3d v1(p1, p2);
  Vec3d v2(p1, p3);
  Vec3d v3(p1, p4);
  Vec3d n;

  Cross (v1, v2, n);
  double vol = n * v3;

  return (vol > 0);
}











/* ************* JacobianPointFunction **************************** */




// class JacobianPointFunction : public MinFunction
// {
// public:
//   Mesh::T_POINTS & points;
//   const Mesh::T_VOLELEMENTS & elements;
//   TABLE<INDEX> elementsonpoint;
//   PointIndex actpind;
  
// public:
//   JacobianPointFunction (Mesh::T_POINTS & apoints, 
// 			 const Mesh::T_VOLELEMENTS & aelements);
  
//   virtual void SetPointIndex (PointIndex aactpind);
//   virtual double Func (const Vector & x) const;
//   virtual double FuncGrad (const Vector & x, Vector & g) const;
//   virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
// };


JacobianPointFunction :: 
JacobianPointFunction (Mesh::T_POINTS & apoints, 
		       const Mesh::T_VOLELEMENTS & aelements)
  : points(apoints), elements(aelements), elementsonpoint(apoints.Size())
{
  INDEX i;
  int j;
  
  for (i = 1; i <= elements.Size(); i++)
    {
      for (j = 1; j <= elements.Get(i).NP(); j++)
	elementsonpoint.Add1 (elements.Get(i).PNum(j), i);  
    }

  onplane = false;
}

void JacobianPointFunction :: SetPointIndex (PointIndex aactpind)
{
  actpind = aactpind; 
}  


double JacobianPointFunction :: Func (const Vector & v) const
{
  int j;
  double badness = 0;

  Point<3> hp = points.Elem(actpind);

  points.Elem(actpind) = hp + Vec<3> (v.Get(1), v.Get(2), v.Get(3));

  if(onplane)
    points.Elem(actpind) -= (v.Get(1)*nv(0)+v.Get(2)*nv(1)+v.Get(3)*nv(2)) * nv;


  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
    {
      int eli = elementsonpoint.Get(actpind, j);
      badness += elements.Get(eli).CalcJacobianBadness (points);
    }
  
  points.Elem(actpind) = hp; 

  return badness;
}





double JacobianPointFunction :: 
FuncGrad (const Vector & x, Vector & g) const
{
  int j, k;
  int lpi;
  double badness = 0;//, hbad;

  Point<3> hp = points.Elem(actpind);
  points.Elem(actpind) = hp + Vec<3> (x.Get(1), x.Get(2), x.Get(3));

  if(onplane)
    points.Elem(actpind) -= (x.Get(1)*nv(0)+x.Get(2)*nv(1)+x.Get(3)*nv(2)) * nv;

  Vec<3> hderiv;
  //Vec3d vdir;
  g.SetSize(3);
  g = 0;

  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
    {
      int eli = elementsonpoint.Get(actpind, j);
      const Element & el = elements.Get(eli);

      lpi = 0;
      for (k = 1; k <= el.GetNP(); k++)
	if (el.PNum(k) == actpind)
	  lpi = k;
      if (!lpi) cerr << "loc point not found" << endl;

      badness += elements.Get(eli).
	CalcJacobianBadnessGradient (points, lpi, hderiv);

      for(k=0; k<3; k++)
	g.Elem(k+1) += hderiv(k);
	
      /*
      for (k = 1; k <= 3; k++)
	{
	  vdir = Vec3d(0,0,0);
	  vdir.X(k) = 1;

	  hbad = elements.Get(eli).
	    CalcJacobianBadnessDirDeriv (points, lpi, vdir, hderiv);
	  //(*testout) << "hderiv " << k << ": " << hderiv << endl;
	  g.Elem(k) += hderiv;
	  if (k == 1)
	    badness += hbad;
	}
      */
    }

  if(onplane)
    {
      double scal = nv(0)*g.Get(1) + nv(1)*g.Get(2) + nv(2)*g.Get(3);
      g.Elem(1) -= scal*nv(0);
      g.Elem(2) -= scal*nv(1);
      g.Elem(3) -= scal*nv(2);
    }

  //(*testout) << "g = " << g << endl;

  
  points.Elem(actpind) = hp; 

  return badness;
}


double JacobianPointFunction :: 
FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
{
  int j, k;
  int lpi;
  double badness = 0;

  Point<3> hp = points.Elem(actpind);
  points.Elem(actpind) = Point<3> (hp + Vec3d (x.Get(1), x.Get(2), x.Get(3)));

  if(onplane)
    points.Elem(actpind) -= (Vec3d (x.Get(1), x.Get(2), x.Get(3))*nv) * nv;

  double hderiv;
  deriv = 0;
  Vec<3> vdir(dir.Get(1), dir.Get(2), dir.Get(3));
 
  if(onplane)
    {
      double scal = vdir * nv;
      vdir -= scal*nv;
    }

  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
    {
      int eli = elementsonpoint.Get(actpind, j);
      const Element & el = elements.Get(eli);

      lpi = 0;
      for (k = 1; k <= el.GetNP(); k++)
	if (el.PNum(k) == actpind)
	  lpi = k;
      if (!lpi) cerr << "loc point not found" << endl;

      badness += elements.Get(eli).
	CalcJacobianBadnessDirDeriv (points, lpi, vdir, hderiv);
      deriv += hderiv;
    }
  
  points.Elem(actpind) = hp; 

  return badness;
  
}










#ifdef SOLIDGEOMxxxx
void Mesh :: ImproveMesh (const CSGeometry & geometry, OPTIMIZEGOAL goal)
{
  INDEX i, eli;
  int j;
  int typ = 1;

  if (!&geometry || geometry.GetNSurf() == 0)
    {
      ImproveMesh (goal);
      return;
    }

  const char * savetask = multithread.task;
  multithread.task = "Smooth Mesh";


  TABLE<INDEX> surfelementsonpoint(points.Size());
  Vector x(3), xsurf(2), xedge(1);
  int surf, surf1, surf2, surf3;

  int uselocalh = mparam.uselocalh;

  (*testout).precision(8);
  (*testout) << "Improve Mesh" << "\n";
  PrintMessage (3, "ImproveMesh");
  //  (*mycout) << "Vol = " << CalcVolume (points, volelements) << endl;


  for (i = 1; i <= surfelements.Size(); i++)
    for (j = 1; j <= 3; j++)
      surfelementsonpoint.Add1 (surfelements.Get(i).PNum(j), i);


  PointFunction * pf;
  if (typ == 1)
    pf = new PointFunction(points, volelements);
  else
    pf = new CheapPointFunction(points, volelements);

  //  pf->SetLocalH (h);
  
  Opti3FreeMinFunction freeminf(*pf);
  Opti3SurfaceMinFunction surfminf(*pf);
  Opti3EdgeMinFunction edgeminf(*pf);
  
  OptiParameters par;
  par.maxit_linsearch = 20;
  par.maxit_bfgs = 20;

  int printmod = 1;
  char printdot = '.';
  if (points.Size() > 1000)
    {
      printmod = 10;
      printdot = '+';
    }
  if (points.Size() > 10000)
    {
      printmod = 100;
      printdot = '*';
    }

  for (i = 1; i <= points.Size(); i++)
    {
      //      if (ptyps.Get(i) == FIXEDPOINT) continue;
      if (ptyps.Get(i) != INNERPOINT) continue;

      if (multithread.terminate)
	throw NgException ("Meshing stopped");
      /*
      if (multithread.terminate)
	break;
      */
      multithread.percent = 100.0 * i /points.Size();

      /*
      if (points.Size() < 1000)
	PrintDot ();
      else
	if (i % 10 == 0)
	  PrintDot ('+');
      */
      if (i % printmod == 0) PrintDot (printdot);

      //    (*testout) << "Now point " << i << "\n";
      //    (*testout) << "Old: " << points.Get(i) << "\n";

      pf->SetPointIndex (i);

      //      if (uselocalh)
      {
	double lh = GetH (points.Get(i));
	pf->SetLocalH (GetH (points.Get(i)));
	par.typx = lh / 10;
	//	  (*testout) << "lh(" << points.Get(i) << ") = " << lh << "\n";
      }

      surf1 = surf2 = surf3 = 0;

      for (j = 1; j <= surfelementsonpoint.EntrySize(i); j++)
	{
	  eli = surfelementsonpoint.Get(i, j);
	  int surfi = surfelements.Get(eli).GetIndex();

	  if (surfi)
	    {
	      surf = GetFaceDescriptor(surfi).SurfNr();
	    
	      if (!surf1)
		surf1 = surf;
	      else if (surf1 != surf)
		{
		  if (!surf2)
		    surf2 = surf;
		  else if (surf2 != surf)
		    surf3 = surf;
		}
	    }
	  else
	    {
	      surf1 = surf2 = surf3 = 1;   // simulates corner point
	    }
	}


      if (surf2 && !surf3)
	{
	  //      (*testout) << "On Edge" << "\n";
	  /*
	    xedge = 0;
	    edgeminf.SetPoint (geometry.GetSurface(surf1),
	    geometry.GetSurface(surf2), 
	    points.Elem(i));
	    BFGS (xedge, edgeminf, par);

	    edgeminf.CalcNewPoint (xedge, points.Elem(i));
	  */
	}

      if (surf1 && !surf2)
	{
	  //      (*testout) << "In Surface" << "\n";
	  /*
	    xsurf = 0;
	    surfminf.SetPoint (geometry.GetSurface(surf1),
	    points.Get(i));
	    BFGS (xsurf, surfminf, par);
   
	    surfminf.CalcNewPoint (xsurf, points.Elem(i));
	  */
	}
 
      if (!surf1)
	{
	  //      (*testout) << "In Volume" << "\n";
	  x = 0;
	  freeminf.SetPoint (points.Elem(i));
	  //	  par.typx = 
	  BFGS (x, freeminf, par);

	  points.Elem(i).X() += x.Get(1);
	  points.Elem(i).Y() += x.Get(2);
	  points.Elem(i).Z() += x.Get(3);
	}
      
      //    (*testout) << "New Point: " << points.Elem(i) << "\n" << "\n";
    
    }
  PrintDot ('\n');
  //  (*mycout) << "Vol = " << CalcVolume (points, volelements) << endl;

  multithread.task = savetask;

}
#endif



  
void Mesh :: ImproveMesh (OPTIMIZEGOAL goal)
{
  int typ = 1;
  
  (*testout) << "Improve Mesh" << "\n";
  PrintMessage (3, "ImproveMesh");

  int np = GetNP();
  int ne = GetNE();


  ARRAY<double,PointIndex::BASE> perrs(np);
  perrs = 1.0;

  double bad1 = 0;
  double badmax = 0;

  if (goal == OPT_QUALITY)
    {
      for (int i = 1; i <= ne; i++)
	{
	  const Element & el = VolumeElement(i);
	  if (el.GetType() != TET)
	    continue;
	  
	  double hbad = CalcBad (points, el, 0);
	  for (int j = 0; j < 4; j++)
	    perrs[el[j]] += hbad;
	  
	  bad1 += hbad;
	}
      
      for (PointIndex i = PointIndex::BASE; i < np+PointIndex::BASE; i++)
	if (perrs[i] > badmax) 
	  badmax = perrs[i];
      badmax = 0;
    }

  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (points, volelements);
      (*testout) << "Total badness = " << bad1 << endl;
      PrintMessage (5, "Total badness = ", bad1);
    }
  
  Vector x(3);
  
  (*testout).precision(8);
  
  //int uselocalh = mparam.uselocalh;


  PointFunction * pf;

  if (typ == 1)
    pf = new PointFunction(points, volelements);
  else
    pf = new CheapPointFunction(points, volelements);

  //  pf->SetLocalH (h);
  
  Opti3FreeMinFunction freeminf(*pf);

  OptiParameters par;
  par.maxit_linsearch = 20;
  par.maxit_bfgs = 20;

  ARRAY<double, PointIndex::BASE> pointh (points.Size());

  if(lochfunc)
    {
      for(int i=1; i<=points.Size(); i++)
	pointh[i] = GetH(points.Get(i));
    }
  else
    {
      pointh = 0;
      for(int i=0; i<GetNE(); i++)
	{
	  const Element & el = VolumeElement(i+1);
	  double h = pow(el.Volume(points),1./3.);
	  for(int j=1; j<=el.GetNV(); j++)
	    if(h > pointh[el.PNum(j)])
	      pointh[el.PNum(j)] = h;
	}
    }
 

  int printmod = 1;
  char printdot = '.';
  if (points.Size() > 1000)
    {
      printmod = 10;
      printdot = '+';
    }
  if (points.Size() > 10000)
    {
      printmod = 100;
      printdot = '*';
    }


  const char * savetask = multithread.task;
  multithread.task = "Smooth Mesh";
  
  for (PointIndex i = PointIndex::BASE; 
       i < points.Size()+PointIndex::BASE; i++)
    if ( (*this)[i].Type() == INNERPOINT && perrs[i] > 0.01 * badmax)
      {
	if (multithread.terminate)
	  throw NgException ("Meshing stopped");

	multithread.percent = 100.0 * (i+1-PointIndex::BASE) / points.Size();
        /*
	if (points.Size() < 1000)
	  PrintDot ();
	else
	  if ( (i+1-PointIndex::BASE) % 10 == 0)
	    PrintDot ('+');
        */
        if (  (i+1-PointIndex::BASE) % printmod == 0) PrintDot (printdot);

	double lh = pointh[i];
	pf->SetLocalH (lh);
	par.typx = lh;

	freeminf.SetPoint (points[i]);
	pf->SetPointIndex (i);

	x = 0;
	int pok;
	pok = freeminf.Func (x) < 1e10; 

	if (!pok)
	  {
	    pok = pf->MovePointToInner ();

	    freeminf.SetPoint (points[i]);
	    pf->SetPointIndex (i);
	  }

	if (pok)
	  {
            //*testout << "start BFGS, pok" << endl;
	    BFGS (x, freeminf, par);
            //*testout << "BFGS complete, pok" << endl;
	    points[i](0) += x.Get(1);
	    points[i](1) += x.Get(2);
	    points[i](2) += x.Get(3);
	  }
      }
  PrintDot ('\n');
  
  
  delete pf;

  multithread.task = savetask;

  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (points, volelements);
      (*testout) << "Total badness = " << bad1 << endl;
      PrintMessage (5, "Total badness = ", bad1);
    }
}




// Improve Condition number of Jacobian, any elements  
void Mesh :: ImproveMeshJacobian (OPTIMIZEGOAL goal, const BitArray * usepoint)
{
  int i, j;
  
  (*testout) << "Improve Mesh Jacobian" << "\n";
  PrintMessage (3, "ImproveMesh Jacobian");

  int np = GetNP();
  int ne = GetNE();

  
  Vector x(3);
  
  (*testout).precision(8);
  
  JacobianPointFunction pf(points, volelements);
  

  OptiParameters par;
  par.maxit_linsearch = 20;
  par.maxit_bfgs = 20;
  
  BitArray badnodes(np);
  badnodes.Clear();

  for (i = 1; i <= ne; i++)
    {
      const Element & el = VolumeElement(i);
      double bad = el.CalcJacobianBadness (Points());
      if (bad > 1)
	for (j = 1; j <= el.GetNP(); j++)
	  badnodes.Set (el.PNum(j));
    }

  ARRAY<double, PointIndex::BASE> pointh (points.Size());

  if(lochfunc)
    {
      for(i = 1; i<=points.Size(); i++)
	pointh[i] = GetH(points.Get(i));
    }
  else
    {
      pointh = 0;
      for(i=0; i<GetNE(); i++)
	{
	  const Element & el = VolumeElement(i+1);
	  double h = pow(el.Volume(points),1./3.);
	  for(j=1; j<=el.GetNV(); j++)
	    if(h > pointh[el.PNum(j)])
	      pointh[el.PNum(j)] = h;
	}
    }
 


  const char * savetask = multithread.task;
  multithread.task = "Smooth Mesh Jacobian";
  
  for (i = 1; i <= points.Size(); i++)
    {
      if ((*this)[PointIndex(i)].Type() != INNERPOINT)
	continue;

      if(usepoint && !usepoint->Test(i))
	continue;

      //(*testout) << "improvejac, p = " << i << endl;

      if (goal == OPT_WORSTCASE && !badnodes.Test(i))
	continue;
      //	(*testout) << "smoot p " << i << endl;

      /*
	if (multithread.terminate)
	break;
      */
      if (multithread.terminate)
	throw NgException ("Meshing stopped");

      multithread.percent = 100.0 * i / points.Size();

      if (points.Size() < 1000)
	PrintDot ();
      else
	if (i % 10 == 0)
	  PrintDot ('+');

      double lh = pointh[i];
      par.typx = lh;

      pf.SetPointIndex (i);

      x = 0;
      int pok = (pf.Func (x) < 1e10); 

      if (pok)
	{
          //*testout << "start BFGS, Jacobian" << endl;
	  BFGS (x, pf, par);
          //*testout << "end BFGS, Jacobian" << endl;
	  points.Elem(i)(0) += x.Get(1);
	  points.Elem(i)(1) += x.Get(2);
	  points.Elem(i)(2) += x.Get(3);
	}
      else
	{
	  cout << "el not ok" << endl;
	}
    }
  PrintDot ('\n');
  

  multithread.task = savetask;
}




// Improve Condition number of Jacobian, any elements  
void Mesh :: ImproveMeshJacobianOnSurface (const BitArray & usepoint, 
					   const ARRAY< Vec<3>* > & nv,
					   OPTIMIZEGOAL goal,
					   const ARRAY< ARRAY<int,PointIndex::BASE>* > * idmaps)
{
  int i, j;
  
  (*testout) << "Improve Mesh Jacobian" << "\n";
  PrintMessage (3, "ImproveMesh Jacobian");

  int np = GetNP();
  int ne = GetNE();

  
  Vector x(3);
  
  (*testout).precision(8);
  
  JacobianPointFunction pf(points, volelements);

  ARRAY< ARRAY<int,PointIndex::BASE>* > locidmaps;
  const ARRAY< ARRAY<int,PointIndex::BASE>* > * used_idmaps;

  if(idmaps)
    used_idmaps = idmaps;
  else
    {
      used_idmaps = &locidmaps;
      
      for(i=1; i<=GetIdentifications().GetMaxNr(); i++)
	{
	  if(GetIdentifications().GetType(i) == Identifications::PERIODIC)
	    {
	      locidmaps.Append(new ARRAY<int,PointIndex::BASE>);
	      GetIdentifications().GetMap(i,*locidmaps.Last(),true);
	    }
	}
    }

  
  bool usesum = (used_idmaps->Size() > 0);
  MinFunctionSum pf_sum;
  
  JacobianPointFunction * pf2ptr = NULL;
  if(usesum)
    {
      pf2ptr = new JacobianPointFunction(points, volelements);
      pf_sum.AddFunction(pf);
      pf_sum.AddFunction(*pf2ptr);
    }
  

  OptiParameters par;
  par.maxit_linsearch = 20;
  par.maxit_bfgs = 20;
  
  BitArray badnodes(np);
  badnodes.Clear();

  for (i = 1; i <= ne; i++)
    {
      const Element & el = VolumeElement(i);
      double bad = el.CalcJacobianBadness (Points());
      if (bad > 1)
	for (j = 1; j <= el.GetNP(); j++)
	  badnodes.Set (el.PNum(j));
    }

  ARRAY<double, PointIndex::BASE> pointh (points.Size());
 
  if(lochfunc)
    {
      for(i=1; i<=points.Size(); i++)
	pointh[i] = GetH(points.Get(i));
    }
  else
    {
      pointh = 0;
      for(i=0; i<GetNE(); i++)
	{
	  const Element & el = VolumeElement(i+1);
	  double h = pow(el.Volume(points),1./3.);
	  for(j=1; j<=el.GetNV(); j++)
	    if(h > pointh[el.PNum(j)])
	      pointh[el.PNum(j)] = h;
	}
    }


  const char * savetask = multithread.task;
  multithread.task = "Smooth Mesh Jacobian";
  
  for (i = 1; i <= points.Size(); i++)
    if ( usepoint.Test(i) )
      {
	//(*testout) << "improvejac, p = " << i << endl;

	if (goal == OPT_WORSTCASE && !badnodes.Test(i))
	  continue;
	//	(*testout) << "smoot p " << i << endl;

	/*
	if (multithread.terminate)
	  break;
	*/
	if (multithread.terminate)
	  throw NgException ("Meshing stopped");

	multithread.percent = 100.0 * i / points.Size();

	if (points.Size() < 1000)
	  PrintDot ();
	else
	  if (i % 10 == 0)
	    PrintDot ('+');

	double lh = pointh[i];//GetH(points.Get(i));
	par.typx = lh;

	pf.SetPointIndex (i);

	int brother = -1;
	if(usesum)
	  {
	    for(j=0; brother == -1 && j<used_idmaps->Size(); j++)
	      {
		if(i < (*used_idmaps)[j]->Size() + PointIndex::BASE)
		  {
		    brother = (*(*used_idmaps)[j])[i];
		    if(brother == i || brother == 0)
		      brother = -1;
		  }
	      }
	    if(brother >= i)
	      {
		pf2ptr->SetPointIndex(brother);
		pf2ptr->SetNV(*nv[brother-1]);
	      }
	  }

	if(usesum && brother < i)
	  continue;

	//pf.UnSetNV(); x = 0;
	//(*testout) << "before " << pf.Func(x);

	pf.SetNV(*nv[i-1]);

	x = 0;
	int pok = (brother == -1) ? (pf.Func (x) < 1e10) : (pf_sum.Func (x) < 1e10);

	if (pok)
	  {
	    
	    if(brother == -1)
	      BFGS (x, pf, par);
	    else
	      BFGS (x, pf_sum, par);


	    for(j=1; j<=3; j++)
	      points.Elem(i)(j-1) += x.Get(j);// - scal*nv[i-1].X(j);

	    if(brother != -1)
	      for(j=1; j<=3; j++)
		points.Elem(brother)(j-1) += x.Get(j);// - scal*nv[brother-1].X(j);


	  }
	else
	  {
	    cout << "el not ok" << endl;
	    (*testout) << "el not ok" << endl
		       << "   func " << ((brother == -1) ? pf.Func(x) : pf_sum.Func (x)) << endl;
	    if(brother != -1)
	      (*testout) << "   func1 " << pf.Func(x) << endl
			 << "   func2 " << pf2ptr->Func(x) << endl;
	  }
      }
  
  PrintDot ('\n');

  delete pf2ptr;
  for(i=0; i<locidmaps.Size(); i++)
    delete locidmaps[i];

  multithread.task = savetask;
}




}
