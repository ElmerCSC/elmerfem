#ifndef FILE_SOLID
#define FILE_SOLID

/**************************************************************************/
/* File:   solid.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   1. Dez. 95                                                     */
/**************************************************************************/

/*

  Constructive Solid Model (csg)
  
*/




class Solid;

class SolidIterator
{
public:
  SolidIterator () { ; }
  virtual ~SolidIterator () { ; }
  virtual void Do (Solid * sol) = 0;
};



class Solid
{
public:
  
  typedef enum optyp1 { TERM, TERM_REF, SECTION, UNION, SUB, ROOT, DUMMY } optyp;
  
private:
  char * name;
  Primitive * prim;
  Solid * s1, * s2;
  
  optyp op;
  bool visited;
  double maxh;

  // static int cntnames;

public:
  Solid (Primitive * aprim);
  Solid (optyp aop, Solid * as1, Solid * as2 = NULL);
  ~Solid ();

  const char * Name () const { return name; }
  void SetName (const char * aname);

  Solid * Copy (class CSGeometry & geom) const;
  void Transform (Transformation<3> & trans);

  
  void IterateSolid (SolidIterator & it, bool only_once = 0);

  
  void Boundaries (const Point<3> & p, ARRAY<int> & bounds) const;
  int NumPrimitives () const;
  void GetSurfaceIndices (ARRAY<int> & surfind) const;
  void GetSurfaceIndices (IndexSet & iset) const;

  void GetTangentialSurfaceIndices (const Point<3> & p, ARRAY<int> & surfids, double eps) const;
  void GetTangentialSurfaceIndices2 (const Point<3> & p, const Vec<3> & v, ARRAY<int> & surfids, double eps) const;
  void GetTangentialSurfaceIndices3 (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2, ARRAY<int> & surfids, double eps) const;


  Primitive * GetPrimitive ()
    { return (op == TERM || op == TERM_REF) ? prim : NULL; }
  const Primitive * GetPrimitive () const
    { return (op == TERM || op == TERM_REF) ? prim : NULL; }

  Solid * S1() { return s1; }
  Solid * S2() { return s2; }

  // geometric tests

  bool IsIn (const Point<3> & p, double eps = 1e-6) const;
  bool IsStrictIn (const Point<3> & p, double eps = 1e-6) const;
  bool VectorIn (const Point<3> & p, const Vec<3> & v, double eps = 1e-6) const;
  bool VectorStrictIn (const Point<3> & p, const Vec<3> & v, double eps = 1e-6) const;
  
  bool VectorIn2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
		  double eps) const;
  bool VectorIn2Rec (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
		     double eps) const;


  /// compute localization in point p
  void TangentialSolid (const Point<3> & p, Solid *& tansol, ARRAY<int> & surfids, double eps) const;

  /// compute localization in point p tangential to vector t
  void TangentialSolid2 (const Point<3> & p, const Vec<3> & t,
			 Solid *& tansol, ARRAY<int> & surfids, double eps) const;

  /** compute localization in point p, with second order approximation to edge
      p + s t + s*s/2 t2 **/
  void TangentialSolid3 (const Point<3> & p, const Vec<3> & t, const Vec<3> & t2, 
			 Solid *& tansol, ARRAY<int> & surfids, double eps) const;



  /** tangential solid, which follows the edge
      p + s t + s*s/2 t2
      with second order, and the neighbouring face
      p + s t + s*s/2 t2 + r m
      with first order
  **/
  void TangentialEdgeSolid (const Point<3> & p, const Vec<3> & t, const Vec<3> & t2, 
			    const Vec<3> & m, 
			    Solid *& tansol, ARRAY<int> & surfids, double eps) const;


  void CalcOnePrimitiveSpecialPoints (const Box<3> & box, ARRAY<Point<3> > & pts) const;

  ///
  int Edge (const Point<3> & p, const Vec<3> & v, double eps) const;
  ///
  int OnFace (const Point<3> & p, const Vec<3> & v, double eps) const;
  ///
  void Print (ostream & str) const;
  ///
  void CalcSurfaceInverse ();
  ///
  Solid * GetReducedSolid (const BoxSphere<3> & box) const;
  

  void SetMaxH (double amaxh)
    { maxh = amaxh; }
  double GetMaxH () const
    { return maxh; }

  void GetSolidData (ostream & ost, int first = 1) const;
  static Solid * CreateSolid (istream & ist, const SYMBOLTABLE<Solid*> & solids);


  static BlockAllocator ball;
  void * operator new(size_t /* s */) 
  {
    return ball.Alloc();
  }

  void operator delete (void * p)
  {
    ball.Free (p);
  }


protected:
  ///

  void RecBoundaries (const Point<3> & p, ARRAY<int> & bounds, 
		      int & in, int & strin) const;
  ///
  void RecTangentialSolid (const Point<3> & p, Solid *& tansol, ARRAY<int> & surfids, 
			   int & in, int & strin, double eps) const;

  void RecTangentialSolid2 (const Point<3> & p, const Vec<3> & vec, 
			    Solid *& tansol, ARRAY<int> & surfids, 
			    int & in, int & strin, double eps) const;
  ///
  void RecTangentialSolid3 (const Point<3> & p, const Vec<3> & vec,const Vec<3> & vec2, 
			    Solid *& tansol, ARRAY<int> & surfids, 
			    int & in, int & strin, double eps) const;
  ///
  void RecTangentialEdgeSolid (const Point<3> & p, const Vec<3> & t, const Vec<3> & t2, 
			       const Vec<3> & m, 
			       Solid *& tansol, ARRAY<int> & surfids, 
			       int & in, int & strin, double eps) const;

  ///
  void RecEdge (const Point<3> & p, const Vec<3> & v,
                int & in, int & strin, int & faces, double eps) const;
  ///
  void CalcSurfaceInverseRec (int inv);
  ///
  Solid * RecGetReducedSolid (const BoxSphere<3> & box, INSOLID_TYPE & in) const;
  ///
  void RecGetSurfaceIndices (ARRAY<int> & surfind) const;
  void RecGetTangentialSurfaceIndices (const Point<3> & p, ARRAY<int> & surfids, double eps) const;
  void RecGetTangentialSurfaceIndices2 (const Point<3> & p, const Vec<3> & v, ARRAY<int> & surfids, double eps) const;
  void RecGetTangentialSurfaceIndices3 (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2, 
					ARRAY<int> & surfids, double eps) const;
  void RecGetTangentialEdgeSurfaceIndices (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2, const Vec<3> & m,
					   ARRAY<int> & surfids, double eps) const;
  void RecGetSurfaceIndices (IndexSet & iset) const;

  void RecCalcOnePrimitiveSpecialPoints (ARRAY<Point<3> > & pts) const;

  friend class SolidIterator;
  friend class ClearVisitedIt;
  friend class RemoveDummyIterator;
  friend class CSGeometry;
};


inline ostream & operator<< (ostream & ost, const Solid & sol)
{
  sol.Print (ost);
  return ost;
}






class ReducePrimitiveIterator : public SolidIterator
{
  const BoxSphere<3> & box;
public:
  ReducePrimitiveIterator (const BoxSphere<3> & abox)
    : SolidIterator(), box(abox) { ; }
  virtual ~ReducePrimitiveIterator () { ; }
  virtual void Do (Solid * sol)
  {
    if (sol -> GetPrimitive())
      sol -> GetPrimitive() -> Reduce (box);
  }
};


class UnReducePrimitiveIterator : public SolidIterator
{
public:
  UnReducePrimitiveIterator () { ; }
  virtual ~UnReducePrimitiveIterator () { ; }
  virtual void Do (Solid * sol)
  {
    if (sol -> GetPrimitive())
      sol -> GetPrimitive() -> UnReduce ();
  }
};


#endif
