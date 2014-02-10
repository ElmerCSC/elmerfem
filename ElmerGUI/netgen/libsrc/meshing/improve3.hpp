#ifndef FILE_IMPROVE3
#define FILE_IMPROVE3




///
class MeshOptimize3d
{
public:
  void CombineImprove (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY);
  void SplitImprove (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY);
  void SwapImprove (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY,
		    const BitArray * working_elements = NULL);
  void SwapImproveSurface (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY,
			   const BitArray * working_elements = NULL,
			   const ARRAY< ARRAY<int,PointIndex::BASE>* > * idmaps = NULL);
  void SwapImprove2 (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY);
};



extern double CalcBad (const Mesh::T_POINTS & points, const Element & elem,
		       double h);

extern double CalcTotalBad (const Mesh::T_POINTS & points, 
			    const Mesh::T_VOLELEMENTS & elements);

extern int WrongOrientation (const Mesh::T_POINTS & points, const Element & el);


/* Functional depending of inner point inside triangular surface */


class MinFunctionSum : public MinFunction
{
protected:
  ARRAY<MinFunction*> functions;
 
public:
  
  virtual double Func (const Vector & x) const;
  virtual void Grad (const Vector & x, Vector & g) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
  virtual double GradStopping (const Vector & x) const;

  void AddFunction(MinFunction & fun);
  
  const MinFunction & Function(int i) const;
  MinFunction & Function(int i);  
};
  


class PointFunction1 : public MinFunction
{
  Mesh::T_POINTS & points;
  const ARRAY<INDEX_3> & faces;
  double h;
public:
  PointFunction1 (Mesh::T_POINTS & apoints, 
		  const ARRAY<INDEX_3> & afaces,
		  double ah);
  
  virtual double Func (const Vector & x) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double GradStopping (const Vector & x) const;
};


class JacobianPointFunction : public MinFunction
{
public:
  Mesh::T_POINTS & points;
  const Mesh::T_VOLELEMENTS & elements;
  TABLE<INDEX> elementsonpoint;
  PointIndex actpind;

  bool onplane;
  Vec<3> nv;
  
public:
  JacobianPointFunction (Mesh::T_POINTS & apoints, 
			 const Mesh::T_VOLELEMENTS & aelements);
  
  virtual void SetPointIndex (PointIndex aactpind);
  virtual double Func (const Vector & x) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;

  inline void SetNV(const Vec<3> & anv) {nv = anv; onplane = true;}
  inline void UnSetNV(void) {onplane = false;}
};



#endif
