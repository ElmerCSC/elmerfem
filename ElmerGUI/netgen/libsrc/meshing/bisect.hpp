#ifndef BISECT
#define BISECT

class BisectionOptions
{
public:
  const char * outfilename;
  const char * mlfilename;
  const char * refinementfilename;
  const char * femcode;
  int maxlevel;
  int usemarkedelements;
  bool refine_hp;
  bool refine_p;
  BisectionOptions ();
};

class ZRefinementOptions
{
public:
  int minref;
  ZRefinementOptions();
};


/*
extern void BisectTets (Mesh &, const CSGeometry *,
			BisectionOptions & opt);
*/

extern void BisectTetsCopyMesh (Mesh &, const class CSGeometry *,
				BisectionOptions & opt);

extern void ZRefinement (Mesh &, const CSGeometry *,
			 ZRefinementOptions & opt);





class Refinement
{
  MeshOptimize2d * optimizer2d;

public:
  Refinement ();
  virtual ~Refinement ();
  
  void Refine (Mesh & mesh);
  void Bisect (Mesh & mesh, class BisectionOptions & opt, ARRAY<double> * quality_loss = NULL);
  void MakeSecondOrder (Mesh & mesh);

  virtual void PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint, 
			     int surfi, 
			     const PointGeomInfo & gi1, 
			     const PointGeomInfo & gi2,
			     Point<3> & newp, PointGeomInfo & newgi);

  virtual void PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
			     int surfi1, int surfi2, 
			     const EdgePointGeomInfo & ap1, 
			     const EdgePointGeomInfo & ap2,
			     Point<3> & newp, EdgePointGeomInfo & newgi);

  virtual Vec<3> GetTangent (const Point<3> & p, int surfi1, int surfi2,
                             const EdgePointGeomInfo & egi) const;

  virtual Vec<3> GetNormal (const Point<3> & p, int surfi1, 
                            const PointGeomInfo & gi) const;


  virtual void ProjectToSurface (Point<3> & p, int surfi);

  virtual void ProjectToSurface (Point<3> & p, int surfi, const PointGeomInfo & /* gi */)
  {
    ProjectToSurface (p, surfi);
  }

  virtual void ProjectToEdge (Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const;


  void ValidateSecondOrder (Mesh & mesh);
  void ValidateRefinedMesh (Mesh & mesh, 
			    ARRAY<INDEX_2> & parents);

  MeshOptimize2d * Get2dOptimizer(void)
  {
    return optimizer2d;
  }
  void Set2dOptimizer(MeshOptimize2d * opti)
  {
    optimizer2d = opti;
  }

  
  virtual void LocalizeEdgePoints(Mesh & mesh) const {;}
};

#endif
