#ifndef FILE_MESHSURF
#define FILE_MESHSURF

///
class Meshing2Surfaces : public Meshing2
{
  ///
  const Surface & surface;
  
public:
  ///
  //  Meshing2Surfaces (const Surface & asurf);
  ///
  Meshing2Surfaces (const Surface & asurf, const Box<3> & aboundingbox);

protected:
  ///
  virtual void DefineTransformation (const Point3d & p1, const Point3d & p2,
				     const PointGeomInfo * geominfo1,
				     const PointGeomInfo * geominfo2);
  ///
  virtual void TransformToPlain (const Point3d & locpoint, 
				 const MultiPointGeomInfo & geominfo,
				 Point2d & plainpoint, 
				 double h, int & zone);
  ///
  virtual int TransformFromPlain (Point2d & plainpoint,
				  Point3d & locpoint, 
				  PointGeomInfo & gi,
				  double h);
  ///
  virtual double CalcLocalH (const Point3d & p, double gh) const;
};



///
class MeshOptimize2dSurfaces : public MeshOptimize2d
  {
  ///
  const CSGeometry & geometry;

public:
    ///
    MeshOptimize2dSurfaces (const CSGeometry & ageometry); 
   
    ///
    virtual void ProjectPoint (INDEX surfind, Point<3> & p) const;
    ///
    virtual void ProjectPoint2 (INDEX surfind, INDEX surfind2, Point<3> & p) const;
    ///
    virtual void GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const;
};





class RefinementSurfaces : public Refinement
{
  const CSGeometry & geometry;

public:
  RefinementSurfaces (const CSGeometry & ageometry);
  virtual ~RefinementSurfaces ();
  
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
                             const EdgePointGeomInfo & ap1) const;

  virtual Vec<3> GetNormal (const Point<3> & p, int surfi1, 
                            const PointGeomInfo & gi) const;


  virtual void ProjectToSurface (Point<3> & p, int surfi);

  virtual void ProjectToEdge (Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const;

};



#endif

