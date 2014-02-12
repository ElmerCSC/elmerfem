#ifndef FILE_GEOM2DMESH
#define FILE_GEOM2DMESH

/**************************************************************************/
/* File:   geom2dmesh.hh                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   22. Jan. 01                                                    */
/**************************************************************************/


class Refinement2d : public Refinement
{
  const SplineGeometry2d & geometry;

public:
  Refinement2d (const SplineGeometry2d & ageometry);
  virtual ~Refinement2d ();
  
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

  virtual void ProjectToSurface (Point<3> & p, int surfi, const PointGeomInfo & /* gi */);

  virtual void ProjectToEdge (Point<3> & p, int surfi1, int surfi2, 
                              const EdgePointGeomInfo & egi) const;			     
};






#endif
