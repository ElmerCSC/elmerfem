#include <mystdlib.h>

#include <csg.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>

namespace netgen
{

  Refinement2d :: Refinement2d (const SplineGeometry2d & ageometry)
    : Refinement(), geometry(ageometry)
  {
    ;
  }

  Refinement2d :: ~Refinement2d ()
  {
    ;
  }
  

  void Refinement2d :: 
  PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
		int surfi, 
		const PointGeomInfo & gi1, 
		const PointGeomInfo & gi2,
		Point<3> & newp, PointGeomInfo & newgi)
  {
    newp = p1+secpoint*(p2-p1);
    newgi.trignum = 1;
  }



  void Refinement2d :: 
  PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint, 
		int surfi1, int surfi2, 
		const EdgePointGeomInfo & ap1, 
		const EdgePointGeomInfo & ap2,
		Point<3> & newp, EdgePointGeomInfo & newgi)
  {
    Point<2> p2d;
  
    p2d = geometry.GetSplines().Get(ap1.edgenr) -> 
      GetPoint (((1-secpoint)*ap1.dist+secpoint*ap2.dist));
  
    //  (*testout) << "refine 2d line, ap1.dist, ap2.dist = " << ap1.dist << ", " << ap2.dist << endl;
    //  (*testout) << "p1, p2 = " << p1 << p2 << ", newp = " << p2d << endl;

    newp = Point3d (p2d(0), p2d(1), 0);
    newgi.edgenr = ap1.edgenr;
    newgi.dist = ((1-secpoint)*ap1.dist+secpoint*ap2.dist);
  };



  Vec<3> Refinement2d :: GetTangent (const Point<3> & p, int surfi1, int surfi2,
                                     const EdgePointGeomInfo & ap1) const
  {
    Vec<2> t2d = geometry.GetSplines().Get(ap1.edgenr) -> GetTangent(ap1.dist);
    return Vec<3> (t2d(0), t2d(1), 0);
  }

  Vec<3> Refinement2d :: GetNormal (const Point<3> & p, int surfi1, 
                                    const PointGeomInfo & gi) const
  {
    return Vec<3> (0,0,1);
  }


  void Refinement2d :: ProjectToSurface (Point<3> & p, int surfi, const PointGeomInfo & /* gi */)
  {
    p(2) = 0;
  }


  void Refinement2d :: ProjectToEdge (Point<3> & p, int surfi1, int surfi2, 
                                      const EdgePointGeomInfo & egi) const
  {
    Point<2> p2d (p(0), p(1)), pp;
    double t;
    geometry.GetSplines().Get(egi.edgenr) -> Project (p2d, pp, t);
    p = Point<3> (pp(0), pp(1), 0);
  }
}
