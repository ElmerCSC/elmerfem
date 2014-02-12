#include <mystdlib.h>

#include <csg.hpp>
#include <meshing.hpp>



namespace netgen
{
  /*
Meshing2Surfaces :: Meshing2Surfaces (const Surface & asurface)
  : surface(asurface)
{
  ;
}
  */
Meshing2Surfaces :: Meshing2Surfaces (const Surface & asurf,
				      const Box<3> & abb)
  : Meshing2(abb), surface(asurf)
{
  ;
}


void Meshing2Surfaces :: DefineTransformation (const Point3d & p1, const Point3d & p2,
					       const PointGeomInfo * geominfo1,
					       const PointGeomInfo * geominfo2)
{
  ((Surface&)surface).DefineTangentialPlane (p1, p2);
}

void Meshing2Surfaces :: TransformToPlain (const Point3d & locpoint, 
					   const MultiPointGeomInfo & geominfo,
					   Point2d & planepoint, 
					   double h, int & zone)
{
  Point<2> hp;
  surface.ToPlane (locpoint, hp, h, zone);
  planepoint.X() = hp(0);
  planepoint.Y() = hp(1);
}

int Meshing2Surfaces :: TransformFromPlain (Point2d & planepoint,
					    Point3d & locpoint, 
					    PointGeomInfo & gi,
					    double h)
{
  Point<3> hp;
  Point<2> hp2 (planepoint.X(), planepoint.Y());
  surface.FromPlane (hp2, hp, h);
  locpoint = hp;
  gi.trignum = 1;
  return 0;
}



double Meshing2Surfaces :: CalcLocalH (const Point3d & p, double gh) const
{
  return surface.LocH (p, 3, 1, gh);
  /*
    double loch = mesh.lochfunc->GetH(p);
    if (gh < loch) loch = gh;
    return loch;
    */
}






MeshOptimize2dSurfaces :: MeshOptimize2dSurfaces (const CSGeometry & ageometry)
  : MeshOptimize2d(), geometry(ageometry)
{
  ;
}


void MeshOptimize2dSurfaces :: ProjectPoint (INDEX surfind, Point<3> & p) const
{
  Point<3> hp = p;
  geometry.GetSurface(surfind)->Project (hp);
  p = hp;
}

void MeshOptimize2dSurfaces :: ProjectPoint2 (INDEX surfind, INDEX surfind2, 
					      Point<3> & p) const
{
  Point<3> hp = p;
  ProjectToEdge ( geometry.GetSurface(surfind), 
		  geometry.GetSurface(surfind2), hp);
  p = hp;
}


void MeshOptimize2dSurfaces :: 
GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const
{
  Vec<3> hn = n;
  geometry.GetSurface(surfind)->CalcGradient (p, hn);
  hn.Normalize();
  n = hn;

  /*
  if (geometry.GetSurface(surfind)->Inverse())
    n *= -1;
  */
}
  






RefinementSurfaces :: RefinementSurfaces (const CSGeometry & ageometry)
  : Refinement(), geometry(ageometry)
{
  if(geometry.GetNSurf() == 0)
    cerr << endl 
	 << "WARNING: Intializing 2D refinement with 0-surface geometry" << endl
	 << "==========================================================" << endl
	 << endl << endl;
}

RefinementSurfaces :: ~RefinementSurfaces ()
{
  ;
}
  
void RefinementSurfaces :: 
PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
	      int surfi, 
	      const PointGeomInfo & gi1, 
	      const PointGeomInfo & gi2,
	      Point<3> & newp, PointGeomInfo & newgi)
{
  Point<3> hnewp;
  hnewp = p1+secpoint*(p2-p1);
  if (surfi != -1)
    {
      geometry.GetSurface (surfi) -> Project (hnewp);
      newgi.trignum = 1;
    }

  newp = hnewp;
}

void RefinementSurfaces :: 
PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
	      int surfi1, int surfi2, 
	      const EdgePointGeomInfo & ap1, 
	      const EdgePointGeomInfo & ap2,
	      Point<3> & newp, EdgePointGeomInfo & newgi)
{
  Point<3> hnewp = p1+secpoint*(p2-p1);
  //(*testout) << "hnewp " << hnewp << " s1 " << surfi1 << " s2 " << surfi2 << endl;
  if (surfi1 != -1 && surfi2 != -1 && surfi1 != surfi2)
    {
      netgen::ProjectToEdge (geometry.GetSurface(surfi1), 
                             geometry.GetSurface(surfi2), 
                             hnewp);
      // (*testout) << "Pointbetween, newp = " << hnewp << endl
      // << ", err = " << sqrt (sqr (hnewp(0))+ sqr(hnewp(1)) + sqr (hnewp(2))) - 1 << endl;
      newgi.edgenr = 1;
      //(*testout) << "hnewp (a1) " << hnewp << endl;
    }
  else if (surfi1 != -1)
    {
      geometry.GetSurface (surfi1) -> Project (hnewp);
      //(*testout) << "hnewp (a2) " << hnewp << endl;
    }

  newp = hnewp;
};

Vec<3> RefinementSurfaces :: GetTangent (const Point<3> & p, int surfi1, int surfi2,
                                         const EdgePointGeomInfo & ap1) const
{
  Vec<3> n1 = geometry.GetSurface (surfi1)->GetNormalVector (p);
  Vec<3> n2 = geometry.GetSurface (surfi2)->GetNormalVector (p);
  Vec<3> tau = Cross (n1, n2).Normalize();
  return tau;
}

Vec<3> RefinementSurfaces :: GetNormal (const Point<3> & p, int surfi1, 
                                        const PointGeomInfo & gi) const
{
  return geometry.GetSurface (surfi1)->GetNormalVector (p);
}



void RefinementSurfaces :: ProjectToSurface (Point<3> & p, int surfi)
{
  if (surfi != -1)
    geometry.GetSurface (surfi) -> Project (p);
};

void RefinementSurfaces :: ProjectToEdge (Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const
{
  netgen::ProjectToEdge (geometry.GetSurface(surfi1), 
                         geometry.GetSurface(surfi2), 
                         p);
  
}


}
