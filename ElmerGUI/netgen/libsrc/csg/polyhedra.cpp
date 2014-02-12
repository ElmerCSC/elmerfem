#include <mystdlib.h>

#include <linalg.hpp>
#include <csg.hpp>

namespace netgen
{

Polyhedra::Face::Face (int pi1, int pi2, int pi3,
		       const ARRAY<Point<3> > & points,
		       int ainputnr)
{
  inputnr = ainputnr;

  pnums[0] = pi1;
  pnums[1] = pi2;
  pnums[2] = pi3;


  bbox.Set (points[pi1]);
  bbox.Add (points[pi2]);
  bbox.Add (points[pi3]);

  v1 = points[pi2] - points[pi1];
  v2 = points[pi3] - points[pi1];

  n = Cross (v1, v2);

  nn = n;
  nn.Normalize();
  //  PseudoInverse (v1, v2, w1, w2);
  
  Mat<2,3> mat;
  Mat<3,2> inv;
  for (int i = 0; i < 3; i++)
    {
      mat(0,i) = v1(i);
      mat(1,i) = v2(i);
    }
  CalcInverse (mat, inv);
  for (int i = 0; i < 3; i++)
    {
      w1(i) = inv(i,0);
      w2(i) = inv(i,1);
    }
}


Polyhedra :: Polyhedra ()
{
  surfaceactive.SetSize(0);
  surfaceids.SetSize(0);
  eps_base1 = 1e-8;
}

Polyhedra :: ~Polyhedra ()
{
  ;
}

Primitive * Polyhedra :: CreateDefault ()
{
  return new Polyhedra();
}

INSOLID_TYPE Polyhedra :: BoxInSolid (const BoxSphere<3> & box) const
{
  /*
  for (i = 1; i <= faces.Size(); i++)
    if (FaceBoxIntersection (i, box))
      return DOES_INTERSECT;
  */
  for (int i = 0; i < faces.Size(); i++)
    {
      if (!faces[i].bbox.Intersect (box))
	continue;
      //(*testout) << "face " << i << endl;

      const Point<3> & p1 = points[faces[i].pnums[0]];
      const Point<3> & p2 = points[faces[i].pnums[1]];
      const Point<3> & p3 = points[faces[i].pnums[2]];

      if (fabs (faces[i].nn * (p1 - box.Center())) > box.Diam()/2)
	continue;

      //(*testout) << "still in loop" << endl;

      double dist2 = MinDistTP2 (p1, p2, p3, box.Center());
      //(*testout) << "p1 " << p1 << " p2 " << p2 << " p3 " << p3 << endl
      //		 << " box.Center " << box.Center() << " box.Diam() " << box.Diam() << endl
      //	 << " dist2 " << dist2 << " sqr(box.Diam()/2) " << sqr(box.Diam()/2) << endl;
      if (dist2 < sqr (box.Diam()/2))
	{
	  //(*testout) << "DOES_INTERSECT" << endl;
	  return DOES_INTERSECT;
	}
    };

  return PointInSolid (box.Center(), 1e-3 * box.Diam());
}


INSOLID_TYPE Polyhedra :: PointInSolid (const Point<3> & p,
					double eps) const
{
  //(*testout) << "PointInSolid p " << p << " eps " << eps << endl;
  //(*testout) << "bbox " << poly_bbox << endl;

  if((p(0) > poly_bbox.PMax()(0) + eps) || (p(0) < poly_bbox.PMin()(0) - eps) ||
     (p(1) > poly_bbox.PMax()(1) + eps) || (p(1) < poly_bbox.PMin()(1) - eps) ||
     (p(2) > poly_bbox.PMax()(2) + eps) || (p(2) < poly_bbox.PMin()(2) - eps))
    {
      //(*testout) << "returning IS_OUTSIDE" << endl;
      return IS_OUTSIDE;
    }

  Vec<3> n, v1, v2;

  // random (?) numbers:
  n(0) = -0.424621;
  n(1) = 0.15432;
  n(2) = 0.89212238;

  int cnt = 0;

  for (int i = 0; i < faces.Size(); i++)
    {
      const Point<3> & p1 = points[faces[i].pnums[0]];
      
      Vec<3> v0 = p - p1;

      double lam3 = faces[i].nn * v0;

      if(fabs(lam3) < eps)
	{
	  double lam1 = (faces[i].w1 * v0);
	  double lam2 = (faces[i].w2 * v0);
	  if (lam1 >= -eps_base1 && lam2 >= -eps_base1 && lam1+lam2 <= 1+eps_base1)
	    {
	      //(*testout) << "returning DOES_INTERSECT" << endl;
	      return DOES_INTERSECT;
	    }
	}
      else
	{
	  lam3 = -(faces[i].n * v0) / (faces[i].n * n);

	  if (lam3 < 0) continue;

	  Vec<3> rs = v0 + lam3 * n;
	  
	  double lam1 = (faces[i].w1 * rs);
	  double lam2 = (faces[i].w2 * rs);
	  if (lam1 >= 0 && lam2 >= 0 && lam1+lam2 <= 1)
	    cnt++;
	}

    }

  //(*testout) << " cnt = " << cnt%2 << endl;
  return (cnt % 2) ? IS_INSIDE : IS_OUTSIDE;
}




void Polyhedra :: GetTangentialSurfaceIndices (const Point<3> & p, 
					       ARRAY<int> & surfind, double eps) const
{
  for (int i = 0; i < faces.Size(); i++)
    {
      const Point<3> & p1 = points[faces[i].pnums[0]];
      
      Vec<3> v0 = p - p1;
      double lam3 = -(faces[i].nn * v0); // n->nn

      if (fabs (lam3) > eps) continue;

      double lam1 = (faces[i].w1 * v0);
      double lam2 = (faces[i].w2 * v0);

      if (lam1 >= -eps_base1 && lam2 >= -eps_base1 && lam1+lam2 <= 1+eps_base1)
	if (!surfind.Contains (GetSurfaceId(i)))
	  surfind.Append (GetSurfaceId(i));
    }

}



INSOLID_TYPE Polyhedra :: VecInSolid (const Point<3> & p,
				      const Vec<3> & v,
				      double eps) const
{
  ARRAY<int> point_on_faces;
  INSOLID_TYPE res(DOES_INTERSECT);

  Vec<3> vn = v;
  vn.Normalize();
  for (int i = 0; i < faces.Size(); i++)
    {
      const Point<3> & p1 = points[faces[i].pnums[0]];
      
      Vec<3> v0 = p - p1;
      double lam3 = -(faces[i].nn * v0); // n->nn 


      if (fabs (lam3) > eps) continue;
      //(*testout) << "lam3 <= eps" << endl;

      double lam1 = (faces[i].w1 * v0);
      double lam2 = (faces[i].w2 * v0);

      if (lam1 >= -eps_base1 && lam2 >= -eps_base1 && lam1+lam2 <= 1+eps_base1)
	{
	  point_on_faces.Append(i);

	  double scal = vn * faces[i].nn; // n->nn
	
	  res = DOES_INTERSECT;
	  if (scal > eps_base1) res = IS_OUTSIDE;
	  if (scal < -eps_base1) res = IS_INSIDE;
	}
    }
  
  //(*testout) << "point_on_faces.Size() " << point_on_faces.Size() 
  //	     << " res " << res << endl;

  if (point_on_faces.Size() == 0)
    return PointInSolid (p, 0);
  if (point_on_faces.Size() == 1)
    return res;



  
  double mindist(0);
  bool first = true;

  for(int i=0; i<point_on_faces.Size(); i++)
    {
      for(int j=0; j<3; j++)
	{
	  double dist = Dist(p,points[faces[point_on_faces[i]].pnums[j]]);
	  if(dist > eps && (first || dist < mindist))
	    {
	      mindist = dist;
	      first = false;
	    }
	}
    }
  
  Point<3> p2 = p + (1e-2*mindist) * vn;
  res = PointInSolid (p2, eps);

  //  (*testout) << "mindist " << mindist << " res " << res << endl;

  return res;
  
 
}


/*
INSOLID_TYPE Polyhedra :: VecInSolid2 (const Point<3> & p,
				       const Vec<3> & v1,
				       const Vec<3> & v2,
				       double eps) const
{
  INSOLID_TYPE res;

  res = VecInSolid(p,v1,eps);
  if(res != DOES_INTERSECT)
    return res;

  int point_on_n_faces = 0;

  Vec<3> v1n = v1;
  v1n.Normalize();
  Vec<3> v2n = v2;
  v2n.Normalize();


  for (int i = 0; i < faces.Size(); i++)
    {
      const Point<3> & p1 = points[faces[i].pnums[0]];
      
      Vec<3> v0 = p - p1;
      double lam3 = -(faces[i].n * v0);

      if (fabs (lam3) > eps) continue;

      double lam1 = (faces[i].w1 * v0);
      double lam2 = (faces[i].w2 * v0);

      if (lam1 >= -eps && lam2 >= -eps && lam1+lam2 <= 1+eps)
	{
	  double scal1 = v1n * faces[i].n;
	  if (fabs (scal1) > eps) continue;


	  point_on_n_faces++;

	  double scal2 = v2n * faces[i].n;
	  res = DOES_INTERSECT;
	  if (scal2 > eps) res = IS_OUTSIDE;
	  if (scal2 < -eps) res = IS_INSIDE;
	}
    }

  if (point_on_n_faces == 1)
    return res;

  cerr << "primitive::vecinsolid2 makes nonsense for polyhedra" << endl;

  return Primitive :: VecInSolid2 (p, v1, v2, eps);
}
*/



INSOLID_TYPE Polyhedra :: VecInSolid2 (const Point<3> & p,
				       const Vec<3> & v1,
				       const Vec<3> & v2,
				       double eps) const
{
  //(*testout) << "VecInSolid2 eps " << eps << endl;
  INSOLID_TYPE res = VecInSolid(p,v1,eps);
  //(*testout) << "VecInSolid = " <<res <<endl;

  if(res != DOES_INTERSECT)
    return res;

  int point_on_n_faces = 0;

  Vec<3> v1n = v1;
  v1n.Normalize();
  Vec<3> v2n = v2 - (v2 * v1n) * v1n;
  v2n.Normalize();

  double cosv2, cosv2max = -1;

  
  for (int i = 0; i < faces.Size(); i++)
    {
      const Point<3> & p1 = points[faces[i].pnums[0]];
      
      Vec<3> v0 = p - p1;
      if (fabs (faces[i].nn * v0) > eps) continue; // n->nn
      if (fabs (v1n * faces[i].nn) > eps_base1) continue; // n->nn

      double lam1 = (faces[i].w1 * v0);
      double lam2 = (faces[i].w2 * v0);

      if (lam1 >= -eps_base1 && lam2 >= -eps_base1 && lam1+lam2 <= 1+eps_base1)
	{
	  // v1 is in face

	  Point<3> fc = Center (points[faces[i].pnums[0]],
				points[faces[i].pnums[1]],
				points[faces[i].pnums[2]]);

	  Vec<3> vpfc = fc - p;
	  cosv2 = (v2n * vpfc) / vpfc.Length();
	  if (cosv2 > cosv2max)
	    {
	      cosv2max = cosv2;
	      point_on_n_faces++;

	      double scal2 = v2n * faces[i].nn; // n->nn
	      res = DOES_INTERSECT;
	      if (scal2 > eps_base1) res = IS_OUTSIDE;
	      if (scal2 < -eps_base1) res = IS_INSIDE;

	    }
	}
    }

  if (point_on_n_faces >= 1)
    return res;

  (*testout) << "primitive::vecinsolid2 makes nonsense for polyhedra" << endl;
  cerr << "primitive::vecinsolid2 makes nonsense for polyhedra" << endl;

  return Primitive :: VecInSolid2 (p, v1, v2, eps);
}



void Polyhedra :: GetTangentialVecSurfaceIndices2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
						   ARRAY<int> & surfind, double eps) const
{
  Vec<3> v1n = v1;
  v1n.Normalize();
  Vec<3> v2n = v2; //  - (v2 * v1n) * v1n;
  v2n.Normalize();


  for (int i = 0; i < faces.Size(); i++)
    {
      const Point<3> & p1 = points[faces[i].pnums[0]];
      
      Vec<3> v0 = p - p1;
      if (fabs (v0 * faces[i].nn) > eps) continue; // n->nn
      if (fabs (v1n * faces[i].nn) > eps_base1) continue; // n->nn
      if (fabs (v2n * faces[i].nn) > eps_base1) continue; // n->nn

      double lam01 = (faces[i].w1 * v0);
      double lam02 = (faces[i].w2 * v0);
      double lam03 = 1-lam01-lam02;
      double lam11 = (faces[i].w1 * v1);
      double lam12 = (faces[i].w2 * v1);
      double lam13 = -lam11-lam12;
      double lam21 = (faces[i].w1 * v2);
      double lam22 = (faces[i].w2 * v2);
      double lam23 = -lam21-lam22;

      bool ok1 = lam01 > eps_base1 ||
	lam01 > -eps_base1 && lam11 > eps_base1 ||
	lam01 > -eps_base1 && lam11 > -eps_base1 && lam21 > eps_base1;

      bool ok2 = lam02 > eps_base1 ||
	lam02 > -eps_base1 && lam12 > eps_base1 ||
	lam02 > -eps_base1 && lam12 > -eps_base1 && lam22 > eps_base1;
      
      bool ok3 = lam03 > eps_base1 ||
	lam03 > -eps_base1 && lam13 > eps_base1 ||
	lam03 > -eps_base1 && lam13 > -eps_base1 && lam23 > eps_base1;

      if (ok1 && ok2 && ok3)
	{
	  if (!surfind.Contains (GetSurfaceId(faces[i].planenr)))
	    surfind.Append (GetSurfaceId(faces[i].planenr));
	}
    }  
}












void Polyhedra :: GetPrimitiveData (const char *& classname, 
				    ARRAY<double> & coeffs) const
{
  classname = "Polyhedra";
  coeffs.SetSize(0);
  coeffs.Append (points.Size());
  coeffs.Append (faces.Size());
  coeffs.Append (planes.Size());

  /*
  int i, j;
  for (i = 1; i <= planes.Size(); i++)
    {
      planes.Elem(i)->Print (*testout);
    }
  for (i = 1; i <= faces.Size(); i++)
    {
      (*testout) << "face " << i << " has plane " << faces.Get(i).planenr << endl;
      for (j = 1; j <= 3; j++)
	(*testout) << points.Get(faces.Get(i).pnums[j-1]);
      (*testout) << endl;
    }
  */
}

void Polyhedra :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  ;
}

void Polyhedra :: Reduce (const BoxSphere<3> & box)
{
  for (int i = 0; i < planes.Size(); i++)
    surfaceactive[i] = 0;

  for (int i = 0; i < faces.Size(); i++)
    if (FaceBoxIntersection (i, box))
      surfaceactive[faces[i].planenr] = 1;
}

void Polyhedra :: UnReduce ()
{
  for (int i = 0; i < planes.Size(); i++)
    surfaceactive[i] = 1;
}

int Polyhedra :: AddPoint (const Point<3> & p)
{
  if(points.Size() == 0)
    poly_bbox.Set(p);
  else
    poly_bbox.Add(p);

  return points.Append (p);
}

int Polyhedra :: AddFace (int pi1, int pi2, int pi3, int inputnum)
{
  (*testout) << "polyhedra, add face " << pi1 << ", " << pi2 << ", " << pi3 << endl;

  if(pi1 == pi2 || pi2 == pi3 || pi3 == pi1)
    {
      ostringstream msg;
      msg << "Illegal point numbers for polyhedron face: " << pi1+1 << ", " << pi2+1 << ", " << pi3+1;
      throw NgException(msg.str());
    }

  faces.Append (Face (pi1, pi2, pi3, points, inputnum));
  
  Point<3> p1 = points[pi1];
  Point<3> p2 = points[pi2];
  Point<3> p3 = points[pi3];

  Vec<3> v1 = p2 - p1;
  Vec<3> v2 = p3 - p1;

  Vec<3> n = Cross (v1, v2); 
  n.Normalize();

  Plane pl (p1, n);
//   int inverse;
//   int identicto = -1;
//   for (int i = 0; i < planes.Size(); i++)
//     if (pl.IsIdentic (*planes[i], inverse, 1e-9*max3(v1.Length(),v2.Length(),Dist(p2,p3))))
//       {
// 	if (!inverse)
// 	  identicto = i;
//       }
//   //  cout << "is identic = " << identicto << endl;
//   identicto = -1;    // changed April 10, JS

//   if (identicto != -1)
//     faces.Last().planenr = identicto;
//   else
    {
      planes.Append (new Plane (p1, n));
      surfaceactive.Append (1);
      surfaceids.Append (0);
      faces.Last().planenr = planes.Size()-1;
    }

//  (*testout) << "is plane nr " << faces.Last().planenr << endl;

  return faces.Size();
}



int Polyhedra :: FaceBoxIntersection (int fnr, const BoxSphere<3> & box) const
{
  /*
  (*testout) << "check face box intersection, fnr = " << fnr << endl;
  (*testout) << "box = " << box << endl;
  (*testout) << "face-box = " << faces[fnr].bbox << endl;
  */

  if (!faces[fnr].bbox.Intersect (box))
    return 0;

  const Point<3> & p1 = points[faces[fnr].pnums[0]];
  const Point<3> & p2 = points[faces[fnr].pnums[1]];
  const Point<3> & p3 = points[faces[fnr].pnums[2]];

  double dist2 = MinDistTP2 (p1, p2, p3, box.Center());
  /*
  (*testout) << "p1 = " << p1 << endl;
  (*testout) << "p2 = " << p2 << endl;
  (*testout) << "p3 = " << p3 << endl;

  (*testout) << "box.Center() = " << box.Center() << endl;
  (*testout) << "center = " << box.Center() << endl;
  (*testout) << "dist2 = " << dist2 << endl;
  (*testout) << "diam = " << box.Diam() << endl;
  */
  if (dist2 < sqr (box.Diam()/2))
    {
      //      (*testout) << "intersect" << endl;
      return 1;
    }
  return 0;
}


void Polyhedra :: GetPolySurfs(ARRAY < ARRAY<int> * > & polysurfs)
{
  int maxnum = -1;
  
  for(int i = 0; i<faces.Size(); i++)
    {
      if(faces[i].inputnr > maxnum)
	maxnum = faces[i].inputnr;
    }
  
  polysurfs.SetSize(maxnum+1);
  for(int i=0; i<polysurfs.Size(); i++)
    polysurfs[i] = new ARRAY<int>;

  for(int i = 0; i<faces.Size(); i++)
    polysurfs[faces[i].inputnr]->Append(faces[i].planenr);
}


void Polyhedra::CalcSpecialPoints (ARRAY<Point<3> > & pts) const
{
  for (int i = 0; i < points.Size(); i++)
    pts.Append (points[i]);
}


void Polyhedra :: AnalyzeSpecialPoint (const Point<3> & pt, 
				       ARRAY<Point<3> > & specpts) const
{
  ;
}

Vec<3> Polyhedra :: SpecialPointTangentialVector (const Point<3> & p, int s1, int s2) const
{
  const double eps = 1e-10*poly_bbox.Diam();

  for (int fi1 = 0; fi1 < faces.Size(); fi1++)
    for (int fi2 = 0; fi2 < faces.Size(); fi2++)
      {
	int si1 = faces[fi1].planenr;
	int si2 = faces[fi2].planenr;

	if (surfaceids[si1] != s1 || surfaceids[si2] != s2) continue;

	//(*testout) << "check pair fi1/fi2 " << fi1 << "/" << fi2 << endl;
	
	Vec<3> n1 = GetSurface(si1) . GetNormalVector (p);
	Vec<3> n2 = GetSurface(si2) . GetNormalVector (p);
	Vec<3> t = Cross (n1, n2);

	//(*testout) << "t = " << t << endl;


	/*
	int samepts = 0;
	for (int j = 0; j < 3; j++)
	  for (int k = 0; k < 3; k++)
	    if (Dist(points[faces[fi1].pnums[j]],
		     points[faces[fi2].pnums[k]]) < eps)
	      samepts++;
	if (samepts < 2) continue;
	*/

	bool shareedge = false;
	for(int j = 0; !shareedge && j < 3; j++)
	  {
	    Vec<3> v1 = points[faces[fi1].pnums[(j+1)%3]] - points[faces[fi1].pnums[j]];
	    double smax = v1.Length();
	    v1 *= 1./smax;
	    
	    int pospos;
	    if(fabs(v1(0)) > 0.5)
	      pospos = 0;
	    else if(fabs(v1(1)) > 0.5)
	      pospos = 1;
	    else
	      pospos = 2;

	    double sp = (p(pospos) - points[faces[fi1].pnums[j]](pospos)) / v1(pospos);
	    if(sp < -eps || sp > smax+eps)
	      continue;

	    for (int k = 0; !shareedge && k < 3; k ++)
	      {
		 Vec<3> v2 = points[faces[fi2].pnums[(k+1)%3]] - points[faces[fi2].pnums[k]];
		 v2.Normalize();
		 if(v2 * v1 > 0)
		   v2 -= v1;
		 else
		   v2 += v1;
		 
		 //(*testout) << "v2.Length2() " << v2.Length2() << endl;

		 if(v2.Length2() > 1e-18)
		   continue;

		 double sa,sb;

		 sa = (points[faces[fi2].pnums[k]](pospos) - points[faces[fi1].pnums[j]](pospos)) / v1(pospos);
		 sb = (points[faces[fi2].pnums[(k+1)%3]](pospos) - points[faces[fi1].pnums[j]](pospos)) / v1(pospos);
		 

		 if(Dist(points[faces[fi1].pnums[j]] + sa*v1, points[faces[fi2].pnums[k]]) > eps)
		   continue;

		 if(sa > sb)
		   {
		     double aux = sa; sa = sb; sb = aux;
		   }

		 //testout->precision(20);
		 //(*testout) << "sa " << sa << " sb " << sb << " smax " << smax << " sp " << sp  << " v1 " << v1 << endl;
		 //testout->precision(8);


		 shareedge = (sa < -eps && sb > eps) ||
		   (sa < smax-eps && sb > smax+eps) ||
		   (sa > -eps && sb < smax+eps);

		 if(!shareedge)
		   continue;

		 sa = max2(sa,0.);
		 sb = min2(sb,smax);

		 if(sp < sa+eps)
		   shareedge = (t * v1 > 0);
		 else if (sp > sb-eps)
		   shareedge = (t * v1 < 0);
		   
	      }
	  }
	if (!shareedge) continue;

	t.Normalize();
	  
	
	return t;
      }

  return Vec<3> (0,0,0);
}


}


