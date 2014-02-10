#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <gprim.hpp>

namespace netgen
{
int
IntersectTriangleLine (const Point<3> ** tri, const Point<3> ** line)
{
  Vec3d vl(*line[0], *line[1]);
  Vec3d vt1(*tri[0], *tri[1]);
  Vec3d vt2(*tri[0], *tri[2]);
  Vec3d vrs(*tri[0], *line[0]);

  static DenseMatrix a(3), ainv(3);
  static Vector rs(3), lami(3);
  int i;

  /*
  (*testout) << "Tri-Line inters: " << endl
	     << "tri = " << *tri[0] << ", " << *tri[1] << ", " << *tri[2] << endl
	     << "line = " << *line[0] << ", " << *line[1] << endl;
  */
  for (i = 1; i <= 3; i++)
    {
      a.Elem(i, 1) = -vl.X(i);
      a.Elem(i, 2) = vt1.X(i);
      a.Elem(i, 3) = vt2.X(i);
      rs.Elem(i) = vrs.X(i);
    }

  double det = a.Det();

  double arel = vl.Length() * vt1.Length() * vt2.Length();
  /*
  double amax = 0;
  for (i = 1; i <= 9; i++)
    if (fabs (a.Get(i)) > amax)
      amax = fabs(a.Get(i));
  */
  // new !!!!
  if (fabs (det) <= 1e-10 * arel)
    {
#ifdef DEVELOP      
      // line parallel to triangle !
      // cout << "ERROR: IntersectTriangleLine degenerated" << endl;
      //      (*testout) << "WARNING: IntersectTriangleLine degenerated\n";
      /*
      (*testout) << "lin-tri intersection: " << endl
		 << "line = " << *line[0] << " - " << *line[1] << endl
		 << "tri = " << *tri[0] << " - " << *tri[1] << " - " << *tri[2] << endl
		 << "lami = " << lami << endl
		 << "pc = " << ( *line[0] + lami.Get(1) * vl ) << endl
		 << "   = " << ( *tri[0] + lami.Get(2) * vt1 + lami.Get(3) * vt2) << endl
		 << " a = " << a << endl
		 << " ainv = " << ainv << endl
		 << " det(a) = " << det << endl
		 << " rs = " << rs << endl;
      */
#endif
      return 0;
    }

  CalcInverse (a, ainv);
  ainv.Mult (rs, lami);

  //  (*testout) << "lami = " << lami << endl;

  double eps = 1e-6;
  if (
      (lami.Get(1) >= -eps && lami.Get(1) <= 1+eps && 
       lami.Get(2) >= -eps && lami.Get(3) >= -eps && 
       lami.Get(2) + lami.Get(3) <= 1+eps)  && !
      (lami.Get(1) >= eps && lami.Get(1) <= 1-eps && 
       lami.Get(2) >= eps && lami.Get(3) >= eps && 
       lami.Get(2) + lami.Get(3) <= 1-eps) )


     {
#ifdef DEVELOP
       //      cout << "WARNING: IntersectTriangleLine degenerated" << endl;
      (*testout) << "WARNING: IntersectTriangleLine numerical inexact" << endl;

      (*testout) << "lin-tri intersection: " << endl
		 << "line = " << *line[0] << " - " << *line[1] << endl
		 << "tri = " << *tri[0] << " - " << *tri[1] << " - " << *tri[2] << endl
		 << "lami = " << lami << endl
		 << "pc = " << ( *line[0] + lami.Get(1) * vl ) << endl
		 << "   = " << ( *tri[0] + lami.Get(2) * vt1 + lami.Get(3) * vt2) << endl
		 << " a = " << a << endl
		 << " ainv = " << ainv << endl
		 << " det(a) = " << det << endl
		 << " rs = " << rs << endl;
#endif
    }
      

  if (lami.Get(1) >= 0 && lami.Get(1) <= 1 && 
      lami.Get(2) >= 0 && lami.Get(3) >= 0 && lami.Get(2) + lami.Get(3) <= 1)
    {

      return 1;
    }

  return 0;
}





int IntersectTetTriangle (const Point<3> ** tet, const Point<3> ** tri,
			  const int * tetpi, const int * tripi)
{
  int i, j;
  double diam = Dist (*tri[0], *tri[1]);
  double epsrel = 1e-8;
  double eps = diam * epsrel;

  double eps2 = eps * eps;
  int cnt = 0;

  int tetp1 = -1, tetp2 = -1;
  int trip1 = -1, trip2 = -1;
  int tetp3, tetp4, trip3;

  /*
  for (i = 0; i < 4; i++)
    loctetpi[i] = -1;
  */


  if (!tetpi)
    {
      for (i = 0; i <= 2; i++)
	{
	  //	  loctripi[i] = -1;
	  for (j = 0; j <= 3; j++)
	    {
	      if (Dist2 (*tet[j], *tri[i]) < eps2)
		{
		  //		  loctripi[i] = j;
		  //		  loctetpi[j] = i;
		  cnt++;
		  tetp2 = tetp1;
		  tetp1 = j;
		  trip2 = trip1;
		  trip1 = i;
		  break;
		}
	    }
	}
    }
  else
    {
      for (i = 0; i <= 2; i++)
	{
	  //	  loctripi[i] = -1;
	  for (j = 0; j <= 3; j++)
	    {
	      if (tetpi[j] == tripi[i])
		{
		  //		  loctripi[i] = j;
		  //		  loctetpi[j] = i;
		  cnt++;
		  tetp2 = tetp1;
		  tetp1 = j;
		  trip2 = trip1;
		  trip1 = i;
		  break;
		}
	    }
	}
    }  
  
  //  (*testout) << "cnt = " << cnt << endl;


  //  (*testout) << "tet-trig inters, cnt = " << cnt << endl;
  
  // cnt .. number of common points
  switch (cnt)
    {
    case 0:
      {
	Vec3d no, n;
	int inpi[3];

	// check, if some trigpoint is in tet:

	for (j = 0; j < 3; j++)
	  inpi[j] = 1;

	for (i = 1; i <= 4; i++)
	  {
	    int pi1 = i % 4;
	    int pi2 = (i+1) % 4;
	    int pi3 = (i+2) % 4;
	    int pi4 = (i+3) % 4;

	    Vec3d v1 (*tet[pi1], *tet[pi2]);
	    Vec3d v2 (*tet[pi1], *tet[pi3]);
	    Vec3d v3 (*tet[pi1], *tet[pi4]);
	    Cross (v1, v2, n);

	    // n /= n.Length();
	    double nl = n.Length();

	    if (v3 * n > 0)
	      n *= -1;

	    int outeri = 1;
	    for (j = 0; j < 3; j++)
	      {
		Vec3d v(*tet[pi1], *tri[j]);
		if ( v * n < eps * nl)
		  outeri = 0;
		else
		  inpi[j] = 0;
	      }

	    if (outeri)
	      return 0;
	  }

	if (inpi[0] || inpi[1] || inpi[2])
	  {
	    return 1;
	  }


	// check, if some tet edge intersects triangle:
	const Point<3> * line[2], *tetf[3];
	for (i = 0; i <= 2; i++)
	  for (j = i+1; j <= 3; j++)
	    {
	      line[0] = tet[i];
	      line[1] = tet[j];

	      if (IntersectTriangleLine (tri, &line[0]))
		return 1;
	    }

	// check, if triangle line intersects tet face:
	for (i = 0; i <= 3; i++)
	  {
	    for (j = 0; j <= 2; j++)
	      tetf[j] = tet[(i+j) % 4];
	    
	    for (j = 0; j <= 2; j++)
	      {
		line[0] = tri[j];
		line[1] = tri[(j+1) % 3];
		
		if (IntersectTriangleLine (&tetf[0], &line[0]))
		  return 1;
	      }
	  }


	return 0;
//GH	break;
      }
    case 1:
      {
	trip2 = 0;
	while (trip2 == trip1)
	  trip2++;
	trip3 = 3 - trip1 - trip2;

	tetp2 = 0;
	while (tetp2 == tetp1)
	  tetp2++;
	tetp3 = 0;
	while (tetp3 == tetp1 || tetp3 == tetp2)
	  tetp3++;
	tetp4 = 6 - tetp1 - tetp2 - tetp3;

	Vec3d vtri1 = *tri[trip2] - *tri[trip1];
	Vec3d vtri2 = *tri[trip3] - *tri[trip1];
	Vec3d ntri;
	Cross (vtri1, vtri2, ntri);

	// tri durch tet ?
	// fehlt noch


	// test 3 tet-faces:
	for (i = 1; i <= 3; i++)
	  {
	    Vec3d vtet1, vtet2;
	    switch (i)
	      {
	      case 1:
		{
		  vtet1 = *tet[tetp2] - *tet[tetp1];
		  vtet2 = *tet[tetp3] - *tet[tetp1];
		  break;
		}
	      case 2:
		{
		  vtet1 = *tet[tetp3] - *tet[tetp1];
		  vtet2 = *tet[tetp4] - *tet[tetp1];
		  break;
		}
	      case 3:
		{
		  vtet1 = *tet[tetp4] - *tet[tetp1];
		  vtet2 = *tet[tetp2] - *tet[tetp1];
		  break;
		}
	      }
	    
	    Vec3d ntet;
	    Cross (vtet1, vtet2, ntet);
	    
	    Vec3d crline = Cross (ntri, ntet);

	    double lcrline = crline.Length();

	    if (lcrline < eps * eps * eps * eps)  // new change !
	      continue;

	    if (vtri1 * crline + vtri2 * crline < 0)
	      crline *= -1;

	    crline /= lcrline;

	    double lam1, lam2, lam3, lam4;
	    LocalCoordinates (vtri1, vtri2, crline, lam1, lam2);
	    LocalCoordinates (vtet1, vtet2, crline, lam3, lam4);
	    
	    if (lam1 > -epsrel && lam2 > -epsrel &&
		lam3 > -epsrel && lam4 > -epsrel)
	      {
		
		/*
		(*testout) << "lcrline = " << lcrline 
			   << " eps = " << eps << " diam = " << diam << endl;
		 
		(*testout) << "hit, cnt == 1 " 
			   << "lam1,2,3,4 = " << lam1 << ", " 
			   << lam2 << ", " << lam3 << ", " << lam4
			   << "\n";
		*/
		return 1;
	      }
	  }
	return 0;
//GH	break;
      }
    case 2:
      {
	// common edge
	tetp3 = 0;
	while (tetp3 == tetp1 || tetp3 == tetp2)
	  tetp3++;
	tetp4 = 6 - tetp1 - tetp2 - tetp3;
	trip3 = 3 - trip1 - trip2;

	//	(*testout) << "trip1,2,3 = " << trip1 << ", " << trip2 << ", " << trip3 << endl;
	//	(*testout) << "tetp1,2,3,4 = " << tetp1 << ", " << tetp2 
	//		   << ", " << tetp3 << ", " << tetp4 << endl;

	Vec3d vtri = *tri[trip3] - *tri[trip1];
	Vec3d vtet1 = *tet[tetp3] - *tri[trip1];
	Vec3d vtet2 = *tet[tetp4] - *tri[trip1];

	Vec3d n = *tri[trip2] - *tri[trip1];
	n /= n.Length();

	vtet1 -= (n * vtet1) * n;
	vtet2 -= (n * vtet2) * n;


	double lam1, lam2;
	LocalCoordinates (vtet1, vtet2, vtri, lam1, lam2);
	
	if (lam1 < -epsrel || lam2 < -epsrel)
	  return 0;
	else
	  {
	    /*

	    (*testout) << "vtet1 = " << vtet1 << endl;
	    (*testout) << "vtet2 = " << vtet2 << endl;
	    (*testout) << "vtri = " << vtri << endl;
	    (*testout) << "lam1 = " << lam1 << " lam2 = " << lam2 << endl;
	    (*testout) << (lam1 * (vtet1 * vtet1) + lam2 * (vtet1 * vtet2))
		       << " = " << (vtet1 * vtri) << endl;
	    (*testout) << (lam1 * (vtet1 * vtet2) + lam2 * (vtet2 * vtet2))
		       << " = " << (vtet2 * vtri) << endl;
	    
	    (*testout) << "tet = ";
	    for (j = 0; j < 4; j++)
	      (*testout) << (*tet[j]) << " ";
	    (*testout) << endl;
	    (*testout) << "tri = ";
	    for (j = 0; j < 3; j++)
	      (*testout) << (*tri[j]) << " ";
	    (*testout) << endl;

	    (*testout) << "hit, cnt == 2" << endl;
	    */
	    
	    return 1;
	  }
	  
	break;
      }
    case 3:
      {
	// common face
	return 0;
      }
    }

  (*testout) << "hit, cnt = " << cnt << endl;
  return 1;
}





int IntersectTetTriangleRef (const Point<3> ** tri, const int * tripi)
{
  int i, j;
  double eps = 1e-8;
  double eps2 = eps * eps;

  static Point<3> rtetp1(0, 0, 0);
  static Point<3> rtetp2(1, 0, 0);  
  static Point<3> rtetp3(0, 1, 0); 
  static Point<3> rtetp4(0, 0, 1);

  static const Point<3> * tet[] = { &rtetp1, &rtetp2, &rtetp3, &rtetp4 };
  static int tetpi[] = { 1, 2, 3, 4 };


  //  return IntersectTetTriangle (tet, tri, tetpi, tripi);

  
  int cnt = 0;

  int tetp1 = -1, tetp2 = -1;
  int trip1 = -1, trip2 = -1;
  int tetp3, tetp4, trip3;

  /*
  if (!tetpi)
    {
      for (i = 0; i <= 2; i++)
	{
	  for (j = 0; j <= 3; j++)
	    {
	      if (Dist2 (*tet[j], *tri[i]) < eps2)
		{
		  cnt++;
		  tetp2 = tetp1;
		  tetp1 = j;
		  trip2 = trip1;
		  trip1 = i;
		  break;
		}
	    }
	}
    }
  else
  */
    {
      for (i = 0; i <= 2; i++)
	{
	  for (j = 0; j <= 3; j++)
	    {
	      if (tetpi[j] == tripi[i])
		{
		  cnt++;
		  tetp2 = tetp1;
		  tetp1 = j;
		  trip2 = trip1;
		  trip1 = i;
		  break;
		}
	    }
	}
    }  
  
  //  (*testout) << "cnt = " << cnt << endl;


  switch (cnt)
    {
    case 0:
      {
	Vec3d no, n;
	//	int inpi[3];
	int pside[3][4];

	for (j = 0; j < 3; j++)
	  {
	    pside[j][0] = (*tri[j])(0) > -eps;
	    pside[j][1] = (*tri[j])(1) > -eps;
	    pside[j][2] = (*tri[j])(2) > -eps;
	    pside[j][3] = (*tri[j])(0) + (*tri[j])(1) + (*tri[j])(2) < 1+eps;
	  }

	
	for (j = 0; j < 4; j++)
	  {
	    if (!pside[0][j] && !pside[1][j] && !pside[2][j])
	      return 0;
	  }

	for (j = 0; j < 3; j++)
	  {
	    if (pside[j][0] && pside[j][1] && pside[j][2] && pside[j][3])
	      return 1;
	  }


	const Point<3> * line[2], *tetf[3];
	for (i = 0; i <= 2; i++)
	  for (j = i+1; j <= 3; j++)
	    {
	      line[0] = tet[i];
	      line[1] = tet[j];

	      if (IntersectTriangleLine (tri, &line[0]))
		return 1;
	    }

	for (i = 0; i <= 3; i++)
	  {
	    for (j = 0; j <= 2; j++)
	      tetf[j] = tet[(i+j) % 4];
	    
	    for (j = 0; j <= 2; j++)
	      {
		line[0] = tri[j];
		line[1] = tri[(j+1) % 3];

	      if (IntersectTriangleLine (&tetf[0], &line[0]))
		return 1;
	      }
	  }


	return 0;
	break;
      }
    case 1:
      {
	trip2 = 0;
	if (trip2 == trip1)
	  trip2++;
	trip3 = 3 - trip1 - trip2;

	tetp2 = 0;
	while (tetp2 == tetp1)
	  tetp2++;
	tetp3 = 0;
	while (tetp3 == tetp1 || tetp3 == tetp2)
	  tetp3++;
	tetp4 = 6 - tetp1 - tetp2 - tetp3;

	Vec3d vtri1 = *tri[trip2] - *tri[trip1];
	Vec3d vtri2 = *tri[trip3] - *tri[trip1];
	Vec3d ntri;
	Cross (vtri1, vtri2, ntri);

	// tri durch tet ?

	/*
	Vec3d vtet1(*tet[tetp1], *tet[tetp2]);
	Vec3d vtet2(*tet[tetp1], *tet[tetp3]);
	Vec3d vtet3(*tet[tetp1], *tet[tetp4]);
	Vec3d sol;
	
	SolveLinearSystem (vtet1, vtet2, vtet3, vtri1, sol);
	if (sol.X() > 0 && sol.Y() > 0 && sol.Z() > 0)
	  return 1;

	SolveLinearSystem (vtet1, vtet2, vtet3, vtri2, sol);
	if (sol.X() > 0 && sol.Y() > 0 && sol.Z() > 0)
	  return 1;
	*/

	// test 3 tet-faces:
	for (i = 1; i <= 3; i++)
	  {
	    Vec3d vtet1, vtet2;
	    switch (i)
	      {
	      case 1:
		{
		  vtet1 = *tet[tetp2] - *tet[tetp1];
		  vtet2 = *tet[tetp3] - *tet[tetp1];
		  break;
		}
	      case 2:
		{
		  vtet1 = *tet[tetp3] - *tet[tetp1];
		  vtet2 = *tet[tetp4] - *tet[tetp1];
		  break;
		}
	      case 3:
		{
		  vtet1 = *tet[tetp4] - *tet[tetp1];
		  vtet2 = *tet[tetp2] - *tet[tetp1];
		  break;
		}
	      }
	    
	    Vec3d ntet;
	    Cross (vtet1, vtet2, ntet);
	    
	    Vec3d crline = Cross (ntri, ntet);

	    double lcrline = crline.Length();
	    if (lcrline < eps * eps)
	      continue;


	    if (vtri1 * crline + vtri2 * crline < 0)
	      crline *= -1;

	    double lam1, lam2, lam3, lam4;
	    LocalCoordinates (vtri1, vtri2, crline, lam1, lam2);
	    LocalCoordinates (vtet1, vtet2, crline, lam3, lam4);
	    
	    if (lam1 > -eps && lam2 > -eps &&
		lam3 > -eps && lam4 > -eps)
	      {
		//		(*testout) << "hit, cnt == 1" << "\n";
		return 1;
	      }
	  }

	return 0;
	break;
      }
    case 2:
      {
	// common edge
	tetp3 = 0;
	while (tetp3 == tetp1 || tetp3 == tetp2)
	  tetp3++;
	tetp4 = 6 - tetp1 - tetp2 - tetp3;
	trip3 = 3 - trip1 - trip2;

	//	(*testout) << "trip1,2,3 = " << trip1 << ", " << trip2 << ", " << trip3 << endl;
	//	(*testout) << "tetp1,2,3,4 = " << tetp1 << ", " << tetp2 
	//		   << ", " << tetp3 << ", " << tetp4 << endl;

	Vec3d vtri = *tri[trip3] - *tri[trip1];
	Vec3d vtet1 = *tet[tetp3] - *tri[trip1];
	Vec3d vtet2 = *tet[tetp4] - *tri[trip1];

	Vec3d n = *tri[trip2] - *tri[trip1];
	n /= n.Length();

	vtet1 -= (n * vtet1) * n;
	vtet2 -= (n * vtet2) * n;


	double lam1, lam2;
	LocalCoordinates (vtet1, vtet2, vtri, lam1, lam2);
	
	if (lam1 < -eps || lam2 < -eps)
	  return 0;
	else
	  {

// 	    (*testout) << "vtet1 = " << vtet1 << endl;
// 	    (*testout) << "vtet2 = " << vtet2 << endl;
// 	    (*testout) << "vtri = " << vtri << endl;
// 	    (*testout) << "lam1 = " << lam1 << " lam2 = " << lam2 << endl;

// 	    (*testout) << (lam1 * (vtet1 * vtet1) + lam2 * (vtet1 * vtet2))
// 		       << " = " << (vtet1 * vtri) << endl;
// 	    (*testout) << (lam1 * (vtet1 * vtet2) + lam2 * (vtet2 * vtet2))
// 		       << " = " << (vtet2 * vtri) << endl;
	    
// 	    (*testout) << "tet = ";
// 	    for (j = 0; j < 4; j++)
// 	      (*testout) << (*tet[j]) << " ";
// 	    (*testout) << endl;
// 	    (*testout) << "tri = ";
// 	    for (j = 0; j < 3; j++)
// 	      (*testout) << (*tri[j]) << " ";
// 	    (*testout) << endl;

// 	    (*testout) << "hit, cnt == 2" << endl;

	    return 1;
	  }
	  
	break;
      }
    case 3:
      {
	// common face
	return 0;
      }
    }

  (*testout) << "hit, cnt = " << cnt << endl;
  return 1;
}











int IntersectTriangleTriangle (const Point<3> ** tri1, const Point<3> ** tri2)
{
  int i, j;
  double diam = Dist (*tri1[0], *tri1[1]);
  double epsrel = 1e-8;
  double eps = diam * epsrel;
  double eps2 = eps * eps;



  int cnt = 0;
  /*
  int tri1pi[3];
  int tri2pi[3];
  */

  //  int tri1p1 = -1; 
  /// int tri1p2 = -1;
  //  int tri2p1 = -1;
  //   int tri2p2 = -1;
  //  int tri1p3, tri2p3;

  /*
  for (i = 0; i < 3; i++)
    tri1pi[i] = -1;
  */
  for (i = 0; i <= 2; i++)
    {
      //      tri2pi[i] = -1;
      for (j = 0; j <= 2; j++)
	{
	  if (Dist2 (*tri1[j], *tri2[i]) < eps2)
	    {
	      //	      tri2pi[i] = j;
	      //	      tri1pi[j] = i;
	      cnt++;
	      //	      tri1p2 = tri1p1;
	      //	      tri1p1 = j;
	      //	      tri2p2 = tri2p1;
	      //	      tri2p1 = i;
	      break;
	    }
	}
    }
  
  switch (cnt)
    {
    case 0:
      {
	const Point<3> * line[2];
	
	for (i = 0; i <= 2; i++)
	  {
	    line[0] = tri2[i];
	    line[1] = tri2[(i+1)%3];

	    if (IntersectTriangleLine (tri1, &line[0]))
	      {
		(*testout) << "int1, line = " << *line[0] << " - " << *line[1] << endl;
		return 1;
	      }
	  }	

	for (i = 0; i <= 2; i++)
	  {
	    line[0] = tri1[i];
	    line[1] = tri1[(i+1)%3];

	    if (IntersectTriangleLine (tri2, &line[0]))
	      {
		(*testout) << "int2, line = " << *line[0] << " - " << *line[1] << endl;
		return 1;
	      }
	  }	
	break;
      }
    default:
      return 0;
    }

  return 0;
}



void
LocalCoordinates (const Vec3d & e1, const Vec3d & e2,
		  const Vec3d & v, double & lam1, double & lam2)
{
  double m11 = e1 * e1;
  double m12 = e1 * e2;
  double m22 = e2 * e2;
  double rs1 = v * e1;
  double rs2 = v * e2;
  
  double det = m11 * m22 - m12 * m12;
  lam1 = (rs1 * m22 - rs2 * m12)/det;
  lam2 = (m11 * rs2 - m12 * rs1)/det;
}





int CalcSphereCenter (const Point<3> ** pts, Point<3> & c)
{
  Vec3d row1 (*pts[0], *pts[1]);
  Vec3d row2 (*pts[0], *pts[2]);
  Vec3d row3 (*pts[0], *pts[3]);

  Vec3d rhs(0.5 * (row1*row1),
	    0.5 * (row2*row2),
	    0.5 * (row3*row3));
  Transpose (row1, row2, row3);
  
  Vec3d sol;
  if (SolveLinearSystem (row1, row2, row3, rhs, sol))
    {
      (*testout) << "CalcSphereCenter: degenerated" << endl;
      return 1;
    }

  c = *pts[0] + sol;
  return 0;
}





int CalcTriangleCenter (const Point3d ** pts, Point3d & c)
{
  static DenseMatrix a(2), inva(2);
  static Vector rs(2), sol(2);
  double h = Dist(*pts[0], *pts[1]);

  Vec3d v1(*pts[0], *pts[1]);
  Vec3d v2(*pts[0], *pts[2]);

  rs.Elem(1) = v1 * v1;
  rs.Elem(2) = v2 * v2;

  a.Elem(1,1) = 2 * rs.Get(1);
  a.Elem(1,2) = a.Elem(2,1) = 2 * (v1 * v2);
  a.Elem(2,2) = 2 * rs.Get(2);

  if (fabs (a.Det()) <= 1e-12 * h * h)
    {
      (*testout) << "CalcTriangleCenter: degenerated" << endl;
      return 1;
    }

  CalcInverse (a, inva);
  inva.Mult (rs, sol);

  c = *pts[0];
  v1 *= sol.Get(1);
  v2 *= sol.Get(2);

  c += v1;
  c += v2;

  return 0;
}



double ComputeCylinderRadius (const Point3d & p1, 
			      const Point3d & p2,
			      const Point3d & p3, 
			      const Point3d & p4)
{
  Vec3d v12(p1, p2);
  Vec3d v13(p1, p3);
  Vec3d v14(p1, p4);

  Vec3d n1 = Cross (v12, v13);
  Vec3d n2 = Cross (v14, v12);
		
  double n1l = n1.Length();
  double n2l = n2.Length();
  n1 /= n1l;
  n2 /= n2l;

  double v12len = v12.Length();
  double h1 = n1l / v12len;
  double h2 = n2l / v12len;
  
  /*
  (*testout) << "n1 = " << n1 << " n2 = " << n2 
	     << "h1 = " << h1 << " h2 = " << h2 << endl;
  */
  return ComputeCylinderRadius (n1, n2, h1, h2);
}




/*
  Two triangles T1 and T2 have normals n1 and n2.
  The height over the common edge is h1, and h2.
 */
double ComputeCylinderRadius (const Vec3d & n1, const Vec3d & n2,
				     double h1, double h2)
{
  Vec3d t1, t2;
  double n11 = n1 * n1;
  double n12 = n1 * n2;
  double n22 = n2 * n2;
  double det = n11 * n22 - n12 * n12;
  
  if (fabs (det) < 1e-14 * n11 * n22)
    return 1e20;

  // a biorthogonal bases   (ti * nj) = delta_ij:
  t1 = (n22/det) * n1 + (-n12/det) * n2;
  t2 = (-n12/det) * n1 + (n11/det) * n2;

  // normalize:
  t1 /= t1.Length();
  t2 /= t2.Length();

  /*
    vector to center point has form
    v = lam1 n1 + lam2 n2
    and fulfills
    t2 v = h1/2
    t1 v = h2/2
  */

  double lam1 = 0.5 * h2 / (n1 * t1);
  double lam2 = 0.5 * h1 / (n2 * t2);
  
  double rad = (lam1 * n1 + lam2 * n2).Length();
  /*
  (*testout) << "n1 = " << n1
	     << " n2 = " << n2
	     << " t1 = " << t1
	     << " t2 = " << t2
	     << " rad = " << rad << endl;
  */
  return rad;
}
    





double MinDistLP2 (const Point2d & lp1, const Point2d & lp2, const Point2d & p)
{
  Vec2d v(lp1, lp2);
  Vec2d vlp(lp1, p);

  // dist(lam) = \| vlp \|^2 - 2 lam (v1p, v) + lam^2 \| v \|^2

  // lam = (v * vlp) / (v * v);
  // if (lam < 0) lam = 0;
  // if (lam > 1) lam = 1;

  double num = v*vlp;
  double den = v*v;

  if (num <= 0) 
    return Dist2 (lp1, p);

  if (num >= den) 
    return Dist2 (lp2, p);
  
  if (den > 0)
    {
      return vlp.Length2() - num * num /den;
    }
  else
    return vlp.Length2();
}




double MinDistLP2 (const Point3d & lp1, const Point3d & lp2, const Point3d & p)
{
  Vec3d v(lp1, lp2);
  Vec3d vlp(lp1, p);

  // dist(lam) = \| vlp \|^2 - 2 lam (v1p, v) + lam^2 \| v \|^2

  // lam = (v * vlp) / (v * v);
  // if (lam < 0) lam = 0;
  // if (lam > 1) lam = 1;

  double num = v*vlp;
  double den = v*v;

  if (num <= 0) 
    return Dist2 (lp1, p);

  if (num >= den) 
    return Dist2 (lp2, p);
  
  if (den > 0)
    {
      return vlp.Length2() - num * num /den;
    }
  else
    return vlp.Length2();
}



double MinDistTP2 (const Point3d & tp1, const Point3d & tp2, 
		   const Point3d & tp3, const Point3d & p)
{
  double lam1, lam2;
  double res;

  LocalCoordinates (Vec3d (tp1, tp2), Vec3d (tp1, tp3),
		    Vec3d (tp1, p), lam1, lam2);
  int in1 = lam1 >= 0;
  int in2 = lam2 >= 0;
  int in3 = lam1+lam2 <= 1;
  
  if (in1 && in2 && in3)
    {
      Point3d pp = tp1 + lam1 * Vec3d(tp1, tp2) + lam2 *  Vec3d (tp1, tp3);
      res = Dist2 (p, pp);
    }
  else
    {
      res = Dist2 (tp1, p);
      if (!in1)
	{
	  double hv = MinDistLP2 (tp1, tp3, p);
	  if (hv < res) res = hv; 
	}
      if (!in2)
	{
	  double hv = MinDistLP2 (tp1, tp2, p);
	  if (hv < res) res = hv; 
	}
      if (!in3)
	{
	  double hv = MinDistLP2 (tp2, tp3, p);
	  if (hv < res) res = hv; 
	}
      /*
      double d1 = MinDistLP2 (tp1, tp2, p);
      double d2 = MinDistLP2 (tp1, tp3, p);
      double d3 = MinDistLP2 (tp2, tp3, p);
      res = min3 (d1, d2, d3);
      */
    }

  return res;

  Vec3d pp1(tp1, p);
  Vec3d v1(tp1, tp2), v2(tp1, tp3);

  double c = pp1.Length2();
  double cx = -2 * (pp1 * v1);
  double cy = -2 * (pp1 * v2);
  double cxx = v1.Length2();
  double cxy = 2 * (v1 * v2);
  double cyy = v2.Length2();

  QuadraticPolynomial2V pol (-c, -cx, -cy, -cxx, -cxy, -cyy);
  double res2 =  - pol.MaxUnitTriangle ();

  if (fabs (res - res2) > 1e-8)
    cout << "res and res2 differ: " << res << " != " << res2 << endl;
  return res2;
}


// 0 checks !!!
double MinDistLL2 (const Point3d & l1p1, const Point3d & l1p2,
		  const Point3d & l2p1, const Point3d & l2p2)
{
  // dist(lam1,lam2) = \| l2p1+lam2v2 - (l1p1+lam1 v1) \|
  // min !

  Vec3d l1l2 (l1p1, l2p1);
  Vec3d v1 (l1p1, l1p2);
  Vec3d v2 (l2p1, l2p2);

  double a11, a12, a22, rs1, rs2;
  double lam1, lam2, det;

  a11 = v1*v1;
  a12 = -(v1*v2);
  a22 = v2*v2;
  rs1 = l1l2 * v1;
  rs2 = - (l1l2 * v2);
  
  det = a11 * a22 - a12 * a12;
  if (det < 1e-14 * a11 * a22) 
    det = 1e-14 * a11 * a22;  // regularization should be stable

  if (det < 1e-20)
    det = 1e-20;


  lam1 = (a22 * rs1 - a12 * rs2) / det;
  lam2 = (-a12 * rs1 + a11 * rs2) / det;

  if (lam1 >= 0 && lam2 >= 0 && lam1 <= 1 && lam2 <= 1)
    {
      Vec3d v = l1l2 + (-lam1) * v1 + lam2 * v2;
      return v.Length2();
    }

  double minv, hv;
  minv = MinDistLP2 (l1p1, l1p2, l2p1);
  hv =  MinDistLP2 (l1p1, l1p2, l2p2);
  if (hv < minv) minv = hv;

  hv =  MinDistLP2 (l2p1, l2p2, l1p1);
  if (hv < minv) minv = hv;
  hv =  MinDistLP2 (l2p1, l2p2, l1p2);
  if (hv < minv) minv = hv;

  return minv;
}
			 
}


