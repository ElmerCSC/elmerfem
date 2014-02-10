#include <mystdlib.h>

#include "meshing.hpp"
#include <csg.hpp>
#include <geometry2d.hpp>

namespace netgen
{

  int CheckSurfaceMesh (const Mesh & mesh)
  {
    PrintMessage (3, "Check Surface mesh");

    int nf = mesh.GetNSE();
    INDEX_2_HASHTABLE<int> edges(nf+2);
    int i, j;
    INDEX_2 i2;
    int cnt1 = 0, cnt2 = 0;

    for (i = 1; i <= nf; i++)
      for (j = 1; j <= 3; j++)
	{
	  i2.I1() = mesh.SurfaceElement(i).PNumMod(j);
	  i2.I2() = mesh.SurfaceElement(i).PNumMod(j+1);
	  if (edges.Used(i2))
	    {
	      int hi;
	      hi = edges.Get(i2);
	      if (hi != 1) 
		PrintSysError ("CheckSurfaceMesh, hi = ", hi);
	      edges.Set(i2, 2);
	      cnt2++;
	    }
	  else
	    {
	      Swap (i2.I1(), i2.I2());
	      edges.Set(i2, 1);
	      cnt1++;
	    }
	}
  

    if (cnt1 != cnt2)
      {
	PrintUserError ("Surface mesh not consistent");
	//      MyBeep(2);
	//      (*mycout) << "cnt1 = " << cnt1 << " cnt2 = " << cnt2 << endl;
	return 0;
      }
    return 1;
  }



  int CheckSurfaceMesh2 (const Mesh & mesh)
  {
    int i, j, k;
    const Point<3> *tri1[3], *tri2[3];

    for (i = 1; i <= mesh.GetNOpenElements(); i++)
      {
	PrintDot ();
	for (j = 1; j < i; j++)
	  {
	    for (k = 1; k <= 3; k++)
	      {
		tri1[k-1] = &mesh.Point (mesh.OpenElement(i).PNum(k));
		tri2[k-1] = &mesh.Point (mesh.OpenElement(j).PNum(k));
	      }
	    if (IntersectTriangleTriangle (&tri1[0], &tri2[0]))
	      {
		PrintSysError ("Surface elements are intersecting");
		(*testout) << "Intersecting: " << endl;
		for (k = 0; k <= 2; k++)
		  (*testout) << *tri1[k] << "   ";
		(*testout) << endl;
		for (k = 0; k <= 2; k++)
		  (*testout) << *tri2[k] << "   ";
		(*testout) << endl;
	      }

	  }
      }
    return 0;
  }





  static double TriangleQualityInst (const Point3d & p1, const Point3d & p2,
				     const Point3d & p3)
  {
    // quality 0 (worst) .. 1 (optimal)

    Vec3d v1, v2, v3;
    double s1, s2, s3;
    double an1, an2, an3;

    v1 = p2 - p1;
    v2 = p3 - p1;
    v3 = p3 - p2;

    an1 = Angle (v1, v2);
    v1 *= -1;
    an2 = Angle (v1, v3);
    an3 = Angle (v2, v3);

    s1 = sin (an1/2);
    s2 = sin (an2/2);
    s3 = sin (an3/2);

    return 8 * s1 * s2 * s3;
  }














  void MeshQuality2d (const Mesh & mesh)
  {
    int ncl = 20, cl;
    ARRAY<INDEX> incl(ncl);
    INDEX i;
    SurfaceElementIndex sei;
    double qual;

    incl = 0;

    for (sei = 0; sei < mesh.GetNSE(); sei++)
      {
	qual = TriangleQualityInst (mesh[mesh[sei][0]],
				    mesh[mesh[sei][1]],
				    mesh[mesh[sei][2]]);

	cl = int ( (ncl-1e-3) * qual ) + 1;
	incl.Elem(cl)++;
      }

    (*testout) << endl << endl;

    (*testout) << "Points:           " << mesh.GetNP() << endl;
    (*testout) << "Surface Elements: " << mesh.GetNSE() << endl;

    (*testout) << endl;
    (*testout) << "Elements in qualityclasses:" << endl;
    (*testout).precision(2);
    for (i = 1; i <= ncl; i++)
      {
	(*testout) << setw(4) << double (i-1)/ncl << " - "
		   << setw(4) << double (i) / ncl << ": "
		   << incl.Get(i) << endl;
      }
  }


  static double TetElementQuality (const Point3d & p1, const Point3d & p2,
				   const Point3d & p3, const Point3d & p4)
  {
    double vol, l, l4, l5, l6;


    Vec3d v1 = p2 - p1;
    Vec3d v2 = p3 - p1;
    Vec3d v3 = p4 - p1;

    vol = fabs ((Cross (v1, v2) * v3)) / 6;
    l4 = Dist (p2, p3);
    l5 = Dist (p2, p4);
    l6 = Dist (p3, p4);

    l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

    if (vol <= 1e-8 * l * l * l) return 1e-10;

    return vol/(l*l*l) * 1832.82;    // 6^4 * sqrt(2)
  }





  double teterrpow = 2;

  double CalcTetBadness (const Point3d & p1, const Point3d & p2,
			 const Point3d & p3, const Point3d & p4, double h)
  {
    double vol, l, ll, lll, ll1, ll2, ll3, ll4, ll5, ll6;
    double err;

    Vec3d v1 (p1, p2);
    Vec3d v2 (p1, p3);
    Vec3d v3 (p1, p4);

    vol = -Determinant (v1, v2, v3) / 6;

    ll1 = v1.Length2();
    ll2 = v2.Length2();
    ll3 = v3.Length2();
    ll4 = Dist2 (p2, p3);
    ll5 = Dist2 (p2, p4);
    ll6 = Dist2 (p3, p4);

    ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
    l = sqrt (ll);
    lll = l * ll;

    if (vol <= 1e-24 * lll)
      return 1e24;

    err = 0.0080187537 * lll / vol;    // sqrt(216) / (6^4 * sqrt(2))

    if (h > 0)
      err += ll / (h * h) + 
	h * h * ( 1 / ll1 + 1 / ll2 + 1 / ll3 + 
		  1 / ll4 + 1 / ll5 + 1 / ll6 ) - 12;
    
    if (teterrpow == 2)
      return err*err;
    return pow (err, teterrpow);
  }


  double CalcTetBadnessGrad (const Point3d & p1, const Point3d & p2,
			     const Point3d & p3, const Point3d & p4, double h,
			     int pi, Vec<3> & grad)
  {
    double vol, l, ll, lll;
    double err;

    const Point3d *pp1, *pp2, *pp3, *pp4;

    pp1 = &p1;
    pp2 = &p2;
    pp3 = &p3;
    pp4 = &p4;
  
    switch (pi)
      {
      case 2:
	{
	  swap (pp1, pp2);
	  swap (pp3, pp4);
	  break;
	}
      case 3:
	{
	  swap (pp1, pp3);
	  swap (pp2, pp4);
	  break;
	}
      case 4:
	{
	  swap (pp1, pp4);
	  swap (pp3, pp2);
	  break;
	}
      }
  

    Vec3d v1 (*pp1, *pp2);
    Vec3d v2 (*pp1, *pp3);
    Vec3d v3 (*pp1, *pp4);

    Vec3d v4 (*pp2, *pp3);
    Vec3d v5 (*pp2, *pp4);
    Vec3d v6 (*pp3, *pp4);

    vol = -Determinant (v1, v2, v3) / 6;  

    Vec3d gradvol;
    Cross (v5, v4, gradvol);
    gradvol *= (-1.0/6.0);


    double ll1 = v1.Length2();
    double ll2 = v2.Length2();
    double ll3 = v3.Length2();
    double ll4 = v4.Length2();
    double ll5 = v5.Length2();
    double ll6 = v6.Length2();

    ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
    l = sqrt (ll);
    lll = l * ll;

    if (vol <= 1e-24 * lll)
      { 
	grad = Vec3d (0, 0, 0);
	return 1e24;
      }



    Vec3d gradll1 (*pp2, *pp1);
    Vec3d gradll2 (*pp3, *pp1);
    Vec3d gradll3 (*pp4, *pp1);
    gradll1 *= 2;
    gradll2 *= 2;
    gradll3 *= 2;

    Vec3d gradll (gradll1);
    gradll += gradll2;
    gradll += gradll3;

    /*
    Vec3d gradll;
    gradll = v1+v2+v3;
    gradll *= -2;
    */

    err = 0.0080187537 * lll / vol; 


    gradll *= (0.0080187537 * 1.5 * l / vol);
    Vec3d graderr(gradll);
    gradvol *= ( -0.0080187537 * lll / (vol * vol) );
    graderr += gradvol;
  
    if (h > 0)
      {
	/*
	Vec3d gradll1 (*pp2, *pp1);
	Vec3d gradll2 (*pp3, *pp1);
	Vec3d gradll3 (*pp4, *pp1);
	gradll1 *= 2;
	gradll2 *= 2;
	gradll3 *= 2;
	*/
	err += ll / (h*h) + 
	  h*h * ( 1 / ll1 + 1 / ll2 + 1 / ll3 + 
		  1 / ll4 + 1 / ll5 + 1 / ll6 ) - 12;

	graderr += (1/(h*h) - h*h/(ll1*ll1)) * gradll1;
	graderr += (1/(h*h) - h*h/(ll2*ll2)) * gradll2;
	graderr += (1/(h*h) - h*h/(ll3*ll3)) * gradll3;
	cout << "?";
      }


    double errpow;
    if (teterrpow == 2)
      {
        errpow = err*err;   
        grad = (2 * err) * graderr;
      }
    else
      {
        errpow = pow (err, teterrpow);
        grad = (teterrpow * errpow / err) * graderr;
      }
    return errpow;
  }
  




  /*

  double CalcTetBadness (const Point3d & p1, const Point3d & p2,
  const Point3d & p3, const Point3d & p4, double h)
  {
  double vol, l;
  double err;


  Vec3d v1 (p1, p2);
  Vec3d v2 (p1, p3);
  Vec3d v3 (p1, p4);

  vol = -Determinant (v1, v2, v3) / 6;

  double l1 = v1.Length();
  double l2 = v2.Length();
  double l3 = v3.Length();
  double l4 = Dist (p2, p3);
  double l5 = Dist (p2, p4);
  double l6 = Dist (p3, p4);

  l = l1 + l2 + l3 + l4 + l5 + l6;

  // just for timing
  // l += 1e-40 * CalcTetBadnessNew (p1, p2, p3, p4, h);

  if (vol <= 1e-24 * l * l * l)
  { 
  return 1e24;
  }

  err = (l*l*l) / (1832.82 * vol);    // 6^4 * sqrt(2)
  
  if (h > 0)
  err += l / h + 
  h * (1 / l1 + 1/l2 + 1/l3 + 1/l4 + 1/l5 + 1/l6) - 12;

  return pow (err, teterrpow);
  }


  
  double CalcTetBadnessGrad (const Point3d & p1, const Point3d & p2,
  const Point3d & p3, const Point3d & p4, double h,
  int pi, Vec3d & grad)
  {
  double vol, l;
  double err;

  const Point3d *pp1, *pp2, *pp3, *pp4;

  pp1 = &p1;
  pp2 = &p2;
  pp3 = &p3;
  pp4 = &p4;
  
  switch (pi)
  {
  case 2:
  {
  swap (pp1, pp2);
  swap (pp3, pp4);
  break;
  }
  case 3:
  {
  swap (pp1, pp3);
  swap (pp2, pp4);
  break;
  }
  case 4:
  {
  swap (pp1, pp4);
  swap (pp3, pp2);
  break;
  }
  }
  

  Vec3d v1 (*pp1, *pp2);
  Vec3d v2 (*pp1, *pp3);
  Vec3d v3 (*pp1, *pp4);

  Vec3d v4 (*pp2, *pp3);
  Vec3d v5 (*pp2, *pp4);
  Vec3d v6 (*pp3, *pp4);


  //   Vec3d n;
  //   Cross (v1, v2, n);
  //   vol = - (n * v3) / 6;


  vol = -Determinant (v1, v2, v3) / 6;  

  Vec3d gradvol;
  Cross (v5, v4, gradvol);
  gradvol *= (-1.0/6.0);


  double l1 = v1.Length();
  double l2 = v2.Length();
  double l3 = v3.Length();
  double l4 = v4.Length();
  double l5 = v5.Length();
  double l6 = v6.Length();

  l = l1 + l2 + l3 +l4 + l5 + l6;

  Vec3d gradl1 (*pp2, *pp1);
  Vec3d gradl2 (*pp3, *pp1);
  Vec3d gradl3 (*pp4, *pp1);
  gradl1 /= l1;
  gradl2 /= l2;
  gradl3 /= l3;

  Vec3d gradl (gradl1);
  gradl += gradl2;
  gradl += gradl3;


  if (vol <= 1e-24 * l * l * l)
  { 
  grad = Vec3d (0, 0, 0);
  return 1e24;
  }


  double c1 = 1.0 / 1832.82;      // 6^4 * sqrt(2)
  err = c1 * (l*l*l) / vol; 


  gradl *= (c1 * 3 * l * l / vol);
  Vec3d graderr(gradl);
  gradvol *= ( -c1 * l * l * l / (vol * vol) );
  graderr+= gradvol;
  
  if (h > 0)
  {
  err += l / h + 
  h * ( 1 / l1 + 1 / l2 + 1 / l3 + 
  1 / l4 + 1 / l5 + 1 / l6 ) - 12;

  graderr += (1/h - h/(l1*l1)) * gradl1;
  graderr += (1/h - h/(l2*l2)) * gradl2;
  graderr += (1/h - h/(l3*l3)) * gradl3;
  cout << "?";
  }

  double errpow = pow (err, teterrpow);
  grad = (teterrpow * errpow / err) * graderr;
  
  return errpow;
  }
  
  */




  
  /*
    double CalcVolume (const ARRAY<Point3d> & points,
    const Element & el)
    {
    Vec3d v1 = points.Get(el.PNum(2)) - 
    points.Get(el.PNum(1));
    Vec3d v2 = points.Get(el.PNum(3)) - 
    points.Get(el.PNum(1));
    Vec3d v3 = points.Get(el.PNum(4)) - 
    points.Get(el.PNum(1)); 
         
    return -(Cross (v1, v2) * v3) / 6;	 
    }  
  */

  double CalcVolume (const ARRAY<Point3d> & points, 
		     const ARRAY<Element> & elements)
  {
    double vol;
    Vec3d v1, v2, v3;
  
    vol = 0;
    for (int i = 0; i < elements.Size(); i++)
      {
	v1 = points.Get(elements[i][1]) - points.Get(elements[i][0]);
	v2 = points.Get(elements[i][2]) - points.Get(elements[i][0]);
	v3 = points.Get(elements[i][3]) - points.Get(elements[i][0]);
	vol -= (Cross (v1, v2) * v3) / 6;	 
      }
    return vol;
  }

  
  

  void MeshQuality3d (const Mesh & mesh, ARRAY<int> * inclass)
  { 
    int ncl = 20;
    signed int cl;
    ARRAY<INDEX> incl(ncl);
    INDEX i;
    double qual;
    double sum = 0;
    int nontet  = 0;

    for (i = 1; i <= incl.Size(); i++)
      incl.Elem(i) = 0;

    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      {
	if (mesh[ei].GetType() != TET)
	  {
	    nontet++;
	    continue;
	  }

	qual = TetElementQuality (mesh.Point(mesh[ei][0]),
				  mesh.Point(mesh[ei][1]),
				  mesh.Point(mesh[ei][2]),
				  mesh.Point(mesh[ei][3]));

	if (qual > 1) qual = 1;
	cl = int (ncl * qual ) + 1;
     
	if (cl < 1) cl = 1; 
	if (cl > ncl) cl = ncl;

	incl.Elem(cl)++;
	if (inclass) (*inclass)[ei] = cl;
	sum += 1/qual;
      }

    (*testout) << endl << endl;
    (*testout) << "Points:           " << mesh.GetNP() << endl;
    (*testout) << "Volume Elements:  " << mesh.GetNE() << endl;
    if (nontet)
      (*testout) << nontet << " non tetrahedral elements" << endl;
    (*testout) << endl;

    (*testout) << "Volume elements in qualityclasses:" << endl;
    (*testout).precision(2);
    for (i = 1; i <= ncl; i++)
      {
	(*testout) << setw(4) << double (i-1)/ncl << " - "
		   << setw(4) << double (i) / ncl << ": "
		   << incl.Get(i) << endl;
      }
    (*testout) << "total error: " << sum << endl;
  }


  void SaveEdges (const Mesh & mesh, const char * geomfile, double h, char * filename)
  {
    ofstream of (filename);
    int i;
    const Segment * seg;
  
    of << "edges" << endl;
    of << geomfile << endl;
    of << h << endl;

    of << mesh.GetNP() << endl;
    for (i = 1; i <= mesh.GetNP(); i++)
      of << mesh.Point(i)(0) << " "
	 << mesh.Point(i)(1) << " "
	 << mesh.Point(i)(2) << "\n";
    
    of << 2 * mesh.GetNSeg() << endl;
    for (i = 1; i <= mesh.GetNSeg(); i++)
      {
	seg = &mesh.LineSegment(i);

	of << seg->p2 << " " << seg->p1 << " " << seg->si << "\n";
      }
   
  }


  void SaveSurfaceMesh (const Mesh & mesh,
			double h,
			char * filename)

  {
    INDEX i;

    ofstream outfile(filename);

    outfile << "surfacemesh" << endl;
    outfile << h << endl;

    outfile << mesh.GetNP() << endl;
    for (i = 1; i <= mesh.GetNP(); i++)
      outfile << mesh.Point(i)(0) << " "
	      << mesh.Point(i)(1) << " "
	      << mesh.Point(i)(2) << endl;

  

    outfile << mesh.GetNSE() << endl;
    for (i = 1; i <= mesh.GetNSE(); i++)
      {
	const Element2d & el = mesh.SurfaceElement(i);

	if (mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
	  outfile << mesh.SurfaceElement(i).PNum(1) << " "
		  << mesh.SurfaceElement(i).PNum(2) << " "
		  << mesh.SurfaceElement(i).PNum(3) << endl;
	if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0)
	  outfile << mesh.SurfaceElement(i).PNum(1) << " "
		  << mesh.SurfaceElement(i).PNum(3) << " "
		  << mesh.SurfaceElement(i).PNum(2) << endl;
      }
  }


#ifdef OLD
  void Save2DMesh (
		   const Mesh & mesh2d,
		   const ARRAY<SplineSegment *> * splines,
		   ostream & outfile)

  {
    int i, j;
    outfile.precision (6);
  
    outfile << "areamesh2" << endl;


    outfile << endl;
    outfile << mesh2d.GetNSeg() << endl;
    for (i = 1; i <= mesh2d.GetNSeg(); i++)
      outfile << mesh2d.LineSegment(i).si << "        "
	      << mesh2d.LineSegment(i).p1 << " "
	      << mesh2d.LineSegment(i).p2 << "  " << endl;
  

    outfile << mesh2d.GetNSE() << endl;
    for (i = 1; i <= mesh2d.GetNSE(); i++)
      {
	outfile << mesh2d.SurfaceElement(i).GetIndex() << "         ";
	outfile << mesh2d.SurfaceElement(i).GetNP() << " ";
	for (j = 1; j <= mesh2d.SurfaceElement(i).GetNP(); j++)
	  outfile << mesh2d.SurfaceElement(i).PNum(j) << " ";
	outfile << endl;
      }

    outfile << mesh2d.GetNP() << endl;
    for (i = 1; i <= mesh2d.GetNP(); i++)
      outfile << mesh2d.Point(i).X() << " "
	      << mesh2d.Point(i).Y() << endl;

    if (splines)
      {
	outfile << splines->Size() << endl;
	for (i = 1; i <= splines->Size(); i++)
	  splines->Get(i) -> PrintCoeff (outfile);
      }
    else
      outfile << "0" << endl;
  }
#endif








  void SaveVolumeMesh (const Mesh & mesh, 
		       const CSGeometry & geometry,
		       char * filename)
  {
    INDEX i;

    ofstream outfile(filename);
    outfile << "volumemesh" << endl;

    outfile << mesh.GetNSE() << endl;
    for (i = 1; i <= mesh.GetNSE(); i++)
      {
	if (mesh.SurfaceElement(i).GetIndex())
	  outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex ()).SurfNr()
		  << "\t";
	else
	  outfile << "0" << "\t";
	outfile << mesh.SurfaceElement(i)[0] << " "
		<< mesh.SurfaceElement(i)[1] << " "
		<< mesh.SurfaceElement(i)[2] << endl;
      }
    outfile << mesh.GetNE() << endl;
    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      outfile << mesh[ei].GetIndex() << "\t"
	      << mesh[ei][0] << " " << mesh[ei][1] << " "
	      << mesh[ei][2] << " " << mesh[ei][3] << endl;

    outfile << mesh.GetNP() << endl;
    for (i = 1; i <= mesh.GetNP(); i++)
      outfile << mesh.Point(i)(0) << " "
	      << mesh.Point(i)(1) << " "
	      << mesh.Point(i)(2) << endl;

#ifdef SOLIDGEOM
    outfile << geometry.GetNSurf() << endl;
    for (i = 1; i <= geometry.GetNSurf(); i++)
      geometry.GetSurface(i) -> Print (outfile);
#endif
  }




  int CheckCode ()
  {
    return 1;

    /*
      char st[100];
      ifstream ist("pw");

      if (!ist.good()) return 0;
      ist >> st;
      if (strcmp (st, "JKULinz") == 0) return 1;
      return 0;
    */
  }



  /* ******************** CheckMesh ******************************* */

  /// Checks, whether mesh contains a valid 3d mesh
  int CheckMesh3D (const Mesh & mesh)
  {
    INDEX_3_HASHTABLE<int> faceused(mesh.GetNE()/3);
    INDEX i;
    int j, k, l;
    INDEX_3 i3;
    int ok = 1;
    ElementIndex ei;

    for (i = 1; i <= mesh.GetNSE(); i++)
      {
	const Element2d & el = mesh.SurfaceElement(i);
      
	if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0 ||
	    mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
	  {
	    for (j = 1; j <= 3; j++)
	      i3.I(j) = el.PNum(j);
	  
	    i3.Sort();
	    faceused.Set (i3, 1);
	  }
      }
  
    for (ei = 0; ei < mesh.GetNE(); ei++)
      {
	const Element & el = mesh[ei];

	for (j = 1; j <= 4; j++)
	  {
	    l = 0;
	    for (k = 1; k <= 4; k++)
	      {
		if (j != k)
		  {
		    l++;
		    i3.I(l) = el.PNum(k);
		  }
	      }

	    i3.Sort();
	    if (faceused.Used(i3))
	      faceused.Set(i3, faceused.Get(i3)+1);
	    else
	      faceused.Set (i3, 1);
	  }
      }


    for (i = 1; i <= mesh.GetNSE(); i++)
      {
	const Element2d & el = mesh.SurfaceElement(i);

	for (j = 1; j <= 3; j++)
	  i3.I(j) = el.PNum(j);
      
	i3.Sort();
	k = faceused.Get (i3);
	if (k != 2)
	  {
	    ok = 0;
	    (*testout) << "face " << i << " with points " 
		       << i3.I1() << "-" << i3.I2() << "-" << i3.I3() 
		       << " has " << k << " elements" << endl;
	  }
      }
  
    for (ei = 0; ei < mesh.GetNE(); ei++)
      {
	const Element & el = mesh[ei];

	for (j = 1; j <= 4; j++)
	  {
	    l = 0;
	    for (k = 1; k <= 4; k++)
	      {
		if (j != k)
		  {
		    l++;
		    i3.I(l) = el.PNum(k);
		  }
	      }

	    i3.Sort();
	    k = faceused.Get(i3);
	    if (k != 2)
	      {
		ok = 0;
		(*testout) << "element " << ei << " with face " 
			   << i3.I1() << "-" << i3.I2() << "-"
			   << i3.I3() 
			   << " has " << k << " elements" << endl;
	      }
	  }
      }





    /*
      for (i = 1; i <= faceused.GetNBags(); i++)
      for (j = 1; j <= faceused.GetBagSize(i); j++)
      {
      faceused.GetData(i, j, i3, k);
      if (k != 2)
      {
      (*testout) << "Face: " << i3.I1() << "-" 
      << i3.I2() << "-" << i3.I3() << " has " 
      << k << " Faces " << endl;
      cerr << "Face Error" << endl;
      ok = 0;
      }
      }
    */


    if (!ok)
      {
	(*testout) << "surfelements: " << endl;
	for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	    const Element2d & el = mesh.SurfaceElement(i);
	    (*testout) << setw(5) << i << ":" 
		       << setw(6) << el.GetIndex() 
		       << setw(6) << el.PNum(1) 
		       << setw(4) << el.PNum(2) 
		       << setw(4) << el.PNum(3)  << endl;
	  }
	(*testout) << "volelements: " << endl;
	for (ei = 0; ei < mesh.GetNE(); ei++)
	  {
	    const Element & el = mesh[ei];
	    (*testout) << setw(5) << i << ":" 
		       << setw(6) << el.GetIndex() 
		       << setw(6) << el[0] << setw(4) << el[1]
		       << setw(4) << el[2] << setw(4) << el[3] << endl;
	  }
      }


    return ok;
  }



  void RemoveProblem (Mesh & mesh, int domainnr)
  {
    int i, j, k;
  
    mesh.FindOpenElements(domainnr);
    int np = mesh.GetNP();

    BitArrayChar<PointIndex::BASE> ppoints(np);
  
    // int ndom = mesh.GetNDomains();

    PrintMessage (3, "Elements before Remove: ", mesh.GetNE());
    // for (k = 1; k <= ndom; k++)
    k = domainnr;
      {
	ppoints.Clear();
      
	for (i = 1; i <= mesh.GetNOpenElements(); i++)
	  {
	    const Element2d & sel = mesh.OpenElement(i);
	    if (sel.GetIndex() == k)
	      {
		for (j = 1; j <= sel.GetNP(); j++)
		  ppoints.Set (sel.PNum(j));
	      }
	  }

	for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	  {
	    const Element & el = mesh[ei];
	    if (el.GetIndex() == k)
	      {
		int todel = 0;
		for (j = 0; j < el.GetNP(); j++)
		  if (ppoints.Test (el[j]))
		    todel = 1;
	      
		if (el.GetNP() != 4)
		  todel = 0;
	      
		if (todel)
		  {
		    mesh[ei].Delete();
		    // ei--;
		  }
	      }
	  }
      }
  
    mesh.Compress();
    PrintMessage (3, "Elements after Remove: ", mesh.GetNE());
  }


}
