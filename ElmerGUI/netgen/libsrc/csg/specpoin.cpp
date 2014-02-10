#include <mystdlib.h>
#include <meshing.hpp>
#include <csg.hpp>


/*
  Special Point calculation uses the global Flags:

  relydegtest       when to rely on degeneration ?
  calccp            calculate points of intersection ?
  cpeps1            eps for degenerated poi
  calcep            calculate points of extreme coordinates ?
  epeps1            eps for degenerated edge
  epeps2            eps for axis parallel pec
  epspointdist      eps for distance of special points 
*/


#undef DEVELOP
//#define DEVELOP


namespace netgen
{
  ARRAY<Box<3> > boxes;


  void ProjectToEdge (const Surface * f1, const Surface * f2, Point<3> & hp);



  SpecialPoint :: SpecialPoint (const SpecialPoint & sp)
  {
    p = sp.p;
    v = sp.v;
    s1 = sp.s1;
    s2 = sp.s2;
    s1_orig = sp.s1_orig;
    s2_orig = sp.s2_orig;
    layer = sp.layer;
    unconditional = sp.unconditional;
  }
  
  SpecialPoint & SpecialPoint :: operator= (const SpecialPoint & sp)
  {
    p = sp.p;
    v = sp.v;
    s1 = sp.s1;
    s2 = sp.s2;
    s1_orig = sp.s1_orig;
    s2_orig = sp.s2_orig;
    layer = sp.layer;
    unconditional = sp.unconditional;
    return *this;
  }


  void SpecialPoint :: Print (ostream & str) const
  {
    str << "p = " << p << "   v = " << v 
	<< " s1/s2 = " << s1 << "/" << s2;
    str << " layer = " << layer
	<< " unconditional = " << unconditional
	<< endl;
  }


  static ARRAY<int> numprim_hist;

  SpecialPointCalculation :: SpecialPointCalculation ()
  {
    ideps = 1e-9;
  }

  void SpecialPointCalculation :: 
  CalcSpecialPoints (const CSGeometry & ageometry, 
		     ARRAY<MeshPoint> & apoints)
  {
    geometry = &ageometry;
    points = &apoints;

    size = geometry->MaxSize();
    (*testout) << "Find Special Points" << endl;
    (*testout) << "maxsize = " << size << endl;

    cpeps1 = 1e-6; 
    epeps1 = 1e-3; 
    epeps2 = 1e-6; 

    epspointdist2 = sqr (size * 1e-8); 
    relydegtest = size * 1e-4; 


    BoxSphere<3> box (Point<3> (-size, -size, -size),
		      Point<3> ( size,  size,  size));

    box.CalcDiamCenter();
    PrintMessage (3, "main-solids: ", geometry->GetNTopLevelObjects());

    numprim_hist.SetSize (geometry->GetNSurf()+1);
    numprim_hist = 0;

    for (int i = 0; i < geometry->GetNTopLevelObjects(); i++)
      {
	const TopLevelObject * tlo = geometry->GetTopLevelObject(i);

	(*testout) << "tlo " << i << ":" << endl
		   << *tlo->GetSolid() << endl;

	if (tlo->GetSolid())
	  {
	    ARRAY<Point<3> > hpts;
	    tlo->GetSolid()->CalcOnePrimitiveSpecialPoints (box, hpts);
            // if (hpts.Size())
            //  cout << "oneprimitivespecialpoints = " << hpts << endl;
	    for (int j = 0; j < hpts.Size(); j++)
	      AddPoint (hpts[j], tlo->GetLayer());
	  }

	CalcSpecialPointsRec (tlo->GetSolid(), tlo->GetLayer(),
			      box, 1, 1, 1);
      }
 
  
    geometry->DeleteIdentPoints();
    for (int i = 0; i < geometry->GetNIdentifications(); i++)
      {
	CloseSurfaceIdentification * ident =
	  dynamic_cast<CloseSurfaceIdentification * >(geometry->identifications[i]);
	
	if(!ident || !ident->IsSkewIdentification())
	  continue;

	for(int j=0; j<points->Size(); j++)
	  {
	    if(fabs(ident->GetSurface1().CalcFunctionValue((*points)[j])) < 1e-15)
	      {
		Point<3> auxpoint = (*points)[j];
		ident->GetSurface2().SkewProject(auxpoint,ident->GetDirection());
		geometry->AddIdentPoint(auxpoint);
		geometry->AddIdentPoint((*points)[j]);
		AddPoint (auxpoint,1);

#ifdef DEVELOP
		(*testout) << "added identpoint " << auxpoint << "; proj. of "
			   <<  (*points)[j] << endl;
#endif
		break;
	      }
	  }
      }
    

    // add user point:
    for (int i = 0; i < geometry->GetNUserPoints(); i++)
      AddPoint (geometry->GetUserPoint(i), 1);
	

    PrintMessage (3, "Found points ", apoints.Size());

    for (int i = 0; i < boxesinlevel.Size(); i++)
      (*testout) << "level " << i << " has " 
		 << boxesinlevel[i] << " boxes" << endl;
    (*testout) << "numprim_histogramm = " << endl << numprim_hist << endl;
  }
  


  void SpecialPointCalculation :: 
  CalcSpecialPointsRec (const Solid * sol, int layer,
			const BoxSphere<3> & box, 
			int level, bool calccp, bool calcep)
  {
    // boxes.Append (box);

#ifdef DEVELOP
    *testout << "lev " << level << ", box = " << box << endl;
    *testout << "calccp = " << calccp << ", calcep = " << calcep << endl;
    *testout << "locsol = " << *sol << endl;
#endif

    if (multithread.terminate)
      {
	*testout << "boxes = " << boxes << endl;
	*testout << "boxesinlevel = " << boxesinlevel << endl;
	throw NgException ("Meshing stopped");
      }


    if (!sol) return;

    if (level >= 100)
      {
	MyStr err =
	  MyStr("Problems in CalcSpecialPoints\nPoint: ") + MyStr (box.Center());
	throw NgException (err.c_str());
      }


    bool decision;
    bool possiblecrossp, possibleexp;  // possible cross or extremalpoint
    bool surecrossp = 0, sureexp = 0;          // sure ...
  
    static ARRAY<int> locsurf;  // attention: array is static

    static int cntbox = 0;
    cntbox++;

    if (level <= boxesinlevel.Size())
      boxesinlevel.Elem(level)++;
    else
      boxesinlevel.Append (1);

    /*
      numprim = sol -> NumPrimitives();
      sol -> GetSurfaceIndices (locsurf);
    */

    geometry -> GetIndependentSurfaceIndices (sol, box, locsurf);
    int numprim = locsurf.Size();

#ifdef DEVELOP
    (*testout) << "numprim = " << numprim << endl;
#endif

    numprim_hist[numprim]++;

    Point<3> p = box.Center();


    // explicit solution for planes only and at most one quadratic
    if (numprim <= 5)
      {
	int nplane = 0, nquad = 0, quadi = -1;
	const QuadraticSurface *qsurf = 0, *qsurfi;

	for (int i = 0; i < numprim; i++)
	  {
	    qsurfi = dynamic_cast<const QuadraticSurface*> 
	      (geometry->GetSurface(locsurf[i]));

	    if (qsurfi) nquad++;
	    if (dynamic_cast<const Plane*> (qsurfi))
	      nplane++;
	    else
	      {
		quadi = i;
		qsurf = qsurfi;
	      }
	  }

	/*
	if (nquad == numprim && nplane == numprim-2)
	  return;
	*/

#ifdef DEVELOP
	(*testout) << "nquad " << nquad << " nplane " << nplane << endl;
#endif

	if (nquad == numprim && nplane >= numprim-1)
	  {
	    ARRAY<Point<3> > pts;
	    ARRAY<int> surfids;

	    for (int k1 = 0; k1 < numprim - 2; k1++)
	      for (int k2 = k1 + 1; k2 < numprim - 1; k2++)
		for (int k3 = k2 + 1; k3 < numprim; k3++)
		  if (k1 != quadi && k2 != quadi && k3 != quadi)
		    {
		      ComputeCrossPoints (dynamic_cast<const Plane*> (geometry->GetSurface(locsurf[k1])),
					  dynamic_cast<const Plane*> (geometry->GetSurface(locsurf[k2])),
					  dynamic_cast<const Plane*> (geometry->GetSurface(locsurf[k3])),
					  pts);
		     
		      for (int j = 0; j < pts.Size(); j++)
			if (Dist (pts[j], box.Center()) < box.Diam()/2)
			  {
			    Solid * tansol;
			    sol -> TangentialSolid (pts[j], tansol, surfids, 1e-9*size);

			    if(!tansol)
			      continue;

			    bool ok1 = false, ok2 = false, ok3 = false;
			    int rep1 = geometry->GetSurfaceClassRepresentant(locsurf[k1]);
			    int rep2 = geometry->GetSurfaceClassRepresentant(locsurf[k2]);
			    int rep3 = geometry->GetSurfaceClassRepresentant(locsurf[k3]);
			    for(int jj=0; jj<surfids.Size(); jj++)
			      {
				int actrep = geometry->GetSurfaceClassRepresentant(surfids[jj]);
				if(actrep == rep1) ok1 = true;
				if(actrep == rep2) ok2 = true;
				if(actrep == rep3) ok3 = true;				  
			      }

			    
			    if (tansol && ok1 && ok2 && ok3)
			    // if (sol -> IsIn (pts[j], 1e-6*size) && !sol->IsStrictIn (pts[j], 1e-6*size))
			      {
				if (AddPoint (pts[j], layer))
				  (*testout) << "cross point found, 1: " << pts[j] << endl;
			      }  
			    delete tansol;
			  }
		    }


	    if (qsurf)
	      {
		for (int k1 = 0; k1 < numprim - 1; k1++)
		  for (int k2 = k1 + 1; k2 < numprim; k2++)
		    if (k1 != quadi && k2 != quadi)
		      {
			ComputeCrossPoints (dynamic_cast<const Plane*> (geometry->GetSurface(locsurf[k1])),
					    dynamic_cast<const Plane*> (geometry->GetSurface(locsurf[k2])),
					    qsurf, pts);
			//(*testout) << "checking pot. crosspoints: " << pts << endl;

			for (int j = 0; j < pts.Size(); j++)
			  if (Dist (pts[j], box.Center()) < box.Diam()/2)
			    {
			      Solid * tansol;
			      sol -> TangentialSolid (pts[j], tansol, surfids, 1e-9*size);

			      if(!tansol)
				continue;
			      			      
			      bool ok1 = false, ok2 = false, ok3 = true;//false;
			      int rep1 = geometry->GetSurfaceClassRepresentant(locsurf[k1]);
			      int rep2 = geometry->GetSurfaceClassRepresentant(locsurf[k2]);
			      //int rep3 = geometry->GetSurfaceClassRepresentant(quadi);
			      for(int jj=0; jj<surfids.Size(); jj++)
				{
				  int actrep = geometry->GetSurfaceClassRepresentant(surfids[jj]);
				  if(actrep == rep1) ok1 = true;
				  if(actrep == rep2) ok2 = true;
				  //if(actrep == rep3) ok3 = true;				  
				}


			      if (tansol && ok1 && ok2 && ok3)
				//if (sol -> IsIn (pts[j], 1e-6*size) && !sol->IsStrictIn (pts[j], 1e-6*size) )
				{
				  if (AddPoint (pts[j], layer))
				    (*testout) << "cross point found, 2: " << pts[j] << endl;
				}  
			      delete tansol;
			    }
		      }


		for (int k1 = 0; k1 < numprim; k1++)
		  if (k1 != quadi)
		    {
		      ComputeExtremalPoints (dynamic_cast<const Plane*> (geometry->GetSurface(locsurf[k1])),
					     qsurf, pts);
		      
		      for (int j = 0; j < pts.Size(); j++)
			if (Dist (pts[j], box.Center()) < box.Diam()/2)
			  {
			    Solid * tansol;
			    sol -> TangentialSolid (pts[j], tansol, surfids, 1e-9*size);
			    if (tansol)
			      // sol -> IsIn (pts[j], 1e-6*size) && !sol->IsStrictIn (pts[j], 1e-6*size) )
			      {
				if (AddPoint (pts[j], layer))
				  (*testout) << "extremal point found, 1: " << pts[j] << endl;
			      }  
			    delete tansol;
			  }
		    }
	      }
	    
	    return;
	  }
      }


    
    possiblecrossp = (numprim >= 3) && calccp;
    surecrossp = 0;

    if (possiblecrossp && (locsurf.Size() <= 5 || level > 50))
      {
	decision = 1;
	surecrossp = 0;

	for (int k1 = 1; k1 <= locsurf.Size() - 2; k1++)
	  for (int k2 = k1 + 1; k2 <= locsurf.Size() - 1; k2++)
	    for (int k3 = k2 + 1; k3 <= locsurf.Size(); k3++)
	      {
		int nc, deg;
		nc = CrossPointNewtonConvergence 
		  (geometry->GetSurface(locsurf.Get(k1)), 
		   geometry->GetSurface(locsurf.Get(k2)), 
		   geometry->GetSurface(locsurf.Get(k3)), box );
	      
		deg = CrossPointDegenerated 
		  (geometry->GetSurface(locsurf.Get(k1)), 
		   geometry->GetSurface(locsurf.Get(k2)), 
		   geometry->GetSurface(locsurf.Get(k3)), box );
	      
		if (!nc && !deg) decision = 0;
		if (nc) surecrossp = 1;
	      }

#ifdef DEVELOP
        (*testout) << "dec = " << decision << ", surcp = " << surecrossp << endl;
#endif

	if (decision && surecrossp)
	  {
	    for (int k1 = 1; k1 <= locsurf.Size() - 2; k1++)
	      for (int k2 = k1 + 1; k2 <= locsurf.Size() - 1; k2++)
		for (int k3 = k2 + 1; k3 <= locsurf.Size(); k3++)
		  {
		    if (CrossPointNewtonConvergence 
			(geometry->GetSurface(locsurf.Get(k1)), 
			 geometry->GetSurface(locsurf.Get(k2)), 
			 geometry->GetSurface(locsurf.Get(k3)), box ) )
		      {
                        
			Point<3> pp = p;
			CrossPointNewton 
			  (geometry->GetSurface(locsurf.Get(k1)), 
			   geometry->GetSurface(locsurf.Get(k2)), 
			   geometry->GetSurface(locsurf.Get(k3)), pp);
              
			BoxSphere<3> hbox (pp, pp);
			hbox.Increase (1e-8*size);

			if (pp(0) > box.PMin()(0) - 1e-5*size && 
			    pp(0) < box.PMax()(0) + 1e-5*size &&
			    pp(1) > box.PMin()(1) - 1e-5*size && 
			    pp(1) < box.PMax()(1) + 1e-5*size &&
			    pp(2) > box.PMin()(2) - 1e-5*size && 
			    pp(2) < box.PMax()(2) + 1e-5*size &&
			    sol -> IsIn (pp, 1e-6*size) && !sol->IsStrictIn (pp, 1e-6*size) &&
			    !CrossPointDegenerated
			    (geometry->GetSurface(locsurf.Get(k1)), 
			     geometry->GetSurface(locsurf.Get(k2)), 
			     geometry->GetSurface(locsurf.Get(k3)), hbox ))

			  { 
			    //                AddCrossPoint (locsurf, sol, p);
			    BoxSphere<3> boxp (pp, pp);
			    boxp.Increase (1e-3*size);
			    boxp.CalcDiamCenter();
			    ARRAY<int> locsurf2;

			    geometry -> GetIndependentSurfaceIndices (sol, boxp, locsurf2);
			  
			    bool found1 = false, found2 = false, found3 = false;
			    for (int i = 0; i < locsurf2.Size(); i++)
			      {
				if (locsurf2[i] == locsurf.Get(k1)) found1 = true;
				if (locsurf2[i] == locsurf.Get(k2)) found2 = true;
				if (locsurf2[i] == locsurf.Get(k3)) found3 = true;
			      }

			    if (found1 && found2 && found3)
			      if (AddPoint (pp, layer))
				{
				  (*testout) << "Crosspoint found: " << pp 
					     << " diam = " << box.Diam()
					     << ",  surfs: " 
					     << locsurf.Get(k1) << "," 
					     << locsurf.Get(k2) << "," 
					     << locsurf.Get(k3) << endl;
				}
			  }
		      }
		  }
	  }
      
	if (decision)
	  possiblecrossp = 0;
      }


    possibleexp = (numprim >= 2) && calcep;

    // (*testout) << "l = " << level << "locsize = " << locsurf.Size() << " possexp = " << possibleexp << "\n";
    if (possibleexp && (numprim <= 5 || level >= 50))
      {
	decision = 1;
	sureexp = 0;

	/*
	(*testout) << "extremal surfs = ";
	for (int k5 = 0; k5 < locsurf.Size(); k5++)
	  (*testout) << typeid(*geometry->GetSurface(locsurf[k5])).name() << " ";
	(*testout) << "\n";
	*/

	for (int k1 = 0; k1 < locsurf.Size() - 1; k1++)
	  for (int k2 = k1+1; k2 < locsurf.Size(); k2++)
	    {
	      const Surface * surf1 = geometry->GetSurface(locsurf[k1]);
	      const Surface * surf2 = geometry->GetSurface(locsurf[k2]);
	      /*
	      (*testout) << "edgecheck, types = " << typeid(*surf1).name() << ", " << typeid(*surf2).name()
			 << "edge-newton-conv = " << EdgeNewtonConvergence (surf1, surf2, p)
			 << "edge-deg = " << EdgeDegenerated (surf1, surf2, box)
			 << "\n";
	      */

	      if (EdgeNewtonConvergence (surf1, surf2, p) ) 
		sureexp = 1;
	      else
		{
		  if (!EdgeDegenerated (surf1, surf2, box)) 
		    decision = 0;
		}
	    }

	// (*testout) << "l = " << level << " dec/sureexp = " << decision << sureexp << endl;

	if (decision && sureexp)
	  {
	    for (int k1 = 0; k1 < locsurf.Size() - 1; k1++)
	      for (int k2 = k1+1; k2 < locsurf.Size(); k2++)
		{
		  const Surface * surf1 = geometry->GetSurface(locsurf[k1]);
		  const Surface * surf2 = geometry->GetSurface(locsurf[k2]);

		  if (EdgeNewtonConvergence (surf1, surf2, p))
		    {
		      EdgeNewton (surf1, surf2, p);
		    
		      Point<3> pp;
		      if (IsEdgeExtremalPoint (surf1, surf2, p, pp, box.Diam()/2))
			{
			  (*testout) << "extremalpoint (nearly) found:" << pp << endl;

			  if (Dist (pp, box.Center()) < box.Diam()/2 &&
			      sol -> IsIn (pp, 1e-6*size) && !sol->IsStrictIn (pp, 1e-6*size) )
			    {
			      if (AddPoint (pp, layer))
				(*testout) << "Extremal point found: " << pp << endl;//"(eps="<<1e-9*size<<")"<< endl;
			    }  
			}            
		    }
		}
	  }
	if (decision)
	  possibleexp = 0;
      }
 

    // (*testout) << "l = " << level << " poss cp/ep sure exp = " << possiblecrossp << " " << possibleexp << " " << sureexp << "\n";
    if (possiblecrossp || possibleexp)
      {
	BoxSphere<3> sbox;
	for (int i = 0; i < 8; i++)
	  {
	    box.GetSubBox (i, sbox);
	    sbox.Increase (1e-4 * sbox.Diam());

	    Solid * redsol = sol -> GetReducedSolid (sbox);

	    if (redsol)
	      {
		CalcSpecialPointsRec (redsol, layer, sbox, level+1, calccp, calcep);
		delete redsol;
	      }
	  }
      }
  }





  /******* Tests for Point of intersection **********************/



  bool SpecialPointCalculation :: 
  CrossPointNewtonConvergence (const Surface * f1, 
			       const Surface * f2, 
			       const Surface * f3,
			       const BoxSphere<3> & box)
  {
    Vec<3> grad, rs, x;
    Mat<3> jacobi, inv;
    Point<3> p = box.Center();

    f1->CalcGradient (p, grad);
    jacobi(0,0) = grad(0);
    jacobi(0,1) = grad(1);
    jacobi(0,2) = grad(2);

    f2->CalcGradient (p, grad);
    jacobi(1,0) = grad(0);
    jacobi(1,1) = grad(1);
    jacobi(1,2) = grad(2);

    f3->CalcGradient (p, grad);
    jacobi(2,0) = grad(0);
    jacobi(2,1) = grad(1);
    jacobi(2,2) = grad(2);

    if (fabs (Det (jacobi)) > 1e-8)
      {
	double gamma = f1 -> HesseNorm() + f2 -> HesseNorm() + f3 -> HesseNorm();
	if (gamma == 0.0) return 1;

	CalcInverse (jacobi, inv);

	rs(0) = f1->CalcFunctionValue (p);
	rs(1) = f2->CalcFunctionValue (p);
	rs(2) = f3->CalcFunctionValue (p);

	x = inv * rs;

	double beta = 0;
	for (int i = 0; i < 3; i++)
	  {
	    double sum = 0;
	    for (int j = 0; j < 3; j++)
	      sum += fabs (inv(i,j));
	    if (sum > beta)  beta = sum;
	  }
	double eta = Abs (x);


#ifdef DEVELOP
        *testout << "check Newton: " << "beta = " << beta << ", gamma = " << gamma << ", eta = " << eta << endl;
        double rad = 1.0 / (beta * gamma);
        *testout << "rad = " << rad << endl;
#endif
        
	return (beta * gamma * eta < 0.1) && (2 > box.Diam()*beta*gamma);
      }
    return 0;

  }




  bool SpecialPointCalculation :: 
  CrossPointDegenerated (const Surface * f1,
			 const Surface * f2, 
			 const Surface * f3, 
			 const BoxSphere<3> & box) const
  {
    Mat<3> mat;
    Vec<3> g1, g2, g3;
    double normprod;

    if (box.Diam() > relydegtest) return 0;

    f1->CalcGradient (box.Center(), g1);
    normprod = Abs2 (g1);

    f2->CalcGradient (box.Center(), g2);
    normprod *= Abs2 (g2);
 
    f3->CalcGradient (box.Center(), g3);
    normprod *= Abs2 (g3);

    for (int i = 0; i < 3; i++)
      {
	mat(i,0) = g1(i);
	mat(i,1) = g2(i);
	mat(i,2) = g3(i);
      }

    return sqr (Det (mat)) < sqr(cpeps1) * normprod;
  }
 




  void SpecialPointCalculation :: CrossPointNewton (const Surface * f1, 
						    const Surface * f2, 
						    const Surface * f3, Point<3> & p)
  {
    Vec<3> g1, g2, g3;
    Vec<3> rs, sol;
    Mat<3> mat;

    int i = 10;
    while (i > 0)
      {
	i--;
	rs(0) = f1->CalcFunctionValue (p);
	rs(1) = f2->CalcFunctionValue (p);
	rs(2) = f3->CalcFunctionValue (p);

	f1->CalcGradient (p, g1);
	f2->CalcGradient (p, g2);
	f3->CalcGradient (p, g3);

	for (int j = 0; j < 3; j++)
	  {
	    mat(0, j) = g1(j);
	    mat(1, j) = g2(j);
	    mat(2, j) = g3(j);
	  }
	mat.Solve (rs, sol);
	if (sol.Length2() < 1e-24 && i > 1) i = 1;

#ifdef DEVELOP
        *testout << "CrossPointNewton, err = " << sol.Length2() << endl;
#endif
	p -= sol;
      }
  }




  /******* Tests for Point on edges **********************/




  bool SpecialPointCalculation :: 
  EdgeNewtonConvergence (const Surface * f1, const Surface * f2, 
			 const Point<3> & p)
  {
    Vec<3> g1, g2, sol;
    Vec<2> vrs;
    Mat<2,3> mat;
    Mat<3,2> inv;

    f1->CalcGradient (p, g1);
    f2->CalcGradient (p, g2);

    if ( sqr(g1 * g2) < (1 - 1e-8) * Abs2 (g1) * Abs2 (g2))
      {
	double gamma = f1 -> HesseNorm() + f2 -> HesseNorm();
	if (gamma < 1e-32) return 1;
	gamma = sqr (gamma);
      
	for (int i = 0; i < 3; i++)
	  {
	    mat(0,i) = g1(i);
	    mat(1,i) = g2(i);
	  }

	CalcInverse (mat, inv);

	vrs(0) = f1->CalcFunctionValue (p);
	vrs(1) = f2->CalcFunctionValue (p);
	sol = inv * vrs;

	double beta = 0;
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 2; j++)
	    beta += inv(i,j) * inv(i,j);
	// beta = sqrt (beta);

	double eta = Abs2 (sol);

	// alpha = beta * gamma * eta;
	return (beta * gamma * eta < 0.01);
      }
    return 0;
  }




  bool SpecialPointCalculation :: 
  EdgeDegenerated (const Surface * f1,
		   const Surface * f2, 
		   const BoxSphere<3> & box) const
  {
    // perform newton steps. normals parallel ?
    // if not decideable: return 0 
  
    Point<3> p = box.Center();
    Vec<3> g1, g2, sol;
    Vec<2> vrs;
    Mat<2,3> mat;

    int i = 20;
    while (i > 0)
      {
	if (Dist2 (p, box.Center()) > sqr(box.Diam()))
	  return 0;

	i--;
	vrs(0) = f1->CalcFunctionValue (p);
	vrs(1) = f2->CalcFunctionValue (p);

	f1->CalcGradient (p, g1);
	f2->CalcGradient (p, g2);

	if ( sqr (g1 * g2) > (1 - 1e-10) * Abs2 (g1) * Abs2 (g2))
	  return 1;

	for (int j = 0; j < 3; j++)
	  {
	    mat(0,j) = g1(j);
	    mat(1,j) = g2(j);
	  }
	mat.Solve (vrs, sol);

	if (Abs2 (sol) < 1e-24 && i > 1) i = 1;
	p -= sol;
      }

    return 0;
  }






  void SpecialPointCalculation :: EdgeNewton (const Surface * f1, 
					      const Surface * f2, Point<3> & p)
  {
    Vec<3> g1, g2, sol;
    Vec<2> vrs;
    Mat<2,3> mat;

    int i = 10;
    while (i > 0)
      {
	i--;
	vrs(0) = f1->CalcFunctionValue (p);
	vrs(1) = f2->CalcFunctionValue (p);

	f1->CalcGradient (p, g1);
	f2->CalcGradient (p, g2);

	//(*testout) << "p " << p << " f1 " << vrs(0) << " f2 " << vrs(1) << " g1 " << g1 << " g2 " << g2 << endl;

	for (int j = 0; j < 3; j++)
	  {
	    mat(0,j) = g1(j);
	    mat(1,j) = g2(j);
	  }
	mat.Solve (vrs, sol);
	
	if (Abs2 (sol) < 1e-24 && i > 1) i = 1;
	p -= sol;
      }
  }



  bool SpecialPointCalculation :: 
  IsEdgeExtremalPoint (const Surface * f1, const Surface * f2, 
		       const Point<3> & p, Point<3> & pp, double rad)
  {
    Vec<3> g1, g2, t, t1, t2;

    f1->CalcGradient (p, g1);
    f2->CalcGradient (p, g2);
  
    t = Cross (g1, g2);
    t.Normalize();

    Point<3> p1 = p + rad * t;
    Point<3> p2 = p - rad * t;

    EdgeNewton (f1, f2, p1);
    EdgeNewton (f1, f2, p2);
  
    f1->CalcGradient (p1, g1);
    f2->CalcGradient (p1, g2);
    t1 = Cross (g1, g2);
    t1.Normalize();

    f1->CalcGradient (p2, g1);
    f2->CalcGradient (p2, g2);
    t2 = Cross (g1, g2);
    t2.Normalize();

    double val = 1e-8 * rad * rad;
    for (int j = 0; j < 3; j++)
      if ( (t1(j) * t2(j) < -val) )
	{
	  pp = p;
	  ExtremalPointNewton (f1, f2, j+1, pp);
	  return 1;
	}

    return 0;
  }









  /********** Tests of Points of extremal coordinates  ****************/


  void SpecialPointCalculation :: ExtremalPointNewton (const Surface * f1, 
						       const Surface * f2, 
						       int dir, Point<3> & p)
  {
    Vec<3> g1, g2, v, curv;
    Vec<3> rs, x, y1, y2, y;
    Mat<3> h1, h2;
    Mat<3> jacobi;

    int i = 50;
    while (i > 0)
      {
	i--;
	rs(0) = f1->CalcFunctionValue (p);
	rs(1) = f2->CalcFunctionValue (p);

	f1 -> CalcGradient (p, g1);
	f2 -> CalcGradient (p, g2);

	f1 -> CalcHesse (p, h1);
	f2 -> CalcHesse (p, h2);

	v = Cross (g1, g2);

	rs(2) = v(dir-1);

	jacobi(0,0) = g1(0);
	jacobi(0,1) = g1(1);
	jacobi(0,2) = g1(2);

	jacobi(1,0) = g2(0);
	jacobi(1,1) = g2(1);
	jacobi(1,2) = g2(2);


	switch (dir)
	  {
	  case 1:
	    {
	      y1(0) = 0;
	      y1(1) = g2(2);
	      y1(2) = -g2(1);
	      y2(0) = 0;
	      y2(1) = -g1(2);
	      y2(2) = g1(1);
	      break;
	    }
	  case 2:
	    {
	      y1(0) = -g2(2);
	      y1(1) = 0;
	      y1(2) = g2(0);
	      y2(0) = g1(2);
	      y2(1) = 0;
	      y2(2) = -g1(0);
	      break;
	    }
	  case 3:
	    {
	      y1(0) = g2(1);
	      y1(1) = -g2(0);
	      y1(2) = 0;
	      y2(0) = -g1(1);
	      y2(1) = g1(0);
	      y2(2) = 0;
	      break;
	    }
	  }

	y = h1 * y1 + h2 * y2;

	jacobi(2,0) = y(0);
	jacobi(2,1) = y(1);
	jacobi(2,2) = y(2);

	/*
	(*testout) << "p " << p << " f1 " << rs(0) << " f2 " << rs(1) << endl
		   << " jacobi " << jacobi << endl
		   << " rhs " << rs << endl;
	*/	

	jacobi.Solve (rs, x);

	if (Abs2 (x) < 1e-24 && i > 1)
	  {
	    i = 1;
	  }

	
	double minval(Abs2(rs)),minfac(1);
	double startval(minval);
	for(double fac = 1; fac > 1e-7; fac *= 0.6)
	  {
	    Point<3> testpoint = p-fac*x;

	    rs(0) = f1->CalcFunctionValue (testpoint);
	    rs(1) = f2->CalcFunctionValue (testpoint);

	    f1 -> CalcGradient (testpoint, g1);
	    f2 -> CalcGradient (testpoint, g2);

	    v = Cross (g1, g2);

	    rs(2) = v(dir-1);

	    double val = Abs2(rs);

	    if(val < minval)
	      {
		minfac = fac;
		if(val < 0.5 * startval)
		  break;
		minval = val;
	      }

	  }
	p -= minfac*x;
	

	//p -= x;
      }


    if (Abs2 (x) > 1e-20)
      {
	(*testout) << "Error: extremum Newton not convergent" << endl;
	(*testout) << "dir = " << dir << endl;
	(*testout) << "p = " << p << endl;
	(*testout) << "x = " << x << endl;
      }
  }

  void SpecialPointCalculation :: 
  ComputeCrossPoints (const Plane * plane1, 
		      const Plane * plane2, 
		      const Plane * plane3, 
		      ARRAY<Point<3> > & pts)
  {
    Mat<3> mat;
    Vec<3> rhs, sol;
    Point<3> p0(0,0,0);

    pts.SetSize (0);
    for (int i = 0; i < 3; i++)
      {
	const Plane * pi(NULL);
	switch (i)
	  {
	  case 0: pi = plane1; break;
	  case 1: pi = plane2; break;
	  case 2: pi = plane3; break;
	  }

	double val;
	Vec<3> hvec;
	val = pi -> CalcFunctionValue(p0);
	pi -> CalcGradient (p0, hvec);

	for (int j = 0; j < 3; j++)
	  mat(i,j) = hvec(j);
	rhs(i) = -val;
      }

    if (fabs (Det (mat)) > 1e-8)
      {
	mat.Solve (rhs, sol);
	pts.Append (Point<3> (sol));
      }
  }





  void SpecialPointCalculation :: 
  ComputeCrossPoints (const Plane * plane1, 
		      const Plane * plane2, 
		      const QuadraticSurface * quadric, 
		      ARRAY<Point<3> > & pts)
  {
    Mat<2,3> mat;
    Mat<3,2> inv;
    Vec<2> rhs;
    Vec<3> sol, t;
    Point<3> p0(0,0,0);

    pts.SetSize (0);
    for (int i = 0; i < 2; i++)
      {
	const Plane * pi(NULL);
	switch (i)
	  {
	  case 0: pi = plane1; break;
	  case 1: pi = plane2; break;
	  }

	double val;
	Vec<3> hvec;
	val = pi -> CalcFunctionValue(p0);
	pi -> CalcGradient (p0, hvec);

	for (int j = 0; j < 3; j++)
	  mat(i,j) = hvec(j);
	rhs(i) = -val;
      }
    CalcInverse (mat, inv);
    sol = inv * rhs;
    t = Cross (mat.Row(0), mat.Row(1));

    if (t.Length() > 1e-8)
      {
	Point<3> p (sol);
	// quadratic on  p + s t = 0
	double quad_a;
	Vec<3> quad_b;
	Mat<3> quad_c;
	
	quad_a = quadric -> CalcFunctionValue(p);
	quadric -> CalcGradient (p, quad_b);
	quadric -> CalcHesse (p, quad_c);
	
	double a, b, c;
	a = quad_a;
	b = quad_b * t;
	c = 0.5 * t * (quad_c * t);

	// a  + s b + s^2 c = 0;
	double disc = b*b-4*a*c;
	if (disc > 1e-10 * fabs (b))
	  {
	    disc = sqrt (disc);
	    double s1 = (-b-disc) / (2*c);
	    double s2 = (-b+disc) / (2*c);

	    pts.Append (p + s1 * t);
	    pts.Append (p + s2 * t);
	  }
      }
  }












  void SpecialPointCalculation :: 
  ComputeExtremalPoints (const Plane * plane, 
			 const QuadraticSurface * quadric, 
			 ARRAY<Point<3> > & pts)
  {
    // 3 equations:
    // surf1 = 0  <===> plane_a + plane_b x = 0;
    // surf2 = 0  <===> quad_a + quad_b x + x^T quad_c x = 0
    // (grad 1 x grad 2)(i) = 0  <====> (grad 1 x e_i) . grad_2 = 0

    pts.SetSize (0);

    Point<3> p0(0,0,0);
    double plane_a, quad_a;
    Vec<3> plane_b, quad_b, ei;
    Mat<3> quad_c;

    plane_a = plane -> CalcFunctionValue(p0);
    plane -> CalcGradient (p0, plane_b);

    quad_a = quadric -> CalcFunctionValue(p0);
    quadric -> CalcGradient (p0, quad_b);
    quadric -> CalcHesse (p0, quad_c);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	quad_c(i,j) *= 0.5;

    for (int dir = 0; dir <= 2; dir++)
      {
	ei = 0.0; ei(dir) = 1;
	Vec<3> v1 = Cross (plane_b, ei);
	
	// grad_2 . v1 ... linear:
	double g2v1_c = v1 * quad_b;
	Vec<3> g2v1_l = 2.0 * (quad_c * v1);

	// find line of two linear equations:
	
	Vec<2> rhs;
	Vec<3> sol;
	Mat<2,3> mat;

	for (int j = 0; j < 3; j++)
	  {
	    mat(0,j) = plane_b(j);
	    mat(1,j) = g2v1_l(j);
	  }
	rhs(0) = -plane_a;
	rhs(1) = -g2v1_c;

	Vec<3> t = Cross (plane_b, g2v1_l);
	if (Abs2(t) > 0)
	  {
	    mat.Solve (rhs, sol);
	    
	    // solve quadratic equation along line  sol + alpha t ....
	    double a = quad_a + quad_b * sol + sol * (quad_c * sol);
	    double b = quad_b * t + 2 * (sol * (quad_c * t));
	    double c = t * (quad_c * t);

	    // solve a + b alpha + c alpha^2:

	    if (fabs (c) > 1e-32)
	      {
		double disc = sqr (0.5*b/c) - a/c;
		if (disc > 0)
		  {
		    disc = sqrt (disc);
		    double alpha1 = -0.5*b/c + disc;
		    double alpha2 = -0.5*b/c - disc;

		    pts.Append (Point<3> (sol+alpha1*t));
		    pts.Append (Point<3> (sol+alpha2*t));
		    /*
		    cout << "sol1 = " << sol + alpha1 * t
			 << ", sol2 = " << sol + alpha2 * t << endl;
		    */
		  }
	      }
	  }
      }
  }






  /*
    bool SpecialPointCalculation :: ExtremalPointPossible (const Surface * f1, 
    const Surface * f2, 
    int dir, 
    const BoxSphere<3> & box)
    {
    double hn1, hn2, gn1, gn2;
    Point<3> p;
    Vec<3> g1, g2, v;
    double f3;
    double r = box.Diam()/2;

    p = box.Center();

    f1 -> CalcGradient (p, g1);
    f2 -> CalcGradient (p, g2);

    gn1 = g1.Length();
    gn2 = g2.Length();

    hn1 = f1 -> HesseNorm ();
    hn2 = f2 -> HesseNorm ();

    v = Cross (g1, g2);
    f3 = fabs (v(dir-1));

    //  (*testout) << "f3 = " << f3 << "  r = " << r 
    //             << "normbound = " 
    //             << (hn1 * (gn2 + r * hn2) + hn2 * (gn1 + r * hn1)) << endl;
 
    return (f3 <= 3 * r * (hn1 * (gn2 + r * hn2) + hn2 * (gn1 + r * hn1)));
    }



    bool SpecialPointCalculation :: 
    ExtremalPointNewtonConvergence (const Surface * f1, const Surface * f2, 
    int dir, 
    const BoxSphere<3> & box)
    {
    return box.Diam() < 1e-8;
    }


    bool SpecialPointCalculation :: 
    ExtremalPointDegenerated (const Surface * f1, const Surface * f2, 
    int dir, const BoxSphere<3> & box)
    {
    double gn1, gn2;
    Point<3> p;
    Vec<3> g1, g2, v;
    double maxderiv;
    double minv;
    Vec<3> curv, t;
    Vec<2> rs, x;
    Mat<3> h1, h2;
    Mat<2> a, inv;
    double leftside;

    if (box.Diam() > relydegtest) return 0;

    p = box.Center();

    f1 -> CalcGradient (p, g1);
    f2 -> CalcGradient (p, g2);
    gn1 = g1.Length();
    gn2 = g2.Length();

    v = Cross (g1, g2);
    if (Abs (v) < epeps1 * gn1 * gn2) return 1;       // irregular edge

    f1 -> CalcHesse (p, h1);
    f2 -> CalcHesse (p, h2);

    //  hn1 = f1 -> HesseNorm ();
    //  hn2 = f2 -> HesseNorm ();

    t = v;
    a(0, 0) = g1 * g1;
    a(0, 1) = 
    a(1, 0) = g1 * g2;
    a(1, 1) = g2 * g2;
  
    rs(0) = g1(dir-1);
    rs(1) = g2(dir-1);

    a.Solve (rs, x);

    //  (*testout) << "g1 = " << g1 << " g2 = " << g2 << endl;
    //  (*testout) << "lam = " << x << endl;
    //  (*testout) << "h2 = " << h2 << endl;

    leftside = fabs (x(0) * ( t * (h1 * t)) + 
    x(1) * ( t * (h2 * t)));

    //  (*testout) << "leftside = " << leftside << endl;

    if (leftside < epeps2 * Abs2 (v)) return 1;  

    return 0;
    }
  */
 

  bool SpecialPointCalculation :: AddPoint (const Point<3> & p, int layer)
  {
    for (int i = 0; i < points->Size(); i++)
      if (Dist2 ( (*points)[i], p) < epspointdist2 &&
	  (*points)[i].GetLayer() == layer)
	return false;

    points->Append (MeshPoint(p, layer));
    PrintMessageCR (3, "Found points ", points->Size());
    return true;
  }







  void SpecialPointCalculation :: 
  AnalyzeSpecialPoints (const CSGeometry & ageometry,
			ARRAY<MeshPoint> & apoints, 
			ARRAY<SpecialPoint> & specpoints)
  {
    ARRAY<int> surfind, rep_surfind, surfind2, rep_surfind2, surfind3;

    ARRAY<Vec<3> > normalvecs;
    Vec<3> nsurf;

    ARRAY<int> specpoint2point;
    specpoints.SetSize (0);

    geometry = &ageometry;

    double geomsize = ageometry.MaxSize();
 
    (*testout) << "AnalyzeSpecialPoints\n";

    if (!apoints.Size()) return;



#ifdef VERTSORT
    for (int i = 0; i < apoints.Size(); i++)
      for (int j = 0; j < apoints.Size()-1; j++)
        if (apoints[j](2) > apoints[j+1](2))
          swap (apoints[j], apoints[j+1]);
#endif







    Box<3> bbox (apoints[0], apoints[0]);
    for (int i = 1; i < apoints.Size(); i++)
      bbox.Add (apoints[i]);
    bbox.Increase (0.1 * bbox.Diam());

    //testout->precision(20);
    (*testout) << "bbox = " << bbox << endl;
    (*testout) << "points = " << apoints << endl;

    Point3dTree searchtree (bbox.PMin(), bbox.PMax());
    ARRAY<int> locsearch;

    for (int si = 0; si < ageometry.GetNTopLevelObjects(); si++)
      {
	const TopLevelObject * tlo = ageometry.GetTopLevelObject(si);

	const Solid * sol = tlo->GetSolid();
	const Surface * surf = tlo->GetSurface();


	for (int i = 0; i < apoints.Size(); i++)
	  {
	    Point<3> p = apoints[i];
	    
#ifdef DEVELOP
	    *testout << "                               test point " << p << endl;
#endif	    

	    if (tlo->GetLayer() != apoints[i].GetLayer())
	      continue;
	    
	    Solid * locsol;
	    sol -> TangentialSolid (p, locsol, surfind, ideps*geomsize);



	    rep_surfind.SetSize (surfind.Size());
	    int num_indep_surfs = 0;
	    
	    for (int j = 0; j < surfind.Size(); j++)
	      {
		rep_surfind[j] = ageometry.GetSurfaceClassRepresentant (surfind[j]);
		bool found = false;
		for (int k = 0; !found && k < j; k++)
		  found = (rep_surfind[k] == rep_surfind[j]);
		if(!found)
		  num_indep_surfs++;
	      }
	    

#ifdef DEVELOP
	    *testout << "surfs = " << surfind << endl;
	    *testout << "rep_surfs = " << rep_surfind << endl;
#endif

	    if (!locsol) continue;
	  
	    // get all surface indices, 
	    if (surf)
	      {
		// locsol -> GetSurfaceIndices (surfind);
		bool hassurf = 0;
		for (int m = 0; m < surfind.Size(); m++)
		  if (ageometry.GetSurface(surfind[m]) == surf)
		    hassurf = 1;

		if (!hassurf)
		  continue;

		nsurf = surf->GetNormalVector (p);
	      }

	    /*
	    // get independent surfaces of tangential solid
	    BoxSphere<3> box(p,p);
	    box.Increase (1e-6*geomsize);
	    box.CalcDiamCenter();
	    ageometry.GetIndependentSurfaceIndices (locsol, box, surfind);
	    */

	    // ageometry.GetIndependentSurfaceIndices (surfind);


	    normalvecs.SetSize(surfind.Size());
	    for (int j = 0; j < surfind.Size(); j++)
	      normalvecs[j] = 
		ageometry.GetSurface(surfind[j]) -> GetNormalVector(apoints[i]);


	    for (int j = 0; j < normalvecs.Size(); j++)
	      for (int k = 0; k < normalvecs.Size(); k++)
		{
		  if (rep_surfind[j] == rep_surfind[k]) continue;
		  //if (j == k) continue;

		  Vec<3> t;

		  if (dynamic_cast<const Polyhedra*> (ageometry.surf2prim[surfind[j]]) && 
		      ageometry.surf2prim[surfind[j]] == 
		      ageometry.surf2prim[surfind[k]])
		    {
		      t = ageometry.surf2prim[surfind[j]] -> 
			SpecialPointTangentialVector (p, surfind[j], surfind[k]);
		    }
		  else
		    {
		      t = Cross (normalvecs[j], normalvecs[k]);
		    }


		  if (Abs2 (t) < 1e-8)
		    continue;

#ifdef DEVELOP
		  *testout << "           tangential vector " << t << endl;
#endif

		  t.Normalize();

		  
		  // try tangential direction t
		  if (surf && fabs (nsurf * t) > 1e-6)
		    continue;

		
#ifdef DEVELOP
		  *testout << "           j " << j << " k " << k << endl;
#endif  

		  if (!surf)
		    {
		      // compute second order approximation
		      // c(s) = p + s t + s*s/2 t2
		      Vec<3> gradj, gradk;
		      Mat<3> hessej, hessek;
		      ageometry.GetSurface (surfind[j]) -> CalcGradient (p, gradj);
		      ageometry.GetSurface (surfind[k]) -> CalcGradient (p, gradk);
		      ageometry.GetSurface (surfind[j]) -> CalcHesse (p, hessej);
		      ageometry.GetSurface (surfind[k]) -> CalcHesse (p, hessek);
		      
		      Vec<2> rhs;
		      Vec<3> t2;
		      Mat<2,3> mat;
		      Mat<3,2> inv;
		      for (int l = 0; l < 3; l++)
			{
			  mat(0,l) = gradj(l);
			  mat(1,l) = gradk(l);
			}
		      rhs(0) = -t * (hessej * t);
		      rhs(1) = -t * (hessek * t);

		      CalcInverse (mat, inv);
		      t2 = inv * rhs;

		      
		      /*
		      ageometry.GetIndependentSurfaceIndices 
			(locsol, p, t, surfind2);
		      */

		      Solid * locsol2;
		      locsol -> TangentialSolid3 (p, t, t2, locsol2, surfind2, ideps*geomsize); 
		      if (!locsol2) continue;
		      
		      // locsol2 -> GetTangentialSurfaceIndices3 (p, t, t2, surfind2, 1e-9*geomsize);

		      rep_surfind2.SetSize (surfind2.Size());
		      for (int j2 = 0; j2 < surfind2.Size(); j2++)
			rep_surfind2[j2] = ageometry.GetSurfaceClassRepresentant (surfind2[j2]);

#ifdef DEVELOP
		      (*testout) << "surfind2 = " << endl << surfind2 << endl;
#endif
		      ARRAY<int> surfind2_aux(surfind2);
		      ageometry.GetIndependentSurfaceIndices (surfind2_aux);
#ifdef DEVELOP
		      (*testout) << "surfind2,rep = " << endl << surfind2_aux << endl;
#endif

		      bool ok = true;

		      // intersecting surfaces must be in second order tangential solid
		      /*
		      if (!surfind2.Contains(surfind[j]) ||
			  !surfind2.Contains(surfind[k]))
			ok = false;
		      */
		      if (!surfind2_aux.Contains(rep_surfind[j]) ||
			  !surfind2_aux.Contains(rep_surfind[k]))
			ok = false;

#ifdef DEVELOP
		      (*testout) << "ok,1 = " << ok << endl;
#endif

		      // there must be 2 different tangential faces to the edge
		      int cnt_tang_faces = 0;
		      for (int l = 0; l < surfind2.Size(); l++)
			{
			  Vec<3> nv =
			    ageometry.GetSurface(surfind2[l]) -> GetNormalVector(p);

			 
			  Vec<3> m1 = Cross (t, nv);
			  Vec<3> m2 = -m1;
			  bool isface1 = 0, isface2 = 0;
			  
			  Solid * locsol3;

			  // locsol2 -> TangentialSolid2 (p, m1, locsol3, surfind3, 1e-9*geomsize);
			  locsol -> TangentialEdgeSolid (p, t, t2, m1, locsol3, surfind3, ideps*geomsize);

			  //ageometry.GetIndependentSurfaceIndices (surfind3);

			  if (surfind3.Contains(surfind2[l]))
			    isface1 = 1;
			  delete locsol3;
			  
			  // locsol2 -> TangentialSolid2 (p, m2, locsol3, surfind3, 1e-9*geomsize);
			  locsol -> TangentialEdgeSolid (p, t, t2, m2, locsol3, surfind3, ideps*geomsize); 

			  // ageometry.GetIndependentSurfaceIndices (surfind3);

			  
			  if (surfind3.Contains(surfind2[l]))
			    isface2 = 1;
			  delete locsol3;

			  if (isface1 != isface2)
			    cnt_tang_faces++;
			}

#ifdef DEVELOP
		      (*testout) << "cnt_tang = " << cnt_tang_faces << endl;
#endif

		      if (cnt_tang_faces < 1)
			ok = false;

		      delete locsol2;
		      if (!ok) continue;
		    }

		  
		  // edge must be on tangential surface
		  bool isedge = 
		    locsol->VectorIn (p, t) &&
		    !locsol->VectorStrictIn (p, t);
		  
#ifdef DEVELOP
		  (*testout) << "isedge,1 = " << isedge << "\n";
#endif		
  
		  // there must exist at least two different faces on edge
		  if (isedge)
		    {
		      // *testout << "succ 1" << endl;
		      int cnts = 0;
		      for (int m = 0; m < surfind.Size(); m++)
			{
			  if (fabs (normalvecs[m] * t) > 1e-6)
			    continue;
			  
			  Vec<3> s = Cross (normalvecs[m], t);
			  Vec<3> t2a = t + 0.01 *s;
			  Vec<3> t2b = t - 0.01 *s;
			  
			  bool isface =
			    (locsol->VectorIn (p, t2a, 1e-6*geomsize) &&
			     !locsol->VectorStrictIn (p, t2a, 1e-6*geomsize))
			    ||
			    (locsol->VectorIn (p, t2b, 1e-6*geomsize) &&
			     !locsol->VectorStrictIn (p, t2b, 1e-6*geomsize));
			  
			  /*
			  bool isface =
			    (locsol->VectorIn (p, t2a) &&
			     !locsol->VectorStrictIn (p, t2a))
			    ||
			    (locsol->VectorIn (p, t2b) &&
			     !locsol->VectorStrictIn (p, t2b));
			  */

			  if (isface)
			    {
			      cnts++;
			    }
			}
		      if (cnts < 2) isedge = 0;
		    }
		  
		  if (isedge)
		    {
#ifdef DEVELOP
		      *testout << "success" << endl;
#endif
		      int spi = -1;
		      
		      const double searchradius = 1e-4*geomsize;//1e-5*geomsize;
		      searchtree.GetIntersecting (apoints[i]-Vec3d(searchradius,searchradius,searchradius), 
						  apoints[i]+Vec3d(searchradius,searchradius,searchradius), 
						  locsearch);
		      
		      for (int m = 0; m < locsearch.Size(); m++)
			{
			  if (Dist2 (specpoints[locsearch[m]].p, apoints[i]) < 1e-10*geomsize
			      && Abs2(specpoints[locsearch[m]].v - t) < 1e-8)
			    {
			      spi = locsearch[m];
			      break;
			    }
			}
		      
		      
		      if (spi == -1)
			{
			  spi = specpoints.Append (SpecialPoint()) - 1;
			  specpoint2point.Append (i);
			  specpoints.Last().unconditional = 0;
			  searchtree.Insert (apoints[i], spi);
			}

		      if(!specpoints[spi].unconditional)
			{
			  specpoints[spi].p = apoints[i];
			  specpoints[spi].v = t;
			  //if (surfind.Size() >= 3)
			  if (num_indep_surfs >= 3)
			    specpoints[spi].unconditional = 1;
			  specpoints[spi].s1 = rep_surfind[j];
			  specpoints[spi].s2 = rep_surfind[k];
			  specpoints[spi].s1_orig = surfind[j];
			  specpoints[spi].s2_orig = surfind[k];
			  specpoints[spi].layer = apoints[i].GetLayer();
			  for (int up = 0; up < geometry->GetNUserPoints(); up++)
			    if (Dist (geometry->GetUserPoint(up), apoints[i]) < 1e-8*geomsize)
			      specpoints[spi].unconditional = 1;
			  for (int ip = 0; ip < geometry->GetNIdentPoints(); ip++)
			    if (Dist (geometry->GetIdentPoint(ip), apoints[i]) < 1e-8*geomsize)
			      specpoints[spi].unconditional = 1;
			}
		    }
		  
		}

	    delete locsol;
	  }
      }

    /*
    BitArray testuncond (specpoints.Size());
    testuncond.Clear();
    for(int i = 0; i<specpoints.Size(); i++)
      {
	if(testuncond.Test(i))
	  continue;
	
	ARRAY<int> same;
	same.Append(i);
	
	for(int j = i+1; j<specpoints.Size(); j++)
	  {
	    if(Dist(specpoints[i].p,specpoints[j].p) < 1e-20)
	      {
		same.Append(j);
		testuncond.Set(j);
	      }
	  }
	
	if(same.Size() < 3)
	  for(int j=0; j<same.Size(); j++)
	    {
	      (*testout) << "setting " << specpoints[same[j]].p << "; " << specpoints[same[j]].v << "; " 
			 <<specpoints[same[j]].unconditional << " to conditional" << endl;
	      specpoints[same[j]].unconditional=0;
	    }
      }
    */


    // if special point is unconditional on some solid,
    // it must be unconditional everywhere:

    BitArray uncond (apoints.Size());
    uncond.Clear();

    for (int i = 0; i < specpoints.Size(); i++)
      if (specpoints[i].unconditional)
	uncond.Set (specpoint2point[i]);
  
    for (int i = 0; i < specpoints.Size(); i++)
      specpoints[i].unconditional = 
	uncond.Test (specpoint2point[i]) ? 1 : 0;
  }
}
