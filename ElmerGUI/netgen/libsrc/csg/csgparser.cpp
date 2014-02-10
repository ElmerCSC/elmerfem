#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{
  //using namespace netgen;


  static kwstruct defkw[] =
    {
      { TOK_RECO,    "algebraic3d" },
      { TOK_SOLID,   "solid" },
      { TOK_TLO,     "tlo" },
      { TOK_CURVE2D, "curve2d" },
      { TOK_CURVE3D, "curve3d" },
      { TOK_BOUNDINGBOX, "boundingbox" },
      { TOK_OR,      "or" },
      { TOK_AND,     "and" },
      { TOK_NOT,     "not" },
      { TOK_SINGULAR, "singular" },
      { TOK_EDGE,     "edge" },
      { TOK_FACE,     "face" },
      { TOK_POINT,    "point" },
      { TOK_IDENTIFY, "identify" },
      { TOK_CLOSESURFACES, "closesurfaces" },
      { TOK_CLOSEEDGES, "closeedges" },
      { TOK_PERIODIC,  "periodic" },
      { TOK_BOUNDARYCONDITION, "boundarycondition" },
      { TOK_BOUNDARYCONDITIONNAME, "boundaryconditionname" },
      { TOK_DEFINE, "define" },
      { TOK_CONSTANT, "constant" },
      { TOKEN_TYPE(0) }
    };

  static primstruct defprim[] =
    {
      { TOK_PLANE,     "plane" },
      { TOK_SPHERE,    "sphere" },
      { TOK_CYLINDER,  "cylinder" },
      { TOK_CONE,      "cone" },
      { TOK_ELLIPTICCYLINDER, "ellipticcylinder" },
      { TOK_ELLIPSOID, "ellipsoid" },
      { TOK_ORTHOBRICK, "orthobrick" },
      { TOK_POLYHEDRON, "polyhedron" },
      { TOK_TORUS, "torus" },

      { TOK_TUBE,      "tube" },
      { TOK_GENCYL,    "gencyl" },
      { TOK_EXTRUSION,  "extrusion" },
      { TOK_REVOLUTION, "revolution" },

      { TOK_TRANSLATE, "translate" },
      { TOK_MULTITRANSLATE, "multitranslate" },
      { TOK_ROTATE,   "rotate" },
      { TOK_MULTIROTATE, "multirotate" },
      { PRIMITIVE_TYPE(0) }
    };

  static CSGeometry * geom;


  CSGScanner :: CSGScanner (istream & ascanin)
  {
    scanin = &ascanin;
    token = TOK_END;
    num_value = 0;
    linenum = 1;
  }


  void CSGScanner :: ReadNext ()
  {
    char ch;
  

    // scan whitespaces
    do
      { 
	scanin->get(ch);

	//if (ch == '\n') 
	//  linenum++;

	// end of file reached
	if (scanin->eof())
	  {
	    token = TOK_END;
	    return;
	  }
	if (ch == '\n') 
	  linenum++;


	// skip comment line
	if (ch == '#')
	  {
	    while (ch != '\n')
	      {
		scanin->get(ch);
		if (scanin->eof())
		  {
		    token = TOK_END;
		    return;
		  }
	      }
	    linenum++;
	  }	
      }
    while (isspace(ch));
  
    switch (ch)
      {
      case '(': case ')': 
      case '[': case ']': 
      case '-':
      case '=': case ',': case ';':
	{
	  token = TOKEN_TYPE (ch);
	  break;
	}
  
      default:
	{
	  if (isdigit (ch) || ch == '.')
	    {
	      scanin->putback (ch);
	      (*scanin) >> num_value;
	      token = TOK_NUM;
	      return;
	    }

	  if (isalpha (ch))
	    {
	      string_value = string (1, ch);
	      scanin->get(ch);
	      while (isalnum(ch) || ch == '_')
		{
		  string_value += ch;
		  scanin->get(ch);
		}
	      scanin->putback (ch);
	    }

	  int nr = 0;
	  while (defkw[nr].kw)
	    {
	      if (string_value == defkw[nr].name)
		{
		  token = defkw[nr].kw;
		  return;
		}
	      nr++;
	    }

	  nr = 0;
	  while (defprim[nr].kw)
	    {
	      if (string_value == defprim[nr].name)
		{
		  token = TOK_PRIMITIVE;
		  prim_token = defprim[nr].kw;
		  return;
		}
	      nr++;
	    }

	  token = TOK_STRING;
	}
      }
  }

  void CSGScanner :: Error (const string & err)
  {
    stringstream errstr;
    errstr << "Parsing error in line " << linenum << ": " << endl << err << endl;
    throw string(errstr.str());
  }


  /*
    Solid = Term { OR Term }
    Term  = Primary { AND Primary }
    Primary = PRIM | IDENT | ( Solid ) | NOT Primary
  */

  void ParseChar (CSGScanner & scan, char ch)
  {
    if (scan.GetToken() != TOKEN_TYPE(ch)) 
      scan.Error (string ("token '") + string(1, ch) + string("' expected"));
    scan.ReadNext();
  }
  
  double ParseNumber(CSGScanner & scan)
  {
    if (scan.GetToken() == '-')
      {
	scan.ReadNext();
	return -ParseNumber (scan);
      }
    if (scan.GetToken() != TOK_NUM) scan.Error ("number expected");
    double val = scan.GetNumValue();
    scan.ReadNext();
    return val;
  }

  Vec<3> ParseVector (CSGScanner & scan)
  {
    Vec<3> v;
    v(0) = ParseNumber (scan);
    ParseChar (scan, ',');
    v(1) = ParseNumber (scan);
    ParseChar (scan, ',');
    v(2) = ParseNumber (scan);
    return v;
  }


  CSGScanner & operator>> (CSGScanner & scan, char ch)
  {
    if (scan.GetToken() != TOKEN_TYPE(ch)) 
      scan.Error (string ("token '") + string(1, ch) + string("' expected"));
    scan.ReadNext();
    return scan;
  }

  CSGScanner & operator>> (CSGScanner & scan, double & d)
  {
    d = ParseNumber (scan);
    return scan;
  }

  CSGScanner & operator>> (CSGScanner & scan, int & i)
  {
    i = int (ParseNumber (scan));
    return scan;
  }

  CSGScanner & operator>> (CSGScanner & scan, Point<3> & p)
  {
    scan >> p(0) >> ',' >> p(1) >> ',' >> p(2);
    return scan;
  }

  CSGScanner & operator>> (CSGScanner & scan, Vec<3> & v)
  {
    scan >> v(0) >> ',' >> v(1) >> ',' >> v(2);
    return scan;
  }


  Solid * ParseSolid (CSGScanner & scan);
  Solid * ParseTerm (CSGScanner & scan);
  Solid * ParsePrimary (CSGScanner & scan);
 

  Solid * ParsePrimary (CSGScanner & scan)
  {
    if (scan.GetToken() == TOK_PRIMITIVE)
      {
	switch (scan.GetPrimitiveToken())
	  {
	  case TOK_PLANE:
	    {
	      Point<3> p;
	      Vec<3> v;
	      
	      scan.ReadNext();
	      scan >> '(' >> p >> ';' >> v >> ')';

	      OneSurfacePrimitive * surf = new Plane ( p, v );
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }

	  case TOK_CYLINDER:
	    {
	      Point<3> pa, pb;
	      double r;
	      
	      scan.ReadNext();
	      scan >> '(' >> pa >> ';' >> pb >> ';' >> r >> ')';

	      OneSurfacePrimitive * surf = new Cylinder ( pa, pb, r );
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }

	  case TOK_ELLIPTICCYLINDER:
	    {
	      Point<3> pa;
	      Vec<3> vl, vs;
	      
	      scan.ReadNext();
	      scan >> '(' >> pa >> ';' >> vl >> ';' >> vs >> ')';

	      OneSurfacePrimitive * surf = new EllipticCylinder ( pa, vl, vs);
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }


	  case TOK_ELLIPSOID:
	    {
	      Point<3> pa;
	      Vec<3> v1, v2, v3;
	      
	      scan.ReadNext();
	      scan >> '(' >> pa >> ';' >> v1 >> ';' >> v2 >> ';' >> v3 >> ')';

	      OneSurfacePrimitive * surf = new Ellipsoid ( pa, v1, v2, v3);
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }


	  case TOK_CONE:
	    {
	      Point<3> pa, pb;
	      double ra, rb;
	      
	      scan.ReadNext();
	      scan >> '(' >> pa >> ';' >> ra >> ';' >> pb >> ';' >> rb >> ')';

	      OneSurfacePrimitive * surf = new Cone ( pa, pb, ra, rb );
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }



	  case TOK_SPHERE:
	    {
	      Point<3> p;
	      double r;
	      
	      scan.ReadNext();
	      scan >> '(' >> p >> ';' >> r >> ')';

	      OneSurfacePrimitive * surf = new Sphere ( p, r );
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }

	  case TOK_ORTHOBRICK:
	    {
	      Point<3> pa, pb;
	      
	      scan.ReadNext();
	      scan >> '(' >> pa >> ';' >> pb >> ')';
	      

	      Primitive * nprim = new OrthoBrick (pa, pb);
	      geom->AddSurfaces (nprim);
	      return new Solid (nprim);
	    } 

	  case TOK_POLYHEDRON:
	    {
	      // Added by Dalibor Lukas, October 15, 2003

	      Point<3> p;
	      //int pi1, pi2, pi3, pi4;
	      
	      scan.ReadNext();
	      ParseChar (scan, '(');
	      
	      Polyhedra * polyhedron = new Polyhedra;

	      // scanning the points
	      while (1)
		{
		  p = Point<3> (ParseVector (scan));
		  ParseChar (scan, ';');

		  polyhedron->AddPoint(p);

		  if (scan.GetToken() == ';')
		    {
		      scan.ReadNext();
		      break;
		    }
		}

	      // scanning the faces
	      int inputface = 0;
	      while (1)
		{
		  ARRAY<int> pnums,cleaned_pnums;
		  for(int i=0; i<3; i++)
		    {
		      pnums.Append((int) (ParseNumber (scan)));
		      if(i<2) 
			ParseChar (scan, ',');
		    }

		  if (scan.GetToken() == TOK_COMMA)
		    {
		      ParseChar (scan, ',');
		      pnums.Append((int) (ParseNumber (scan)));	
		    }

		  for(int i=0; i<pnums.Size(); i++)
		    if(!cleaned_pnums.Contains(pnums[i]))
		      cleaned_pnums.Append(pnums[i]);

		  if(cleaned_pnums.Size() == 3)
		    {
		      polyhedron->AddFace(cleaned_pnums[0]-1,
					  cleaned_pnums[1]-1,
					  cleaned_pnums[2]-1,
					  inputface);
		    }
		  else if(cleaned_pnums.Size() == 4)
		    {
		      polyhedron->AddFace(cleaned_pnums[0]-1,
					  cleaned_pnums[1]-1,
					  cleaned_pnums[2]-1,
					  inputface);
		      polyhedron->AddFace(cleaned_pnums[0]-1,
					  cleaned_pnums[2]-1,
					  cleaned_pnums[3]-1,
					  inputface);
		    }
		  else
		    {
		      ostringstream msg;
		      msg << "Something wrong with polyhedron face:";
		      for(int i=0; i<pnums.Size(); i++)
			msg << " " << pnums[i];
		      throw NgException(msg.str());
		    }
		  
		      

		  if (scan.GetToken() == ')')
		    {
		      scan.ReadNext();
		      break;
		    }
		  scan.ReadNext();
		  inputface++;
		}

	      geom->AddSurfaces (polyhedron);
	      return new Solid (polyhedron);
	    }


	  case TOK_REVOLUTION:
	    {
	      Point<3> p0,p1;

	      scan.ReadNext();
	      scan >> '(' >> p0 >> ';' >> p1 >> ';';

	      string spline = scan.GetStringValue();
	      
	      scan.ReadNext();
	      scan >> ')';
	      
	      if(!geom->GetSplineCurve2d(spline))
		{
		  scan.Error ( string("2D Spline curve not found: ") + spline );
		  break;
		}

	      Primitive * nprim = new Revolution(p0,p1,
						 *(geom->GetSplineCurve2d(spline)));

	      geom->AddSurfaces (nprim);
	      return new Solid(nprim);
	    }


	  case TOK_EXTRUSION: 
	    {   
	      scan.ReadNext();
	      scan >> '(';
	      string epath = scan.GetStringValue();
	      scan.ReadNext();
	      scan >> ';';
	      string profile = scan.GetStringValue();

	      
	      scan.ReadNext();
	      Vec<3> z_dir;
	      scan >> ';' >> z_dir(0) >> ',' >> z_dir(1) >> ',' >> z_dir(2) >> ')';
	      
	      if(!geom->GetSplineCurve2d(profile))
		{
		  scan.Error ( string("2D Spline curve not found: ") + profile );
		  break;
		}
	      if(!geom->GetSplineCurve3d(epath))
		{
		  scan.Error ( string("2D Spline curve not found: ") + epath );
		  break;
		}
	      
	      Primitive * nprim = new Extrusion(*(geom->GetSplineCurve3d(epath)),
						*(geom->GetSplineCurve2d(profile)),
						z_dir);
	      geom->AddSurfaces (nprim);
	      return new Solid(nprim);
	    }


	  /// Torus 
    	  /// Lorenzo Codecasa (codecasa@elet.polimi.it)
    	  /// April 27th, 2005 
	  ///
	  /// begin...
	  case TOK_TORUS:
	    {     
	      Point<3> pc; 
	      Vec<3> vn;
	      double R, r;
	      
	      scan.ReadNext();
	      scan >> '(' >> pc >> ';' >> vn >> ';' >> R >> ';' >> r >> ')';

	      OneSurfacePrimitive * surf = new Torus ( pc, vn, R, r );
	      geom->AddSurfaces (surf);
	      return new Solid (surf);
	    }
	  /// ..end




	  case TOK_TRANSLATE: 
	    {
	      Vec<3> v;
	      scan.ReadNext();

	      ParseChar (scan, '(');
	      v = ParseVector (scan);
	      ParseChar (scan, ';');
	      
	      Solid * sol1 = ParseSolid (scan);

	      ParseChar (scan, ')');

	      Solid * nsol = sol1 -> Copy(*geom);
	      Transformation<3> trans(v);
	      nsol -> Transform (trans);
	      return nsol;
	    }


	  case TOK_ROTATE: 
	    {
	      Point<3> c;
	      Vec<3> v;
	      scan.ReadNext();

	      scan >> '(' >> c >> ';' >> v >> ';';

	      Solid * sol1 = ParseSolid (scan);

	      ParseChar (scan, ')');

	      Solid * nsol = sol1 -> Copy(*geom);
	      Transformation<3> trans(c,v(0),v(1),v(2));
	      nsol -> Transform (trans);
	      return nsol;
	    }


	  case TOK_MULTITRANSLATE: 
	    {
	      Vec<3> v;
	      int n;
	      
	      scan.ReadNext();

	      scan >> '(' >> v >> ';' >> n >> ';';

	      Solid * sol1 = ParseSolid (scan);
	      
	      scan >> ')';
	      
	      Solid * hsol = sol1;
	      for (int i = 1; i <= n; i++)
		{
		  Solid * nsol = sol1 -> Copy(*geom);
		  Transformation<3> trans(double(i) * v);
		  
		  nsol -> Transform (trans);
		  hsol = new Solid (Solid::UNION, hsol, nsol); 
		}
	      return hsol;
	    }


	  case TOK_MULTIROTATE: 
	    {
	      Point<3> c;
	      Vec<3> v;
	      int n;
	      
	      scan.ReadNext();

	      scan >> '(' >> c >> ';' >> v >> ';' >> n >> ';';
	      Solid * sol1 = ParseSolid (scan);
	      scan >> ')';

	      Transformation<3> trans(c, v(0), v(1), v(2));
	      Transformation<3> multi(Vec<3>(0,0,0));
	      Transformation<3> ht;

	      Solid * hsol = sol1;
	      for (int i = 1; i <= n; i++)
		{
		  Solid * nsol = sol1 -> Copy(*geom);

		  nsol -> Transform (multi);
		  hsol = new Solid (Solid::UNION, hsol, nsol); 

		  ht=multi;
		  multi.Combine (trans, ht);
		}
	      return hsol;
	    }


	  default:
	    {
	      scan.Error (string ("unknown primary ") + scan.GetStringValue());
	    }

	  }
      }

    else if (scan.GetToken() == TOK_STRING &&
	     geom->GetSolid(scan.GetStringValue()))

      {
	Solid * sol = const_cast<Solid*> (geom->GetSolid(scan.GetStringValue()));
	scan.ReadNext();
	return sol;
      }

    else if (scan.GetToken() == TOK_NOT)

      {
	scan.ReadNext();
	Solid * sol1 = ParsePrimary (scan);
	return new Solid (Solid::SUB, sol1);
      }

    else if (scan.GetToken() == '(')

      {
	scan.ReadNext();
	Solid * sol1 = ParseSolid (scan);
	scan.ReadNext();
	return sol1;
      }

    scan.Error (string ("not a primary, name = ")+
		scan.GetStringValue());
    return 0;
  }



  Solid * ParseTerm (CSGScanner & scan)
  {
    Solid * sol = ParsePrimary(scan);
    while (scan.GetToken() == TOK_AND)
      {
	scan.ReadNext();
	Solid * sol2 = ParsePrimary(scan);
	sol = new Solid (Solid::SECTION, sol, sol2);
      }
    return sol;
  }


  Solid * ParseSolid (CSGScanner & scan)
  {
    Solid * sol = ParseTerm(scan);
    while (scan.GetToken() == TOK_OR)
      {
	scan.ReadNext();
	Solid * sol2 = ParseTerm(scan);
	sol = new Solid (Solid::UNION, sol, sol2);
      }
    return sol;
  }



  void ParseFlags (CSGScanner & scan, Flags & flags)
  {
    while (scan.GetToken() == '-')
      {
	scan.ReadNext();
	string name = scan.GetStringValue();
	scan.ReadNext();
	if (scan.GetToken() == '=')
	  {
	    scan.ReadNext();
	    if (scan.GetToken() == TOK_STRING)
	      {
		flags.SetFlag (name.c_str(), scan.GetStringValue().c_str());
		scan.ReadNext();
	      }
	    else if (scan.GetToken() == '[')
	      {
		scan.ReadNext();

		if(scan.GetToken() == '-' || scan.GetToken() == TOK_NUM)
		  {
		    ARRAY<double> vals;
		    vals.Append (ParseNumber(scan));
		    while (scan.GetToken() == ',')
		      {
			scan.ReadNext();
			vals.Append (ParseNumber(scan));
		      }
		    ParseChar (scan, ']');
		    flags.SetFlag (name.c_str(), vals);
		  }
		else
		  { // string list
		    ARRAY<char*> vals;
		    string val = scan.GetStringValue();
		    vals.Append(new char[val.size()+1]);
		    strcpy(vals.Last(),val.c_str());
		    scan.ReadNext();

		    while (scan.GetToken() == ',')
		      {
			scan.ReadNext();
			val = scan.GetStringValue();
			vals.Append(new char[val.size()+1]);
			strcpy(vals.Last(),val.c_str());
			scan.ReadNext();
		      }
		    ParseChar (scan, ']');
		    flags.SetFlag (name.c_str(), vals);
		    for(int i=0; i<vals.Size(); i++)
		      delete [] vals[i];
		  }
	      }
	    else if (scan.GetToken() == TOK_NUM)
	      {
		flags.SetFlag (name.c_str(), scan.GetNumValue());
		scan.ReadNext();
	      }
	  }     
	else
	  {
	    flags.SetFlag (name.c_str());
	  }
      }
  }


  /*
    Main parsing function for CSG geometry
  */
  CSGeometry * ParseCSG (istream & istr)
  {
    CSGScanner scan(istr);
    
    geom = new CSGeometry;

    scan.ReadNext();
    if (scan.GetToken() != TOK_RECO)  // keyword 'algebraic3d'
      return 0;

    scan.ReadNext();

    try
      {
	while (1)
	  {
	    if (scan.GetToken() == TOK_END) break;
	    
	    if (scan.GetToken() == TOK_SOLID)
	      {
		scan.ReadNext();
		if (scan.GetToken() != TOK_STRING)
		  scan.Error ("name identifier expected");
		string solidname = scan.GetStringValue();

		scan.ReadNext();

		ParseChar (scan, '=');
		Solid * solid = ParseSolid (scan);

		Flags flags;
		ParseFlags (scan, flags);

		geom->SetSolid (solidname.c_str(), new Solid (Solid::ROOT, solid)); 
		geom->SetFlags (solidname.c_str(), flags); 
		
		ParseChar (scan, ';');
		
		PrintMessage (4, "define solid ", solidname);
	      }

	    else if (scan.GetToken() == TOK_TLO)

	      { // a TopLevelObject definition

		scan.ReadNext();
		
		string name = scan.GetStringValue();
		scan.ReadNext();

		if (scan.GetToken() != TOK_STRING)

		  { // a solid TLO

		    Flags flags;
		    ParseFlags (scan, flags);
		    
		    ParseChar (scan, ';');
		    if (!geom->GetSolid (name))
		      scan.Error ("Top-Level-Object "+name+" not defined");

		    int tlonr = 
		      geom->SetTopLevelObject ((Solid*)geom->GetSolid(name));
		    TopLevelObject * tlo = geom->GetTopLevelObject (tlonr);

		    if (flags.NumListFlagDefined ("col"))
		      {
			const ARRAY<double> & col =
			  flags.GetNumListFlag ("col");
			tlo->SetRGB (col[0], col[1], col[2]);
		      }

		    if (flags.GetDefineFlag ("transparent"))
		      tlo->SetTransparent (1);

		    tlo->SetMaterial (flags.GetStringFlag ("material", ""));
		    tlo->SetLayer (int(flags.GetNumFlag ("layer", 1)));
		    if (flags.NumFlagDefined ("maxh"))
		      tlo->SetMaxH (flags.GetNumFlag("maxh", 1e10));
		  }

		else
		  
		  { // a surface TLO

		    string surfname = scan.GetStringValue();
		    scan.ReadNext();

		    Flags flags;
		    ParseFlags (scan, flags);
		    
		    ParseChar (scan, ';');

		    ARRAY<int> si;
		    geom->GetSolid(surfname)->GetSurfaceIndices(si);
		    int tlonr = 
		      geom->SetTopLevelObject ((Solid*)geom->GetSolid(name),
					       (Surface*)geom->GetSurface(si.Get(1)));
		    TopLevelObject * tlo = geom->GetTopLevelObject (tlonr);
		    if (flags.NumListFlagDefined ("col"))
		      {
			const ARRAY<double> & col = flags.GetNumListFlag ("col");
			tlo->SetRGB (col.Get(1), col.Get(2), col.Get(3));
		      }
		    if (flags.GetDefineFlag ("transparent"))
		      tlo->SetTransparent (1);

		    if (flags.NumFlagDefined ("maxh"))
		      tlo->SetMaxH (flags.GetNumFlag("maxh", 1e10));
		    tlo->SetLayer (int(flags.GetNumFlag ("layer", 1)));
		    tlo->SetBCProp (int(flags.GetNumFlag ("bc", -1)));
		    if ( flags.StringFlagDefined("bcname") )
		      tlo->SetBCName ( flags.GetStringFlag ("bcname", "default") );
		  }
	      }
	    
	    else if (scan.GetToken() == TOK_IDENTIFY)

	      {
		
		scan.ReadNext();
		switch (scan.GetToken())
		  {
		  case TOK_CLOSESURFACES:
		    {
		      scan.ReadNext();
		      
		      string name1 = scan.GetStringValue();
		      scan.ReadNext();
		      
		      string name2 = scan.GetStringValue();
		      scan.ReadNext();

		      Flags flags;
		      ParseFlags (scan, flags);
		      
		      ParseChar (scan, ';');
		      
		      
		      ARRAY<int> si1, si2;
		      geom->GetSolid(name1)->GetSurfaceIndices(si1);
		      geom->GetSolid(name2)->GetSurfaceIndices(si2);

		      const TopLevelObject * domain = 0;
		      if (flags.StringFlagDefined ("tlo"))
			{
			  domain = 
			    geom->GetTopLevelObject (geom->GetSolid(flags.GetStringFlag ("tlo","")));
			  if (!domain) 
			    scan.Error ("identification needs undefined tlo");
			}

		      geom->AddIdentification 
			(new CloseSurfaceIdentification 
			 (geom->GetNIdentifications()+1, *geom, 
			  geom->GetSurface (si1[0]), geom->GetSurface (si2[0]),
			  domain,
			  flags));

		      break;
		    }
		    
		  case TOK_PERIODIC:
		    {
		      scan.ReadNext();
		      
		      string name1 = scan.GetStringValue();
		      scan.ReadNext();

		      string name2 = scan.GetStringValue();
		      scan.ReadNext();

		      ParseChar (scan, ';');

		      
		      ARRAY<int> si1, si2;
		      geom->GetSolid(name1)->GetSurfaceIndices(si1);
		      geom->GetSolid(name2)->GetSurfaceIndices(si2);
		      
		      geom->AddIdentification 
			(new PeriodicIdentification 
			 (geom->GetNIdentifications()+1,
			  *geom,
			  geom->GetSurface (si1.Get(1)),
			  geom->GetSurface (si2.Get(1))));
		      break;
		    }

		  default:
		    scan.Error ("keyword 'closesurfaces' or 'periodic' expected");
		  }
		
	      }

	    else if (scan.GetToken() == TOK_SINGULAR)

	      {
		
		scan.ReadNext();
		switch (scan.GetToken())
		  {
		  case TOK_FACE:
		    {
		      scan.ReadNext();
		      
		      string name1 = scan.GetStringValue();  // tlo
		      scan.ReadNext();
		      
		      string name2 = scan.GetStringValue();
		      scan.ReadNext();
		      
		      Flags flags;
		      ParseFlags (scan, flags);
		      int factor = int(flags.GetNumFlag("factor",1)); 
		      // cout << "Singular Face with factor " << factor << endl; 
		      PrintMessageCR (3, "Singular Face  with factor ", factor);

		      ParseChar (scan, ';');
		      
		      const Solid * sol = geom->GetSolid(name2);

		      if(!sol)
			scan.Error ("unknown solid in singular face definition");
		      else
			for (int i = 0; i < geom->GetNTopLevelObjects(); i++)
			  if (name1 == geom->GetTopLevelObject (i)->GetSolid()->Name())
			    geom->singfaces.Append (new SingularFace (i+1, sol,factor));

		      break;
		    }

		  case TOK_EDGE:
		    {
		      scan.ReadNext();
		      
		      string name1 = scan.GetStringValue();
		      scan.ReadNext();
		      
		      string name2 = scan.GetStringValue();
		      scan.ReadNext();
		      
		      Flags flags;
		      ParseFlags (scan, flags);
		      int factor = int(flags.GetNumFlag("factor",1));
		      double maxhinit = flags.GetNumFlag("maxh",-1);
		      ParseChar (scan, ';');
		      
		      const Solid * s1 = geom->GetSolid(name1);
		      const Solid * s2 = geom->GetSolid(name2);
		      PrintMessageCR (3, "Singular Edge  with factor ", factor);

		      int domnr = -1;
		      if (flags.StringFlagDefined ("tlo"))
			{
			  const Solid * sol =
			    geom->GetSolid(flags.GetStringFlag ("tlo",""));
			  
			  for (int i = 0; i < geom->GetNTopLevelObjects(); i++)
			    if (geom->GetTopLevelObject(i)->GetSolid() == sol)
			      domnr = i;
			  
			  // cout << "domnr = " << domnr;
			}

		      if(!s1 || !s2)
			scan.Error ("unknown solid ins singular edge definition");
		      else
			geom->singedges.Append (new SingularEdge (1, domnr, 
								  *geom, s1, s2, factor,
								  maxhinit));
		      break;
		    }

		  case TOK_POINT:
		    {
		      scan.ReadNext();
		      
		      string name1 = scan.GetStringValue();
		      scan.ReadNext();
		      string name2 = scan.GetStringValue();
		      scan.ReadNext();
		      string name3 = scan.GetStringValue();
		      scan.ReadNext();
		      
		      Flags flags;
		      ParseFlags (scan, flags);
		      int factor = int(flags.GetNumFlag("factor",1)); 
		      ParseChar (scan, ';');
		      
		      const Solid * s1 = geom->GetSolid(name1);
		      const Solid * s2 = geom->GetSolid(name2);
		      const Solid * s3 = geom->GetSolid(name3);
		      // cout << "Singular Point with factor " << factor << endl; 
		      PrintMessageCR (3, "Singular Point  with factor ", factor);
		      geom->singpoints.Append (new SingularPoint (1, s1, s2, s3, factor));
		      break;
		    }
		  default:
		    scan.Error ("keyword 'face' or 'edge' or 'point' expected");
		  }
	      }

	    
	    else if (scan.GetToken() == TOK_POINT)
	      {
		Point<3> p;

		scan.ReadNext();
		ParseChar (scan, '(');
		p = Point<3> (ParseVector (scan));
		ParseChar (scan, ')');


		Flags flags;
		ParseFlags (scan, flags);
		int factor = int(flags.GetNumFlag("factor",0)); 

		ParseChar (scan, ';');

		geom->AddUserPoint (p, factor);
	      }

	    else if (scan.GetToken() == TOK_BOUNDINGBOX)
	      {
		Point<3> p1, p2;
		
		scan.ReadNext();
		ParseChar (scan, '(');
		p1 = Point<3> (ParseVector (scan));
		ParseChar (scan, ';');
		p2 = Point<3> (ParseVector (scan));
		ParseChar (scan, ')');
		ParseChar (scan, ';');

		geom->SetBoundingBox (Box<3> (p1, p2));
	      }

	    else if (scan.GetToken() == TOK_CURVE2D)
	      {
		scan.ReadNext();

		
		if (scan.GetToken() != TOK_STRING)
		  scan.Error ("name identifier expected");
		string curvename = scan.GetStringValue();

		scan.ReadNext();

		ParseChar (scan, '=');
		ParseChar (scan, '(');
		
		SplineGeometry<2> * newspline = new SplineGeometry<2>;
		newspline->CSGLoad(scan);

		ParseChar (scan, ')');
		ParseChar (scan, ';');

		geom->SetSplineCurve(curvename.c_str(),newspline);

		PrintMessage (4, "define 2d curve ", curvename);
	      }

	    else if (scan.GetToken() == TOK_CURVE3D)
	      {
		scan.ReadNext();

		
		if (scan.GetToken() != TOK_STRING)
		  scan.Error ("name identifier expected");
		string curvename = scan.GetStringValue();

		scan.ReadNext();

		ParseChar (scan, '=');
		ParseChar (scan, '(');
		
		SplineGeometry<3> * newspline = new SplineGeometry<3>;
		newspline->CSGLoad(scan);

		ParseChar (scan, ')');
		ParseChar (scan, ';');

		geom->SetSplineCurve(curvename.c_str(),newspline);

		PrintMessage (4, "define 3d curve ", curvename);
	      }

	    else if (scan.GetToken() == TOK_BOUNDARYCONDITION)
	      {
		scan.ReadNext();
		
		string name1 = scan.GetStringValue();
		scan.ReadNext();
		
		string name2 = scan.GetStringValue();
		scan.ReadNext();
		
		int num = int (ParseNumber (scan));
		ParseChar (scan, ';');


		CSGeometry::BCModification bcm;
		bcm.bcname = NULL;
		ARRAY<int> si;
		
		geom->GetSolid(name1)->GetSurfaceIndices(si);
		if(si.Size() == 0)
		  {
		    string errstring = "solid \""; errstring += name1; errstring += "\" has no surfaces";
		    scan.Error (errstring);
		  }
	
		bcm.tlonr = -1;
		int i;	
		for (i = 0; i < geom->GetNTopLevelObjects(); i++)
		  if (string (geom->GetTopLevelObject(i)->GetSolid()->Name())
		      == name2)
		    {
		      bcm.tlonr = i;
		      break;
		    }
		if(bcm.tlonr == -1)
		  {
		    string errstring = "tlo \""; errstring += name2; errstring += "\" not found";
		    scan.Error(errstring);
		  }
		
		
		bcm.bcnr = num;
		for (i = 0; i < si.Size(); i++)
		  {
		    bcm.si = si[i];
		    geom->bcmodifications.Append (bcm);
		  }
	      }
	    
	    else if (scan.GetToken() == TOK_BOUNDARYCONDITIONNAME)
	      {
		scan.ReadNext();
		
		string name1 = scan.GetStringValue();
		scan.ReadNext();
		
		string name2 = scan.GetStringValue();
		scan.ReadNext();

		string bcname = scan.GetStringValue();
		scan.ReadNext();
		ParseChar(scan, ';');
		

		CSGeometry::BCModification bcm;
		bcm.bcname = NULL;


		ARRAY<int> si;
		
		geom->GetSolid(name1)->GetSurfaceIndices(si);
		if(si.Size() == 0)
		  {
		    string errstring = "solid \""; errstring += name1; errstring += "\" has no surfaces";
		    scan.Error (errstring);
		  }
	
		bcm.tlonr = -1;
		int i;	
		for (i = 0; i < geom->GetNTopLevelObjects(); i++)
		  if (string (geom->GetTopLevelObject(i)->GetSolid()->Name())
		      == name2)
		    {
		      bcm.tlonr = i;
		      break;
		    }
		if(bcm.tlonr == -1)
		  {
		    string errstring = "tlo \""; errstring += name2; errstring += "\" not found";
		    scan.Error(errstring);
		  }
		
		
		bcm.bcnr = -1;
		for (i = 0; i < si.Size(); i++)
		  {
		    bcm.si = si[i];
		    geom->bcmodifications.Append (bcm);
		    geom->bcmodifications.Last().bcname = new string(bcname);
		  }
	      }
	    
	    else if (scan.GetToken() == TOK_DEFINE)
	      {
		scan.ReadNext();
		string name;
		double val;
		
		switch (scan.GetToken())
		  {
		  case TOK_CONSTANT:
		    scan.ReadNext();
		      
		    name = scan.GetStringValue();
		    scan.ReadNext();

		    ParseChar(scan, '=');
		    val = ParseNumber(scan);

		    if(name == "identprec")
		      geom->SetIdEps(val);
		    
		    

		    break;
		  default:
		    scan.Error ("keyword 'constant' expected");
		  }
	      }


	    else
	      {
		cout << "read unidentified token " << scan.GetToken() 
		     << " (as char: \"" << char(scan.GetToken()) << "\")"
		     << " string = " << scan.GetStringValue() << endl;
		scan.ReadNext();
	      }
	  }
      }
    catch (string errstr)
      {
	cout << "caught error " << errstr << endl;
	throw NgException (errstr);
      }



    (*testout) << geom->GetNTopLevelObjects() << " TLOs:" << endl;
    for (int i = 0; i < geom->GetNTopLevelObjects(); i++)
      {
	const TopLevelObject * tlo = geom->GetTopLevelObject(i);
	if (tlo->GetSolid())
	  (*testout) << i << ": " << *tlo->GetSolid() << endl;
      }

    (*testout) << geom->GetNSurf() << " Surfaces" << endl;
    for (int i = 0; i < geom->GetNSurf(); i++)
      (*testout) << i << ": " << *geom->GetSurface(i) << endl;

    return geom;
    /*
      do
      {
      scan.ReadNext();
      if (scan.GetToken() == TOK_STRING)
      cout << "found string " << scan.GetStringValue() << endl;
      else
      cout << "token = " << int(scan.GetToken()) << endl;
      }
      while (scan.GetToken() != TOK_END);
    */
  }


};

