%{
//define YYDEBUG 1

extern int yylex ();

#include <mystdlib.h>

#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>

namespace netgen
{
netgen::CSGeometry * parsegeom;
}

using namespace netgen;

// extern ARRAY<Surface*> surfaces;
// extern SYMBOLTABLE<Solid*> solids;


int yyerror (char * s);
splinetube * tube;
spline3d * middlecurve;
Point<3> splinep1;
BSplineCurve2d *bspline;
Flags parseflags;
extern int linenum;
ARRAY<double> doublearray;
ARRAY<char*> stringarray;

Polyhedra * polyhedron;
// Revolution * revolution;
%}

%union {
double val;
char * chptr;
Solid * solidtype;
}

%token <val> NUM
%token TOK_SOLID, TOK_RECO, TOK_TLO, TOK_BOUNDINGBOX
%token <chptr> IDENT IDENTSOLID
%token <solidtype> TOK_SPHERE TOK_CYLINDER TOK_CONE TOK_PLAIN TOK_TUBE TOK_GENCYL TOK_ORTHOBRICK TOK_POLYHEDRON TOK_REVOLUTION
%left <solidtype> TOK_OR TOK_AND TOK_NOT
%token <solidtype> TOK_TRANSLATE TOK_MULTITRANSLATE TOK_ROTATE TOK_MULTIROTATE
%type <solidtype> solid solidprimitive 
%type <void> splinesegmentlist splinesegment readbspline bsplinepointlist
%type <chptr> anyident
%token TOK_SINGULAR TOK_EDGE TOK_POINT
%token TOK_IDENTIFY TOK_CLOSESURFACES TOK_CLOSEEDGES TOK_PERIODIC
%token TOK_BOUNDARYCONDITION
%type <void> polyhedronpoints polyhedronfaces polyhedronpoint polyhedronface
%type <void> revolutionpoints revolutionpoint


%%
input: 	
		{
		  linenum = 1;
		}
	TOK_RECO
	recsoliddef
	recadddef
	
		{
		int i;
		extern ARRAY<char*> parsestrings;
		for (i = 0; i < parsestrings.Size(); i++)
 		  delete [] parsestrings[i];
		parsestrings.SetSize(0);
		}
	;


recsoliddef:
	/* empty */
        | recsoliddef soliddef ';'
        ;

soliddef:
          TOK_SOLID IDENT '=' solid
                  { parsegeom->SetSolid($2, new Solid (Solid::ROOT, $4)); 
		  }
	  flaglist
		  { 
		  parsegeom->SetFlags($2, parseflags); 
		  }		
        ;

solid:
          solidprimitive
        | IDENTSOLID              { $$ = (Solid*)parsegeom->GetSolid($1);  }
        | solid TOK_OR solid      { $$ = new Solid (Solid::UNION, $1, $3); }
        | solid TOK_AND solid     { $$ = new Solid (Solid::SECTION, $1, $3); }
        | TOK_NOT solid           { $$ = new Solid (Solid::SUB, $2); }
        | '(' solid ')'           { $$ = $2; }
        ;

solidprimitive:
          TOK_SPHERE '(' NUM ',' NUM ',' NUM ';'
                         NUM ')'
                 { 
		   OneSurfacePrimitive * surf = new Sphere (Point<3> ($3, $5, $7), $9); 
		   parsegeom -> AddSurface (surf);
		   surf->SetSurfaceId (0, parsegeom->GetNSurf()-1);
                   $$ = new Solid (surf);
                 }
        | TOK_CYLINDER '(' NUM ',' NUM ',' NUM ';'
                           NUM ',' NUM ',' NUM ';'
                           NUM ')'

                 {
		   OneSurfacePrimitive * surf = new Cylinder (Point<3> ($3, $5, $7),
	                                       Point<3> ($9, $11, $13), $15);
		   parsegeom->AddSurface (surf);
		   surf->SetSurfaceId (0, parsegeom->GetNSurf()-1);
                   $$ = new Solid (surf);
                 }
        | TOK_CONE '(' NUM ',' NUM ',' NUM ';'
                       NUM ';' 
                       NUM ',' NUM ',' NUM ';'
                       NUM ')'

                 {
		   OneSurfacePrimitive * surf = new Cone (Point<3> ($3, $5, $7),
	                                     	 Point<3> ($11, $13, $15), $9, $17);
		   parsegeom->AddSurface (surf);
		   surf->SetSurfaceId (0, parsegeom->GetNSurf()-1);
                   $$ = new Solid (surf);
                 }
        | TOK_PLAIN '(' NUM ',' NUM ',' NUM ';'
                           NUM ',' NUM ',' NUM ')'
                 {
		   OneSurfacePrimitive * surf = new Plane ( Point<3> ($3, $5, $7),
	                                        Vec<3> ($9, $11, $13) );
		   parsegeom->AddSurface (surf);
		   surf->SetSurfaceId (0, parsegeom->GetNSurf()-1);
                   $$ = new Solid (surf);
                 }
	| TOK_ORTHOBRICK  '(' NUM ',' NUM ',' NUM ';'
	                      NUM ',' NUM ',' NUM ')'
		{
		  Primitive * nprim = new OrthoBrick (Point<3> ($3, $5, $7),
		                                      Point<3> ($9, $11, $13));
                  for (int j = 0; j < nprim->GetNSurfaces(); j++)
		    {
		      parsegeom->AddSurface (&nprim->GetSurface(j));
		      nprim->SetSurfaceId (j, parsegeom->GetNSurf()-1);
		      $$ = new Solid (nprim);
		    }
		} 


	| TOK_POLYHEDRON '(' 
		{
		  polyhedron = new Polyhedra ();
		}
		polyhedronpoints ';' ';'
		polyhedronfaces ')'
		{
		  int j;
                  for (j = 0; j < polyhedron->GetNSurfaces(); j++)
		    {
		      parsegeom->AddSurface (&polyhedron->GetSurface(j));
		      polyhedron->SetSurfaceId (j, parsegeom->GetNSurf()-1);
		      $$ = new Solid (polyhedron);
		    }	
		}


/*
	| TOK_REVOLUTION '(' NUM ',' NUM ',' NUM ';'
	                     NUM ',' NUM ',' NUM ';' 
		{
                    revolution = new Revolution (Point<3> ($3, $5, $7),
						 Point<3> ($9, $11, $13));
		}
 		    revolutionpoints
			  ')'
	        {
		  revolution -> Finish ();
		  int j;
                  for (j = 0; j < revolution->GetNSurfaces(); j++)
		    {
		      parsegeom->AddSurface (&revolution->GetSurface(j));
		      revolution->SetSurfaceId (j, parsegeom->GetNSurf()-1);
		      $$ = new Solid (revolution);
		    }	
	        }
*/


	| TOK_TRANSLATE  '(' NUM ',' NUM ',' NUM ';' solid ')'
		{
                  Solid * nsol = $9 -> Copy(*parsegeom);
                  Vec<3> v($3, $5, $7);
                  Transformation<3> trans(v);
		  nsol -> Transform (trans);
                  $$ = nsol;
		}

	| TOK_MULTITRANSLATE  '(' NUM ',' NUM ',' NUM ';' NUM ';' solid ')'
		{
		  int i;
                  Solid * hsol = $11;
	          for (i = 1; i <= $9; i++)
                    {
                      Solid * nsol = $11 -> Copy(*parsegeom);
                      Vec<3> v($3, $5, $7);
		      v *= i;
                      Transformation<3> trans(v);
		      nsol -> Transform (trans);
                      hsol = new Solid (Solid::UNION, hsol, nsol); 
                    }
		  $$ = hsol;
		}

	| TOK_ROTATE  '(' NUM ',' NUM ',' NUM ';' NUM ',' NUM ',' NUM ';' solid ')'
		{
                  Solid * nsol = $15 -> Copy(*parsegeom);
                  Point<3> c($3, $5, $7);
                  Transformation<3> rot(c, $9, $11, $13);
		  nsol -> Transform (rot);
                  $$ = nsol;
		}

	| TOK_MULTIROTATE  '(' NUM ',' NUM ',' NUM ';'
			       NUM ',' NUM ',' NUM ';'
			       NUM ';' solid ')'
		{
		  int i;
                  Solid * hsol = $17;                      

		  Point<3> c($3, $5, $7);
                  Transformation<3> trans(c, $9, $11, $13);
		  Transformation<3> multi(Vec<3>(0,0,0));
		  Transformation<3> ht;

	          for (i = 1; i <= $15; i++)
                    {
                      Solid * nsol = $17 -> Copy(*parsegeom);
		      nsol -> Transform (multi);
                      hsol = new Solid (Solid::UNION, hsol, nsol); 

		      ht=multi;
                      multi.Combine (trans, ht);
                    }
		  $$ = hsol;
		}


/*
	| TOK_TUBE '(' NUM ',' NUM ',' NUM ','
		{
		middlecurve = new spline3d;
		splinep1.X() = $3; splinep1.Y() = $5; splinep1.Z() = $7;
		}
	 splinesegmentlist ';' NUM ')'
		{
		   Surface * surf = new splinetube (*middlecurve, $12);
		   parsegeom->AddSurface (surf);
                   $$ = new Solid (surf, parsegeom->GetNSurf());
		}
		
	| TOK_GENCYL '(' NUM ',' NUM ',' NUM ';'
	                 NUM ',' NUM ',' NUM ';' 
	                 NUM ',' NUM ',' NUM ';'
			readbspline
			')'
	        {
		   Surface * surf = new GeneralizedCylinder
	           (*bspline, Point<3> ($3, $5, $7), Vec<3> ($9, $11, $13),
	             Vec<3> ($15, $17, $19) );
		   parsegeom->AddSurface (surf);
                   $$ = new Solid (surf, parsegeom->GetNSurf());
	        }
*/	             
	;


polyhedronpoints:
	  polyhedronpoints ';' polyhedronpoint
	| polyhedronpoint
	;
polyhedronfaces:
	  polyhedronfaces ';' polyhedronface
	| polyhedronface
	;
polyhedronpoint:
	  NUM ',' NUM ',' NUM
	{
		polyhedron->AddPoint (Point<3> ($1, $3, $5));
 		cout << " " << $1 << " " << $3 << " " << $5 << endl;
	}
	;
polyhedronface:
	  NUM ',' NUM ',' NUM
 	{ 
		polyhedron->AddFace (int($1)-1, int($3)-1, int($5)-1);
		cout << $1 << " " << $3 << " " << $5 << endl; 
	}
	| NUM ',' NUM ',' NUM ',' NUM
 	{ 	
		cout << "face, 1 = " << $1 << " " << $3 << " "  << $5 << " " << $7 << endl; 
		polyhedron->AddFace (int($1)-1, int($3)-1, int($5)-1);
		polyhedron->AddFace (int($1)-1, int($5)-1, int($7)-1);
		cout << $1 << $3 << $5 << $7 << endl; 
	}
	;

revolutionpoints:
	  revolutionpoints ';' revolutionpoint
	| revolutionpoint
	;
revolutionpoint:
	  NUM ',' NUM 
	{
//		revolution->AddPoint (Point<2> ($1, $3));
 		cout << " " << $1 << " " << $3 << endl;
	}
	;




	
splinesegmentlist:
	  splinesegment
        | splinesegment ',' splinesegmentlist
	;
splinesegment:
	  NUM ',' NUM ',' NUM ','
	  NUM ',' NUM ',' NUM
	 	{
	        middlecurve->AddSegment (splinep1, Point<3> ($1, $3, $5), Point<3> ($7, $9, $11));
		splinep1(0) = $7; splinep1(1) = $9; splinep1(2) = $11;
		}
	;
	
	
readbspline:
          NUM ',' NUM
                {
                bspline = new BSplineCurve2d;
                bspline -> AddPoint (Point<2> ($1, $3));
		cout << "first point" << endl;
                }
          bsplinepointlist
        ;
bsplinepointlist: 
          /* empty */ { } 
        | ',' NUM ',' NUM
               {
               bspline -> AddPoint (Point<2> ($2, $4));
		cout << "Add Point: " << $2 << "-" << $4 << endl;
               }
          bsplinepointlist
        ;









recadddef:
	/* empty */
        | recadddef adddef ';'
        ;

adddef:
	TOK_SINGULAR TOK_EDGE  NUM IDENTSOLID IDENTSOLID
	{ cout << "singular edge:" << $3 << " between "
		 << $4 << " and " << $5 << endl; 
	  parsegeom->singedges.Append 
	    (new SingularEdge ($3, parsegeom->GetSolid($4), 
				parsegeom->GetSolid($5)));
	}
	|
	TOK_SINGULAR TOK_POINT  NUM IDENTSOLID IDENTSOLID IDENTSOLID
	{ cout << "singular point:" << $3 << " between "
		 << $4 << ", " << $5 << " and " << $6 << endl; 
	  parsegeom->singpoints.Append 
	    (new SingularPoint ($3, parsegeom->GetSolid($4), 
				parsegeom->GetSolid($5),
				parsegeom->GetSolid($6)));
	}
	|
	TOK_IDENTIFY TOK_CLOSESURFACES IDENTSOLID IDENTSOLID flaglist
	{ 	
	  ARRAY<int> si1, si2;
	  parsegeom->GetSolid($3)->GetSurfaceIndices(si1);
	  parsegeom->GetSolid($4)->GetSurfaceIndices(si2);

	  parsegeom->AddIdentification (
		new CloseSurfaceIdentification (
			parsegeom->GetNIdentifications()+1,
			*parsegeom, 
			parsegeom->GetSurface (si1[0]),
			parsegeom->GetSurface (si2[0]),
			parseflags));
	}
	|
	TOK_IDENTIFY TOK_CLOSEEDGES IDENTSOLID IDENTSOLID IDENTSOLID
	{ 	
	  ARRAY<int> si1, si2, si3;
	  parsegeom->GetSolid($3)->GetSurfaceIndices(si1);
	  parsegeom->GetSolid($4)->GetSurfaceIndices(si2);
	  parsegeom->GetSolid($5)->GetSurfaceIndices(si3);

	  parsegeom->AddIdentification (
		new CloseEdgesIdentification (
			parsegeom->GetNIdentifications()+1,
			*parsegeom, 
			parsegeom->GetSurface (si1.Get(1)),
			parsegeom->GetSurface (si2.Get(1)),
			parsegeom->GetSurface (si3.Get(1))));
	}
	|
	TOK_IDENTIFY TOK_PERIODIC IDENTSOLID IDENTSOLID 
	{ 
	  ARRAY<int> si1, si2;
	  parsegeom->GetSolid($3)->GetSurfaceIndices(si1);
	  parsegeom->GetSolid($4)->GetSurfaceIndices(si2);

	  parsegeom->AddIdentification (
		new PeriodicIdentification (
			parsegeom->GetNIdentifications()+1,
			*parsegeom,
			parsegeom->GetSurface (si1.Get(1)),
			parsegeom->GetSurface (si2.Get(1))));
	}
	| 
	TOK_TLO IDENTSOLID flaglist
	{
	  int tlonr = 
            parsegeom->SetTopLevelObject ((Solid*)parsegeom->GetSolid($2));
	  TopLevelObject * tlo = parsegeom->GetTopLevelObject (tlonr);
	  if (parseflags.NumListFlagDefined ("col"))
	    {
	      const ARRAY<double> & col = parseflags.GetNumListFlag ("col");
	      tlo->SetRGB (col[0], col[1], col[2]);
	    }
	  if (parseflags.GetDefineFlag ("transparent"))
	    tlo->SetTransparent (1);
	}
	|
	TOK_TLO IDENTSOLID IDENTSOLID flaglist
	{
	  ARRAY<int> si;
	  parsegeom->GetSolid($3)->GetSurfaceIndices(si);
	  int tlonr = 
	    parsegeom->SetTopLevelObject ((Solid*)parsegeom->GetSolid($2),
					  (Surface*)parsegeom->GetSurface(si.Get(1)));
	  TopLevelObject * tlo = parsegeom->GetTopLevelObject (tlonr);
	  if (parseflags.NumListFlagDefined ("col"))
	    {
	      const ARRAY<double> & col = parseflags.GetNumListFlag ("col");
	      tlo->SetRGB (col.Get(1), col.Get(2), col.Get(3));
	    }
	  if (parseflags.GetDefineFlag ("transparent"))
	    tlo->SetTransparent (1);
	}
	|
	TOK_BOUNDINGBOX '(' NUM ',' NUM ',' NUM ';' 
			    NUM ',' NUM ',' NUM ')'
	{
	  parsegeom->SetBoundingBox (Box<3> (Point<3> ($3, $5, $7), 
                                             Point<3> ($9, $11, $13)));
	}
	| 
	TOK_POINT '(' NUM ',' NUM ',' NUM ')' 
	{
	  parsegeom->AddUserPoint (Point<3> ($3, $5, $7));
	}
	|
	TOK_BOUNDARYCONDITION IDENTSOLID IDENTSOLID NUM
	{
  	  CSGeometry::BCModification bcm;
	  ARRAY<int> si;

	  parsegeom->GetSolid($2)->GetSurfaceIndices(si);
	
          bcm.tlonr = -1;
	  int i;	
	  for (i = 0; i < parsegeom->GetNTopLevelObjects(); i++)
	    if (strcmp (parsegeom->GetTopLevelObject(i)->GetSolid()->Name(), $3) == 0)
	      {
	        bcm.tlonr = i;
	        break;
              }

	  bcm.bcnr = int($4);
	  for (i = 0; i < si.Size(); i++)
            {
	      bcm.si = si[i];
              parsegeom->bcmodifications.Append (bcm);
            }
	}
	;	




flaglist:
	      { parseflags.DeleteFlags (); }
	  recflaglist
	;

recflaglist:
	  /* empty */
	| flag recflaglist
	;
	
flag:
	  '-' anyident
	                   { parseflags.SetFlag ($2); }
	| '-' anyident '=' anyident
	                   { parseflags.SetFlag ($2, $4); }
	| '-' anyident '=' NUM
	                   { parseflags.SetFlag ($2, $4); }
	| '-' anyident '=' numlistbrack
			   { parseflags.SetFlag ($2, doublearray); }
	| '-' anyident '=' stringlistbrack
			   { parseflags.SetFlag ($2, stringarray); }
	;     	


numlistbrack:
	  '['   			{ doublearray.SetSize (0); }
	   numlist ']'
	;

numlist:
          NUM                        { doublearray.Append ($1); }
        | numlist ',' NUM            { doublearray.Append ($3); }
        ;

stringlistbrack:
	'['
	{
	  int i;
	  for (i = 0; i < stringarray.Size(); i++)
	     delete stringarray[i];
	  stringarray.SetSize (0);
	}
	  stringlist  ']'
	;

stringlist:
	  anyident
	        {
			stringarray.Append (new char[strlen($1)+1]);
			strcpy (stringarray.Last(), $1);
		}
					
	| stringlist ',' anyident
           	{
			stringarray.Append (new char[strlen($3)+1]);
			strcpy (stringarray.Last(), $3);
		}
	;





anyident:
	  IDENT            { $$ = $1; }
	| IDENTSOLID       { $$ = $1; }
%%


int yyerror (char * s)
{	
  cerr << s << " in line " << linenum << endl;
  return 0;
}

