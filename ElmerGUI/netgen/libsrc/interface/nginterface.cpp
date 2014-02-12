#include <mystdlib.h>


#include <meshing.hpp>
#include <csg.hpp>
#include <geometry2d.hpp>
#include <stlgeom.hpp>


#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif

#ifdef ACIS
#include <acisgeom.hpp>
#endif

#ifdef SOCKETS
#include "../sockets/sockets.hpp"
#endif

#ifndef NOTCL
#include <visual.hpp>
#endif


#include "nginterface.h"

// #include <FlexLexer.h>


// #include <mystdlib.h>


namespace netgen
{
#include "writeuser.hpp"

	extern AutoPtr<Mesh> mesh;
#ifndef NOTCL
  extern VisualSceneMesh vsmesh;
  extern Tcl_Interp * tcl_interp;
#endif

  extern AutoPtr<SplineGeometry2d> geometry2d;
  extern AutoPtr<CSGeometry> geometry;
  extern STLGeometry * stlgeometry;

#ifdef OCCGEOMETRY
  extern OCCGeometry * occgeometry;
#endif
#ifdef ACIS
  extern ACISGeometry * acisgeometry;
#endif

#ifdef OPENGL
  extern VisualSceneSolution vssolution;
#endif
  extern CSGeometry * ParseCSG (istream & istr);

#ifdef SOCKETS
  extern AutoPtr<ClientSocket> clientsocket;
  //extern ARRAY< AutoPtr < ServerInfo > > servers;
  extern ARRAY< ServerInfo* > servers;
#endif
}


using namespace netgen;

/*
  extern void * operator new (size_t s);
  extern void * operator new [] (size_t s);
  extern void operator delete (void * p);
  extern void operator delete [] (void * p);
*/

// extern FlexLexer * lexer;



void Ng_LoadGeometry (const char * filename)
{

  
  if (printmessage_importance>0)
    cout << "CALLED NG LOAD GEOMETRY" << endl; 
  
  geometry.Reset (new CSGeometry ());
  geometry2d.Reset ();

#ifdef OCCGEOMETRY
  delete occgeometry;
  occgeometry = 0;
#endif
#ifdef ACIS
  delete acisgeometry;
  acisgeometry = 0;
#endif
  
  // he: if filename is empty, return
  // can be used to reset geometry
  if (strcmp(filename,"")==0) 
    return;
     
  ifstream infile (filename);
  
  if ((strcmp (&filename[strlen(filename)-3], "geo") == 0) ||
      (strcmp (&filename[strlen(filename)-3], "GEO") == 0) ||
      (strcmp (&filename[strlen(filename)-3], "Geo") == 0))
    {
 

      geometry.Reset( netgen::ParseCSG(infile) );   

      if (!geometry)
	{
	  geometry.Reset (new CSGeometry ());
	  //throw NgException ("input file not found");
	  cerr << "Error: input file \"" << filename << "\" not found" << endl;
	}

      geometry -> FindIdenticSurfaces(1e-6);

#ifdef PARALLEL
      int id, rc, ntasks;
      MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      if ( id > 0 )
	{
	  geometry->CalcTriangleApproximation ( geometry->BoundingBox(), 0.001, 20 );
	  return;
	}
#endif

      Box<3> box (geometry->BoundingBox());
#ifdef NOTCL
      double detail = 0.001;
      double facets = 20;
      geometry->CalcTriangleApproximation(box, detail, facets);
      
#else
      double detail = atof (Tcl_GetVar (tcl_interp, "::geooptions.detail", 0));
      double facets = atof (Tcl_GetVar (tcl_interp, "::geooptions.facets", 0));
      
      if (atoi (Tcl_GetVar (tcl_interp, "::geooptions.drawcsg", 0)))
	geometry->CalcTriangleApproximation(box, detail, facets);
#endif
    }

  else if (strcmp (&filename[strlen(filename)-4], "in2d") == 0)
    {
      geometry2d.Reset (new SplineGeometry2d());
      geometry2d -> Load (filename);
    }

  else if ((strcmp (&filename[strlen(filename)-3], "stl") == 0) ||
	   (strcmp (&filename[strlen(filename)-3], "STL") == 0) ||
	   (strcmp (&filename[strlen(filename)-3], "Stl") == 0))
    {
      stlgeometry = STLGeometry :: Load (infile);
      stlgeometry->edgesfound = 0;
      Mesh meshdummy;
      stlgeometry->Clear();
      stlgeometry->BuildEdges();
      stlgeometry->MakeAtlas(meshdummy);
      stlgeometry->CalcFaceNums();
      stlgeometry->AddFaceEdges();
      stlgeometry->LinkEdges();
    }

#ifdef OCCGEOMETRY
  else if ((strcmp (&filename[strlen(filename)-4], "iges") == 0) ||
	   (strcmp (&filename[strlen(filename)-3], "igs") == 0) ||
	   (strcmp (&filename[strlen(filename)-3], "IGS") == 0) ||
	   (strcmp (&filename[strlen(filename)-4], "IGES") == 0))
    {
      PrintMessage (1, "Load IGES geometry file ", filename);
      occgeometry = LoadOCC_IGES (filename);
    }
  else if ((strcmp (&filename[strlen(filename)-4], "step") == 0) ||
	   (strcmp (&filename[strlen(filename)-3], "stp") == 0) ||
	   (strcmp (&filename[strlen(filename)-3], "STP") == 0) ||
	   (strcmp (&filename[strlen(filename)-4], "STEP") == 0))
    {
      PrintMessage (1, "Load STEP geometry file ", filename);
      occgeometry = LoadOCC_STEP (filename);
    }
#endif

#ifdef ACIS
  else if (
		  strcmp (&filename[strlen(filename)-3], "sat") == 0 ||
		  ( strlen(filename) >= 7 && strcmp ( &filename[ strlen( filename)-7 ], "sat.tet" ) == 0 )
		)
    {
      PrintMessage (1, "Load ACIS geometry file ", filename);
      acisgeometry = LoadACIS_SAT (filename);
    }
#endif
  else
    {
      //throw NgException("Unknown geometry extension");
		cerr << "Error: Unknown geometry extension \"" << filename[strlen(filename)-3] << "\"" << endl;
    }


}                          


void Ng_LoadMeshFromStream ( istream & input )
{
  mesh.Reset (new Mesh());
  mesh -> Load(input);
  if(input.good())
    {
      string auxstring;
      input >> auxstring;
      if(auxstring == "csgsurfaces")
	{
	  if (geometry)
            {
              geometry.Reset (new CSGeometry (""));
            }
	  if (stlgeometry)
	    {
	      delete stlgeometry;
	      stlgeometry = NULL;
	    }
#ifdef OCCGEOMETRY
	  if (occgeometry)
	    {
 	      delete occgeometry;
 	      occgeometry = NULL;
	    }
#endif
#ifdef ACIS
	  if (acisgeometry)
	    {
 	      delete acisgeometry;
 	      acisgeometry = NULL;
	    }
#endif
 	  geometry2d.Reset (0);
	  
 	  geometry -> LoadSurfaces(input);
	}
    }
 
}


void Ng_LoadMesh (const char * filename)
{
  if ( (strlen (filename) > 4) &&
       strcmp (filename + (strlen (filename)-4), ".vol") != 0 )
    {
      mesh.Reset (new Mesh());
      ReadFile(*mesh,filename);

      //mesh->SetGlobalH (mparam.maxh);
      //mesh->CalcLocalH();
      return;
    }

  ifstream infile(filename);
  Ng_LoadMeshFromStream(infile);
}

void Ng_LoadMeshFromString (char * mesh_as_string)
{
  istringstream instream(mesh_as_string);
  Ng_LoadMeshFromStream(instream);
}
  



int Ng_GetDimension ()
{
  return (mesh) ? mesh->GetDimension() : -1;
}

int Ng_GetNP ()
{
  return (mesh) ? mesh->GetNP() : 0;
}

int Ng_GetNV ()
{
  return (mesh) ? mesh->GetNV() : 0;
}

int Ng_GetNE ()
{
  if(!mesh) return 0;
  if (mesh->GetDimension() == 3)
    return mesh->GetNE();
  else
    return mesh->GetNSE();
}

int Ng_GetNSE ()
{
  if(!mesh) return 0;
  if (mesh->GetDimension() == 3)
    return mesh->GetNSE();
  else
    return mesh->GetNSeg();
}

void Ng_GetPoint (int pi, double * p)
{
  if (pi < 1 || pi > mesh->GetNP())
    {
      if (printmessage_importance>0)
        cout << "Ng_GetPoint: illegal point " << pi << endl;
      return;
    }

  const Point3d & hp = mesh->Point (pi);
  p[0] = hp.X();
  p[1] = hp.Y();
  if (mesh->GetDimension() == 3)
    p[2] = hp.Z();
}


NG_ELEMENT_TYPE Ng_GetElement (int ei, int * epi, int * np)
{
  if (mesh->GetDimension() == 3)
    {
      int i;
      const Element & el = mesh->VolumeElement (ei);
      for (i = 0; i < el.GetNP(); i++)
	epi[i] = el.PNum(i+1);
      
      if (np)
	*np = el.GetNP();

      if (el.GetType() == PRISM)
	{
	  // degenerated prism, (should be obsolete)
	  const int map1[] = { 3, 2, 5, 6, 1 };
	  const int map2[] = { 1, 3, 6, 4, 2 };
	  const int map3[] = { 2, 1, 4, 5, 3 };
	  
	  const int * map = NULL;
	  int deg1 = 0, deg2 = 0, deg3 = 0;
	  //int deg = 0;
	  if (el.PNum(1) == el.PNum(4)) { map = map1; deg1 = 1; }
	  if (el.PNum(2) == el.PNum(5)) { map = map2; deg2 = 1; }
	  if (el.PNum(3) == el.PNum(6)) { map = map3; deg3 = 1; }
	  
	  switch (deg1+deg2+deg3)
	    {
	      {
	      case 1:
                if (printmessage_importance>0)
                  cout << "degenerated prism found, deg = 1" << endl;
		for (i = 0; i < 5; i++)
		  epi[i] = el.PNum (map[i]);
		
		if (np) *np = 5;
		return NG_PYRAMID;
		break;
	      }
	    case 2:
	      {
                if (printmessage_importance>0)
                  cout << "degenerated prism found, deg = 2" << endl;
		if (!deg1) epi[3] = el.PNum(4);
		if (!deg2) epi[3] = el.PNum(5);
		if (!deg3) epi[3] = el.PNum(6);
		
		if (np) *np = 4;
		return NG_TET;
		break;
	      }
	    default:
	      ;
	    }
	  
	}

      return NG_ELEMENT_TYPE (el.GetType());
    }
  else
    {
      int i;
      const Element2d & el = mesh->SurfaceElement (ei);
      for (i = 0; i < el.GetNP(); i++)
	epi[i] = el.PNum(i+1);      

      if (np) *np = el.GetNP();
      return NG_ELEMENT_TYPE (el.GetType());
      /*
      switch (el.GetNP())
	{
	case 3: return NG_TRIG; 
	case 4: return NG_QUAD; 
	case 6: return NG_TRIG6; 
	}
      */
    }

  // should not occur
  return NG_TET;
}


NG_ELEMENT_TYPE Ng_GetElementType (int ei)
{
  if (mesh->GetDimension() == 3)
    {
      return NG_ELEMENT_TYPE (mesh->VolumeElement (ei).GetType());
    }
  else
    {
      const Element2d & el = mesh->SurfaceElement (ei);
      switch (el.GetNP())
	{
	case 3: return NG_TRIG; 
	case 4: return NG_QUAD; 
	case 6: return NG_TRIG6; 
	}
    }

  // should not occur
  return NG_TET;
}



int Ng_GetElementIndex (int ei)
{
  if (mesh->GetDimension() == 3)
    return mesh->VolumeElement(ei).GetIndex();
  else
    {
      int ind = mesh->SurfaceElement(ei).GetIndex(); 
      ind = mesh->GetFaceDescriptor(ind).BCProperty();
      return ind;
    }
}

void Ng_SetElementIndex(const int ei, const int index)
{
  mesh->VolumeElement(ei).SetIndex(index);
}

char * Ng_GetElementMaterial (int ei)
{
  static char empty[] = "";
  if (mesh->GetDimension() == 3)
    {
      int ind = mesh->VolumeElement(ei).GetIndex();
      // cout << "ind = " << ind << endl;
      const char * mat = mesh->GetMaterial (ind);
      if (mat)
	return const_cast<char*> (mat);
      else 
	return empty;
    }
  // add astrid
  else
    {
      int ind = mesh->SurfaceElement(ei).GetIndex();
      ind = mesh->GetFaceDescriptor(ind).BCProperty();
      const char * mat = mesh->GetMaterial ( ind );
      if (mat)
	return const_cast<char*> (mat);
      else
	return empty;
    }
  return 0;
}

char * Ng_GetDomainMaterial (int dom)
{
  static char empty[] = "";
  // astrid
  if ( 1 ) // mesh->GetDimension() == 3)
    {
      const char * mat = mesh->GetMaterial(dom);
      if (mat)
	return const_cast<char*> (mat);
      else 
	return empty;      
    }

  return 0;
}


NG_ELEMENT_TYPE Ng_GetSurfaceElement (int ei, int * epi, int * np)
{
  if (mesh->GetDimension() == 3)
    {
      const Element2d & el = mesh->SurfaceElement (ei);
      for (int i = 0; i < el.GetNP(); i++)
	epi[i] = el[i];
      
      if (np) *np = el.GetNP();
      
      return NG_ELEMENT_TYPE (el.GetType());
    }
  else
    {
      const Segment & seg = mesh->LineSegment (ei);

      if (seg.pmid < 0)
	{
	  epi[0] = seg.p1;
	  epi[1] = seg.p2;
	  
	  if (np) *np = 2;
	  return NG_SEGM;
	}
      else
	{
	  epi[0] = seg.p1;
	  epi[1] = seg.p2;
	  epi[2] = seg.pmid;

	  if (np) *np = 3;
	  return NG_SEGM3;
	}
    }

  return NG_TRIG;
}

int Ng_GetSurfaceElementIndex (int ei)
{
  if (mesh->GetDimension() == 3)
    return mesh->GetFaceDescriptor(mesh->SurfaceElement(ei).GetIndex()).BCProperty();
  else
    return mesh->LineSegment(ei).si;
}

int Ng_GetSurfaceElementSurfaceNumber (int ei)
{
  if (mesh->GetDimension() == 3)
    return mesh->GetFaceDescriptor(mesh->SurfaceElement(ei).GetIndex()).SurfNr();
  else
    return mesh->LineSegment(ei).si;
}
int Ng_GetSurfaceElementFDNumber (int ei)
{
  if (mesh->GetDimension() == 3)
    return mesh->SurfaceElement(ei).GetIndex();
  else
    return -1;
}


char * Ng_GetSurfaceElementBCName (int ei)
{
  if ( mesh->GetDimension() == 3 )
    return const_cast<char *>(mesh->GetFaceDescriptor(mesh->SurfaceElement(ei).GetIndex()).GetBCName().c_str());
  else
    return const_cast<char *>(mesh->LineSegment(ei).GetBCName().c_str());
}


// Inefficient (but maybe safer) version:
//void Ng_GetSurfaceElementBCName (int ei, char * name)
//{
//  if ( mesh->GetDimension() == 3 )
//      strcpy(name,mesh->GetFaceDescriptor(mesh->SurfaceElement(ei).GetIndex()).GetBCName().c_str());
//  else
//      strcpy(name,mesh->LineSegment(ei).GetBCName().c_str());
//}

char * Ng_GetBCNumBCName (int bcnr)
{
  return const_cast<char *>(mesh->GetBCName(bcnr).c_str());
}


// Inefficient (but maybe safer) version:
//void Ng_GetBCNumBCName (int bcnr, char * name)
//{
//    strcpy(name,mesh->GetBCName(bcnr).c_str());
//}

void Ng_GetNormalVector (int sei, int locpi, double * nv)
{
  nv[0] = 0; 
  nv[1] = 0;
  nv[2] = 1;

  (*testout) << "Ng_GetNormalVector (sei = " << sei << ", locpi = " << locpi << ")" << endl;
  
  if (mesh->GetDimension() == 3)
    {
      Vec<3> n;
      Point<3> p;
      p = mesh->Point (mesh->SurfaceElement(sei).PNum(locpi));

      int surfi = mesh->GetFaceDescriptor(mesh->SurfaceElement(sei).GetIndex()).SurfNr();

      (*testout) << "surfi = " << surfi << endl;
#ifdef OCCGEOMETRY
      if (occgeometry)
	{
	  PointGeomInfo gi = mesh->SurfaceElement(sei).GeomInfoPi(locpi);
	  occgeometry->GetSurface (surfi).GetNormalVector(p, gi, n);
	  nv[0] = n(0);
	  nv[1] = n(1);
	  nv[2] = n(2);
	}
      else
#endif
      if (geometry)
	{
	  (*testout) << "geometry defined" << endl;
	  n = geometry->GetSurface (surfi) -> GetNormalVector(p);
	  (*testout) << "aus is" << endl;
	  nv[0] = n(0);
	  nv[1] = n(1);
	  nv[2] = n(2);
	}
    }
}



void Ng_SetPointSearchStartElement(const int el)
{
  mesh->SetPointSearchStartElement(el);
}


int Ng_FindElementOfPoint (double * p, double * lami, int build_searchtree, 
			   const int * const indices, const int numind)
  
{
  ARRAY<int> * dummy(NULL);
  int ind = -1;

  if(indices != NULL)
    {
      dummy = new ARRAY<int>(numind);
      for(int i=0; i<numind; i++) (*dummy)[i] = indices[i];
    }

  if (mesh->GetDimension() == 3)
    {
      Point3d p3d(p[0], p[1], p[2]);
      ind = 
	mesh->GetElementOfPoint(p3d, lami, dummy, build_searchtree != 0);
    }
  else
    {
      double lam3[3];
      Point3d p2d(p[0], p[1], 0);
      ind = 
	mesh->GetElementOfPoint(p2d, lam3, dummy, build_searchtree != 0);
      

      if(mesh->SurfaceElement(ind).GetType()==QUAD)
	{
	  lami[0] = lam3[0];
	  lami[1] = lam3[1];
	}
      else 
	{
	  lami[0] = 1-lam3[0]-lam3[1];
	  lami[1] = lam3[0];
	}
    }

  delete dummy;

  return ind;
}

int Ng_FindSurfaceElementOfPoint (double * p, double * lami, int build_searchtree, 
				  const int * const indices, const int numind)
  
{
  ARRAY<int> * dummy(NULL);
  int ind = -1;

  if(indices != NULL)
    {
      dummy = new ARRAY<int>(numind);
      for(int i=0; i<numind; i++) (*dummy)[i] = indices[i];
    }

  if (mesh->GetDimension() == 3)
    {
      Point3d p3d(p[0], p[1], p[2]);
      ind = 
	mesh->GetSurfaceElementOfPoint(p3d, lami, dummy, build_searchtree != 0);
    }
  else
    {
      //throw NgException("FindSurfaceElementOfPoint for 2D meshes not yet implemented");
		cerr << "FindSurfaceElementOfPoint for 2D meshes not yet implemented" << endl;
    }

  delete dummy;

  return ind;
}


int Ng_IsElementCurved (int ei)
{
  if (mesh->GetDimension() == 2)
    return mesh->GetCurvedElements().IsSurfaceElementCurved (ei-1);
  else
    return mesh->GetCurvedElements().IsElementCurved (ei-1);
}


int Ng_IsSurfaceElementCurved (int sei)
{
  if (mesh->GetDimension() == 2)
    return mesh->GetCurvedElements().IsSegmentCurved (sei-1);
  else
    return mesh->GetCurvedElements().IsSurfaceElementCurved (sei-1);
}




void Ng_GetElementTransformation (int ei, const double * xi, 
				  double * x, double * dxdxi)
{
  if (mesh->GetDimension() == 2)
    {
      Point<2> xl(xi[0], xi[1]);
      Point<3> xg;
      Mat<3,2> dx;

      mesh->GetCurvedElements().CalcSurfaceTransformation (xl, ei-1, xg, dx);

      if (x)
	{
	  for (int i = 0; i < 2; i++)
	    x[i] = xg(i);
	}
	  
      if (dxdxi)
	{
	  for (int i=0; i<2; i++)
	    {
	      dxdxi[2*i] = dx(i,0);
	      dxdxi[2*i+1] = dx(i,1);
	    }
	}
    }
  else
    {
      Point<3> xl(xi[0], xi[1], xi[2]);
      Point<3> xg;
      Mat<3,3> dx;

      mesh->GetCurvedElements().CalcElementTransformation (xl, ei-1, xg, dx);

      // still 1-based arrays
      if (x)
	{
	  for (int i = 0; i < 3; i++)
	    x[i] = xg(i);
	}

      if (dxdxi)
	{
	  for (int i=0; i<3; i++)
	    {
	      dxdxi[3*i] = dx(i,0);
	      dxdxi[3*i+1] = dx(i,1);
              dxdxi[3*i+2] = dx(i,2);
	    }
	}
    }
}



void Ng_GetBufferedElementTransformation (int ei, const double * xi, 
                                          double * x, double * dxdxi,
                                          void * buffer, int buffervalid)
{
  // buffer = 0;
  // buffervalid = 0;
  if (mesh->GetDimension() == 2)
    {
      return Ng_GetElementTransformation (ei, xi, x, dxdxi);
    }
  else
    {
      mesh->GetCurvedElements().CalcElementTransformation (reinterpret_cast<const Point<3> &> (*xi), 
                                                           ei-1, 
                                                           reinterpret_cast<Point<3> &> (*x), 
                                                           reinterpret_cast<Mat<3,3> &> (*dxdxi), 
                                                           buffer, (buffervalid != 0));

      /*
      Point<3> xl(xi[0], xi[1], xi[2]);
      Point<3> xg;
      Mat<3,3> dx;
      // buffervalid = 0;
      mesh->GetCurvedElements().CalcElementTransformation (xl, ei-1, xg, dx, buffer, buffervalid);

      // still 1-based arrays
      if (x)
	{
	  for (int i = 0; i < 3; i++)
	    x[i] = xg(i);
	}

      if (dxdxi)
	{
	  for (int i=0; i<3; i++)
	    {
	      dxdxi[3*i] = dx(i,0);
	      dxdxi[3*i+1] = dx(i,1);
              dxdxi[3*i+2] = dx(i,2);
	    }
	}
      */
    }
}








void Ng_GetMultiElementTransformation (int ei, int n,
                                       const double * xi, int sxi,
                                       double * x, int sx,
                                       double * dxdxi, int sdxdxi)
{
  if (mesh->GetDimension() == 2)
    {
      for (int i = 0; i < n; i++)
        {
          Point<2> xl(xi[i*sxi], xi[i*sxi+1]);
          Point<3> xg;
          Mat<3,2> dx;

          mesh->GetCurvedElements().CalcSurfaceTransformation (xl, ei-1, xg, dx);

          if (x)
            {
              x[i*sx  ] = xg(0);
              x[i*sx+1] = xg(1);
            }
	  
          if (dxdxi)
            {
              dxdxi[i*sdxdxi  ] = dx(0,0);
              dxdxi[i*sdxdxi+1] = dx(0,1);
              dxdxi[i*sdxdxi+2] = dx(1,0);
              dxdxi[i*sdxdxi+3] = dx(1,1);
            }
        }
    }
  else
    {
      mesh->GetCurvedElements().CalcMultiPointElementTransformation (ei-1, n, xi, sxi, x, sx, dxdxi, sdxdxi);
    }
}








void Ng_GetSurfaceElementTransformation (int sei, const double * xi, 
					 double * x, double * dxdxi)
{
  if (mesh->GetDimension() == 2)
    {
      Point<3> xg;
      Vec<3> dx;

      mesh->GetCurvedElements().CalcSegmentTransformation (xi[0], sei-1, xg, dx);

      if (x)
        for (int i = 0; i < 2; i++)
	  x[i] = xg(i);
	  
      if (dxdxi)
        for (int i=0; i<2; i++)
	  dxdxi[i] = dx(i);

    }
  else
    {
      Point<2> xl(xi[0], xi[1]);
      Point<3> xg;
      Mat<3,2> dx;
      
      mesh->GetCurvedElements().CalcSurfaceTransformation (xl, sei-1, xg, dx);
      
      for (int i=0; i<3; i++)
	{
	  if (x)
	    x[i] = xg(i);
	  if (dxdxi)
	    {
	      dxdxi[2*i] = dx(i,0);
	      dxdxi[2*i+1] = dx(i,1);
	    }
	}
    }
}



void Ng_GetSurfaceElementNeighbouringDomains(const int selnr, int & in, int & out)
{
  if ( mesh->GetDimension() == 3 )
    {
      in = mesh->GetFaceDescriptor(mesh->SurfaceElement(selnr).GetIndex()).DomainIn();
      out = mesh->GetFaceDescriptor(mesh->SurfaceElement(selnr).GetIndex()).DomainOut();
    }
  else
    {
      in = mesh -> LineSegment(selnr) . domin;
      out = mesh -> LineSegment(selnr) . domout;
    }
}


#ifdef PARALLEL
// Is Element ei an element of this processor ??
bool Ng_IsGhostEl (int ei)
{
  if ( mesh->GetDimension() == 3 )
    return mesh->VolumeElement(ei).IsGhost();
  else
    return false;
}

void Ng_SetGhostEl(const int ei, const bool aisghost )
{
  if ( mesh -> GetDimension () == 3 )
    mesh -> VolumeElement(ei).SetGhost (aisghost);
}

bool Ng_IsGhostSEl (int ei)
{
  if ( mesh -> GetDimension () == 3 )
    return mesh->SurfaceElement(ei).IsGhost();
  else
    return false;
}

void Ng_SetGhostSEl(const int ei, const bool aisghost )
{
  if ( mesh -> GetDimension () == 3 )
    mesh -> SurfaceElement(ei).SetGhost (aisghost);
}


bool Ng_IsGhostVert ( int pnum )
{
  return mesh -> Point ( pnum ).IsGhost() ;  
}
bool Ng_IsGhostEdge ( int ednum )
{
  return mesh -> GetParallelTopology() . IsGhostEdge ( ednum ); 
}

bool Ng_IsGhostFace ( int fanum )
{
  return mesh -> GetParallelTopology() . IsGhostFace ( fanum ); 
}

// void Ng_SetGhostVert ( const int pnum, const bool aisghost );
// void Ng_SetGhostEdge ( const int ednum, const bool aisghost );
// void Ng_SetGhostFace ( const int fanum, const bool aisghost );


bool Ng_IsExchangeEl ( int elnum )
{ return mesh -> GetParallelTopology() . IsExchangeElement ( elnum ); }

bool Ng_IsExchangeSEl ( int selnum )
{ return mesh -> GetParallelTopology() . IsExchangeSEl ( selnum ); }

void Ng_UpdateOverlap()
{ mesh->UpdateOverlap(); }

int Ng_Overlap ()
{ return mesh->GetParallelTopology() . Overlap(); }

#endif

void Ng_SetRefinementFlag (int ei, int flag)
{
  if (mesh->GetDimension() == 3)
    {
      mesh->VolumeElement(ei).SetRefinementFlag (flag != 0);
      mesh->VolumeElement(ei).SetStrongRefinementFlag (flag >= 10);
    }
  else
    {
      mesh->SurfaceElement(ei).SetRefinementFlag (flag != 0);
      mesh->SurfaceElement(ei).SetStrongRefinementFlag (flag >= 10);
    }
}

void Ng_SetSurfaceRefinementFlag (int ei, int flag)
{
  if (mesh->GetDimension() == 3)
    {
      mesh->SurfaceElement(ei).SetRefinementFlag (flag != 0);
      mesh->SurfaceElement(ei).SetStrongRefinementFlag (flag >= 10);
    }
}


void Ng_Refine (NG_REFINEMENT_TYPE reftype)
{
  NgLock meshlock (mesh->MajorMutex(), 1);

  BisectionOptions biopt;
  biopt.usemarkedelements = 1;
  biopt.refine_p = 0;
  biopt.refine_hp = 0;
  if (reftype == NG_REFINE_P)
    biopt.refine_p = 1;
  if (reftype == NG_REFINE_HP)
    biopt.refine_hp = 1;
  Refinement * ref;
  MeshOptimize2d * opt = NULL;

  if (geometry2d)
    ref = new Refinement2d(*geometry2d);
  else if (stlgeometry)
    ref = new RefinementSTLGeometry(*stlgeometry);
#ifdef OCCGEOMETRY
  else if (occgeometry)
    ref = new OCCRefinementSurfaces (*occgeometry);
#endif
#ifdef ACIS
  else if (acisgeometry)
    {
      ref = new ACISRefinementSurfaces (*acisgeometry);
      opt = new ACISMeshOptimize2dSurfaces(*acisgeometry);
      ref->Set2dOptimizer(opt);
    }
#endif
  else if (geometry && mesh->GetDimension() == 3)
    {
      ref = new RefinementSurfaces(*geometry);
      opt = new MeshOptimize2dSurfaces(*geometry);
      ref->Set2dOptimizer(opt);
    }
  else
    {
      ref = new Refinement();
    }


  ref -> Bisect (*mesh, biopt);

  mesh -> UpdateTopology();

  // mesh -> GetCurvedElements().BuildCurvedElements (ref, mparam.elementorder);
  delete ref;
  delete opt;
}

void Ng_SecondOrder ()
{
  if (stlgeometry)
    {
      RefinementSTLGeometry ref (*stlgeometry);
      ref.MakeSecondOrder (*mesh);
    }

  else if (geometry2d)
    {
      Refinement2d ref (*geometry2d);
      ref.MakeSecondOrder (*mesh);
    }

  else if (geometry && mesh->GetDimension() == 3)

    {
      RefinementSurfaces ref (*geometry);
      ref.MakeSecondOrder (*mesh);
    }
  else
    {
      if (printmessage_importance>0)
        cout << "no geom" << endl;
      Refinement ref;
      ref.MakeSecondOrder (*mesh);
    }

  mesh -> UpdateTopology();
}

/*
void Ng_HPRefinement (int levels)
{
  Refinement * ref;

  if (stlgeometry)
    ref = new RefinementSTLGeometry (*stlgeometry);
  else if (geometry2d)
    ref = new Refinement2d (*geometry2d);
  else
    ref = new RefinementSurfaces (*geometry);


  HPRefinement (*mesh, ref, levels);
}

void Ng_HPRefinement (int levels, double parameter)
{
  Refinement * ref;

  if (stlgeometry)
    ref = new RefinementSTLGeometry (*stlgeometry);
  else if (geometry2d)
    ref = new Refinement2d (*geometry2d);
  else
    ref = new RefinementSurfaces (*geometry);


  HPRefinement (*mesh, ref, levels, parameter);
}
*/

void Ng_HPRefinement (int levels, double parameter, bool setorders,
                      bool ref_level)
{
  Refinement * ref;

  if (stlgeometry)
    ref = new RefinementSTLGeometry (*stlgeometry);
  else if (geometry2d)
    ref = new Refinement2d (*geometry2d);
  else
    ref = new RefinementSurfaces (*geometry);


  HPRefinement (*mesh, ref, levels, parameter, setorders, ref_level);
}


void Ng_HighOrder (int order, bool rational)
{
  NgLock meshlock (mesh->MajorMutex(), true);

  Refinement * ref;

  if (stlgeometry)
    ref = new RefinementSTLGeometry (*stlgeometry);
#ifdef OCCGEOMETRY
  else if (occgeometry)
    ref = new OCCRefinementSurfaces (*occgeometry);
#endif
#ifdef ACIS
  else if (acisgeometry)
    {
      ref = new ACISRefinementSurfaces (*acisgeometry);
    }
#endif
  else if (geometry2d)
    ref = new Refinement2d (*geometry2d);
  else
    {
       ref = new RefinementSurfaces (*geometry);
    }
 
  // cout << "parameter 1: " << argv[1] << " (conversion to int = " << atoi(argv[1]) << ")" << endl;
 

  mesh -> GetCurvedElements().BuildCurvedElements (ref, order, rational);


  /*
  if(mesh)
    mesh -> GetCurvedElements().BuildCurvedElements (ref, order, rational);
  */

  delete ref;
}












int Ng_ME_GetNVertices (NG_ELEMENT_TYPE et)
{
  switch (et)
    {
    case NG_SEGM:
    case NG_SEGM3:
      return 2;

    case NG_TRIG:
    case NG_TRIG6:
      return 3;

    case NG_QUAD:
      return 4;

    case NG_TET:
    case NG_TET10:
      return 4;

    case NG_PYRAMID:
      return 5;

    case NG_PRISM:
   case NG_PRISM12:
      return 6;

    case NG_HEX:
      return 8;

    default:
      cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return 0;
}

int Ng_ME_GetNEdges (NG_ELEMENT_TYPE et)
{
  switch (et)
    {
    case NG_SEGM:
    case NG_SEGM3:
      return 1;

    case NG_TRIG:
    case NG_TRIG6:
      return 3;

    case NG_QUAD:
      return 4;

    case NG_TET:
    case NG_TET10:
      return 6;

    case NG_PYRAMID:
      return 8;

    case NG_PRISM:
    case NG_PRISM12:
      return 9;

    case NG_HEX:
      return 12;

    default:
      cerr << "Ng_ME_GetNEdges, illegal element type " << et << endl;
    }
  return 0;
}


int Ng_ME_GetNFaces (NG_ELEMENT_TYPE et)
{
  switch (et)
    {
    case NG_SEGM:
    case NG_SEGM3:
      return 0;

    case NG_TRIG:
    case NG_TRIG6:
      return 1;

    case NG_QUAD:
    case NG_QUAD6:
      return 1;

    case NG_TET:
    case NG_TET10:
      return 4;

    case NG_PYRAMID:
      return 5;

    case NG_PRISM:
    case NG_PRISM12:
      return 5;

    case NG_HEX:
      return 6;

    default:
      cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return 0;
}


const NG_POINT * Ng_ME_GetVertices (NG_ELEMENT_TYPE et)
{
  static double segm_points [][3] = 
    { { 1, 0, 0 },
      { 0, 0, 0 } };

  static double trig_points [][3] = 
    { { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 } };

  static double quad_points [][3] = 
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 } };

  static double tet_points [][3] = 
    { { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 0, 0, 0 } };

  static double pyramid_points [][3] =
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 },
      { 0, 0, 1-1e-7 },
    };    
  
  static double prism_points[][3] = 
    {
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 0, 0, 1 }
    };

  switch (et)
    {
    case NG_SEGM:
    case NG_SEGM3:
      return segm_points;

    case NG_TRIG:
    case NG_TRIG6:
      return trig_points;

    case NG_QUAD:
    case NG_QUAD6:
      return quad_points;

    case NG_TET:
    case NG_TET10:
      return tet_points;

    case NG_PYRAMID:
      return pyramid_points;

    case NG_PRISM:
    case NG_PRISM12:
      return prism_points;

    case NG_HEX:
    default:
      cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
    }
  return 0;
}



const NG_EDGE * Ng_ME_GetEdges (NG_ELEMENT_TYPE et)
{
  static int segm_edges[1][2] =
    { { 1, 2 }};

  static int trig_edges[3][2] =
    { { 3, 1 },
      { 3, 2 },
      { 1, 2 }};

  static int quad_edges[4][2] =
    { { 1, 2 },
      { 4, 3 },
      { 1, 4 },
      { 2, 3 }};


  static int tet_edges[6][2] =
    { { 4, 1 },
      { 4, 2 },
      { 4, 3 }, 
      { 1, 2 },
      { 1, 3 },
      { 2, 3 }};

  static int prism_edges[9][2] =
    { { 3, 1 },
      { 1, 2 },
      { 3, 2 },
      { 6, 4 },
      { 4, 5 },
      { 6, 5 },
      { 3, 6 },
      { 1, 4 },
      { 2, 5 }};

  static int pyramid_edges[8][2] =
    { { 1, 2 },
      { 2, 3 },
      { 1, 4 },
      { 4, 3 },
      { 1, 5 },
      { 2, 5 },
      { 3, 5 },
      { 4, 5 }};



  switch (et)
    {
    case NG_SEGM:
    case NG_SEGM3:
      return segm_edges;

    case NG_TRIG:
    case NG_TRIG6:
      return trig_edges;

    case NG_QUAD:
    case NG_QUAD6:
      return quad_edges;

    case NG_TET:
    case NG_TET10:
      return tet_edges;

    case NG_PYRAMID:
      return pyramid_edges;

    case NG_PRISM:
    case NG_PRISM12:
      return prism_edges;

    case NG_HEX:
    default:
      cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
    }
  return 0;  
}


const NG_FACE * Ng_ME_GetFaces (NG_ELEMENT_TYPE et)
{
  static int tet_faces[4][4] =
    { { 4, 2, 3, 0 },
      { 4, 1, 3, 0 },
      { 4, 1, 2, 0 },
      { 1, 2, 3, 0 } };
  
  static int prism_faces[5][4] =
    {
      { 1, 2, 3, 0 },
      { 4, 5, 6, 0 },
      { 3, 1, 4, 6 },
      { 1, 2, 5, 4 },
      { 2, 3, 6, 5 } 
    };
  
  static int pyramid_faces[5][4] =
    {
      { 1, 2, 5, 0 },
      { 2, 3, 5, 0 },
      { 3, 4, 5, 0 },
      { 4, 1, 5, 0 },
      { 1, 2, 3, 4 } 
    };
  
  static int trig_faces[1][4] = 
    {
      { 1, 2, 3, 0 },
    };

  switch (et)
    {
    case NG_TET:
    case NG_TET10:
      return tet_faces;

    case NG_PRISM:
    case NG_PRISM12:
      return prism_faces;

    case NG_PYRAMID:
      return pyramid_faces;


    case NG_SEGM:
    case NG_SEGM3:

    case NG_TRIG:
    case NG_TRIG6:
      return trig_faces;
    case NG_QUAD:


    case NG_HEX:

    default:
      cerr << "Ng_ME_GetFaces, illegal element type " << et << endl;
    }
  return 0;
}


int Ng_GetNEdges()
{
  return mesh->GetTopology().GetNEdges();
}
int Ng_GetNFaces()
{
  return mesh->GetTopology().GetNFaces();
}



int Ng_GetElement_Edges (int elnr, int * edges, int * orient)
{
  const MeshTopology & topology = mesh->GetTopology();
  if (mesh->GetDimension() == 3)
    return topology.GetElementEdges (elnr, edges, orient);
  else
    return topology.GetSurfaceElementEdges (elnr, edges, orient);
}

int Ng_GetElement_Faces (int elnr, int * faces, int * orient)
{
  const MeshTopology & topology = mesh->GetTopology();
  if (mesh->GetDimension() == 3)
    return topology.GetElementFaces (elnr, faces, orient);
  else
    {
      faces[0] = elnr;
      if (orient) orient[0] = 0;
      return 1;
    }
}

int Ng_GetSurfaceElement_Edges (int elnr, int * edges, int * orient)
{
  const MeshTopology & topology = mesh->GetTopology();
  if (mesh->GetDimension() == 3)
    return topology.GetSurfaceElementEdges (elnr, edges, orient);
  else
    {
      if (orient)
	topology.GetSegmentEdge(elnr, edges[0], orient[0]);
      else
	edges[0] = topology.GetSegmentEdge(elnr);
    }
  return 1;
  /*
    int i, ned;
    const MeshTopology & topology = mesh->GetTopology();
    ARRAY<int> ia;
    topology.GetSurfaceElementEdges (elnr, ia);
    ned = ia.Size();
    for (i = 1; i <= ned; i++)
    edges[i-1] = ia.Get(i);

    if (orient)
    {
    topology.GetSurfaceElementEdgeOrientations (elnr, ia);
    for (i = 1; i <= ned; i++)
    orient[i-1] = ia.Get(i);
    }
    return ned;
  */
}

int Ng_GetSurfaceElement_Face (int selnr, int * orient)
{
  if (mesh->GetDimension() == 3)
    {
      const MeshTopology & topology = mesh->GetTopology();
      if (orient)
	*orient = topology.GetSurfaceElementFaceOrientation (selnr);
      return topology.GetSurfaceElementFace (selnr);
    }
  return -1;
}

int Ng_GetFace_Vertices (int fnr, int * vert)
{
  const MeshTopology & topology = mesh->GetTopology();
  ArrayMem<int,4> ia;
  topology.GetFaceVertices (fnr, ia);
  for (int i = 0; i < ia.Size(); i++)
    vert[i] = ia[i];
  //  cout << "face verts = " << ia << endl;
  return ia.Size();
}


int Ng_GetFace_Edges (int fnr, int * edge)
{
  const MeshTopology & topology = mesh->GetTopology();
  ArrayMem<int,4> ia;
  topology.GetFaceEdges (fnr, ia);
  for (int i = 0; i < ia.Size(); i++)
    edge[i] = ia[i];
  return ia.Size();
}

void Ng_GetEdge_Vertices (int ednr, int * vert)
{
  const MeshTopology & topology = mesh->GetTopology();
  topology.GetEdgeVertices (ednr, vert[0], vert[1]);
}


int Ng_GetNVertexElements (int vnr)
{
  if (mesh->GetDimension() == 3)
    return mesh->GetTopology().GetVertexElements(vnr).Size();
  else
    return mesh->GetTopology().GetVertexSurfaceElements(vnr).Size();
}

void Ng_GetVertexElements (int vnr, int * els)
{
  FlatArray<int> ia(0,0);
  if (mesh->GetDimension() == 3)
    ia = mesh->GetTopology().GetVertexElements(vnr);
  else
    ia = mesh->GetTopology().GetVertexSurfaceElements(vnr);
  for (int i = 0; i < ia.Size(); i++)
    els[i] = ia[i];
}


int Ng_GetElementOrder (int enr)
{
  if (mesh->GetDimension() == 3)
    return mesh->VolumeElement(enr).GetOrder();
  else
    return mesh->SurfaceElement(enr).GetOrder();
}

void Ng_GetElementOrders (int enr, int * ox, int * oy, int * oz)
{
  if (mesh->GetDimension() == 3)
    mesh->VolumeElement(enr).GetOrder(*ox, *oy, *oz);
  else
    mesh->SurfaceElement(enr).GetOrder(*ox, *oy, *oz);
}

void Ng_SetElementOrder (int enr, int order)
{
  if (mesh->GetDimension() == 3)
    return mesh->VolumeElement(enr).SetOrder(order);
  else
    return mesh->SurfaceElement(enr).SetOrder(order);
}

void Ng_SetElementOrders (int enr, int ox, int oy, int oz)
{
  if (mesh->GetDimension() == 3)
    mesh->VolumeElement(enr).SetOrder(ox, oy, oz);
  else
    mesh->SurfaceElement(enr).SetOrder(ox, oy);
}


int Ng_GetSurfaceElementOrder (int enr)
{
  return mesh->SurfaceElement(enr).GetOrder();
}

//HERBERT: falsche Anzahl von Argumenten
//void Ng_GetSurfaceElementOrders (int enr, int * ox, int * oy, int * oz)
void Ng_GetSurfaceElementOrders (int enr, int * ox, int * oy)
{
  int d; 
  mesh->SurfaceElement(enr).GetOrder(*ox, *oy, d);
}

void Ng_SetSurfaceElementOrder (int enr, int order)
{
  return mesh->SurfaceElement(enr).SetOrder(order);
}

void Ng_SetSurfaceElementOrders (int enr, int ox, int oy)
{
  mesh->SurfaceElement(enr).SetOrder(ox, oy);
}


int Ng_GetNLevels ()
{
  return (mesh) ? mesh->mglevels : 0;
}


void Ng_GetParentNodes (int ni, int * parents)
{
  if (ni <= mesh->mlbetweennodes.Size())
    {
      parents[0] = mesh->mlbetweennodes.Get(ni).I1();
      parents[1] = mesh->mlbetweennodes.Get(ni).I2();
    }
  else
    parents[0] = parents[1] = 0;
}


int Ng_GetParentElement (int ei)
{
  if (mesh->GetDimension() == 3)
    {
      if (ei <= mesh->mlparentelement.Size())
	return mesh->mlparentelement.Get(ei);
    }
  else
    {
      if (ei <= mesh->mlparentsurfaceelement.Size())
	return mesh->mlparentsurfaceelement.Get(ei);
    }
  return 0;
}


int Ng_GetParentSElement (int ei)
{
  if (mesh->GetDimension() == 3)
    {
      if (ei <= mesh->mlparentsurfaceelement.Size())
	return mesh->mlparentsurfaceelement.Get(ei);
    }
  else
    {
      return 0;
    }
  return 0;
}





int Ng_GetClusterRepVertex (int pi)
{
  return mesh->GetClusters().GetVertexRepresentant(pi);
}

int Ng_GetClusterRepEdge (int pi)
{
  return mesh->GetClusters().GetEdgeRepresentant(pi);
}

int Ng_GetClusterRepFace (int pi)
{
  return mesh->GetClusters().GetFaceRepresentant(pi);
}

int Ng_GetClusterRepElement (int pi)
{
  return mesh->GetClusters().GetElementRepresentant(pi);
}






void Ng_InitSolutionData (Ng_SolutionData * soldata)
{
  soldata -> name = NULL;
  soldata -> data = NULL;
  soldata -> components = 1;
  soldata -> dist = 1;
  soldata -> order = 1;
  soldata -> iscomplex = 0;
  soldata -> draw_surface = 1;
  soldata -> draw_volume = 1;
  soldata -> soltype = NG_SOLUTION_NODAL;
  soldata -> solclass = 0;
}

void Ng_SetSolutionData (Ng_SolutionData * soldata)
{
#ifdef OPENGL
  //   vssolution.ClearSolutionData ();
  VisualSceneSolution::SolData * vss = new VisualSceneSolution::SolData;

  //  cout << "Add solution " << soldata->name << ", type = " << soldata->soltype << endl;

  vss->name = new char[strlen (soldata->name)+1];
  strcpy (vss->name, soldata->name);

  vss->data = soldata->data;
  vss->components = soldata->components;
  vss->dist = soldata->dist;
  vss->order = soldata->order;
  vss->iscomplex = bool(soldata->iscomplex);
  vss->draw_surface = soldata->draw_surface;
  vss->draw_volume = soldata->draw_volume;
  vss->soltype = VisualSceneSolution::SolType (soldata->soltype);
  vss->solclass = soldata->solclass;
  vssolution.AddSolutionData (vss);
#endif
}

void Ng_ClearSolutionData ()
{
#ifdef OPENGL
  vssolution.ClearSolutionData();
#endif
}



void Ng_Redraw ()
{
#ifdef OPENGL
  extern bool nodisplay; // he: global in ngappinit.cpp
  if (!nodisplay)
  {
    vssolution.UpdateSolutionTimeStamp();
    Render();
  }
#endif
}


void Ng_SetVisualizationParameter (const char * name, const char * value)
{
#ifdef OPENGL
#ifndef NOTCL
  char buf[100];
  sprintf (buf, "visoptions.%s", name);
  if (printmessage_importance>0)
  {
    cout << "name = " << name << ", value = " << value << endl;
    cout << "set tcl-variable " << buf << " to " << value << endl;
  }
  Tcl_SetVar (tcl_interp, buf, const_cast<char*> (value), 0);
  Tcl_Eval (tcl_interp, "Ng_Vis_Set parameters;");
#endif
#endif
}




int firsttime = 1;
int animcnt = 0;
void PlayAnimFile(const char* name, int speed, int maxcnt)
{
  //extern Mesh * mesh;

  /*
  if (mesh.Ptr()) mesh->DeleteMesh();
  if (!mesh.Ptr()) mesh = new Mesh();
  */
  mesh.Reset (new Mesh());

  int ne, np, i;

  char str[80];
  char str2[80];

  //int tend = 5000;
  //  for (ti = 1; ti <= tend; ti++)
  //{
  int rti = (animcnt%(maxcnt-1)) + 1;
  animcnt+=speed;
  
  sprintf(str2,"%05i.sol",rti);
  strcpy(str,"mbssol/");
  strcat(str,name);
  strcat(str,str2);

  if (printmessage_importance>0)
    cout << "read file '" << str << "'" << endl;
  
  ifstream infile(str);
  infile >> ne;
  for (i = 1; i <= ne; i++)
    {
      int j;
      Element2d tri(TRIG);
      tri.SetIndex(1); //faceind
      
      for (j = 1; j <= 3; j++)
	infile >> tri.PNum(j);

      infile >> np;
      for (i = 1; i <= np; i++)
	{
	  Point3d p;
	  infile >> p.X() >> p.Y() >> p.Z();
	  if (firsttime)
	    mesh->AddPoint (p);
	  else
	    mesh->Point(i) = Point<3> (p);
	}

      //firsttime = 0;
      Ng_Redraw();
   }
}

		
int Ng_GetNPeriodicVertices (int idnr)
{
  ARRAY<INDEX_2> apairs;
  mesh->GetIdentifications().GetPairs (idnr, apairs);
  return apairs.Size();
}


// pairs should be an integer array of 2*npairs
void Ng_GetPeriodicVertices (int idnr, int * pairs)
{
  ARRAY<INDEX_2> apairs;
  mesh->GetIdentifications().GetPairs (idnr, apairs);
  for (int i = 0; i < apairs.Size(); i++)
    {
      pairs[2*i] = apairs[i].I1();
      pairs[2*i+1] = apairs[i].I2();
    }
      
}



int Ng_GetNPeriodicEdges (int idnr)
{
  ARRAY<INDEX,PointIndex::BASE> map;
  //const MeshTopology & top = mesh->GetTopology();
  int nse = mesh->GetNSeg();

  int cnt = 0;
  //  for (int id = 1; id <= mesh->GetIdentifications().GetMaxNr(); id++)
    {
      mesh->GetIdentifications().GetMap(idnr, map);
      //(*testout) << "ident-map " << id << ":" << endl << map << endl;

      for (SegmentIndex si = 0; si < nse; si++)
	{
	  PointIndex other1 = map[(*mesh)[si].p1];
	  PointIndex other2 = map[(*mesh)[si].p2];
	  //  (*testout) << "seg = " << (*mesh)[si] << "; other = " 
	  //     << other1 << "-" << other2 << endl;
	  if (other1 && other2 && mesh->IsSegment (other1, other2))
	    {
	      cnt++;
	    }
	}
    }
  return cnt;
}

void Ng_GetPeriodicEdges (int idnr, int * pairs)
{
  ARRAY<INDEX,PointIndex::BASE> map;
  const MeshTopology & top = mesh->GetTopology();
  int nse = mesh->GetNSeg();

  int cnt = 0;
  //  for (int id = 1; id <= mesh->GetIdentifications().GetMaxNr(); id++)
    {
      mesh->GetIdentifications().GetMap(idnr, map);
      
      //(*testout) << "map = " << map << endl;

      for (SegmentIndex si = 0; si < nse; si++)
	{
	  PointIndex other1 = map[(*mesh)[si].p1];
	  PointIndex other2 = map[(*mesh)[si].p2];
	  if (other1 && other2 && mesh->IsSegment (other1, other2))
	    {
	      SegmentIndex otherseg = mesh->SegmentNr (other1, other2);
	      pairs[cnt++] = top.GetSegmentEdge (si+1);
	      pairs[cnt++] = top.GetSegmentEdge (otherseg+1);
	    }
	}
    }
}



void Ng_PushStatus (const char * str)
{
  PushStatus (MyStr (str));
}

void Ng_PopStatus ()
{
  PopStatus ();
}

void Ng_SetThreadPercentage (double percent)
{
  SetThreadPercent (percent);
}

void Ng_GetStatus (char ** str, double & percent)
{
  MyStr s;
  GetStatus(s,percent);
  *str = new char[s.Length()+1];
  strcpy(*str,s.c_str());  
}


void Ng_SetTerminate(void)
{
  multithread.terminate = 1;
}
void Ng_UnSetTerminate(void)
{
  multithread.terminate = 0;
}

int Ng_ShouldTerminate(void)
{
  return multithread.terminate;
}

///// Added by Roman Stainko ....
int Ng_GetVertex_Elements( int vnr, int* elems )
{
  const MeshTopology& topology = mesh->GetTopology();
  ArrayMem<int,4> indexArray;
  topology.GetVertexElements( vnr, indexArray );
  
  for( int i=0; i<indexArray.Size(); i++ )
    elems[i] = indexArray[i];
  
  return indexArray.Size();
}

///// Added by Roman Stainko ....
int Ng_GetVertex_SurfaceElements( int vnr, int* elems )
{
  const MeshTopology& topology = mesh->GetTopology();
  ArrayMem<int,4> indexArray;
  topology.GetVertexSurfaceElements( vnr, indexArray );
  
  for( int i=0; i<indexArray.Size(); i++ )
    elems[i] = indexArray[i];
  
  return indexArray.Size();
}

///// Added by Roman Stainko ....
int Ng_GetVertex_NElements( int vnr )
{
  const MeshTopology& topology = mesh->GetTopology();
  ArrayMem<int,4> indexArray;
  topology.GetVertexElements( vnr, indexArray );
  
  return indexArray.Size();
}

///// Added by Roman Stainko ....
int Ng_GetVertex_NSurfaceElements( int vnr )
{
  const MeshTopology& topology = mesh->GetTopology();
  ArrayMem<int,4> indexArray;
  topology.GetVertexSurfaceElements( vnr, indexArray );

  return indexArray.Size();
}



#ifdef SOCKETS
int Ng_SocketClientOpen( const int port, const char * host )
{
  try
    {
      if(host)
	clientsocket.Reset(new ClientSocket(port,host));
      else
	clientsocket.Reset(new ClientSocket(port));
    }
  catch( SocketException e)
    {
      cerr << e.Description() << endl;
      return 0;
    }
  return 1;
}
 
void Ng_SocketClientWrite( const char * write, char** reply)
{
  string output = write;
  (*clientsocket) << output;
  string sreply;
  (*clientsocket) >> sreply;
  *reply = new char[sreply.size()+1];
  strcpy(*reply,sreply.c_str());
}


void Ng_SocketClientClose ( void )
{
  clientsocket.Reset(NULL);
}


void Ng_SocketClientGetServerHost ( const int number, char ** host )
{
  *host = new char[servers[number]->host.size()+1];
  strcpy(*host,servers[number]->host.c_str());
}

void Ng_SocketClientGetServerPort ( const int number, int * port )
{
  *port = servers[number]->port;
}

void Ng_SocketClientGetServerClientID ( const int number, int * id )
{
  *id = servers[number]->clientid;
}

#endif // SOCKETS




#ifdef PARALLEL
void Ng_SetElementPartition ( const int elnr, const int part )
{
  mesh->VolumeElement(elnr+1).SetPartition(part);

}
int Ng_GetElementPartition ( const int elnr )
{
  return mesh->VolumeElement(elnr+1).GetPartition();
}
#endif


void Ng_InitPointCurve(double red, double green, double blue)
{
  mesh->InitPointCurve(red, green, blue);
}

void Ng_AddPointCurvePoint(const double * point)
{
  Point3d pt;
  pt.X() = point[0];
  pt.Y() = point[1];
  pt.Z() = point[2];
  mesh->AddPointCurvePoint(pt);
}


void Ng_SaveMesh ( const char * meshfile )
{
  mesh -> Save(string(meshfile));
}


int Ng_Bisect_WithInfo ( const char * refinementfile, double ** qualityloss, int * qualityloss_size )
{
  BisectionOptions biopt;
  biopt.outfilename = NULL; // "ngfepp.vol";
  biopt.femcode = "fepp";
  biopt.refinementfilename = refinementfile;
  
  Refinement * ref;
  MeshOptimize2d * opt = NULL;
  if (stlgeometry)
    ref = new RefinementSTLGeometry(*stlgeometry);
#ifdef OCCGEOMETRY
  else if (occgeometry)
    ref = new OCCRefinementSurfaces (*occgeometry);
#endif
#ifdef ACIS
  else if (acisgeometry)
    {
      ref = new ACISRefinementSurfaces(*acisgeometry);
      opt = new ACISMeshOptimize2dSurfaces(*acisgeometry);
      ref->Set2dOptimizer(opt);
    }
#endif
  else
    {
      ref = new RefinementSurfaces(*geometry);
      opt = new MeshOptimize2dSurfaces(*geometry);
      ref->Set2dOptimizer(opt);
    }
  
  if(!mesh->LocalHFunctionGenerated())
    mesh->CalcLocalH();
  
  mesh->LocalHFunction().SetGrading (mparam.grading);

  ARRAY<double> * qualityloss_arr = NULL;
  if(qualityloss != NULL)
    qualityloss_arr = new ARRAY<double>;

  ref -> Bisect (*mesh, biopt, qualityloss_arr);

  int retval = 0;

  if(qualityloss != NULL)
    {
      *qualityloss = new double[qualityloss_arr->Size()+1];

      for(int i = 0; i<qualityloss_arr->Size(); i++)
	(*qualityloss)[i+1] = (*qualityloss_arr)[i];

      retval = qualityloss_arr->Size();

      delete qualityloss_arr;
    }

  mesh -> UpdateTopology();
  mesh -> GetCurvedElements().BuildCurvedElements (ref, mparam.elementorder);
  
  multithread.running = 0;
  delete ref;
  delete opt;

  return retval;
}

void Ng_Bisect ( const char * refinementfile )
{
  Ng_Bisect_WithInfo( refinementfile, NULL, NULL );
}





/*
  number of nodes of type nt
  nt = 0 is Vertex
  nt = 1 is Edge
  nt = 2 is Face
  nt = 3 is Cell
*/
int Ng_GetNNodes (int nt)
{
  switch (nt)
    {
    case 0: return mesh -> GetNV();
    case 1: return mesh->GetTopology().GetNEdges();
    case 2: return mesh->GetTopology().GetNFaces();
    case 3: return mesh -> GetNE();
    }
  return -1;
}


int Ng_GetClosureNodes (int nt, int nodenr, int nodeset, int * nodes)
{
  switch (nt)
    {
    case 3:  // The closure of a cell
      {
        int cnt = 0;
        if (nodeset & 1)  // Vertices
          {
            const Element & el = (*mesh)[ElementIndex(nodenr)];
            for (int i = 0; i < el.GetNP(); i++)
              { 
                nodes[cnt++] = 0;
                nodes[cnt++] = el[i] - PointIndex::BASE;
              }
          }

        if (nodeset & 2)  // Edges
          {
            int edges[12];
            int ned;
            ned = mesh->GetTopology().GetElementEdges (nodenr+1, edges, 0);
            for (int i = 0; i < ned; i++)
              {
                nodes[cnt++] = 1;
                nodes[cnt++] = edges[i]-1;
              }
          }

        if (nodeset & 4)  // Faces
          {
            int faces[12];
            int nfa;
            nfa = mesh->GetTopology().GetElementFaces (nodenr+1, faces, 0);
            for (int i = 0; i < nfa; i++)
              {
                nodes[cnt++] = 2;
                nodes[cnt++] = faces[i]-1;
              }
          }

        if (nodeset & 8)  // Cell
          {
            nodes[cnt++] = 3;
            nodes[cnt++] = nodenr;
          }

        return cnt/2;
      }
    default:
      {
        cerr << "GetClosureNodes not implemented for Nodetype " << nt << endl;
      }
    }
  return 0;
}



int Ng_GetNElements (int dim)
{
  switch (dim)
    {
    case 0: return mesh -> GetNV();
    case 1: return mesh -> GetNSeg();
    case 2: return mesh -> GetNSE();
    case 3: return mesh -> GetNE();
    }
  return -1;
}



  /*
    closure nodes of element
    nodeset is bit-coded, bit 0 includes Vertices, bit 1 edges, etc
    E.g., nodeset = 6 includes edge and face nodes
    nodes is pair of integers (nodetype, nodenr) 
    return value is number of nodes
   */
int Ng_GetElementClosureNodes (int dim, int elementnr, int nodeset, int * nodes)
{
  switch (dim)
    {
    case 3:  // The closure of a volume element = CELL
      {
        return Ng_GetClosureNodes (3, elementnr, nodeset, nodes);
      }
    case 2:
      {
        int cnt = 0;
        if (nodeset & 1)  // Vertices
          {
            const Element2d & el = (*mesh)[SurfaceElementIndex(elementnr)];
            for (int i = 0; i < el.GetNP(); i++)
              { 
                nodes[cnt++] = 0;
                nodes[cnt++] = el[i] - PointIndex::BASE;
              }
          }

        if (nodeset & 2)  // Edges
          {
            int edges[12];
            int ned;
            ned = mesh->GetTopology().GetSurfaceElementEdges (elementnr+1, edges, 0);
            for (int i = 0; i < ned; i++)
              {
                nodes[cnt++] = 1;
                nodes[cnt++] = edges[i]-1;
              }
          }

        if (nodeset & 4)  // Faces
          {
            int face = mesh->GetTopology().GetSurfaceElementFace (elementnr+1);
            nodes[cnt++] = 2;
            nodes[cnt++] = face-1;
          }

        return cnt/2;
      }
    default:
      {
        cerr << "GetClosureNodes not implemented for Element of dimension " << dim << endl;
      }
    }
  return 0;
}
