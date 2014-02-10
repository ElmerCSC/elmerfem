//
//  Write user dependent output file
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "writeuser.hpp"


void RegisterUserFormats (ARRAY<const char*> & names)
{
  const char *types[] =
    {
      "Neutral Format",
      "Surface Mesh Format" ,
      "DIFFPACK Format",
      "TecPlot Format",     
      "Tochnog Format",
      "Abaqus Format",
      "Fluent Format",
      "Permas Format",
      "FEAP Format",
      "Elmer Format",
      "STL Format",
      "VRML Format",
      "Gmsh Format",
      "JCMwave Format",
      "TET Format",
      //      { "Chemnitz Format" },
      0
    };

  for (int i = 0; types[i]; i++)
    names.Append (types[i]);
}



bool WriteUserFormat (const string & format,
		      const Mesh & mesh,
		      const CSGeometry & geom, 
		      const string & filename)
{
  PrintMessage (1, "Export mesh to file ", filename, 
		", format is ", format);

  if (format == "Neutral Format")
    WriteNeutralFormat (mesh, geom, filename);

  else if (format == "Surface Mesh Format")
    WriteSurfaceFormat (mesh, filename);

  else if (format == "DIFFPACK Format")
    WriteDiffPackFormat (mesh, geom, filename);

  else if (format == "Tochnog Format")
    WriteTochnogFormat (mesh, filename);

  else if (format == "TecPlot Format")
    cerr << "ERROR: TecPlot format currently out of order" << endl;
      // WriteTecPlotFormat (mesh, geom, filename);

  else if (format == "Abaqus Format")
    WriteAbaqusFormat (mesh, filename);

  else if (format == "Fluent Format")
    WriteFluentFormat (mesh, filename);

  else if (format == "Permas Format")
    WritePermasFormat (mesh, filename);

  else if (format == "FEAP Format")
    WriteFEAPFormat (mesh, filename);

  else if (format == "Elmer Format")
    WriteElmerFormat (mesh, filename);

  else if (format == "STL Format")
    WriteSTLFormat (mesh, filename);

  else if (format == "VRML Format")
    WriteVRMLFormat (mesh, 1, filename);

  else if (format == "Fepp Format")
    WriteFEPPFormat (mesh, geom, filename);

  else if (format ==  "EdgeElement Format")
    WriteEdgeElementFormat (mesh, geom, filename);

  else if (format == "Chemnitz Format")
    WriteUserChemnitz (mesh, filename);

  else if (format == "Gmsh Format")
    WriteGmshFormat (mesh, geom, filename);
 
  else if (format == "JCMwave Format")
    WriteJCMFormat (mesh, geom, filename);

#ifdef OLIVER
  else if (format == "TET Format")
    WriteTETFormat( mesh, filename);//, "High Frequency" );
#endif

  else 
    {
      return 1;
    }

  return 0;
}




/*
 *  Neutral mesh format
 *  points, elements, surface elements
 */

void WriteNeutralFormat (const Mesh & mesh,
			 const CSGeometry & geom,
			 const string & filename)
{
  cout << "write neutral, new" << endl;
  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int nseg = mesh.GetNSeg();
  int i, j;
  
  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;

  ofstream outfile (filename.c_str());

  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);
  
  outfile << np << "\n";
  
  for (i = 1; i <= np; i++)
    {
      const Point3d & p = mesh.Point(i);
      
      outfile.width(10);
      outfile << p.X() << " ";
      outfile.width(9);
      outfile << p.Y() << " ";
      if (mesh.GetDimension() == 3)
	{
	  outfile.width(9);
	  outfile << p.Z();
	  }
      outfile << "\n";
    }

  if (mesh.GetDimension() == 3)
    {
      outfile << ne << "\n";
      for (i = 1; i <= ne; i++)
	{
	  Element el = mesh.VolumeElement(i);
	  if (inverttets)
	    el.Invert();
	  outfile.width(4);
	  outfile << el.GetIndex() << "  ";
	  for (j = 1; j <= el.GetNP(); j++)
	    {
	      outfile << " ";
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << "\n";
	}
    }

  outfile << nse << "\n";
  for (i = 1; i <= nse; i++)
    {
      Element2d el = mesh.SurfaceElement(i);
      if (invertsurf)
	el.Invert();
      outfile.width(4);
      outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << "    ";
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << el.PNum(j);
	}
      outfile << "\n";
    }


  if (mesh.GetDimension() == 2)
    {
      outfile << nseg << "\n";
      for (i = 1; i <= nseg; i++)
	{
	  const Segment & seg = mesh.LineSegment(i);
	  outfile.width(4);
	  outfile << seg.si << "    ";

	  outfile << " ";
	  outfile.width(8);
	  outfile << seg.p1;
	  outfile << " ";
	  outfile.width(8);
	  outfile << seg.p2;

	  outfile << "\n";
	}
    }
}









void WriteSurfaceFormat (const Mesh & mesh,
			 const string & filename)
{
  // surface mesh
  int i, j;
  
  cout << "Write Surface Mesh" << endl;
  
  ofstream outfile (filename.c_str());
  
  outfile << "surfacemesh" << endl;
  
  outfile << mesh.GetNP() << endl;
  for (i = 1; i <= mesh.GetNP(); i++)
    {
      for (j = 0; j < 3; j++)
	{
	  outfile.width(10);
	  outfile << mesh.Point(i)(j) << " ";
	}
      outfile << endl;
    }
  outfile << mesh.GetNSE() << endl;
  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      for (j = 1; j <= 3; j++)
	{
	  outfile.width(8);
	  outfile << mesh.SurfaceElement(i).PNum(j);
	}
      outfile << endl;
    }
}





/*
 *  save surface mesh as STL file
 */

void WriteSTLFormat (const Mesh & mesh,
		     const string & filename)
{
  cout << "\nWrite STL Surface Mesh" << endl;
  
  ofstream outfile (filename.c_str());
  
  int i;
  
  outfile.precision(10);
  
  outfile << "solid" << endl;
  
  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      outfile << "facet normal ";
      const Point3d& p1 = mesh.Point(mesh.SurfaceElement(i).PNum(1));
      const Point3d& p2 = mesh.Point(mesh.SurfaceElement(i).PNum(2));
      const Point3d& p3 = mesh.Point(mesh.SurfaceElement(i).PNum(3));
      
      Vec3d normal = Cross(p2-p1,p3-p1);
      if (normal.Length() != 0)
	{
	  normal /= (normal.Length());		  
	}
      
      outfile << normal.X() << " " << normal.Y() << " " << normal.Z() << "\n";
      outfile << "outer loop\n";
      
      outfile << "vertex " << p1.X() << " " << p1.Y() << " " << p1.Z() << "\n";
      outfile << "vertex " << p2.X() << " " << p2.Y() << " " << p2.Z() << "\n";
      outfile << "vertex " << p3.X() << " " << p3.Y() << " " << p3.Z() << "\n";
      
      outfile << "endloop\n";
      outfile << "endfacet\n"; 
    }
  outfile << "endsolid" << endl;
}





/*
 *
 *  write surface mesh as VRML file
 *
 */

void WriteVRMLFormat (const Mesh & mesh,
		      bool faces,
		      const string & filename)
{

  if (faces)

    {
      // Output in VRML, IndexedFaceSet is used 
      // Bartosz Sawicki <sawickib@ee.pw.edu.pl>

      int np = mesh.GetNP();
      int nse = mesh.GetNSE();
      int i, j;

      ofstream outfile (filename.c_str());

      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      outfile << "#VRML V2.0 utf8 \n"
	         "Background {\n"
		 "    skyColor [1 1 1]\n"
     		 "    groundColor [1 1 1]\n"
		 "}\n"
		 "Group{ children [\n"
		 "Shape{ \n"
		 "appearance Appearance { material Material { }} \n"
                 "geometry IndexedFaceSet { \n"
                 "coord Coordinate { point [ \n";  
	        

      for (i = 1; i <= np; i++)
        {
          const Point3d & p = mesh.Point(i);
          outfile.width(10);
          outfile << p.X() << " ";
          outfile << p.Y() << " ";
          outfile << p.Z() << " \n";
	}

      outfile << "  ] } \n"
                 "coordIndex [ \n";               
	
      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);
      
	  for (j = 1; j <= 3; j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j)-1;
	    }
	  outfile << " -1 \n";
	}
      
      outfile << "  ] \n";

      //define number and RGB definitions of colors
      outfile << "color Color { color [1 0 0, 0 1 0, 0 0 1, 1 1 0]} \n"
                 "colorIndex [\n";
      
      for (i = 1; i <= nse; i++)
	{
	  outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex ()).BCProperty();      
          outfile << endl;
	}
      
      outfile << " ] \n"
                 "colorPerVertex FALSE \n"
                 "creaseAngle 0 \n"
		 "solid FALSE \n"
                 "ccw FALSE \n"
		 "convex TRUE \n"
                 "} } # end of Shape\n"
		 "] }\n";
                         
    } /* end of VRMLFACES */


  else

    {
        // Output in VRML, IndexedLineSet is used
	// Bartosz Sawicki <sawickib@ee.pw.edu.pl>

      int np = mesh.GetNP();
      int nse = mesh.GetNSE();
      int i, j;

      ofstream outfile (filename.c_str());

      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      outfile << "#VRML V2.0 utf8 \n"
	         "Background {\n"
		 "    skyColor [1 1 1]\n"
     		 "    groundColor [1 1 1]\n"
		 "}\n"
		 "Group{ children [\n"
	         "Shape{ \n"
		 "appearance Appearance { material Material { }} \n"
                 "geometry IndexedLineSet { \n"
                 "coord Coordinate { point [ \n";  
	        

      for (i = 1; i <= np; i++)
        {
          const Point3d & p = mesh.Point(i);
          outfile.width(10);
          outfile << p.X() << " ";
          outfile << p.Y() << " ";
          outfile << p.Z() << " \n";
	}

      outfile << "  ] } \n"
                 "coordIndex [ \n";               
	
      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);
      
	  for (j = 1; j <= 3; j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j)-1;
	    }
	  outfile.width(8);  
	  outfile << el.PNum(1)-1; 
	  outfile << " -1 \n";
	}
      
      outfile << "  ] \n";

/* Uncomment if you want color mesh    
      outfile << "color Color { color [1 1 1, 0 1 0, 0 0 1, 1 1 0]} \n"
                 "colorIndex [\n";
      
      for (i = 1; i <= nse; i++)
	{
	  outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex ()).BCProperty();      
          outfile << endl;
	}
      
      outfile << " ] \n"
*/ 
      outfile << "colorPerVertex FALSE \n"
                 "} } #end of Shape\n"
		 "] } \n";
                         
    }

}






/*
 * FEPP .. a finite element package developed at University Linz, Austria
 */
void WriteFEPPFormat (const Mesh & mesh,
		      const CSGeometry & geom,
		      const string & filename)
{
  
  ofstream outfile (filename.c_str());

  if (mesh.GetDimension() == 3)

    {

      // output for FEPP
      
      int np = mesh.GetNP();
      int ne = mesh.GetNE();
      int nse = mesh.GetNSE();
      int ns = mesh.GetNFD();
      int i, j;

      outfile.precision(5);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);
      
      outfile << "volumemesh4" << endl;
      outfile << nse << endl;
      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);

	  //	  int facenr = mesh.facedecoding.Get(el.GetIndex()).surfnr;
	  outfile.width(4);
	  outfile << el.GetIndex() << " ";
	  outfile.width(4);
	  //	  outfile << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
	  outfile << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
	  outfile.width(4);
	  outfile << el.GetNP() << "    ";
	  for (j = 1; j <= el.GetNP(); j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << "\n";
	}


      outfile << ne << "\n";
      for (i = 1; i <= ne; i++)
	{
	  const Element & el = mesh.VolumeElement(i);
	  outfile.width(4);
	  outfile << el.GetIndex() << " ";
	  outfile.width(4);
	  outfile << el.GetNP() << " ";
	  for (j = 1; j <= el.GetNP(); j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << "\n";
	}

      outfile << np << "\n";
      for (i = 1; i <= np; i++)
	{
	  const Point3d & p = mesh.Point(i);

	  outfile.width(10);
	  outfile << p.X() << " ";
	  outfile.width(9);
	  outfile << p.Y() << " ";
	  outfile.width(9);
	  outfile << p.Z() << "\n";
	}

      /*      
      if (typ == WRITE_FEPPML)
	{
	  int nbn =  mesh.mlbetweennodes.Size();
	  outfile << nbn << "\n";
	  for (i = 1; i <= nbn; i++)
	    outfile << mesh.mlbetweennodes.Get(i).I1() << " "
		    << mesh.mlbetweennodes.Get(i).I2() << "\n";
	  

	  //	  int ncon = mesh.connectedtonode.Size();
	  //	  outfile << ncon << "\n";
	  //	  for (i = 1; i <= ncon; i++)
	  //	    outfile << i << " " << mesh.connectedtonode.Get(i) << endl;
	}
      */


      // write CSG surfaces
      if (&geom && geom.GetNSurf() >= ns)
	{
	  outfile << ns << endl;
	  for (i = 1; i <= ns; i++)
	    geom.GetSurface(mesh.GetFaceDescriptor(i).SurfNr())->Print(outfile);
	}
      else 
	outfile << "0" << endl;
    }

  
  else
    
    { // 2D fepp format
      
      ;
      /*
      extern SplineGeometry2d * geometry2d;
      if (geometry2d)
	Save2DMesh (mesh, &geometry2d->GetSplines(), outfile);
      else
	Save2DMesh (mesh, 0, outfile);
      */
    }
}






/*
 *  Edge element mesh format
 *  points, elements, edges
 */

void WriteEdgeElementFormat (const Mesh & mesh,
			     const CSGeometry & geom,
			     const string & filename)
{
  cout << "write edge element format" << endl;

  const MeshTopology * top = &mesh.GetTopology();
  int npoints = mesh.GetNP();
  int nelements = mesh.GetNE();
  int nsurfelem = mesh.GetNSE();
  int nedges = top->GetNEdges();
  int i, j;
  
  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;
  ARRAY<int> edges;

  ofstream outfile (filename.c_str());

  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);


  // vertices with coordinates  
  outfile << npoints << "\n";
  for (i = 1; i <= npoints; i++)
    {
      const Point3d & p = mesh.Point(i);
      
      outfile.width(10);
      outfile << p.X() << " ";
      outfile.width(9);
      outfile << p.Y() << " ";
      outfile.width(9);
      outfile << p.Z() << "\n";
    }

  // element - edge - list
  outfile << nelements << " " << nedges << "\n";
  for (i = 1; i <= nelements; i++)
    {
      Element el = mesh.VolumeElement(i);
      if (inverttets)
      	el.Invert();
      outfile.width(4);
      outfile << el.GetIndex() << "  ";
      outfile.width(8);
      outfile << el.GetNP();
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << el.PNum(j);
	}

      top->GetElementEdges(i,edges);
      outfile << endl << "      ";
      outfile.width(8);
      outfile << edges.Size();
      for (j=1; j <= edges.Size(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << edges[j-1];
	}
      outfile << "\n";

      // orientation:
      top->GetElementEdgeOrientations(i,edges);
      outfile << "              ";
      for (j=1; j <= edges.Size(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << edges[j-1];
	}
      outfile << "\n";
    }

  // surface element - edge - list (with boundary conditions)
  outfile << nsurfelem << "\n";
  for (i = 1; i <= nsurfelem; i++)
    {
      Element2d el = mesh.SurfaceElement(i);
      if (invertsurf)
	el.Invert();
      outfile.width(4);
      outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << "  ";
      outfile.width(8);
      outfile << el.GetNP();
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << el.PNum(j);
	}

      top->GetSurfaceElementEdges(i,edges);
      outfile << endl << "      ";
      outfile.width(8);
      outfile << edges.Size();
      for (j=1; j <= edges.Size(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << edges[j-1];
	}
      outfile << "\n";
    }


  int v1, v2;
  // edge - vertex - list
  outfile << nedges << "\n";
  for (i=1; i <= nedges; i++)
    {
      top->GetEdgeVertices(i,v1,v2);
      outfile.width(4);
      outfile << v1;
      outfile << " ";
      outfile.width(8);
      outfile << v2 << endl;
    }
}









#ifdef OLDSTYLE_WRITE


void WriteFile (int typ,
		const Mesh & mesh,
		const CSGeometry & geom,
		const char * filename,
		const char * geomfile,
		double h)
{

  
  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;








  if (typ == WRITE_EDGEELEMENT)
    {
      // write edge element file
      // Peter Harscher, ETHZ

      cout << "Write Edge-Element Format" << endl;

      ofstream outfile (filename);

      int i, j;
      int ned;

      // hash table representing edges;
      INDEX_2_HASHTABLE<int> edgeht(mesh.GetNP());

      // list of edges
      ARRAY<INDEX_2> edgelist;

      // edge (point) on boundary ?
      BitArray bedge, bpoint(mesh.GetNP());
      
      static int eledges[6][2] = { { 1, 2 } , { 1, 3 } , { 1, 4 },
				   { 2, 3 } , { 2, 4 } , { 3, 4 } };

      // fill hashtable   (point1, point2)  ---->  edgenr
      for (i = 1; i <= mesh.GetNE(); i++)
	{
	  const Element & el = mesh.VolumeElement (i);
	  INDEX_2 edge;
	  for (j = 1; j <= 6; j++)
	    {
	      edge.I1() = el.PNum (eledges[j-1][0]);
	      edge.I2() = el.PNum (eledges[j-1][1]);
	      edge.Sort();

	      if (!edgeht.Used (edge))
		{
		  edgelist.Append (edge);
		  edgeht.Set (edge, edgelist.Size());
		}
	    }
	}

      
      // set bedges, bpoints
      bedge.SetSize (edgelist.Size());
      bedge.Clear();
      bpoint.Clear();

      for (i = 1; i <= mesh.GetNSE(); i++)
	{
	  const Element2d & sel = mesh.SurfaceElement(i);
	  for (j = 1; j <= 3; j++)
	    {
	      bpoint.Set (sel.PNum(j));

	      INDEX_2 edge;
	      edge.I1() = sel.PNum(j);
	      edge.I2() = sel.PNum(j%3+1);
	      edge.Sort();

	      bedge.Set (edgeht.Get (edge));
	    }
	}



      outfile << mesh.GetNE() << endl;
      // write element ---> point
      for (i = 1; i <= mesh.GetNE(); i++)
	{
	  const Element & el = mesh.VolumeElement(i);
	  
	  outfile.width(8);
	  outfile << i;
	  for (j = 1; j <= 4; j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << endl;
	}

      // write element ---> edge
      for (i = 1; i <= mesh.GetNE(); i++)
	{
	  const Element & el = mesh.VolumeElement (i);
	  INDEX_2 edge;
	  for (j = 1; j <= 6; j++)
	    {
	      edge.I1() = el.PNum (eledges[j-1][0]);
	      edge.I2() = el.PNum (eledges[j-1][1]);
	      edge.Sort();

	      outfile.width(8);
	      outfile << edgeht.Get (edge);
	    }
	  outfile << endl;
	}

      // write points
      outfile << mesh.GetNP() << endl;
      outfile.precision (6);
      for (i = 1; i <= mesh.GetNP(); i++)
	{
	  const Point3d & p = mesh.Point(i);
	  
	  for (j = 1; j <= 3; j++)
	    {
	      outfile.width(8);
	      outfile << p.X(j);
	    }
	  outfile << "       "
		  << (bpoint.Test(i) ? "1" : 0) << endl;
	}

      // write edges
      outfile << edgelist.Size() << endl;
      for (i = 1; i <= edgelist.Size(); i++)
	{
	  outfile.width(8);
	  outfile << edgelist.Get(i).I1();
	  outfile.width(8);
	  outfile << edgelist.Get(i).I2();
	  outfile << "       "
		  << (bedge.Test(i) ? "1" : "0") << endl;
	}
    }




}
#endif
}

