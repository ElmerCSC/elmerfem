//
//  Write Fluent file
//  Johannes Gerstmayr, University Linz
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{

#include "writeuser.hpp"



void WriteFluentFormat (const Mesh & mesh,
			const string & filename)

{
  cout << "start writing fluent export" << endl;
      
  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int i, j;

  ofstream outfile (filename.c_str());
  char str[100];

  outfile.precision(6);
  //outfile.setf (ios::fixed, ios::floatfield);
  //outfile.setf (ios::showpoint);
      
  outfile << "(0 \"Exported file from NETGEN \")" << endl;
  outfile << "(0 \"Dimension:\")" << endl;
  outfile << "(2 3)" << endl << endl;

  outfile << "(0 \"Nodes:\")" << endl;

  //number of nodes:
  sprintf(str,"(10 (0 1 %x 1))",np); //hexadecimal!!!
  outfile << str << endl;

  //nodes of zone 1:
  sprintf(str,"(10 (7 1 %x 1)(",np); //hexadecimal!!!
  outfile << str << endl;
  for (i = 1; i <= np; i++)
    {
      const Point3d & p = mesh.Point(i);

      //outfile.width(10);
      outfile << p.X() << " ";
      outfile << p.Y() << " ";
      outfile << p.Z() << "\n";
    }
  outfile << "))" << endl << endl;

  //write faces with elements

  outfile << "(0 \"Faces:\")" << endl;

  Element2d face, face2;
  int i2, j2;
  ARRAY<INDEX_3> surfaceelp;
  ARRAY<int> surfaceeli;
  ARRAY<int> locels;

  //no cells=no tets
  //no faces=2*tets

  int noverbface = 2*ne-nse/2;
      
  sprintf(str,"(13 (0 1 %x 0))",(noverbface+nse)); //hexadecimal!!!
  outfile << str << endl;
      
  sprintf(str,"(13 (4 1 %x 2 3)(",noverbface); //hexadecimal!!!
  outfile << str << endl;

  const_cast<Mesh&> (mesh).BuildElementSearchTree();

  for (i = 1; i <= ne; i++)
    {
      if (ne > 2000)
	{
	  if (i%2000 == 0)
	    {
	      cout << (double)i/(double)ne*100. << "%" << endl;
	    }
	}

      Element el = mesh.VolumeElement(i);
      //if (inverttets)
      //  el.Invert();
	  
      //outfile << el.GetIndex() << "    ";
      if (el.GetNP() != 4) {cout << "only tet-meshes supported in write fluent!" << endl;}
	  
      //faces:
	  
      Box3d box;
      el.GetBox(mesh.Points(), box);
      box.IncreaseRel(1e-6);

      mesh.GetIntersectingVolEls(box.PMin(),box.PMax(),locels);
      int nel = locels.Size();
      int locind;

      //cout << "nel=" << nel << endl;

      for (j = 1; j <= el.GetNFaces(); j++)
	{
	  el.GetFace(j, face);
	  face.Invert();
	  int eli2 = 0;
	  int stopsig = 0;
	      
	  for (i2 = 1; i2 <= nel; i2++)
	    {
	      locind = locels.Get(i2);
	      //cout << "  locind=" << locind << endl;

	      Element el2 = mesh.VolumeElement(locind);
	      //if (inverttets)
	      //  el2.Invert();

	      for (j2 = 1; j2 <= el2.GetNFaces(); j2++)
		{
		  el2.GetFace(j2, face2);

		  if (face2.HasFace(face)) {eli2 = locind; stopsig = 1; break;}
		}
	      if (stopsig) break;
	    }
	      
	  if (eli2==i) cout << "error in WRITE_FLUENT!!!" << endl;
	      
	  if (eli2 > i) //dont write faces two times!
	    {
	      //i: left cell, eli: right cell
	      outfile << hex << face.PNum(2) << " "
		<< hex << face.PNum(1) << " "
		<< hex << face.PNum(3) << " "
		<< hex << i  << " "
		<< hex << eli2 << "\n";
	    }
	  if (eli2 == 0) 
	    {
	      surfaceelp.Append(INDEX_3(face.PNum(2),face.PNum(1),face.PNum(3)));
	      surfaceeli.Append(i);
	    }
	}
    }
  outfile << "))" << endl;
      
  sprintf(str,"(13 (2 %x %x 3 3)(",(noverbface+1),noverbface+nse); //hexadecimal!!!
  outfile << str << endl;

  for (i = 1; i <= surfaceelp.Size(); i++)
    {
      outfile << hex << surfaceelp.Get(i).I1() << " "
	      << hex << surfaceelp.Get(i).I2() << " "
	      << hex << surfaceelp.Get(i).I3() << " "
	      << hex << surfaceeli.Get(i) << " " << 0 << "\n";
    }

  outfile << "))" << endl << endl;

  outfile << "(0 \"Cells:\")" << endl;
      
  sprintf(str,"(12 (0 1 %x 0))",ne); //hexadecimal!!!
  outfile << str << endl;

  sprintf(str,"(12 (1 1 %x 1 2))",ne); //hexadecimal!!!
  outfile << str << endl << endl;




  outfile << "(0 \"Zones:\")\n"
	  << "(45 (1 fluid fluid)())\n"
    //      << "(45 (2 velocity-inlet velocity_inlet.1)())\n"
    //      << "(45 (3 pressure-outlet pressure_outlet.2)())\n"
	  << "(45 (2 wall wall)())\n"
	  << "(45 (4 interior default-interior)())\n" << endl;

  cout << "done" << endl;
}

}
