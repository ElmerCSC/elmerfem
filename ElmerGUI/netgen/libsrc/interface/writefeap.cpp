//
// Write FEAP file
// FEAP by Bob Taylor, Berkely
//
// contact Peter Wriggers or Albrecht Rieger, Hannover
// rieger@ibnm.uni-hannover.de
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{

#include "writeuser.hpp"


void WriteFEAPFormat (const Mesh & mesh,
		      const string & filename)
  
{
  // Feap format by A. Rieger 
  // rieger@ibnm.uni-hannover.de

  int inverttets = mparam.inverttets;
  //int invertsurf = mparam.inverttrigs;

  int i, j;

  double scale = 1;   // globflags.GetNumFlag ("scale", 1);
  
  ofstream outfile(filename.c_str());
 
  outfile << "feap" << "\n";
  outfile << mesh.GetNP();
  outfile << ",";
  outfile << mesh.GetNE();
  outfile << ",";
  outfile << "1,3,3,4" << "\n" << "\n"; 
  outfile << "!numnp,numel,nummat,ndm,ndf,nen";
  outfile << "\n";
      
  outfile << "\n" << "\n";
  outfile << "!node,,         X           Y           Z" << "\n";
  outfile << "COOR" << "\n";
  outfile.precision(4);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);

  for (i = 1; i <= mesh.GetNP(); i++)
    {
      outfile.width(5);
      outfile << i;
      outfile << ",,";
      outfile.width(10);
      outfile << mesh.Point(i)(0)/scale << "  ";
      outfile.width(10);
      outfile << mesh.Point(i)(1)/scale << "  ";
      outfile.width(10);
      outfile << mesh.Point(i)(2)/scale << "\n";
    }   
      
  outfile << "\n" << "\n";
  outfile << "!elm,,mat,     n1      n2      n3      n4" << "\n";
  outfile << "ELEM" << "\n";

  for (i = 1; i <= mesh.GetNE(); i++)
    {
      Element el = mesh.VolumeElement(i);
      if (inverttets)
	el.Invert();


      outfile.width(5);
      outfile << i;
      outfile << ",,";
      outfile << el.GetIndex();
      outfile << ",";


      for (j = 1; j <= el.NP(); j++)
	{
	  outfile.width(8);
	  outfile << el.PNum(j);
	}
      outfile << "\n";
    }
      
  outfile << "\n" << "\n";
      
      
  /*

  //outfile << "SLOA" << "\n";
  //outfile << "2,3,3" << "\n";
  //outfile << GetNSE() << "\n";
  outfile << "selm" << "\n" << GetNSE() << "\n";
  for (i = 1; i <= GetNSE(); i++)
  {
  if (SurfaceElement(i).GetIndex())
  {
  outfile.width(8);
  outfile << facedecoding.Get(SurfaceElement(i).GetIndex ()).surfnr;
  //outfile.width(8);	  
  //outfile << facedecoding.Get(SurfaceElement(i).GetIndex ()).domin;
  //outfile.width(8);	  
  //outfile << facedecoding.Get(SurfaceElement(i).GetIndex ()).domout;
  }
  else
  outfile << "       0       0       0";
  
  
  Element2d sel = SurfaceElement(i);
  if (invertsurf)
  sel.Invert();
  //outfile.width(8);
  //outfile << sel.GetNP();
  //if (facedecoding.Get(SurfaceElement(i).GetIndex ()).surfnr == 4)
  //{
  for (j = 1; j <= sel.GetNP(); j++)
  {
  outfile.width(8);	  
  outfile << sel.PNum(j);
  }
  //outfile.width(8);	
  //outfile << "0.0";
  //outfile.width(8);	
  //outfile << "0.0";
  //outfile.width(8);	
  //outfile << "1.0" << "\n";
  //}
  outfile << "\n";
  //outfile << endl;
  }
  */



  // BEGIN CONTACT OUTPUT
  /*      
	  int masterindex, slaveindex;
	  cout << "Master Surface index = ";
	  cin >> masterindex;
	  cout << "Slave Surface index  = ";
	  cin >> slaveindex;


	  // CONTACT SURFACE 1
	  outfile << "\n";
	  outfile << "\n";
	  outfile << "surface,1" << "\n";;
	  outfile.width(6);
	  outfile << "tria" << "\n";;
	  outfile.width(13);
	  outfile << "facet" << "\n";;
	  zz = 0;
	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	  Element2d sel = mesh.SurfaceElement(i);
	  if (invertsurf)
	  sel.Invert();
	  if (mesh.GetFaceDescriptor(sel.GetIndex ()).BCProperty() == masterindex)
	  {
	  zz++;
	  outfile.width(14);
	  outfile << zz;
	  outfile << ",,";
	  for (j = 1; j <= sel.GetNP(); j++)
	  {
	  outfile << sel.PNum(j);
	  outfile << ",";
	  }
	  outfile << "\n";
	  }
	  }


	  // CONTACT SURFACE 2
	  outfile << "\n";
	  outfile << "\n";
	  outfile << "surface,2" << "\n";;
	  outfile.width(6);
	  outfile << "tria" << "\n";;
	  outfile.width(13);
	  outfile << "facet" << "\n";;
	  zz = 0;
	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	  
	  Element2d sel = mesh.SurfaceElement(i);
	  if (invertsurf)
	  sel.Invert();
	  if (mesh.GetFaceDescriptor(sel.GetIndex ()).BCProperty() == slaveindex)
	  {
	  zz++;
	  outfile.width(14);
	  outfile << zz;
	  outfile << ",,";
	  for (j = 1; j <= sel.GetNP(); j++)
	  {
	  outfile << sel.PNum(j);
	  outfile << ",";
	  }
	  outfile << "\n";
	  }
	  }
      
	  outfile << "\n";
	  outfile << "\n";
  */      
      
  // END CONTACT OUTPUT

  cout << "done" << endl;
}
}
