
//
//  Write Elmer file
//
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>
#include <sys/stat.h>


namespace netgen
{
#include "writeuser.hpp"



void WriteElmerFormat (const Mesh &mesh,
			 const string &filename)
{
  cout << "write elmer mesh files" << endl;
  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int i, j;
  char str[200];
  
  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;

#ifdef WIN32
  char a[256];
  sprintf( a, "mkdir %s", filename.c_str() );
  system( a );
#else
  int rc = mkdir(filename.c_str(), S_IRWXU|S_IRWXG);
#endif

  sprintf( str, "%s/mesh.header", filename.c_str() );
  ofstream outfile_h(str);
  sprintf( str, "%s/mesh.nodes", filename.c_str() );
  ofstream outfile_n(str);
  sprintf( str, "%s/mesh.elements", filename.c_str() );
  ofstream outfile_e(str);
  sprintf( str, "%s/mesh.boundary", filename.c_str() );
  ofstream outfile_b(str);

  // fill hashtable

  INDEX_3_HASHTABLE<int> face2volelement(ne);

  for (i = 1; i <= ne; i++)
    {
      const Element & el = mesh.VolumeElement(i);
      INDEX_3 i3;
      int k, l;
      for (j = 1; j <= 4; j++)   // loop over faces of tet
	{
	  l = 0;
	  for (k = 1; k <= 4; k++)
	    if (k != j)
	      {
		l++;
		i3.I(l) = el.PNum(k);
	      }
	  i3.Sort();
	  face2volelement.Set (i3, i);
	}
    }

//  outfile.precision(6);
//  outfile.setf (ios::fixed, ios::floatfield);
//  outfile.setf (ios::showpoint);
  
  outfile_h << np << " " << ne << " " << nse << "\n";
  outfile_h << "2"     << "\n";
  outfile_h << "303 " << nse << "\n";
  outfile_h << "504 " << ne << "\n";
  
  for (i = 1; i <= np; i++)
    {
      const Point3d & p = mesh.Point(i);
      
      outfile_n << i << " -1 ";
      outfile_n << p.X() << " ";
      outfile_n << p.Y() << " ";
      outfile_n << p.Z() << "\n";
    }

  for (i = 1; i <= ne; i++)
    {
      Element el = mesh.VolumeElement(i);
      if (inverttets) el.Invert();
      sprintf( str, "5%02d", (int)el.GetNP() );
      outfile_e << i << " " << el.GetIndex() << " " << str <<  "  ";
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile_e << " ";
	  outfile_e << el.PNum(j);
	}
      outfile_e << "\n";
    }

  for (i = 1; i <= nse; i++)
    {
      Element2d el = mesh.SurfaceElement(i);
      if (invertsurf) el.Invert();
      sprintf( str, "3%02d", (int)el.GetNP() );
      {
	  INDEX_3 i3;
	  for (j = 1; j <= 3; j++) i3.I(j) = el.PNum(j);
	  i3.Sort();
	  
	  int elind = face2volelement.Get(i3);
          outfile_b << i << " " << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << 
         " " << elind << " 0 "  << str << "    ";
      }
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile_b << " ";
	  outfile_b << el.PNum(j);
	}
      outfile_b << "\n";
    }
}

}
