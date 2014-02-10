//
//  Write Tochnog file
//
//  by
//
//  Andreas Seltmann
//  email:  A.Seltmann@lsw.uni-heidelberg.de
//
#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>


namespace netgen
{
#include "writeuser.hpp"


void WriteTochnogFormat (const Mesh & mesh,
			 const string & filename)
{
  cout << "\nWrite Tochnog Volume Mesh" << endl;

  ofstream outfile (filename.c_str());

  outfile << "(Nodes and Elements generated with NETGEN" << endl;
  outfile << " " << filename << ")" << endl;

  outfile.precision(8);

  outfile << "(Nodes)" << endl;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int i, j;

  for (i = 1; i <= np; i++)
    {
      outfile << "node " << " " << i << " ";
      outfile << mesh.Point(i)(0) << " ";
      outfile << mesh.Point(i)(1) << " ";
      outfile << mesh.Point(i)(2) << "\n";
    }

  int elemcnt = 0; //element counter
  int finished = 0;
  int indcnt = 1; //index counter

  while (!finished)
    {
      int actcnt = 0;
      const Element & el1 = mesh.VolumeElement(1);
      int non = el1.GetNP();
      if (non == 4)
	{
	  outfile << "(Elements, type=-tet4)" << endl;
	} 
      else
	{
	  cout << "unsupported Element type!!!" << endl;	  
	}

      for (i = 1; i <= ne; i++)
	{
	  const Element & el = mesh.VolumeElement(i);
	      
	  if (el.GetIndex() == indcnt)
	    {
	      actcnt++;
	      if (el.GetNP() != non) 
		{
		  cout << "different element-types in a subdomain are not possible!!!" << endl;
		  continue;
		}
		  
	      elemcnt++;
	      outfile << "element " << elemcnt << " -tet4 ";
	      if (non == 4)
		{
		  outfile << el.PNum(1) << " ";
		  outfile << el.PNum(2) << " ";
		  outfile << el.PNum(4) << " ";
		  outfile << el.PNum(3) << "\n";
		}
	      else
		{
		  cout << "unsupported Element type!!!" << endl;
		  for (j = 1; j <= el.GetNP(); j++)
		    {
		      outfile << el.PNum(j);
		      if (j != el.GetNP()) outfile << ", ";
		    }
		  outfile << "\n";
		}
	    }
	}	  
      indcnt++;
      if (elemcnt == ne) {finished = 1; cout << "all elements found by Index!" << endl;}
      if (actcnt == 0) {finished = 1;}
    }

  cout << "done" << endl;
}

}
