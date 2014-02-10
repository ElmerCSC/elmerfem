//
//  Write diffpack file
//
//  by
//  Bartosz Sawicki <sawickib@ee.pw.edu.pl>
//  extended by
//  Jacques Lechelle <jacques.lechelle@wanadoo.fr>
//
#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>


namespace netgen
{
#include "writeuser.hpp"


void WriteDiffPackFormat (const Mesh & mesh,
			  const CSGeometry & geom,
			  const string & filename)
{
  //   double scale = globflags.GetNumFlag ("scale", 1);
  double scale = 1;

  ofstream outfile(filename.c_str());

  if (mesh.GetDimension() == 3)

    {
      // Output compatible to Diffpack grid format
      // Bartosz Sawicki <sawickib@ee.pw.edu.pl>

      int np = mesh.GetNP();
      int ne = mesh.GetNE();
      int nse = mesh.GetNSE();
      ARRAY <int> BIname;
      ARRAY <int> BCsinpoint;
      int i, j, k, l;


      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      const Element & eldummy = mesh.VolumeElement((int)1);
      outfile << "\n\n"
	"Finite element mesh (GridFE):\n\n"
	"  Number of space dim. =   3\n"
	"  Number of elements   =  " << ne << "\n"
	"  Number of nodes      =  " << np << "\n\n"
	"  All elements are of the same type : dpTRUE\n"
	"  Max number of nodes in an element: "<< eldummy.GetNP() << "\n"
	"  Only one subdomain               : dpFALSE\n"
	"  Lattice data                     ? 0\n\n\n\n";
      
      for (i = 1; i <= nse; i++) 
	{
	  int BI=mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex()).BCProperty();
	  int nbi=BIname.Size();
	  int found=0;
	  for (j = 1; j <= nbi; j++)
	    if(BI == BIname.Get(j)) found = 1;
	  if( ! found ) BIname.Append(BI);	    	     
	}
      
      outfile << "  " << BIname.Size() <<  " Boundary indicators:  ";
      for (i =1 ; i <= BIname.Size(); i++)
	outfile << BIname.Get(i) << " ";
      outfile << "\n\n\n";
      
      outfile << "  Nodal coordinates and nodal boundary indicators,\n"
	"  the columns contain:\n"
	"   - node number\n"
	"   - coordinates\n"
	"   - no of boundary indicators that are set (ON)\n"
	"   - the boundary indicators that are set (ON) if any.\n"
	"#\n";

      for (i = 1; i <= np; i++)
        {
          const Point3d & p = mesh.Point(i);

          outfile.width(4);
          outfile << i << "  (";
          outfile.width(10);
          outfile << p.X()/scale << ", ";
          outfile.width(9);
          outfile << p.Y()/scale << ", ";
          outfile.width(9);
          outfile << p.Z()/scale << ") ";
	 
	  if(mesh[PointIndex(i)].Type() != INNERPOINT) 
	    {
	      BCsinpoint.DeleteAll();
	      for (j = 1; j <= nse; j++) 
		{
		  for (k = 1; k <= mesh.SurfaceElement(j).GetNP(); k++) 
		    {
		      if(mesh.SurfaceElement(j).PNum(k)==i) 
			{
			  int BC=mesh.GetFaceDescriptor(mesh.SurfaceElement(j).GetIndex()).BCProperty();
			  int nbcsp=BCsinpoint.Size();
			  int found = 0;
			  for (l = 1; l <= nbcsp; l++)
			    if(BC == BCsinpoint.Get(l)) found = 1;
			  if( ! found ) BCsinpoint.Append(BC); 	    	     
			}
		    }
		}
	      int nbcsp = BCsinpoint.Size();
	      outfile << "[" << nbcsp << "] ";
	      for (j = 1; j <= nbcsp; j++)
		outfile << BCsinpoint.Get(j) << " ";
	      outfile << "\n";
            }
          else outfile << "[0]\n";

        }

      outfile << "\n"
	"  Element types and connectivity\n"
	"  the columns contain:\n"
	"   - element number\n"
	"   - element type\n"
	"   - subdomain number\n"
	"   - the global node numbers of the nodes in the element.\n"
	"#\n";

      for (i = 1; i <= ne; i++)
        {
          const Element & el = mesh.VolumeElement(i);
          outfile.width(5);
          if(el.GetNP()==4)
            outfile << i << "  ElmT4n3D ";
          else
            outfile << i << "  ElmT10n3D ";
          outfile.width(4);
          outfile << el.GetIndex() << "    ";
          if(el.GetNP()==10)
            {
	      outfile.width(8);
	      outfile << el.PNum(1);
	      outfile.width(8);
	      outfile << el.PNum(3);
	      outfile.width(8);
	      outfile << el.PNum(2);
	      outfile.width(8);
	      outfile << el.PNum(4);
	      outfile.width(8);
	      outfile << el.PNum(6);
	      outfile.width(8);
	      outfile << el.PNum(8);
	      outfile.width(8);
	      outfile << el.PNum(5);
	      outfile.width(8);
	      outfile << el.PNum(7);
	      outfile.width(8);
	      outfile << el.PNum(10);
	      outfile.width(8);
	      outfile << el.PNum(9);
            }
          else
            {
	      outfile.width(8);
	      outfile << el.PNum(1);
	      outfile.width(8);
	      outfile << el.PNum(3);
	      outfile.width(8);
	      outfile << el.PNum(2);
	      outfile.width(8);
	      outfile << el.PNum(4);
            }
          outfile << "\n";
        }
    } /* Diffpack */

  else

    {
      // Output compatible to Diffpack grid format 2D

      int np = mesh.GetNP();
      //int ne = mesh.GetNE();
      int nse = mesh.GetNSE();
      ARRAY <int> BIname;
      ARRAY <int> BCsinpoint;
      int i, j, k, l;


      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      outfile << "\n\n"
	"Finite element mesh (GridFE):\n\n"
	"  Number of space dim. =  2\n"
	"  Number of elements   =  " << nse << "\n"
	"  Number of nodes      =  " << np << "\n\n"
	"  All elements are of the same type : dpTRUE\n"
	"  Max number of nodes in an element: 3\n"
	"  Only one subdomain               : dpFALSE\n"
	"  Lattice data                     ? 0\n\n\n\n";
      
      for (i = 1; i <= nse; i++) 
	{
	  int BI=mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex()).BCProperty();
	  int nbi=BIname.Size();
	  int found=0;
	  for (j = 1; j <= nbi; j++)
	    if(BI == BIname.Get(j)) found = 1;
	  if( ! found ) BIname.Append(BI);	    	     
	}
      
      outfile << "  " << BIname.Size() <<  " Boundary indicators:  ";
      for (i =1 ; i <= BIname.Size(); i++)
	outfile << BIname.Get(i) << " ";
      outfile << "\n\n\n";
      
      outfile << "  Nodal coordinates and nodal boundary indicators,\n"
	"  the columns contain:\n"
	"   - node number\n"
	"   - coordinates\n"
	"   - no of boundary indicators that are set (ON)\n"
	"   - the boundary indicators that are set (ON) if any.\n"
	"#\n";

      for (i = 1; i <= np; i++)
        {
          const Point3d & p = mesh.Point(i);

          outfile.width(4);
          outfile << i << "  (";
          outfile.width(10);
          outfile << p.X()/scale << ", ";
          outfile.width(9);
          outfile << p.Y()/scale << ", ";
	 
	  if(mesh[PointIndex(i)].Type() != INNERPOINT) 
	    {
	      BCsinpoint.DeleteAll();
	      for (j = 1; j <= nse; j++) 
		{
		  for (k = 1; k <= 2; k++) 
		    {
		      if(mesh.SurfaceElement(j).PNum(k)==i) 
			{
			  int BC=mesh.GetFaceDescriptor(mesh.SurfaceElement(j).GetIndex()).BCProperty();
			  int nbcsp=BCsinpoint.Size();
			  int found = 0;
			  for (l = 1; l <= nbcsp; l++)
			    if(BC == BCsinpoint.Get(l)) found = 1;
			  if( ! found ) BCsinpoint.Append(BC); 	    	     
			}
		    }
		}
	      int nbcsp = BCsinpoint.Size();
	      outfile << "[" << nbcsp << "] ";
	      for (j = 1; j <= nbcsp; j++)
		outfile << BCsinpoint.Get(j) << " ";
	      outfile << "\n";
            }
          else outfile << "[0]\n";

        }

      outfile << "\n"
	"  Element types and connectivity\n"
	"  the columns contain:\n"
	"   - element number\n"
	"   - element type\n"
	"   - subdomain number\n"
	"   - the global node numbers of the nodes in the element.\n"
	"#\n";

      for (i = 1; i <= nse; i++)
        {
          const Element2d & el = mesh.SurfaceElement(i);
          outfile.width(5);
          outfile << i << "  ElmT3n2D ";
          outfile.width(4);
          outfile << el.GetIndex() << "    ";
	  outfile.width(8);
	  outfile << el.PNum(1);
	  outfile.width(8);
	  outfile << el.PNum(3);
	  outfile.width(8);
	  outfile << el.PNum(2);
          outfile << "\n";
        }
    }
}
}
