/*************************************
 * Write Gmsh file
 * First issue the 04/26/2004 by Paul CARRICO (paul.carrico@free.fr)
 * At the moment, the GMSH format is available for
 * linear tetrahedron elements i.e. in 3D
 * (based on Neutral Format)
 *
 * Second issue the 05/05/2004 by Paul CARRICO
 * Thanks to Joachim Schoeberl for the correction of a minor bug
 * the 2 initial Gmsh Format (i.e. volume format and surface format) are group together)
 * in only one file
 **************************************/

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "writeuser.hpp"



/*
 *  GMSH mesh format
 *  points, elements, surface elements and physical entities
 */

void WriteGmshFormat (const Mesh & mesh,
			 const CSGeometry & geom,
			 const string & filename)
{
  ofstream outfile (filename.c_str());
  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);

  int np = mesh.GetNP();  /// number of point
  int ne = mesh.GetNE();  /// number of element
  int nse = mesh.GetNSE();  /// number of surface element (BC)
  int i, j, k, l;


  /*
   * 3D section : Linear volume elements (only tetrahedra)
   */

   if (ne > 0 && mesh.VolumeElement(1).GetNP() == 4)
      {
      cout << "Write GMSH Format \n";
      cout << "The GMSH format is available for linear tetrahedron elements only in 3D\n" << endl;

      int inverttets = mparam.inverttets;
      int invertsurf = mparam.inverttrigs;


      /// Write nodes
      outfile << "$NOD\n";
      outfile << np << "\n";
  
      for (i = 1; i <= np; i++)
          {
          const Point3d & p = mesh.Point(i);
          outfile << i << " "; /// node number
          outfile << p.X() << " ";
          outfile << p.Y() << " ";
          outfile << p.Z() << "\n";
          }
      outfile << "$ENDNOD\n";

      /// write elements
      outfile << "$ELM\n";
      outfile << ne + nse << "\n";  ////  number of elements + number of surfaces BC

     for (i = 1; i <= nse; i++)
         {
         Element2d el = mesh.SurfaceElement(i);
         if (invertsurf) el.Invert();
         outfile << i;
         outfile << " ";
         outfile << "2";
         outfile << " ";
         outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
         /// that means that physical entity = elementary entity (arbitrary approach)
         outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
         outfile << "3";
         outfile << " ";
                 for (j = 1; j <= el.GetNP(); j++)
                     {
                     outfile << " ";
                     outfile << el.PNum(j);
                     }
                     outfile << "\n";
         }


         for (i = 1; i <= ne; i++)
             {
             Element el = mesh.VolumeElement(i);
             if (inverttets) el.Invert();
             outfile << nse + i; /// element number
             outfile << " ";
             outfile << "4"; /// element type i.e. Tetraedron == 4
             outfile << " ";
             outfile << 100000 + el.GetIndex();
             /// that means that physical entity = elementary entity (arbitrary approach)
             outfile << " ";
             outfile << 100000 + el.GetIndex();   /// volume number
             outfile << " ";
             outfile << "4";  /// number of nodes i.e. 4 for a tetrahedron
                                                                                                        
             for (j = 1; j <= el.GetNP(); j++)
                 {
                 outfile << " ";
                 outfile << el.PNum(j);
                 }
             outfile << "\n";
             }


             outfile << "$ENDELM\n";
   }

   /*
    * End of 3D section
    */





  /*
   * 2D section : available for triangles and quadrangles
   */
      else if (ne == 0)   /// means that there's no 3D element
              {
              cout << "\n Write Gmsh Surface Mesh (triangle and/or quadrangles)" << endl;

              /// Write nodes
              outfile << "$NOD\n";
              outfile << np << "\n";

              for (i = 1; i <= np; i++)
              {
              const Point3d & p = mesh.Point(i);
              outfile << i << " "; /// node number
              outfile << p.X() << " ";
              outfile << p.Y() << " ";
              outfile << p.Z() << "\n";
              }
              outfile << "$ENDNOD\n";


              /// write triangles & quadrangles
              outfile << "$ELM\n";
              outfile << nse << "\n";

              for (k = 1; k <= nse; k++)
              {
              const Element2d & el = mesh.SurfaceElement(k);


              outfile << k;
              outfile << " ";
              outfile << (el.GetNP()-1);   // 2 for a triangle and 3 for a quadrangle
              outfile << " ";
              outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
              /// that means that physical entity = elementary entity (arbitrary approach)
              outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
              outfile << (el.GetNP());    // number of node per surfacic element
              outfile << " ";

              for (l = 1; l <= el.GetNP(); l++)
                  {
                  outfile << " ";
                  outfile << el.PNum(l);
                  }
	              outfile << "\n";
		  
               }
               outfile << "$ENDELM$ \n";
    }

   /*
    * End of 2D section
    */

     else
    {
    cout << " Invalide element type for Gmsh volume Format !\n";
    }


}
}


