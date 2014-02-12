//
//  Write dolfin file
//
//  by
//  Kent-Andre Mardal <kent-and@simula.no>


#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{

#include "writeuser.hpp"



  void WriteDolfinFormat (const Mesh & mesh, const string & filename)
  {
    cout << "start writing dolfin export" << endl;

    int np = mesh.GetNP();
    int ne = mesh.GetNE();
    int nse = mesh.GetNSE();
    int nsd = mesh.GetDimension(); 
    int invertsurf = mparam.inverttrigs;
    int i, j;

    ofstream outfile (filename.c_str());

    char str[100];
    outfile.precision(8);
    outfile.setf (ios::fixed, ios::floatfield);
    outfile.setf (ios::showpoint);

    if ( nsd == 3) {

      outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" <<endl; 
      outfile << ""<<endl; 

      outfile << "<dolfin xmlns:dolfin=\"http://www.phi.chalmers.se/dolfin/\">"<<endl;
      outfile << "  <mesh celltype=\"tetrahedron\" dim=\"3\">" <<endl; 
      outfile << "      <vertices size=\""<<np<<"\">"<<endl; 
      for (i = 1; i <= np; i++) { 
        const Point3d & p = mesh.Point(i);
        outfile << "      <vertex index=\""<<i-1<<"\" x=\""<<p.X()<<"\" y=\""<<p.Y()<<"\" z=\""<<p.Z()<<"\"/>"<<endl; 
      }
      outfile << "      </vertices>"<<endl; 



      outfile << "      <cells size=\""<<ne<<"\">"<<endl; 
      for (i = 1; i <= ne; i++) {
        const Element & el = mesh.VolumeElement(i);

        outfile << "      <tetrahedron index=\""<<i-1<<"\" v0=\""<<el.PNum(1)-1<<"\" v1=\""<<el.PNum(2)-1<<"\" v2=\""<<el.PNum(3)-1<<"\" v3=\""<<el.PNum(4)-1<<"\"/>"<<endl; 
      }
      outfile << "      </cells>"<<endl; 
    }
    outfile << "   </mesh>"<<endl; 
    outfile << "</dolfin>"<<endl; 

    cout << "done writing dolfin export" << endl;
  }
}
