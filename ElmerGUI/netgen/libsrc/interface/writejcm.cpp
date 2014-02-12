//
//  Write JCMwave file
//  07.07.2005, Sven Burger, ZIB Berlin
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

void WriteJCMFormat (const Mesh & mesh,
                     const CSGeometry & geom,
                     const string & filename)
{
  if (mesh.GetDimension() != 3)
  {
    cout <<"\n Error: Dimension 3 only supported by this output format!"<<endl;
    return;
  }

  int bc_at_infinity = 0;
  int i, j, jj, ct(0), counter;
  double dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3, vol;

  // number of points
  int np = mesh.GetNP();

  // Identic points
  ARRAY<int,1> identmap1, identmap2, identmap3;
  mesh.GetIdentifications().GetMap(1, identmap1);
  mesh.GetIdentifications().GetMap(2, identmap2);
  mesh.GetIdentifications().GetMap(3, identmap3);

  // number of volume elements
  int ne = mesh.GetNE();
  int ntets = 0;
  int nprisms = 0;
  for (i = 1; i <= ne; i++)
  {
    Element el = mesh.VolumeElement(i);
    if (el.GetNP() == 4)
    {
      ntets++;
      // Check that no two points on a tetrahedron are identified with each other
      for (j = 1; j <= 4; j++)
        for (jj = 1; jj <=4; jj++)
        {
          if (identmap1.Elem(el.PNum(j)) == el.PNum(jj))
          {
            cout << "\n Error: two points on a tetrahedron identified (1) with each other"
                 << "\n REFINE MESH !" << endl;
            return;
          }
          if (identmap2.Elem(el.PNum(j)) == el.PNum(jj))
          {
            cout << "\n Error: two points on a tetrahedron identified (2) with each other"
                 << "\n REFINE MESH !" << endl;
            return;
          }
          if (identmap3.Elem(el.PNum(j)) == el.PNum(jj))
          {
            cout << "\n Error: two points on a tetrahedron identified (3) with each other"
                 << "\n REFINE MESH !" << endl;
            return;
          }
        }      
      
    }
    else if (el.GetNP() == 6)
      nprisms++;
  }
  if ( ne != (ntets+nprisms))
  {
    cout<< "\n Error in determining number of volume elements!\n"
        << "\n Prisms and tetrahedra only implemented in the JCMwave format!\n"<<endl;
    return;
  }

  if (nprisms > 0)
    cout << " Please note: Boundaries at infinity have to carry the bc-attribute '-bc="
         << bc_at_infinity <<"'."<<endl; 

  // number of surface elements
  int nse = mesh.GetNSE();
  // number of boundary triangles
  int nbtri = 0;
  // number of boundary quadrilaterals
  int nbquad = 0;
  // array with 1 if point on any tetra, 0 else 
  // this is needed in order to arrange the prism points in the right order
  ARRAY<int,1> pointsOnTetras;
  pointsOnTetras.SetSize (mesh.GetNP());
  pointsOnTetras = 0;
  for (i = 1; i <= ne; i++)
  {
    Element el = mesh.VolumeElement(i);
    if (el.GetNP() == 4)
    {
      for (j = 1; j <= 4; j++)
        pointsOnTetras.Set(el.PNum(j).GetInt(),1);     
    }
  }

  // number of boundary triangles and boundary quadrilaterals
  for (i = 1; i <= nse; i++)
  {
    Element2d el = mesh.SurfaceElement(i);
    if (el.GetNP() == 3 &&
        ( mesh.GetFaceDescriptor (el.GetIndex()).DomainIn()==0  ||
          mesh.GetFaceDescriptor (el.GetIndex()).DomainOut()==0 ) )
      nbtri++;
    else if (el.GetNP() == 4 &&
             ( mesh.GetFaceDescriptor (el.GetIndex()).DomainIn()==0 ||
               mesh.GetFaceDescriptor (el.GetIndex()).DomainOut()==0 ) )
      nbquad++;
  }
  
  ofstream outfile (filename.c_str());
  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);
  
  outfile << "/* <BLOBHead>\n";
  outfile << "__BLOBTYPE__=Grid\n";
  outfile << "__OWNER__=JCMwave\n";
  outfile << "<I>SpaceDim=3\n";
  outfile << "<I>ManifoldDim=3\n";
  outfile << "<I>NRefinementSteps=0\n";
  outfile << "<I>NPoints="<<np<<"\n";
  outfile << "<I>NTetrahedra="<<ntets<<"\n";
  outfile << "<I>NPrisms="<<nprisms<<"\n";
  outfile << "<I>NBoundaryTriangles="<<nbtri<<"\n";
  outfile << "<I>NBoundaryQuadrilaterals="<<nbquad<<"\n";
  outfile << "*/\n";
  outfile << "\n";
  outfile << "# output from Netgen\n\n";
  int nDomains=mesh.GetNDomains();
  for (i=1; i<=nDomains; i++)
  {
    if (mesh.GetMaterial(i))
      outfile << "#" << mesh.GetMaterial(i) 
              << ": Material ID = " 
              << i << "\n";
  }

  outfile << "# Points\n";
  cout << " Please note: The unit of length in the .geo file is assumed to be 'microns'."<<endl; 
  for (i = 1; i <= np; i++)
  {
    const Point<3> & p = mesh.Point(i);
    outfile << i << "\n";
    outfile << p(0) << "e-6\n";
    outfile << p(1) << "e-6\n";
    outfile << p(2) << "e-6\n\n";
  }

  outfile << "\n";
  outfile << "# Tetrahedra\n";
  counter = 0;
  for (i = 1; i <= ne; i++)
  {
    Element el = mesh.VolumeElement(i);
    if (el.GetNP() == 4)
    {
      counter++;
      dx1 = mesh.Point(el.PNum(2))(0) - mesh.Point(el.PNum(1))(0);
      dx2 = mesh.Point(el.PNum(3))(0) - mesh.Point(el.PNum(1))(0);
      dx3 = mesh.Point(el.PNum(4))(0) - mesh.Point(el.PNum(1))(0);
      dy1 = mesh.Point(el.PNum(2))(1) - mesh.Point(el.PNum(1))(1);
      dy2 = mesh.Point(el.PNum(3))(1) - mesh.Point(el.PNum(1))(1);
      dy3 = mesh.Point(el.PNum(4))(1) - mesh.Point(el.PNum(1))(1);
      dz1 = mesh.Point(el.PNum(2))(2) - mesh.Point(el.PNum(1))(2);
      dz2 = mesh.Point(el.PNum(3))(2) - mesh.Point(el.PNum(1))(2);
      dz3 = mesh.Point(el.PNum(4))(2) - mesh.Point(el.PNum(1))(2);
      vol = (dy1*dz2-dz1*dy2)*dx3 + (dz1*dx2-dx1*dz2)*dy3 + (dx1*dy2-dy1*dx2)*dz3;

      if ( vol > 0 )
        for (j = 1; j <= 4; j++)
          outfile << el.PNum(j)<<"\n";
      else
      {
        for (j = 2; j >= 1; j--)
          outfile << el.PNum(j)<<"\n";
        for (j = 3; j <= 4; j++)
          outfile << el.PNum(j)<<"\n";
      }  
      outfile << el.GetIndex() << "\n\n";
    }
  }
  if ( counter != ntets)
  {
    cout<< "\n Error in determining number of tetras!\n"<<endl;
    return;
  }

  outfile << "\n";
  outfile << "# Prisms\n";
  counter = 0;
  for (i = 1; i <= ne; i++)
  {
    Element el = mesh.VolumeElement(i);
    if (el.GetNP() == 6)
    {
      counter++;
      dx1 = mesh.Point(el.PNum(2))(0) - mesh.Point(el.PNum(1))(0);
      dx2 = mesh.Point(el.PNum(3))(0) - mesh.Point(el.PNum(1))(0);
      dx3 = mesh.Point(el.PNum(4))(0) - mesh.Point(el.PNum(1))(0);
      dy1 = mesh.Point(el.PNum(2))(1) - mesh.Point(el.PNum(1))(1);
      dy2 = mesh.Point(el.PNum(3))(1) - mesh.Point(el.PNum(1))(1);
      dy3 = mesh.Point(el.PNum(4))(1) - mesh.Point(el.PNum(1))(1);
      dz1 = mesh.Point(el.PNum(2))(2) - mesh.Point(el.PNum(1))(2);
      dz2 = mesh.Point(el.PNum(3))(2) - mesh.Point(el.PNum(1))(2);
      dz3 = mesh.Point(el.PNum(4))(2) - mesh.Point(el.PNum(1))(2);
      vol = (dy1*dz2-dz1*dy2)*dx3 + (dz1*dx2-dx1*dz2)*dy3 + (dx1*dy2-dy1*dx2)*dz3;

      if (pointsOnTetras.Get(el.PNum(1)) &&
          pointsOnTetras.Get(el.PNum(2)) &&
          pointsOnTetras.Get(el.PNum(3)))
      {
        if (vol > 0)
          for (j = 1; j <= 6; j++)
            outfile << el.PNum(j)<<"\n";
        else
        {
          for (j = 3; j >= 1; j--)
            outfile << el.PNum(j)<<"\n";
          for (j = 6; j >= 4; j--)
            outfile << el.PNum(j)<<"\n";
        }
      }
      else if ( pointsOnTetras.Get(el.PNum(4)) &&
                pointsOnTetras.Get(el.PNum(5)) &&
                pointsOnTetras.Get(el.PNum(6))    )
      {
        if ( vol < 0 )
        {
          for (j = 4; j <= 6; j++)
            outfile << el.PNum(j)<<"\n";
          for (j = 1; j <= 3; j++)
            outfile << el.PNum(j)<<"\n";
        }
        else
        {
          for (j = 6; j >= 4; j--)
            outfile << el.PNum(j)<<"\n";
          for (j = 3; j >= 1; j--)
            outfile << el.PNum(j)<<"\n";
        }
      }
      else 
      {
        cout << "\n Error in determining prism point numbering!\n"<<endl;
        return;
      }
      outfile << el.GetIndex() << "\n\n";
    }
  }
  if ( counter != nprisms)
  {
    cout<< "\n Error in determining number of prisms!\n"<<endl;
    return;
  }

  int npid1 = 0;
  int npid2 = 0;
  int npid3 = 0;
  for (i=1; i<=np; i++)
  {
    if (identmap1.Elem(i))
      npid1++;
    if (identmap2.Elem(i))
      npid2++;
    if (identmap3.Elem(i))
      npid3++;
  }

  outfile << "\n";
  outfile << "# Boundary triangles\n";  
  outfile << "# Number of identified points in 1-direction: " << npid1 << "\n";
  outfile << "# Number of identified points in 2-direction: " << npid2 << "\n";
  outfile << "# Number of identified points in 3-direction: " << npid3 << "\n";
  for (i = 1; i <= nse; i++)
  {
    Element2d el = mesh.SurfaceElement(i);
    if (el.GetNP() == 3
        && (mesh.GetFaceDescriptor (el.GetIndex()).DomainIn()==0
            || mesh.GetFaceDescriptor (el.GetIndex()).DomainOut()==0))
    {
      outfile <<"# T\n";
      for (j = 1; j <= 3; j++)
        outfile << el.PNum(j)<<"\n";
      if (mesh.GetFaceDescriptor (el.GetIndex()).BCProperty()==bc_at_infinity)
        outfile << 1000 << "\n";      
      else
        outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << "\n";      
      if (mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() == bc_at_infinity)
        outfile << "-2\n\n";
      else if (identmap1.Elem(el.PNum(1))
               &&identmap1.Elem(el.PNum(2))
               &&identmap1.Elem(el.PNum(3)))
      {
        outfile << "-1\n";
        for (j = 1; j <= 3; j++)
          outfile << identmap1.Elem(el.PNum(j))<<"\n";
        outfile << "\n";
      }
      else if (identmap2.Elem(el.PNum(1))
               &&identmap2.Elem(el.PNum(2))
               &&identmap2.Elem(el.PNum(3)))
      {
        outfile << "-1\n";
        for (j = 1; j <= 3; j++)
          outfile << identmap2.Elem(el.PNum(j))<<"\n";
        outfile << "\n";
      }
      else if (identmap3.Elem(el.PNum(1))
               &&identmap3.Elem(el.PNum(2))
               &&identmap3.Elem(el.PNum(3)))
      {
        outfile << "-1\n";
        for (j = 1; j <= 3; j++)
          outfile << identmap3.Elem(el.PNum(j))<<"\n";
        outfile << "\n";
      }
      else
        outfile << "1\n\n";
        
    }
  }

  outfile << "\n";
  outfile << "# Boundary quadrilaterals\n";
  for (i = 1; i <= nse; i++)
  {
    Element2d el = mesh.SurfaceElement(i);

    if (el.GetNP() == 4
        && (mesh.GetFaceDescriptor (el.GetIndex()).DomainIn()==0
            || mesh.GetFaceDescriptor (el.GetIndex()).DomainOut()==0))
    {
      if      (pointsOnTetras.Get(el.PNum(1)) &&
               pointsOnTetras.Get(el.PNum(2)))
        ct = 0;
      else if (pointsOnTetras.Get(el.PNum(2)) &&
               pointsOnTetras.Get(el.PNum(3)))
        ct = 1;
      else if (pointsOnTetras.Get(el.PNum(3)) &&
               pointsOnTetras.Get(el.PNum(4)))
        ct = 2;
      else if (pointsOnTetras.Get(el.PNum(4)) &&
               pointsOnTetras.Get(el.PNum(1)))
        ct = 3;
      else
        cout << "\nWarning: Quadrilateral with inconsistent points found!"<<endl;
      
      for (j = 1; j <= 4; j++)
      {
        jj = j + ct;
        if ( jj >= 5 )
          jj = jj - 4;
        outfile << el.PNum(jj)<<"\n";
      }
      outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << "\n";      
      if (mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() == bc_at_infinity)
      {
        outfile << "-2\n\n";
        cout << "\nWarning: Quadrilateral at infinity found (this should not occur)!"<<endl;
      }
      else if ( identmap1.Elem(el.PNum(1)) &&
                identmap1.Elem(el.PNum(2)) &&
                identmap1.Elem(el.PNum(3)) &&
                identmap1.Elem(el.PNum(4))    )
      {
        outfile << "-1\n";
        for (j = 1; j <= 4; j++)
        {
          jj = j + ct;
          if ( jj >= 5 )
            jj = jj - 4;
          outfile << identmap1.Elem(el.PNum(jj))<<"\n";
        }
        outfile << "\n";
      }
      else if ( identmap2.Elem(el.PNum(1)) &&
                identmap2.Elem(el.PNum(2)) &&
                identmap2.Elem(el.PNum(3)) &&
                identmap2.Elem(el.PNum(4))    )
      {
        outfile << "-1\n";
        for (j = 1; j <= 4; j++)
        {
          jj = j + ct;
          if ( jj >= 5 )
            jj = jj - 4;
          outfile << identmap2.Elem(el.PNum(jj))<<"\n";
        }
        outfile << "\n";
      }
      else if ( identmap3.Elem(el.PNum(1)) &&
                identmap3.Elem(el.PNum(2)) &&
                identmap3.Elem(el.PNum(3)) &&
                identmap3.Elem(el.PNum(4))    )
      {
        outfile << "-1\n";
        for (j = 1; j <= 4; j++)
        {
          jj = j + ct;
          if ( jj >= 5 )
            jj = jj - 4;
          outfile << identmap3.Elem(el.PNum(jj))<<"\n";
        }
        outfile << "\n";
      }
      else
        outfile << "1\n\n";
    }
  }

  cout << " JCMwave grid file written." << endl;
}

}

