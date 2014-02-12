//
//  Read user dependent output file
//


#include <mystdlib.h>


#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "writeuser.hpp"

  void ReadFile (Mesh & mesh,
                 const string & hfilename)
  {
    cout << "Read User File" << endl;

    const char * filename = hfilename.c_str();

    int i, j;

    char reco[100];
    int np, nbe;



    // ".surf" - mesh
  
    if ( (strlen (filename) > 5) &&
         strcmp (&filename[strlen (filename)-5], ".surf") == 0 )
    
      {
        cout << "Surface file" << endl;
      
        ifstream in (filename);
      
        in >> reco;
        in >> np;
        for (i = 1; i <= np; i++)
          {
            Point3d p;
            in >> p.X() >> p.Y() >> p.Z();
            mesh.AddPoint (p);
          }

        mesh.ClearFaceDescriptors();
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
      
        in >> nbe;
        //      int invert = globflags.GetDefineFlag ("invertsurfacemesh");
        for (i = 1; i <= nbe; i++)
          {
            Element2d el;
            int hi;
	  
            el.SetIndex(1);
	  
            // in >> hi; 
            for (j = 1; j <= 3; j++)
              {
                in >> el.PNum(j);
                // el.PNum(j)++;
                if (el.PNum(j) < PointIndex(1) || 
                    el.PNum(j) > PointIndex(np))
                  {
                    cerr << "Point Number " << el.PNum(j) << " out of range 1..."
                         << np << endl;
                    return;
                  }
              }
	  
            /*
              if (invert)
              swap (el.PNum(2), el.PNum(3));
            */
	  
            mesh.AddSurfaceElement (el);
          }
      
      
        cout << "points: " << np << " faces: " << nbe << endl;
      }
  
  
  

  
    // Universal mesh (AVL)
    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".unv") == 0 )
      {  
        int i, j, k;
      
        double h;
        char reco[100];
        int np, nbe;
        int invert;
      
      
        ifstream in(filename);

        invert = 0;    // globflags.GetDefineFlag ("invertsurfacemesh");
        double scale = 1;  // globflags.GetNumFlag ("scale", 1);

        mesh.ClearFaceDescriptors();
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));


        while (in.good())
          {
            in >> reco;
            if (strcmp (reco, "NODES") == 0)
              {
                cout << "nodes found" << endl;
                for (j = 1; j <= 4; j++)
                  in >> reco;  // read dummy

                while (1)
                  {
                    int pi, hi;
                    double x, y, z;
                    Point3d p;

                    in >> pi;
                    if (pi == -1)
                      break;

                    in >> hi >> hi >> hi;
                    in >> p.X() >> p.Y() >> p.Z();

                    p.X() *= scale;
                    p.Y() *= scale;
                    p.Z() *= scale;


                    mesh.AddPoint (p);
                  }
              }

            if (strcmp (reco, "ELEMENTS") == 0)
              {
                cout << "elements found" << endl;
                for (j = 1; j <= 4; j++)
                  in >> reco;  // read dummy

                while (1)
                  {
                    int hi;
                    in >> hi;
                    if (hi == -1) break;
                    for (j = 1; j <= 7; j++)
                      in >> hi;
	      
                    Element2d el;
                    el.SetIndex(1);
                    in >> el.PNum(1) >> el.PNum(2) >> el.PNum(3);
	      
                    if (invert)
                      swap (el.PNum(2), el.PNum(3));
                    mesh.AddSurfaceElement (el);	  
	      
                    for (j = 1; j <= 5; j++)
                      in >> hi;
                  }
              }
          }
      

        Point3d pmin, pmax;
        mesh.GetBox (pmin, pmax);
        cout << "bounding-box = " << pmin << "-" << pmax << endl;
      }



    // fepp format2d:
  
    if ( (strlen (filename) > 7) &&
         strcmp (&filename[strlen (filename)-7], ".mesh2d") == 0 )
      {
        cout << "Reading FEPP2D Mesh" << endl;
      
        char buf[100];
        int np, ne, nseg, i, j;

        ifstream in (filename);

        in >> buf;

        in >> nseg;
        for (i = 1; i <= nseg; i++)
          {
            int bound, p1, p2;
            in >> bound >> p1 >> p2;
            // forget them
          }

        in >> ne;
        for (i = 1; i <= ne; i++)
          {
            int mat, nelp;
            in >> mat >> nelp;
            Element2d el (nelp == 3 ? TRIG : QUAD);
            el.SetIndex (mat);
            for (j = 1; j <= nelp; j++)
              in >> el.PNum(j);
            mesh.AddSurfaceElement (el);
          }

        in >> np;
        for (i = 1; i <= np; i++)
          {
            Point3d p(0,0,0);
            in >> p.X() >> p.Y();
            mesh.AddPoint (p);
          }
      }

  
    else if ( (strlen (filename) > 5) &&
              strcmp (&filename[strlen (filename)-5], ".mesh") == 0 )
      {
        cout << "Reading Neutral Format" << endl;
      
        int np, ne, nse, i, j;

        ifstream in (filename);

        in >> np;

        if (in.good())
          {
            // file starts with an integer

            for (i = 1; i <= np; i++)
              {
                Point3d p(0,0,0);
                in >> p.X() >> p.Y() >> p.Z();
                mesh.AddPoint (p);
              }
	  
            in >> ne;
            for (i = 1; i <= ne; i++)
              {
                int mat;
                in >> mat;
                Element el (4);
                el.SetIndex (mat);
                for (j = 1; j <= 4; j++)
                  in >> el.PNum(j);
                mesh.AddVolumeElement (el);
              }

            mesh.AddFaceDescriptor (FaceDescriptor (1, 1, 0, 0));
	  
            in >> nse;
            for (i = 1; i <= nse; i++)
              {
                int mat, nelp;
                in >> mat;
                Element2d el (TRIG);
                el.SetIndex (mat);
                for (j = 1; j <= 3; j++)
                  in >> el.PNum(j);
                mesh.AddSurfaceElement (el);
              }
          }
        else
          {
            char buf[100];
            in.clear();
            do
              {
                in >> buf;
                cout << "buf = " << buf << endl;
                if (strcmp (buf, "points") == 0)
                  {
                    in >> np;
                    cout << "np = " << np << endl;
                  }
              }
            while (in.good());
          }
      }


    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".emt") == 0 )
      {
        ifstream inemt (filename);
      
        string pktfile = filename;
        int len = strlen (filename);
        pktfile[len-3] = 'p';
        pktfile[len-2] = 'k';
        pktfile[len-1] = 't';
        cout << "pktfile = " << pktfile << endl;

        int np, nse, i;
        int num, bcprop;
        ifstream inpkt (pktfile.c_str());
        inpkt >> np;
        ARRAY<double> values(np);
        for (i = 1; i <= np; i++)
          {
            Point3d p(0,0,0);
            inpkt >> p.X() >> p.Y() >> p.Z()
                  >> bcprop >> values.Elem(i);
            mesh.AddPoint (p);
          }      

        mesh.ClearFaceDescriptors();
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(1).SetBCProperty (1);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(2).SetBCProperty (2);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(3).SetBCProperty (3);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(4).SetBCProperty (4);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(5).SetBCProperty (5);

        int p1, p2, p3;
        double value;
        inemt >> nse;
        for (i = 1; i <= nse; i++)
          {
            inemt >> p1 >> p2 >> p3 >> bcprop >> value;

            if (bcprop < 1 || bcprop > 4)
              cerr << "bcprop out of range, bcprop = " << bcprop << endl;
            p1++;
            p2++;
            p3++;
            if (p1 < 1 || p1 > np || p2 < 1 || p2 > np || p3 < 1 || p3 > np)
              {
                cout << "p1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << endl;
              }

            if (i > 110354) Swap (p2, p3);
            if (mesh.Point(p1)(0) < 0.25)
              Swap (p2,p3);

            Element2d el(TRIG);

            if (bcprop == 1)
              {
                if (values.Get(p1) < -69999)
                  el.SetIndex(1);
                else
                  el.SetIndex(2);
              }
            else
              el.SetIndex(3);


            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            mesh.AddSurfaceElement (el);
          }


        ifstream incyl ("ngusers/guenter/cylinder.surf");
        int npcyl, nsecyl; 
        incyl >> npcyl;
        cout << "npcyl = " << npcyl << endl;
        for (i = 1; i <= npcyl; i++)
          {
            Point3d p(0,0,0);
            incyl >> p.X() >> p.Y() >> p.Z();
            mesh.AddPoint (p);
          }
        incyl >> nsecyl;
        cout << "nsecyl = " << nsecyl << endl;
        for (i = 1; i <= nsecyl; i++)
          {
            incyl >> p1 >> p2 >> p3;
            p1 += np;
            p2 += np;
            p3 += np;
            Element2d el(TRIG);
            el.SetIndex(5);
            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            mesh.AddSurfaceElement (el);
          }
      }


    // .tet mesh
    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".tet") == 0 )
      {
        ReadTETFormat (mesh, filename);
      }
  }
  
}

