// Write Chemnitz file format


#include <mystdlib.h>

#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{

class POINT3D
  {
  public:
  POINT3D () { };
  double x, y, z;
  };

class VOLELEMENT
  {
  public:
  VOLELEMENT () {};
  int domnr, p1, p2, p3, p4;
  int faces[4];
  };
  
class SURFELEMENT
  {
  public:
  SURFELEMENT () { };
  int snr, p1, p2, p3;
  };
  

class FACE
  {
  public:
  FACE () { };
  int p1, p2, p3;
  int edges[3];
  };

class EDGE
  {
  public:
  EDGE () { };
  int p1, p2;
  };

static ARRAY<POINT3D> points;
static ARRAY<VOLELEMENT> volelements;
static ARRAY<SURFELEMENT> surfelements;

static ARRAY<FACE> faces;
static ARRAY<EDGE> edges;


void ReadFile (char * filename)
  {
  int i, n;
  ifstream infile(filename);
  char reco[100];
  
  
  infile >> reco;  // file format recognition
  
  infile >> n;   // number of surface elements
  cout << n << " Surface elements" << endl;
  
  for (i = 1; i <= n; i++)
    {
    SURFELEMENT sel;
    infile >> sel.snr >> sel.p1 >> sel.p2 >> sel.p3;
    surfelements.Append (sel);
    }
    
  infile >> n;   // number of volume elements
  cout << n << " Volume elements" << endl;
  
  for (i = 1; i <= n; i++)
    {
    VOLELEMENT el;
    infile >> el.p1 >> el.p2 >> el.p3 >> el.p4;
    volelements.Append (el);
    }
    
  infile >> n;   // number of points 
  cout << n << " Points" << endl;
  
  for (i = 1; i <= n; i++)
    {
    POINT3D p;
    infile >> p.x >> p.y >> p.z;
    points.Append (p);
    }
  }
  
  

void ReadFileMesh (const Mesh & mesh)
{
  int i, n;
  
  n = mesh.GetNSE();   // number of surface elements
  cout << n << " Surface elements" << endl;
  
  for (i = 1; i <= n; i++)
    {
      SURFELEMENT sel;
      const Element2d & el = mesh.SurfaceElement(i);
      sel.snr = el.GetIndex();
      sel.p1 = el.PNum(1);
      sel.p2 = el.PNum(2);
      sel.p3 = el.PNum(3);
      surfelements.Append (sel);
    }
    
  n = mesh.GetNE();   // number of volume elements
  cout << n << " Volume elements" << endl;
  
  for (i = 1; i <= n; i++)
    {
      VOLELEMENT el;
      const Element & nel = mesh.VolumeElement(i);
      el.p1 = nel.PNum(1);
      el.p2 = nel.PNum(2);
      el.p3 = nel.PNum(3);
      el.p4 = nel.PNum(4);
      //      infile >> el.p1 >> el.p2 >> el.p3 >> el.p4;
      volelements.Append (el);
    }
    
  n = mesh.GetNP();   // number of points 
  cout << n << " Points" << endl;
  
  for (i = 1; i <= n; i++)
    {
      POINT3D p;
      Point3d mp = mesh.Point(i);
      p.x = mp.X();
      p.y = mp.Y();
      p.z = mp.Z();
      //      infile >> p.x >> p.y >> p.z;
      points.Append (p);
    }
  }
  



void Convert ()
  {
  int i, j, facei, edgei;
  INDEX_3 i3;
  INDEX_2 i2;

  INDEX_3_HASHTABLE<int> faceindex(volelements.Size()/5 + 1);
  INDEX_2_HASHTABLE<int> edgeindex(volelements.Size()/5 + 1);
  
  for (i = 1; i <= volelements.Size(); i++)
    {
    for (j = 1; j <= 4; j++)
      {
      switch (j)
        {
        case 1:
          i3.I1() = volelements.Get(i).p2;
          i3.I2() = volelements.Get(i).p3;
          i3.I3() = volelements.Get(i).p4;
          break;
        case 2:
          i3.I1() = volelements.Get(i).p1;
          i3.I2() = volelements.Get(i).p3;
          i3.I3() = volelements.Get(i).p4;
          break;
         case 3:
          i3.I1() = volelements.Get(i).p1;
          i3.I2() = volelements.Get(i).p2;
          i3.I3() = volelements.Get(i).p4;
          break;
         case 4:
          i3.I1() = volelements.Get(i).p1;
          i3.I2() = volelements.Get(i).p2;
          i3.I3() = volelements.Get(i).p3;
          break;
		 default:
			 i3.I1()=i3.I2()=i3.I3()=0;
        }
      i3.Sort();
      if (faceindex.Used (i3)) 
        facei = faceindex.Get(i3);
      else
        {
        FACE fa;
        fa.p1 = i3.I1();
        fa.p2 = i3.I2();
        fa.p3 = i3.I3();
        facei = faces.Append (fa);
        faceindex.Set (i3, facei);
        } 
        
      volelements.Elem(i).faces[j-1] = facei;  
      }    
    
    } 
 

  for (i = 1; i <= faces.Size(); i++)
    {
    for (j = 1; j <= 3; j++)
      {
      switch (j)
        {
        case 1:
          i2.I1() = faces.Get(i).p2;
          i2.I2() = faces.Get(i).p3;
          break;
        case 2:
          i2.I1() = faces.Get(i).p1;
          i2.I2() = faces.Get(i).p3;
          break;
         case 3:
          i2.I1() = faces.Get(i).p1;
          i2.I2() = faces.Get(i).p2;
          break;
		 default:
			 i2.I1()=i2.I2()=0;
        }
      if (i2.I1() > i2.I2()) swap (i2.I1(), i2.I2());
      if (edgeindex.Used (i2)) 
        edgei = edgeindex.Get(i2);
      else
        {
        EDGE ed;
        ed.p1 = i2.I1();
        ed.p2 = i2.I2();
        edgei = edges.Append (ed);
        edgeindex.Set (i2, edgei);
        } 
        
      faces.Elem(i).edges[j-1] = edgei;  
      }    
    
    }  
 
  }  
  
  
void WriteFile (ostream & outfile)
  {
  int i;
  
  outfile 
  	<< "#VERSION: 1.0" << endl
  	<< "#PROGRAM: NETGEN" << endl
  	<< "#EQN_TYPE: POISSON" << endl
  	<< "#DIMENSION: 3D" << endl
  	<< "#DEG_OF_FREE: 1" << endl
  	<< "#DESCRIPTION: I don't know" << endl
  	<< "##RENUM: not done" << endl
  	<< "#USER: Kleinzen" << endl
  	<< "DATE: 10.06.1996" << endl;
  
  outfile << "#HEADER:   8" << endl
  	<< points.Size() << "  " << edges.Size() << "  " 
  	<< faces.Size() << "  " << volelements.Size() << "  0  0  0  0" << endl;
  
  outfile << "#VERTEX:   " << points.Size() << endl;
  for (i = 1; i <= points.Size(); i++)
    outfile << "  " << i << "  " << points.Get(i).x << "  " << points.Get(i).y 
    	<< "  " << points.Get(i).z << endl;
    	
  outfile << "#EDGE:  " << edges.Size() << endl;
  for (i = 1; i <= edges.Size(); i++)
    outfile << "  " << i << "  1  " 
    	<< edges.Get(i).p1 << "  " 
    	<< edges.Get(i).p2 
    	<< "  0" << endl;
    
  outfile << "#FACE:  " << faces.Size() << endl;  
  for (i = 1; i <= faces.Size(); i++)
    outfile << "  " << i << "  1  3  " 
    	<< faces.Get(i).edges[0] << "  " 
    	<< faces.Get(i).edges[1] << "  " 
    	<< faces.Get(i).edges[2] << endl;
    	
  outfile << "#SOLID:  " << volelements.Size() << endl;
  for (i = 1; i <= volelements.Size(); i++)
    outfile << "  " << i << "  1  4  " 
    	<< volelements.Get(i).faces[0] << "  "
    	<< volelements.Get(i).faces[1] << "  "
    	<< volelements.Get(i).faces[2] << "  "
    	<< volelements.Get(i).faces[3] << endl;
    	
  outfile << "#END_OF_DATA" << endl;
  }
    

void WriteUserChemnitz (const Mesh & mesh,
			const string & filename)
{
  ofstream outfile (filename.c_str());

  ReadFileMesh (mesh);
  Convert ();
  
  WriteFile (outfile);
  cout << "Wrote Chemnitz standard file" << endl;
}
}
