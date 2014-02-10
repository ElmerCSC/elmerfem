#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{
  GeomSearch3d :: GeomSearch3d() 
  {
    size.i1 = 0; size.i2 = 0; size.i3 = 0; 
  };

  GeomSearch3d :: ~GeomSearch3d()
  {
    //delete old Hashtable:
    if (size.i1 != 0)
      {
	for (int i = 0; i < size.i1*size.i2*size.i3; i++)
	  delete hashtable[i];
      } 
  }

  void GeomSearch3d :: Init (ARRAY <FrontPoint3,PointIndex::BASE> *pointsi, ARRAY <FrontFace> *facesi)
  {
    points = pointsi;
    faces = facesi;
    size.i1 = 0; size.i2 = 0; size.i3 = 0; 
    reset = 1;
    hashcount = 1;
  }

  void GeomSearch3d :: ElemMaxExt(Point3d& minp, Point3d& maxp, const MiniElement2d& elem)
  {
    maxp.X()=(*points)[elem.PNum(1)].P()(0);
    maxp.Y()=(*points)[elem.PNum(1)].P()(1);
    maxp.Z()=(*points)[elem.PNum(1)].P()(2);
    minp.X()=(*points)[elem.PNum(1)].P()(0);
    minp.Y()=(*points)[elem.PNum(1)].P()(1);
    minp.Z()=(*points)[elem.PNum(1)].P()(2);
  
    for (int i=2; i <= 3; i++)
      {
	maxp.X()=max2((*points)[elem.PNum(i)].P()(0),maxp.X());
	maxp.Y()=max2((*points)[elem.PNum(i)].P()(1),maxp.Y());
	maxp.Z()=max2((*points)[elem.PNum(i)].P()(2),maxp.Z());
	minp.X()=min2((*points)[elem.PNum(i)].P()(0),minp.X());
	minp.Y()=min2((*points)[elem.PNum(i)].P()(1),minp.Y());
	minp.Z()=min2((*points)[elem.PNum(i)].P()(2),minp.Z());
      }
  }

  void GeomSearch3d :: MinCoords(const Point3d& p1, Point3d& p2)
  {
    p2.X()=min2(p1.X(),p2.X());
    p2.Y()=min2(p1.Y(),p2.Y());
    p2.Z()=min2(p1.Z(),p2.Z());
  }

  void GeomSearch3d :: MaxCoords(const Point3d& p1, Point3d& p2)
  {
    p2.X()=max2(p1.X(),p2.X());
    p2.Y()=max2(p1.Y(),p2.Y());
    p2.Z()=max2(p1.Z(),p2.Z());
  }

  void GeomSearch3d :: Create()
  {
    INDEX i,j,k;
    if (reset)
      {
	const double hashelemsizefactor = 4;
	reset = 0;
	/*
	  minext=Point3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
	  maxext=Point3d(MINDOUBLE, MINDOUBLE, MINDOUBLE);
	*/
	ElemMaxExt(minext, maxext, faces->Get(1).Face());
	Point3d maxp, minp;
	Vec3d midext(0,0,0);
      
	//get max Extension of Frontfaces
	for (i = 1; i <= faces->Size(); i++)
	  {
	    ElemMaxExt(minp, maxp, faces->Get(i).Face());
	    MinCoords(minp, minext);
	    MaxCoords(maxp, maxext);
	    midext+=maxp-minp;
	  }


	maxextreal = maxext;
	maxext = maxext + 1e-4 * (maxext - minext);

	midext*=1./faces->Size();
	Vec3d boxext = maxext - minext;
      
	//delete old Hashtable:
	if (size.i1 != 0)
	  {
	    for (i = 1; i <= size.i1*size.i2*size.i3; i++)
	      {
		delete hashtable.Get(i);
	      }
	  } 
      
	size.i1 = int (boxext.X()/midext.X()/hashelemsizefactor+1);
	size.i2 = int (boxext.Y()/midext.Y()/hashelemsizefactor+1);
	size.i3 = int (boxext.Z()/midext.Z()/hashelemsizefactor+1);
	// PrintMessage (5, "hashsizes = ", size.i1, ", ", size.i2, ", ", size.i3);
      
	elemsize.X()=boxext.X()/size.i1;
	elemsize.Y()=boxext.Y()/size.i2;
	elemsize.Z()=boxext.Z()/size.i3;

	//create Hasharrays:
	hashtable.SetSize(size.i1*size.i2*size.i3);
	for (i = 1; i <= size.i1; i++)
	  {
	    for (j = 1; j <= size.i2; j++)
	      {
		for (k = 1; k <= size.i3; k++)
		  {
		    INDEX ind=i+(j-1)*size.i1+(k-1)*size.i2*size.i1;
		    hashtable.Elem(ind) = new ARRAY <int> ();
		  }
	      }
	  }
      }
    else
      {
	//Clear all Hash-Arrays
	for (i = 1; i <= size.i1; i++)
	  {
	    for (j = 1; j <= size.i2; j++)
	      {
		for (k = 1; k <= size.i3; k++)
		  {
		    INDEX ind=i+(j-1)*size.i1+(k-1)*size.i2*size.i1;
		    hashtable.Elem(ind)->SetSize(0);
		  }
	      }
	  }	  
      }
  
    //Faces in Hashtable einfuegen:
    for (i = 1; i <= faces->Size(); i++)
      {
	AddElem(faces->Get(i).Face(),i);
      }
  
  }

  void GeomSearch3d :: AddElem(const MiniElement2d& elem, INDEX elemnum)
  {
    Point3d minp, maxp;
    ElemMaxExt(minp, maxp, elem);
    int sx = int ((minp.X()-minext.X())/elemsize.X()+1.);
    int ex = int ((maxp.X()-minext.X())/elemsize.X()+1.);
    int sy = int ((minp.Y()-minext.Y())/elemsize.Y()+1.);
    int ey = int ((maxp.Y()-minext.Y())/elemsize.Y()+1.);
    int sz = int ((minp.Z()-minext.Z())/elemsize.Z()+1.);
    int ez = int ((maxp.Z()-minext.Z())/elemsize.Z()+1.);
  
    for (int ix = sx; ix <= ex; ix++)
      for (int iy = sy; iy <= ey; iy++)
        for (int iz = sz; iz <= ez; iz++)
          {
            INDEX ind=ix+(iy-1)*size.i1+(iz-1)*size.i2*size.i1;
            if (ind < 1 || ind > size.i1 * size.i2 * size.i3)
              {
                cerr << "Illegal hash-position";
                cerr << "Position: " << ix << "," << iy << "," << iz << endl;
		    throw NgException ("Illegal position in Geomsearch");
              }
            hashtable.Elem(ind)->Append(elemnum);		      
          }
  }

  void GeomSearch3d :: GetLocals(ARRAY<MiniElement2d> & locfaces,  ARRAY<INDEX> & findex,
				 INDEX fstind, const Point3d& p0, double xh)
  {
    hashcount++;
  
    Point3d minp, maxp, midp; 

    minp=p0-Vec3d(xh,xh,xh); //lay cube over sphere
    maxp=p0+Vec3d(xh,xh,xh);

    MaxCoords(minext,minp); //cube may not be out of hash-region
    MinCoords(maxextreal,maxp);


    int cluster = faces->Get(fstind).Cluster();
  
    int sx = int((minp.X()-minext.X())/elemsize.X()+1.);
    int ex = int((maxp.X()-minext.X())/elemsize.X()+1.);
    int sy = int((minp.Y()-minext.Y())/elemsize.Y()+1.);
    int ey = int((maxp.Y()-minext.Y())/elemsize.Y()+1.);
    int sz = int((minp.Z()-minext.Z())/elemsize.Z()+1.);
    int ez = int((maxp.Z()-minext.Z())/elemsize.Z()+1.);
    int ix,iy,iz,i,k;

    int cnt1 = 0;  // test, how efficient hashtable is
    int cnt2 = 0;
    int cnt3 = 0;
  
    for (ix = sx; ix <= ex; ix++)
      {
	for (iy = sy; iy <= ey; iy++)
	  {
	    for (iz = sz; iz <= ez; iz++)
	      {
		INDEX ind=ix+(iy-1)*size.i1+(iz-1)*size.i2*size.i1;
	      
		//go through all elements in one hash area
		const ARRAY <int> & area = *hashtable.Elem(ind);
		for (k = 1; k <= area.Size(); k++)
		  {
		    cnt2++;
		    i = area.Get(k);
		    if (faces->Get(i).Cluster() == cluster && 
			faces->Get(i).Valid() &&
			faces->Get(i).HashValue() != hashcount && 
			i != fstind)
		      {
			cnt1++;
			const MiniElement2d & face = faces->Get(i).Face();
		      
			const Point3d & p1 = (*points)[face.PNum(1)].P();
			const Point3d & p2 = (*points)[face.PNum(2)].P();
			const Point3d & p3 = (*points)[face.PNum(3)].P();
		      
			midp = Center (p1, p2, p3);
		      
			// if (Dist2 (midp, p0) <= xh*xh)  
                        if((Dist2 (p1, p0) <= xh*xh) ||
                           (Dist2 (p2, p0) <= xh*xh) ||
                           (Dist2 (p3, p0) <= xh*xh) ||
                           (Dist2 (midp, p0) <= xh*xh) )  // by Jochen Wild
			  {
			    cnt3++;
			    locfaces.Append(faces->Get(i).Face());
			    findex.Append(i);
			    faces->Elem(i).SetHashValue(hashcount);
			  }
		      }
		  }
	      }
	  }
      }
    /*
      if (faces->Size() != 0 && hashcount % 200 == 0)
      {
      (*mycout) << "n.o.f= " << faces->Size();
      (*mycout) << ", n.o.lf= " << locfaces.Size();
      (*mycout) << ", hashf= " << (double)cnt2/(double)faces->Size();
      (*mycout) << " (" << (double)cnt1/(double)faces->Size();
      (*mycout) << ", " << (double)cnt3/(double)faces->Size() << ")" << endl;
      }
    */

  }

}
