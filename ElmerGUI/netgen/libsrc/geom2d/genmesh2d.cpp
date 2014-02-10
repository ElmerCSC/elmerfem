#include <mystdlib.h>
#include <csg.hpp>
#include <geometry2d.hpp>
#include "meshing.hpp"

namespace netgen
{

  // static ARRAY<Point<2> > points2;
  //  static ARRAY<int> lp1, lp2;


  extern void Optimize2d (Mesh & mesh, MeshingParameters & mp);







  void MeshFromSpline2D (SplineGeometry2d & geometry,
			 Mesh *& mesh, 
			 MeshingParameters & mp)
  {
    PrintMessage (1, "Generate Mesh from spline geometry");

    double h = mp.maxh;

    Box<2> bbox = geometry.GetBoundingBox ();

    if (bbox.Diam() < h) 
      {
	h = bbox.Diam();
	mp.maxh = h;
      }

    mesh = new Mesh;
    mesh->SetDimension (2);

    geometry.PartitionBoundary (h, *mesh);

    // marks mesh points for hp-refinement
    for (int i = 0; i < geometry.GetNP(); i++)
      if (geometry.GetPoint(i).hpref)
	{
	  double mindist = 1e99;
	  PointIndex mpi(0);
	  Point<2> gp = geometry.GetPoint(i);
	  Point<3> gp3(gp(0), gp(1), 0);
	  for (PointIndex pi = PointIndex::BASE; 
	       pi < mesh->GetNP()+PointIndex::BASE; pi++)
	    if (Dist2(gp3, (*mesh)[pi]) < mindist)
	      {
		mpi = pi;
		mindist = Dist2(gp3, (*mesh)[pi]);
	      }
	  (*mesh)[mpi].Singularity(1.);
	}


    int maxdomnr = 0;
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
	if ( (*mesh)[si].domin > maxdomnr) maxdomnr = (*mesh)[si].domin;
	if ( (*mesh)[si].domout > maxdomnr) maxdomnr = (*mesh)[si].domout;
      }

    mesh->ClearFaceDescriptors();
    for (int i = 1; i <= maxdomnr; i++)
      mesh->AddFaceDescriptor (FaceDescriptor (i, 0, 0, i));

    // set ARRAY<string*> bcnames... 
    // number of bcnames
    int maxsegmentindex = 0;
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
	if ( (*mesh)[si].si > maxsegmentindex) maxsegmentindex = (*mesh)[si].si;
      }

    mesh->SetNBCNames(maxsegmentindex);

    for ( int sindex = 0; sindex < maxsegmentindex; sindex++ )
      {
	mesh->SetBCName ( sindex, geometry.GetBCName( sindex+1 ) );
      }

    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
 	(*mesh)[si].SetBCName ( (*mesh).GetBCNamePtr( (*mesh)[si].si-1 ) );
      }
    Point3d pmin(bbox.PMin()(0), bbox.PMin()(1), -bbox.Diam());
    Point3d pmax(bbox.PMax()(0), bbox.PMax()(1), bbox.Diam());

    mesh->SetLocalH (pmin, pmax, mparam.grading);
    mesh->SetGlobalH (h);
  
    mesh->CalcLocalH();

    int bnp = mesh->GetNP(); // boundary points

    int hquad = mparam.quad;


    for (int domnr = 1; domnr <= maxdomnr; domnr++)
      if (geometry.GetDomainTensorMeshing (domnr))
        { // tensor product mesh
          
          ARRAY<PointIndex, PointIndex::BASE> nextpi(bnp);
          ARRAY<int, PointIndex::BASE> si1(bnp), si2(bnp);
          PointIndex firstpi;
          
          nextpi = -1;
          si1 = -1;
          si2 = -1;
          for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
            {
              int p1 = -1, p2 = -2;

              if ( (*mesh)[si].domin == domnr)
                { p1 = (*mesh)[si].p1; p2 = (*mesh)[si].p2; }
              if ( (*mesh)[si].domout == domnr)
                { p1 = (*mesh)[si].p2; p2 = (*mesh)[si].p1; }
              
              if (p1 == -1) continue;

              nextpi[p1] = p2;       // counter-clockwise
              
              int index = (*mesh)[si].si;
              if (si1[p1] != index && si2[p1] != index)
                { si2[p1] = si1[p1]; si1[p1] = index; }
              if (si1[p2] != index && si2[p2] != index)
                { si2[p2] = si1[p2]; si1[p2] = index; }
            }

          PointIndex c1(0), c2, c3, c4;  // 4 corner points
          int nex = 1, ney = 1;

          for (PointIndex pi = 1; pi <= si2.Size(); pi++)
            if (si2[pi] != -1)
              { c1 = pi; break; }      

          for (c2 = nextpi[c1]; si2[c2] == -1; c2 = nextpi[c2], nex++);
          for (c3 = nextpi[c2]; si2[c3] == -1; c3 = nextpi[c3], ney++);
          for (c4 = nextpi[c3]; si2[c4] == -1; c4 = nextpi[c4]);



          ARRAY<PointIndex> pts ( (nex+1) * (ney+1) );   // x ... inner loop
          pts = -1;

          for (PointIndex pi = c1, i = 0; pi != c2; pi = nextpi[pi], i++)
            pts[i] = pi;
          for (PointIndex pi = c2, i = 0; pi != c3; pi = nextpi[pi], i++)
            pts[(nex+1)*i+nex] = pi;
          for (PointIndex pi = c3, i = 0; pi != c4; pi = nextpi[pi], i++)
            pts[(nex+1)*(ney+1)-i-1] = pi;
          for (PointIndex pi = c4, i = 0; pi != c1; pi = nextpi[pi], i++)
            pts[(nex+1)*(ney-i)] = pi;


          for (PointIndex pix = nextpi[c1], ix = 0; pix != c2; pix = nextpi[pix], ix++)
            for (PointIndex piy = nextpi[c2], iy = 0; piy != c3; piy = nextpi[piy], iy++)
              {
                Point<3> p = (*mesh)[pix] + ( (*mesh)[piy] - (*mesh)[c2] );
                pts[(nex+1)*(iy+1) + ix+1] = mesh -> AddPoint (p , 1, FIXEDPOINT);
              }

          for (int i = 0; i < ney; i++)
            for (int j = 0; j < nex; j++)
              {
                Element2d el(QUAD);
                el[0] = pts[i*(nex+1)+j];
                el[1] = pts[i*(nex+1)+j+1];
                el[2] = pts[(i+1)*(nex+1)+j+1];
                el[3] = pts[(i+1)*(nex+1)+j];
                el.SetIndex (domnr);

                mesh -> AddSurfaceElement (el);
              }
        }




    for (int domnr = 1; domnr <= maxdomnr; domnr++)
      {
        if (geometry.GetDomainTensorMeshing (domnr)) continue;

	if ( geometry.GetDomainMaxh ( domnr ) > 0 )
	  h = geometry.GetDomainMaxh(domnr);


	PrintMessage (3, "Meshing domain ", domnr, " / ", maxdomnr);

	int oldnf = mesh->GetNSE();

        mparam.quad = hquad || geometry.GetDomainQuadMeshing (domnr);

	Meshing2 meshing (Box<3> (pmin, pmax));

	for (PointIndex pi = PointIndex::BASE; pi < bnp+PointIndex::BASE; pi++)
	  meshing.AddPoint ( (*mesh)[pi], pi);
      

	PointGeomInfo gi;
	gi.trignum = 1;
	for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
	  {
	    if ( (*mesh)[si].domin == domnr)
	      meshing.AddBoundaryElement ( (*mesh)[si].p1 + 1 - PointIndex::BASE, 
					   (*mesh)[si].p2 + 1 - PointIndex::BASE, gi, gi);
	    if ( (*mesh)[si].domout == domnr)
	      meshing.AddBoundaryElement ( (*mesh)[si].p2 + 1 - PointIndex::BASE, 
					   (*mesh)[si].p1 + 1 - PointIndex::BASE, gi, gi);
	  }


	mparam.checkoverlap = 0;

	meshing.GenerateMesh (*mesh, h, domnr);

	for (SurfaceElementIndex sei = oldnf; sei < mesh->GetNSE(); sei++)
	  (*mesh)[sei].SetIndex (domnr);


	// astrid
	char * material;
	geometry.GetMaterial( domnr, material );
	if ( material )
	  {
	    (*mesh).SetMaterial ( domnr,  material );
	  }

      }

    mparam.quad = hquad;



    int hsteps = mp.optsteps2d;

    mp.optimize2d = "smcm"; 
    mp.optsteps2d = hsteps/2;
    Optimize2d (*mesh, mp);

    mp.optimize2d = "Smcm"; 
    mp.optsteps2d = (hsteps+1)/2;
    Optimize2d (*mesh, mp);

    mp.optsteps2d = hsteps;

    mesh->Compress();
    mesh -> SetNextMajorTimeStamp();


#ifdef OPENGL
    extern void Render();
    Render();
#endif

  }








}
