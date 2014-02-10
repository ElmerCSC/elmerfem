
#include <mystdlib.h>

#include "meshing.hpp"
#ifndef CURVEDELEMS_NEW

namespace netgen
{
    

// ----------------------------------------------------------------------------
//      CurvedElements
// ----------------------------------------------------------------------------

    CurvedElements :: CurvedElements (const Mesh & amesh)
	: mesh(amesh), top(mesh.GetTopology())
    {
	isHighOrder = 0;
	nvisualsubsecs = 2;
	nIntegrationPoints = 10;
    }


    CurvedElements :: ~CurvedElements ()
    {
      ;
    }


  void CurvedElements :: BuildCurvedElements(Refinement * ref, int polydeg, bool rational)
    {
      if (mesh.coarsemesh)
	{
	  mesh.coarsemesh->GetCurvedElements().BuildCurvedElements (ref, polydeg, rational);
	  SetHighOrder();
	  return;
	}

      PrintMessage (2, "Build curved elements, order = ", polydeg);

      NgLock lock(const_cast<Mesh&>(mesh).Mutex(), 1);
      isHighOrder = 0;
      lock.UnLock();

	const_cast<Mesh &>(mesh).UpdateTopology();

	// set order of edges and faces

	BaseFiniteElement2D * fe2d;

	FEEdge edge (*this);
	FESegm segm (*this);
	FETrig trig (*this);
	FEQuad quad (*this);

	int i, k, e, f;

	ARRAY<bool> edgedone;

	edgedone.SetSize (top.GetNEdges());

	edgeorder.SetSize (top.GetNEdges());
	faceorder.SetSize (top.GetNFaces());

	int nedgestocurve = mesh.GetNSeg();

	edgedone = 0;
	edgeorder = 1;
	faceorder = 1;

	if (polydeg == 1)
	{
	    isHighOrder = 0;
	    return;
	}

	
	/*
	for (e = 0; e < top.GetNEdges(); e++)
	  {
	    edgedone = 0;
	    edgeorder[e] = 1;
	  }

	for (f = 0; f < top.GetNFaces(); f++)
	    faceorder[f] = 1;
	*/

	for (i = 1; i <= mesh.GetNSeg(); i++) 
	    edgeorder[top.GetSegmentEdge(i)-1] = polydeg;


	if (mesh.GetDimension() == 3)
	  {
	    for (i = 1; i <= mesh.GetNSE(); i++)
	      {
		faceorder[top.GetSurfaceElementFace(i)-1] = polydeg;
		
		Element2d elem = mesh[(SurfaceElementIndex) (i-1)];
		
		ARRAY<int> edgenrs;
		top.GetSurfaceElementEdges(i, edgenrs);
		
		nedgestocurve += top.GetNEdges(elem.GetType());
		
		for (int e = 0; e < top.GetNEdges(elem.GetType()); e++)
		  edgeorder[edgenrs[e]-1] = polydeg;
	      }
	  }


	PrintMessage (1, "Building curved elements, order = ", polydeg);
	PushStatusF ("curving edges");



        // set edgecoeffs and facecoeffs arrays index and size

	edgecoeffsindex.SetSize (top.GetNEdges()+1);
	facecoeffsindex.SetSize (top.GetNFaces()+1);

	edgecoeffsindex[0] = 0;
	for (e = 2; e <= top.GetNEdges()+1; e++)
	    edgecoeffsindex[e-1] = edgecoeffsindex[e-2] + edgeorder[e-2]-1;

	facecoeffsindex[0] = 0;
	for (f = 2; f <= top.GetNFaces()+1; f++)
	{
	    switch (top.GetFaceType (f-1))
	    {
		case TRIG:
		    facecoeffsindex[f-1] = facecoeffsindex[f-2] + 
			(faceorder[f-2]-1)*(faceorder[f-2]-2)/2;
		    break;
		case QUAD:
		    facecoeffsindex[f-1] = facecoeffsindex[f-2] +
			(faceorder[f-2]-1)*(faceorder[f-2]-1);
		    break;
	    }
	}

	edgecoeffs.SetSize(edgecoeffsindex[top.GetNEdges()]);
	facecoeffs.SetSize(facecoeffsindex[top.GetNFaces()]);


        
	// evaluate edge points

	PointGeomInfo newgi;          // dummy variable, only needed for function call
	EdgePointGeomInfo newepgi;    // dummy variable, only needed for function call
	Point3d xexact;               // new point to be stored in ARRAY edgepts

	ARRAY<double> xi, wi;
	ComputeGaussRule(nIntegrationPoints, xi, wi);

	for (i=0; i<edgecoeffsindex[top.GetNEdges()]; i++)
	    edgecoeffs[i] = Vec<3>(0.,0.,0.);




	// all edges belonging to segments

	for (i=0; i<mesh.GetNSeg(); i++) 
	{
	  if (multithread.terminate) return;

	  SetThreadPercent( double(100*i/nedgestocurve) );

	  int edgenr = top.GetSegmentEdge(i+1);

	  if (edgedone[edgenr-1]) continue;

	  edgedone[edgenr-1] = 1;

            Segment s = mesh.LineSegment(i+1); 

	    segm.SetElementNumber (i+1);
	
	    for (k = 2; k <= segm.GetEdgeOrder(); k++)
	      edgecoeffs[edgecoeffsindex[edgenr-1]+k-2] = Vec<3>(0.,0.,0.);

	    for (int l = 0; l < nIntegrationPoints; l++)
	      {
		segm.SetReferencePoint (Point<1>(xi[l]));
		segm.CalcVertexShapes ();
		segm.CalcEdgeLaplaceShapes ();
		
		Point<3> xv(0,0,0);

		for (int v = 0; v < 2; v++)
		  xv = xv + segm.GetVertexShape(v) * mesh.Point(segm.GetVertexNr(v));
		  		
		double secpoint = xi[l];

		ref->PointBetween (mesh.Point(segm.GetVertexNr(1)),
				   mesh.Point(segm.GetVertexNr(0)), secpoint,
				   s.surfnr2, s.surfnr1,
				   s.epgeominfo[1], s.epgeominfo[0],
				   xexact, newepgi);
		
		for (int k = 2; k <= segm.GetEdgeOrder(); k++)
		  edgecoeffs[edgecoeffsindex[edgenr-1]+k-2] -=
		    wi[l] * segm.GetEdgeLaplaceShape(k-2) * Vec<3>(xexact - xv);
		
	      }
	    
	    for (k = 2; k <= segm.GetEdgeOrder(); k++)
	      edgecoeffs[edgecoeffsindex[edgenr-1]+k-2] =
		(2.0*(k-1.0)+1.0)*edgecoeffs[edgecoeffsindex[edgenr-1]+k-2];
	
	}





	// all edges belonging to surface elements
	
	if (mesh.GetDimension() == 3)
	  {
	    int nedgescurved = mesh.GetNSeg();
	    for (int i=0; i<mesh.GetNSE(); i++) 
	      {
		if (multithread.terminate) return;
		
		//		SetThreadPercent( double(100*(mesh.GetNSeg()+i)/nedgestocurve) );
		Element2d elem = mesh[(SurfaceElementIndex) i];
		const ELEMENT_EDGE * eledges = MeshTopology::GetEdges(elem.GetType());
		
		ARRAY<int> edgenrs;
		ARRAY<int> orient;
		top.GetSurfaceElementEdges(i+1, edgenrs);
		top.GetSurfaceElementEdgeOrientations(i+1, orient);

		for (int e = 0; e < top.GetNEdges(elem.GetType()); e++)
		  {
//		    cout << "e = " << e << "/" << top.GetNEdges(elem.GetType()) <<  endl;

		    nedgescurved++;

		    if (edgedone[edgenrs[e]-1]) continue;
		    
		    edgedone[edgenrs[e]-1] = 1;

		    SetThreadPercent( double(100*(nedgescurved)/nedgestocurve) );

		    edge.SetElementNumber (edgenrs[e]);

		    for (k = 2; k <= edge.GetEdgeOrder(); k++)
		      edgecoeffs[edgecoeffsindex[edgenrs[e]-1]+k-2] = Vec<3>(0.,0.,0.);

		    for (int l = 0; l < nIntegrationPoints; l++)
		      {
//			cout << "." << flush;
			edge.SetReferencePoint (Point<1>(xi[l]));
			edge.CalcVertexShapes ();
			edge.CalcEdgeLaplaceShapes ();
			
			Point<3> xv(0,0,0);
			for (int v = 0; v < 2; v++)
			  xv = xv + edge.GetVertexShape(v) * mesh.Point(edge.GetVertexNr(v));

			double secpoint = xi[l];

			if (orient[e] == 1)
			  ref->PointBetween (mesh.Point(edge.GetVertexNr(1)),
					     mesh.Point(edge.GetVertexNr(0)), secpoint,
					     mesh.GetFaceDescriptor(elem.GetIndex()).SurfNr(),
					     elem.GeomInfoPi(eledges[e][1]),
					     elem.GeomInfoPi(eledges[e][0]),
					     xexact, newgi);
			else
			  ref->PointBetween (mesh.Point(edge.GetVertexNr(1)),
					     mesh.Point(edge.GetVertexNr(0)), secpoint,
					     mesh.GetFaceDescriptor(elem.GetIndex()).SurfNr(),
					     elem.GeomInfoPi(eledges[e][0]),
					     elem.GeomInfoPi(eledges[e][1]),
					     xexact, newgi);

			for (k = 2; k <= edge.GetEdgeOrder(); k++)
			  edgecoeffs[edgecoeffsindex[edgenrs[e]-1]+k-2] -=
			    wi[l] * edge.GetEdgeLaplaceShape(k-2) * Vec<3>(xexact - xv);
		      }	
//		    cout << endl;
		    for (k = 2; k <= edge.GetEdgeOrder(); k++)
		      edgecoeffs[edgecoeffsindex[edgenrs[e]-1]+k-2] =
			(2.0*(k-1.0)+1.0)*edgecoeffs[edgecoeffsindex[edgenrs[e]-1]+k-2];
		    
		}
	      }
	  }




/*

	// L2-Projection for edges


	cout << "WARNING: L2-Projection for edges" << endl;

	if (mesh.GetDimension() == 3)
	{
	    for (int i=0; i<mesh.GetNSE(); i++) 
	    {
		Element2d elem = mesh[(SurfaceElementIndex) i];
		const ELEMENT_EDGE * eledges = MeshTopology::GetEdges(elem.GetType());
		
		ARRAY<int> edgenrs;
		ARRAY<int> orient;
		top.GetSurfaceElementEdges(i+1, edgenrs);
		top.GetSurfaceElementEdgeOrientations(i+1, orient);
		
		for (int e = 0; e < top.GetNEdges(elem.GetType()); e++)
		{
		    edge.SetElementNumber (edgenrs[e]);

		    int npoints = edge.GetEdgeOrder()-1;

		    if (npoints == 0) continue;

		    DenseMatrix mat(npoints);
		    DenseMatrix inv(npoints);
		    Vector vec[3];
	    
		    for (int k = 0; k < 3; k++)
		    {
			vec[k].SetSize(npoints);
			for (int n = 1; n <= npoints; n++) vec[k].Set(n, 0.);
		    }
		    
		    for (int l = 0; l < nIntegrationPoints; l++)
		    {
			double w = wi[l];
			
			edge.SetReferencePoint (Point<1>(xi[l]));
			edge.CalcVertexShapes ();
			edge.CalcEdgeShapes ();
		
			for (int n = 0; n < npoints; n++)
			    for (int m = 0; m < npoints; m++)
				mat.Set(n+1, m+1, mat.Get(n+1,m+1) +
					edge.GetEdgeShape(n) * edge.GetEdgeShape(m) * w);
		
			Point<3> xv(0,0,0);
			for (int v = 0; v < 2; v++)
			    xv = xv + edge.GetVertexShape(v) * mesh.Point(edge.GetVertexNr(v));
			
			double secpoint = xi[l];
			
			ref->PointBetween (mesh.Point(edge.GetVertexNr(1)),
					   mesh.Point(edge.GetVertexNr(0)), secpoint,
					   mesh.GetFaceDescriptor(elem.GetIndex()).SurfNr(),
					   elem.GeomInfoPi(eledges[e][1]),
					   elem.GeomInfoPi(eledges[e][0]),
					   xexact, newgi);
		
			for (int k = 2; k <= edge.GetEdgeOrder(); k++)
			{
			    vec[0].Set(k-1, vec[0].Get(k-1) + Vec<3>(xexact - xv)(0)*edge.GetEdgeShape(k-2)*w );
			    vec[1].Set(k-1, vec[1].Get(k-1) + Vec<3>(xexact - xv)(1)*edge.GetEdgeShape(k-2)*w );
			    vec[2].Set(k-1, vec[2].Get(k-1) + Vec<3>(xexact - xv)(2)*edge.GetEdgeShape(k-2)*w );
			}
		
		    }


		    CalcInverse(mat,inv);
	    
		    Vector a0, a1, a2;
		    
		    a0 = inv*vec[0];
		    a1 = inv*vec[1];
		    a2 = inv*vec[2];

		    int index = edgecoeffsindex[edge.GetEdgeNr()-1];

		    for (int n = 0; n < npoints; n++, index++)
			edgecoeffs[index] =  Vec<3>(a0(n+1), a1(n+1), a2(n+1));
		}
	    }
	}


	for (int i=0; i<mesh.GetNSeg(); i++) 
	{
	    int edgenr = top.GetSegmentEdge(i+1);

            Segment s = mesh.LineSegment(i+1); 

	    segm.SetElementNumber (i+1);

	    int npoints = segm.GetEdgeOrder()-1;

	    if (npoints == 0) continue;

	    DenseMatrix mat(npoints);
	    DenseMatrix inv(npoints);
	    Vector vec[3];

	    for (int k = 0; k < 3; k++)
	    {
		vec[k].SetSize(npoints);
		for (int n = 1; n <= npoints; n++) vec[k].Set(n, 0.);
	    }
	
	    for (int l = 0; l < nIntegrationPoints; l++)
	    {
		double w = wi[l];

		segm.SetReferencePoint (Point<1>(xi[l]));
		segm.CalcVertexShapes ();
		segm.CalcEdgeShapes ();
		
		for (int n = 0; n < npoints; n++)
		    for (int m = 0; m < npoints; m++)
			mat.Set(n+1, m+1, mat.Get(n+1,m+1) +
				segm.GetEdgeShape(n) * segm.GetEdgeShape(m) * w);
		
		Point<3> xv(0,0,0);
		for (int v = 0; v < 2; v++)
		    xv = xv + segm.GetVertexShape(v) * mesh.Point(segm.GetVertexNr(v));
		
		double secpoint = xi[l];
		
		if (segm.GetEdgeOrientation() == -1) secpoint = 1. - secpoint; // reverse orientation
		
		ref->PointBetween (mesh.Point(segm.GetVertexNr(1)),
				   mesh.Point(segm.GetVertexNr(0)), secpoint,
				   s.surfnr2, s.surfnr1,
				   s.epgeominfo[1], s.epgeominfo[0],
				   xexact, newepgi);
		
		for (int k = 2; k <= segm.GetEdgeOrder(); k++)
		{
		    vec[0].Set(k-1, vec[0].Get(k-1) + Vec<3>(xexact - xv)(0)*segm.GetEdgeShape(k-2)*w );
		    vec[1].Set(k-1, vec[1].Get(k-1) + Vec<3>(xexact - xv)(1)*segm.GetEdgeShape(k-2)*w );
		    vec[2].Set(k-1, vec[2].Get(k-1) + Vec<3>(xexact - xv)(2)*segm.GetEdgeShape(k-2)*w );
		}
		
	    }


	    CalcInverse(mat,inv);
	    
	    Vector a0, a1, a2;

	    a0 = inv*vec[0];
	    a1 = inv*vec[1];
	    a2 = inv*vec[2];

	    int index = edgecoeffsindex[segm.GetEdgeNr()-1];

	    for (int n = 0; n < npoints; n++, index++)
		edgecoeffs[index] =  Vec<3>(a0(n+1), a1(n+1), a2(n+1));



	}

*/





	// evaluate face points

	if (mesh.GetDimension() == 3)
	  {
	    PopStatus ();
	    PushStatusF ("curving faces");
	    
	    for (int j=0; j<facecoeffsindex[top.GetNFaces()]; j++)
	      facecoeffs[j] = Vec<3>(0.,0.,0.);
	    
	    for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)   // for all surface elements
	      {
		if (multithread.terminate) return;

	        SetThreadPercent( double(100*i/mesh.GetNSE()) );

		Element2d elem = mesh[i];
		
		if (elem.GetType() == TRIG)
		    fe2d = &trig;
                else
		    fe2d = &quad;

		fe2d->SetElementNumber (i+1);

		int npoints = fe2d->GetNFaceShapes();

		if (npoints == 0) continue;

		DenseMatrix mat(npoints);
		DenseMatrix inv(npoints);
		Vector vec[3];

		for (int k = 0; k < 3; k++)
		{
		    vec[k].SetSize(npoints);
		    for (int n = 1; n <= npoints; n++) vec[k].Set(n, 0.);
		}

		for (int j = 0; j < nIntegrationPoints; j++)
		{
		    for (int k = 0; k < nIntegrationPoints; k++)
		    {
			double w;
			Point<2> xr;

			if (elem.GetType() == TRIG)
			  {
			    w = wi[j]*wi[k]*(1-xi[j]);
			    xr = Point<2> (xi[j], xi[k]*(1-xi[j]));
			  }
			else
			  {
			    w = wi[j]*wi[k];
			    xr = Point<2> (xi[j], xi[k]);
			  }

			fe2d->SetReferencePoint (xr);
			fe2d->CalcFaceShapes ();
			fe2d->CalcVertexShapes ();
			fe2d->CalcEdgeShapes ();
			fe2d->CalcFaceLaplaceShapes ();

			// integration over the product of the gradients of the face shapes

			for (int n = 0; n < npoints; n++)
			  for (int m = 0; m < npoints; m++)
			    mat.Set(n+1, m+1,
				    mat.Get(n+1,m+1) +
				    fe2d->GetFaceDShape(n)*fe2d->GetFaceDShape(m)*w);

			// integration over the difference between the exact geometry and the one
			// defined by vertex and edge shape functions times face shape

			// double giu = 0, giv = 0;
			PointGeomInfo gi;
			gi.trignum = elem.GeomInfoPi(1).trignum;
			gi.u = 0.0;
			gi.v = 0.0;
			Point<3> xve(0.,0.,0.);

			// vertex shape functions
			for (int v = 0; v < fe2d->GetNVertices(); v++)
			  {
			    xve = xve + fe2d->GetVertexShape(v) * mesh.Point(fe2d->GetVertexNr(v));
			    gi.u += fe2d->GetVertexShape(v) * elem.GeomInfoPi(v+1).u;
			    gi.v += fe2d->GetVertexShape(v) * elem.GeomInfoPi(v+1).v;
			  }

			// edge shape functions
			int index = 0;
			for (int e = 0; e < fe2d->GetNEdges(); e++)
			  {
			    int gindex = edgecoeffsindex[fe2d->GetEdgeNr(e)-1];
			    for (int k = 2; k <= fe2d->GetEdgeOrder(e); k++, index++, gindex++)
			      xve = xve + fe2d->GetEdgeShape(index) * edgecoeffs[gindex];
			  }

			// exact point

			Point<3> xexact = xve;
			ref->ProjectToSurface (xexact, mesh.GetFaceDescriptor(elem.GetIndex()).SurfNr(), gi);

			Vec<3> v2 = w*(Vec<3>(xexact)-Vec<3>(xve));

			for (int k = 0; k < 3; k++)
			  for (int n = 0; n < npoints; n++)
			    vec[k].Set(n+1, vec[k].Get(n+1) - fe2d->GetFaceLaplaceShape(n)*v2(k));
		    }
		}

		CalcInverse(mat,inv);
		
		Vector a0(npoints), a1(npoints), a2(npoints);
		
		/*
		a0 = inv*vec[0];
		a1 = inv*vec[1];
		a2 = inv*vec[2];
		*/
		inv.Mult (vec[0], a0);
		inv.Mult (vec[1], a1);
		inv.Mult (vec[2], a2);

		int index = facecoeffsindex[fe2d->GetFaceNr()-1];
		
		for (int n = 0; n < npoints; n++, index++)
		  facecoeffs[index] =  Vec<3>(a0.Elem(n+1), a1.Elem(n+1), a2.Elem(n+1));
	      }
	  }
	


	
/*
	cout << "WARNING: L2-Projection for faces" << endl;

	// evaluate face points

	if (mesh.GetDimension() == 3)
	{
	    for (int i=0; i<facecoeffsindex[top.GetNFaces()]; i++)
		facecoeffs[i] = Vec<3>(0.,0.,0.);

	    for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)   // for all surface elements
	    {
		Element2d elem = mesh[i];
		
		if (elem.GetType() == TRIG)
		    fe2d = &trig;
                else
		    fe2d = &quad;

		fe2d->SetElementNumber (i+1);

		int npoints = fe2d->GetNFaceShapes();

		if (npoints == 0) continue;

		DenseMatrix mat(npoints);
		DenseMatrix inv(npoints);
		Vector vec[3];

		for (int k = 0; k < 3; k++)
		{
		    vec[k].SetSize(npoints);
		    for (int n = 1; n <= npoints; n++) vec[k].Set(n, 0.);
		}

		for (int j = 0; j < nIntegrationPoints; j++)
		{
		    for (int k = 0; k < nIntegrationPoints; k++)
		    {
			double w;
			Point<2> xr;

			if (elem.GetType() == TRIG)
			{
			    w = wi[j]*wi[k]*(1-xi[j]);
			    xr = Point<2> (xi[j], xi[k]*(1-xi[j]));
			}
			else
			{
			    w = wi[j]*wi[k];
			    xr = Point<2> (xi[j], xi[k]);
			}

			fe2d->SetReferencePoint (xr);
//			fe2d->CalcFaceDShape (false, true);
			fe2d->CalcFaceShapes ();

			// integration over the product of the gradients of the face shapes

			for (int n = 0; n < npoints; n++)
			    for (int m = 0; m < npoints; m++)
				    mat.Set(n+1, m+1, mat.Get(n+1,m+1) +
					    fe2d->GetFaceShape(n)*fe2d->GetFaceShape(m)*w);

			// integration over the difference between the exact geometry and the one
			// defined by vertex and edge shape functions times face shape

			Point<3> xve(0.,0.,0.);

			// vertex shape functions
			fe2d->CalcVertexShapes ();
//			fe2d->CalcVertexShape (true, false);
			for (int v = 0; v < fe2d->GetNVertices(); v++)
			    xve = xve + fe2d->GetVertexShape(v) * mesh.Point(fe2d->GetVertexNr(v));

			// edge shape functions
//			fe2d->CalcEdgeShape (true, false);
			fe2d->CalcEdgeShapes ();

			int index = 0;
			for (int e = 0; e < fe2d->GetNEdges(); e++)
			{
			    int gindex = edgecoeffsindex[fe2d->GetEdgeNr(e)-1];

			    for (int k = 2; k <= fe2d->GetEdgeOrder(e); k++, index++, gindex++)
				xve = xve + fe2d->GetEdgeShape(index) * edgecoeffs[gindex];
			}

			// exact point

			Point<3> xexact = xve;
			ref->ProjectToSurface (xexact, mesh.GetFaceDescriptor(elem.GetIndex()).SurfNr());

			Vec<3> v = w*(Vec<3>(xexact)-Vec<3>(xve));

			fe2d->CalcFaceLaplaceShapes ();

			for (int k = 0; k < 3; k++)
			    for (int n = 0; n < npoints; n++)
				vec[k].Set(n+1, vec[k].Get(n+1) + fe2d->GetFaceShape(n)*v(k));
			}
		    }

 		    CalcInverse(mat,inv);

		    Vector a0, a1, a2;

		    a0 = inv*vec[0];
		    a1 = inv*vec[1];
		    a2 = inv*vec[2];

		    int index = facecoeffsindex[fe2d->GetFaceNr()-1];

		    for (int n = 0; n < npoints; n++, index++)
			facecoeffs[index] =  Vec<3>(a0(n+1), a1(n+1), a2(n+1));
	    }
	}
*/


    PrintMessage (5, "reducing order");
   
    for (e = 0; e < top.GetNEdges(); e++)
      if (edgeorder[e] > 1)
	{
	  int i;
	  double maxcoeff = 0.;

	  for (i = edgecoeffsindex[e]; i < edgecoeffsindex[e+1]; i++)
	    maxcoeff = max2 (maxcoeff, edgecoeffs[i].Length());

	  if (maxcoeff < 1e-12) edgeorder[e] = 1;
	}

    for (f = 0; f < top.GetNFaces(); f++)
      if (faceorder[f] > 1)
	{
	  int i;
	  double maxcoeff = 0.;

	  for (i = facecoeffsindex[f]; i < facecoeffsindex[f+1]; i++)
	    maxcoeff = max (maxcoeff, facecoeffs[i].Length());

	  if (maxcoeff < 1e-12) faceorder[f] = 1;
	}
    
    isHighOrder = 1;              // true

    PrintMessage(1, "done");
    PopStatus();
    //	cout << "finished" << endl;
    }

} // namespace netgen

#endif
