#include <mystdlib.h>

#include "meshing.hpp"
#ifndef CURVEDELEMS_NEW

namespace netgen
{
    
    
  // computes Gaussean integration formula on (0,1) with n points
  // in: Numerical algs in C (or, was it the Fortran book ?)
  void ComputeGaussRule (int n, ARRAY<double> & xi, ARRAY<double> & wi)
  {
    xi.SetSize (n);
    wi.SetSize (n);
	
    int m = (n+1)/2;
    double p1, p2, p3;
    double pp, z, z1;
    for (int i = 1; i <= m; i++)
      {
	z = cos ( M_PI * (i - 0.25) / (n + 0.5));
	while(1)
	  {
	    p1 = 1; p2 = 0;
	    for (int j = 1; j <= n; j++)
	      {
		p3 = p2; p2 = p1;
		p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
	      }
	    // p1 is legendre polynomial

	    pp = n * (z*p1-p2) / (z*z - 1);
	    z1 = z;
	    z = z1-p1/pp;

	    if (fabs (z - z1) < 1e-14) break;
	  }

	xi[i-1] = 0.5 * (1 - z);
	xi[n-i] = 0.5 * (1 + z);
	wi[i-1] = wi[n-i] = 1.0 / ( (1  - z * z) * pp * pp);
      }
  }



  // ----------------------------------------------------------------------------
  //      PolynomialBasis
  // ----------------------------------------------------------------------------


  void PolynomialBasis :: CalcLegendrePolynomials (double x)
  {
    double p1 = 1.0, p2 = 0.0, p3;

    lp[0] = 1.0;

    for (int j=1; j<=order; j++)
      {
	p3=p2; p2=p1;
	p1=((2.0*j-1.0)*(2*x-1)*p2-(j-1.0)*p3)/j;
	lp[j] = p1;
      }
  }


  void PolynomialBasis :: CalcScaledLegendrePolynomials (double x, double t)
  {
    double p1 = 1.0, p2 = 0.0, p3;

    lp[0] = 1.0;
	
    double x2mt = 2*x-t;
    double t2 = t*t;

    for (int j=1; j<=order; j++)
      {
	p3=p2; p2=p1;
	p1=((2.0*j-1.0)*x2mt*p2-t2*(j-1.0)*p3)/j;
	lp[j] = p1;
      }
  }


  void PolynomialBasis :: CalcDLegendrePolynomials (double x)
  {
    double p1 = 0., p2 = 0., p3;

    dlp[0] = 0.;

    for (int j = 1; j <= order-1; j++)
      {
	p3=p2; p2=p1;
	p1=((2.*j-1.)*(2*lp[j-1]+(2*x-1)*p2)-(j-1.)*p3)/j;
	dlp[j] = p1;
      }
  }


  void PolynomialBasis :: CalcF (double x)
  {
    CalcLegendrePolynomials (x);

    for (int j = 0; j<=order-2; j++)
      f[j] = (lp[j+2]-lp[j])/(2.0*(j+1)+1)/2.0;
  }


  void PolynomialBasis :: CalcDf (double x)
  {
    CalcLegendrePolynomials (x);

    for (int j = 0; j <= order-2; j++)
      df[j] = lp[j+1];
  }


  void PolynomialBasis :: CalcFDf (double x)
  {
    CalcLegendrePolynomials (x);

    for (int j = 0; j<=order-2; j++)
      {
	f[j] = (lp[j+2]-lp[j])/(2.0*(j+1)+1)/2.0;
	df[j] = lp[j+1];
      }
  }


  void PolynomialBasis :: CalcDDf (double x)
  {
    CalcLegendrePolynomials (x);
    CalcDLegendrePolynomials (x);

    for (int j = 0; j <= order-2; j++) ddf[j] = dlp[j+1];
  }


  void PolynomialBasis :: CalcFScaled (double x, double t)
  {
    double tt = t*t;
    double x2mt = 2*x-t;


      

    double p1 = 0.5*x2mt, p2 = -0.5, p3 = 0.0;
    for (int j=2; j<=order; j++)
      {
	p3=p2; p2=p1;
	p1=((2.0*j-3.0) * x2mt * p2 - tt*(j-3.0)*p3)/j;
	f[j-2] = p1;
      }

    /*
      CalcF (x/t);
      double fac = t*t;
      for (int j = 0; j<=order-2; j++, fac*=t) 
      f[j] *= fac;
    */
  }

  void PolynomialBasis :: CalcFScaled (int p, double x, double t, double * values)
  {
    double tt = t*t;
    double x2mt = 2*x-t;
    
    double p1 = 0.5*x2mt, p2 = -0.5, p3 = 0.0;
    for (int j=2; j<=p; j++)
      {
	p3=p2; p2=p1;
	p1=((2.0*j-3.0) * x2mt * p2 - tt*(j-3.0)*p3)/j;
	values[j-2] = p1;
      }
    
    /*
      CalcF (x/t);
      double fac = t*t;
      for (int j = 0; j<=order-2; j++, fac*=t) 
      f[j] *= fac;
    */
  }




  // ----------------------------------------------------------------------------
  //      BaseFiniteElement1D
  // ----------------------------------------------------------------------------


  void BaseFiniteElement1D :: CalcVertexShapes ()
  {
    vshape[0] = xi(0);
    vshape[1] = 1-xi(0);

    vdshape[0] = 1;
    vdshape[1] = -1;

    /*
      if (edgeorient == -1)
      {
      Swap (vshape[0], vshape[1]);
      Swap (vdshape[0], vdshape[1]);
      }
    */
	
  }


  void BaseFiniteElement1D :: CalcEdgeShapes ()
  {
    b.SetOrder (edgeorder);
    if (edgeorient == 1)
      b.CalcFDf( 1-xi(0) );
    else
      b.CalcFDf( xi(0) );

    for (int k = 2; k <= edgeorder; k++)
      {
	eshape[k-2] = b.GetF(k);
	edshape[k-2] = -b.GetDf(k);
      }
  }


  void BaseFiniteElement1D :: CalcEdgeLaplaceShapes ()
  {
    b.SetOrder (edgeorder);
    if (edgeorient == 1)
      b.CalcDDf( 1-xi(0) );
    else
      b.CalcDDf( xi(0) );

    for (int k = 2; k <= edgeorder; k++)
      eddshape[k-2] = b.GetDDf(k);
  }




  // ----------------------------------------------------------------------------
  //      BaseFiniteElement2D
  // ----------------------------------------------------------------------------


  void BaseFiniteElement2D :: SetElementNumber (int aelnr)
  {
    int locmaxedgeorder = -1;
	
    BaseFiniteElement<2> :: SetElementNumber (aelnr);
    const Element2d & elem = mesh[(SurfaceElementIndex) (elnr-1)]; 
    top.GetSurfaceElementEdges (elnr, &(edgenr[0]), &(edgeorient[0]));
    facenr = top.GetSurfaceElementFace (elnr);
    faceorient = top.GetSurfaceElementFaceOrientation (elnr);
	
    for (int v = 0; v < nvertices; v++)
      vertexnr[v] = elem[v];
	
    for (int e = 0; e < nedges; e++)
      {
	edgeorder[e] = curv.GetEdgeOrder (edgenr[e]-1); // 1-based
	locmaxedgeorder = max2 (edgeorder[e], locmaxedgeorder);
      }
	
    faceorder = curv.GetFaceOrder (facenr-1); // 1-based
    CalcNFaceShapes ();
	
    if (locmaxedgeorder > maxedgeorder)
      {
	maxedgeorder = locmaxedgeorder;
	eshape.SetSize(nedges * (maxedgeorder-1));
	edshape.SetSize(nedges * (maxedgeorder-1));
      }
	
    if (faceorder > maxfaceorder)
      {
	maxfaceorder = faceorder;
	fshape.SetSize( nfaceshapes ); // number of face shape functions
	fdshape.SetSize( nfaceshapes );
	fddshape.SetSize ( nfaceshapes );
      }
  };
    

  // ----------------------------------------------------------------------------
  //      BaseFiniteElement3D
  // ----------------------------------------------------------------------------


  void BaseFiniteElement3D :: SetElementNumber (int aelnr)
  {
    int locmaxedgeorder = -1;
    int locmaxfaceorder = -1;
    int v, f, e;

    BaseFiniteElement<3> :: SetElementNumber (aelnr);
    Element elem = mesh[(ElementIndex) (elnr-1)];
    top.GetElementEdges (elnr, &(edgenr[0]), &(edgeorient[0]));
    top.GetElementFaces (elnr, &(facenr[0]), &(faceorient[0]));
	
    for (v = 0; v < nvertices; v++)
      vertexnr[v] = elem[v];
	
    for (f = 0; f < nfaces; f++)
      {
	surfacenr[f] = top.GetFace2SurfaceElement (facenr[f]);
	// surfaceorient[f] = top.GetSurfaceElementFaceOrientation (surfacenr[f]);
      }
	
    for (e = 0; e < nedges; e++)
      {
	edgeorder[e] = curv.GetEdgeOrder (edgenr[e]-1); // 1-based
	locmaxedgeorder = max (edgeorder[e], locmaxedgeorder);
      }
	
    for (f = 0; f < nfaces; f++)
      {
	faceorder[f] = curv.GetFaceOrder (facenr[f]-1); // 1-based
	locmaxfaceorder = max (faceorder[f], locmaxfaceorder);
      }
	
    CalcNFaceShapes ();
	
    if (locmaxedgeorder > maxedgeorder)
      {
	maxedgeorder = locmaxedgeorder;
	eshape.SetSize(nedges * (maxedgeorder-1));
	edshape.SetSize(nedges * (maxedgeorder-1));
      }
	
    if (locmaxfaceorder > maxfaceorder)
      {
	maxfaceorder = locmaxfaceorder;
	fshape.SetSize( nfaces * (maxfaceorder-1) * (maxfaceorder-1)); // number of face shape functions
	fdshape.SetSize( nfaces * (maxfaceorder-1) * (maxfaceorder-1));
      }
  };
    
    


  // ----------------------------------------------------------------------------
  //      FETrig
  // ----------------------------------------------------------------------------


  void FETrig :: SetReferencePoint (Point<2> axi)
  {
    BaseFiniteElement2D :: SetReferencePoint (axi);
    lambda(0) = xi(0);
    lambda(1) = xi(1);
    lambda(2) = 1-xi(0)-xi(1);

    dlambda(0,0) =  1; dlambda(0,1) =  0;
    dlambda(1,0) =  0; dlambda(1,1) =  1;
    dlambda(2,0) = -1; dlambda(2,1) = -1;
  }


  void FETrig :: SetVertexSingularity (int v, int exponent)
  {
    int i;
    if (1-lambda(v) < EPSILON) return;

    Point<3> lamold = lambda;

    Vec<2> dfac;

    double fac = pow(1-lambda(v),exponent-1);

    for (i = 0; i < 2; i++)
      {
	dfac(i) = -(exponent-1)*pow(1-lambda(v),exponent-2)*dlambda(v,i);
	dlambda(v,i) *= exponent * pow(1-lambda(v),exponent-1);
      }

    lambda(v) = 1-pow(1-lambda(v),exponent);

    for (i = 0; i < nvertices; i++)
      {
	if (i == v) continue;
	for (int j = 0; j < 2; j++)
	  dlambda(i,j) = dlambda(i,j) * fac + lamold(i) * dfac(j);

	lambda(i) *= fac;
      }
  }



  void FETrig :: CalcVertexShapes ()
  {
    for (int v = 0; v < nvertices; v++)
      {
	vshape[v] = lambda(v);
	vdshape[v](0) = dlambda(v,0);
	vdshape[v](1) = dlambda(v,1);
      }
  }


  void FETrig :: CalcEdgeShapes ()
  {
    int index = 0;
    for (int e = 0; e < nedges; e++)
      {
	if (edgeorder[e] <= 1) continue;

	int i0 = eledge[e][0]-1;
	int i1 = eledge[e][1]-1;

	double arg = lambda(i0) + lambda(i1); // = 1-lambda[i2];

	if (fabs(arg) < EPSILON) // division by 0
	  {
	    for (int k = 2; k <= edgeorder[e]; k++)
	      {
		eshape[index] = 0;
		edshape[index] = Vec<2>(0,0);
		index++;
	      }
	    continue;
	  }

	if (edgeorient[e] == -1) Swap (i0, i1); // reverse orientation

	double xi = lambda(i1)/arg;

	b1.SetOrder (edgeorder[e]);
	b1.CalcFDf (xi);

	double decay = arg;
	double ddecay;
	    
	double l0 = lambda(i0);
	double l0x = dlambda(i0,0);
	double l0y = dlambda(i0,1);

	double l1 = lambda(i1);
	double l1x = dlambda(i1,0);
	double l1y = dlambda(i1,1);

	for (int k = 2; k <= edgeorder[e]; k++)
	  {        
	    ddecay = k*decay;
	    decay *= arg;
		
	    eshape[index] = b1.GetF (k) * decay;
	    edshape[index] = Vec<2>
	      (b1.GetDf(k) * (l1x*arg - l1*(l0x+l1x)) * 
	       decay / (arg * arg) + b1.GetF(k) * ddecay * (l0x+l1x),
	       b1.GetDf(k) * (l1y*arg - l1*(l0y+l1y)) *
	       decay / (arg * arg) + b1.GetF(k) * ddecay * (l0y+l1y));
	    index++;
	  }
      }
    // (*testout) << "index = " << index << endl;
    // (*testout) << "eshape = " << eshape << ", edshape = " << edshape << endl;


    // eshape = 0.0;
    //	edshape = 0.0;

    /*
      int index = 0;
      for (int e = 0; e < nedges; e++)
      {
      if (edgeorder[e] <= 1) continue;

      int i0 = eledge[e][0]-1;
      int i1 = eledge[e][1]-1;

      if (edgeorient[e] == -1) Swap (i0, i1); // reverse orientation
	    
      double x = lambda(i1)-lambda(i0);
      double y = 1-lambda(i0)-lambda(i1);
      double fy = (1-y)*(1-y);

      // double p3 = 0, p3x = 0, p3y = 0;
      // double p2 = -1, p2x = 0, p2y = 0;
      // double p1 = x, p1x = 1, p1y = 0;

      double p3 = 0, p3x = 0, p3y = 0;
      double p2 = -0.5, p2x = 0, p2y = 0;
      double p1 = 0.5*x, p1x = 0.5, p1y = 0;

      for (int j=2; j<= edgeorder[e]; j++)
      {
      p3=p2; p3x = p2x; p3y = p2y;
      p2=p1; p2x = p1x; p2y = p1y;
      double c1 = (2.0*j-3) / j;
      double c2 = (j-3.0) / j;
		
      p1  = c1 * x * p2 - c2 * fy * p3;
      p1x = c1 * p2 + c1 * x * p2x - c2 * fy * p3x;
      p1y = c1 * x * p2y - (c2 * 2 * (y-1) * p3 + c2 * fy * p3y);

      eshape[index] = p1;
      for (int k = 0; k < 2; k++)
      edshape[index](k) = 
      p1x * ( dlambda(i1,k) - dlambda(i0,k) ) + 
      p1y * (-dlambda(i1,k) - dlambda(i0,k) );
      index++;
      }    

      }
      // (*testout) << "eshape = " << eshape << ", edshape = " << edshape << endl;
      */
  }


  void FETrig :: CalcFaceShapes ()
  {
    int index = 0;

    int i0 = elface[0][0]-1;
    int i1 = elface[0][1]-1;
    int i2 = elface[0][2]-1;

    // sort lambda_i's by the corresponing vertex numbers

    if (vertexnr[i1] < vertexnr[i0]) Swap (i0, i1);
    if (vertexnr[i2] < vertexnr[i0]) Swap (i0, i2);
    if (vertexnr[i2] < vertexnr[i1]) Swap (i1, i2);

    double arg = lambda(i1) + lambda(i2);

    if (fabs(arg) < EPSILON) // division by 0
      {
	for (int k = 0; k < nfaceshapes; k++)
	  {
	    fshape[index] = 0;
	    fdshape[index] = Vec<2>(0,0);
	    index++;
	  }
	return;
      }

    b1.SetOrder (faceorder);
    b2.SetOrder (faceorder);

    b1.CalcFDf (lambda(i0));
    b2.CalcFDf (lambda(i2)/arg);

    double decay = 1;
    double ddecay;

    double l0 = lambda(i0);
    double l1 = lambda(i1);
    double l2 = lambda(i2);
    double dl0x = dlambda(i0,0);
    double dl0y = dlambda(i0,1);
    double dl1x = dlambda(i1,0);
    double dl1y = dlambda(i1,1);
    double dl2x = dlambda(i2,0);
    double dl2y = dlambda(i2,1);

    double dargx = dl1x + dl2x;
    double dargy = dl1y + dl2y;

    for (int n1 = 2; n1 < faceorder; n1++)
      {
	ddecay = (n1-1)*decay;
	decay *= arg;
	    
	for (int n0 = 2; n0 < faceorder-n1+2; n0++)
	  {
	    fshape[index] = b1.GetF(n0) * b2.GetF(n1) * decay;
	    fdshape[index] = Vec<2>
	      (b1.GetDf(n0) * dl0x * b2.GetF(n1) * decay +
	       b1.GetF(n0) * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay +
	       b1.GetF(n0) * b2.GetF(n1) * ddecay * dargx,
	       b1.GetDf(n0) * dl0y * b2.GetF(n1) * decay +
	       b1.GetF(n0) * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay +
	       b1.GetF(n0) * b2.GetF(n1) * ddecay * dargy);
	    index++;
	  }
      }
  }



  void FETrig :: CalcFaceLaplaceShapes ()
  {
    int index = 0;

    int i0 = elface[0][0]-1;
    int i1 = elface[0][1]-1;
    int i2 = elface[0][2]-1;

    if (vertexnr[i1] < vertexnr[i0]) Swap (i0, i1);
    if (vertexnr[i2] < vertexnr[i0]) Swap (i0, i2);
    if (vertexnr[i2] < vertexnr[i1]) Swap (i1, i2);

    double arg = lambda(i1) + lambda(i2);

    if (fabs(arg) < EPSILON) // division by 0
      {
	for (int k = 0; k < nfaceshapes; k++)
	  fddshape[k] = 0;
	return;
      }

    b1.SetOrder (faceorder);
    b2.SetOrder (faceorder);

    b1.CalcFDf (lambda(i0));
    b1.CalcDDf (lambda(i0));
    b2.CalcFDf (lambda(i2)/arg);
    b2.CalcDDf (lambda(i2)/arg);

    double decay = 1;
    double ddecay = 0;
    double dddecay;

    double l0 = lambda(i0);
    double l1 = lambda(i1);
    double l2 = lambda(i2);
    double dl0x = dlambda(i0,0);
    double dl0y = dlambda(i0,1);
    double dl1x = dlambda(i1,0);
    double dl1y = dlambda(i1,1);
    double dl2x = dlambda(i2,0);
    double dl2y = dlambda(i2,1);

    double dargx = dl1x + dl2x;
    double dargy = dl1y + dl2y;

    for (int n1 = 2; n1 < faceorder; n1++)
      {
	dddecay = (n1-1)*ddecay;
	ddecay = (n1-1)*decay;
	decay *= arg;
	    
	for (int n0 = 2; n0 < faceorder-n1+2; n0++)
	  {
	    fddshape[index] = 

	      //  b1.GetDf(n0) * dl0x * b2.GetF(n1) * decay .... derived

	      b1.GetDDf(n0) * dl0x * dl0x * b2.GetF(n1) * decay +
	      b1.GetDf(n0) * dl0x * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay +
	      b1.GetDf(n0) * dl0x * b2.GetF(n1) * ddecay * dargx +

		    
	      //  b1.GetF(n0) * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay ... derived

	      b1.GetDf(n0) * dl0x * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay +
	      b1.GetF(n0) * b2.GetDDf(n1) * pow((dl2x * arg - l2 * dargx)/(arg*arg),2) * decay +
	      b1.GetF(n0) * b2.GetDf(n1) * (-2*dargx/arg) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay +
	      b1.GetF(n0) * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * ddecay * dargx +

		    
	      //  b1.GetF(n0) * b2.GetF(n1) * ddecay * dargx ... derived

	      b1.GetDf(n0) * dl0x * b2.GetF(n1) * ddecay * dargx +
	      b1.GetF(n0) * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * ddecay * dargx +
	      b1.GetF(n0) * b2.GetF(n1) * dddecay * dargx * dargx +

		    
	      //  b1.GetDf(n0) * dl0y * b2.GetF(n1) * decay ... derived

	      b1.GetDDf(n0) * dl0y * dl0y * b2.GetF(n1) * decay +
	      b1.GetDf(n0) * dl0y * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay +
	      b1.GetDf(n0) * dl0y * b2.GetF(n1) * ddecay * dargy +
		    

	      //  b1.GetF(n0) * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay ... derived

	      b1.GetDf(n0) * dl0y * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay +
	      b1.GetF(n0) * b2.GetDDf(n1) * pow((dl2y * arg - l2 * dargy)/(arg*arg),2) * decay +
	      b1.GetF(n0) * b2.GetDf(n1) * (-2*dargy/arg) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay +
	      b1.GetF(n0) * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * ddecay * dargy +
		    

	      //  b1.GetF(n0) * b2.GetF(n1) * ddecay * dargy ... derived

	      b1.GetDf(n0) * dl0y * b2.GetF(n1) * ddecay * dargy +
	      b1.GetF(n0) * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * ddecay * dargy +
	      b1.GetF(n0) * b2.GetF(n1) * dddecay * dargy * dargy;

	    index++;
	  }
      }
  }



  // ----------------------------------------------------------------------------
  //      FEQuad
  // ----------------------------------------------------------------------------


  void FEQuad :: CalcVertexShapes ()
  {
    vshape[0] = (1-xi(0)) * (1-xi(1));
    vshape[1] = (  xi(0)) * (1-xi(1));
    vshape[2] = (  xi(0)) * (  xi(1));
    vshape[3] = (1-xi(0)) * (  xi(1));

    vdshape[0] = Vec<2> ( -(1-xi(1)), -(1-xi(0)) );
    vdshape[1] = Vec<2> (  (1-xi(1)), -(  xi(0)) );
    vdshape[2] = Vec<2> (  (  xi(1)),  (  xi(0)) );
    vdshape[3] = Vec<2> ( -(  xi(1)),  (1-xi(0)) );
  }


  void FEQuad :: CalcEdgeShapes ()
  {
    int index = 0;

    double arg0[4] = { xi(0), 1-xi(0), 1-xi(1), xi(1) };
    double arg1[4] = { 1-xi(1), xi(1), 1-xi(0), xi(0) };
    double darg0[4] = {  1, -1, -1,  1 };
    double darg1[4] = { -1,  1, -1,  1 };
	
    for (int e = 0; e < nedges; e++)
      {
	b1.SetOrder (edgeorder[e]);
	b1.CalcFDf (edgeorient[e] == 1 ? arg0[e] : 1-arg0[e]);

	double decay = arg1[e];
	double ddecay;

	for (int k = 2; k <= edgeorder[e]; k++, index++)
	  {
	    ddecay = k*decay;
	    decay *= arg1[e];
	    // linear decay
	    eshape[index] = b1.GetF(k) * arg1[e];

	    if (e < 2)
	      edshape[index] = Vec<2>
		(darg0[e] * edgeorient[e] * b1.GetDf(k) * arg1[e],
		 b1.GetF(k) * darg1[e]);
	    else
	      edshape[index] = Vec<2>
		(b1.GetF(k) * darg1[e],
		 darg0[e] * edgeorient[e] * b1.GetDf(k) * arg1[e]);

	    // exponential decay
	    /*		eshape[index] = b1.GetF(k) * decay;

	    if (e < 2)
	    edshape[index] = Vec<2>
	    (darg0[e] * edgeorient[e] * b1.GetDf(k) * decay,
	    b1.GetF(k) * ddecay * darg1[e]);
	    else
	    edshape[index] = Vec<2>
	    (b1.GetF(k) * ddecay * darg1[e],
	    darg0[e] * edgeorient[e] * b1.GetDf(k) * decay);
	    */
	  }
      }
  }


  void FEQuad :: CalcFaceShapes ()
  {
    int index = 0;

    // find index of point with smallest number

    int i0 = 0;
    for (int i = 1; i < 4; i++)
      if (vertexnr[elface[0][i]-1] < vertexnr[elface[0][i0]-1]) i0 = i;

    double x;
    double y;
    double dxx;
    double dxy;
    double dyx;
    double dyy;

    switch (i0)
      {
      case 0:
	x = xi(0); y = xi(1);
	dxx = 1; dxy = 0;
	dyx = 0; dyy = 1;
	break;
      case 1:
	x = xi(1); y = 1-xi(0);
	dxx = 0; dxy = 1;
	dyx = -1; dyy = 0;
	break;
      case 2: 
	x = 1-xi(0); y = 1-xi(1);
	dxx = -1; dxy = 0;
	dyx = 0; dyy = -1;
	break;
      case 3:
	x = 1-xi(1); y = xi(0);
	dxx = 0; dxy =-1;
	dyx = 1; dyy = 0;
	break;
      }

    if (vertexnr[elface[0][(i0+1) % 4]-1] > vertexnr[elface[0][(i0+3) % 4]-1]) 
      {
	Swap (x,y);
	Swap (dxx, dyx);
	Swap (dxy, dyy);
      }

    b1.SetOrder (faceorder);
    b2.SetOrder (faceorder);

    b1.CalcFDf (x);
    b2.CalcFDf (y);

    for (int n0 = 2; n0 <= faceorder; n0++)
      for (int n1 = 2; n1 <= faceorder; n1++)
	{
	  fshape[index] = b1.GetF(n0) * b2.GetF(n1);
	  fdshape[index] = Vec<2> (b1.GetDf(n0) * dxx * b2.GetF(n1) + b1.GetF(n0) * b2.GetDf(n1) * dyx,
				   b1.GetDf(n0) * dxy * b2.GetF(n1) + b1.GetF(n0) * b2.GetDf(n1) * dyy);
	  index++;
	}
  }


  void FEQuad :: CalcFaceLaplaceShapes ()
  {
    int index = 0;

    // find index of point with smallest number

    int i0 = 0;
    for (int i = 1; i < 4; i++)
      if (vertexnr[elface[0][i]-1] < vertexnr[elface[0][i0]-1]) i0 = i;

    double x;
    double y;
    double dxx;
    double dxy;
    double dyx;
    double dyy;

    switch (i0)
      {
      case 0:
	x = xi(0); y = xi(1);
	dxx = 1; dxy = 0;
	dyx = 0; dyy = 1;
	break;
      case 1:
	x = xi(1); y = 1-xi(0);
	dxx = 0; dxy = 1;
	dyx = -1; dyy = 0;
	break;
      case 2: 
	x = 1-xi(0); y = 1-xi(1);
	dxx = -1; dxy = 0;
	dyx = 0; dyy = -1;
	break;
      case 3:
	x = 1-xi(1); y = xi(0);
	dxx = 0; dxy =-1;
	dyx = 1; dyy = 0;
	break;
      }

    if (vertexnr[elface[0][(i0+1) % 4]-1] > vertexnr[elface[0][(i0+3) % 4]-1]) 
      {
	Swap (x,y);
	Swap (dxx, dyx);
	Swap (dxy, dyy);
      }

    b1.SetOrder (faceorder);
    b2.SetOrder (faceorder);

    b1.CalcFDf (x);
    b1.CalcDDf (x);
    b2.CalcFDf (y);
    b2.CalcDDf (y);

    for (int n0 = 2; n0 <= faceorder; n0++)
      for (int n1 = 2; n1 <= faceorder; n1++)
	{
	  fddshape[index] =
	    b1.GetDDf(n0) * dxx * dxx * b2.GetF(n1) +
	    2*  b1.GetDf(n0) * dxx * b2.GetDf(n1) * dyx +
	    b1.GetF(n0) * b2.GetDDf(n1) * dyx * dyx +

	    b1.GetDDf(n0) * dxy * dxy * b2.GetF(n1) +
	    2*  b1.GetDf(n0) * dxy * b2.GetDf(n1) * dyy +
	    b1.GetF(n0) * b2.GetDDf(n1) * dyy * dyy;

	  index++;
	}
  }



  // ----------------------------------------------------------------------------
  //      FETet
  // ----------------------------------------------------------------------------


  void FETet :: SetReferencePoint (Point<3> axi)
  {
    BaseFiniteElement3D :: SetReferencePoint (axi);
	
    lambda(0) = xi(0);
    lambda(1) = xi(1);
    lambda(2) = xi(2);
    lambda(3) = 1-xi(0)-xi(1)-xi(2);

    dlambda(0,0) =  1; dlambda(0,1) =  0; dlambda(0,2) =   0;
    dlambda(1,0) =  0; dlambda(1,1) =  1; dlambda(1,2) =   0;
    dlambda(2,0) =  0; dlambda(2,1) =  0; dlambda(2,2) =   1;
    dlambda(3,0) = -1; dlambda(3,1) = -1; dlambda(3,2) =  -1;
  }


  void FETet :: CalcVertexShapes ()
  {
    for (int v = 0; v < nvertices; v++)
      {
	vshape[v] = lambda(v);
	vdshape[v](0) = dlambda(v,0);
	vdshape[v](1) = dlambda(v,1);
	vdshape[v](2) = dlambda(v,2);
      }
  }

  void FETet :: CalcVertexShapesOnly ()
  {
    for (int v = 0; v < nvertices; v++)
      vshape[v] = lambda(v);
  }


  void FETet :: CalcEdgeShapes ()
  {
    int index = 0;

    for (int e = 0; e < nedges; e++)
      {
	int i0 = eledge[e][0]-1;
	int i1 = eledge[e][1]-1;

	double arg = lambda(i0)+lambda(i1);

	if (fabs(arg) < EPSILON) // division by 0
	  {
	    for (int k = 2; k <= edgeorder[e]; k++)
	      {
		eshape[index] = 0;
		edshape[index] = Vec<3>(0,0,0);
		index++;
	      }
	    continue;
	  }

	if (edgeorient[e] == -1) Swap (i0, i1);

	double xi = lambda[i1]/arg;

	b1.SetOrder (edgeorder[e]);
	b1.CalcFDf (xi);

	double decay = arg;
	double ddecay;
	    
	double l0 = lambda(i0);
	double dl0x = dlambda(i0,0);
	double dl0y = dlambda(i0,1);
	double dl0z = dlambda(i0,2);

	double l1 = lambda(i1);
	double dl1x = dlambda(i1,0);
	double dl1y = dlambda(i1,1);
	double dl1z = dlambda(i1,2);

	double dargx = dl0x + dl1x;
	double dargy = dl0y + dl1y;
	double dargz = dl0z + dl1z;
                           
	for (int k = 2; k <= edgeorder[e]; k++)
	  {        
	    ddecay = k*decay;
	    decay *= arg;

	    eshape[index] = b1.GetF (k) * decay;
	    edshape[index] = Vec<3>
	      (b1.GetDf(k) * (dl1x*arg - l1*dargx) * 
	       decay / (arg * arg) + b1.GetF(k) * ddecay * dargx,
	       b1.GetDf(k) * (dl1y*arg - l1*dargy) * 
	       decay / (arg * arg) + b1.GetF(k) * ddecay * dargy,
	       b1.GetDf(k) * (dl1z*arg - l1*dargz) * 
	       decay / (arg * arg) + b1.GetF(k) * ddecay * dargz);

	    index++;
	  }
      }



    /*
      int index = 0;
      for (int e = 0; e < nedges; e++)
      {
      if (edgeorder[e] <= 1) continue;

      int i0 = eledge[e][0]-1;
      int i1 = eledge[e][1]-1;

      if (edgeorient[e] == -1) Swap (i0, i1); // reverse orientation
	    
      double x = lambda(i1)-lambda(i0);
      double y = 1-lambda(i0)-lambda(i1);
      double fy = (1-y)*(1-y);

      // double p3 = 0, p3x = 0, p3y = 0;
      // double p2 = -1, p2x = 0, p2y = 0;
      // double p1 = x, p1x = 1, p1y = 0;

      double p3 = 0, p3x = 0, p3y = 0;
      double p2 = -0.5, p2x = 0, p2y = 0;
      double p1 = 0.5*x, p1x = 0.5, p1y = 0;

      for (int j=2; j<= edgeorder[e]; j++)
      {
      p3=p2; p3x = p2x; p3y = p2y;
      p2=p1; p2x = p1x; p2y = p1y;
      double c1 = (2.0*j-3) / j;
      double c2 = (j-3.0) / j;
		
      p1  = c1 * x * p2 - c2 * fy * p3;
      p1x = c1 * p2 + c1 * x * p2x - c2 * fy * p3x;
      p1y = c1 * x * p2y - (c2 * 2 * (y-1) * p3 + c2 * fy * p3y);

      eshape[index] = p1;
      for (int k = 0; k < 3; k++)
      edshape[index](k) = 
      p1x * ( dlambda(i1,k) - dlambda(i0,k) ) + 
      p1y * (-dlambda(i1,k) - dlambda(i0,k) );
      index++;
      }    

      }

    */




  }




  void FETet :: CalcEdgeShapesOnly ()
  {
    int index = 0;

    for (int e = 0; e < nedges; e++)
      {
	int i0 = eledge[e][0]-1;
	int i1 = eledge[e][1]-1;

	if (edgeorient[e] == -1) Swap (i0, i1);

	double arg = lambda(i0)+lambda(i1);
	double xi = lambda[i1];

	/*
	  b1.SetOrder (edgeorder[e]);
	  b1.CalcFScaled (xi, arg);
	    
	  for (int k = 2; k <= edgeorder[e]; k++, index++)
	  eshape[index] = b1.GetF (k); 
	*/
	b1.CalcFScaled (edgeorder[e], xi, arg, &eshape[index]);
	index += edgeorder[e]-1;
      }
  }






  void FETet :: CalcFaceShapes ()
  {
    int index = 0;

    for (int f = 0; f < nfaces; f++)
      {
	if (faceorder[f] <= 2) continue;

	int i0 = elface[f][0]-1;
	int i1 = elface[f][1]-1;
	int i2 = elface[f][2]-1;

	if (vertexnr[i1] < vertexnr[i0]) Swap (i0, i1);
	if (vertexnr[i2] < vertexnr[i0]) Swap (i0, i2);
	if (vertexnr[i2] < vertexnr[i1]) Swap (i1, i2);

	double arg = lambda(i1) + lambda(i2);
	double arg2 = lambda(i0) + lambda(i1) + lambda(i2);

	if (fabs(arg) < EPSILON || fabs(arg2) < EPSILON) // division by 0
	  {
	    for (int k = 0; k < nfaceshapes[f]; k++)
	      {
		fshape[index] = 0;
		fdshape[index] = Vec<3>(0,0,0);
		index++;
	      }
	    continue;
	  }

	b1.SetOrder (faceorder[f]);
	b2.SetOrder (faceorder[f]);
	    
	b1.CalcFDf (lambda(i0)/arg2);
	b2.CalcFDf (lambda(i2)/arg);
	    
	double decay = 1;
	double ddecay;
	    
	double l0 = lambda(i0);
	double l1 = lambda(i1);
	double l2 = lambda(i2);
	double dl0x = dlambda(i0,0);
	double dl0y = dlambda(i0,1);
	double dl0z = dlambda(i0,2);
	double dl1x = dlambda(i1,0);
	double dl1y = dlambda(i1,1);
	double dl1z = dlambda(i1,2);
	double dl2x = dlambda(i2,0);
	double dl2y = dlambda(i2,1);
	double dl2z = dlambda(i2,2);
	    
	double dargx = dl1x + dl2x;
	double dargy = dl1y + dl2y;
	double dargz = dl1z + dl2z;
	double darg2x = dl0x + dl1x + dl2x;
	double darg2y = dl0y + dl1y + dl2y;
	double darg2z = dl0z + dl1z + dl2z;

	for (int n1 = 2; n1 < faceorder[f]; n1++)
	  {
	    ddecay = (n1-1)*decay;
	    decay *= arg;
		
	    double decay2 = arg2;
	    double ddecay2;

	    for (int n0 = 2; n0 < faceorder[f]-n1+2; n0++)
	      {
		ddecay2 = n0*decay2;
		decay2 *= arg2;

		fshape[index] = b1.GetF(n0) * b2.GetF(n1) * decay * decay2;
		fdshape[index] = Vec<3>
		  (b1.GetDf(n0) * (dl0x * arg2 - l0 * darg2x)/(arg2*arg2) * b2.GetF(n1) * decay * decay2 +
		   b1.GetF(n0) * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargx * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * decay * ddecay2 * darg2x,
			
		   b1.GetDf(n0) * (dl0y * arg2 - l0 * darg2y)/(arg2*arg2) * b2.GetF(n1) * decay * decay2 +
		   b1.GetF(n0) * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargy * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * decay * ddecay2 * darg2y,

		   b1.GetDf(n0) * (dl0z * arg2 - l0 * darg2z)/(arg2*arg2) * b2.GetF(n1) * decay * decay2 +
		   b1.GetF(n0) * b2.GetDf(n1) * (dl2z * arg - l2 * dargz)/(arg*arg) * decay * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargz * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * decay * ddecay2 * darg2z);

		index++;
	      }
	  }
      }
  }




  // ----------------------------------------------------------------------------
  //      FEPrism
  // ----------------------------------------------------------------------------


  void FEPrism :: SetReferencePoint (Point<3> axi)
  {
    BaseFiniteElement3D :: SetReferencePoint (axi);
	
    lambda(0) = xi(0);
    lambda(1) = xi(1);
    lambda(2) = 1-xi(0)-xi(1);
    lambda(3) = xi(2);

    dlambda(0,0) =  1; dlambda(0,1) =  0; dlambda(0,2) =   0;
    dlambda(1,0) =  0; dlambda(1,1) =  1; dlambda(1,2) =   0;
    dlambda(2,0) = -1; dlambda(2,1) = -1; dlambda(2,2) =   0;
    dlambda(3,0) =  0; dlambda(3,1) =  0; dlambda(3,2) =   1;
  }


  void FEPrism :: CalcVertexShapes ()
  {
    double zcomp = 1-lambda(3);

    for (int v = 0; v < nvertices; v++)
      {
	if (v == 3) zcomp = 1-zcomp;

	vshape[v] = lambda(v % 3) * zcomp;
	vdshape[v](0) = dlambda(v % 3,0) * zcomp;
	vdshape[v](1) = dlambda(v % 3,1) * zcomp;
	vdshape[v](2) = lambda(v % 3) * (-dlambda(3,2));

	if (v >= 3) vdshape[v](2) *= -1;
      }
  }


  void FEPrism :: CalcEdgeShapes ()
  {
    int index = 0;
    int e;
    // triangle edge shapes
	
    for (e = 0; e < 6; e++)
      {
	int i0 = (eledge[e][0]-1) % 3;
	int i1 = (eledge[e][1]-1) % 3;

	double arg = lambda[i0]+lambda[i1];

	if (fabs(arg) < EPSILON) // division by 0
	  {
	    for (int k = 2; k <= edgeorder[e]; k++)
	      {
		eshape[index] = 0;
		edshape[index] = Vec<3>(0,0,0);
		index++;
	      }
	    continue;
	  }

	if (edgeorient[e] == -1) Swap (i0, i1);

	double xi = lambda[i1]/arg;

	b1.SetOrder (edgeorder[e]);
	b1.CalcFDf (xi);

	double decay = arg;
	double ddecay;

	double zarg = e < 3 ? (1-lambda(3)) : lambda(3);
	double zcomp = zarg;
	double dzcomp;
	    
	double l0 = lambda(i0);
	double dl0x = dlambda(i0,0);
	double dl0y = dlambda(i0,1);

	double l1 = lambda(i1);
	double dl1x = dlambda(i1,0);
	double dl1y = dlambda(i1,1);

	double dargx = dl0x + dl1x;
	double dargy = dl0y + dl1y;
                           

	/*

	for (int k = 2; k <= edgeorder[e]; k++)
	{        
	ddecay = k*decay;
	decay *= arg;

	dzcomp = k*zcomp;
	zcomp *= zarg;

	eshape[index] = b1.GetF (k) * decay * zcomp;
	edshape[index] = Vec<3>
	((b1.GetDf(k) * (dl1x*arg - l1*dargx) * 
	decay / (arg * arg) + b1.GetF(k) * ddecay * dargx) * zcomp,
	(b1.GetDf(k) * (dl1y*arg - l1*dargy) * 
	decay / (arg * arg) + b1.GetF(k) * ddecay * dargy) * zcomp,
	b1.GetF(k) * decay * dzcomp * (e < 3 ? -1 : 1));
	index++;
	}
	*/

	// JS linear in z-direction (for thin plates !)
	for (int k = 2; k <= edgeorder[e]; k++)
	  {        
	    ddecay = k*decay;
	    decay *= arg;

	    // dzcomp = k*zcomp;
	    // zcomp *= zarg;
	    dzcomp = 1;
	    zcomp = zarg;

	    eshape[index] = b1.GetF (k) * decay * zcomp;
	    edshape[index] = Vec<3>
	      ((b1.GetDf(k) * (dl1x*arg - l1*dargx) * 
		decay / (arg * arg) + b1.GetF(k) * ddecay * dargx) * zcomp,
	       (b1.GetDf(k) * (dl1y*arg - l1*dargy) * 
		decay / (arg * arg) + b1.GetF(k) * ddecay * dargy) * zcomp,
	       b1.GetF(k) * decay * dzcomp * (e < 3 ? -1 : 1));
	    index++;
	  }



      }

    // rectangle edge shapes
	
    for (e = 6; e < nedges; e++)
      {
	int i0 = eledge[e][0]-1;

	double arg = lambda[i0]; 

	if (fabs(arg) < EPSILON) // division by 0
	  {
	    for (int k = 2; k <= edgeorder[e]; k++)
	      {
		eshape[index] = 0.;
		edshape[index] = Vec<3>(0.,0.,0.);
		index++;
	      }
	    continue;
	  }

	double xi = lambda[3];

	if (edgeorient[e] == -1) xi = 1-xi;

	b1.SetOrder (edgeorder[e]);
	b1.CalcFDf (xi);

	double decay = arg;
	double ddecay;
	    
	double l0 = lambda(i0);
	double l0x = dlambda(i0,0);
	double l0y = dlambda(i0,1);

	for (int k = 2; k <= edgeorder[e]; k++)
	  {        
	    ddecay = k*decay;
	    decay *= arg;
		
	    eshape[index] = b1.GetF (k) * decay;
	    edshape[index] = Vec<3>
	      (b1.GetF(k) * ddecay * l0x,
	       b1.GetF(k) * ddecay * l0y,
	       b1.GetDf(k) * edgeorient[e] * decay);
	    index++;
	  }
      }
  }


  void FEPrism :: CalcFaceShapes ()
  {
    int index = 0;
    int f;

    // triangle face parts

    for (f = 0; f < 2; f++)
      {
	int i0 = elface[f][0]-1;
	int i1 = elface[f][1]-1;
	int i2 = elface[f][2]-1;

	if (vertexnr[i1] < vertexnr[i0]) Swap (i0, i1);
	if (vertexnr[i2] < vertexnr[i0]) Swap (i0, i2);
	if (vertexnr[i2] < vertexnr[i1]) Swap (i1, i2);

	i0 = i0 % 3;
	i1 = i1 % 3;
	i2 = i2 % 3;

	double arg = lambda(i1) + lambda(i2);

	if (fabs(arg) < EPSILON) // division by 0
	  {
	    for (int k = 0; k < nfaceshapes[f]; k++)
	      {
		fshape[index] = 0;
		fdshape[index] = Vec<3>(0,0,0);
		index++;
	      }
	    continue;
	  }

	b1.SetOrder (faceorder[f]);
	b2.SetOrder (faceorder[f]);
	    
	b1.CalcFDf (lambda(i0));
	b2.CalcFDf (lambda(i2)/arg);
	    
	double decay = 1;
	double ddecay;
	    
	double l0 = lambda(i0);
	double l1 = lambda(i1);
	double l2 = lambda(i2);
	double dl0x = dlambda(i0,0);
	double dl0y = dlambda(i0,1);
	double dl0z = dlambda(i0,2);
	double dl1x = dlambda(i1,0);
	double dl1y = dlambda(i1,1);
	double dl1z = dlambda(i1,2);
	double dl2x = dlambda(i2,0);
	double dl2y = dlambda(i2,1);
	double dl2z = dlambda(i2,2);
	    
	double dargx = dl1x + dl2x;
	double dargy = dl1y + dl2y;
	double dargz = dl1z + dl2z;

	double arg2 = f == 0 ? 1-xi(2) : xi(2);
	double darg2z = f == 0 ? -1 : 1;
	    
	for (int n1 = 2; n1 < faceorder[f]; n1++)
	  {
	    ddecay = (n1-1)*decay;
	    decay *= arg;
		
	    double decay2 = arg2;
	    double ddecay2;

	    for (int n0 = 2; n0 < faceorder[f]-n1+2; n0++)
	      {
		ddecay2 = n0*decay2;
		decay2 *= arg2;

		fshape[index] = b1.GetF(n0) * b2.GetF(n1) * decay * decay2;
		fdshape[index] = Vec<3>
		  (b1.GetDf(n0) * dl0x * b2.GetF(n1) * decay * decay2 +
		   b1.GetF(n0) * b2.GetDf(n1) * (dl2x * arg - l2 * dargx)/(arg*arg) * decay * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargx * decay2,
			
		   b1.GetDf(n0) * dl0y * b2.GetF(n1) * decay * decay2 +
		   b1.GetF(n0) * b2.GetDf(n1) * (dl2y * arg - l2 * dargy)/(arg*arg) * decay * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargy * decay2,

		   b1.GetDf(n0) * dl0z * b2.GetF(n1) * decay * decay2 +
		   b1.GetF(n0) * b2.GetDf(n1) * (dl2z * arg - l2 * dargz)/(arg*arg) * decay * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargz * decay2 +
		   b1.GetF(n0) * b2.GetF(n1) * decay * ddecay2 * darg2z);

		index++;
	      }
	  }
      }


    // quad parts

    for (f = 2; f < nfaces; f++)
      {
	// find index of point with smallest number
	  
	int i, i0 = 0;
	for (i = 1; i < 4; i++)
	  if (vertexnr[elface[f][i]-1] < vertexnr[elface[f][i0]-1]) i0 = i;
	    
	double arg = 0;
	double dargx = 0;
	double dargy = 0;
	double dargz = 0;
	for (i = 0; i < 4; i++)
	  {
	    arg += lambda((elface[f][i]-1) % 3)/2.0;
	    dargx += dlambda((elface[f][i]-1) % 3,0)/2.0;
	    dargy += dlambda((elface[f][i]-1) % 3,1)/2.0;
	    dargz += dlambda((elface[f][i]-1) % 3,2)/2.0;
	  }
	    
	if (fabs(arg) < EPSILON) // division by 0
	  {
	    for (int k = 0; k < nfaceshapes[f]; k++)
	      {
		fshape[index] = 0;
		fdshape[index] = Vec<3>(0,0,0);
		index++;
	      }
	    continue;
	  }

	int i1 = (i0+3) % 4;
	int j = (elface[f][i0]-1) % 3;

	double lam = lambda(j)/arg;
	double dlamx = (dlambda(j,0)*arg-lambda(j)*dargx)/(arg*arg);
	double dlamy = (dlambda(j,1)*arg-lambda(j)*dargy)/(arg*arg);
	double dlamz = (dlambda(j,2)*arg-lambda(j)*dargz)/(arg*arg);
			    
	double x;
	double z;
	double dxx;
	double dxy;
	double dxz;
	double dzx;
	double dzy;
	double dzz;

	int ratvar;
	/*
	  switch (i0)
	  {
	  case 0:
	  x = xi(2); z = lam;

	  dxx = 0;     dxy = 0;     dxz = 1;
	  dzx = dlamx; dzy = dlamy; dzz = dlamz;
	  ratvar = 1;
	  break;
	  case 1:
	  x = 1-lam; z = xi(2);
	  dxx = -dlamx; dxy = -dlamy; dxz = -dlamz;
	  dzx = 0;      dzy = 0;      dzz = 1;
	  ratvar = 0;
	  break;
	  case 2: 
	  x = 1-xi(2); z = 1-lam;
	  dxx = 0;      dxy = 0;      dxz = -1;
	  dzx = -dlamx; dzy = -dlamy; dzz = -dlamz;
	  ratvar = 1;
	  break;
	  case 3:
	  x = lam; z = 1-xi(2);
	  dxx = dlamx; dxy = dlamy; dxz = dlamz;
	  dzx = 0;     dzy = 0;     dzz = -1;
	  ratvar = 0;
	  break;
	  }
	*/

	ratvar = 0;
	x = 1-lam;
	dxx = -dlamx; dxy = -dlamy; dxz = -dlamz;
	if (i0 == 0 || i0 == 1)
	  {
	    z = xi(2);
	    dzx = 0; dzy = 0; dzz = 1;
	  }
	else
	  {
	    z = 1-xi(2);
	    dzx = 0; dzy = 0; dzz = -1;
	  }

	int ix = i0 ^ 1;
	int iz = 3-i0;

	if (vertexnr[elface[f][ix]-1] > vertexnr[elface[f][iz]-1])
	  {
	    Swap (x,z);
	    Swap (dxx, dzx);
	    Swap (dxy, dzy);
	    Swap (dxz, dzz);
	    ratvar = 1-ratvar;
	  }

	b1.SetOrder (faceorder[f]);
	b2.SetOrder (faceorder[f]);
	    
	b1.CalcFDf (x);
	b2.CalcFDf (z);
	    
	double decay = arg;
	double ddecay;
	    
	for (int n0 = 2; n0 <= faceorder[f]; n0++)
	  {
	    ddecay = n0*decay;
	    decay *= arg;
		
	    if (ratvar == 1) decay = arg;

	    for (int n1 = 2; n1 <= faceorder[f]; n1++)
	      {
		if (ratvar == 1)
		  {
		    ddecay = n1*decay;
		    decay *= arg;
		  }

		fshape[index] = b1.GetF(n0) * b2.GetF(n1) * decay;
		fdshape[index] = Vec<3>
		  (b1.GetDf(n0) * dxx * b2.GetF(n1) * decay +
		   b1.GetF(n0) * b2.GetDf(n1) * dzx * decay +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargx,

		   b1.GetDf(n0) * dxy * b2.GetF(n1) * decay +
		   b1.GetF(n0) * b2.GetDf(n1) * dzy * decay +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargy,
			
		   b1.GetDf(n0) * dxz * b2.GetF(n1) * decay +
		   b1.GetF(n0) * b2.GetDf(n1) * dzz * decay +
		   b1.GetF(n0) * b2.GetF(n1) * ddecay * dargz);

		index++;
	      }
	  }
      }
  }
    


  // ----------------------------------------------------------------------------
  //      FEPyramid
  // ----------------------------------------------------------------------------


  void FEPyramid :: SetReferencePoint (Point<3> axi)
  {
    BaseFiniteElement3D :: SetReferencePoint (axi);
  }


  void FEPyramid :: CalcVertexShapes ()
  {
    double x = xi(0);
    double y = xi(1);
    double z = xi(2);

    if (z == 1.) z = 1-1e-10;
    vshape[0] = (1-z-x)*(1-z-y) / (1-z);
    vshape[1] = x*(1-z-y) / (1-z);
    vshape[2] = x*y / (1-z);
    vshape[3] = (1-z-x)*y / (1-z);
    vshape[4] = z;

    double z1 = 1-z;
    double z2 = z1*z1;

    vdshape[0] = Vec<3>( -(z1-y)/z1, -(z1-x)/z1, ((x+y+2*z-2)*z1+(z1-y)*(z1-x))/z2 );
    vdshape[1] = Vec<3>( (z1-y)/z1,  -x/z1,      (-x*z1+x*(z1-y))/z2 );
    vdshape[2] = Vec<3>( y/z1,       x/z1,       x*y/z2 );
    vdshape[3] = Vec<3>( -y/z1,      (z1-x)/z1,  (-y*z1+y*(z1-x))/z2 );
    vdshape[4] = Vec<3>( 0, 0, 1 );
  }


  void FEPyramid :: CalcEdgeShapes ()
  {
    int index = 0;

    for (int e = 0; e < GetNEdges(); e++)
      {
	for (int k = 2; k <= edgeorder[e]; k++)
	  {        
	    eshape[index] = 0;
	    edshape[index] = Vec<3>(0,0,0);
	    index++;
	  }
      }
  }


  void FEPyramid :: CalcFaceShapes ()
  {
    int index = 0;

    for (int f = 0; f < GetNFaces(); f++)
      {
	for (int k = 0; k < nfaceshapes[f]; k++)
	  {
	    fshape[index] = 0;
	    fdshape[index] = Vec<3>(0,0,0);
	    index++;
	  }
      }
  }
    




  // ----------------------------------------------------------------------------
  //      FEHex
  // ----------------------------------------------------------------------------


  void FEHex :: SetReferencePoint (Point<3> axi)
  {
    BaseFiniteElement3D :: SetReferencePoint (axi);
  }


  void FEHex :: CalcVertexShapes ()
  {
    double x = xi(0);
    double y = xi(1);
    double z = xi(2);

    vshape[0] = (1-x)*(1-y)*(1-z);
    vshape[1] =    x *(1-y)*(1-z); 
    vshape[2] =    x *   y *(1-z);
    vshape[3] = (1-x)*   y *(1-z);
    vshape[4] = (1-x)*(1-y)* z;
    vshape[5] =    x *(1-y)* z;
    vshape[6] =    x *   y * z;
    vshape[7] = (1-x)*   y * z;

    vdshape[0] = Vec<3>(-(1-y)*(1-z), -(1-x)*(1-z), -(1-x)*(1-y));
    vdshape[1] = Vec<3>( (1-y)*(1-z),    -x *(1-z),    -x *(1-y));
    vdshape[2] = Vec<3>(    y *(1-z),     x *(1-z),    -x * y);
    vdshape[3] = Vec<3>(   -y *(1-z),  (1-x)*(1-z), -(1-x)*y);
    vdshape[4] = Vec<3>(-(1-y)*   z, -(1-x)* z,  (1-x)*(1-y));
    vdshape[5] = Vec<3>( (1-y)*   z,    -x * z,     x *(1-y));
    vdshape[6] = Vec<3>(    y *   z,     x * z,     x * y);
    vdshape[7] = Vec<3>(   -y *   z,  (1-x)* z,  (1-x)*y);

  }


  void FEHex :: CalcEdgeShapes ()
  {
    int index = 0;

    for (int e = 0; e < GetNEdges(); e++)
      {
	for (int k = 2; k <= edgeorder[e]; k++)
	  {        
	    eshape[index] = 0;
	    edshape[index] = Vec<3>(0,0,0);
	    index++;
	  }
      }
  }


  void FEHex :: CalcFaceShapes ()
  {
    int index = 0;

    for (int f = 0; f < GetNFaces(); f++)
      {
	for (int k = 0; k < nfaceshapes[f]; k++)
	  {
	    fshape[index] = 0;
	    fdshape[index] = Vec<3>(0,0,0);
	    index++;
	  }
      }
  }
    








  int CurvedElements :: IsSurfaceElementCurved (int elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(SurfaceElementIndex) elnr].hp_elnr];

	return mesh.coarsemesh->GetCurvedElements().IsSurfaceElementCurved (hpref_el.coarse_elnr);
      }




    Element2d elem = mesh[(SurfaceElementIndex) elnr];

    switch (elem.GetType())
      {
      case TRIG:
	{
	  FETrig fe2d(*this);
	  fe2d.SetElementNumber (elnr+1);
	  return (fe2d.IsCurved());
	  break;
	}

      case QUAD:
	{
	  FEQuad fe2d(*this);
	  fe2d.SetElementNumber (elnr+1);
	  return (fe2d.IsCurved());
	  break;
	}

      }
    return 0;
  }



  int CurvedElements :: IsElementCurved (int elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(ElementIndex) elnr].hp_elnr];

	return mesh.coarsemesh->GetCurvedElements().IsElementCurved (hpref_el.coarse_elnr);
      }



    Element elem = mesh[(ElementIndex) elnr];

    switch (elem.GetType())
      {
      case TET:
	{
	  FETet fe3d(*this);
	  fe3d.SetElementNumber (elnr+1);
	  return (fe3d.IsCurved());
	  break;
	}

      case PRISM:
	{
	  FEPrism fe3d(*this);
	  fe3d.SetElementNumber (elnr+1);
	  return (fe3d.IsCurved());
	  break;
	}

      case PYRAMID:
	{
	  FEPyramid fe3d(*this);
	  fe3d.SetElementNumber (elnr+1);
	  return (fe3d.IsCurved());
	  break;
	}

      case HEX:
	{
	  FEHex fe3d(*this);
	  fe3d.SetElementNumber (elnr+1);
	  return (fe3d.IsCurved());
	  break;
	}

      }

    return 0;
  }


  void CurvedElements :: CalcSegmentTransformation (double xi, int segnr,
						    Point<3> * x, Vec<3> * dxdxi)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(SegmentIndex) segnr].hp_elnr];
	  
	// xi umrechnen
	double lami[2];
	lami[0] = xi;
	lami[1] = 1-xi;

	double coarse_xi = 0;
	for (int i = 0; i < 2; i++)
	  coarse_xi += hpref_el.param[i][0] * lami[i];

	mesh.coarsemesh->GetCurvedElements().CalcSegmentTransformation (coarse_xi, hpref_el.coarse_elnr, x, dxdxi);
	return;
      }


    // xi = 1-xi;  // or, this way  ????? JS
    FESegm segm (*this);
    segm.SetElementNumber (segnr+1);
    segm.SetReferencePoint (Point<1>(xi));

    //	segm.CalcVertexShapes (x != NULL, dxdxi != NULL);
    segm.CalcVertexShapes ();

    if (x)
      {
	(*x) = Point<3>(0,0,0);
	for (int v = 0; v < 2; v++)
	  (*x) = (*x) + segm.GetVertexShape(v) * mesh.Point(segm.GetVertexNr(v));
      }
	
    if (dxdxi)
      {
	(*dxdxi) = Vec<3>(0,0,0);
	for (int v = 0; v < 2; v++)
	  (*dxdxi) = (*dxdxi) + segm.GetVertexDShape(v) * mesh.Point(segm.GetVertexNr(v));
      }
	
    if (segm.GetEdgeOrder() > 1)
      {
	//	    segm.CalcEdgeShapes (x != NULL, dxdxi != NULL);
	segm.CalcEdgeShapes ();
	    
	if (x)
	  {
	    int gindex = edgecoeffsindex[segm.GetEdgeNr()-1];
	    for (int k = 2; k <= segm.GetEdgeOrder(); k++, gindex++)
	      (*x) = (*x) + segm.GetEdgeShape(k-2) * edgecoeffs[gindex];
	  }
	    
	if (dxdxi)
	  {
	    int gindex = edgecoeffsindex[segm.GetEdgeNr()-1];
	    for (int k = 2; k <= segm.GetEdgeOrder(); k++, gindex++)
	      (*dxdxi) = (*dxdxi) + segm.GetEdgeDShape(k-2) * edgecoeffs[gindex];
	  }
      }
  }
  


  void CurvedElements :: CalcMultiPointSegmentTransformation (ARRAY<double> * xi, int segnr,
							      ARRAY<Point<3> > * x,
							      ARRAY<Vec<3> > * dxdxi)
  {
    int size = xi->Size();

    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(SegmentIndex) segnr].hp_elnr];
	  
	// xi umrechnen
	ARRAY< Point<2> > lami;
	lami.SetSize (size);
	for (int i = 0; i<size; i++)
	  {
	    lami[i](0) = (*xi)[i];
	    lami[i](1) = 1-(*xi)[i];
	  }

	ARRAY<double> coarse_xi;
	coarse_xi.SetSize (size);
	coarse_xi = 0;

	  
	for (int j = 0; j < 2; j++)
	  for (int i = 0; i<size; i++)
	    coarse_xi[i] += hpref_el.param[j][0] * lami[i](j);

	mesh.coarsemesh->GetCurvedElements().CalcMultiPointSegmentTransformation
	  (&coarse_xi, hpref_el.coarse_elnr, x, dxdxi);
	return;
      }



    FESegm segm (*this);
    segm.SetElementNumber (segnr+1);
    x->SetSize (size);
    dxdxi->SetSize (size);

    for (int i = 0; i < size; i++)
      {
	segm.SetReferencePoint (Point<1>((*xi)[i]));

	//	segm.CalcVertexShapes (x != NULL, dxdxi != NULL);
	segm.CalcVertexShapes ();

	(*x)[i] = Point<3>(0,0,0);
	(*dxdxi)[i] = Vec<3>(0,0,0);

	for (int v = 0; v < 2; v++)
	  {
	    (*x)[i] = (*x)[i] + segm.GetVertexShape(v) * mesh.Point(segm.GetVertexNr(v));
	    (*dxdxi)[i] = (*dxdxi)[i] + segm.GetVertexDShape(v) * mesh.Point(segm.GetVertexNr(v));
	  }

	if (segm.GetEdgeOrder() > 1)
	  {
	    //	    segm.CalcEdgeShapes (x != NULL, dxdxi != NULL);
	    segm.CalcEdgeShapes ();

	    int gindex = edgecoeffsindex[segm.GetEdgeNr()-1];
	    for (int k = 2; k <= segm.GetEdgeOrder(); k++, gindex++)
	      {
		(*x)[i] = (*x)[i] + segm.GetEdgeShape(k-2) * edgecoeffs[gindex];
		(*dxdxi)[i] = (*dxdxi)[i] + segm.GetEdgeDShape(k-2) * edgecoeffs[gindex];
	      }
	  }
      }
  }





  void CurvedElements :: CalcSurfaceTransformation (Point<2> xi, int elnr,
						    Point<3> * x, Mat<3,2> * dxdxi)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(SurfaceElementIndex) elnr].hp_elnr];
	  
	// xi umrechnen
	double lami[4];
	FlatVector vlami(4, lami);
	vlami = 0;
	mesh[(SurfaceElementIndex) elnr] . GetShapeNew (xi, vlami);

	Mat<2,2> trans;
	Mat<3,2> dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<2> dlami(4);
	    dlami = 0;
	    mesh[(SurfaceElementIndex) elnr] . GetDShapeNew (xi, dlami);	  
	      
	    trans = 0;
	    for (int k = 0; k < 2; k++)
	      for (int l = 0; l < 2; l++)
		{
		  double sum = 0;
		  for (int i = 0; i < 4; i++)
		    sum += hpref_el.param[i][l] * dlami(i, k);
		  trans(l,k) = sum;
		}
	  }

	Point<2> coarse_xi(0,0);
	for (int i = 0; i < 4; i++)
	  {
	    coarse_xi(0) += hpref_el.param[i][0] * lami[i];
	    coarse_xi(1) += hpref_el.param[i][1] * lami[i];
	  }
	mesh.coarsemesh->GetCurvedElements().CalcSurfaceTransformation (coarse_xi, hpref_el.coarse_elnr, x, &dxdxic);

	if (dxdxi)
	  *dxdxi = dxdxic * trans;

	return;
      }




    const Element2d & elem = mesh[(SurfaceElementIndex) elnr];

    BaseFiniteElement2D * fe2d;

    // char locmem[max2(sizeof(FEQuad), sizeof(FETrig))];
    char locmemtrig[sizeof(FETrig)];
    char locmemquad[sizeof(FEQuad)];
    switch (elem.GetType())
      {
      case TRIG: fe2d = new (locmemtrig) FETrig (*this); break;
      case QUAD: fe2d = new (locmemquad) FEQuad (*this); break;
      }

    //	fe2d->SetElementNumber (elnr+1);
    fe2d->SetReferencePoint (xi);
    fe2d->CalcVertexShapes ();

    if (x)
      {
	(*x) = Point<3>(0,0,0);
	for (int v = 0; v < fe2d->GetNVertices(); v++)
	  (*x) = (*x) + fe2d->GetVertexShape(v) * mesh.Point(elem[v]);
	// (*x) = (*x) + fe2d->GetVertexShape(v) * mesh.Point(fe2d->GetVertexNr(v));
      }

    if (dxdxi)
      {
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 2; j++)
	    (*dxdxi)(i,j) = 0;
                  
	for (int v = 0; v < elem.GetNV(); v++)
	  for (int i = 0; i < 3; i++)
	    for (int j = 0; j < 2; j++)
	      (*dxdxi)(i,j) += fe2d->GetVertexDShape(v)(j) * mesh.Point(elem[v]).X(i+1);
	// (*dxdxi)(i,j) += fe2d->GetVertexDShape(v)(j) * mesh.Point(fe2d->GetVertexNr(v)).X(i+1);
      }

    if (IsHighOrder())
      {
	fe2d->SetElementNumber (elnr+1);

	//	    fe2d->CalcEdgeShapes (x != NULL, dxdxi != NULL);
	fe2d->CalcEdgeShapes ();
	if (x)
	  {
	    int index = 0;
	    for (int e = 0; e < fe2d->GetNEdges(); e++)
	      {
		int gindex = edgecoeffsindex[fe2d->GetEdgeNr(e)-1];

		for (int k = 2; k <= fe2d->GetEdgeOrder(e); k++, index++, gindex++)
		  (*x) = (*x) + fe2d->GetEdgeShape(index) * edgecoeffs[gindex];
	      }
	  }

	if (dxdxi)
	  {
	    int index = 0;
	    for (int e = 0; e < fe2d->GetNEdges(); e++)
	      {
		int gindex = edgecoeffsindex[fe2d->GetEdgeNr(e)-1];
		for (int k = 2; k <= fe2d->GetEdgeOrder(e); k++, index++, gindex++)
		  for (int i = 0; i < 3; i++)
		    for (int j = 0; j < 2; j++)
		      (*dxdxi)(i,j) += fe2d->GetEdgeDShape(index)(j) * edgecoeffs[gindex](i);
	      }
	  }

	if (mesh.GetDimension() == 3)
	  {
	    //		fe2d->CalcFaceShapes (x != NULL, dxdxi != NULL);
	    fe2d->CalcFaceShapes ();

	    if (x)
	      {
		int gindex = facecoeffsindex[fe2d->GetFaceNr()-1];
		for (int index = 0; index < fe2d->GetNFaceShapes(); index++, gindex++)
		  {
		    (*x) = (*x) + fe2d->GetFaceShape(index) * facecoeffs[gindex];
		  }
	      }

	    if (dxdxi)
	      {
		int gindex = facecoeffsindex[fe2d->GetFaceNr()-1];
		for (int index = 0; index < fe2d->GetNFaceShapes(); index++, gindex++)
		  for (int i = 0; i < 3; i++)
		    for (int j = 0; j < 2; j++)
		      (*dxdxi)(i,j) += fe2d->GetFaceDShape(index)(j) * facecoeffs[gindex](i);
	      }
	  }
      } 

    fe2d -> ~BaseFiniteElement2D();
  }	








  void CurvedElements :: CalcMultiPointSurfaceTransformation (ARRAY< Point<2> > * xi, int elnr,
							      ARRAY< Point<3> > * x,
							      ARRAY< Mat<3,2> > * dxdxi)
  {

    int size = xi->Size();
      
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(SurfaceElementIndex) elnr].hp_elnr];
	  
	// xi umrechnen
	ARRAY< Point<2> > coarse_xi;
	ARRAY< Mat<2,2> > trans;
	ARRAY< Mat<3,2> > dxdxic;

	coarse_xi.SetSize (size);
	coarse_xi = 0;
	trans.SetSize (size);
	dxdxic.SetSize (size);

	for (int mp = 0; mp<size; mp++)
	  {
	    double lami[4];
	    FlatVector vlami(4, lami);
	    vlami = 0;
	    mesh[(SurfaceElementIndex) elnr] . GetShapeNew ((*xi)[mp], vlami);
	      
	    MatrixFixWidth<2> dlami(4);
	    dlami = 0;
	    mesh[(SurfaceElementIndex) elnr] . GetDShapeNew ((*xi)[mp], dlami);	  
	      
	    trans[mp] = 0;
	    for (int k = 0; k < 2; k++)
	      for (int l = 0; l < 2; l++)
		{
		  double sum = 0;
		  for (int i = 0; i < 4; i++)
		    sum += hpref_el.param[i][l] * dlami(i, k);
		  trans[mp](l,k) = sum;
		}
	      
	    for (int i = 0; i < 4; i++)
	      {
		coarse_xi[mp](0) += hpref_el.param[i][0] * lami[i];
		coarse_xi[mp](1) += hpref_el.param[i][1] * lami[i];
	      }
	  }
	  
	mesh.coarsemesh->GetCurvedElements().CalcMultiPointSurfaceTransformation (&coarse_xi, hpref_el.coarse_elnr, x, &dxdxic);
	  
	for (int mp = 0; mp < size; mp++)
	  {
	    (*dxdxi)[mp] = dxdxic[mp] * trans[mp];
	  }
	  
	return;
      }
    

    x->SetSize (size);
    dxdxi->SetSize (size);

    const Element2d & elem = mesh[(SurfaceElementIndex) elnr];

    BaseFiniteElement2D * fe2d;

    // char locmem[max2(sizeof(FEQuad), sizeof(FETrig))];
    char locmemtrig[sizeof(FETrig)];
    char locmemquad[sizeof(FEQuad)];
    switch (elem.GetType())
      {
      case TRIG: fe2d = new (locmemtrig) FETrig (*this); break;
      case QUAD: fe2d = new (locmemquad) FEQuad (*this); break;
      }
      
    fe2d->SetElementNumber (elnr+1);

    for (int mp = 0; mp < size; mp++)
      {
	fe2d->SetReferencePoint ((*xi)[mp]);
	fe2d->CalcVertexShapes ();
      
	(*x)[mp] = Point<3>(0,0,0);
      
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 2; j++)
	    (*dxdxi)[mp](i,j) = 0;
      
	for (int v = 0; v < fe2d->GetNVertices(); v++)
	  {
	    (*x)[mp] = (*x)[mp] + fe2d->GetVertexShape(v) * mesh.Point(fe2d->GetVertexNr(v));
	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 2; j++)
		(*dxdxi)[mp](i,j) += fe2d->GetVertexDShape(v)(j) * mesh.Point(fe2d->GetVertexNr(v)).X(i+1);
	  }
	  
	if (IsHighOrder())
	  {
	    //	    fe2d->CalcEdgeShapes (x != NULL, dxdxi != NULL);
	    fe2d->CalcEdgeShapes ();

	    int index = 0;
	    for (int e = 0; e < fe2d->GetNEdges(); e++)
	      {
		int gindex = edgecoeffsindex[fe2d->GetEdgeNr(e)-1];

		for (int k = 2; k <= fe2d->GetEdgeOrder(e); k++, index++, gindex++)
		  {
		    (*x)[mp] = (*x)[mp] + fe2d->GetEdgeShape(index) * edgecoeffs[gindex];
		    for (int i = 0; i < 3; i++)
		      for (int j = 0; j < 2; j++)
			(*dxdxi)[mp](i,j) += fe2d->GetEdgeDShape(index)(j) * edgecoeffs[gindex](i);	    
		  }
	      }


	    if (mesh.GetDimension() == 3)
	      {
		//		fe2d->CalcFaceShapes (x != NULL, dxdxi != NULL);
		fe2d->CalcFaceShapes ();

		int gindex = facecoeffsindex[fe2d->GetFaceNr()-1];
		for (int index = 0; index < fe2d->GetNFaceShapes(); index++, gindex++)
		  {
		    (*x)[mp] = (*x)[mp] + fe2d->GetFaceShape(index) * facecoeffs[gindex];
		    for (int i = 0; i < 3; i++)
		      for (int j = 0; j < 2; j++)
			(*dxdxi)[mp](i,j) += fe2d->GetFaceDShape(index)(j) * facecoeffs[gindex](i);
		  }
	      }
	  } 
      }
      
    fe2d -> ~BaseFiniteElement2D();
  }	



















  void CurvedElements :: CalcElementTransformation (Point<3> xi, int elnr,
						    Point<3> * x, Mat<3,3> * dxdxi)
  {

    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(ElementIndex) elnr].hp_elnr];
	

	/*
	 *testout << " Curved Element " << elnr << endl; 
	 *testout << " hpref_el.coarse_elnr " << hpref_el.coarse_elnr << endl; 
	 *testout << " hpref_el.param = " << endl; 
	 for(int j=0;j< hpref_el.np; j++) 
	 { 
	 for(int k=0;k<3;k++)
	 *testout << hpref_el.param[j][k] << "\t"; 
	 *testout << endl; 
	 } 
	*/

	double lami[8];
	FlatVector vlami(8, lami);
	vlami = 0;
	mesh[(ElementIndex) elnr] . GetShapeNew (xi, vlami);
	
	Point<3> coarse_xi(0,0,0);
	for (int i = 0; i < hpref_el.np; i++)
	  for (int l = 0; l < 3; l++)
	    coarse_xi(l) += hpref_el.param[i][l] * lami[i];
	  
	Mat<3,3> trans, dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<3> dlami(8);
	    dlami = 0;
	    mesh[(ElementIndex) elnr] . GetDShapeNew (xi, dlami);	  
	      
	    trans = 0;
	    for (int k = 0; k < 3; k++)
	      for (int l = 0; l < 3; l++)
		{
		  double sum = 0;
		  for (int i = 0; i < hpref_el.np; i++)
		    sum += hpref_el.param[i][l] * dlami(i, k);
		  trans(l,k) = sum;
		}
	  }
	/*	  
	 *testout << " x " << x << endl;
	 // << "\t" << x.X(2) << "\t" << x.X(3) << endl;  
	 *testout << " Element Trafo " << coarse_xi(0) << " \t " << coarse_xi(1) << " \t " << coarse_xi(2) << endl;
	 */
	 


	mesh.coarsemesh->GetCurvedElements().CalcElementTransformation (coarse_xi, hpref_el.coarse_elnr, x, &dxdxic);

	if (dxdxi)
	  *dxdxi = dxdxic * trans;
	return;
      }









    Element elem = mesh[(ElementIndex) elnr];
    BaseFiniteElement3D * fe3d;

    // char locmem[max2(sizeof(FETet), sizeof(FEPrism))];
    char locmemtet[sizeof(FETet)];
    char locmemprism[sizeof(FEPrism)];
    char locmempyramid[sizeof(FEPyramid)];
    char locmemhex[sizeof(FEHex)];
    switch (elem.GetType())
      {
      case TET: fe3d = new (locmemtet) FETet (*this); break;
      case PRISM: fe3d = new (locmemprism) FEPrism (*this); break;
      case PYRAMID: fe3d = new (locmempyramid) FEPyramid (*this); break;
      case HEX: fe3d = new (locmemhex) FEHex (*this); break;
      }
	
    // fe3d->SetElementNumber (elnr+1);
    fe3d->SetReferencePoint (xi);

    fe3d->CalcVertexShapes ();
    //	fe3d->CalcVertexShapes (x != NULL, dxdxi != NULL);

	
    if (x)
      {
	(*x) = Point<3>(0,0,0);
	for (int v = 0; v < fe3d->GetNVertices(); v++)
	  (*x) += fe3d->GetVertexShape(v) * Vec<3> (mesh.Point(elem[v]));
      }

    if (dxdxi)
      {
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    (*dxdxi)(i,j) = 0;
	    
	for (int v = 0; v < fe3d->GetNVertices(); v++)
	  for (int i = 0; i < 3; i++)
	    for (int j = 0; j < 3; j++)
	      (*dxdxi)(i,j) += fe3d->GetVertexDShape(v)(j) * mesh.Point(elem[v]).X(i+1);
      }

    if (IsHighOrder())
      {
	fe3d->SetElementNumber (elnr+1);
	//	    fe3d->CalcEdgeShapes (x != NULL, dxdxi != NULL);
	fe3d->CalcEdgeShapes ();


	if (x)
	  {
	    int index = 0;
	    for (int e = 0; e < fe3d->GetNEdges(); e++)
	      {
		int gindex = edgecoeffsindex[fe3d->GetEdgeNr(e)-1];
		for (int k = 2; k <= fe3d->GetEdgeOrder(e); k++, index++, gindex++)
		  (*x) += fe3d->GetEdgeShape(index) * edgecoeffs[gindex];
	      }
	  }

	if (dxdxi)
	  {
	    int index = 0;
	    for (int e = 0; e < fe3d->GetNEdges(); e++)   
	      {
		int gindex = edgecoeffsindex[fe3d->GetEdgeNr(e)-1];
		for (int k = 2; k <= fe3d->GetEdgeOrder(e); k++, index++, gindex++)
		  for (int i = 0; i < 3; i++)
		    for (int j = 0; j < 3; j++)
		      (*dxdxi)(i,j) += fe3d->GetEdgeDShape(index)(j) * edgecoeffs[gindex](i);
	      }
	  }


	// cout << "switched off face mapping" << endl;
	if (mesh.GetDimension() == 3)   // JS
	  {
	    fe3d->CalcFaceShapes ();
	    //		fe3d->CalcFaceShapes (x != NULL, dxdxi != NULL);

	    if (x)
	      {
		int index = 0;
		for (int f = 0; f < fe3d->GetNFaces(); f++)
		  {
		    int gindex = facecoeffsindex[fe3d->GetFaceNr(f)-1];
		    for (int k = 0; k < fe3d->GetNFaceShapes(f); k++, index++, gindex++)
		      (*x) += fe3d->GetFaceShape(index) * facecoeffs[gindex];
		  }
	      }

	    if (dxdxi)
	      {
		int index = 0;
		for (int f = 0; f < fe3d->GetNFaces(); f++)
		  {
		    int gindex = facecoeffsindex[fe3d->GetFaceNr(f)-1];
		    for (int k = 0; k < fe3d->GetNFaceShapes(f); k++, index++, gindex++)
		      for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			  (*dxdxi)(i,j) += fe3d->GetFaceDShape(index)(j) * facecoeffs[gindex](i);
		  }
	      }
	  } 
      }
	
    fe3d -> ~BaseFiniteElement3D();
  }









#ifdef ROBERT

  void CurvedElements :: CalcMultiPointElementTransformation (ARRAY< Point<3> > * xi, int elnr,
							      ARRAY< Point<3> > * x,
							      ARRAY< Mat<3,3> > * dxdxi)
  {
    int size = (*xi).Size();
    x->SetSize (size);

    if (dxdxi) dxdxi->SetSize (size);

    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(ElementIndex) elnr].hp_elnr];
	
	ARRAY< Mat<3,3> > trans, dxdxic;
	if (dxdxi)
	  {
	    trans.SetSize (size);
	    dxdxic.SetSize (size);
	  }

	ARRAY<Point<3> > coarse_xi(size);

	double lami[8];
	FlatVector vlami(8, lami);
	const Element & el = mesh[(ElementIndex) elnr];
	int el_np = el.GetNP();

	for (int mp = 0; mp < size; mp++)
	  {
	    el.GetShapeNew ((*xi)[mp], vlami);
	    
	    Point<3> pc(0,0,0);
	    for (int i = 0; i < el_np; i++)
	      for (int l = 0; l < 3; l++)
		pc(l) += hpref_el.param[i][l] * lami[i];
	
	    coarse_xi[mp] = pc;

	    if (dxdxi)
	      {
		MatrixFixWidth<3> dlami(8);
		dlami = 0;
		mesh[(ElementIndex) elnr] . GetDShapeNew ((*xi)[mp], dlami);	  
		
		trans[mp] = 0;
		for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 3; l++)
		    {
		      double sum = 0;
		      for (int i = 0; i < 8; i++)
			sum += hpref_el.param[i][l] * dlami(i, k);
		      trans[mp](l,k) = sum;
		    }
	      }
	  }

	if (dxdxi)
	  mesh.coarsemesh->GetCurvedElements().CalcMultiPointElementTransformation (&coarse_xi, hpref_el.coarse_elnr, x, &dxdxic);
	else
	  mesh.coarsemesh->GetCurvedElements().CalcMultiPointElementTransformation (&coarse_xi, hpref_el.coarse_elnr, x, 0);
	
	if (dxdxi)
	  for (int mp = 0; mp < size; mp++)
	    (*dxdxi)[mp] = dxdxic[mp] * trans[mp];
	
	return;
      }
    








    Element elem = mesh[(ElementIndex) elnr];
    BaseFiniteElement3D * fe3d;
    
    // char locmem[max2(sizeof(FETet), sizeof(FEPrism))];
    char locmemtet[sizeof(FETet)];
    char locmemprism[sizeof(FEPrism)];
    char locmempyramid[sizeof(FEPyramid)];
    char locmemhex[sizeof(FEHex)];
    switch (elem.GetType())
      {
      case TET: fe3d = new (locmemtet) FETet (*this); break;
      case PRISM: fe3d = new (locmemprism) FEPrism (*this); break;
      case PYRAMID: fe3d = new (locmempyramid) FEPyramid (*this); break;
      case HEX: fe3d = new (locmemhex) FEHex (*this); break;
      }
    
    fe3d->SetElementNumber (elnr+1);


    for (int mp = 0; mp < size; mp++)
      {
	fe3d->SetReferencePoint ((*xi)[mp]);
    
	fe3d->CalcVertexShapes ();
	//	fe3d->CalcVertexShapes (x != NULL, dxdxi != NULL);
    
	(*x)[mp] = Point<3>(0,0,0);
	if (dxdxi)
	  for (int i = 0; i < 3; i++)
	    for (int j = 0; j < 3; j++)
	      (*dxdxi)[mp](i,j) = 0;
	
	for (int v = 0; v < fe3d->GetNVertices(); v++)
	  {
	    (*x)[mp] += fe3d->GetVertexShape(v) * Vec<3> (mesh.Point(fe3d->GetVertexNr(v)));
	    if (dxdxi)
	      for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		  (*dxdxi)[mp](i,j) += fe3d->GetVertexDShape(v)(j) * mesh.Point(fe3d->GetVertexNr(v)).X(i+1);
	  }
    

	if (IsHighOrder())
	  {
	    //	    fe3d->CalcEdgeShapes (x != NULL, dxdxi != NULL);
	    fe3d->CalcEdgeShapes ();
	
	    int index = 0;
	    for (int e = 0; e < fe3d->GetNEdges(); e++)
	      {
		int gindex = edgecoeffsindex[fe3d->GetEdgeNr(e)-1];
		for (int k = 2; k <= fe3d->GetEdgeOrder(e); k++, index++, gindex++)
		  {
		    (*x)[mp] += fe3d->GetEdgeShape(index) * edgecoeffs[gindex];
		    if (dxdxi)
		      for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			  (*dxdxi)[mp](i,j) += fe3d->GetEdgeDShape(index)(j) * edgecoeffs[gindex](i);
		  }
	      }

	    if (mesh.GetDimension() == 3)
	      {
		fe3d->CalcFaceShapes ();
		//		fe3d->CalcFaceShapes (x != NULL, dxdxi != NULL);
	    
		int index = 0;
		for (int f = 0; f < fe3d->GetNFaces(); f++)
		  {
		    int gindex = facecoeffsindex[fe3d->GetFaceNr(f)-1];
		    for (int k = 0; k < fe3d->GetNFaceShapes(f); k++, index++, gindex++)
		      {
			(*x)[mp] += fe3d->GetFaceShape(index) * facecoeffs[gindex];
			if (dxdxi)
			  for (int i = 0; i < 3; i++)
			    for (int j = 0; j < 3; j++)
			      (*dxdxi)[mp](i,j) += fe3d->GetFaceDShape(index)(j) * facecoeffs[gindex](i);
		      }
		  }
	      }
	  }
      }
	
    fe3d -> ~BaseFiniteElement3D();
  }

#endif






  void CurvedElements :: CalcMultiPointElementTransformation (ARRAY< Point<3> > * xi, int elnr,
							      ARRAY< Point<3> > * x,
							      ARRAY< Mat<3,3> > * dxdxi)
  {
    int size = (*xi).Size();
    x->SetSize (size);

    if (dxdxi) dxdxi->SetSize (size);

    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [ mesh[(ElementIndex) elnr].hp_elnr];
	
	ARRAY< Mat<3,3> > trans, dxdxic;
	if (dxdxi)
	  {
	    trans.SetSize (size);
	    dxdxic.SetSize (size);
	  }

	ARRAY<Point<3> > coarse_xi(size);

	double lami[8];
	FlatVector vlami(8, lami);
	const Element & el = mesh[(ElementIndex) elnr];
	int el_np = el.GetNP();

	for (int mp = 0; mp < size; mp++)
	  {
	    el.GetShapeNew ((*xi)[mp], vlami);
	    
	    Point<3> pc(0,0,0);
	    for (int i = 0; i < el_np; i++)
	      for (int l = 0; l < 3; l++)
		pc(l) += hpref_el.param[i][l] * lami[i];
	
	    coarse_xi[mp] = pc;

	    if (dxdxi)
	      {
		MatrixFixWidth<3> dlami(8);
		dlami = 0;
		mesh[(ElementIndex) elnr] . GetDShapeNew ((*xi)[mp], dlami);	  
		
		trans[mp] = 0;
		for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 3; l++)
		    {
		      double sum = 0;
		      for (int i = 0; i < 8; i++)
			sum += hpref_el.param[i][l] * dlami(i, k);
		      trans[mp](l,k) = sum;
		    }
	      }
	  }

	if (dxdxi)
	  mesh.coarsemesh->GetCurvedElements().CalcMultiPointElementTransformation (&coarse_xi, hpref_el.coarse_elnr, x, &dxdxic);
	else
	  mesh.coarsemesh->GetCurvedElements().CalcMultiPointElementTransformation (&coarse_xi, hpref_el.coarse_elnr, x, 0);
	
	if (dxdxi)
	  for (int mp = 0; mp < size; mp++)
	    (*dxdxi)[mp] = dxdxic[mp] * trans[mp];
	
	return;
      }
    








    Element elem = mesh[(ElementIndex) elnr];
    BaseFiniteElement3D * fe3d;
    
    // char locmem[max2(sizeof(FETet), sizeof(FEPrism))];
    char locmemtet[sizeof(FETet)];
    char locmemprism[sizeof(FEPrism)];
    char locmempyramid[sizeof(FEPyramid)];
    char locmemhex[sizeof(FEHex)];
    switch (elem.GetType())
      {
      case TET: fe3d = new (locmemtet) FETet (*this); break;
      case PRISM: fe3d = new (locmemprism) FEPrism (*this); break;
      case PYRAMID: fe3d = new (locmempyramid) FEPyramid (*this); break;
      case HEX: fe3d = new (locmemhex) FEHex (*this); break;
      }
    
    fe3d->SetElementNumber (elnr+1);

    if (dxdxi)
      for (int mp = 0; mp < size; mp++)
	(*dxdxi)[mp] = 0.0;
    
    Vec<3> vertices[8];
    ArrayMem<Vec<3>, 100> loc_edgecoefs(0);
    ArrayMem<Vec<3>, 100> loc_facecoefs(0);

    for (int v = 0; v < fe3d->GetNVertices(); v++)
      vertices[v] = Vec<3> (mesh.Point(fe3d->GetVertexNr(v)));

    if (IsHighOrder())
      {
	for (int e = 0; e < fe3d->GetNEdges(); e++)
	  {
	    int gindex = edgecoeffsindex[fe3d->GetEdgeNr(e)-1];
	    for (int k = 2; k <= fe3d->GetEdgeOrder(e); k++, gindex++)
	      loc_edgecoefs.Append (edgecoeffs[gindex]);
	  }    
	
	if (mesh.GetDimension() == 3)
	  for (int f = 0; f < fe3d->GetNFaces(); f++)
	    {
	      int gindex = facecoeffsindex[fe3d->GetFaceNr(f)-1];
	      for (int k = 0; k < fe3d->GetNFaceShapes(f); k++, gindex++)
		loc_facecoefs.Append (facecoeffs[gindex]);
	    }    
      }



    if (!dxdxi)

      for (int mp = 0; mp < size; mp++)
	{
	  fe3d->SetReferencePoint ((*xi)[mp]);
	  fe3d->CalcVertexShapesOnly ();
	  
	  if (IsHighOrder())
	    {
	      fe3d->CalcEdgeShapesOnly ();
	      if (mesh.GetDimension() == 3)
		fe3d->CalcFaceShapes ();
	    }
	  
	  Point<3> hp (0,0,0);
	  for (int v = 0; v < fe3d->GetNVertices(); v++)
	    hp += fe3d->GetVertexShape(v) * vertices[v];
	  for (int index = 0; index < loc_edgecoefs.Size(); index++)
	    hp += fe3d->GetEdgeShape(index) * loc_edgecoefs[index]; 
	  for (int index = 0; index < loc_facecoefs.Size(); index++)
	    hp += fe3d->GetFaceShape(index) * loc_facecoefs[index]; // edgecoeffs[gindex];
	  (*x)[mp] = hp;
	}
    
    else
      
      for (int mp = 0; mp < size; mp++)
	{
	  fe3d->SetReferencePoint ((*xi)[mp]);
	  fe3d->CalcVertexShapes ();
	  
	  if (IsHighOrder())
	    {
	      fe3d->CalcEdgeShapes ();
	      if (mesh.GetDimension() == 3)
		fe3d->CalcFaceShapes ();
	    }
	  
	  Point<3> hp (0,0,0);
	  for (int v = 0; v < fe3d->GetNVertices(); v++)
	    hp += fe3d->GetVertexShape(v) * vertices[v];
	  for (int index = 0; index < loc_edgecoefs.Size(); index++)
	    hp += fe3d->GetEdgeShape(index) * loc_edgecoefs[index]; 
	  for (int index = 0; index < loc_facecoefs.Size(); index++)
	    hp += fe3d->GetFaceShape(index) * loc_facecoefs[index]; // edgecoeffs[gindex];
	  (*x)[mp] = hp;
	  

	  for (int v = 0; v < fe3d->GetNVertices(); v++)
	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 3; j++)
		(*dxdxi)[mp](i,j) += fe3d->GetVertexDShape(v)(j) * vertices[v](i);
	  
	  for (int index = 0; index < loc_edgecoefs.Size(); index++)
	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 3; j++)
		(*dxdxi)[mp](i,j) += fe3d->GetEdgeDShape(index)(j) * loc_edgecoefs[index](i);
	  
	  for (int index = 0; index < loc_facecoefs.Size(); index++)
	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 3; j++)
		(*dxdxi)[mp](i,j) += fe3d->GetFaceDShape(index)(j) * loc_facecoefs[index](i);
	}
    
    fe3d -> ~BaseFiniteElement3D();
  }



  /*
    class Init
    {
    PolynomialBasis b;
    public:
    Init ();
    };

    Init :: Init()
    {
    ofstream edge("edge.out");
    b.SetOrder(10);
    for (double x = 0; x <= 1; x += 0.01)
    {
    b.CalcFScaled(x, 1);
    edge << x;
    for (int j = 2; j <= 5; j++)
    edge << " " << b.GetF(j);
    edge << endl;
    }
    }

    Init in;
  */

} // namespace netgen


#endif
