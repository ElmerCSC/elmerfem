#include <mystdlib.h>

#include "meshing.hpp"
#ifdef CURVEDELEMS_NEW

#include "../general/autodiff.hpp"


namespace netgen
{
  
  //   bool rational = true;


  
  static void ComputeGaussRule (int n, ARRAY<double> & xi, ARRAY<double> & wi)
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
  
  

  // compute edge bubbles up to order n, x \in (-1, 1)
  static void CalcEdgeShape (int n, double x, double * shape)
  {
    double p1 = x, p2 = -1, p3 = 0;
    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	shape[j-2] = p1;
      } 
  }

  static void CalcEdgeDx (int n, double x, double * dshape)
  {
    double p1 = x, p2 = -1, p3 = 0;
    double p1dx = 1, p2dx = 0, p3dx = 0;

    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p3dx = p2dx; p2dx = p1dx;

	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	p1dx = ( (2*j-3) * (x * p2dx + p2) - (j-3) * p3dx) / j;

	dshape[j-2] = p1dx;
      }    
  }

  static void CalcEdgeShapeDx (int n, double x, double * shape, double * dshape)
  {
    double p1 = x, p2 = -1, p3 = 0;
    double p1dx = 1, p2dx = 0, p3dx = 0;

    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p3dx = p2dx; p2dx = p1dx;

	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	p1dx = ( (2*j-3) * (x * p2dx + p2) - (j-3) * p3dx) / j;

	shape[j-2] = p1;
	dshape[j-2] = p1dx;
      }    
  }

  // compute L_i(x/t) * t^i
  static void CalcScaledEdgeShape (int n, double x, double t, double * shape)
  {
    double p1 = x, p2 = -1, p3 = 0;
    for (int j=0; j<=n-2; j++)
      {
	p3=p2; p2=p1;
	p1=( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
	shape[j] = p1;
      }    
  }

  template <int DIST>
  static void CalcScaledEdgeShapeDxDt (int n, double x, double t, double * dshape)
  {
    double p1 = x, p2 = -1, p3 = 0;
    double p1dx = 1, p1dt = 0;
    double p2dx = 0, p2dt = 0;
    double p3dx = 0, p3dt = 0;
     
    for (int j=0; j<=n-2; j++)
      {
	p3=p2; p3dx=p2dx; p3dt = p2dt;
	p2=p1; p2dx=p1dx; p2dt = p1dt;

	p1   = ( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
	p1dx = ( (2*j+1) * (x * p2dx + p2) - t*t*(j-1) * p3dx) / (j+2);
	p1dt = ( (2*j+1) * x * p2dt - (j-1)* (t*t*p3dt+2*t*p3)) / (j+2);

	// shape[j] = p1;
	dshape[DIST*j  ] = p1dx;
	dshape[DIST*j+1] = p1dt;
      }    
  }


  template <class Tx, class Tres>
  static void LegendrePolynomial (int n, Tx x, Tres * values)
  {
    switch (n)
      {
      case 0:
	values[0] = 1;
	break;
      case 1:
	values[0] = 1;
	values[1] = x;
	break;

      default:

	if (n < 0) return;

	Tx p1 = 1.0, p2 = 0.0, p3;
	
	values[0] = 1.0;
	for (int j=1; j<=n; j++)
	  {
	    p3 = p2; p2 = p1;
	    p1 = ((2.0*j-1.0)*x*p2 - (j-1.0)*p3) / j;
	    values[j] = p1;
	  }
      }
  }

  template <class Tx, class Tt, class Tres>
  static void ScaledLegendrePolynomial (int n, Tx x, Tt t, Tres * values)
  {
    switch (n)
      {
      case 0:
	values[0] = 1.0;
	break;

      case 1:
	values[0] = 1.0;
	values[1] = x;
	break;

      default:

	if (n < 0) return;

	Tx p1 = 1.0, p2 = 0.0, p3;
	values[0] = 1.0;
	for (int j=1; j<=n; j++)
	  {
	    p3 = p2; p2 = p1;
	    p1 = ((2.0*j-1.0)*x*p2 - t*t*(j-1.0)*p3) / j;
	    values[j] = p1;
	  }
      }
  }



template <class S, class T>
inline void JacobiPolynomial (int n, S x, double alpha, double beta, T * values)
{
  S p1 = 1.0, p2 = 0.0, p3;

  if (n >= 0) 
    p2 = values[0] = 1.0;
  if (n >= 1) 
    p1 = values[1] = 0.5 * (2*(alpha+1)+(alpha+beta+2)*(x-1));

  for (int i  = 1; i < n; i++)
    {
      p3 = p2; p2=p1;
      p1 =
	1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
	( 
	 ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) + 
	   (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
	 * p2
	 - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * p3
	 );
      values[i+1] = p1;
    }
}




template <class S, class St, class T>
inline void ScaledJacobiPolynomial (int n, S x, St t, double alpha, double beta, T * values)
{
  /*
  S p1 = 1.0, p2 = 0.0, p3;

  if (n >= 0) values[0] = 1.0;
  */

  S p1 = 1.0, p2 = 0.0, p3;

  if (n >= 0) 
    p2 = values[0] = 1.0;
  if (n >= 1) 
    p1 = values[1] = 0.5 * (2*(alpha+1)*t+(alpha+beta+2)*(x-t));

  for (int i=1; i < n; i++)
    {
      p3 = p2; p2=p1;
      p1 =
	1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
	( 
	 ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) * t + 
	   (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
	 * p2
	 - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * t * t * p3
	 );
      values[i+1] = p1;
    }
}





  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
  template <class Tx, class Ty, class Ts>
  static void CalcTrigShape (int n, Tx x, Ty y, Ts * shape)
  { 
    if (n < 3) return;
    Tx hx[50], hy[50*50];
    // ScaledLegendrePolynomial (n-3, 2*x-1, 1-y, hx);
    ScaledJacobiPolynomial (n-3, x, 1-y, 2, 2, hx);


    // LegendrePolynomial (n-3, 2*y-1, hy);
    for (int ix = 0; ix <= n-3; ix++)
      //  JacobiPolynomial (n-3, 2*y-1, 0, 0, hy+50*ix);
      JacobiPolynomial (n-3, 2*y-1, 2*ix+5, 2, hy+50*ix);

    int ii = 0;
    Tx bub = (1+x-y)*y*(1-x-y);
    for (int iy = 0; iy <= n-3; iy++)
      for (int ix = 0; ix <= n-3-iy; ix++)
	shape[ii++] = bub * hx[ix]*hy[iy+50*ix];
  }



  static void CalcTrigShapeDxDy (int n, double x, double y, double * dshape)
  { 
    if (n < 3) return;

    AutoDiff<2> adx(x, 0);
    AutoDiff<2> ady(y, 1);
    AutoDiff<2> res[2000];
    CalcTrigShape (n, adx, ady, &res[0]);
    int ndof = (n-1)*(n-2)/2;
    for (int i = 0; i < ndof; i++)
      {
	dshape[2*i] = res[i].DValue(0);
	dshape[2*i+1] = res[i].DValue(1);
      }

    /*    
	  if (n < 3) return;
    
    int ndof = (n-1)*(n-2)/2;
    double h1[1000], h2[1000];
    double eps = 1e-4;
  
    CalcTrigShape (n, x+eps, y, h1);
    CalcTrigShape (n, x-eps, y, h2);

    for (int i = 0; i < ndof; i++)
      dshape[2*i] = (h1[i]-h2[i])/(2*eps);

    CalcTrigShape (n, x, y+eps, h1);
    CalcTrigShape (n, x, y-eps, h2);

    for (int i = 0; i < ndof; i++)
      dshape[2*i+1] = (h1[i]-h2[i])/(2*eps);
    */
  }








  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
  template <class Tx, class Ty, class Tt, class Tr>
  static void CalcScaledTrigShape (int n, Tx x, Ty y, Tt t, Tr * shape)
  {
    if (n < 3) return;

    Tx hx[50], hy[50*50];
    // ScaledLegendrePolynomial (n-3, (2*x-1), t-y, hx);
    ScaledJacobiPolynomial (n-3, x, t-y, 2, 2, hx);

    // ScaledLegendrePolynomial (n-3, (2*y-1), t, hy);
    for (int ix = 0; ix <= n-3; ix++)
      ScaledJacobiPolynomial (n-3, 2*y-1, t, 2*ix+5, 2, hy+50*ix);


    int ii = 0;
    Tx bub = (t+x-y)*y*(t-x-y);
    for (int iy = 0; iy <= n-3; iy++)
      for (int ix = 0; ix <= n-3-iy; ix++)
	shape[ii++] = bub * hx[ix]*hy[iy+50*ix];
  }


  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
  static void CalcScaledTrigShapeDxDyDt (int n, double x, double y, double t, double * dshape)
  {
    if (n < 3) return;
    AutoDiff<3> adx(x, 0);
    AutoDiff<3> ady(y, 1);
    AutoDiff<3> adt(t, 2);
    AutoDiff<3> res[2000];
    CalcScaledTrigShape (n, adx, ady, adt, &res[0]);
    double dshape1[6000];
    int ndof = (n-1)*(n-2)/2;
    for (int i = 0; i < ndof; i++)
      {
	dshape[3*i] = res[i].DValue(0);
	dshape[3*i+1] = res[i].DValue(1);
	dshape[3*i+2] = res[i].DValue(2);
      }

    /*
    if (n < 3) return;
    double hvl[1000], hvr[1000];
    int nd = (n-1)*(n-2)/2;
    
    double eps = 1e-6;

    CalcScaledTrigShape (n, x-eps, y, t, hvl);
    CalcScaledTrigShape (n, x+eps, y, t, hvr);
    for (int i = 0; i < nd; i++)
      dshape[3*i] = (hvr[i]-hvl[i])/(2*eps);

    CalcScaledTrigShape (n, x, y-eps, t, hvl);
    CalcScaledTrigShape (n, x, y+eps, t, hvr);
    for (int i = 0; i < nd; i++)
      dshape[3*i+1] = (hvr[i]-hvl[i])/(2*eps);

    CalcScaledTrigShape (n, x, y, t-eps, hvl);
    CalcScaledTrigShape (n, x, y, t+eps, hvr);
    for (int i = 0; i < nd; i++)
      dshape[3*i+2] = (hvr[i]-hvl[i])/(2*eps);
    */

    /*
    for (int i = 0; i < 3*nd; i++)
      if (fabs (dshape[i]-dshape1[i]) > 1e-8 * fabs(dshape[i]) + 1e-6)
	{
	  cerr
	  cerr << "big difference: " << dshape1[i] << " != " << dshape[i] << endl;
	  }
    */
  }

    

  

  CurvedElements :: CurvedElements (const Mesh & amesh)
    : mesh (amesh)
  {
    order = 1;
    rational = 0;
  }


  CurvedElements :: ~CurvedElements()
  {
    ;
  }


  void CurvedElements :: BuildCurvedElements(Refinement * ref, int aorder,
                                             bool arational)
  {
    order = aorder;
    ishighorder = 0;

    if (mesh.coarsemesh)
      {
	mesh.coarsemesh->GetCurvedElements().BuildCurvedElements (ref, aorder, arational);
        order = aorder;
        rational = arational;
        ishighorder = (order > 1);
	return;
      }
    
    PrintMessage (1, "Curve elements, order = ", aorder);
    if (rational) PrintMessage (1, "curved elements with rational splines");

    const_cast<Mesh&> (mesh).UpdateTopology();
    const MeshTopology & top = mesh.GetTopology();


    rational = arational;

    ARRAY<int> edgenrs;

    edgeorder.SetSize (top.GetNEdges());
    faceorder.SetSize (top.GetNFaces());

    edgeorder = 1;
    faceorder = 1;

    if (rational)
      {
        edgeweight.SetSize (top.GetNEdges());
        edgeweight = 1.0;
      }

    
    if (aorder <= 1) return;


    if (rational) aorder = 2;


    if (mesh.GetDimension() == 3)
      for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	{
	  top.GetSurfaceElementEdges (i+1, edgenrs);
	  for (int j = 0; j < edgenrs.Size(); j++)
	    edgeorder[edgenrs[j]-1] = aorder;
	  faceorder[top.GetSurfaceElementFace (i+1)-1] = aorder;
	}
    for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
      edgeorder[top.GetSegmentEdge (i+1)-1] = aorder;

    if (rational)
      {
        edgeorder = 2;
        faceorder = 1;
      }


    edgecoeffsindex.SetSize (top.GetNEdges()+1);
    int nd = 0;
    for (int i = 0; i < top.GetNEdges(); i++)
      {
	edgecoeffsindex[i] = nd;
	nd += max (0, edgeorder[i]-1);
      }
    edgecoeffsindex[top.GetNEdges()] = nd;

    edgecoeffs.SetSize (nd);
    edgecoeffs = Vec<3> (0,0,0);
    

    facecoeffsindex.SetSize (top.GetNFaces()+1);
    nd = 0;
    for (int i = 0; i < top.GetNFaces(); i++)
      {
	facecoeffsindex[i] = nd;
	if (top.GetFaceType(i+1) == TRIG)
	  nd += max (0, (faceorder[i]-1)*(faceorder[i]-2)/2);
	else
	  nd += max (0, sqr(faceorder[i]-1));
      }
    facecoeffsindex[top.GetNFaces()] = nd;

    facecoeffs.SetSize (nd);
    facecoeffs = Vec<3> (0,0,0);


    if (!ref || aorder <= 1) 
      {
        order = aorder;
        return;
      }
    
    ARRAY<double> xi, weight;

    ComputeGaussRule (aorder+4, xi, weight);  // on (0,1)

    PrintMessage (3, "Curving edges");

    if (mesh.GetDimension() == 3 || rational)
      for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
        {
	  SetThreadPercent(double(i)/mesh.GetNSE()*100.);
          const Element2d & el = mesh[i];
          top.GetSurfaceElementEdges (i+1, edgenrs);
          for (int j = 0; j < edgenrs.Size(); j++)
            edgenrs[j]--;
          const ELEMENT_EDGE * edges = MeshTopology::GetEdges (el.GetType());

          for (int i2 = 0; i2 < edgenrs.Size(); i2++)
            {
              PointIndex pi1 = edges[i2][0]-1;
              PointIndex pi2 = edges[i2][1]-1;

              bool swap = el[pi1] > el[pi2];

              Point<3> p1 = mesh[el[pi1]];
              Point<3> p2 = mesh[el[pi2]];

              int order1 = edgeorder[edgenrs[i2]];
              int ndof = max (0, order1-1);

              if (rational && order1 >= 2)
                {
                  Point<3> pm = Center (p1, p2);

                  int surfnr = mesh.GetFaceDescriptor(el.GetIndex()).SurfNr();

                  Vec<3> n1 = ref -> GetNormal (p1, surfnr, el.GeomInfoPi(edges[i2][0]));
                  Vec<3> n2 = ref -> GetNormal (p2, surfnr, el.GeomInfoPi(edges[i2][1]));

                  // p3 = pm + alpha1 n1 + alpha2 n2
                  
                  Mat<2> mat, inv;
                  Vec<2> rhs, sol;

                  mat(0,0) = n1*n1;
                  mat(0,1) = mat(1,0) = n1*n2;
                  mat(1,1) = n2*n2;
                  
                  rhs(0) = n1 * (p1-pm);
                  rhs(1) = n2 * (p2-pm);
                  

                  Point<3> p3;

                  if (fabs (Det (mat)) > 1e-10)
                    {
                      CalcInverse (mat, inv);
                      sol = inv * rhs;

                      p3 = pm + sol(0) * n1 + sol(1) * n2;
                    }
                  else
                    p3 = pm;

                  edgecoeffs[edgecoeffsindex[edgenrs[i2]]] = Vec<3> (p3);


                  double wold = 1, w = 1, dw = 0.1;
                  double dold = 1e99;
                  while (fabs (dw) > 1e-12)
                    {
                      Vec<3> v05 = 0.25 * Vec<3> (p1) + 0.5*w* Vec<3>(p3) + 0.25 * Vec<3> (p2);
                      v05 /= 1 + (w-1) * 0.5;
                      Point<3> p05 (v05), pp05(v05);
                      ref -> ProjectToSurface (pp05, surfnr,  el.GeomInfoPi(edges[i2][0]));
                      double d = Dist (pp05, p05);
                      
                      if (d < dold)
                        {
                          dold = d;
                          wold = w;
                          w += dw;
                        }
                      else
                        {
                          dw *= -0.7;
                          w = wold + dw;
                        }
                    }
                  
                  edgeweight[edgenrs[i2]] = w;
                  continue;
                }
	    
              Vector shape(ndof);
              DenseMatrix mat(ndof, ndof), inv(ndof, ndof),
                rhs(ndof, 3), sol(ndof, 3);
	    
              rhs = 0.0;
              mat = 0.0;
              for (int j = 0; j < xi.Size(); j++)
                {
                  Point<3> p;
                  Point<3> pp;
                  PointGeomInfo ppgi;
		
                  if (swap)
                    {
                      p = p1 + xi[j] * (p2-p1);
                      ref -> PointBetween (p1, p2, xi[j], 
                                           mesh.GetFaceDescriptor(el.GetIndex()).SurfNr(),
                                           el.GeomInfoPi(edges[i2][0]),
                                           el.GeomInfoPi(edges[i2][1]),
                                           pp, ppgi);
                    }
                  else
                    {
                      p = p2 + xi[j] * (p1-p2);
                      ref -> PointBetween (p2, p1, xi[j], 
                                           mesh.GetFaceDescriptor(el.GetIndex()).SurfNr(),
                                           el.GeomInfoPi(edges[i2][1]),
                                           el.GeomInfoPi(edges[i2][0]),
                                           pp, ppgi);
                    }
		
                  Vec<3> dist = pp - p;
		
                  CalcEdgeShape (order1, 2*xi[j]-1, &shape(0));
		
                  for (int k = 0; k < ndof; k++)
                    for (int l = 0; l < ndof; l++)
                      mat(k,l) += weight[j] * shape(k) * shape(l);
		
                  for (int k = 0; k < ndof; k++)
                    for (int l = 0; l < 3; l++)
                      rhs(k,l) += weight[j] * shape(k) * dist(l);
                }
	    
              CalcInverse (mat, inv);
              Mult (inv, rhs, sol);

	      
	    
              int first = edgecoeffsindex[edgenrs[i2]];
              for (int j = 0; j < ndof; j++)
                for (int k = 0; k < 3; k++)
                  edgecoeffs[first+j](k) = sol(j,k);
            }
        }


    
    for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
      {
	SetThreadPercent(double(i)/mesh.GetNSeg()*100.);
	const Segment & seg = mesh[i];
	PointIndex pi1 = mesh[i].p1;
	PointIndex pi2 = mesh[i].p2;

	bool swap = (pi1 > pi2);

	Point<3> p1 = mesh[pi1];
	Point<3> p2 = mesh[pi2];

	int segnr = top.GetSegmentEdge (i+1)-1;

	int order1 = edgeorder[segnr];
	int ndof = max (0, order1-1);


        if (rational)
          {
            Vec<3> tau1 = ref -> GetTangent (p1, seg.surfnr2, seg.surfnr1, seg.epgeominfo[0]);
            Vec<3> tau2 = ref -> GetTangent (p2, seg.surfnr2, seg.surfnr1, seg.epgeominfo[1]);
            // p1 + alpha1 tau1 = p2 + alpha2 tau2;

            Mat<3,2> mat;
            Mat<2,3> inv;
            Vec<3> rhs;
            Vec<2> sol;
            for (int j = 0; j < 3; j++)
              {
                mat(j,0) = tau1(j); 
                mat(j,1) = -tau2(j); 
                rhs(j) = p2(j)-p1(j); 
              }
            CalcInverse (mat, inv);
            sol = inv * rhs;

            Point<3> p3 = p1+sol(0) * tau1;
            edgecoeffs[edgecoeffsindex[segnr]] = Vec<3> (p3);

            
//             double dopt = 1e99;
//             double wopt = 0;
//             for (double w = 0; w <= 2; w += 0.0001)
//               {
//                 Vec<3> v05 = 0.25 * Vec<3> (p1) + 0.5*w* Vec<3>(p3) + 0.25 * Vec<3> (p2);
//                 v05 /= 1 + (w-1) * 0.5;
//                 Point<3> p05 (v05), pp05(v05);
//                 ref -> ProjectToEdge (pp05, seg.surfnr1, seg.surfnr2, seg.epgeominfo[0]);
//                 double d = Dist (pp05, p05);
//                 if (d < dopt)
//                   {
//                     wopt = w;
//                     dopt = d;
//                   }
//               }
            
            double wold = 1, w = 1, dw = 0.1;
            double dold = 1e99;
            while (fabs (dw) > 1e-12)
              {
                Vec<3> v05 = 0.25 * Vec<3> (p1) + 0.5*w* Vec<3>(p3) + 0.25 * Vec<3> (p2);
                v05 /= 1 + (w-1) * 0.5;
                Point<3> p05 (v05), pp05(v05);
                ref -> ProjectToEdge (pp05, seg.surfnr1, seg.surfnr2, seg.epgeominfo[0]);
                double d = Dist (pp05, p05);

                if (d < dold)
                  {
                    dold = d;
                    wold = w;
                    w += dw;
                  }
                else
                  {
                    dw *= -0.7;
                    w = wold + dw;
                  }
                // *testout << "w = " << w << ", dw = " << dw << endl;
              }

            // cout << "wopt = " << w << ", dopt = " << dold << endl;
            edgeweight[segnr] = w;
            
//             cout << "p1 = " << p1 << ", tau1 = " << tau1 << ", alpha1 = " << sol(0) << endl;
//             cout << "p2 = " << p2 << ", tau2 = " << tau2 << ", alpha2 = " << -sol(1) << endl;
//             cout << "p+alpha tau = " << p1 + sol(0) * tau1 
//                  << " =?= " << p2 +sol(1) * tau2 << endl;
            
          }

        else
          
          {

            Vector shape(ndof);
            DenseMatrix mat(ndof, ndof), inv(ndof, ndof),
              rhs(ndof, 3), sol(ndof, 3);

            rhs = 0.0;
            mat = 0.0;
            for (int j = 0; j < xi.Size(); j++)
              {
                Point<3> p;

                Point<3> pp;
                EdgePointGeomInfo ppgi;
	    
                if (swap)
                  {
                    p = p1 + xi[j] * (p2-p1);
                    ref -> PointBetween (p1, p2, xi[j], 
                                         seg.surfnr2, seg.surfnr1, 
                                         seg.epgeominfo[0], seg.epgeominfo[1],
                                         pp, ppgi);
                  }
                else
                  {
                    p = p2 + xi[j] * (p1-p2);
                    ref -> PointBetween (p2, p1, xi[j], 
                                         seg.surfnr2, seg.surfnr1, 
                                         seg.epgeominfo[1], seg.epgeominfo[0],
                                         pp, ppgi);
                  }
	    
                Vec<3> dist = pp - p;

                CalcEdgeShape (order1, 2*xi[j]-1, &shape(0));

                for (int k = 0; k < ndof; k++)
                  for (int l = 0; l < ndof; l++)
                    mat(k,l) += weight[j] * shape(k) * shape(l);

                for (int k = 0; k < ndof; k++)
                  for (int l = 0; l < 3; l++)
                    rhs(k,l) += weight[j] * shape(k) * dist(l);
              }


            CalcInverse (mat, inv);
            Mult (inv, rhs, sol);

            int first = edgecoeffsindex[segnr];
            for (int j = 0; j < ndof; j++)
              for (int k = 0; k < 3; k++)
                edgecoeffs[first+j](k) = sol(j,k);
          }
      }

   

    
    PrintMessage (3, "Curving faces");

    if (mesh.GetDimension() == 3)
    for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
      {
	SetThreadPercent(double(i)/mesh.GetNSE()*100.);
	const Element2d & el = mesh[i];
	int facenr = top.GetSurfaceElementFace (i+1)-1;

	if (el.GetType() == TRIG && order >= 3)
	  {
	    int fnums[] = { 0, 1, 2 };
	    if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	    if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
	    if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);

	    int order1 = faceorder[facenr];
	    int ndof = max (0, (order1-1)*(order1-2)/2);
	    
	    Vector shape(ndof);
	    DenseMatrix mat(ndof, ndof), inv(ndof, ndof),
	      rhs(ndof, 3), sol(ndof, 3);
	    
	    rhs = 0.0;
	    mat = 0.0;

	    for (int jx = 0; jx < xi.Size(); jx++)
	      for (int jy = 0; jy < xi.Size(); jy++)
		{
		  double y = xi[jy];
		  double x = (1-y) * xi[jx];
		  double lami[] = { x, y, 1-x-y };
		  double wi = weight[jx]*weight[jy]*(1-y);

		  Point<2> xii (x, y);
		  Point<3> p1, p2;
		  CalcSurfaceTransformation (xii, i, p1);
		  p2 = p1;
		  ref -> ProjectToSurface (p2, mesh.GetFaceDescriptor(el.GetIndex()).SurfNr());

		  Vec<3> dist = p2-p1;
		
		  CalcTrigShape (order1, lami[fnums[1]]-lami[fnums[0]],
				 1-lami[fnums[1]]-lami[fnums[0]], &shape(0));

		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < ndof; l++)
		      mat(k,l) += wi * shape(k) * shape(l);
		  
		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < 3; l++)
		      rhs(k,l) += wi * shape(k) * dist(l);
		}

	    // *testout << "mat = " << endl << mat << endl;
	    // CalcInverse (mat, inv);
	    // Mult (inv, rhs, sol);

	    for (int i = 0; i < ndof; i++)
	      for (int j = 0; j < 3; j++)
		sol(i,j) = rhs(i,j) / mat(i,i);   // Orthogonal basis !

	    int first = facecoeffsindex[facenr];
	    for (int j = 0; j < ndof; j++)
	      for (int k = 0; k < 3; k++)
		facecoeffs[first+j](k) = sol(j,k);
	  }
      }
    
    PrintMessage (3, "Complete");


    // compress edge and face tables
    int newbase = 0;
    for (int i = 0; i < edgeorder.Size(); i++)
      {
	bool curved = 0;
	int oldbase = edgecoeffsindex[i];
	nd = edgecoeffsindex[i+1] - edgecoeffsindex[i];

	for (int j = 0; j < nd; j++)
	  if (edgecoeffs[oldbase+j].Length() > 1e-12)
	    curved = 1;
        if (rational) curved = 1;

	if (curved && newbase != oldbase)
	  for (int j = 0; j < nd; j++)
	    edgecoeffs[newbase+j] = edgecoeffs[oldbase+j];

	edgecoeffsindex[i] = newbase;
	if (!curved) edgeorder[i] = 1;
	if (curved) newbase += nd;
      }
    edgecoeffsindex.Last() = newbase;


    newbase = 0;
    for (int i = 0; i < faceorder.Size(); i++)
      {
	bool curved = 0;
	int oldbase = facecoeffsindex[i];
	nd = facecoeffsindex[i+1] - facecoeffsindex[i];

	for (int j = 0; j < nd; j++)
	  if (facecoeffs[oldbase+j].Length() > 1e-12)
	    curved = 1;

	if (curved && newbase != oldbase)
	  for (int j = 0; j < nd; j++)
	    facecoeffs[newbase+j] = facecoeffs[oldbase+j];

	facecoeffsindex[i] = newbase;
	if (!curved) faceorder[i] = 1;
	if (curved) newbase += nd;
      }
    facecoeffsindex.Last() = newbase;

    ishighorder = (order > 1);
    // (*testout) << "edgecoeffs = " << endl << edgecoeffs << endl;
    // (*testout) << "facecoeffs = " << endl << facecoeffs << endl;
  }










  // ***********************  Transform edges *****************************

  
  bool CurvedElements ::  IsSegmentCurved (SegmentIndex elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsSegmentCurved (hpref_el.coarse_elnr);
      }

    SegmentInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = 2;
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	info.edgenr = top.GetSegmentEdge (elnr+1)-1;	
	info.ndof += edgeorder[info.edgenr]-1;
      }

    return (info.ndof > info.nv);
  }


 
  

  void CurvedElements :: 
  CalcSegmentTransformation (double xi, SegmentIndex elnr,
			     Point<3> * x, Vec<3> * dxdxi, bool * curved)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	// xi umrechnen
	double lami[2] = { xi, 1-xi };
	double dlami[2] = { 1, -1 };

	double coarse_xi = 0;
	double trans = 0;
	for (int i = 0; i < 2; i++)
	  {
	    coarse_xi += hpref_el.param[i][0] * lami[i];
	    trans += hpref_el.param[i][0] * dlami[i];
	  }

	mesh.coarsemesh->GetCurvedElements().CalcSegmentTransformation (coarse_xi, hpref_el.coarse_elnr, x, dxdxi, curved);
	if (dxdxi) *dxdxi *= trans;
	
	return;
      }
    

    Vector shapes, dshapes;
    ARRAY<Vec<3> > coefs;

    SegmentInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = 2;

    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	info.edgenr = top.GetSegmentEdge (elnr+1)-1;	
	info.ndof += edgeorder[info.edgenr]-1;
      }

    CalcElementShapes (info, xi, shapes);
    GetCoefficients (info, coefs);

    *x = 0;
    for (int i = 0; i < shapes.Size(); i++)
      *x += shapes(i) * coefs[i];


    if (dxdxi)
      {
	CalcElementDShapes (info, xi, dshapes);
	
	*dxdxi = 0;
	for (int i = 0; i < shapes.Size(); i++)
	  for (int j = 0; j < 3; j++)
	    (*dxdxi)(j) += dshapes(i) * coefs[i](j);
      }

    if (curved)
      *curved = (info.order > 1);

    // cout << "Segment, |x| = " << Abs2(Vec<3> (*x) ) << endl;
  }




  void CurvedElements :: 
  CalcElementShapes (SegmentInfo & info, double xi, Vector & shapes) const
  {
    if (rational && info.order == 2)
      {
        shapes.SetSize(3);
        double w = edgeweight[info.edgenr];
        shapes(0) = xi*xi;
        shapes(1) = (1-xi)*(1-xi);
        shapes(2) = 2*w*xi*(1-xi);
        shapes *= 1.0 / (1 + (w-1) *2*xi*(1-xi));
        return;
      }


    shapes.SetSize(info.ndof);
    shapes(0) = xi;
    shapes(1) = 1-xi;

    if (info.order >= 2)
      {
	if (mesh[info.elnr].p1 > mesh[info.elnr].p2)
	  xi = 1-xi;
	CalcEdgeShape (edgeorder[info.edgenr], 2*xi-1, &shapes(2));
      }
  }

  void CurvedElements :: 
  CalcElementDShapes (SegmentInfo & info, double xi, Vector & dshapes) const
  {
    if (rational && info.order == 2)
      {
        dshapes.SetSize(3);
        double wi = edgeweight[info.edgenr];
        double shapes[3];
        shapes[0] = xi*xi;
        shapes[1] = (1-xi)*(1-xi);
        shapes[2] = 2*wi*xi*(1-xi);
        double w = 1 + (wi-1) *2*xi*(1-xi);
        double dw = (wi-1) * (2 - 4*xi);
        
        dshapes(0) = 2*xi;
        dshapes(1) = 2*(xi-1);
        dshapes(2) = 2*wi*(1-2*xi);

        for (int j = 0;j < 3; j++)
          dshapes(j) = dshapes(j) / w - shapes[j] * dw / (w*w);
        return;
      }






    dshapes.SetSize(info.ndof);
    dshapes = 0;
    dshapes(0) = 1;
    dshapes(1) = -1;

    // int order = edgeorder[info.edgenr];

    if (info.order >= 2)
      {
        double fac = 2;
	if (mesh[info.elnr].p1 > mesh[info.elnr].p2)
          {
            xi = 1-xi; 
            fac *= -1;
          }
	CalcEdgeDx (edgeorder[info.edgenr], 2*xi-1, &dshapes(2));
        for (int i = 2; i < dshapes.Size(); i++)
          dshapes(i) *= fac;
      }

    // ??? not implemented ????
  }

  void CurvedElements :: 
  GetCoefficients (SegmentInfo & info, ARRAY<Vec<3> > & coefs) const
  {
    const Segment & el = mesh[info.elnr];

    coefs.SetSize(info.ndof);

    coefs[0] = Vec<3> (mesh[el.p1]);
    coefs[1] = Vec<3> (mesh[el.p2]);

    if (info.order >= 2)
      {
	int first = edgecoeffsindex[info.edgenr]; 
	int next = edgecoeffsindex[info.edgenr+1]; 
	for (int i = 0; i < next-first; i++)
	  coefs[i+2] = edgecoeffs[first+i];
      }
  }













  // ********************** Transform surface elements *******************


  bool CurvedElements :: IsSurfaceElementCurved (SurfaceElementIndex elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsSurfaceElementCurved (hpref_el.coarse_elnr);
      }

    const Element2d & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();
    
    SurfaceElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = (type == TRIG) ? 3 : 4;
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.edgenrs[i]--;
	info.facenr = top.GetSurfaceElementFace (elnr+1)-1;

	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    return (info.ndof > info.nv);
  }
  
  void CurvedElements :: 
  CalcSurfaceTransformation (Point<2> xi, SurfaceElementIndex elnr,
			     Point<3> * x, Mat<3,2> * dxdxi, bool * curved)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	  // xi umrechnen
	double lami[4];
	FlatVector vlami(4, lami);
	vlami = 0;
	mesh[elnr].GetShapeNew (xi, vlami);
	
	Mat<2,2> trans;
	Mat<3,2> dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<2> dlami(4);
	    dlami = 0;
	    mesh[elnr].GetDShapeNew (xi, dlami);	  
	    
	    trans = 0;
	    for (int k = 0; k < 2; k++)
	      for (int l = 0; l < 2; l++)
		for (int i = 0; i < hpref_el.np; i++)
		  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
	    }
	
	Point<2> coarse_xi(0,0);
	for (int i = 0; i < hpref_el.np; i++)
	  for (int j = 0; j < 2; j++)
	    coarse_xi(j) += hpref_el.param[i][j] * lami[i];
	
	mesh.coarsemesh->GetCurvedElements().CalcSurfaceTransformation (coarse_xi, hpref_el.coarse_elnr, x, &dxdxic, curved);
	
	if (dxdxi)
	  *dxdxi = dxdxic * trans;
	
	return;
      }
    


    Vector shapes;
    DenseMatrix dshapes;
    ARRAY<Vec<3> > coefs;

    const Element2d & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    SurfaceElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = (type == TRIG) ? 3 : 4;
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.edgenrs[i]--;
	info.facenr = top.GetSurfaceElementFace (elnr+1)-1;


	bool firsttry = true;
	bool problem = false;

	while(firsttry || problem)
	  {
	    problem = false;

	    for (int i = 0; !problem && i < info.edgenrs.Size(); i++)
	      {
		if(info.edgenrs[i]+1 >= edgecoeffsindex.Size())
		  problem = true;
		else
		  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	      }
	    if(info.facenr+1 >= facecoeffsindex.Size())
	      problem = true;
	    else
	      info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];

	    if(problem && !firsttry)
	      throw NgException("something wrong with curved elements");
	    
	    if(problem)
	      BuildCurvedElements(NULL,order,rational);

	    firsttry = false;
	  }
      }

    CalcElementShapes (info, xi, shapes);
    GetCoefficients (info, coefs);

    *x = 0;
    for (int i = 0; i < coefs.Size(); i++)
      *x += shapes(i) * coefs[i];

    if (dxdxi)
      {
	CalcElementDShapes (info, xi, dshapes);
	
	*dxdxi = 0;
	for (int i = 0; i < coefs.Size(); i++)
	  for (int j = 0; j < 3; j++)
	    for (int k = 0; k < 2; k++)
	      (*dxdxi)(j,k) += dshapes(i,k) * coefs[i](j);
      }

    if (curved)
      *curved = (info.ndof > info.nv);

    
    // cout.precision(12);
    // cout << "face, |x| = " << Abs2(Vec<3> (*x) )/sqr(M_PI) << endl;
  }




  void CurvedElements :: 
  CalcElementShapes (SurfaceElementInfo & info, const Point<2> & xi, Vector & shapes) const
  {
    const Element2d & el = mesh[info.elnr];

    shapes.SetSize(info.ndof);
    shapes = 0;	  

    if (rational && info.order >= 2)
      {
        shapes.SetSize(6);
        double w = 1;
        double lami[3] = { xi(0), xi(1), 1-xi(0)-xi(1) };
        for (int j = 0; j < 3; j++)
          shapes(j) = lami[j] * lami[j];

        const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TRIG);
        for (int j = 0; j < 3; j++)
          {
            double wi = edgeweight[info.edgenrs[j]];
            shapes(j+3) = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
            w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
          }

        shapes *= 1.0 / w;
        return;
      }

    switch (el.GetType())
      {
      case TRIG:
	{
	  shapes(0) = xi(0);
	  shapes(1) = xi(1);
	  shapes(2) = 1-xi(0)-xi(1);

	  if (info.order == 1) return;

	  int ii = 3;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TRIG);
	  
	  for (int i = 0; i < 3; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, shapes(vi1)-shapes(vi2), shapes(vi1)+shapes(vi2), &shapes(ii));
		  ii += eorder-1;
		}
	    }

	  int forder = faceorder[info.facenr];
	  if (forder >= 3)
	    {
	      int fnums[] = { 0, 1, 2 };
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      
	      CalcTrigShape (forder, 
			     shapes(fnums[1])-shapes(fnums[0]),
			     1-shapes(fnums[1])-shapes(fnums[0]), &shapes(ii));
	    }
	  break;
	}

      case QUAD:
	{
	  shapes(0) = (1-xi(0))*(1-xi(1));
	  shapes(1) =    xi(0) *(1-xi(1));
	  shapes(2) =    xi(0) *   xi(1) ;
	  shapes(3) = (1-xi(0))*   xi(1) ;

	  if (info.order == 1) return;
	  
	  double mu[4] = { 
	    1 - xi(0) + 1 - xi(1), 
	        xi(0) + 1 - xi(1), 
	        xi(0) +     xi(1), 
	    1 - xi(0) +     xi(1), 
	  };
	    
	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (QUAD);
	  
	  for (int i = 0; i < 4; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcEdgeShape (eorder, mu[vi1]-mu[vi2], &shapes(ii));
		  double lame = shapes(vi1)+shapes(vi2);
		  for (int j = 0; j < order-1; j++)
		    shapes(ii+j) *= lame;
		  ii += eorder-1;
		}
	    }
	  
	  for (int i = ii; i < info.ndof; i++)
	    shapes(i) = 0;

	  break;
	}
      };
  }


  void CurvedElements :: 
  CalcElementDShapes (SurfaceElementInfo & info, const Point<2> & xi, DenseMatrix & dshapes) const
  {
    const Element2d & el = mesh[info.elnr];
    ELEMENT_TYPE type = el.GetType();

    double lami[4];

    dshapes.SetSize(info.ndof,2);
    dshapes = 0;	  

    // *testout << "calcelementdshapes, info.ndof = " << info.ndof << endl;

    if (rational && info.order >= 2)
      {
        double w = 1;
        double dw[2] = { 0, 0 };


        lami[0] = xi(0); lami[1] = xi(1); lami[2] = 1-xi(0)-xi(1);
        double dlami[3][2] = { { 1, 0 }, { 0, 1 }, { -1, -1 }};
        double shapes[6];

        for (int j = 0; j < 3; j++)
          {
            shapes[j] = lami[j] * lami[j];
            dshapes(j,0) = 2 * lami[j] * dlami[j][0];
            dshapes(j,1) = 2 * lami[j] * dlami[j][1];
          }

        const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TRIG);
        for (int j = 0; j < 3; j++)
          {
            double wi = edgeweight[info.edgenrs[j]];

            shapes[j+3] = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
            for (int k = 0; k < 2; k++)
              dshapes(j+3,k) = 2*wi* (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
                                      lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);

            w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
            for (int k = 0; k < 2; k++)
              dw[k] += 2*(wi-1) * (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
                                  lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);
          }
        // shapes *= 1.0 / w;
        dshapes *= 1.0 / w;
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 2; j++)
            dshapes(i,j) -= shapes[i] * dw[j] / (w*w);
        return;
      }





    switch (type)
      {
      case TRIG:
	{
	  dshapes(0,0) = 1;
	  dshapes(1,1) = 1;
	  dshapes(2,0) = -1;
	  dshapes(2,1) = -1;
	  
	  if (info.order == 1) return;

          // *testout << "info.order = " << info.order << endl;


	  lami[0] = xi(0);
	  lami[1] = xi(1);
	  lami[2] = 1-xi(0)-xi(1);

	  int ii = 3;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TRIG);
	  
	  for (int i = 0; i < 3; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShapeDxDt<2> (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &dshapes(ii,0));

		  Mat<2,2> trans;
		  for (int j = 0; j < 2; j++)
		    {
		      trans(0,j) = dshapes(vi1,j)-dshapes(vi2,j);
		      trans(1,j) = dshapes(vi1,j)+dshapes(vi2,j);
		    }
		  
		  for (int j = 0; j < eorder-1; j++)
		    {
		      double ddx = dshapes(ii+j,0);
		      double ddt = dshapes(ii+j,1);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		    }

		  ii += eorder-1;
		}
	    }

	  int forder = faceorder[info.facenr];
          // *testout << "forder = " << forder << endl;
	  if (forder >= 3)
	    {
	      int fnums[] = { 0, 1, 2 };
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      
	      CalcTrigShapeDxDy (forder, 
				 lami[fnums[1]]-lami[fnums[0]],
				 1-lami[fnums[1]]-lami[fnums[0]], &dshapes(ii,0));

	      int nd = (forder-1)*(forder-2)/2;
	      Mat<2,2> trans;
	      for (int j = 0; j < 2; j++)
		{
		  trans(0,j) = dshapes(fnums[1],j)-dshapes(fnums[0],j);
		  trans(1,j) = -dshapes(fnums[1],j)-dshapes(fnums[0],j);
		}

	      for (int j = 0; j < nd; j++)
		{
		  double ddx = dshapes(ii+j,0);
		  double ddt = dshapes(ii+j,1);
		  dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		  dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		}
	    }

	  break;
	}
      case QUAD:
	{
	  dshapes(0,0) = -(1-xi(1));
	  dshapes(0,1) = -(1-xi(0));
	  dshapes(1,0) =  (1-xi(1));
	  dshapes(1,1) =    -xi(0);
	  dshapes(2,0) =     xi(1);
	  dshapes(2,1) =     xi(0);
	  dshapes(3,0) =    -xi(1);
	  dshapes(3,1) =  (1-xi(0));

	  if (info.order == 1) return;

	  double shapes[4] = {
	    (1-xi(0))*(1-xi(1)),
	       xi(0) *(1-xi(1)),
	       xi(0) *   xi(1) ,
	    (1-xi(0))*   xi(1) 
	  };

	  double mu[4] = { 
	    1 - xi(0) + 1 - xi(1), 
	        xi(0) + 1 - xi(1), 
	        xi(0) +     xi(1), 
	    1 - xi(0) +     xi(1), 
	  };

	  double dmu[4][2] = {
	    { -1, -1 },
	    { 1, -1 },
	    { 1, 1 },
	    { -1, 1 } };
	    
	  // double hshapes[20], hdshapes[20];
          ArrayMem<double, 20> hshapes(order+1), hdshapes(order+1);

	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (QUAD);
	  
	  for (int i = 0; i < 4; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcEdgeShapeDx (eorder, mu[vi1]-mu[vi2], &hshapes[0], &hdshapes[0]);

		  double lame = shapes[vi1]+shapes[vi2];
		  double dlame[2] = {
		    dshapes(vi1, 0) + dshapes(vi2, 0),
		    dshapes(vi1, 1) + dshapes(vi2, 1) };
		    
		  for (int j = 0; j < eorder-1; j++)
		    for (int k = 0; k < 2; k++)
		      dshapes(ii+j, k) = 
			lame * hdshapes[j] * (dmu[vi1][k]-dmu[vi2][k])
			+ dlame[k] * hshapes[j];

		  ii += eorder-1;
		}
	    }

	  /*	  
	  *testout << "quad, dshape = " << endl << dshapes << endl;
	  for (int i = 0; i < 2; i++)
	    {
	      Point<2> xil = xi, xir = xi;
	      Vector shapesl(dshapes.Height()), shapesr(dshapes.Height());
	      xil(i) -= 1e-6;
	      xir(i) += 1e-6;
	      CalcElementShapes (info, xil, shapesl);
	      CalcElementShapes (info, xir, shapesr);
	      
	      for (int j = 0; j < dshapes.Height(); j++)
		dshapes(j,i) = 1.0 / 2e-6 * (shapesr(j)-shapesl(j));
	    }
	  
	  *testout << "quad, num dshape = " << endl << dshapes << endl;
	  */
	  break;
	}
      };
  }



  void CurvedElements :: 
  GetCoefficients (SurfaceElementInfo & info, ARRAY<Vec<3> > & coefs) const
  {
    const Element2d & el = mesh[info.elnr];

    coefs.SetSize (info.ndof);
    // coefs = Vec<3> (0,0,0);
    
    for (int i = 0; i < info.nv; i++)
      coefs[i] = Vec<3> (mesh[el[i]]);
    
    if (info.order == 1) return;

    int ii = info.nv;
	  
    for (int i = 0; i < info.edgenrs.Size(); i++)
      {
	int first = edgecoeffsindex[info.edgenrs[i]];
	int next = edgecoeffsindex[info.edgenrs[i]+1];
	for (int j = first; j < next; j++, ii++)
	  coefs[ii] = edgecoeffs[j];
      }
    
    int first = facecoeffsindex[info.facenr];
    int next = facecoeffsindex[info.facenr+1];
    for (int j = first; j < next; j++, ii++)
      coefs[ii] = facecoeffs[j];
  }







  // ********************** Transform volume elements *******************


  bool CurvedElements :: IsElementCurved (ElementIndex elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsElementCurved (hpref_el.coarse_elnr);
      }

    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();
    
    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNVertices (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	  info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	  info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	for (int i = 0; i < info.nfaces; i++)
	  info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
      }

    return (info.ndof > info.nv);
  }






  void CurvedElements :: 
  CalcElementTransformation (Point<3> xi, ElementIndex elnr,
			     Point<3> * x, Mat<3,3> * dxdxi, //  bool * curved,
                             void * buffer, bool valid)
  {
    if (mesh.coarsemesh)
      {
	  const HPRefElement & hpref_el =
	    (*mesh.hpelements) [mesh[elnr].hp_elnr];
	  
	  // xi umrechnen
	  double lami[8];
	  FlatVector vlami(8, lami);
	  vlami = 0;
	  mesh[elnr].GetShapeNew (xi, vlami);

	  Mat<3,3> trans, dxdxic;
	  if (dxdxi)
	    {
	      MatrixFixWidth<3> dlami(8);
	      dlami = 0;
	      mesh[elnr].GetDShapeNew (xi, dlami);	  
	      
	      trans = 0;
	      for (int k = 0; k < 3; k++)
		for (int l = 0; l < 3; l++)
		  for (int i = 0; i < hpref_el.np; i++)
		    trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
	    }

	  Point<3> coarse_xi(0,0,0);
	  for (int i = 0; i < hpref_el.np; i++)
	    for (int j = 0; j < 3; j++)
	      coarse_xi(j) += hpref_el.param[i][j] * lami[i];

	  mesh.coarsemesh->GetCurvedElements().CalcElementTransformation (coarse_xi, hpref_el.coarse_elnr, x, &dxdxic /* , curved */);

	  if (dxdxi)
	    *dxdxi = dxdxic * trans;

	  return;
	}


    Vector shapes;
    MatrixFixWidth<3> dshapes;

    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    ElementInfo hinfo;
    ElementInfo & info = (buffer) ? *static_cast<ElementInfo*> (buffer) : hinfo;
    

    if (!valid)
      {
        info.elnr = elnr;
        info.order = order;
        info.ndof = info.nv = MeshTopology::GetNVertices (type);
        if (info.order > 1)
          {
            const MeshTopology & top = mesh.GetTopology();
            
            info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
            for (int i = 0; i < info.nedges; i++)
              info.edgenrs[i]--;
            
            info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
            for (int i = 0; i < info.nfaces; i++)
              info.facenrs[i]--;
            
            for (int i = 0; i < info.nedges; i++)
              info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
            for (int i = 0; i < info.nfaces; i++)
              info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
          }
      }

    CalcElementShapes (info, xi, shapes);

    Vec<3> * coefs =  (info.ndof <= 10) ? 
      &info.hcoefs[0] : new Vec<3> [info.ndof];

    if (info.ndof > 10 || !valid)
      GetCoefficients (info, coefs);

    if (x)
      {
        *x = 0;
        for (int i = 0; i < shapes.Size(); i++)
          *x += shapes(i) * coefs[i];
      }

    if (dxdxi)
      {
        if (valid && info.order == 1 && info.nv == 4)   // a linear tet
          {
            *dxdxi = info.hdxdxi;
          }
        else
          {
            CalcElementDShapes (info, xi, dshapes);
            
            *dxdxi = 0;
            for (int i = 0; i < shapes.Size(); i++)
              for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                  (*dxdxi)(j,k) += dshapes(i,k) * coefs[i](j);
            
            info.hdxdxi = *dxdxi;
          }
      }

    // *testout << "curved_elements, dshapes = " << endl << dshapes << endl;

    //    if (curved) *curved = (info.ndof > info.nv);

    if (info.ndof > 10) delete [] coefs;
  }




  void CurvedElements ::   CalcElementShapes (ElementInfo & info, const Point<3> & xi, Vector & shapes) const
  {
    const Element & el = mesh[info.elnr];

    if (rational && info.order >= 2)
      {
        shapes.SetSize(10);
        double w = 1;
        double lami[4] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
        for (int j = 0; j < 4; j++)
          shapes(j) = lami[j] * lami[j];

        const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TET);
        for (int j = 0; j < 6; j++)
          {
            double wi = edgeweight[info.edgenrs[j]];
            shapes(j+4) = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
            w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
          }

        shapes *= 1.0 / w;
        return;
      }

    shapes.SetSize(info.ndof);
    
    switch (el.GetType())
      {
      case TET:
	{
	  shapes(0) = xi(0);
	  shapes(1) = xi(1);
	  shapes(2) = xi(2);
	  shapes(3) = 1-xi(0)-xi(1)-xi(2);

	  if (info.order == 1) return;

	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TET);
	  for (int i = 0; i < 6; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, shapes(vi1)-shapes(vi2), shapes(vi1)+shapes(vi2), &shapes(ii));
		  ii += eorder-1;
		}
	    }
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces (TET);
	  for (int i = 0; i < 4; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if (forder >= 3)
		{
		  int fnums[] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 }; 
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		  if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);

		  CalcScaledTrigShape (forder, 
				       shapes(fnums[1])-shapes(fnums[0]), shapes(fnums[2]), 
				       shapes(fnums[0])+shapes(fnums[1])+shapes(fnums[2]), &shapes(ii));
		  ii += (forder-1)*(forder-2)/2;
		  // CalcScaledEdgeShape (forder, shapes(vi1)-shapes(vi2), shapes(vi1)+shapes(vi2), &shapes(ii));
		  // ii += forder-1;
		}
	    }


	  break;
	}
      case PRISM:
	{
	  double lami[6] = { xi(0), xi(1), 1-xi(0)-xi(1), xi(0), xi(1), 1-xi(0)-xi(1) };
	  double lamiz[6] = { 1-xi(2), 1-xi(2), 1-xi(2), xi(2), xi(2), xi(2) };
	  for (int i = 0; i < 6; i++)
	    shapes(i) = lami[i%3] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	  for (int i = 6; i < info.ndof; i++)
	    shapes(i) = 0;

          if (info.order == 1) return;


	  int ii = 6;
 	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (PRISM);
	  for (int i = 0; i < 6; i++)    // horizontal edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1) % 3, vi2 = (edges[i][1]-1) % 3;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &shapes(ii));
		  double facz = (i < 3) ? (1-xi(2)) : xi(2);
		  for (int j = 0; j < eorder-1; j++)
		    shapes(ii+j) *= facz;

		  ii += eorder-1;
		}
	    }

	  for (int i = 6; i < 9; i++)    // vertical edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  double bubz = lamiz[vi1]*lamiz[vi2];
		  double polyz = lamiz[vi1] - lamiz[vi2];
		  double bubxy = lami[(vi1)%3];

		  for (int j = 0; j < eorder-1; j++)
		    {
		      shapes(ii+j) = bubxy * bubz;
		      bubz *= polyz;
		    }
		  ii += eorder-1;
		}
	    }

	  // FACE SHAPES
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces (PRISM);
	  for (int i = 0; i < 2; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if ( forder < 3 ) continue;
	      int fav[3] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 };
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 
	      if(el[fav[1]] > el[fav[2]]) swap(fav[1],fav[2]);
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 	

	      CalcTrigShape (forder, 
			     lami[fav[2]]-lami[fav[1]], lami[fav[0]],
			     &shapes(ii));
	      
	      int ndf = (forder+1)*(forder+2)/2 - 3 - 3*(forder-1);
	      for ( int j = 0; j < ndf; j++ )
		shapes(ii+j) *= lamiz[fav[1]];
	      ii += ndf;
	    }
	  break;
	}

      case PYRAMID:
	{
	  shapes = 0.0;
	  double x = xi(0);
	  double y = xi(1);
	  double z = xi(2);
	  
	  if (z == 1.) z = 1-1e-10;
	  shapes[0] = (1-z-x)*(1-z-y) / (1-z);
	  shapes[1] = x*(1-z-y) / (1-z);
	  shapes[2] = x*y / (1-z);
	  shapes[3] = (1-z-x)*y / (1-z);
	  shapes[4] = z;
          
          if (info.order == 1) return;

	  int ii = 5;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (PYRAMID);
	  for (int i = 0; i < 4; i++)    // horizontal edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, shapes[vi1]-shapes[vi2], shapes[vi1]+shapes[vi2], &shapes(ii));
		  double fac = (shapes[vi1]+shapes[vi2]) / (1-z);
		  for (int j = 0; j < eorder-1; j++)
		    shapes(ii+j) *= fac;

		  ii += eorder-1;
		}
	    }



	  break;
	}

      case HEX:
        {
	  shapes = 0.0;
	  double x = xi(0);
	  double y = xi(1);
	  double z = xi(2);
	  
	  shapes[0] = (1-x)*(1-y)*(1-z);
	  shapes[1] =    x *(1-y)*(1-z);
	  shapes[2] =    x *   y *(1-z);
	  shapes[3] = (1-x)*   y *(1-z);
	  shapes[4] = (1-x)*(1-y)*(z);
	  shapes[5] =    x *(1-y)*(z);
	  shapes[6] =    x *   y *(z);
	  shapes[7] = (1-x)*   y *(z);
          break;
        }
      };
  }


  void CurvedElements :: 
  CalcElementDShapes (ElementInfo & info, const Point<3> & xi, MatrixFixWidth<3> & dshapes) const
  {
    const Element & el = mesh[info.elnr];

    dshapes.SetSize(info.ndof);
    dshapes = 0.0;



    if (rational && info.order >= 2)
      {
        double w = 1;
        double dw[3] = { 0, 0, 0 };

        double lami[4] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
        double dlami[4][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { -1, -1, -1 }};
        double shapes[10];

        for (int j = 0; j < 4; j++)
          {
            shapes[j] = lami[j] * lami[j];
            dshapes(j,0) = 2 * lami[j] * dlami[j][0];
            dshapes(j,1) = 2 * lami[j] * dlami[j][1];
            dshapes(j,2) = 2 * lami[j] * dlami[j][2];
          }

        const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TET);
        for (int j = 0; j < 6; j++)
          {
            double wi = edgeweight[info.edgenrs[j]];

            shapes[j+4] = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
            for (int k = 0; k < 3; k++)
              dshapes(j+4,k) = 2*wi* (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
                                      lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);

            w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
            for (int k = 0; k < 3; k++)
              dw[k] += 2*(wi-1) * (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
                                   lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);
          }
        // shapes *= 1.0 / w;
        dshapes *= 1.0 / w;
        for (int i = 0; i < 10; i++)
          for (int j = 0; j < 3; j++)
            dshapes(i,j) -= shapes[i] * dw[j] / (w*w);
        return;
      }

    switch (el.GetType())
      {
      case TET:
	{
	  dshapes(0,0) = 1;
	  dshapes(1,1) = 1;
	  dshapes(2,2) = 1;
	  dshapes(3,0) = -1;
	  dshapes(3,1) = -1;
	  dshapes(3,2) = -1;

	  if (info.order == 1) return;

	  double lami[] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (TET);
	  for (int i = 0; i < 6; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShapeDxDt<3> (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &dshapes(ii,0));

		  Mat<2,3> trans;
		  for (int j = 0; j < 3; j++)
		    {
		      trans(0,j) = dshapes(vi1,j)-dshapes(vi2,j);
		      trans(1,j) = dshapes(vi1,j)+dshapes(vi2,j);
		    }
		  
		  for (int j = 0; j < order-1; j++)
		    {
		      double ddx = dshapes(ii+j,0);
		      double ddt = dshapes(ii+j,1);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		      dshapes(ii+j,2) = ddx * trans(0,2) + ddt * trans(1,2);
		    }

		  ii += eorder-1;
		}
	    }

	  const ELEMENT_FACE * faces = MeshTopology::GetFaces (TET);
	  for (int i = 0; i < 4; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if (forder >= 3)
		{
		  int fnums[] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 }; 
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		  if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);

		  CalcScaledTrigShapeDxDyDt (forder, 
					     lami[fnums[1]]-lami[fnums[0]], 
					     lami[fnums[2]], lami[fnums[0]]+lami[fnums[1]]+lami[fnums[2]],
					     &dshapes(ii,0));

		  Mat<3,3> trans;
		  for (int j = 0; j < 3; j++)
		    {
		      trans(0,j) = dshapes(fnums[1],j)-dshapes(fnums[0],j);
		      trans(1,j) = dshapes(fnums[2],j);
		      trans(2,j) = dshapes(fnums[0],j)+dshapes(fnums[1],j)+dshapes(fnums[2],j);
		    }
		  
		  int nfd = (forder-1)*(forder-2)/2;
		  for (int j = 0; j < nfd; j++)
		    {
		      double ddx = dshapes(ii+j,0);
		      double ddy = dshapes(ii+j,1);
		      double ddt = dshapes(ii+j,2);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddy * trans(1,0) + ddt * trans(2,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddy * trans(1,1) + ddt * trans(2,1);
		      dshapes(ii+j,2) = ddx * trans(0,2) + ddy * trans(1,2) + ddt * trans(2,2);
		    }

		  ii += nfd;
		}
	    }

	  break;
	}
      case PRISM:
	{
	  double lami[6] = { xi(0), xi(1), 1-xi(0)-xi(1), xi(0), xi(1), 1-xi(0)-xi(1)  };
	  double lamiz[6] = { 1-xi(2), 1-xi(2), 1-xi(2), xi(2), xi(2), xi(2) };
	  double dlamiz[6] = { -1, -1, -1, 1, 1, 1 };
	  double dlami[6][2] = 
	    { { 1, 0, },
	      { 0, 1, },
	      { -1, -1 },
	      { 1, 0, },
	      { 0, 1, },
	      { -1, -1 } };
	  for (int i = 0; i < 6; i++)
	    {
	      // shapes(i) = lami[i%3] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	      dshapes(i,0) = dlami[i%3][0] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	      dshapes(i,1) = dlami[i%3][1] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	      dshapes(i,2) = lami[i%3] * ( (i < 3) ? -1 : 1 );
	    }

	  int ii = 6;

	  if (info.order == 1) return;


	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges (PRISM);
	  for (int i = 0; i < 6; i++)    // horizontal edges
	    {
	      int order = edgeorder[info.edgenrs[i]];
	      if (order >= 2)
		{
		  int vi1 = (edges[i][0]-1) % 3, vi2 = (edges[i][1]-1) % 3;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  Vector shapei(order+1);
		  CalcScaledEdgeShapeDxDt<3> (order, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &dshapes(ii,0) );
		  CalcScaledEdgeShape(order, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &shapei(0) );

		  Mat<2,2> trans;
		  for (int j = 0; j < 2; j++)
		    {
		      trans(0,j) = dlami[vi1][j]-dlami[vi2][j];
		      trans(1,j) = dlami[vi1][j]+dlami[vi2][j];
		    }
		  
		  for (int j = 0; j < order-1; j++)
		    {
		      double ddx = dshapes(ii+j,0);
		      double ddt = dshapes(ii+j,1);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		    }



		  double facz = (i < 3) ? (1-xi(2)) : xi(2);
		  double dfacz = (i < 3) ? (-1) : 1;
		  for (int j = 0; j < order-1; j++)
		    {
		      dshapes(ii+j,0) *= facz;
		      dshapes(ii+j,1) *= facz;
		      dshapes(ii+j,2) = shapei(j) * dfacz;
		    }

		  ii += order-1;
		}
	    }
	  for (int i = 6; i < 9; i++)    // vertical edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  double bubz = lamiz[vi1] * lamiz[vi2];
		  double dbubz = dlamiz[vi1]*lamiz[vi2] + lamiz[vi1]*dlamiz[vi2];
		  double polyz = lamiz[vi1] - lamiz[vi2];
		  double dpolyz = dlamiz[vi1] - dlamiz[vi2];
		  double bubxy = lami[(vi1)%3];
		  double dbubxydx = dlami[(vi1)%3][0];
		  double dbubxydy = dlami[(vi1)%3][1];

		  for (int j = 0; j < eorder-1; j++)
		    {
		      dshapes(ii+j,0) = dbubxydx * bubz;
		      dshapes(ii+j,1) = dbubxydy * bubz;
		      dshapes(ii+j,2) = bubxy * dbubz;

		      dbubz = bubz * dpolyz + dbubz * polyz;
		      bubz *= polyz;
		    }
		  ii += eorder-1;
		}
	    }


	  if (info.order == 2) return;
	  // FACE SHAPES
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces (PRISM);
	  for (int i = 0; i < 2; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if ( forder < 3 ) continue;
	      int ndf = (forder+1)*(forder+2)/2 - 3 - 3*(forder-1);

	      int fav[3] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 };
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 
	      if(el[fav[1]] > el[fav[2]]) swap(fav[1],fav[2]);
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 	

	      MatrixFixWidth<2> dshapei(ndf);
	      Vector shapei(ndf);

	      CalcTrigShapeDxDy (forder, 
				 lami[fav[2]]-lami[fav[1]], lami[fav[0]],
				 &dshapei(0,0));
	      CalcTrigShape (forder, lami[fav[2]]-lami[fav[1]], lami[fav[0]],
				 &shapei(0));
	      
	      Mat<2,2> trans;
	      for (int j = 0; j < 2; j++)
		{
		  trans(0,j) = dlami[fav[2]][j]-dlami[fav[1]][j];
		  trans(1,j) = dlami[fav[0]][j];
		}
		  
	      for (int j = 0; j < order-1; j++)
		{
		  double ddx = dshapes(ii+j,0);
		  double ddt = dshapes(ii+j,1);
		  dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		  dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		}

	      for ( int j = 0; j < ndf; j++ )
		{
		  dshapes(ii+j,0) *= lamiz[fav[1]];
		  dshapes(ii+j,1) *= lamiz[fav[1]];
		  dshapes(ii+j,2) = shapei(j) * dlamiz[fav[1]];
		}
	      ii += ndf;
	    }

	  break;

	}

      case PYRAMID:
	{
	  dshapes = 0.0;
	  double x = xi(0);
	  double y = xi(1);
	  double z = xi(2);
	  
	  if (z == 1.) z = 1-1e-10;
	  double z1 = 1-z;
	  double z2 = z1*z1;
	  
	  dshapes(0,0) = -(z1-y)/z1;
	  dshapes(0,1) = -(z1-x)/z1;
	  dshapes(0,2) = ((x+y+2*z-2)*z1+(z1-y)*(z1-x))/z2;

	  dshapes(1,0) = (z1-y)/z1;
	  dshapes(1,1) = -x/z1;
	  dshapes(1,2) = (-x*z1+x*(z1-y))/z2;

	  dshapes(2,0) = y/z1;
	  dshapes(2,1) = x/z1;
	  dshapes(2,2) = x*y/z2;

	  dshapes(3,0) = -y/z1;
	  dshapes(3,1) = (z1-x)/z1;
	  dshapes(3,2) = (-y*z1+y*(z1-x))/z2;

	  dshapes(4,0) = 0;
	  dshapes(4,1) = 0;
	  dshapes(4,2) = 1;
		  /* old:
	  vdshape[0] = Vec<3>( -(z1-y)/z1, -(z1-x)/z1, ((x+y+2*z-2)*z1+(z1-y)*(z1-x))/z2 );
	  vdshape[1] = Vec<3>( (z1-y)/z1,  -x/z1,      (-x*z1+x*(z1-y))/z2 );
	  vdshape[2] = Vec<3>( y/z1,       x/z1,       x*y/z2 );
	  vdshape[3] = Vec<3>( -y/z1,      (z1-x)/z1,  (-y*z1+y*(z1-x))/z2 );
	  vdshape[4] = Vec<3>( 0, 0, 1 );
		  */
	  break;
	}

      case HEX:
        {
          dshapes = 0.0;

	  double x = xi(0);
	  double y = xi(1);
	  double z = xi(2);

	  // shapes[0] = (1-x)*(1-y)*(1-z);
          dshapes(0,0) = - (1-y)*(1-z);
          dshapes(0,1) = (1-x) * (-1) * (1-z);
          dshapes(0,2) = (1-x) * (1-y) * (-1);

	  // shapes[1] =    x *(1-y)*(1-z);
          dshapes(1,0) = (1-y)*(1-z);
          dshapes(1,1) = -x * (1-z);
          dshapes(1,2) = -x * (1-y);

	  // shapes[2] =    x *   y *(1-z);
          dshapes(2,0) = y * (1-z);
          dshapes(2,1) = x * (1-z);
          dshapes(2,2) = -x * y;

	  // shapes[3] = (1-x)*   y *(1-z);
          dshapes(3,0) = -y * (1-z);
          dshapes(3,1) = (1-x) * (1-z);
          dshapes(3,2) = -(1-x) * y;

	  // shapes[4] = (1-x)*(1-y)*z;
          dshapes(4,0) = - (1-y)*z;
          dshapes(4,1) = (1-x) * (-1) * z;
          dshapes(4,2) = (1-x) * (1-y) * 1;

	  // shapes[5] =    x *(1-y)*z;
          dshapes(5,0) = (1-y)*z;
          dshapes(5,1) = -x * z;
          dshapes(5,2) = x * (1-y);

	  // shapes[6] =    x *   y *z;
          dshapes(6,0) = y * z;
          dshapes(6,1) = x * z;
          dshapes(6,2) = x * y;

	  // shapes[7] = (1-x)*   y *z;
          dshapes(7,0) = -y * z;
          dshapes(7,1) = (1-x) * z;
          dshapes(7,2) = (1-x) * y;

          break;
        }
      }
    
    /*
    DenseMatrix dshapes2 (info.ndof, 3);
    Vector shapesl(info.ndof); 
    Vector shapesr(info.ndof); 
    
    double eps = 1e-6;
    for (int i = 0; i < 3; i++)
      {
	Point<3> xl = xi;
	Point<3> xr = xi;
	
	xl(i) -= eps;
	xr(i) += eps;
	CalcElementShapes (info, xl, shapesl);
	CalcElementShapes (info, xr, shapesr);
	
	for (int j = 0; j < info.ndof; j++)
	  dshapes2(j,i) = (shapesr(j)-shapesl(j)) / (2*eps);
      }
    (*testout) << "dshapes = " << endl << dshapes << endl;
    (*testout) << "dshapes2 = " << endl << dshapes2 << endl;
    dshapes2 -= dshapes;
    (*testout) << "diff = " << endl << dshapes2 << endl;
    */
  }



  void CurvedElements :: 
  GetCoefficients (ElementInfo & info, Vec<3> * coefs) const
  {
    // cout << "getcoeffs, info.elnr = " << info.elnr << endl;
    // cout << "getcoeffs, info.nv = " << info.nv << endl;

    const Element & el = mesh[info.elnr];

    /*
    coefs.SetSize (info.ndof);
    coefs = Vec<3> (0,0,0);
    */
    /*
    for (int i = 0; i < info.ndof; i++)
      coefs[i] = Vec<3> (0,0,0);
    */
    for (int i = 0; i < info.nv; i++)
      coefs[i] = Vec<3> (mesh[el[i]]);

    if (info.order == 1) return;

    int ii = info.nv;
	  
    for (int i = 0; i < info.nedges; i++)
      {
	int first = edgecoeffsindex[info.edgenrs[i]];
	int next = edgecoeffsindex[info.edgenrs[i]+1];
	for (int j = first; j < next; j++, ii++)
	  coefs[ii] = edgecoeffs[j];
      }
    for (int i = 0; i < info.nfaces; i++)
      {
	int first = facecoeffsindex[info.facenrs[i]];
	int next = facecoeffsindex[info.facenrs[i]+1];
	for (int j = first; j < next; j++, ii++)
	  coefs[ii] = facecoeffs[j];
      }
  }































  void CurvedElements :: 
  CalcMultiPointSegmentTransformation (ARRAY<double> * xi, SegmentIndex segnr,
				       ARRAY<Point<3> > * x,
				       ARRAY<Vec<3> > * dxdxi)
  {
    ;
  }

  void CurvedElements :: 
  CalcMultiPointSurfaceTransformation (ARRAY< Point<2> > * xi, SurfaceElementIndex elnr,
				       ARRAY< Point<3> > * x,
				       ARRAY< Mat<3,2> > * dxdxi)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	  // xi umrechnen
	double lami[4];
	FlatVector vlami(4, lami);

	ArrayMem<Point<2>, 50> coarse_xi (xi->Size());
	
	for (int pi = 0; pi < xi->Size(); pi++)
	  {
	    vlami = 0;
	    mesh[elnr].GetShapeNew ( (*xi)[pi], vlami);
	    
	    Point<2> cxi(0,0);
	    for (int i = 0; i < hpref_el.np; i++)
	      for (int j = 0; j < 2; j++)
		cxi(j) += hpref_el.param[i][j] * lami[i];

	    coarse_xi[pi] = cxi;
	  }

	mesh.coarsemesh->GetCurvedElements().
	  CalcMultiPointSurfaceTransformation (&coarse_xi, hpref_el.coarse_elnr, x, dxdxi);


	Mat<2,2> trans;
        Mat<3,2> dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<2> dlami(4);
	    dlami = 0;

	    for (int pi = 0; pi < xi->Size(); pi++)
	      {
		mesh[elnr].GetDShapeNew ( (*xi)[pi], dlami);	  
		
		trans = 0;
		for (int k = 0; k < 2; k++)
		  for (int l = 0; l < 2; l++)
		    for (int i = 0; i < hpref_el.np; i++)
		      trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
		
		dxdxic = (*dxdxi)[pi];
		(*dxdxi)[pi] = dxdxic * trans;
	      }
	  }	

	return;
      }





    Vector shapes;
    DenseMatrix dshapes;
    ARRAY<Vec<3> > coefs;


    const Element2d & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    SurfaceElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = (type == TRIG) ? 3 : 4;
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.edgenrs[i]--;
	info.facenr = top.GetSurfaceElementFace (elnr+1)-1;

	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    GetCoefficients (info, coefs);

    if (x)
      {
	for (int j = 0; j < xi->Size(); j++)
	  {
	    CalcElementShapes (info, (*xi)[j], shapes);
	    (*x)[j] = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      (*x)[j] += shapes(i) * coefs[i];
	  }
      }

    if (dxdxi)
      {
	for (int ip = 0; ip < xi->Size(); ip++)
	  {
	    CalcElementDShapes (info, (*xi)[ip], dshapes);
	
	    (*dxdxi)[ip] = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 2; k++)
		  (*dxdxi)[ip](j,k) += dshapes(i,k) * coefs[i](j);
	  }
      }
  }

  void CurvedElements :: 
  CalcMultiPointElementTransformation (ARRAY< Point<3> > * xi, ElementIndex elnr,
				       ARRAY< Point<3> > * x,
				       ARRAY< Mat<3,3> > * dxdxi)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	  // xi umrechnen
	double lami[8];
	FlatVector vlami(8, lami);


	ArrayMem<Point<3>, 50> coarse_xi (xi->Size());
	
	for (int pi = 0; pi < xi->Size(); pi++)
	  {
	    vlami = 0;
	    mesh[elnr].GetShapeNew ( (*xi)[pi], vlami);
	    
	    Point<3> cxi(0,0,0);
	    for (int i = 0; i < hpref_el.np; i++)
	      for (int j = 0; j < 3; j++)
		cxi(j) += hpref_el.param[i][j] * lami[i];

	    coarse_xi[pi] = cxi;
	  }

	mesh.coarsemesh->GetCurvedElements().
	  CalcMultiPointElementTransformation (&coarse_xi, hpref_el.coarse_elnr, x, dxdxi);


	Mat<3,3> trans, dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<3> dlami(8);
	    dlami = 0;

	    for (int pi = 0; pi < xi->Size(); pi++)
	      {
		mesh[elnr].GetDShapeNew ( (*xi)[pi], dlami);	  
		
		trans = 0;
		for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 3; l++)
		    for (int i = 0; i < hpref_el.np; i++)
		      trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
		
		dxdxic = (*dxdxi)[pi];
		(*dxdxi)[pi] = dxdxic * trans;
	      }
	  }	

	return;
      }








    Vector shapes;
    MatrixFixWidth<3> dshapes;


    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNVertices (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	  info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	  info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	for (int i = 0; i < info.nfaces; i++)
	  info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
	// info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    ARRAY<Vec<3> > coefs(info.ndof);
    GetCoefficients (info, &coefs[0]);
    if (x)
      {
	for (int j = 0; j < xi->Size(); j++)
	  {
	    CalcElementShapes (info, (*xi)[j], shapes);
	    (*x)[j] = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      (*x)[j] += shapes(i) * coefs[i];
	  }
      }
    if (dxdxi)
      {
	for (int ip = 0; ip < xi->Size(); ip++)
	  {
	    CalcElementDShapes (info, (*xi)[ip], dshapes);
	
	    (*dxdxi)[ip] = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		  (*dxdxi)[ip](j,k) += dshapes(i,k) * coefs[i](j);
	  }
      }
  }




  void  CurvedElements :: 
  CalcMultiPointElementTransformation (ElementIndex elnr, int n,
                                       const double * xi, int sxi,
                                       double * x, int sx,
                                       double * dxdxi, int sdxdxi)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	  // xi umrechnen
	double lami[8];
	FlatVector vlami(8, lami);


	ArrayMem<double, 100> coarse_xi (3*n);
	
	for (int pi = 0; pi < n; pi++)
	  {
	    vlami = 0;
            Point<3> pxi;
            for (int j = 0; j < 3; j++)
              pxi(j) = xi[pi*sxi+j];

	    mesh[elnr].GetShapeNew ( pxi, vlami);
	    
	    Point<3> cxi(0,0,0);
	    for (int i = 0; i < hpref_el.np; i++)
	      for (int j = 0; j < 3; j++)
		cxi(j) += hpref_el.param[i][j] * lami[i];

            for (int j = 0; j < 3; j++)
              coarse_xi[3*pi+j] = cxi(j);
	  }

	mesh.coarsemesh->GetCurvedElements().
	  CalcMultiPointElementTransformation (hpref_el.coarse_elnr, n, 
                                               &coarse_xi[0], 3, 
                                               x, sx, 
                                               dxdxi, sdxdxi);

	Mat<3,3> trans, dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<3> dlami(8);
	    dlami = 0;

	    for (int pi = 0; pi < n; pi++)
	      {
                Point<3> pxi;
                for (int j = 0; j < 3; j++)
                  pxi(j) = xi[pi*sxi+j];

		mesh[elnr].GetDShapeNew (pxi, dlami);	  
		
		trans = 0;
		for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 3; l++)
		    for (int i = 0; i < hpref_el.np; i++)
		      trans(l,k) += hpref_el.param[i][l] * dlami(i, k);

                Mat<3> mat_dxdxic, mat_dxdxi;
                for (int j = 0; j < 3; j++)
                  for (int k = 0; k < 3; k++)
                    mat_dxdxic(j,k) = dxdxi[pi*sdxdxi+3*j+k];
		
                mat_dxdxi = mat_dxdxic * trans;

                for (int j = 0; j < 3; j++)
                  for (int k = 0; k < 3; k++)
                    dxdxi[pi*sdxdxi+3*j+k] = mat_dxdxi(j,k);

		// dxdxic = (*dxdxi)[pi];
		// (*dxdxi)[pi] = dxdxic * trans;
	      }
	  }	
	return;
      }








    Vector shapes;
    MatrixFixWidth<3> dshapes;


    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNVertices (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	  info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	  info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	for (int i = 0; i < info.nfaces; i++)
	  info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
	// info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    ARRAY<Vec<3> > coefs(info.ndof);
    GetCoefficients (info, &coefs[0]);
    if (x)
      {
	for (int j = 0; j < n; j++)
	  {
            Point<3> xij, xj;
            for (int k = 0; k < 3; k++)
              xij(k) = xi[j*sxi+k];

	    CalcElementShapes (info, xij, shapes);
	    xj = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      xj += shapes(i) * coefs[i];

            for (int k = 0; k < 3; k++)
              x[j*sx+k] = xj(k);
	  }
      }
    if (dxdxi)
      {
	for (int ip = 0; ip < n; ip++)
	  {
            Point<3> xij;
            for (int k = 0; k < 3; k++)
              xij(k) = xi[ip*sxi+k];

	    CalcElementDShapes (info, xij, dshapes);
            
            Mat<3> dxdxij;
            dxdxij = 0.0;
	    for (int i = 0; i < coefs.Size(); i++)
	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		  dxdxij(j,k) += dshapes(i,k) * coefs[i](j);


	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
                  dxdxi[ip*sdxdxi+3*j+k] = dxdxij(j,k);
	  }
      }
  }









};


#endif
