/* femtools.h */
/* Basic subroutines that are common to all FEM programs:
   shape functions, Gaussian quadratures and isoparametric 
   transformations etc. Only linear elements have been 
   implemented. */

void Squad404(Real *xi,Real *eta,Real *sfun, Real *lder);
void Squad303(Real *xi,Real *eta,Real *sfun, Real *lder);
void Squad408(Real *xi,Real *eta,Real *sfun, Real *lder);
void Squad409(Real *xi,Real *eta,Real *sfun, Real *lder);
void Squad202(Real *xi, Real *sfun, Real *lder);
void Squad203(Real *xi, Real *sfun, Real *lder);
int LocalToGlobalD2(Real *globalcoord,Real *shapeder,
		    Real *shapefunc,Real *xgauss, Real *ygauss, 
		    Real *det,Real *globalder,int nodesd2);
void LocalToGlobalD1(Real *globalcoord,Real *shapeder,
		     Real *shapefunc,Real *xgauss, Real *ygauss, 
		     Real *ratio,int nodesd1);
void SurfaceNormalD1(Real *coord,Real *normal,int nodesd1);
int GlobalToLocalD2(Real *coord,Real xglobal,Real yglobal,
		    Real *xlocal,Real *ylocal);
void GaussD1ToD2(int *sideind,int *elemind,Real xigaussd1,
		 Real *xigaussd2,Real *etagaussd2,int *outward); 
