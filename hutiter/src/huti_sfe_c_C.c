










































#include  <stdlib.h> 
#include  <stdio.h> 
#include  "huti_defs.h" 
#include  "../config.h" 

extern void FC_FUNC_(huti_ccgsolv, HUTI_CCGSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_ccgssolv, HUTI_CCGSSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_cbicgstabsolv, HUTI_CBICGSTABSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_cqmrsolv, HUTI_CQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_ctfqmrsolv, HUTI_CTFQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_cgmressolv, HUTI_CGMRESSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_cbicgstab_2solv, HUTI_CBICGSTAB_2SOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern int huti_num_of_procs;

extern void FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN) (void *u, void *v, int *ipar);



extern void FC_FUNC(cdotu,CDOTU) (int *N, void *x, int *xind, void *y, int *yind, void *blah);
extern void FC_FUNC(cdotc,CDOTC) (int *N, void *x, int *xind, void *y, int *yind, void *blah);
extern void FC_FUNC(scnrm2,SCNRM2) (int *N, void *x, int *xind);








void STDCALLBULL FC_FUNC_(huti_c_cg, HUTI_C_CG) ( void *xvec, void *rhsvec,
		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == NULL)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(cdotu,CDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_ccgsolv, HUTI_CCGSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void FC_FUNC_(huti_c_tfqmr, HUTI_C_TFQMR) ( void *xvec, void *rhsvec,

		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == NULL)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(cdotu,CDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_ctfqmrsolv, HUTI_CTFQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void FC_FUNC_(huti_c_cgs, HUTI_C_CGS) ( void *xvec, void *rhsvec,
		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == NULL)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(cdotu,CDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_ccgssolv, HUTI_CCGSSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}








void FC_FUNC_(huti_c_qmr, HUTI_C_QMR) ( void *xvec, void *rhsvec,
		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == NULL)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == 0)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == 0)
    dotprodfun = FC_FUNC(cdotu,CDOTU);
  if (normfun == 0)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_cqmrsolv, HUTI_CQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );
  return;
}







void FC_FUNC_(huti_c_bicgstab, HUTI_C_BICGSTAB) ( void *xvec, void *rhsvec,
		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == 0)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == 0)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == 0)
    dotprodfun = FC_FUNC(cdotu,CDOTU);
  if (normfun == 0)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_cbicgstabsolv, HUTI_CBICGSTABSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void FC_FUNC_(huti_c_gmres, HUTI_C_GMRES) ( void *xvec, void *rhsvec,
		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == NULL)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(cdotc,CDOTC);
  if (normfun == NULL)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_cgmressolv, HUTI_CGMRESSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  






void FC_FUNC_(huti_c_bicgstab_2, HUTI_C_BICGSTAB_2) ( void *xvec, void *rhsvec,
		int *ipar, double *dpar, void *work,
		void (*matvecsubr)(),
		void (*pcondlsubr)(),
		void (*pcondrsubr)(),
		void (*dotprodfun)(),
		void (*normfun)(),
		void (*mstopfun)() )
{
  HUTI_Init();

  

  if (pcondrsubr == NULL)
    pcondrsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(cdotu,CDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(scnrm2,SCNRM2);

  FC_FUNC_(huti_cbicgstab_2solv, HUTI_CBICGSTAB_2SOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
