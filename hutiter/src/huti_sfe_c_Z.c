










































#include  <stdlib.h> 
#include  <stdio.h> 
#include  "huti_defs.h" 
#include  "../config.h" 

extern void FC_FUNC_(huti_zcgsolv,HUTI_ZCGSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_zcgssolv,HUTI_ZCGSSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_zbicgstabsolv,HUTI_ZBICGSTABSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_zqmrsolv,HUTI_ZQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_ztfqmrsolv,HUTI_ZTFQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_zgmressolv,HUTI_ZGMRESSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_zbicgstab_2solv,HUTI_ZBICGSTAB_2SOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern int huti_num_of_procs;

extern void FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN) (void *u, void *v, int *ipar);



extern void FC_FUNC(zdotu,ZDOTU) (int *N, void *x, int *xind, void *y, int *yind, void *a);
extern void FC_FUNC(zdotc,ZDOTC) (int *N, void *x, int *xind, void *y, int *yind, void *a);
extern void FC_FUNC(dznrm2,DZNRM2) (int *N, void *x, int *xind);









void  STDCALLBULL FC_FUNC_(huti_z_cg,HUTI_Z_CG) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotu,ZDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_zcgsolv,HUTI_ZCGSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void  STDCALLBULL FC_FUNC_(huti_z_tfqmr,HUTI_Z_TFQMR) ( void *xvec, void *rhsvec,

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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotu,ZDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_ztfqmrsolv,HUTI_ZTFQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void  STDCALLBULL FC_FUNC_(huti_z_cgs,HUTI_Z_CGS) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotu,ZDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_zcgssolv,HUTI_ZCGSSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}








void  STDCALLBULL FC_FUNC_(huti_z_qmr,HUTI_Z_QMR) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotu,ZDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_zqmrsolv,HUTI_ZQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );
  return;
}







void  STDCALLBULL FC_FUNC_(huti_z_bicgstab,HUTI_Z_BICGSTAB) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotu,ZDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_zbicgstabsolv,HUTI_ZBICGSTABSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void  STDCALLBULL FC_FUNC_(huti_z_gmres,HUTI_Z_GMRES) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotc,ZDOTC);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_zgmressolv,HUTI_ZGMRESSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  






void  STDCALLBULL FC_FUNC_(huti_z_bicgstab_2,HUTI_Z_BICGSTAB_2) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(zdotu,ZDOTU);
  if (normfun == NULL)
    normfun = FC_FUNC(dznrm2,DZNRM2);

  FC_FUNC_(huti_zbicgstab_2solv,HUTI_ZBICGSTAB_2SOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
