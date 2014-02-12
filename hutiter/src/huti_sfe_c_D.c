










































#include  <stdlib.h> 
#include  <stdio.h> 
#include  "huti_defs.h" 
#include  "../config.h" 

extern void FC_FUNC_(huti_dcgsolv, HUTI_DCGSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_dcgssolv, HUTI_DCGSSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_dbicgstabsolv, HUTI_DBICGSTABSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_dqmrsolv, HUTI_DQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_dtfqmrsolv, HUTI_DTFQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_dgmressolv, HUTI_DGMRESSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_dbicgstab_2solv, HUTI_DBICGSTAB_2SOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern int huti_num_of_procs;

extern void FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN) (void *u, void *v, int *ipar);



extern void FC_FUNC(ddot,DDOT) (int *N, void *x, int *xind, void *y, int *yind);
extern void FC_FUNC(ddot,DDOT) (int *N, void *x, int *xind, void *y, int *yind);
extern void FC_FUNC(dnrm2,DNRM2) (int *N, void *x, int *xind);








void  STDCALLBULL FC_FUNC_(huti_d_cg, HUTI_D_CG) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dcgsolv, HUTI_DCGSOLV)
      ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun,
                  normfun,
                 mstopfun );

  return;
}
  







void  STDCALLBULL FC_FUNC_(huti_d_tfqmr, HUTI_D_TFQMR) ( void *xvec, void *rhsvec,

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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dtfqmrsolv, HUTI_DTFQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void  STDCALLBULL FC_FUNC_(huti_d_cgs, HUTI_D_CGS) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dcgssolv, HUTI_DCGSSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}








void  STDCALLBULL FC_FUNC_(huti_d_qmr, HUTI_D_QMR) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dqmrsolv, HUTI_DQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );
  return;
}







void  STDCALLBULL FC_FUNC_(huti_d_bicgstab, HUTI_D_BICGSTAB) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == 0)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == 0)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == 0)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dbicgstabsolv, HUTI_DBICGSTABSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void STDCALLBULL FC_FUNC_(huti_d_gmres, HUTI_D_GMRES) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dgmressolv, HUTI_DGMRESSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  






void STDCALLBULL FC_FUNC_(huti_d_bicgstab_2, HUTI_D_BICGSTAB_2) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(ddot,DDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(dnrm2,DNRM2);

  FC_FUNC_(huti_dbicgstab_2solv, HUTI_DBICGSTAB_2SOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
