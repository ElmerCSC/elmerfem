










































#include  <stdlib.h> 
#include  <stdio.h> 
#include  "huti_defs.h" 
#include  "../config.h" 

extern void FC_FUNC_(huti_scgsolv,HUTI_SCGSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_scgssolv,HUTI_SCGSSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_sbicgstabsolv,HUTI_SBICGSTABSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_sqmrsolv,HUTI_SQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_stfqmrsolv,HUTI_STFQMRSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_sgmressolv,HUTI_SGMRESSOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern void FC_FUNC_(huti_sbicgstab_2solv,HUTI_SBICGSTAB_2SOLV) ( int *ndim, int *wrkdim,
			void *xvec, void *rhsvec, int *ipar, double *dpar,
			void *work, void (*matvecsubr)(),
			void (*pcondlsubr)(), void (*pcondrsubr)(),
			void (*dotprodfun)(), void   (*normfun)(),
			void (*stopcfun)() );

extern int huti_num_of_procs;

extern void FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN) (void *u, void *v, int *ipar);



extern void FC_FUNC(sdot,SDOT) (int *N, void *x, int *xind, void *y, int *yind);
extern void FC_FUNC(sdot,SDOT) (int *N, void *x, int *xind, void *y, int *yind);
extern void FC_FUNC(snrm2,SNRM2) (int *N, void *x, int *xind);








void FC_FUNC_(huti_s_cg,HUTI_S_CG) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_scgsolv,HUTI_SCGSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void FC_FUNC_(huti_s_tfqmr,HUTI_S_TFQMR) ( void *xvec, void *rhsvec,

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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_stfqmrsolv,HUTI_STFQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void FC_FUNC_(huti_s_cgs,HUTI_S_CGS) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_scgssolv,HUTI_SCGSSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}








void FC_FUNC_(huti_s_qmr,HUTI_S_QMR) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_sqmrsolv,HUTI_SQMRSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );
  return;
}







void FC_FUNC_(huti_s_bicgstab,HUTI_S_BICGSTAB) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_sbicgstabsolv,HUTI_SBICGSTABSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  







void FC_FUNC_(huti_s_gmres,HUTI_S_GMRES) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_sgmressolv,HUTI_SGMRESSOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
  






void FC_FUNC_(huti_s_bicgstab_2,HUTI_S_BICGSTAB_2) ( void *xvec, void *rhsvec,
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
    pcondrsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (pcondlsubr == NULL)
    pcondlsubr = FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN);
  if (dotprodfun == NULL)
    dotprodfun = FC_FUNC(sdot,SDOT);
  if (normfun == NULL)
    normfun = FC_FUNC(snrm2,SNRM2);

  FC_FUNC_(huti_sbicgstab_2solv,HUTI_SBICGSTAB_2SOLV) ( &HUTI_NDIM, &HUTI_WRKDIM, xvec, rhsvec,
                 ipar, dpar, work, matvecsubr, pcondlsubr, pcondrsubr,
                 dotprodfun, normfun, mstopfun );

  return;
}
