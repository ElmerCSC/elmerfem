
/*
  huti_defs.h - HUTIter libarary definitions 

  $Id: huti_defs.h,v 1.1.1.1 2005/04/15 10:31:18 vierinen Exp $

 */

#ifndef _H_HUTI_DEFS
#define _H_HUTI_DEFS

#define NUMBER_OF_PROCESSORS "HUTI_NP"

#define MAX_ERRMSGS 100              /* Maximum number of error messages */


/* HUTI defaults */

#define HUTI_DFLTMAXIT 5000
#define HUTI_DFLTTOLERANCE 0.000001
#define HUTI_DFLTDBUGLVL 0
#define HUTI_DFLTPROCS 1
#define HUTI_EPSILON 1.17549435E-38

/* HUTI status flags */

#define HUTI_OK 0
#define HUTI_CONVERGENCE 1
#define HUTI_MAXITER 2
#define HUTI_DIVERGENCE 3

/* QMR method */

#define HUTI_QMR_RHO_PSI 10
#define HUTI_QMR_DELTA 11
#define HUTI_QMR_EPSILON 12
#define HUTI_QMR_BETA 13
#define HUTI_QMR_GAMMA 14

/* CG method */

#define HUTI_CG_RHO 20

/* CGS method */

#define HUTI_CGS_RHO 25

/* TFQMR method */

#define HUTI_TFQMR_RHO 30

/* BiCGSTAB method */

#define HUTI_BICGSTAB_RHO 35
#define HUTI_BICGSTAB_SNORM 36
#define HUTI_BICGSTAB_OMEGA 37

/* GMRES method */

#define HUTI_GMRES_ALPHA 40
#define HUTI_GMRES_BETA 41

/* BiCGSTAB(2) method */

#define HUTI_BICGSTAB_2_RHO 45


/* HUTI debug levels */

#define HUTI_NO_DEBUG 0
#define HUTI_ITEROUTPUT 1

/* Initial X for solvers */

#define HUTI_RANDOMX 0
#define HUTI_USERSUPPLIEDX 1

/* Matrix type in external operations */

#define HUTI_MAT_NOTTRPSED 0
#define HUTI_MAT_TRPSED 1

/* Storage allocation for various methods */

#define HUTI_CG_WORKSIZE 4
#define HUTI_CGS_WORKSIZE 7
#define HUTI_BICGSTAB_WORKSIZE 8
#define HUTI_QMR_WORKSIZE 14
#define HUTI_TFQMR_WORKSIZE 10
#define HUTI_GMRES_WORKSIZE 7
#define HUTI_BICGSTAB_2_WORKSIZE 8

/* Different stopping criteria */

#define HUTI_TRUERESIDUAL 0
#define HUTI_TRESID_SCALED_BYB 1
#define HUTI_PSEUDORESIDUAL 2
#define HUTI_PRESID_SCALED_BYB 3
#define HUTI_PRESID_SCALED_BYPRECB 4
#define HUTI_XDIFF_NORM 5
#define HUTI_UPPERB_STOPC 6
#define HUTI_USUPPLIED_STOPC 10

/*
 * HUTI ipar structure (used for various parameters)
 */

#define HUTI_IPAR_DFLTSIZE 50

/*
 * Input parameters supplied by user or by initialization
 *
 * General parameters (0-8)
 */

#define HUTI_IPAR_LEN ipar[0]
#define HUTI_DPAR_LEN ipar[1]
#define HUTI_NDIM ipar[2]
#define HUTI_WRKDIM ipar[3]
#define HUTI_DBUGLVL ipar[4]
#define HUTI_EXTOP_MATTYPE ipar[5]

/* Iteration parameters (9-18) */

#define HUTI_MAXIT ipar[9]
#define HUTI_STOPC ipar[11]
#define HUTI_PCOND ipar[12]
#define HUTI_INITIALX ipar[13]
#define HUTI_GMRES_RESTART ipar[14]

/* Parallel environment parameters (19-28) */

#define HUTI_MYPROC ipar[19]
#define HUTI_PROCS ipar[20]

/*
 * Output parameters (30-39)
 */

#define HUTI_INFO ipar[29]
#define HUTI_ITERS ipar[30]

/*
 * HUTI dpar structure (used for various parameters)
 */

#define HUTI_DPAR_DFLTSIZE 10

/* Input parameters supplied by user */

#define HUTI_TOLERANCE dpar[0]
#define HUTI_MAXTOLERANCE dpar[1]


/* Parameters related to other solver in Elmer */
#define HUTI_SGSPARAM dpar[2]
#define HUTI_BICGSTABL_L ipar[15]


/*  
 * End of definitions
 */

#endif  /* _H_HUTI_DEFS */
