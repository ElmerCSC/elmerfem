
/*
 Internal definitions for HUTIter library

 $Id: huti_intdefs.h,v 1.5 2005/06/02 15:35:27 vierinen Exp $
*/

#define AUTOC "!!! This file is generated automatically, do not edit !!!"

#define _CONCAT(x,y)x/**/y

#define MAKE_INCLUDE(nsign,inc,incfile)nsign/**/inc incfile
#define MAKE_DEFINE(nsign,def,dname,dvalue)nsign/**/def dname dvalue

# ifdef S_PRE
#  define PRECISIONC s
#  define PRECISION_COMMENT Single precision
#  define F_PRECISION_TYPE real
#  define C_PRECISION_TYPE float
#  define NORMFUN_PREC_TYPE real
#  define MAKE_F_SUBRN( bn1, bn2, bn3, bn4 ) FC_FUNC_(bn2/**/s/**/bn4,bn1/**/S/**/bn3)
#  define MAKE_SUBRN( bn1, bn2 )  bn1/**/s/**/bn2
#  define PRECD_DUMMY_PCONDFUN  FC_FUNC_(huti_sdummy_pcondfun, HUTI_SDUMMY_PCONDFUN)
#  define PRECD_BLAS_DOTPRODFUN  FC_FUNC(sdot,SDOT)
#  define PRECD_BLAS_DOTPRODFUN_N FC_FUNC(sdot,SDOT)
#  define PRECD_BLAS_NORMFUN FC_FUNC(snrm2,SNRM2)
# endif

# ifdef D_PRE
#  define PRECISIONC d
#  define PRECISION_COMMENT Double precision
#  define F_PRECISION_TYPE double precision
#  define C_PRECISION_TYPE double
#  define NORMFUN_PREC_TYPE double precision
#  define MAKE_SUBRN( bn1, bn2 )  bn1/**/d/**/bn2
#  define MAKE_F_SUBRN( bn1, bn2, bn3, bn4 ) FC_FUNC_(bn2/**/d/**/bn4, bn1/**/D/**/bn3)
#  define PRECD_DUMMY_PCONDFUN FC_FUNC_(huti_ddummy_pcondfun, HUTI_DDUMMY_PCONDFUN)
#  define PRECD_BLAS_DOTPRODFUN FC_FUNC(ddot,DDOT)
#  define PRECD_BLAS_DOTPRODFUN_N FC_FUNC(ddot,DDOT)
#  define PRECD_BLAS_NORMFUN FC_FUNC(dnrm2,DNRM2)
# endif

# ifdef C_PRE
#  define PRECISIONC c
#  define PRECISION_COMMENT Complex
#  define F_PRECISION_TYPE complex
#  define C_PRECISION_TYPE float
#  define NORMFUN_PREC_TYPE real
#  define MAKE_SUBRN( bn1, bn2 )  bn1/**/c/**/bn2
#  define MAKE_F_SUBRN( bn1, bn2, bn3, bn4 ) FC_FUNC_(bn2/**/c/**/bn4, bn1/**/C/**/bn3)
#  define PRECD_DUMMY_PCONDFUN FC_FUNC_(huti_cdummy_pcondfun, HUTI_CDUMMY_PCONDFUN)
#  define PRECD_BLAS_DOTPRODFUN FC_FUNC(cdotu,CDOTU)
#  define PRECD_BLAS_DOTPRODFUN_N FC_FUNC(cdotc,CDOTC)
#  define PRECD_BLAS_NORMFUN FC_FUNC(scnrm2,SCNRM2)
# endif

# ifdef Z_PRE
#  define PRECISIONC z
#  define PRECISION_COMMENT Double complex
#  define F_PRECISION_TYPE double complex
#  define C_PRECISION_TYPE double
#  define NORMFUN_PREC_TYPE double precision
#  define MAKE_SUBRN( bn1, bn2 )  bn1/**/z/**/bn2
#  define MAKE_F_SUBRN( bn1, bn2, bn3, bn4 ) FC_FUNC_(bn2/**/z/**/bn4,bn1/**/Z/**/bn3)
#  define PRECD_DUMMY_PCONDFUN FC_FUNC_(huti_zdummy_pcondfun,HUTI_ZDUMMY_PCONDFUN)
#  define PRECD_BLAS_DOTPRODFUN FC_FUNC(zdotu,ZDOTU)
#  define PRECD_BLAS_DOTPRODFUN_N FC_FUNC(zdotc,ZDOTC)
#  define PRECD_BLAS_NORMFUN FC_FUNC(dznrm2,DZNRM2)
# endif

