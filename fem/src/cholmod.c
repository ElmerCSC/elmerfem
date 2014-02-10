#include "../config.h"

#ifdef HAVE_CHOLMOD
#include "cholmod.h"
#include "SuiteSparseQR_C.h"

typedef struct {
  cholmod_common c;
  cholmod_sparse a;
  cholmod_factor *l;
  int nz;
  double *z;
  SuiteSparseQR_C_factorization *qr;
} cholmod;


cholmod STDCALLBULL *FC_FUNC_(cholmod_ffactorize,CHOLMOD_FFACTORIZE)(int *n,int *rows,int *cols,double *vals)
{
  static double *y;
  int i,j,nz,ok,*p,*q;
  double *xx,*bx;

  cholmod *handle;

  handle = (cholmod *)calloc(sizeof(cholmod),1);
  cholmod_start(&handle->c);

  handle->a.nrow=*n;
  handle->a.ncol=*n;
  handle->a.p=rows;
  handle->a.i=cols;
  handle->a.x=vals;
  handle->a.packed=1;
  handle->a.sorted=1;
  handle->a.stype=-1;
  handle->a.nzmax=rows[*n];
  handle->a.xtype=CHOLMOD_REAL;

  handle->l=cholmod_analyze(&handle->a,&handle->c);
  cholmod_factorize(&handle->a,handle->l,&handle->c);

  return handle;
}

void STDCALLBULL FC_FUNC_(cholmod_fsolve,CHOLMOD_FSOLVE)(cholmod **handle,int *n,double *x, double *b)
{
  double *xx,*bb;
  int i;
  cholmod_dense *dx, *db;

  db = cholmod_zeros((*handle)->a.nrow, 1, (*handle)->a.xtype, &(*handle)->c);
  bb=db->x;
  for(i=0;i<*n;i++) bb[i]=b[i];

  dx=cholmod_solve(CHOLMOD_A,(*handle)->l,db,&(*handle)->c);
  xx=dx->x;
  for(i=0;i<*n;i++) x[i]=xx[i];

  cholmod_free_dense(&dx, &(*handle)->c);
  cholmod_free_dense(&db, &(*handle)->c);
}


void STDCALLBULL FC_FUNC_(cholmod_ffree,CHOLMOD_FFREE)(cholmod **handle)
{
  cholmod_free_factor(&(*handle)->l,&(*handle)->c);
  cholmod_finish(&(*handle)->c);
  free(*handle);
}


cholmod STDCALLBULL *FC_FUNC_(spqr_ffactorize,SPQR_FFACTORIZE)(int *n,int *rows,int *cols,double *vals)
{
  cholmod *handle;
  long int *pp,*ii;

  double *xx,*bb;
  int i,j,k,rank,nsize;
  cholmod_dense *dx, *db, *dr;

  handle = (cholmod *)calloc(sizeof(cholmod),1);
  cholmod_l_start(&handle->c);

  pp=(long int*)malloc(sizeof(long int)*(*n+1));
  ii=(long int*)malloc(sizeof(long int)*rows[*n]);
  for(j=0;j<=*n;j++) pp[j]=rows[j];
  for(j=0;j<rows[*n];j++) ii[j]=cols[j];

  handle->a.nrow=*n;
  handle->a.ncol=*n;
  handle->a.nz = NULL;
  handle->a.x=vals;
  handle->a.p=pp;
  handle->a.i=ii;
  handle->a.packed=1;
  handle->a.sorted=1;
  handle->a.stype= 0;
  handle->a.nzmax=rows[*n];
  handle->a.xtype=CHOLMOD_REAL;

  handle->qr=SuiteSparseQR_C_factorize(SPQR_ORDERING_DEFAULT,SPQR_DEFAULT_TOL,&handle->a,&handle->c);

  rank=handle->c.SPQR_istat[4];
  fprintf (stderr,"rank %ld %d\n", rank, *n) ;

  db = cholmod_l_zeros(handle->a.nrow, 1, handle->a.xtype, &handle->c);
  bb=db->x;

  nsize=*n-rank;
  handle->z=(double*)malloc(sizeof(double)*(*n)*nsize);

  j=0;
  for(i=*n-nsize;i<*n;i++)
  {
     bb[i]=1;
     dr=SuiteSparseQR_C_qmult(SPQR_QX,handle->qr,db,&handle->c);
     bb[i]=0;
     xx=dr->x;
     for(k=0;k<*n;k++) handle->z[j++]=xx[k];
     cholmod_l_free_dense(&dr, &handle->c);
  }
  handle->nz=nsize;
  cholmod_l_free_dense(&db, &handle->c);
  
  return handle;
}

void STDCALLBULL FC_FUNC_(spqr_fsolve,SPQR_FSOLVE)(cholmod **handle,int *n,double *x, double *b)
{
  double *xx,*bb;
  int i,j,k,rank,nsize;
  cholmod_dense *dx, *db, *dr;

  db = cholmod_l_zeros((*handle)->a.nrow, 1, (*handle)->a.xtype, &(*handle)->c);
  bb=db->x;
  for(i=0;i<*n;i++) bb[i]=b[i];

  dr=SuiteSparseQR_C_qmult(SPQR_QTX,(*handle)->qr,db,&(*handle)->c);
  dx=SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B,(*handle)->qr,dr,&(*handle)->c);
  xx=dx->x;
  for(i=0;i<*n;i++) x[i]=xx[i];

  cholmod_l_free_dense(&dr, &(*handle)->c);
  cholmod_l_free_dense(&dx, &(*handle)->c);
  cholmod_l_free_dense(&db, &(*handle)->c);
}

void STDCALLBULL FC_FUNC_(spqr_nz,SPQR_NZ)(cholmod **handle,int *nz)
{
   int i,j,k;
   *nz=(*handle)->nz;
}

void STDCALLBULL FC_FUNC_(spqr_nullspace,SPQR_NULLSPACE)(cholmod **handle,int *n,int *nz,double *z)
{
   int i,j,k;
   *nz=(*handle)->nz;
    k=0;
    for(i=0;i<*nz;i++)
    for(j=0;j<*n;j++)
    {
      z[j*(*nz)+i]=(*handle)->z[k++];
    }
}



int STDCALLBULL FC_FUNC_(spqr_ffree,SPQR_FFREE)(cholmod **handle)
{
  if(!(*handle)->qr) return -1;

  SuiteSparseQR_C_free(&(*handle)->qr,&(*handle)->c);
  if((*handle)->a.i) free((*handle)->a.i);
  if((*handle)->a.p) free((*handle)->a.p);
  if((*handle)->z) free((*handle)->z);
  cholmod_l_finish(&(*handle)->c);
  free(*handle);
  return 0;
}
#endif
