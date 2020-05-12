#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "nrutil.h"

#define FREE_ARG char*
#define SHOWMEM 0

#if SHOWMEM
static int nfloat=0, nint=0, cumnfloat=0, cumnint=0, locnint, locnfloat;

int MemoryUsage()
{
  if(locnint > 1000 || locnfloat > 1000) 
    printf("Memory used real %d %d %d int %d %d %d\n",
	   nfloat,cumnfloat,locnfloat,nint,cumnint,locnint);
  locnfloat = locnint = 0;
}
#endif


/* The following routines are copied from the book
   "Numerical Recipes in C, The art of scientific computing" 
   by Cambridge University Press and include the following

   Non-Copyright Notice: This appendix and its utility routines are 
   herewith placed into the public domain. Anyone may copy them freely
   for any purpose. We of course accept no liability whatsoever for 
   any such use. */
   


void nrerror(char error_text[])
/* standerd error handler */
{
  fprintf(stderr,"run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


/* Vector initialization */

float *vector(int nl,int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v = (float*)malloc((size_t) (nh-nl+1+1)*sizeof(float));
  if (!v) nrerror("allocation failure in vector()");
  return(v-nl+1);
}



int *ivector(int nl,int nh)
/* Allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  if( nh < nl ) {
    printf("Allocation impossible in ivector: %d %d\n",nl,nh);
    exit(1);
  }

  v=(int*) malloc((size_t) (nh-nl+1+1)*sizeof(int));
  if (!v) nrerror("allocation failure in ivector()");

 #if SHOWMEM
  locnint = (nh-nl+1+1);
  nint += locnint;
  cumnint += locnint;
  MemoryUsage();
#endif

  return(v-nl+1);
}


unsigned char *cvector(int nl,int nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;
  
  v=(unsigned char *)malloc((size_t) (nh-nl+1+1)*sizeof(unsigned char));
  if (!v) nrerror("allocation failure in cvector()");
  return(v-nl+1);
}


unsigned long *lvector(int nl,int nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;
  
  v=(unsigned long *)malloc((size_t) (nh-nl+1+1)*sizeof(unsigned long));
  if (!v) nrerror("allocation failure in lvector()");
  return(v-nl+1);
}



double *dvector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  if( nh < nl ) {
    printf("Allocation impossible in dvector: %d %d\n",nl,nh);
    exit(1);
  }

  v=(double *)malloc((size_t) (nh-nl+1+1)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");

#if SHOWMEM
  locnfloat = nh-nl+1+1;
  nfloat += locnfloat;
  cumnfloat += locnfloat;
  MemoryUsage();
#endif

  return(v-nl+1);
}


/* Matrix initialization */

float **matrix(int nrl,int nrh,int ncl,int nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  float **m;
  
  /* allocate pointers to rows */
  m=(float **) malloc((size_t) (nrow+1)*sizeof(float*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+1)*sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++)
    m[i]=m[i-1]+ncol;
  
  return(m);
}


double **dmatrix(int nrl,int nrh,int ncl,int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  if( nrh < nrl || nch < ncl ) {
    printf("Allocation impossible in dmatrix: %d %d %d %d\n",nrl,nrh,ncl,nch);
    exit(1);
  }


  
  /* allocate pointers to rows */
  m=(double **) malloc((size_t) (nrow+1)*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m += 1;
  m -= nrl;

  
  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++)
    m[i]=m[i-1]+ncol;

#if SHOWMEM
  locnfloat = (nrow+1 + (nrow*ncol+1));
  nfloat += locnfloat;
  cumnfloat += locnfloat;
  MemoryUsage();
#endif
  
  return(m);
} 


int **imatrix(int nrl,int nrh,int ncl,int nch)
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  int **m;

  if( nrh < nrl || nch < ncl ) {
    printf("Allocation impossible in dmatrix: %d %d %d %d\n",nrl,nrh,ncl,nch);
    exit(1);
  }
  
  /* allocate pointers to rows */
  m=(int **) malloc((size_t) (nrow+1)*sizeof(int*));
  if (!m) nrerror("allocation failure 1 in imatrix()");
  m += 1;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+1)*sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++)
    m[i]=m[i-1]+ncol;

#if SHOWMEM
  locnint = (nrow+1 + (nrow*ncol+1));
  nint += locnint;
  cumnint += locnint;
  MemoryUsage();
#endif
  
  return(m);
} 



float **submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  int i,j, nrow=oldrh-oldrl+1, ncol=oldcl-newcl;
  float **m;
  
  /* allocate array of pointers to rows */
  m=(float **) malloc((size_t) ((nrow+1)*sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += 1;
  m -= newrl;
  
  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) 
    m[j]=a[i]+ncol;
  
  return(m);
}

/* Tensor initialization */

double ***f3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  t=(double***) malloc((size_t)((nrow+1)*sizeof(double***)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += 1;
  t -= nrl;

  t[nrl]=(double**) malloc((size_t)((nrow*ncol+1)*sizeof(double*)));
  if(!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += 1;
  t[nrl] -= ncl;

  t[nrl][ncl]=(double*) malloc((size_t)((nrow*ncol*ndep+1)*sizeof(double)));
  if(!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += 1;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j] = t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i] = t[i-1]+ncol;
    t[i][ncl] = t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j] = t[i][j-1]+ndep;
  }
  return(t);
}


/* Deallocation routines */

void free_vector(float *v,int nl,int nh)
{
  free((FREE_ARG) (v+nl-1));
}

void free_ivector(int *v,int nl,int nh)
{
#if SHOWMEM
  cumnint -= (nh-nl+1+1);
#endif

  free((FREE_ARG) (v+nl-1));
} 

void free_cvector(unsigned char *v,int nl,int nh)
{
  free((FREE_ARG) (v+nl-1));
}

void free_lvector(unsigned long *v,int nl,int nh)
{
  free((FREE_ARG) (v+nl-1));
}


void free_dvector(double *v,int nl,int nh)
{
#if SHOWMEM
  cumnfloat -= (nh-nl+1+1);
#endif

  free((FREE_ARG) (v+nl-1));
}

void free_matrix(float **m,int nrl,int nrh,int ncl,int nch)
{
  free((FREE_ARG) (m[nrl]+ncl-1));
  free((FREE_ARG) (m+nrl-1));
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
#if SHOWMEM
  int nrow=nrh-nrl+1;
  int ncol=nch-ncl+1;
  cumnfloat -= (nrow+1 + (nrow*ncol+1));
#endif

  free((FREE_ARG) (m[nrl]+ncl-1));
  free((FREE_ARG) (m+nrl-1));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
#if SHOWMEM
  int nrow=nrh-nrl+1;
  int ncol=nch-ncl+1;
  cumnint -= (nrow+1 + (nrow*ncol+1));
#endif

  free((FREE_ARG) (m[nrl]+ncl-1));
  free((FREE_ARG) (m+nrl-1));
}

void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch)
{
  free((FREE_ARG) (b+nrl-1));
}

void free_f3tensor(double ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-1));
  free((FREE_ARG) (t[nrl]+ncl-1));
  free((FREE_ARG) (t+nrl-1));
}
