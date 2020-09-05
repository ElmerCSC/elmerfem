/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter Raback
   Email: elmeradm@csc.fi
   Address: CSC - IT Center for Science Ltd.
            Keilaranta 14
            02101 Espoo, Finland

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <time.h> 

/* Possible monitoring of memory usage, if supported */
#define MEM_USAGE 0
#if MEM_USAGE
#include <sys/resource.h>
#endif

#include "egutils.h"

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
    printf("Allocation impossible in imatrix: %d %d %d %d\n",nrl,nrh,ncl,nch);
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


/* -------------------------------:  common.c  :----------------------------
   Includes common operations for operating vectors and such. */


static Real timer_t0, timer_dt;
static int timer_active = FALSE;
static char timer_filename[600];

void timer_init()
{
  timer_active = FALSE;
}


void timer_activate(const char *prefix)
{
  Real time;
  timer_active = TRUE; 

  time = clock() / (double)CLOCKS_PER_SEC;

  AddExtension(prefix,timer_filename,"time");

  printf("Activating timer (s): %.2f\n",time);
  printf("Saving timer info to file %s\n",timer_filename);

  timer_dt = time;
  timer_t0 = time;
}


void timer_show()
{
  static int visited = 0;
  Real time;
  FILE *out;
#if MEM_USAGE
  int who,ret;
  Real memusage;
  static struct rusage usage;
#endif

  if(!timer_active) return;

  time = clock() / (double)CLOCKS_PER_SEC;
  printf("Elapsed time (s): %.2f %.2f\n",time-timer_t0,time-timer_dt);

#if MEM_USAGE
  who = RUSAGE_SELF;
  ret = getrusage( who, &usage );
  if( !ret ) {
    printf("maxrss %ld\n",usage.ru_maxrss);
    printf("ixrss %ld\n",usage.ru_ixrss);
    printf("idrss %ld\n",usage.ru_idrss);
    printf("isrss %ld\n",usage.ru_isrss);

    memusage = (double) 1.0 * usage.ru_maxrss;
  }
  else {
    printf("Failed to obtain resource usage!\n");
    memusage = 0.0;
  }
#endif

  visited = visited + 1;
  if( visited == 1 ) {
    out = fopen(timer_filename,"w");
  }
  else {
    out = fopen(timer_filename,"a");    
  }

#if MEM_USAGE
  fprintf(out,"%3d %12.4le %12.4le %12.4le\n",visited,time-timer_t0,time-timer_dt,memusage);
#else
  fprintf(out,"%3d %12.4le %12.4le\n",visited,time-timer_t0,time-timer_dt);
#endif

  fclose(out);

  timer_dt = time; 
}




void bigerror(char error_text[])
{
  fprintf(stderr,"The program encountered a major error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


void smallerror(char error_text[])
{
  fprintf(stderr,"The program encountered a minor error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...we'll try to continue...\n");
}



int FileExists(char *filename)
{
  FILE *in;

  if((in = fopen(filename,"r")) == NULL) 
    return(FALSE);  
  else {
    fclose(in);
    return(TRUE);
  }
}


Real Minimum(Real *vector,int first,int last)
/* Returns the smallest value of vector in range [first,last]. */
{
  Real min;
  int i;

  min=vector[first];
  for(i=first+1;i<=last;i++)
    if(min>vector[i]) min=vector[i];

  return(min);
}


int Minimi(Real *vector,int first,int last)
/* Returns the position of the smallest value of vector in range [first,last]. */
{
  Real min;
  int i,mini;

  mini=-1;
  min=vector[first];
  for(i=first+1;i<=last;i++)
    if(min>vector[i]) 
      min=vector[(mini=i)];

  return(mini);
}


Real Maximum(Real *vector,int first,int last)
/* Returns the largest value of vector in range [first,last]. */
{
  Real max;
  int i;

  max=vector[first];
  for(i=first+1;i<=last;i++)
    if(max<vector[i]) max=vector[i];

  return(max);
}


int Maximi(Real *vector,int first,int last)
/* Returns the position of the largest value of vector in range [first,last]. */
{
  Real max;
  int i,maxi;

  maxi=-1;
  max=vector[first];
  for(i=first+1;i<=last;i++) 
    if(max<vector[i]) 
      max=vector[(maxi=i)];

  return(maxi);
}



void AddExtension(const char *fname1,char *fname2,const char *newext)
/* Changes the extension of a filename.
   'fname1' - the original filename
   'fname2' - the new filename
   'newext' - the new extension
   If there is originally no extension it's appended. In this case 
   there has to be room for the extension. 
   */
{
  char *ptr1,*ptr2;

  strcpy(fname2,fname1); 
  ptr1 = strrchr(fname2, '.');

  if (ptr1) {
    int badpoint=FALSE;
    ptr2 = strrchr(fname2, '/');
    if(ptr2 && ptr2 > ptr1) badpoint = TRUE;
    if(!badpoint) *ptr1 = '\0';
  }
  strcat(fname2, ".");
  strcat(fname2,newext);
}


int StringToReal(const char *buf,Real *dest,int maxcnt,char separator)
/*  Finds real numbers separated by a certain separator from a string.
    'buf'       - input string ending to a EOF
    'dest'      - a vector of real numbers
    'maxcnt'    - maximum number of real numbers to be read
    'separator'	- the separator of real numbers
    The number of numbers found is returned in the function value.
    */
{
  int cnt = 0;
  char *ptr1 = (char *)buf, *ptr2;
  
  if (!buf[0]) return 0;
  do {
    ptr2 = strchr(ptr1,separator);
    if (ptr2) ptr2[0] = '\0';
    dest[cnt++] = atof(ptr1);
    if (ptr2) ptr1 = ptr2+1;
  } while (cnt < maxcnt && ptr2 != NULL);

  return cnt;
}


int StringToInteger(const char *buf,int *dest,int maxcnt,char separator)
{
  int cnt = 0;
  char *ptr1 = (char *)buf, *ptr2;
  int ival;
  
  if (!buf[0]) return 0;
  do {
    
    ptr2 = strchr(ptr1,separator);
    if (ptr2) ptr2[0] = '\0';
    ival = atoi(ptr1);

    dest[cnt++] = ival;
    
    if (ptr2) ptr1 = ptr2+1;
  } while (cnt < maxcnt && ptr2 != NULL);

  return cnt;
}

int StringToIntegerNoZero(const char *buf,int *dest,int maxcnt,char separator)
{
  int cnt = 0;
  char *ptr1 = (char *)buf, *ptr2;
  int ival;
  
  if (!buf[0]) return 0;
  do {
    
    ptr2 = strchr(ptr1,separator);
    if (ptr2) ptr2[0] = '\0';
    ival = atoi(ptr1);

    if(ival == 0) break;
      
    dest[cnt++] = ival;      

    if (ptr2) ptr1 = ptr2+1;
  } while (cnt < maxcnt && ptr2 != NULL);

  return cnt;
}


int EchoFile(char *filename)
#define LINELENGTH 100
{
  FILE *in;
  char line[LINELENGTH];

  if((in = fopen(filename,"r")) == NULL) {
    printf("Could not open file '%s'.\n",filename);
    return(1);
  }

  while(fgets(line,LINELENGTH,in) != NULL)
    printf("%s",line);

  fclose(in);
  return(0);
}


int next_int(char **start)
{
  int i;
  char *end;

  i = strtol(*start,&end,10);
  *start = end;
  return(i);
}

int next_int_n(char **start, int n)
{
  int i;
  char *end = *start+n;
  char saved = *end;

  *end = '\0';
  i = strtol(*start,NULL,10);
  *end = saved;
  *start = end;
  return(i);
}

Real next_real(char **start)
{
  Real r;
  char *end;

  r = strtod(*start,&end);

  *start = end;
  return(r);
}


/*
 * sort: sort an (double) array to ascending order, and move the elements of
 *       another (integer) array accordingly. the latter can be used as track
 *       keeper of where an element in the sorted order at position (k) was in
 *       in the original order (Ord[k]), if it is initialized to contain
 *       numbers (0..N-1) before calling sort. 
 *
 * Parameters:
 *
 * N:      int                  / number of entries in the arrays.
 * Key:    double[N]             / array to be sorted.
 * Ord:    int[N]               / change this accordingly.
 */
void SortIndex( int N, double *Key, int *Ord )
{
  double CurrentKey;

  int CurrentOrd;

  int CurLastPos;
  int CurHalfPos;

  int i;
  int j; 
 
  /* Initialize order */
  for(i=1;i<=N;i++)
    Ord[i] = i;

  CurHalfPos = N / 2 + 1;
  CurLastPos = N;
  while( 1 ) {
    if ( CurHalfPos > 1 ) {
      CurHalfPos--;
      CurrentKey = Key[CurHalfPos];
      CurrentOrd = Ord[CurHalfPos];
    } else {
      CurrentKey = Key[CurLastPos];
      CurrentOrd = Ord[CurLastPos];
      Key[CurLastPos] = Key[1];
      Ord[CurLastPos] = Ord[1];
      CurLastPos--;
      if ( CurLastPos == 1 ) {
	Key[1] = CurrentKey;
	Ord[1] = CurrentOrd;
	return;
      }
    }
    i = CurHalfPos;
    j = 2 * CurHalfPos;
    while( j <= CurLastPos ) {
      if ( j < CurLastPos && Key[j] < Key[j + 1] ) {
	j++;
      }
      if ( CurrentKey < Key[j] ) {
	Key[i] = Key[j];
	Ord[i] = Ord[j];
	i  = j;
	j += i;
      } else {
	j = CurLastPos + 1;
      }
    }
    Key[i] = CurrentKey;
    Ord[i] = CurrentOrd;
  }

}



