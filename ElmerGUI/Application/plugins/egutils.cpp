/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter Råback
   Email: Peter.Raback@csc.fi
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


/* -------------------------------:  egutils.c  :----------------------------
   Includes common operations for operating vectors and such.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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


void nrerror(const char error_text[])
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


void bigerror(const char error_text[])
{
  fprintf(stderr,"The program encountered a major error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


void smallerror(const char error_text[])
{
  fprintf(stderr,"The program encountered a minor error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...we'll try to continue...\n");
  exit(1);
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
  int i,mini = 0;

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
  int i,maxi = 0;

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
  char *ptr1;

  strcpy(fname2,fname1); 

  //ML:
  return;

  ptr1 = strchr(fname2, '.');
  if (ptr1) *ptr1 = '\0';
  strcat(fname2, ".");
  strcat(fname2,newext);
}

int StringToStrings(const char *buf,char args[10][10],int maxcnt,char separator)
/*  Finds real numbers separated by a certain separator from a string.
    'buf'       - input string ending to a EOF
    'dest'      - a vector of real numbers
    'maxcnt'    - maximum number of real numbers to be read
    'separator'	- the separator of real numbers
    The number of numbers found is returned in the function value.
    */
{
  int i,cnt,totlen,finish;
  char *ptr1 = (char *)buf, *ptr2;
  

  totlen = strlen(buf);
  finish = 0;
  cnt = 0;

  if (!buf[0]) return 0;

  do {
    ptr2 = strchr(ptr1,separator);
    if(ptr2) {
      for(i=0;i<10;i++) {
	args[cnt][i] = ptr1[i];
	if(ptr1 + i >= ptr2) break;
      }
      args[cnt][i] = '\0';
      ptr1 = ptr2+1;
    }
    else {
      for(i=0;i<10;i++) {
	if(ptr1 + i >= buf+totlen) break;
	args[cnt][i] = ptr1[i];
      }
      args[cnt][i] = '\0';
      finish = 1;
    }
    
    cnt++;
  } while (cnt < maxcnt && !finish);
  
  return cnt;
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

  if (!buf[0]) return 0;
  do {
    ptr2 = strchr(ptr1,separator);
    if (ptr2) ptr2[0] = '\0';
    dest[cnt++] = atoi(ptr1);
    if (ptr2) ptr1 = ptr2+1;
  } while (cnt < maxcnt && ptr2 != NULL);

  return cnt;
}


int next_int(char **start)
{
  int i;
  char *end;

  i = strtol(*start,&end,10);
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



/* Indexing algorithm, Creates an index table */
#define SWAPI(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void SortIndex(int n,double *arr,int *indx)
{
  int i,indxt,ir,itemp,j,k,l;
  int jstack,*istack;
  double a;

  ir = n;
  l = 1;
  jstack = 0;  
  istack = ivector(1,NSTACK);

  for(j=1;j<=n;j++) 
    indx[j] = j;

  for(;;) {
    if (ir-l < M) {
      for(j=l+1;j<=ir;j++) {
	indxt = indx[j];
	a = arr[indxt];
	for(i=j-1;i>=1;i--) {
	  if(arr[indx[i]] <= a) break;
	  indx[i+1] = indx[i];
	}
	indx[i+1] = indxt;
      }
      if(jstack == 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    } 
    else {
      k = (l+ir) >>  1;
      SWAPI(indx[k],indx[l+1]);
      if(arr[indx[l+1]] > arr[indx[ir]]) {
	SWAPI(indx[l+1],indx[ir]);
      }
      if(arr[indx[l]] > arr[indx[ir]]) {
	SWAPI(indx[l],indx[ir]);
      }
      if(arr[indx[l+1]] > arr[indx[l]]) {
	SWAPI(indx[l+1],indx[l]);
      }
      i = l+1;
      j = ir;
      indxt = indx[l];
      a = arr[indxt];
      for(;;) {
	do i++; while(arr[indx[i]] < a);
	do j--; while(arr[indx[j]] > a);
	if(j < i) break;
	SWAPI(indx[i],indx[j]);
      }
      indx[l] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if(jstack > NSTACK) printf("NSTACK too small in SortIndex.");
      if(ir-i+1 >= j-l) {
	istack[jstack]   = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}




