/* common.h */
/* Common subroutines that operate on vectors, matrices and other basic 
   data types: Find the minimum or maximum place or value, find 
   the mean, calculate the mean difference, save to or load from 
   an external file etc. */

#ifndef _COMMON_H_
#define _COMMON_H_

typedef double Real;
#define Rvector       dvector
#define Ivector       ivector
#define Rmatrix       dmatrix
#define Imatrix       imatrix
#define free_Rvector  free_dvector
#define free_Ivector  free_ivector  
#define free_Rmatrix  free_dmatrix
#define free_Imatrix  free_imatrix
#define TRUE 1
#define FALSE 0


/* Numerical Recipes' uncopyrighted vector and matrix allocation 
   and deallocation routines. */
int MemoryUsage();
void nrerror(const char error_text[]);

float *vector(int,int);
int  *ivector(int,int);
unsigned char *cvector(int,int);
unsigned long *lvector(int,int);
double *dvector(int,int);

float **matrix(int,int,int,int);
double **dmatrix(int,int,int,int);
int **imatrix(int,int,int,int);

void free_vector(float *,int,int);
void free_ivector(int *,int,int);
void free_cvector(unsigned char *,int,int);
void free_lvector(unsigned long *,int,int);
void free_dvector(double *,int,int);
 
void free_matrix(float **,int,int,int,int);
void free_dmatrix(double **,int,int,int,int);
void free_imatrix(int **,int,int,int,int);

void bigerror(const char error_text[]);
void smallerror(const char error_text[]);
int  FileExists(char *filename);
Real Minimum(Real *vector,int first,int last);
int  Minimi(Real *vector,int first,int last);
Real Maximum(Real *vector,int first,int last);
int  Maximi(Real *vector,int first,int last);
void AddExtension(const char *fname1,char *fname2,const char *newext);
int StringToStrings(const char *buf,char argv[10][10],int argc,char separator);
int StringToReal(const char *buf,Real *dest,int maxcnt,char separator);
int StringToInteger(const char *buf,int *dest,int maxcnt,char separator);
int next_int(char **start);
Real next_real(char **start);
void SortIndex(int n,double *arr,int *indx);
#endif
