/* nrutil.h */
/* Numerical Recipes' uncopyrighted vector and matrix allocation 
   and deallocation routines. */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

int MemoryUsage();
void nrerror(char error_text[]);

float *vector(int,int);
int  *ivector(int,int);
unsigned char *cvector(int,int);
unsigned long *lvector(int,int);
double *dvector(int,int);

float **matrix(int,int,int,int);
double **dmatrix(int,int,int,int);
int **imatrix(int,int,int,int);
float **submatrix(float **,int,int,int,int,int,int);

double ***f3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh);

void free_vector(float *,int,int);
void free_ivector(int *,int,int);
void free_cvector(unsigned char *,int,int);
void free_lvector(unsigned long *,int,int);
void free_dvector(double *,int,int);
 
void free_matrix(float **,int,int,int,int);
void free_dmatrix(double **,int,int,int,int);
void free_imatrix(int **,int,int,int,int);
void free_submatrix(float **,int,int,int,int);

void free_f3tensor(double ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh);

#endif
