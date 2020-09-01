/* nrutil.h + common.h -> egutils.h */

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

void timer_init();
void timer_activate(const char *prefix);
void timer_show();

void bigerror(char error_text[]);
void smallerror(char error_text[]);
int  FileExists(char *filename);
Real Minimum(Real *vector,int first,int last);
int  Minimi(Real *vector,int first,int last);
Real Maximum(Real *vector,int first,int last);
int  Maximi(Real *vector,int first,int last);
void InspectVector(Real *vector,int first,int last,Real *min,
		   Real *max,int *mini,int *maxi);
int  Steepest(Real *vector,int first,int last);
Real MeanVector(Real *vector,int first,int last);
Real AbsMeanVector(Real *vector,int first,int last);
Real DifferVector(Real *vector1,Real *vector2,int first,int last);
void ReformVector(Real *vector1,int n1,Real *vector2,int n2);
void AdjustVector(Real max,Real min,Real *vector,int first,int last);
int  ReadRealVector(Real *vector,int first,int last,char *filename);
void SaveRealVector(Real *vector,int first,int last,char *filename);
int  ReadRealMatrix(Real **matrix,int row_first,int row_last,
		    int col_first,int col_last,char *filename);
void SaveRealMatrix(Real **matrix,int row_first,int row_last,
		    int col_first,int col_last,char *filename);
int  ReadIntegerVector(int *vector,int first,int last,char *filename);
void SaveIntegerVector(int *vector,int first,int last,char *filename);
int  ReadIntegerMatrix(int **matrix,int row_first,int row_last,
		       int col_first,int col_last,char *filename);
void SaveIntegerMatrix(int **matrix,int row_first,int row_last,
		       int col_first,int col_last,char *filename);
void SaveNonZeros(Real **matrix,int row_first,int row_last,
		  int col_first,int col_last,char *filename);
void AddExtension(const char *fname1,char *fname2,const char *newext);
int StringToReal(const char *buf,Real *dest,int maxcnt,char separator);
int StringToInteger(const char *buf,int *dest,int maxcnt,char separator);
int StringToIntegerNoZero(const char *buf,int *dest,int maxcnt,char separator);
int EchoFile(char *filename);
int next_int(char **start);
int next_int_n(char **start, int n);
Real next_real(char **start);
void SortIndex( int N, double *Key, int *Ord );
#endif
