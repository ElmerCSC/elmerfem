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
Real next_real(char **start);
void SortIndex( int N, double *Key, int *Ord );
#endif
