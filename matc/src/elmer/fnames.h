/***********************************************************************
|
|  FNAMES.H - Last Edited 6. 8. 1988
|
***********************************************************************/

/* matc.c  */

char *doread( void );
VARIABLE *com_quit( void );

/*
 * $Id: fnames.h,v 1.2 2007/06/08 08:12:19 jpr Exp $ 
 *
 * $Log: fnames.h,v $
 * Revision 1.2  2007/06/08 08:12:19  jpr
 * *** empty log message ***
 *
 * Revision 1.1  2005/05/31 09:43:48  vierinen
 * oops
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:37  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

void error( char *fmt, ... );
void sig_trap(int);
int dogets(char *, char *);

void com_init( char *, int, int, VARIABLE *(*)(), int, int, char * );
void com_free( void ) ;
COMMAND *com_check(  char *);
VARIABLE *com_help( VARIABLE * );
VARIABLE *com_source( VARIABLE *);
VARIABLE *com_apply(  VARIABLE *);
VARIABLE *com_pointw(  double (*)(), VARIABLE *);
VARIABLE *com_el(  VARIABLE * );

void mem_free(void *);
void mem_free_all(void);
void *mem_alloc(size_t);

/* eval.c  */

VARIABLE *evaltree(  TREE *);
VARIABLE *evaltreelist(  TREE *);
VARIABLE *evalclause(  CLAUSE *);

VARIABLE *put_values( VARIABLE *, char *, VARIABLE *);
VARIABLE *put_result( VARIABLE *, char *, VARIABLE *, int, int);

/* files.c */

void fil_com_init(void);
VARIABLE *fil_fscanf( VARIABLE *);
VARIABLE *fil_fprintf( VARIABLE *);
VARIABLE *fil_fputs( VARIABLE *);
VARIABLE *fil_fopen( VARIABLE *);
VARIABLE *fil_freopen( VARIABLE *);
VARIABLE *fil_fclose( VARIABLE *);
VARIABLE *fil_save( VARIABLE *);
VARIABLE *fil_load( VARIABLE *);

/* funcs.c  */

FUNCTION *fnc_check( char *);
void fnc_free_entry( FUNCTION *);
VARIABLE *fnc_delete( VARIABLE *);
void fnc_free( void);
VARIABLE *fnc_exec( FUNCTION *, VARIABLE *);
void fnc_com_init( void);

/* jacobi.c */

VARIABLE *mtr_jacob(VARIABLE *);
int jacobi( double *, double *, double *, double *, double *, int, double);

/* lists.c */

void lst_addtail( int, LIST *);
void lst_addhead( int, LIST *);
void lst_add( int, LIST *);
void lst_unlink( int, LIST *);
void lst_free( int, LIST *);
LIST *lst_find( int, char *);
void lst_purge( int);
VARIABLE *lst_print( int);


/* matrix.c */

void mtr_com_init( void);

double func_abs( double);
VARIABLE *mtr_min( VARIABLE *);
VARIABLE *mtr_max( VARIABLE *);
VARIABLE *mtr_sum( VARIABLE *);
VARIABLE *mtr_trace( VARIABLE *);
VARIABLE *mtr_zeros( VARIABLE *);
VARIABLE *mtr_ones( VARIABLE *);
VARIABLE *mtr_rand( VARIABLE *);
VARIABLE *mtr_resize( VARIABLE *);
VARIABLE *mtr_vector( VARIABLE *);
VARIABLE *mtr_eye( VARIABLE *);
VARIABLE *mtr_size( VARIABLE *);

VARIABLE *mtr_LUD( VARIABLE *);
VARIABLE *mtr_det( VARIABLE *);
VARIABLE *mtr_inv( VARIABLE *);
void LUDecomp( double *, int, int *);

VARIABLE *mtr_eig( VARIABLE *);
VARIABLE *mtr_hesse( VARIABLE *);
void vbcalc( double *, double *, double *,int, int);
void hesse( double *, int, int);
void francis( double *, int, int);

/* oper.c */

MATRIX *mat_new( int, int, int);
MATRIX *mat_copy( MATRIX *);
void mat_free( MATRIX *);
MATRIX *opr_vector( MATRIX *, MATRIX * );
MATRIX *opr_resize( MATRIX *, MATRIX * );
MATRIX *opr_apply( MATRIX * );
MATRIX *opr_add( MATRIX *, MATRIX *);
MATRIX *opr_minus( MATRIX *);
MATRIX *opr_subs( MATRIX *, MATRIX *);
MATRIX *opr_mul( MATRIX *, MATRIX *);
MATRIX *opr_pmul( MATRIX *, MATRIX *);
MATRIX *opr_div( MATRIX *, MATRIX *);
MATRIX *opr_pow( MATRIX *, MATRIX *);
MATRIX *opr_trans( MATRIX *);
MATRIX *opr_reduction( MATRIX *, MATRIX *);
MATRIX *opr_lt( MATRIX *, MATRIX *);
MATRIX *opr_le( MATRIX *, MATRIX *);
MATRIX *opr_gt( MATRIX *, MATRIX *);
MATRIX *opr_ge( MATRIX *, MATRIX *);
MATRIX *opr_eq( MATRIX *, MATRIX *);
MATRIX *opr_neq( MATRIX *, MATRIX *);
MATRIX *opr_and( MATRIX *, MATRIX *);
MATRIX *opr_or( MATRIX *, MATRIX *);
MATRIX *opr_not( MATRIX *);


/* optimclause.c */

TREE *optimtree( TREE *);
CLAUSE *optimclause( CLAUSE *);

/* parser.c */

int char_in_list( int, char *);
void scan(void);
TREE *newtree(void);
TREE *args( int, int);
TREE *nameorvar( void);

TREE *par_trans(TREE *);
TREE *par_pow( TREE *);
TREE *par_timesdivide( TREE *);
TREE *par_plusminus(  TREE *);
TREE *par_compare( TREE *);
TREE *par_reduction( TREE *);

TREE *equation( void);
CLAUSE *statement( void);
CLAUSE *blockparse( void);
CLAUSE *funcparse( void);
CLAUSE *ifparse( void);
CLAUSE *whileparse( void);
CLAUSE *parse( void);

void free_treeentry( TREEENTRY *);
void free_tree( TREE *);
void free_clause( CLAUSE *);

VARIABLE *doit( char *);

/* printclause.c */

void printtree( TREE *, FILE *);
void printtreelist( TREE *, FILE *);
int printclause( CLAUSE *, FILE *, int);

/* urand.c */

double urand( int *);


/* VARIABLE.c */

void var_com_init( void);

VARIABLE *var_check( char *);
VARIABLE *var_varlist( void);
void var_print( VARIABLE *);

VARIABLE *var_temp_copy( VARIABLE *);
VARIABLE *var_temp_new( int, int, int);
void var_delete_temp( VARIABLE *);
void var_delete_temp_el( VARIABLE *);

VARIABLE *const_new( char *, int, int, int);
void const_free( void);
VARIABLE *var_new( char *, int, int, int);
VARIABLE *var_rename( VARIABLE *, char *);
void var_free( void);
void var_free_el( VARIABLE *);
void var_delete( char *);
char *var_to_string( VARIABLE *);

/* str.c */
void str_com_init(void);
VARIABLE *str_sprintf(VARIABLE *);
VARIABLE *str_sscanf(VARIABLE *);
VARIABLE *str_matcvt(VARIABLE *);
VARIABLE *str_cvtmat(VARIABLE *);

/* gra.c */
void gra_init(void);
void gra_close_sys(void);
void gra_set_viewport(double, double, double, double);
void gra_set_window(double, double, double, double, double, double);
void gra_perspective(double);
void gra_error(void);
void gra_window_to_viewport(double, double, double, double *, double *);
void gra_translate(double, double, double);
void gra_rotate(double, double, double);
void gra_scale(double, double, double);
void gra_viewpoint(double, double, double, double, double, double);
void gra_dbuffer_null(void);
void gra_com_init(void);
void gra_init_matc(int devtype, char *name);
void gra_mtrans(double,double,double,double *,double *,double *);

/*******************************************************************
                  graphics package definitions
*******************************************************************/ 
#include "gra.h"

/* gra.c */
void gra_getmatrix(GMATRIX);
void gra_setmatrix(GMATRIX);

/* clip.c */
int clip_poly(int *,double *,double *);
int clip_line(int *,double *,double *);

/* c3d.c */
VARIABLE *c3d_gc3d(VARIABLE *);
VARIABLE *c3d_gc3dlevels(VARIABLE *);

/* dri/dri_ps.c */

void gra_ps_open(int);
void gra_ps_close(void);
void gra_ps_clear(void);
void gra_ps_defcolor(int, double, double, double);
void gra_ps_color(int);
void gra_ps_polyline(int, Point *);
void gra_ps_draw(Point *);
void gra_ps_move(Point *);
void gra_ps_polymarker(int, int, Point *);
void gra_ps_marker(int, Point *);
void gra_ps_areafill(int, Point *);
void gra_ps_image(int, int, int, unsigned char *);
void gra_ps_text(double, double, char *);
void gra_ps_flush(void);
void gra_ps_reset(void);
