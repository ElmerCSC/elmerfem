/*
Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif
#include <stdlib.h>
/*
Constants
*/
#define NSCAL  7
#define NVEC3  1
#define NVEC2  1
#define MINHEIGHT 1.0e-6
#define LINESIZE 30
#define BLANK 1
/*
Prototypes
*/
void postgrid_(float  *xi,
	       float  *eta,
	       float  *z_c,
	       float  *z_t,
	       float  *deltaX,
	       int    *imax_in,
	       int    *jmax_in,
	       int    *kcmax_in,
	       int    *ktmax_in,
	       char   *runname,
	       char   *ergnum,
	       int    *maske,
	       int    *flag);
void pregrid_(float  *xi,
	      float  *eta,
	      float  *z_c,
	      float  *z_t,
	      int    *imax_in,
	      int    *jmax_in,
	      int    *kcmax_in,
	      int    *ktmax_in,
	      char   *runname,
	      char   *ergnum,
	      int    *maske,
	      float  *deltaX,
	      int    *flag);
void elmerdata_(int   *imax_in,
		int   *jmax_in,
		int   *kcmax_in,
		int   *ktmax_in,
		float *z_c,
		float *z_t,
		float *vx_c,
		float *vy_c,
		float *vz_c,
		float *age_c,
		float *temp_c,
		float *vx_t,
		float *vy_t,
		float *vz_t,
		float *temp_t_m,
		float *age_t,
		float *omega_t,
		float *Q_bm,
		float *Q_tld,
		float *am_perp,
		float *qx,
		float *qy,
		int   *n_cts,
		int   *maske,
		char  *runname,
		char  *ergnum,
		int   *flag);
void asciidata_(float  *xi,
		float  *eta,
		int   *imax_in,
		int   *jmax_in,
		int   *kcmax_in,
		int   *ktmax_in,
		float *z_c,
		float *z_t,

		float *vx_c,
		float *vy_c,
		float *vz_c,
		float *age_c,
		float *temp_c,
		float *vx_t,
		float *vy_t,
		float *vz_t,
		float *temp_t_m,
		float *age_t,
		float *omega_t,
		float *Q_bm,
		float *Q_tld,
		float *am_perp,
		float *qx,
		float *qy,
		int   *n_cts,
		int   *maske,
		char  *runname,
		char  *ergnum,
		int   *flag);
int get_staggered_grid(float  *xi,
		       float  *eta,
		       float  *z_in,
		       int    imax,
		       int    jmax,
		       int    kmax,
		       float  *deltaX,
		       float  *staggered_grid);
int get_interpolated_property_on_staggered_grid(int imax,
						int jmax,
						int kmax,
						float *property_in,
						float *property_out);
int get_glaciation_info(int imax,
			int jmax,
			int *iced,
			int *mask);
int get_glaciation_boundary_info(int imax,
				 int jmax,
				 int *iced,
				 int *boundary);
void  make_float_from_integer_scalar_field(int   *input_property,
					   float *output_property, 
					   int   number_of_nodes,
					   int   reorder_ice_land_sea_mask);
void readlog_c_(char   *runname,
		int    *imax,
		int    *jmax,
		int    *kcmax,
		int    *ktmax,
		int    *krmax,
		float  *deform,
		float  *deltaX,
		int    *gotit);
