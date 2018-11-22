/* c routine to compute cubic spline interpolation from csa-c library (http://code.google.com/p/csa-c/) */
/*   adapetd from the stand alone program in the sources by F. Gillet-Chaulet */
#include <stdlib.h>

#include "csa.h"

void csa_interpolate_points(int nin, point pin[], int nout, point pout[], int nppc, int k)
{
double* std = NULL;
csa* a = csa_create();

csa_addpoints(a, nin, pin);
csa_addstd(a, nin, std);
if (nppc > 0)
      csa_setnppc(a, nppc);
if (k > 0)
      csa_setk(a, k);

csa_calculatespline(a);
csa_approximatepoints(a, nout, pout);
csa_destroy(a);
}
