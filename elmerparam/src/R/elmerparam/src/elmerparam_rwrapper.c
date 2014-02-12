/*

  ElmerParam - A simple system for parametrized computing
 
  Copyright (C) 2006  CSC - IT Center for Science Ltd.

  Authors: Erik Edelmann <Erik.Edelmann@csc.fi>
           Peter Råback <Peter.Raback@csc.fi>
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
  along with this program (in file elmerparam/COPYING); if not, write to
  the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
  Boston, MA 02110-1301, USA.

 */

#include <R.h>
#include <Rdefines.h>
#include <ctype.h>
#include <elmerparam.h>

SEXP elmerparam_rwrapper (SEXP nfun, SEXP xr, SEXP xi, SEXP tag)
{
  int nr, ni;
  int *ip, *nfunp;
  double *rp, *yp, *tp;
  char *tagp;
  SEXP y;
  int i;
  FILE *fd;

  PROTECT(nfun = AS_INTEGER(nfun));
  PROTECT(xi = AS_INTEGER(xi));
  ni = LENGTH(xi);
  PROTECT(xr = AS_NUMERIC(xr));
  nr = LENGTH(xr);
  PROTECT(tag = AS_CHARACTER(tag));

  nfunp = INTEGER_POINTER(nfun);
  ip = INTEGER_POINTER(xi);
  rp = NUMERIC_POINTER(xr);
  tagp = CHAR(tag);

  PROTECT(y = NEW_NUMERIC(*nfunp));
  yp = NUMERIC_POINTER(y);

  #if 0
  fd = fopen("fuling", "w");
  for (i = 0; i <= 255; i++) {
    if (isprint(tagp[i]))
      fprintf(fd,"%i: %c\n", i, tagp[i]);
    else
      fprintf(fd, "%i: :%i:\n", i, tagp[i]);
  }
  fclose(fd);
  #endif

  /* I have no idea if this way of passing tagp is safe, but it seems to work.*/
  elmer_param_vec(*nfunp, yp, nr, rp, ni, ip, &tagp[32]);

  UNPROTECT(5);

  return y;
}
