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
#include <elmerparam.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i;
  char *tag;
  int *xi, ni, nr;
  double *xr, *tmp, *y;
  int n, m, len, nfun;

  tag = NULL;
  xi = NULL;
  xr = NULL;
  ni = nr = 0;
  nfun = 1;

  switch (nrhs) {

  case 4:
    if (mxIsDouble(prhs[2]) != 1)
      mexErrMsgTxt("4:th argument must be a number.");
    n = mxGetN(prhs[2]);
    m = mxGetM(prhs[2]);
    /*
    if (n != 1 || m != 1)
      mexErrMsgTxt("4:rd argument must be a scalar");
    */

    tmp = mxGetPr(prhs[3]);
    nfun = (int) *tmp;

    /* Fall through.  */
  case 3:
    n = mxGetN(prhs[2]);
    m = mxGetM(prhs[2]);
    len = n*m+1;
    if (len > 1) { 
      if (mxIsChar(prhs[2]) != 1)
        mexErrMsgTxt("3:th argument must be a string.");
      tag = mxCalloc(len, sizeof(char));
      mxGetString(prhs[2], tag, len);
    }

    /* Fall through.  */
  case 2:
    if (!mxIsDouble(prhs[1]))
      mexErrMsgTxt("2:rd argument must numerical.");
    n = mxGetN(prhs[1]);
    m = mxGetM(prhs[1]);
    if (n > 1 && m > 1)
      mexErrMsgTxt("2:rd argument must be a vector.");

    ni = n*m;
    xi = mxCalloc(ni, sizeof(int));
    tmp = mxGetPr(prhs[1]);
    for (i = 0; i < ni; i++)
      xi[i] = (int)(tmp[i] + 0.5);

    /* Fall through.  */
  case 1:
    if (!mxIsDouble(prhs[0]))
      mexErrMsgTxt("1:st argument must be numerical.");
    n = mxGetN(prhs[0]);
    m = mxGetM(prhs[0]);
    if (n > 1 && m > 1)
      mexErrMsgTxt("1:st argument must be a vector.");

    xr = mxGetPr(prhs[0]);
    nr = n*m;
  }

  plhs[0] = mxCreateDoubleMatrix(1,nfun, mxREAL);
  y = mxGetPr(plhs[0]);
  elmer_param_vec (nfun, y, nr, xr, ni, xi, tag);
}
