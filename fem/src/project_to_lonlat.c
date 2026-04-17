/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/

/* Authors: Moritz Hanke (DKRZ), adapted by Thomas Zwinger
    original code provided by Moritz Hanke (DKRZ) in the frame of the TerraDT project
  Email:   thomas.zwinger@csc.fi
  Web:     http://www.csc.fi/elmer
  Address: CSC - IT Center for Science Ltd.
           Keilaranta 14
           02101 Espoo, Finland 

  Original Date: 10.6.2025
*/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include <proj.h>

#include "project_to_lonlat.h"

void convert_to_lonlat(double * x, double * y, int n, const char * crs) {
  // define transformation
  PJ * P =
    proj_create_crs_to_crs(
      PJ_DEFAULT_CTX, crs, "+proj=longlat +datum=WGS84", NULL);

  if (!P) {
    fputs("failed to create transformation", stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // transform all vertices
  for (int i = 0; i < n; ++i) {
    PJ_COORD src_coord = proj_coord(x[i], y[i], 0, 0);
    PJ_COORD tgt_coord = proj_trans(P, PJ_FWD, src_coord);
    x[i] = proj_torad(tgt_coord.lp.lam);
    y[i] = proj_torad(tgt_coord.lp.phi);
  }

  // clean up
  proj_destroy(P);
}
