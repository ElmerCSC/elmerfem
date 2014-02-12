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

double elmer_param(int nr, const double *xr,
                   int ni, const int *xi, const char *tag);

void elmer_param_vec(int nfun, double *fun, int nr, const double *xr,
                     int ni, const int *xi, const char *tag);
