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

#ifndef PROJECT_TO_LONLAT_H
#define PROJECT_TO_LONLAT_H

/**
 * Convert coordinates from EPSG:3413 to longitude/latitude in radians.
 * 
 * @param x Array of x coordinates (EPSG:3413), modified in place to longitude in radians.
 * @param y Array of y coordinates (EPSG:3413), modified in place to latitude in radians.
 * @param n Number of coordinates in the arrays.
 */
void convert_epsg3413_to_lonlat(double * x, double * y, const int n);

#endif // PROJECT_TO_LONLAT_H
