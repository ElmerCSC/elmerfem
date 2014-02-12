/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_body.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Body 2D class, derived from Body base class. 
        
************************************************************************/

#ifndef _ECIF_BODY2D_
#define _ECIF_BODY2D_

#include "ecif_body.h"

class Body2D : public Body
{     
public:
  Body2D();
  Body2D(bodyGmtrType body_type, int ext_id, char* name, colorIndices color = DEFAULT_COLOR_INDEX);
  Body2D(ecif_Body_X& trx_body, bool add_default_layer = false);
  Body2D(bodyGmtrType body_type, int int_id, int ext_id, int nof_fem_elements, int* fem_elem_ids);
  bool acceptsStructuredMesh(int layer);
  int addAllPendingVertices(int layer);

protected:
  bool check();
} ;

#endif
