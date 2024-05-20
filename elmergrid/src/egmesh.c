/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter Raback
   Email: elmeradm@csc.fi
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
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/* --------------------------:  egmesh.c  :----------------------------

   This module includes subroutines that formulate the mesh into structures 
   more useful for the user. These are the functions that should be used for
   assembling. The routines usually operate on structures FemType and 
   BoundaryType.
   */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egnative.h"
#include "egmesh.h"

#define DEBUG 0


void GetElementInfo(int element,struct FemType *data,
		    Real *globalcoord,int *ind,int *material)
/* For a given element gives the coordinates and index for the knot.
   This subroutine uses the standard formulation which uses up more
   memory, but is easier to understand. It requires that the knots
   must first be stored in struct FemType.
   */
{
  int i,indi,nodesd2;
  nodesd2 = data->elementtypes[element]%100; 

  for(i=0;i<nodesd2;i++) {
    indi = ind[i] = data->topology[element][i];
    globalcoord[i] = data->x[indi];
    globalcoord[i+nodesd2] = data->y[indi];
  }
  (*material) = data->material[element];
}

int GetElementDimension(int elementtype)
{
  int elemdim;

  elemdim = 0;

  switch (elementtype / 100) {
  case 1:
    elemdim = 0;
    break;
  case 2: 
    elemdim = 1;
    break;
  case 3:
  case 4:
    elemdim = 2;
    break;
  case 5:
  case 6:
  case 7:
  case 8:
    elemdim = 3;
    break;
  default:
    printf("GetElementDimension: unknown elementtype %d\n",elementtype); 
    
  }
  return(elemdim);
}


int GetMaxElementType(struct FemType *data)
{
  int i,maxelementtype;

  maxelementtype = data->elementtypes[1];
  for(i=1;i<=data->noelements;i++)
    if(data->elementtypes[i] > maxelementtype)
      maxelementtype = data->elementtypes[i];

  return(maxelementtype);
}


int GetMinElementType(struct FemType *data)
{
  int i,minelementtype;

  minelementtype = data->elementtypes[1];
  for(i=1;i<=data->noelements;i++)
    if(data->elementtypes[i] < minelementtype)
      minelementtype = data->elementtypes[i];

  return(minelementtype);
}


int GetMaxElementDimension(struct FemType *data)
{
  int maxelementtype,elemdim;

  maxelementtype = GetMaxElementType(data);
  elemdim = GetElementDimension(maxelementtype);
  return(elemdim);
}


int GetCoordinateDimension(struct FemType *data,int info)
{
  int i,j,noknots,coorddim;
  int coordis;
  Real *coord;
  Real epsilon = 1.0e-20;

  noknots = data->noknots;
  coorddim = 0;

  for(j=3;j>=1;j--) {
    coordis = FALSE;
    if( j==1 ) 
      coord = data->x;
    else if( j==2 ) 
      coord = data->y;
    else
      coord = data->z;
    
    for(i=1;i<=noknots;i++) 
      if( fabs( coord[i] ) > epsilon ) {
	coordis = TRUE; 
	break;
      }
    if( coordis ) coorddim = MAX( coorddim, j ); 
  }
  if(info) printf("Coordinates defined in %d dimensions\n",coorddim);

  return(coorddim);
}


void GetElementSide(int element,int side,int normal,
		    struct FemType *data,int *ind,int *sideelemtype)
/* Give the indices of a given side of a given element. 
   The subroutine is valid for 4, 5, 8, 9, 12 and 16  
   node rectangular elements, and 3 and 6 node triangular 
   elements.
   */
{
  int i,j,elemtype,*elemind=NULL,sides,ind2[MAXNODESD2];

  /* if(element < 1 || element > data->noelements ) {
    printf("Invalid index for element: %d\n",element);
    bigerror("Cannot continue");
    } */
  
  elemtype = data->elementtypes[element];
  elemind = data->topology[element];
  sides = elemtype/100;
  *sideelemtype = 0;

  if(side < 0 && sides > 4) 
    side = -(side+1);
  
  switch (elemtype) {
  case 202:
  case 203:
  case 204:
    *sideelemtype = 101;    
    ind[0] = elemind[side];
    break;

  case 303: /* Linear triangle */
    if(side < 3) {
      *sideelemtype = 202;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
    }
    else if( side < 6 ) {
      *sideelemtype = 101;
      ind[0] = elemind[side-3];
    }
    break;

  case 306: /* 2nd order triangle */
    if(side < 3) {
      *sideelemtype = 203;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
      ind[2] = elemind[side+3];
    }
    else if( side < 9 ) {
      *sideelemtype = 101;
      ind[0] = elemind[side-3];
    }
    break;

  case 310: /* 3rd order triangle */
    if(side < 3) {
      *sideelemtype = 204;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
      ind[2] = elemind[2*side+3];
      ind[3] = elemind[2*side+4];            
    }
    else if( side < 13) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-3];
    }
    break;

  case 404: /* Linear quadrilateral */
    if(side < 4) {
      *sideelemtype = 202;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
    }
    else if(side < 8) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-4];
    }
    break;

  case 405:
    if(side < 4) {
      *sideelemtype = 202;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
    }
    else if(side < 9) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-4];
    }
    break;


  case 408: /* 2nd order quadrilateral */
    if(side < 4) {
      *sideelemtype = 203;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[side+4];
    }
    else if(side < 12) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-4];
    }
    break;

  case 409:
    if(side < 4) {
      *sideelemtype = 203;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[side+4];
    }
    else if(side < 13) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-4];
    }
    break;

  case 412: /* 3rd order quadrilateral */
    if(side < 4) {
      *sideelemtype = 204;      
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[2*side+4];
      ind[3] = elemind[2*side+5];      
    }
    else if(side < 16) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-4];
    }
    break;

  case 416:
    if(side < 4) {
      *sideelemtype = 204;      
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[2*side+4];
      ind[3] = elemind[2*side+5];      
    }
    else if(side < 20) {
      *sideelemtype = 101;      
      ind[0] = elemind[side-4];
    }
    break;

  case 504: /* Linear tetrahedron */
    if(side < 4) {
      *sideelemtype = 303;
      if(side < 3) {
	ind[0] = elemind[side];
	ind[1] = elemind[(side+1)%3];
	ind[2] = elemind[3];
      }
      if(side == 3) {
	ind[0] = elemind[0];
	ind[1] = elemind[2];
	ind[2] = elemind[1];	
      }
    }
    else if(side < 10) {
      *sideelemtype = 202;
      if(side < 7) {
	ind[0] = elemind[side-4];
	ind[1] = elemind[3];
      }
      else {
	ind[0] = elemind[side-7];
	ind[1] = elemind[(side-6)%3];	
      }
    }
    else if(side < 14) {
      *sideelemtype = 101;
      ind[0] = elemind[side-10];
    }      
    break;

  case 510: /* 2nd order tetrahedron */

    if(side < 4) {
      *sideelemtype = 306;
      if(side < 3) {
	ind[0] = elemind[side];
	ind[1] = elemind[(side+1)%3];
	ind[2] = elemind[3];
	ind[3] = elemind[4+side];
	ind[4] = elemind[7+(side+1)%3];
	ind[5] = elemind[7+side];
      }	
      else if(side == 3) {
	ind[0] = elemind[0];
	ind[1] = elemind[1];
	ind[2] = elemind[2];	
	ind[3] = elemind[4];
	ind[4] = elemind[5];
	ind[5] = elemind[6];	
      }
    }
    else if(side < 10) {
      *sideelemtype = 203;
      if(side < 7) {
	ind[0] = elemind[side-4];
	ind[1] = elemind[3];
	ind[2] = elemind[side+3];
      }
      else {
	ind[0] = elemind[side-7];
	ind[1] = elemind[(side-6)%3];	
	ind[2] = elemind[side-3];
      }
    }    
    else if(side < 20) {
      *sideelemtype = 101;
      ind[0] = elemind[side-10];
    }      

    break;

  case 706: /* Linear wedge element */
    if(side < 3) {
      *sideelemtype = 404;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
      ind[2] = elemind[(side+1)%3+3];
      ind[3] = elemind[side+3];  
    }
    else if (side < 5) {
      *sideelemtype = 303;          
      for(i=0;i<3;i++)
	ind[i] = elemind[3*(side-3)+i];
    }
    else if(side < 14) {
      *sideelemtype = 202;
      if(side < 8) {
	ind[0] = elemind[side-5];
	ind[1] = elemind[(side-4)%3];
      }
      if(side < 11) {
	ind[0] = elemind[3+side-8];
	ind[1] = elemind[3+(side-7)%3];
      }
      else {
	ind[0] = elemind[side-11];
	ind[1] = elemind[3+side-11];	
      }
    }
    else if (side < 20) {
      *sideelemtype = 101;
      ind[0] = elemind[side-14];      
    }
    break;

  case 715: /* Quadratic wedge element */
    if(side < 3) {
      *sideelemtype = 408;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
      ind[2] = elemind[(side+1)%3+3];
      ind[3] = elemind[side+3];  
      ind[4] = elemind[6+side];
      ind[5] = elemind[9+(side+1)%3];
      ind[6] = elemind[12+side];
      ind[7] = elemind[9+side];      
    }
    else if (side < 5) {
      *sideelemtype = 306;          
      for(i=0;i<3;i++) {
	ind[i] = elemind[3*(side-3)+i];
	ind[i+3] = elemind[3*(side-3)+6+i];
      }      
    }
    else if(side < 14) {
      *sideelemtype = 202;
      if(side < 8) {
	ind[0] = elemind[side-5];
	ind[1] = elemind[(side-4)%3];
      }
      if(side < 11) {
	ind[0] = elemind[3+side-8];
	ind[1] = elemind[3+(side-7)%3];
      }
      else {
	ind[0] = elemind[side-11];
	ind[1] = elemind[3+side-11];	
      }
    }
    else if (side < 20) {
      *sideelemtype = 101;
      ind[0] = elemind[side-14];      
    }
    break;

  case 605: /* Linear pyramid */
    if(side < 4) {
      *sideelemtype = 303;
      ind[0] = elemind[side];
      ind[1] = elemind[4];
      ind[2] = elemind[(side+1)%4];
    }
    else if (side < 5) {
      *sideelemtype = 404;     
      for(i=0;i<4;i++)
	ind[i] = elemind[i];
    }
    else if(side < 13) {
      *sideelemtype = 202;
      if(side < 9) {
	ind[0] = elemind[side-5];
	ind[1] = elemind[(side-4)%4];
      }
      else {
	ind[0] = elemind[side-9];
	ind[1] = elemind[4];
      }
    }
    else if(side < 18) {
      *sideelemtype = 101;
      ind[0] = elemind[side-13];            
    }
    break;

  case 613: /* 2nd order pyramid */
    if(side < 4) {
      *sideelemtype = 306;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[4];

      ind[3] = elemind[side+5];
      ind[4] = elemind[(side+1)%4+9];
      ind[5] = elemind[side%4+9];
    }
    else if (side == 4) {
      *sideelemtype = 408;     
      for(i=0;i<4;i++)
	ind[i] = elemind[i];
      for(i=0;i<4;i++)
	ind[i+4] = elemind[i+5];
    }
    else if(side < 13) {
      *sideelemtype = 203;
      if(side < 9) {
	ind[0] = elemind[(side-5)];
	ind[1] = elemind[(side-4)%4];
	ind[2] = elemind[side];
      }
      else {
	ind[0] = elemind[side-9];
	ind[1] = elemind[4];
	ind[2] = elemind[side];
      }
    }
    else if(side < 26) {
      *sideelemtype = 101;
      ind[0] = elemind[side-13];            
    }
    break;

  case 808: /* Linear brick */
    if(side < 6) {
      *sideelemtype = 404;
      if(side < 4) {
	ind[0] = elemind[side];
	ind[1] = elemind[(side+1)%4];
	ind[2] = elemind[(side+1)%4+4];
	ind[3] = elemind[side+4];      
      }
      else if(side < 6) {
	for(i=0;i<4;i++)
	  ind[i] = elemind[4*(side-4)+i];
      }
    }
    else if(side < 18) {
      *sideelemtype = 202;
      if(side < 10) {
	ind[0] = elemind[side-6];
	ind[1] = elemind[(side-5)%4];
      }
      else if(side < 14) {
	ind[0] = elemind[side-6];
	ind[1] = elemind[(side-9)%4+4];      
      }
      else {
	ind[0] = elemind[side-14];
	ind[1] = elemind[side-14+4];      
      }
    }
    else if(side < 26) {
      *sideelemtype = 101;
      ind[0] = elemind[side-18];
    }
    break;

  case 820: /* 2nd order brick */
    if(side < 4) {
      *sideelemtype = 408;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[(side+1)%4+4];
      ind[3] = elemind[side+4];      
      ind[4] = elemind[8+side];
      ind[5] = elemind[12+(side+1)%4];
      ind[6] = elemind[16+side];
      ind[7] = elemind[12+side];      
    }
    else if(side < 6) {
      *sideelemtype = 408;
      for(i=0;i<4;i++)
	ind[i] = elemind[4*(side-4)+i];
      for(i=0;i<4;i++)
	ind[i+4] = elemind[8*(side-4)+8+i];
    }
    break;

  case 827: 
    if(side < 4) {
      *sideelemtype = 409;
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
      ind[2] = elemind[(side+1)%4+4];
      ind[3] = elemind[side+4];      
      ind[4] = elemind[8+side];
      ind[5] = elemind[12+(side+1)%4];
      ind[6] = elemind[16+side];
      ind[7] = elemind[12+side];      
      ind[8] = elemind[20+side];
    }
    else {
      *sideelemtype = 409;
      for(i=0;i<4;i++)
	ind[i] = elemind[4*(side-4)+i];
      for(i=0;i<4;i++)
	ind[i+4] = elemind[8*(side-4)+8+i];
      ind[8] = elemind[20+side];
    }
    break;

  default:
    printf("GetElementSide: unknown elementtype %d (elem=%d,side=%d)\n",elemtype,element,side); 
    bigerror("Cannot continue");
  }

  if(normal == -1) {
    if(*sideelemtype == 202 || *sideelemtype == 203 || *sideelemtype == 303 || *sideelemtype == 404) {
      j = *sideelemtype/100-1;
      for(i=0;i<=j;i++)
	ind2[i] = ind[i];
      for(i=0;i<=j;i++)
	ind[i] = ind2[j-i];
    }
  }
}



void GetBoundaryElement(int sideind,struct BoundaryType *bound,struct FemType *data,int *ind,int *sideelemtype)
{
  int element,side,normal,i,n;

  if( sideind > bound->nosides ) {
    *sideelemtype = 0;
    printf("Side element index %d exceeds size of boundary (%d)\n",sideind,bound->nosides);
    return;
  }

  element = bound->parent[sideind];


  /*GetElementSide(elemind2,side,1,data,&sideind2[0],&sideelemtype2); */

  if(element) {
    side = bound->side[sideind];
    normal = bound->normal[sideind];
    GetElementSide(element,side,normal,data,ind,sideelemtype);
  }
  else {
    *sideelemtype = bound->elementtypes[sideind];

    n = *sideelemtype % 100;
    
    for(i=0;i<n;i++)
      ind[i] = bound->topology[sideind][i];

    if(0) {
      printf("sidelemtype = %d\n",*sideelemtype);
      printf("ind = ");
      for(i=0;i<n;i++) printf("%d ",ind[i]);
      printf("\n");
    }
  }
}



int GetElementFaces(int elemtype) 
{
  int basetype=0,elemfaces=0;

  basetype = elemtype / 100;

  switch (basetype) {
  case 1:
    elemfaces = 0;
    break;
  case 2:
    elemfaces = 2;
    break;
  case 3:
    elemfaces = 3;
    break;
  case 4:
    elemfaces = 4;
    break;
  case 5:
    elemfaces = 4;
    break;
  case 6:
    elemfaces = 5;
    break;
  case 7:
    elemfaces = 5;
    break;
  case 8:
    elemfaces = 6;
    break;

  default:
    printf("GetElementFaces: Unknown elementtype %d\n",elemfaces);
  }

  return(elemfaces);
}






int GetElementGraph(int element,int edge,struct FemType *data,int *ind)
{
  int elemtype,basetype,elemnodes;
  int hit,evenodd,quadratic,side;
  int *elemind=NULL;

  elemtype = data->elementtypes[element];
  basetype = elemtype / 100;
  elemnodes = elemtype % 100;
  quadratic = (elemnodes > basetype);
  elemind = data->topology[element];

  ind[0] = ind[1] = 0;
  
  if(quadratic) 
    side = edge / 2;
  else
    side = edge;
  
  
  switch (basetype) {
  case 2:
    if(side == 0) {
      ind[0] = elemind[0];
      ind[1] = elemind[1];
    }
    break;    
  case 3:
    if(side < 3) {
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
    }
    break;
  case 4:
    if(side < 4) {
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
    }
    break;
  case 5:
    if(side < 3) {
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%3];
    }
    else if(side < 6) {
      ind[0] = elemind[side-3];
      ind[1] = elemind[3];
    }
    break;
  case 6:
    if(side < 4) {
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
    }
    else if(side < 8) {
      ind[0] = elemind[side-4];
      ind[1] = elemind[4];
    }
    break;
  case 7:
    switch(side) {
       case 0: ind[0]=elemind[0]; ind[1]=elemind[1]; break;
       case 1: ind[0]=elemind[1]; ind[1]=elemind[2]; break;
       case 2: ind[0]=elemind[2]; ind[1]=elemind[0]; break;
       case 3: ind[0]=elemind[3]; ind[1]=elemind[4]; break;
       case 4: ind[0]=elemind[4]; ind[1]=elemind[5]; break;
       case 5: ind[0]=elemind[5]; ind[1]=elemind[3]; break;
       case 6: ind[0]=elemind[0]; ind[1]=elemind[3]; break;
       case 7: ind[0]=elemind[1]; ind[1]=elemind[4]; break;
       case 8: ind[0]=elemind[2]; ind[1]=elemind[5]; break;
    }
    break;

  case 8:
    if(side < 4) {
      ind[0] = elemind[side];
      ind[1] = elemind[(side+1)%4];
    }
    else if(side < 8) {
      ind[0] = elemind[side-4];
      ind[1] = elemind[side];
    }
    else if(side < 12) {
      ind[0] = elemind[side-4];
      ind[1] = elemind[4+(side+1)%4];
    }
    break;
  }

  hit = (ind[0] || ind[1]);
  
  
  if(hit && quadratic) {
    evenodd = edge - 2*side;
    
    switch (basetype) {
    case 2:
      ind[evenodd] = elemind[2];
      break;
      
    case 3:
      ind[evenodd] = elemind[side+3];
      break;
      
    case 4:
      ind[evenodd] = elemind[side+4];
      break;
      
    case 5:
      ind[evenodd] = elemind[side+4];
      break;
      
    case 6:
      ind[evenodd] = elemind[side+5];
      break;

    case 7:
      ind[evenodd] = elemind[side+6];
      break;
      
    case 8:
      ind[evenodd] = elemind[side+8];
      break;
      
    }
  }

  return(hit);
}




int CalculateIndexwidth(struct FemType *data,int indxis,int *indx)
{
  int i,ind,nonodes,indexwidth;
  int imax,imin,element;

  /* Calculate the maximum bandwidth */
  
  indexwidth = 0;

  for(element=1; element <= data->noelements; element++) {
    imin = data->noknots;
    imax = 0;
    nonodes = data->elementtypes[element]%100;
    for(i=0;i<nonodes;i++) {
      ind = data->topology[element][i];
      if(indxis) ind = indx[ind];
      if(ind == 0) continue;
      if(ind > imax) imax = ind;
      if(ind < imin) imin = ind;
    }
    if(imax-imin > indexwidth) 
      indexwidth = imax-imin;
  }

  if(!indxis) data->indexwidth = indexwidth;
  return(indexwidth);
}


void InitializeKnots(struct FemType *data) 
{
  int i;

  data->timesteps = 0;
  data->noknots = 0;
  data->noelements = 0;
  data->coordsystem = COORD_CART2;
  data->numbering = NUMBER_XY;
  data->created = FALSE;
  data->variables = 0;
  data->maxnodes = 0;
  data->indexwidth = 0;
  data->noboundaries = 0;
  data->mapgeo = 1;
  data->nocorners = 0;

  data->boundarynamesexist = FALSE;
  data->bodynamesexist = FALSE;

  data->nodepermexist = FALSE;  
  
  data->nopartitions = 1;
  data->partitionexist = FALSE;
  data->periodicexist = FALSE;
  data->nodeconnectexist = FALSE;
  data->elemconnectexist = FALSE;

  data->nodalexists = FALSE;
  /* data->invtopoexists = FALSE; */
  data->partitiontableexists = FALSE;

  data->invtopo.created = FALSE;
  data->nodalgraph2.created = FALSE;
  data->dualgraph.created = FALSE;

  
  for(i=0;i<MAXDOFS;i++) {
    data->edofs[i] = 0;
    data->bandwidth[i] = 0;
    data->iterdofs[i] = 0;
    strcpy(data->dofname[i],""); 
  }

  for(i=0;i<MAXBODIES;i++) {
    data->bodyname[i] = NULL;
  }
  for(i=0;i<MAXBCS;i++) {
    data->boundaryname[i] = NULL;
  }
}


void AllocateKnots(struct FemType *data)
{
  int i;

  data->topology = Imatrix(1,data->noelements,0,data->maxnodes-1);
  data->material = Ivector(1,data->noelements);
  data->elementtypes = Ivector(1,data->noelements);

  for(i=1;i<=data->noelements;i++)
    data->material[i] = 0;

  for(i=1;i<=data->noelements;i++)
    data->elementtypes[i] = 0;

  data->x = Rvector(1,data->noknots);
  data->y = Rvector(1,data->noknots);
  data->z = Rvector(1,data->noknots);
  for(i=1;i<=data->noknots;i++)
    data->x[i] = data->y[i] = data->z[i] = 0.0;

  data->created = TRUE;

#if DEBUG
  printf("Allocated for %d %d-node elements resulting to %d nodes\n",
	 data->noelements,data->maxnodes,data->noknots);
#endif
}


static void MovePointCircle(Real *lim,int points,Real *coords,
			    Real x,Real y,Real *dx,Real *dy)
{
  int i;
  Real x0,y0,r,r1,r2,p;

  r2 = fabs(lim[1]-lim[0]);
  if(r2 > fabs(lim[2]-lim[1]))
    r2 = fabs(lim[2]-lim[1]);

  for(i=0;i<points/2;i++) {
    x0 = x-coords[2*i];
    r1 = fabs(coords[2*i+1]);

    if(fabs(x0) >= r2) continue;
    y0 = y-(lim[1]+coords[2*i+1]);
    if(y0 < 0 && lim[0] > lim[1]) continue;
    if(y0 > 0 && lim[2] < lim[1]) continue;
    if(fabs(y0) >= r2) continue;
    r = sqrt(x0*x0+y0*y0);
    if(r < 1.0e-50) continue;

    if(fabs(x0) > fabs(y0)) {
      p = fabs(x0)/r - 1.0;
      if(fabs(x0) <= r1) {
	*dx += p*x0;
	*dy += p*y0;
      }
      else if(fabs(x0) <= r2) {
	*dx += p*x0*(r2-fabs(x0))/(r2-r1);
	*dy += p*y0*(r2-fabs(x0))/(r2-r1);
      } 
    } 
    else {
      p = fabs(y0)/r - 1.0;
      if(fabs(y0) <= r1) {
	*dx += p*x0;
	*dy += p*y0;
      }
      else if(fabs(y0) <= r2) {
	*dx += p*x0*(r2-fabs(y0))/(r2-r1);
	*dy += p*y0*(r2-fabs(y0))/(r2-r1);
      }       
    }
  }
}



static void MovePointLinear(Real *lim,int points,Real *coords,
			    Real x,Real y,Real *dx,Real *dy)
{
  static int i=0;
  Real c,d;

  if(y > lim[0]  &&  y < lim[2]) {

    
    if(x <= coords[0]) {
      d = coords[1];
    }
    else if(x >= coords[points-2]) {
      d = coords[points-1];
    }
    else {
      i = 1;
      while(x > coords[2*i] && i < points/2-1) i++;
      c = (coords[2*i+1]-coords[2*i-1])/(coords[2*i]-coords[2*i-2]);
      d = coords[2*i-1] + c*(x-coords[2*i-2]);
    }

    if(y < lim[1])
      *dy += d*(y-lim[0])/(lim[1]-lim[0]);
    else
      *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
  }
}


static void MovePointAngle(Real *lim,int points,Real *coords,
			   Real x,Real y,Real *dx,Real *dz)
{
  static int i=0;
  Real x1,z1,degs;

  degs = FM_PI/180.0;
  x1 = z1 = 0.0;

  if(y > lim[0]  &&  y < lim[2]) {
    if(x <= coords[0]) {
      x1 = x;
    }
    else {
      i = 1;
      while(x > coords[2*i] && i <= points/2-1) {
	x1 = x1 + cos(degs*coords[2*i-1])*(coords[2*i]-coords[2*i-2]);
	z1 = z1 + sin(degs*coords[2*i-1])*(coords[2*i]-coords[2*i-2]);	
	i++;
      }
      x1 = x1 + cos(degs*coords[2*i-1])*(x-coords[2*i-2]);
      z1 = z1 + sin(degs*coords[2*i-1])*(x-coords[2*i-2]);
    }

    if(y < lim[1]) {
      *dx += (x1-x)*(y-lim[0])/(lim[1]-lim[0]);
      *dz += z1*(y-lim[0])/(lim[1]-lim[0]);
    }
    else {
      *dx += (x1-x)*(lim[2]-y)/(lim[2]-lim[1]);      
      *dz += z1*(lim[2]-y)/(lim[2]-lim[1]);      
    }
  }
}


static void MovePointSinus(Real *lim,int points,Real *coords,
			   Real x,Real y,Real *dx,Real *dy)
{
  Real c,d;

  if(y > lim[0]  &&  y < lim[2]) {

    if(x <= coords[0]) {
      d = 0.0;
    }
    else if(x >= coords[1]) {
      d = coords[3]*sin(coords[2]*2.*FM_PI);
    }
    else {
      c = coords[2]*2.*FM_PI/(coords[1]-coords[0]);
      d = coords[3]*sin(c*(x-coords[0]));
    }

    if(y < lim[1])
      *dy += d*(y-lim[0])/(lim[1]-lim[0]);
    else
      *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
  }
}



static void MovePointPowerSeries(Real *lim,int points,Real *coords,
				 Real x,Real y,Real *dx,Real *dy)
{
  int i,n;
  Real d,t,u;

  if(y > lim[0]  &&  y < lim[2]) {
    t = x;
    if(coords[1] > coords[0]) {
      if(t<coords[0]) t = coords[0];
      if(t>coords[1]) t = coords[1];
    }
    else {
      if(t>coords[0]) t = coords[0];
      if(t<coords[1]) t = coords[1];
    }      
    
    n = points-2;
    u = (t - coords[0])/(coords[1]-coords[0]);
    
    d = 0.0;
    for(i=0;i<n;i++) { 
      d += coords[i+2] * pow(u,i);
    }
    
    if(y < lim[1]) {
      *dy += d*(y-lim[0])/(lim[1]-lim[0]);
    }
    else {
      *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
    }
  }
}


static void MovePointPowerSeries2(Real *lim,int points,Real *coords,
				 Real x,Real y,Real *dx,Real *dy)
{
  int i,j,n;
  Real d,e,t,u,h;

  if(y > lim[0]  &&  y < lim[2]) {
    t = x;
    if(coords[1] > coords[0]) {
      if(t<coords[0]) t = coords[0];
      if(t>coords[1]) t = coords[1];
    }
    else {
      if(t>coords[0]) t = coords[0];
      if(t<coords[1]) t = coords[1];
    }      
    
    n = points-2;
    u = (t - coords[0])/(coords[1]-coords[0]);
    
    d = 0.0;
    
    d = coords[2];
    if(n>=1) d += u * coords[3];

    for(i=2;i<n;i++) { 
      h = 1.0/(i-1);
      e = 1.0;
      for(j=0;j<i;j++)
	e *= (u-j*h);
      d += coords[i+2] * e;
    }
    
    if(y < lim[1]) {
      *dy += d*(y-lim[0])/(lim[1]-lim[0]);
    }
    else {
      *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
    }
  }
}


  
static void MovePointPower(Real *lim,int points,Real *coords,
			   Real x,Real y,Real *dx,Real *dy)
{
  static int i=0;
  Real c,d;
  
  if(y > lim[0]  &&  y < lim[2]) {
    if(x <= coords[0]) {
      d = coords[1];
    }
    else if(x >= coords[points-2]) {
      d = coords[points-1];
    }
    else {
      i = 1;
      while(x > coords[3*i] && i < points/3-1) i++;
      c = (coords[3*i+1]-coords[3*i-2])/pow((coords[3*i]-coords[3*i-3]),coords[3*i-1]);
      d = coords[3*i-2] + c*pow((x-coords[3*i-3]),coords[3*i-1]);
    }
    
    if(y < lim[1])
      *dy += d*(y-lim[0])/(lim[1]-lim[0]);
    else
      *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
  }
}

/* Creates airfoil shapes */
static void MovePointNACAairfoil(Real *lim,int points,Real *coords,
				 Real x,Real y,Real *dx,Real *dy)
{
  Real p,d,t,u;

  if(y < lim[0]  ||  y > lim[2]) return;
  if(x < coords[0]  ||  x > coords[1]) return;

  if(0) {
    printf("x=%.3e y=%.3e lim0=%.3e lim2=%.3e\n",x,y,lim[0],lim[2]);
    printf("naca: %.3e %.3e %.3e\n",coords[0],coords[1],coords[2]);
  }

  t = x;
  if(coords[1] > coords[0]) {
    if(t<coords[0]) t = coords[0];
    if(t>coords[1]) t = coords[1];
  }
  else {
    if(t>coords[0]) t = coords[0];
    if(t<coords[1]) t = coords[1];
  }      
    
  u = (t - coords[0])/(coords[1]-coords[0]);    
  p = 0.2969*sqrt(u) - 0.1260*u - 0.3537*u*u + 0.2843*u*u*u - 0.1015*u*u*u*u;
  
  d = coords[2] * (coords[1]-coords[0]) * p / 0.2;
  
  if(y < lim[1]) {
    *dy += d*(y-lim[0])/(lim[1]-lim[0]);
  }
  else {
    *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
  }


  if(0) printf("d=%.3e p=%.3e u=%.3e dy=%.3e\n",d,p,u,*dy);
}



static void MovePointArc(Real *lim,int points,Real *coords,
			 Real x,Real y,Real *dx,Real *dy)
{
  static int i=0;
  Real sx,sy,ss,r,rat,d,x0,y0;

  if(y > lim[0]  &&  y < lim[2]) {
    if(x <= coords[0]) {
      d = coords[1];
    }
    else if(x >= coords[points-2]) {
      d = coords[points-1];
    }
    else {
      i = 1;
      while(x > coords[3*i] && i < points/3-1) i++;
      sx = 0.5*(coords[3*i]-coords[3*i-3]);
      sy = 0.5*(coords[3*i+1]-coords[3*i-2]);
      r = coords[3*i-1];
      ss = sx*sx+sy*sy;
      rat = sqrt(1.0/ss-1.0/(r*r))*r; 
      x0 = coords[3*i-3] + sx - sy * rat;
      y0 = coords[3*i-2] + sy + sx * rat; 
      d  = y0-sqrt(r*r-(x-x0)*(x-x0))*r/fabs(r);
    }

    if(y < lim[1])
      *dy += d*(y-lim[0])/(lim[1]-lim[0]);
    else
      *dy += d*(lim[2]-y)/(lim[2]-lim[1]);      
  }
}



void CreateKnots(struct GridType *grid,struct CellType *cell,
		 struct FemType *data,int noknots,int info)
/* Saves information concerning the knots to a special structure to avoid 
   repetitious calculations. This should be used unless there is a severe
   lack of memory. GridType includes only rectangular 2D elements.
   */ 
{
  Real globalcoord[DIM*MAXNODESD2];
  Real maplim[3*MAXMAPPINGS];
  int material,nonodes,elemind,elemtype;
  int mode,level,maplevel,dim;
  int celli,cellj,i,j,k,l,ind[MAXNODESD2];
  Real x,y,dx,dy,dz,size,minsize,maxsize;

  InitializeKnots(data);

  dim = data->dim = grid->dimension;
  nonodes = grid->nonodes;
  data->maxnodes = grid->nonodes;
  data->nocells = grid->nocells;
  data->noelements = grid->noelements;
  data->coordsystem = grid->coordsystem;
  data->numbering = grid->numbering;
  data->indexwidth = grid->maxwidth;
  data->noknots = MAX(noknots,grid->noknots);
  data->noboundaries = grid->noboundaries;
  
  AllocateKnots(data);
  minsize = 1.0e20;

  if(dim == 1) 
    elemtype = grid->nonodes + 200;
  else 
    elemtype = grid->nonodes + 400;
    
  for(i=1;i<=data->noelements;i++) 
    data->elementtypes[i] = elemtype;

  /* This numbers the elements the same way the knots are numbered. */
  for(cellj=1;cellj<= grid->ycells ;cellj++) {        /* cells direction up     */
    for(j=1; j<=grid->yelems[cellj]; j++)             /* lines inside cells     */
      for(celli=1;celli<= grid->xcells; celli++)      /* cells direction right  */
	if(k=grid->numbered[cellj][celli]) {
	  material = cell[k].material;
	  for(i=1; i<=grid->xelems[celli]; i++) {

	    elemind = GetElementCoordinates(&(cell)[k],i,j,globalcoord,ind);

	    if(data->noknots == grid->noknots) {
	      for(l=0;l<nonodes;l++) {
		data->topology[elemind][l] = ind[l];
		data->x[ind[l]] = globalcoord[l];
		data->y[ind[l]] = globalcoord[l+nonodes];
	      }
	    data->material[elemind] = material;

	    }
	  }	  
	}
  }    

  /* Map the knots as defined in structures grid */
  for(k=0;k<grid->mappings;k++) {
    j = grid->mappingline[k];
    if(grid->mappingtype[k] > 0) 
      maplim[3*k+1] = grid->y[j];
    else if(grid->mappingtype[k] < 0)
      maplim[3*k+1] = grid->x[j];
    else 
      continue;
    maplim[3*k]   = maplim[3*k+1] - grid->mappinglimits[2*k];
    maplim[3*k+2] = maplim[3*k+1] + grid->mappinglimits[2*k+1];
  }

  mode = 0;
  if(grid->mappings) 
    
    for(level=0;level<10;level++) {
      maplevel = FALSE;    
      for(k=0;k<grid->mappings;k++) 
	if(abs(grid->mappingtype[k]/10) == level)
	  maplevel = TRUE;
      if(maplevel == FALSE) continue;

      if(level >= 5) data->dim = 3;

      for(i=1;i<=data->noknots;i++) {
	x = data->x[i];
	y = data->y[i];
	dx = 0.0;
	dy = 0.0;
	dz = 0.0;

	for(k=0;k<grid->mappings;k++) {
	  mode = grid->mappingtype[k]%10;
	  switch (mode) {
	  case 1:
	    MovePointLinear(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			    x,y,&dx,&dy);
	    break;
	  case 2:
	    MovePointPower(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   x,y,&dx,&dy);
	    break;
	  case 3:
	    MovePointArc(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			 x,y,&dx,&dy);
	    break;
	  case 4:
	    MovePointCircle(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			    x,y,&dx,&dy);
	    break;
	  case 5:
	    MovePointSinus(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   x,y,&dx,&dy);
	    break;
	  case 6:
	    MovePointPowerSeries(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   x,y,&dx,&dy);
	    break;
	  case 7:
	    MovePointPowerSeries2(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
				  x,y,&dx,&dy);
	    break;
	  case 8:
	    MovePointAngle(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   x,y,&dx,&dz);
	    break;	
	  case 9:
	    MovePointNACAairfoil(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
				 x,y,&dx,&dy);
	    break;	


	  case -1:
	    MovePointLinear(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			    y,x,&dy,&dx);
	    break;
	  case -2:
	    MovePointPower(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   y,x,&dy,&dx);
	    break;
	  case -3:
	    MovePointArc(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			 y,x,&dy,&dx);
	    break;
	  case -4:
	    MovePointCircle(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			    y,x,&dy,&dx);
	    break;
	  case -5:
	    MovePointSinus(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   y,x,&dy,&dx);
	    break;
	  case -6:
	    MovePointPowerSeries(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
				 y,x,&dy,&dx);
	    break;
	  case -7:
	    MovePointPowerSeries2(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
				  y,x,&dy,&dx);
	    break;
	  case -8:
	    MovePointAngle(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   y,x,&dy,&dz);
	    break;
	  case -9:
	    MovePointNACAairfoil(&maplim[3*k],grid->mappingpoints[k],grid->mappingparams[k],
			   y,x,&dy,&dx);
	    break;


	  }  
	}

	if(mode == 8 || mode == -8) {
	  data->x[i] += dx;
	  data->y[i] += dy;	  
	  data->z[i] += dz;
	}
	else if(level >= 5) {
	  data->z[i] += dx + dy;
	}
	else {
	  data->x[i] += dx; 
	  data->y[i] += dy;
	}
      }
    }

  minsize = 1.0e20;
  maxsize = 0.0;

  for(i=1;i<=data->noelements;i++) {
    GetElementInfo(i,data,globalcoord,ind,&material);

    dx = globalcoord[0]-globalcoord[1];
    dy = globalcoord[nonodes]-globalcoord[nonodes+1];
    size = dx*dx+dy*dy;
    if(size < minsize) minsize = size;
    if(size > maxsize) maxsize = size;
    dx = globalcoord[0]-globalcoord[nonodes-1];
    dy = globalcoord[nonodes]-globalcoord[2*nonodes-1];
    size = dx*dx+dy*dy;
    if(size < minsize) minsize = size;
    if(size > maxsize) maxsize = size;
  }  

  data->maxsize = sqrt(maxsize);
  data->minsize = sqrt(minsize);

  if(info) printf("Maximum elementsize is %.3e and minimum %.3e.\n",
		  data->maxsize,data->minsize);
}




int CreateVariable(struct FemType *data,int variable,int unknowns,
		   Real value,const char *dofname,int eorder)
/* Create variables for the given data structure */
{
  int i,info=FALSE;
  int timesteps;

  if(variable == 0) return(0);

  if(data->created == FALSE) {
    if(info) printf("CreateVariable: Knots must first be created!\n");
    return(1);
  }
  timesteps = data->timesteps;
  if(timesteps < 1) timesteps = 1;

  if(data->edofs[variable] == 0) {

    data->variables += 1;
    data->edofs[variable] = unknowns;
    data->alldofs[variable] = unknowns * data->noknots;
    data->bandwidth[variable] = unknowns * data->indexwidth;
    data->dofs[variable]  = Rvector(1,timesteps * data->alldofs[variable]);
    if(info) printf("Created variable %s with %d dofs.\n",
		    dofname,data->alldofs[variable]);
    for(i=1;i<=data->alldofs[variable]*timesteps;i++)
      data->dofs[variable][i] = value;  
    data->iterdofs[variable] = 1;
  }
  else if (data->edofs[variable] == unknowns) {
    if(info) printf("CreateVariable: Variable %d exists with correct number of dofs!\n",
		    variable);  
  } 
  else {
    if(info) printf("CreateVariable: Variable %d exists with wrong number of dofs!\n",
		    variable);
    return(2);
  }  


  if(eorder) {
    if (data->eorder[variable] == FALSE) {
      data->eorder[variable] = TRUE;
      data->order[variable] = Ivector(1,data->alldofs[variable]);
      for(i=1;i<=data->alldofs[variable];i++)
	data->order[variable][i] = i;
    }
    if(info) printf("Created index for variable %s.\n",dofname);
  }

  strcpy(data->dofname[variable],dofname);

  return(0);
}



void DestroyKnots(struct FemType *data)
{
  int i;

  if(!data->created) return;
  data->created = FALSE;

  for(i=0;i<MAXDOFS;i++) 
    if(data->edofs[i] != 0) {
      if(data->edofs[i] > 0) 
	free_Rvector(data->dofs[i],1,data->alldofs[i]);
      data->edofs[i] = 0;
    }

  free_Imatrix(data->topology,1,data->noelements,0,data->maxnodes-1);
  free_Ivector(data->material,1,data->noelements);
  free_Ivector(data->elementtypes,1,data->noelements);

  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);

  data->noknots = 0;
  data->noelements = 0;
  data->maxnodes = 0;

  if(data->nocorners > 0)
    free_Ivector(data->corners,1,2*data->nocorners);
}



int CreateBoundary(struct CellType *cell,struct FemType *data,
		   struct BoundaryType *bound,int material1,int material2,
		   int solidmat,int boundarytype,int info)
/* This subroutine makes a boundary which includes all sides that separate
   two materials that fulfill the conditions in the function call. If both
   materials are positive only the sides for which both of the materials
   coincide are accepted. In other cases the negative argument tells which
   conditions the positive argument should fulfill. Note that on a boundary
   where knots are created only for the other material, this material
   should be the latter one in the function call (material). The physical
   properties (emissivity) of the materials are taken from the one given
   by the flag 'solidmat'.
   */
{
  int i,side,more,elem,elemind[2],nosides,no,times;
  int sidemat,thismat,size,setpoint,dimsides,cellside;

  if(data->dim == 1) 
    dimsides = 2;
  else 
    dimsides = 4;

  if(bound->created == TRUE) {
    if(info) printf("CreateBoundary: You tried to recreate the boundary!\n");
    return(1);
  }
  if(!data->created) {
    if(info) printf("CreateBoundary: You tried to create a boundary before the knots were made.");
    return(2);
  }
  if(material1 < 0  &&  material2 < 0) {
    if(info) printf("CreateBoundary: the material arguments are both negative");
    return(3);
  }

  times = 0;

  bound->created = FALSE;
  bound->nosides = 0;
  if(solidmat >= 2) solidmat -= 2;

startpoint:

  /* Go through all elements which have a boundary with the given material, but 
     are not themselves of that material. First only calculate their amount, then
     allocate space and tabulate them. */ 
  nosides = 0;


  for(no=1; no <= data->nocells; no++) 
    for(side=0; side < dimsides; side++) { 

      if(data->dim == 1) 
	cellside = 3-2*side;
      else
	cellside = side;

      setpoint = FALSE;
      sidemat = cell[no].boundary[cellside];
      thismat = cell[no].material;

      /* The free boundary conditions are not allowed if the negative
	 keywords are used. */

      /* Either material must be the one defined. */
      if( material1 >= 0  &&  material1 != sidemat) continue;
      if( material2 >= 0  &&  material2 != thismat) continue;
#if 0
      printf("mat=[%d %d]  sidemat=%d  thismat=%d  side=%d\n",
	     material1,material2,sidemat,thismat,side);
#endif

      if( material2 == -((side+2)%4+1) &&  sidemat ==  material1 &&  
	  sidemat != thismat) setpoint = TRUE;
      if( material1 == -(side+1)       &&  thismat ==  material2 &&
	  sidemat != thismat) setpoint = TRUE;

      if( material1 == MAT_BIGGER    &&  sidemat >  material2 ) setpoint = TRUE;
      if( material1 == MAT_SMALLER   &&  sidemat <  material2 ) setpoint = TRUE;
      if( material1 == MAT_ANYTHING  &&  sidemat != material2 ) setpoint = TRUE;
      if( material2 == MAT_BIGGER    &&  thismat >  material1 ) setpoint = TRUE;
      if( material2 == MAT_SMALLER   &&  thismat <  material1 ) setpoint = TRUE;
      if( material2 == MAT_ANYTHING  &&  thismat != material1 ) setpoint = TRUE;
      if( sidemat == material1   &&  thismat == material2 )     setpoint = TRUE; 

      if(setpoint == TRUE) {
#if 0
	printf("going through boundary %d vs. %d in cell %d\n",material1,material2,no);
#endif
	elem = 0; 
	do  { 
	  elem++;
	  nosides++;
	  more = GetSideInfo(cell,no,side,elem,elemind);

#if 0
	  printf("elem=%d nosides=%d no=%d side=%d elemind=%d %d\n",
		 elem,nosides, no, side, elemind[0], elemind[1]);
#endif	
  
	  /* In the second round the values are tabulated. */
	  if(times) {
	    /* It is assumed that the material pointed by solidmat
	       determines the surface properties. */

	    bound->parent[nosides]  = elemind[0];
	    bound->parent2[nosides] = elemind[1];

	    bound->side[nosides]  = side;
	    bound->side2[nosides] = (side+dimsides/2)%dimsides;

	    bound->types[nosides] = boundarytype;

	    /* The direction of the surface normal must be included */
	    if(solidmat==FIRST) {
	      bound->material[nosides] = sidemat;
	      bound->normal[nosides] = 1;
	    }
	    if(solidmat==SECOND){
	      bound->material[nosides]  = thismat;
	      bound->normal[nosides] = -1;
	    }	    
	  }
	  
	  
	} while(more); 
      } 
    }  

  if(nosides == 0) {
    if(info) printf("No boundary between materials %d and %d exists.\n",
		    material1,material2); 
    return(0);
  }

  if(times == 0) {
    times++;

    /* Allocate space. This has sometimes led to strange errors. 
       The allocation takes place only in the first loop. */
    
    bound->created = TRUE;
    bound->nosides = size = nosides;
    bound->coordsystem = data->coordsystem;
    bound->types = Ivector(1,nosides);
    bound->side = Ivector(1,nosides);
    bound->side2 = Ivector(1,nosides);
    bound->material = Ivector(1,nosides);    
    bound->parent = Ivector(1,nosides);
    bound->parent2 = Ivector(1,nosides);
    bound->normal = Ivector(1,nosides);

    bound->echain = FALSE;
    bound->ediscont = FALSE;

    goto startpoint;
  }
  
  if(info) printf("%d element sides between materials %d and %d were located to type %d.\n",
		  nosides,material1,material2,boundarytype); 

  return(0);
}


int AllocateBoundary(struct BoundaryType *bound,int size)
{
  int i;

  if(bound->created == TRUE) {
    printf("AllocateBoundary: You tried to recreate the boundary!\n");
    return(1);
  }
    
  bound->created = TRUE;
  bound->nosides = size;
  bound->echain = FALSE;
  bound->ediscont = FALSE;

  bound->material = Ivector(1,size);    
  bound->side = Ivector(1,size);
  bound->side2 = Ivector(1,size);
  bound->parent = Ivector(1,size);
  bound->parent2 = Ivector(1,size);
  bound->types = Ivector(1,size);
  bound->normal = Ivector(1,size);

  for(i=1;i<=size;i++) {
    bound->material[i] = 0;
    bound->side[i] = 0;
    bound->side2[i] = 0;
    bound->parent[i] = 0;
    bound->parent2[i] = 0;
    bound->types[i] = 0;
    bound->normal[i] = 1;
  }

  return(0);
}



int DestroyBoundary(struct BoundaryType *bound)
/* Destroys boundaries of various types. */
{
  int i,nosides;

  if(!bound->created) {
    return(1);
  }
  nosides = bound->nosides;
  if(!nosides) {
    bound->created = FALSE;
    return(2);
  }

  free_Ivector(bound->material,1,nosides);
  free_Ivector(bound->side,1,nosides);
  free_Ivector(bound->side2,1,nosides);
  free_Ivector(bound->parent,1,nosides);
  free_Ivector(bound->parent2,1,nosides);
  free_Ivector(bound->types,1,nosides);
  free_Ivector(bound->normal,1,nosides);

  bound->nosides = 0;
  bound->created = FALSE;

#if DEBUG
  printf("%d element sides were destroyed.\n",nosides); 
#endif
  return(0);
}



int CreatePoints(struct CellType *cell,struct FemType *data,
		 struct BoundaryType *bound,
		 int param1,int param2,int pointmode,int pointtype,int info)
{
  int size,i,no,corner,times,elem,node;

  bound->created = FALSE;
  bound->nosides = 0;    
  times = 0;
     
omstart:
  i = 0;

  /* Create nodes that are divided by the two materials specified */
  if(pointmode == 4) {
    for(no=1; no <= data->nocells; no++) 
      if(cell[no].material == param2) {
	
	for(corner=0; corner < 4; corner++) 
	  if(cell[no].boundary[4+corner] == param1) {
	    
	    i++;

	    if(times) {
	      bound->material[i] = param2;
	      bound->types[i] = pointtype;
	      bound->side[i] = 4 + corner;
	      
	      if(corner == BOTLEFT) 
		elem = GetElementIndex(&cell[no],1,1);
	      else if(corner == BOTRIGHT)	 
		elem = GetElementIndex(&cell[no],cell[no].xelem,1);
	      else if(corner == TOPRIGHT)
		elem = GetElementIndex(&cell[no],cell[no].xelem,cell[no].yelem);
	      else if(corner == TOPLEFT)
		elem = GetElementIndex(&cell[no],1,cell[no].yelem);
	      
	      bound->parent[i] = elem;
	    }
	  }
      }
  }  

  if(pointmode == 5) {
    corner = -1;
    for(no=1; no <= data->nocells && corner <0; no++) {
      if(cell[no].xind-1 == param1 && cell[no].yind-1 == param2) 
	corner = BOTLEFT;
      else if(cell[no].xind == param1 && cell[no].yind-1 == param2) 
	corner = BOTRIGHT;
      else if(cell[no].xind == param1 && cell[no].yind == param2) 
	corner = TOPRIGHT;
      else if(cell[no].xind-1 == param1 && cell[no].yind == param2) 
	corner = TOPLEFT;
    }
    if(corner >= 0) {
      i++;
      no--;
    }

    if(times) {
      bound->types[i] = pointtype;
      bound->side[i] = 4 + corner;

      if(corner == BOTLEFT) 
	elem = GetElementIndex(&cell[no],1,1);
      else if(corner == BOTRIGHT)	 
	elem = GetElementIndex(&cell[no],cell[no].xelem,1);
      else if(corner == TOPRIGHT)
	elem = GetElementIndex(&cell[no],cell[no].xelem,cell[no].yelem);
      else if(corner == TOPLEFT)
	elem = GetElementIndex(&cell[no],1,cell[no].yelem);      

      bound->parent[i] = elem;
      node = data->topology[elem][corner];
      if(info) printf("Found node %d at (%.3lg, %.3lg)\n",node,data->x[node],data->y[node]);
    }
  }  

  size = i;
  if(times == 0 && size > 0) {
    AllocateBoundary(bound,size);
    times = 1;
    goto omstart;
  }

  if(info) printf("Created %d new points of type %d in the corner of materials %d and %d.\n",
		  size,pointtype,param1,param2);
  
  return(FALSE);
}



int CreateNewNodes(struct FemType *data,int *order,int material,int newknots)
{
  int i,j,k,l,lmax,ind;
  int newsize,noknots,nonodes;
  int *neworder;
  Real *newx=NULL,*newy=NULL,*newz=NULL;
  Real *newdofs[MAXDOFS];

  noknots = data->noknots;

  printf("Creating %d new nodes for discontinuous boundary.\n",newknots);

  /* Allocate for the new nodes */
  newsize = noknots+newknots;
  newx = Rvector(1,newsize);
  newy = Rvector(1,newsize);
  newz = Rvector(1,newsize);

  neworder = Ivector(1,newsize);

  for(j=1;j<MAXDOFS;j++)
    if(data->edofs[j])
      newdofs[j] = Rvector(1,data->edofs[j]*newsize);
  
  /* Set the new coordinates and dofs */
  j = 0;
  for(i=1;i<=noknots;i++) {
    j++;
    neworder[j] = i;
    newx[j] = data->x[i];
    newy[j] = data->y[i];
    newz[j] = data->z[i];

    for(k=1;k<MAXDOFS;k++) {
      if(lmax = data->edofs[k])
	for(l=1;l<=lmax;l++) 
	  newdofs[k][lmax*(j-1)+l] = data->dofs[k][lmax*(i-1)+l];
    }
    
    if(order[i] < 0) {
      j++;
      neworder[j] = -i;
      newx[j] = data->x[i];
      newy[j] = data->y[i];
      newz[j] = data->z[i];

      for(k=1;k<MAXDOFS;k++) {
	if(lmax = data->edofs[k])
	  for(l=1;l<=lmax;l++) 
	    newdofs[k][lmax*(j-1)+l] = data->dofs[k][lmax*(i-1)+l];
      }
    }
  }
  
  /* Find the old index corresponding to the new one. */
  for(i=1;i<=newsize;i++) 
    if(neworder[i] > 0) {
      if(order[neworder[i]] < 0)  
	order[neworder[i]] = -i; 
      else 
	order[neworder[i]] = i; 
    }
  
  /* Set the new element topology */
  for(i=1;i<=data->noelements;i++) {
    nonodes = data->elementtypes[i]%100;
    for(j=0;j<nonodes;j++) {
      ind = data->topology[i][j];
      if(data->material[i] == material && order[ind] < 0)
	data->topology[i][j] = abs(order[ind])+1;
      else
	data->topology[i][j] = abs(order[ind]);
    }
  }    

  /* Destroy old vectors and set the pointers to the new vectors. */
  free_Rvector(data->x,1,noknots);
  free_Rvector(data->y,1,noknots);
  free_Rvector(data->z,1,noknots);

  for(k=1;k<MAXDOFS;k++)
    if(data->edofs[k]) free_Rvector(data->dofs[k],1,noknots);
  
  data->noknots = newsize;
  data->x = newx;
  data->y = newy;
  data->z = newz;

  for(k=1;k<MAXDOFS;k++) {
    if(data->edofs[k]) {
      data->dofs[k] = newdofs[k];
      data->alldofs[k] = data->edofs[k] * data->noknots;
    }  
  }
  
  return(0);
}


int SetDiscontinuousBoundary(struct FemType *data,struct BoundaryType *bound,
			     int boundtype,int endnodes,int info)
/* Create secondary points for a given boundary.
   This feature is handy when one wants to solve problems with discontinuous
   field variables.
   */
{
  int i,j,bc,ind,sideind[MAXNODESD1];
  int side,parent,newnodes,doublesides,maxtype,newbc;
  int newsuccess,noelements,nonodes,sideelemtype,sidenodes,disconttype;
  int *order=NULL;
  int mat1,mat2,par1,par2,mat1old,mat2old,material;
  static int hitsexist=FALSE,hitslength,*hits=NULL;


  if(boundtype < 0) {
    newbc = TRUE;
    boundtype = -boundtype;
  }
  else {
    newbc = FALSE;
  }

  mat1old = mat2old = 0;
  doublesides = 0;

  /* Compute the number of duplicate boundary elements */
  for(bc=0;bc<MAXBOUNDARIES;bc++) {

    if(bound[bc].created == FALSE) continue;
    if(!bound->nosides) continue;
    
    for(i=1;i<=bound[bc].nosides;i++) {
      if(bound[bc].types[i] == boundtype) {
	par1 = bound[bc].parent[i];
	par2 = bound[bc].parent2[i];
	if(par1 && par2) {
	  doublesides++;
	  mat1 = data->material[par1];
	  mat2 = data->material[par2];
	  if(!mat1old) mat1old = mat1;
	  else if(mat1old != mat1) mat1old = -mat1;
	  if(!mat2old) mat2old = mat2;
	  else if(mat2old != mat2) mat2old = -mat2;
	}
      }
    }
  }  
   
  if(!doublesides) return(1);
    
  if( mat1old > 0 && mat2old > 0 ) 
    material = MIN( mat1old, mat2old );
  else if(mat1old > 0) 
    material = mat1old;
  else if(mat2old > 0) 
    material = mat2old;
  else {
    printf("SetDiscontinuousBoundary: impossible to make the boundary of several materials\n");
    return(2);
  }

  if(info) {
    printf("Creating discontinuous boundary between materials %d and %d\n",mat1old,mat2old);
    printf("New set of nodes will be created for material %d\n",material);
  }


  noelements = data->noelements;    
  order = Ivector(1,data->noknots);
  for(i=1;i<=data->noknots;i++)
    order[i] = i;
    
  /* Compute the endnodes by the fact that they have different number of
     hits */
  if(endnodes == 1) {
    if(!hitsexist) {
      hitslength = 1.1*data->noknots;
      hits = Ivector(1,hitslength);
      hitsexist = TRUE;
    }
    else if(hitslength <= data->noknots) {
      free_Ivector(hits,1,hitslength);
      hitslength = 1.1*data->noknots;
      hits = Ivector(1,hitslength);
    }
    
    for(i=1;i<=data->noknots;i++)
      hits[i] = 0;
    
    for(j=1;j<=noelements;j++) {
      nonodes = data->elementtypes[j]%100;
      for(i=0;i<nonodes;i++)
	hits[data->topology[j][i]] += 1;
    }
  }
    
  /* If requested create a secondary boundary at the other side */
  if(newbc) {
    maxtype = 0;
    for(bc=0;bc<MAXBOUNDARIES;bc++) {
      for(i=1;i<=bound[bc].nosides;i++) {
	maxtype = MAX(maxtype, bound[bc].types[i]);
	if(bound[bc].ediscont) maxtype = MAX(maxtype, bound[bc].discont[i]);
      }
    }
    disconttype = maxtype + 1;
    if(info) printf("Type of the other side of discontinuous boundary set to %d.\n",disconttype);
  }
  else {
    disconttype = boundtype;
  }


  
  /* Find the number of new nodes */
  newnodes = 0;

  for(bc=0;bc<MAXBOUNDARIES;bc++) {

    for(i=1;i<=bound[bc].nosides;i++) {
      
      if(bound[bc].types[i] != boundtype) continue;
      
      if(!bound[bc].ediscont) {
	bound[bc].discont = Ivector(1,bound[bc].nosides);
	for(j=1;j<=bound[bc].nosides;j++) 
	  bound[bc].discont[j] = 0;
	bound[bc].ediscont = TRUE;
      }

      parent = bound[bc].parent2[i];
      if(parent == 0) continue;
      side = bound[bc].side2[i];
      GetElementSide(parent,side,1,data,sideind,&sideelemtype);
      sidenodes = sideelemtype%100;
      
      bound[bc].discont[i] = disconttype;
      
      for(j=0;j<sidenodes;j++) {      
	ind = abs(sideind[j]);
	
	if(endnodes == 2) {
	  if(order[ind] > 0) {
	    newnodes++;
	    order[ind] = -newnodes;
	  }
	}
	else if(endnodes == 0) {
	  if(order[ind] > 0)
	    order[ind] = 0;
	  else if(order[ind] == 0) {
	    newnodes++;
	    order[ind] = -newnodes;	
	  }
	}
	else if(endnodes == 1) {
	  if(order[ind] > 0) {
	    if(hits[ind] < 4) {
	      newnodes++;
	      order[ind] = -newnodes;	    
	    }
	    else
	      order[ind] = 0;
	  }
	  else if(order[ind] == 0) {
	    newnodes++;
	    order[ind] = -newnodes;	
	  }
	}
	
      }      
    }
    
    if(endnodes == 0 || endnodes == 1) {
      for(i=1;i<=data->noknots;i++) 
	if(order[i] == 0) 
	  order[i] = i;
    }
  }

  if(newnodes == 0) return(3);
    
  newsuccess = CreateNewNodes(data,order,material,newnodes);
  return(newsuccess);
}



int FindPeriodicBoundary(struct FemType *data,struct BoundaryType *bound,
                         int boundary1,int boundary2,int info)
/* Create periodic boundary conditions for a given boundary pair
   boundary1, boundary2. 
   */
{
  int i,j,k,l;
  int parent,elemtype;
  int minp[2],maxp[2],bounds[2],dp[2],sumsides[2];

  if(bound->created == FALSE) {
    printf("SetDiscontinuousBoundary: The boundary does not exist!\n");
    return(1);
  }
  if(!bound->nosides) return(0);


  bounds[0] = boundary1;
  bounds[1] = boundary2;
  minp[0] = minp[1] = data->noknots+1;
  maxp[0] = maxp[1] = 0;

  sumsides[0] = sumsides[1] = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound->created == FALSE) continue;

    for(i=1; i <= bound[j].nosides; i++) {

      for(k=0;k<=1;k++) {
	if(bound[j].types[i] == bounds[k]) {

	  sumsides[k] += 1;
	  if(bound[j].parent[i] > maxp[k]) maxp[k] = bound[j].parent[i];
	  if(bound[j].parent[i] < minp[k]) minp[k] = bound[j].parent[i];
	}
      } 
    }
  }

  for (k=0;k<=1;k++) {
    dp[k] = maxp[k] - minp[k];
    if(info) printf("Parents of boundary %d are on the interval [%d, %d]\n",
		    bounds[k],minp[k],maxp[k]);
  }  

  if(dp[0] != dp[1] || sumsides[0] != sumsides[1]) {
    printf("FindPeriodicBoundary: The simple scheme cannot handle these boundaries!\n");
    printf("dp=[%d %d]  sumsides=[%d %d]\n",dp[0],dp[1],sumsides[0],sumsides[1]);
    return(1);
  }

  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound->created == FALSE) continue;

    for(i=1; i <= bound[j].nosides; i++) {

      for(k=0;k<=1;k++) {
	if(bound[j].types[i] == bounds[k]) {
	  parent = bound[j].parent[i];
	  bound[j].parent2[i] = bound[j].parent[i] - minp[k] + minp[(k+1)%2];	  

	  if(!bound[j].ediscont) {
	    bound[j].discont = Ivector(1,bound[j].nosides);
	    for(l=1; l <= bound[j].nosides; l++)
	      bound[j].discont[l] = 0;
	    bound[j].ediscont = TRUE;
	  }

	  bound[j].discont[i] = 2+k;
	  elemtype = data->elementtypes[parent];
	  if(elemtype%100 == 4) {
	    bound[j].side2[i] = (bound[j].side[i] + 2) % 4;
	  }
	  else if(elemtype%100 == 8) {
	    if(bound[j].side[i] < 4) bound[j].side2[i] = (bound[j].side[i] + 2) % 4;
	    if(bound[j].side[i] >= 4) bound[j].side2[i] = 4 + (5 - bound[j].side[i]);
	  }
	}
      } 
    }
  }

  if(info) printf("Periodic boundaries were set with a simple scheme\n");

  return(2);
}



int SetConnectedNodes(struct FemType *data,struct BoundaryType *bound,
		      int bctype,int connecttype,int info)
/* Mark node that are related to a boundary condition of a given bctype.
   This may be used to create strong connections in the partitioning process. */
{
  int i,j,k,bc,sideelemtype,sidenodes,nodesset;
  int sideind[MAXNODESD1],conflicts;

  conflicts = 0;
  nodesset = 0;

  for(bc=0;bc<MAXBOUNDARIES;bc++) {    
    if(bound[bc].created == FALSE) continue;
    if(bound[bc].nosides == 0) continue;
    
    for(i=1;i<=bound[bc].nosides;i++) {
      if( bctype > 0 ) {
	if(bound[bc].types[i] != bctype) continue;
      } 
      else if( bctype == -1 ) {
	if( !bound[bc].parent[i] ) continue;
      }
      else if( bctype == -2 ) {
	if( !bound[bc].parent[i] ) continue;
	if( !bound[bc].parent2[i] ) continue;
      }
      else if( bctype == -3 ) {
	if( !bound[bc].parent[i] ) continue;
	if( bound[bc].parent2[i] ) continue;
      }

      /* If the table pointing the connected nodes does not exist, create it */
      if(!data->nodeconnectexist) {
	data->nodeconnect = Ivector(1,data->noknots);
	for(k=1;k<=data->noknots;k++)
	  data->nodeconnect[k] = 0;
	data->nodeconnectexist = TRUE;
      }
       
      GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
		     data,sideind,&sideelemtype);
      sidenodes = sideelemtype%100;
      
      for(j=0;j<sidenodes;j++) {
	k = sideind[j];
	if( data->nodeconnect[k] != connecttype ) {
	  if( data->nodeconnect[k] ) conflicts += 1;
	  data->nodeconnect[k] = connecttype;
	  nodesset += 1;	  
	}
      }
    }
  }
  if(info) printf("Setting connectivity group %d for %d nodes on boundary %d\n",
		  connecttype,nodesset,bctype);

  if(conflicts) printf("The were %d conflicts in the connectivity set %d\n",
		       conflicts,connecttype);

  return(0);
}


int SetConnectedElements(struct FemType *data,int info)
/* Create connected boundary conditions for a given bctype */
{
  int i,j,k,nonodes,hit,nohits,con;
  int *nodeconnect;

  if(!data->nodeconnectexist) {
    printf("Cannot create connected elements without connected nodes!\n");
    return(1);
  }
  nodeconnect = data->nodeconnect;

  /* Allocated space for the connected elements */
  if(!data->elemconnectexist) {
    printf("Created table for connected elements\n");
    data->elemconnect = Ivector(1,data->noelements);
    for(k=1;k<=data->noelements;k++)
      data->elemconnect[k] = 0;
    data->elemconnectexist = TRUE;
    
    /* Go through all the elements and check which of the elements have 
       nodes that are related to a connected node */
    nohits = 0;  
    for(i=1;i<=data->noelements;i++) {
      nonodes = data->elementtypes[i] % 100;
      hit = FALSE;
      for(j=0;j<nonodes;j++) {
	k = data->topology[i][j];
	con = nodeconnect[k];
	if( con ) {
	  data->elemconnect[i] = MAX( con, data->elemconnect[i] );
	  hit = TRUE;
	}
      }
      if(hit) nohits++;
    }
  
    if(info) printf("Number of connected elements is %d (out of %d)\n",nohits,data->noelements);
    data->elemconnectexist = nohits;
  }

  /* This is a little bit dirty. We set the connections to negative and use the unconnected 
     as a permutation. */
  if( data->elemconnectexist ) {
    if(info) printf("Use connected table as a permutation for creating dual graph!\n");
    j = 0;
    for(i=1;i<=data->noelements;i++) {
      if( data->elemconnect[i] ) {
	data->elemconnect[i] = -abs(data->elemconnect[i]);
      }
      else {
	j++;
	data->elemconnect[i] = j;
      }
    }
  }

  return(0);
}



int FindCorners(struct GridType *grid,struct CellType *cell,
		struct FemType *data,int info)
/* Find the nodes in the mesh that are at material corners. 
   These nodes are often of special interest.
   */
{
  int i,j,k,ind,cellno,elem;
  int allocated,nocorners;

  nocorners = 0;
  allocated = FALSE;
  
omstart:
  
  if(nocorners > 0) {
    data->corners = Ivector(1,2*nocorners);
    data->nocorners = nocorners;
    allocated  = TRUE;
  }
  
  k = 0;
  
    for(i=1;i<=grid->xcells+1;i++) 
      for(j=1;j<=grid->ycells+1;j++) {
      if(grid->structure[j][i] == grid->structure[j][i-1] &&
	 grid->structure[j-1][i] == grid->structure[j-1][i-1])
	continue;
      if(grid->structure[j][i] == grid->structure[j-1][i] &&
	 grid->structure[j][i-1] == grid->structure[j-1][i-1])
	continue;

      /* point (i,j) must now be a corner */
      if(cellno = grid->numbered[j][i]) {
	elem = GetElementIndex(&(cell)[cellno],1,1);
	ind  = BOTLEFT;
      } 
      else if(cellno = grid->numbered[j][i-1]) {
	elem = GetElementIndex(&(cell)[cellno],cell[cellno].xelem,1);
	ind  = BOTRIGHT;
      } 
      else if(cellno = grid->numbered[j-1][i]) {
	elem = GetElementIndex(&(cell)[cellno],1,cell[cellno].yelem);
	ind  = TOPLEFT;
      } 
      else if(cellno = grid->numbered[j-1][i-1]) {
	elem = GetElementIndex(&(cell)[cellno],cell[cellno].xelem,cell[cellno].yelem);
	ind  = TOPRIGHT;
      }
      else continue;

      /* ind is now the index of the corner knot */
      k++;
      
      if(allocated == FALSE) continue;      
      data->corners[2*k-1] = elem;
      data->corners[2*k]   = ind;

    }

  nocorners = k;
  
  if(nocorners == 0) return(0);  
  if(allocated == FALSE) goto omstart;

  if(info) printf("Found %d material corners.\n",nocorners);
  return(0);
}



int ElementsToTriangles(struct FemType *data,struct BoundaryType *bound,
			Real critangle,int info)
/* Make triangles out of rectangular elements */
{
  int i,j,k,l,side,elem,i1,isum,sideelemtype;
  int noelements,elementtype,triangles,noknots,nonodes,newelements,newtype,newmaxnodes;
  int **newtopo=NULL,*newmaterial=NULL,*newelementtypes=NULL;
  int newnodes,*needed=NULL,*divisions=NULL,*division1=NULL;
  int sideind[MAXNODESD1], sideind2[MAXNODESD1];
  int allocated,maxanglej,evenodd,newelem;
  Real dx1,dx2,dy1,dy2,ds1,ds2;
  Real angles[4],maxangle;
  struct FemType data2;

  noelements  = data->noelements;
  noknots = data->noknots;
  allocated = FALSE;

  needed = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    needed[i] = 0;
  for(i=1;i<=noelements;i++) {
    nonodes = data->elementtypes[i] / 100;
    for(j=0;j<nonodes;j++)
      needed[data->topology[i][j]] += 1;
  }

  divisions = Ivector(1,noelements);
  division1 = Ivector(1,noelements);
  for(i=1;i<=noelements;i++)
    divisions[i] = division1[i] = 0;

  /* First divide the elements along the shorter diameter */

  newelements = 0;
  newmaxnodes = 0;

 omstart:

  for(i=1;i<=noelements;i++) {

    elementtype = data->elementtypes[i];

    /* compute the four angles and divide the rectangle so that the largest angle is split */
    maxangle = 0.0;
    maxanglej = 0;
    for(j=0;j<4;j++) {
      dx1 = data->x[data->topology[i][(j+3)%4]] - data->x[data->topology[i][j]];
      dy1 = data->y[data->topology[i][(j+3)%4]] - data->y[data->topology[i][j]];
      dx2 = data->x[data->topology[i][(j+1)%4]] - data->x[data->topology[i][j]];
      dy2 = data->y[data->topology[i][(j+1)%4]] - data->y[data->topology[i][j]];
      ds1 = sqrt(dx1*dx1+dy1*dy1);
      ds2 = sqrt(dx2*dx2+dy2*dy2);
      angles[j] = (180.0/FM_PI) * acos((dx1*dx2+dy1*dy2)/(ds1*ds2));
      
      /* Slightly favor divisions where corner is split */
      if(needed[data->topology[i][j]] == 1) angles[j] *= 1.001;

      if( abs(angles[j] > maxangle)) {
	maxangle = abs(angles[j]);
	maxanglej = j;
      }
    }    
    evenodd = maxanglej % 2;


    /* No triangularization is performed unless the critical angle is exceeded. */
    if( maxangle < critangle ) {
      triangles = 1;
      newtype = elementtype;
      newnodes = elementtype % 100;
    }
    else {
      switch(elementtype) {
      case 404:
	newtype = 303;
	newnodes = 3;
	triangles = 2;
	break;
      case 405:
	newtype = 303;
	newnodes = 3;
	triangles = 4;
	break;
      case 409:
	newtype = 306;
	newnodes = 6;
	triangles = 2;
	break;
      case 416:
	newtype = 310;
	newnodes = 10;
	triangles = 2;
	break;
      default:
	printf("ElementsToTriangles: not implemented for elementtype %d\n",elementtype);
	return(1);
      }
    }

    newmaxnodes = MAX( newnodes, newmaxnodes );
    

    if(!allocated) {
      divisions[i] = triangles;
      division1[i] = newelements;
      newelements += triangles; 
      continue;
    }

    for(j=division1[i]+1;j<=division1[i]+divisions[i];j++) {
      newelementtypes[j] = newtype;
      newmaterial[j] = data->material[i];
    }

    newelem = division1[i]+1;
    if(triangles == 1) {
      for(j=0;j<newnodes;j++)
	newtopo[newelem][j] = data->topology[i][j];
    }
    else {
      if(elementtype == 404 || elementtype == 409 || elementtype == 416) {
	if(evenodd) {
	  newtopo[newelem][0] = data->topology[i][0];
	  newtopo[newelem][1] = data->topology[i][1];
	  newtopo[newelem][2] = data->topology[i][3];
	  newtopo[newelem+1][0] = data->topology[i][2];
	  newtopo[newelem+1][1] = data->topology[i][3];
	  newtopo[newelem+1][2] = data->topology[i][1];
	}
	else {
	  newtopo[newelem][0] = data->topology[i][1];
	  newtopo[newelem][1] = data->topology[i][2];
	  newtopo[newelem][2] = data->topology[i][0];
	  newtopo[newelem+1][0] = data->topology[i][3];
	  newtopo[newelem+1][1] = data->topology[i][0];
	  newtopo[newelem+1][2] = data->topology[i][2];
	}
      }
      if(elementtype == 409) {
	if(evenodd) {
	  newtopo[newelem][3] = data->topology[i][4];
	  newtopo[newelem][4] = data->topology[i][8];
	  newtopo[newelem][5] = data->topology[i][7];
	  newtopo[newelem+1][3] = data->topology[i][6];
	  newtopo[newelem+1][4] = data->topology[i][8];
	  newtopo[newelem+1][5] = data->topology[i][5];      
	}
	else {
	  newtopo[newelem][3] = data->topology[i][5];
	  newtopo[newelem][4] = data->topology[i][8];
	  newtopo[newelem][5] = data->topology[i][4];
	  newtopo[newelem+1][3] = data->topology[i][7];
	  newtopo[newelem+1][4] = data->topology[i][8];
	  newtopo[newelem+1][5] = data->topology[i][6];            	
	}
      }
      if(elementtype == 416) {
	if(evenodd) {
	  newtopo[newelem][3] = data->topology[i][4];
	  newtopo[newelem][4] = data->topology[i][5];
	  newtopo[newelem][5] = data->topology[i][13];
	  newtopo[newelem][6] = data->topology[i][15];
	  newtopo[newelem][7] = data->topology[i][10];
	  newtopo[newelem][8] = data->topology[i][11];
	  newtopo[newelem][9] = data->topology[i][12];
	  
	  newtopo[newelem+1][3] = data->topology[i][8];
	  newtopo[newelem+1][4] = data->topology[i][9];
	  newtopo[newelem+1][5] = data->topology[i][15];
	  newtopo[newelem+1][6] = data->topology[i][13];
	  newtopo[newelem+1][7] = data->topology[i][6];      
	  newtopo[newelem+1][8] = data->topology[i][7];      
	  newtopo[newelem+1][9] = data->topology[i][14];      
	}
	else {
	  newtopo[newelem][3] = data->topology[i][6];
	  newtopo[newelem][4] = data->topology[i][7];
	  newtopo[newelem][5] = data->topology[i][14];
	  newtopo[newelem][6] = data->topology[i][12];
	  newtopo[newelem][7] = data->topology[i][4];
	  newtopo[newelem][8] = data->topology[i][5];
	  newtopo[newelem][9] = data->topology[i][13];
	  
	  newtopo[newelem+1][3] = data->topology[i][10];
	  newtopo[newelem+1][4] = data->topology[i][11];
	  newtopo[newelem+1][5] = data->topology[i][12];
	  newtopo[newelem+1][6] = data->topology[i][14];
	  newtopo[newelem+1][7] = data->topology[i][8];            	
	  newtopo[newelem+1][8] = data->topology[i][9];            	
	  newtopo[newelem+1][9] = data->topology[i][15];            	
	}
      }
      else if(elementtype == 405) {
	newtopo[newelem][0] = data->topology[i][0];
	newtopo[newelem][1] = data->topology[i][1];
	newtopo[newelem][2] = data->topology[i][4];
	newtopo[newelem+1][0] = data->topology[i][1];
	newtopo[newelem+1][1] = data->topology[i][2];
	newtopo[newelem+1][2] = data->topology[i][4];
	newtopo[newelem+2][0] = data->topology[i][2];
	newtopo[newelem+2][1] = data->topology[i][3];
	newtopo[newelem+2][2] = data->topology[i][4];
	newtopo[newelem+3][0] = data->topology[i][3];
	newtopo[newelem+3][1] = data->topology[i][0];
	newtopo[newelem+3][2] = data->topology[i][4];
      }
    }
  }


  if(!allocated)  {

    newtopo = Imatrix(1,newelements,0,newmaxnodes-1);
    newmaterial = Ivector(1,newelements);
    newelementtypes = Ivector(1,newelements);
    allocated = TRUE;
    
    data2 = *data;
    data2.topology = newtopo;
    data2.material = newmaterial;
    data2.elementtypes = newelementtypes;

    goto omstart;
  }


  /* Then make the corresponding mapping for the BCs. 
     This is done in a brute-force way where all the 
     possible new elements are checked. */

  for(j=0;j < MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    
    for(i=1; i <= bound[j].nosides; i++) {

      for(l=1;l<=2;l++) {
	
	if(l==1) 
	  k = bound[j].parent[i];
	else 
	  k = bound[j].parent2[i];
	if(k == 0) continue;


	if(l == 1) 
	  GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
			 data,sideind,&sideelemtype);
	else 
	  GetElementSide(bound[j].parent2[i],bound[j].side2[i],bound[j].normal[i],
			 data,sideind,&sideelemtype);	
	
	if(sideelemtype/100 != 2) {
	  printf("ElementsToTriangles: implemented only for BCs 202 and 203\n");
	  continue;
	}
	
	isum = 0;
	side = 0;
	if(divisions[k] == 1) {
	  elem = division1[k]+1;
	  side = bound[j].side[i];
	  isum = 2;
	  goto nextparent;
	}

	/* Test for all possible elements that could be parents */       
	for(elem=division1[k]+1;elem<=division1[k]+divisions[k];elem++) {
	  isum = 0;
	  for(i1=0;i1<3;i1++) {
	    if(newtopo[elem][i1] == sideind[0]) isum++;
	    if(newtopo[elem][i1] == sideind[1]) isum++;
	  }
	  
	  if(isum != 2) continue;
	  
	  for(side=0;side<3;side++) { 
	    GetElementSide(elem,side,bound[j].normal[i],
			   &data2,sideind2,&sideelemtype);
	    isum = 0;
	    for(i1=0;i1<2;i1++) {
	      if(sideind2[i1] == sideind[0]) isum++;
	      if(sideind2[i1] == sideind[1]) isum++;
	    }
	    
	    if(isum == 2) goto nextparent;
	  }
	}
	
      nextparent:
	if(isum == 2) {
	    if(l == 1) {
	      bound[j].parent[i] = elem;
	      bound[j].side[i] = side;
	    }
	    if(l == 2) {
	      bound[j].parent2[i] = elem;
	      bound[j].side2[i] = side;
	    }
	}     
	else {
	  printf("Failed to find parent for side %d of %d (parent %d)\n",i,j,k);
	} 
      }
    }
  }


  free_Imatrix(data->topology,1,noelements,0,data->maxnodes-1);
  free_Ivector(data->material,1,noelements);
  free_Ivector(data->elementtypes,1,noelements);
  free_Ivector(needed,1,noknots);
  free_Ivector(divisions,1,noelements);
  free_Ivector(division1,1,noelements);

  data->topology = newtopo;
  data->elementtypes = newelementtypes;
  data->material = newmaterial;
  data->noelements = newelements;
  data->maxnodes = newmaxnodes;

  if(info) printf("There are %d elements after triangularization (was %d)\n",
		  newelements,noelements);

  return(0);
}


int PolarCoordinates(struct FemType *data,Real rad,int info)
{
  int i;
  Real fii,zet,dr;

  for(i=1;i<=data->noknots;i++) {
    zet = data->x[i];
    fii = FM_PI/180.0 * data->y[i];
    dr = data->z[i];

    data->z[i] = zet;
    data->x[i] = (rad+dr) * cos(fii);
    data->y[i] = (rad+dr) * sin(fii);
  }

  if(info) printf("Making coordinate transformation from polar to cartesian\n");

  return(0);
}


int CylinderCoordinates(struct FemType *data,int info)
{
  int i;
  Real x,y,z,rad,fii;
  
  if( data->dim == 3 ) {    
    for(i=1;i<=data->noknots;i++) {
      x = data->x[i];
      y = data->y[i];
      fii = FM_PI/180.0 * data->z[i];
      rad = data->x[i];
      
      data->x[i] = rad * cos(fii);
      data->y[i] = rad * sin(fii);
      data->z[i] = y;
    }
  }
  else {
    for(i=1;i<=data->noknots;i++) {
      rad = data->x[i];
      fii = FM_PI/180.0 * data->y[i];
      
      data->x[i] = rad * cos(fii);
      data->y[i] = rad * sin(fii);
    }
  }
    
  if(info) printf("Making coordinate transformation from cylindrical to cartesian\n");

  return(0);
}


int UniteMeshes(struct FemType *data1,struct FemType *data2,
		struct BoundaryType *bound1,struct BoundaryType *bound2,
		int nooverlap, int info)
/* Unites two meshes for one larger mesh */
{
  int i,j,k;
  int noelements,noknots,nonodes,maxnodes;
  int **newtopo=NULL,*newmaterial=NULL,*newelementtypes=NULL;
  Real *newx=NULL,*newy=NULL,*newz=NULL;
  int mat,usenames,*bodynameis,*boundarynameis,*bodyused,*boundaryused;
  int bcmax1,bcmin2,bcoffset;
  int bodymax1,bodymin2,bodyoffset;

  noknots = data1->noknots + data2->noknots;
  noelements  = data1->noelements + data2->noelements;
  maxnodes = MAX(data1->maxnodes,data2->maxnodes);

  if(data2->dim > data1->dim) data1->dim = data2->dim;

  if(0) printf("Uniting two meshes to %d nodes and %d elements.\n",noknots,noelements);

  usenames = data1->bodynamesexist || data1->boundarynamesexist; 
  bcoffset = 0; bodyoffset = 0;

  if( usenames ) {
    bodynameis = Ivector(1,MAXBODIES);
    boundarynameis = Ivector(1,MAXBCS);
    bodyused = Ivector(1,MAXBODIES);
    boundaryused = Ivector(1,MAXBCS);

    for(i=1;i<=MAXBODIES;i++)
      bodynameis[i] = bodyused[i] = FALSE;
    for(i=1;i<=MAXBCS;i++)
      boundarynameis[i] = boundaryused[i] = FALSE;

    /* First mark the original bodies and boundaries that maintain their index */
    for(i=1;i<=data1->noelements;i++) {
      mat = data1->material[i];
      if( mat < MAXBODIES ) {
	if(!bodynameis[mat]) {
	  bodynameis[mat] = -1;
	  bodyused[mat] = TRUE;
	}
      }
    }

    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound1[j].created) continue;     
      for(i=1; i <= bound1[k].nosides; i++) {
	mat = bound1[j].types[i];
	if( mat < MAXBCS ) {
	  if(!boundarynameis[mat]) {
	    boundarynameis[mat] = -1;
	    boundaryused[mat] = TRUE;
	  }
	}
      }
    }

    /* Then mark the joined bodies and boundaries that are not conflicting */ 
    for(i=1;i<=data2->noelements;i++) {
      mat = data2->material[i];
      if( mat < MAXBODIES ) {
	if( !bodynameis[mat] ) {
	  bodynameis[mat] = mat;
	  bodyused[mat] = TRUE;
	  if(!data1->bodyname[mat]) data1->bodyname[mat] = Cvector(0,MAXNAMESIZE);
	  strcpy(data1->bodyname[mat],data2->bodyname[mat]);
	}
      }
    }

    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound2[j].created) continue;     

      for(i=1; i <= bound2[j].nosides; i++) {
	mat = bound2[j].types[i];
	if( mat < MAXBCS ) {
	  if( !boundarynameis[mat] ) {
	    boundarynameis[mat] = mat;
	    boundaryused[mat] = TRUE;
	    if(!data1->boundaryname[mat]) data1->boundaryname[mat] = Cvector(0,MAXNAMESIZE);
	    strcpy(data1->boundaryname[mat],data2->boundaryname[mat]);
	  }
	}
      }
    }


    /* And finally number the conflicting joined bodies and BCs */
    for(i=1;i<=data2->noelements;i++) {
      mat = data2->material[i];
      if( mat < MAXBODIES ) {
	if( bodynameis[mat] == -1) {
	  for(k=1;k<MAXBODIES;k++)
	    if(!bodyused[k]) break;
	  if(info) printf("Renumbering body %d to %d\n",mat,k);
	  bodynameis[mat] = k;	  
	  bodyused[k] = TRUE;
	  if(!data1->bodyname[k]) data1->bodyname[k] = Cvector(0,MAXNAMESIZE);
  	  strcpy(data1->bodyname[k],data2->bodyname[mat]);
	}
      }
    }

    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound2[j].created) continue;     
      for(i=1; i <= bound2[j].nosides; i++) {
	mat = bound2[j].types[i];

	if( mat < MAXBCS ) {
	  if( boundarynameis[mat] == -1) {
	    for(k=1;k<MAXBCS;k++) 
	      if(!boundaryused[k]) break;
	    if(info) printf("Renumbering boundary %d to %d\n",mat,k);
	    boundarynameis[mat] = k;
	    boundaryused[k] = TRUE;
	    if(!data1->boundaryname[k]) data1->boundaryname[k] = Cvector(0,MAXNAMESIZE);
	    strcpy(data1->boundaryname[k],data2->boundaryname[mat]);
	  }
	}
      }
    }

  }
  else if (nooverlap ) {
    bcmax1 = 0;
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound1[j].created) continue;     
      for(i=1; i <= bound1[k].nosides; i++) {
	mat = bound1[j].types[i];
	bcmax1 = MAX( bcmax1, mat );	
      }
    }

    bcmin2 = 1000;
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound2[j].created) continue;     

      for(i=1; i <= bound2[j].nosides; i++) {
	mat = bound2[j].types[i];
	bcmin2 = MIN( bcmin2, mat );
      }
    }
    bcoffset = MAX(0, bcmax1 - bcmin2 + 1);
    if( info ) {
      printf("Max(bc1) is %d and Min(bc2) is %d, using BC offset %d for mesh 2!\n",bcmax1,bcmin2,bcoffset);
    }
        
    bodymax1 = 0;
    for(i=1;i<=data1->noelements;i++) {
      mat = data1->material[i];
      bodymax1 = MAX( bodymax1, mat );	
    }

    bodymin2 = 1000;
    for(i=1;i<=data2->noelements;i++) {
      mat = data2->material[i];
      bodymin2 = MIN( bodymin2, mat );	
    }
    bodyoffset = MAX(0, bodymax1 - bodymin2 + 1);
    if( info ) {
      printf("Max(body1) is %d and Min(body2) is %d, using body offset %d for mesh 2!\n",bodymax1,bodymin2,bodyoffset);
    }
  }
  


  for(j=0;j < MAXBOUNDARIES;j++) {
    if(!bound2[j].created) continue;

    for(k=j;k < MAXBOUNDARIES;k++)
      if(!bound1[k].created) break;
 
    bound1[k].created = bound2[j].created;
    bound1[k].nosides = bound2[j].nosides;
    bound1[k].coordsystem = bound2[j].coordsystem;
    bound1[k].side = bound2[j].side;
    bound1[k].side2 = bound2[j].side2;
    bound1[k].parent = bound2[j].parent;
    bound1[k].parent2 = bound2[j].parent2;
    bound1[k].material = bound2[j].material;
    bound1[k].echain = bound2[j].echain;
    bound1[k].types = bound2[j].types;
    bound1[k].normal = bound2[j].normal;

    for(i=1; i <= bound1[k].nosides; i++) {
      bound1[k].parent[i] += data1->noelements;
      if(bound1[k].parent2[i])
	bound1[k].parent2[i] += data1->noelements;	 
      
      mat = bound2[j].types[i];
      if( usenames ) {
	if( mat < MAXBCS ) {
	  bound1[k].types[i] = boundarynameis[mat];
	}
      } else {
	bound1[k].types[i] = bcoffset + mat;
      }
    }
  }

  data1->maxnodes = maxnodes;
  newtopo = Imatrix(1,noelements,0,maxnodes-1);
  newmaterial = Ivector(1,noelements);
  newelementtypes = Ivector(1,noelements);
  newx = Rvector(1,noknots);
  newy = Rvector(1,noknots);
  newz = Rvector(1,noknots);

  for(i=1;i<=data1->noknots;i++) {
    newx[i] = data1->x[i];
    newy[i] = data1->y[i];
    newz[i] = data1->z[i];
  }
  for(i=1;i<=data2->noknots;i++) {
    newx[i+data1->noknots] = data2->x[i];
    newy[i+data1->noknots] = data2->y[i];
    newz[i+data1->noknots] = data2->z[i];
  }

  for(i=1;i<=data1->noelements;i++) {
    mat = data1->material[i];
    newmaterial[i] = mat;
    newelementtypes[i] = data1->elementtypes[i]; 
    nonodes = newelementtypes[i]%100;
    for(j=0;j<nonodes;j++)
      newtopo[i][j] = data1->topology[i][j];
  }
  for(i=1;i<=data2->noelements;i++) {
    mat = data2->material[i];
    newelementtypes[i+data1->noelements] = data2->elementtypes[i]; 
    nonodes = newelementtypes[i+data1->noelements]%100;
    for(j=0;j<nonodes;j++)
      newtopo[i+data1->noelements][j] = data2->topology[i][j] + data1->noknots;

    if( usenames ) {
      if( mat < MAXBODIES ) {
        newmaterial[i+data1->noelements] = bodynameis[mat];
      }
    }
    else {
      newmaterial[i+data1->noelements] = bodyoffset + mat;
    }
  }

  free_Imatrix(data1->topology,1,data1->noelements,0,data1->maxnodes-1);
  free_Ivector(data1->material,1,data1->noelements);
  free_Rvector(data1->x,1,data1->noknots);
  free_Rvector(data1->y,1,data1->noknots);
  free_Rvector(data1->z,1,data1->noknots);

  free_Imatrix(data2->topology,1,data2->noelements,0,data2->maxnodes-1);
  free_Ivector(data2->material,1,data2->noelements);
  free_Rvector(data2->x,1,data2->noknots);
  free_Rvector(data2->y,1,data2->noknots);
  free_Rvector(data2->z,1,data2->noknots);
  
  data1->noelements = noelements;
  data1->noknots  = noknots;
  data1->topology = newtopo;
  data1->material = newmaterial;
  data1->elementtypes = newelementtypes; 
  data1->x = newx;
  data1->y = newy;
  data1->z = newz;

  if(info) printf("Two meshes were united to one with %d nodes and %d elements.\n",
		  noknots,noelements);

  return(0);
}


int CloneMeshes(struct FemType *data,struct BoundaryType *bound,
		int *ncopies,Real *meshsize,int diffmats,int info)
/* Unites two meshes for one larger mesh */
{
  int i,j,k,l,m;
  int noelements,noknots,nonodes,totcopies,ind,origdim;
  int **newtopo=NULL,*newmaterial=NULL,*newelementtypes=NULL,maxnodes;
  int maxmaterial,maxtype,ncopy,bndr,nosides;
  Real *newx=NULL,*newy=NULL,*newz=NULL;
  Real maxcoord[3],mincoord[3];

  int *vparent=NULL,*vparent2=NULL,*vside=NULL,*vside2=NULL;
  int *vtypes=NULL,*vmaterial=NULL,*vnormal=NULL,*vdiscont=NULL;
  
  if(info) printf("CloneMeshes: copying the mesh to a matrix\n");
  if(diffmats) {
    if(info) printf("CloneMeshes: giving each new entity new material and bc indexes\n");
  }
  
  origdim = data->dim;
  totcopies = 1;
  if( ncopies[2] > 1 ) {
    data->dim = 3;
  }
  else {
    ncopies[2] = 1;
  }

  for(i=0;i<data->dim;i++) {
    if(ncopies[i] > 1) totcopies *= ncopies[i];
  }

  maxcoord[0] = mincoord[0] = data->x[1];
  maxcoord[1] = mincoord[1] = data->y[1];
  maxcoord[2] = mincoord[2] = data->z[1];

  for(i=1;i<=data->noknots;i++) {
    if(data->x[i] > maxcoord[0]) maxcoord[0] = data->x[i]; 
    if(data->x[i] < mincoord[0]) mincoord[0] = data->x[i]; 
    if(data->y[i] > maxcoord[1]) maxcoord[1] = data->y[i]; 
    if(data->y[i] < mincoord[1]) mincoord[1] = data->y[i]; 
    if(data->z[i] > maxcoord[2]) maxcoord[2] = data->z[i]; 
    if(data->z[i] < mincoord[2]) mincoord[2] = data->z[i]; 
  }

  for(i=0;i<origdim;i++) {
    if(maxcoord[i]-mincoord[i] > meshsize[i]) meshsize[i] = maxcoord[i]-mincoord[i];
  }
  if(info) printf("Meshsize to be copied: %lg %lg %lg\n",meshsize[0],meshsize[1],meshsize[2]);

  noknots = totcopies * data->noknots;
  noelements  = totcopies * data->noelements;
  maxnodes = data->maxnodes;

  if(info) printf("Copying the mesh to %d identical domains in %d-dim.\n",totcopies,data->dim);

  data->maxnodes = maxnodes;
  newtopo = Imatrix(1,noelements,0,maxnodes-1);
  newmaterial = Ivector(1,noelements);
  newelementtypes = Ivector(1,noelements);
  newx = Rvector(1,noknots);
  newy = Rvector(1,noknots);
  newz = Rvector(1,noknots);

  for(l=0;l<ncopies[2];l++) {
    for(k=0;k<ncopies[1];k++) {
      for(j=0;j<ncopies[0];j++) {
	for(i=1;i<=data->noknots;i++) {
	  ncopy = j+k*ncopies[0]+l*ncopies[0]*ncopies[1];
	  ind = i + ncopy*data->noknots;

	  newx[ind] = data->x[i] + j*meshsize[0];
	  newy[ind] = data->y[i] + k*meshsize[1];
	  newz[ind] = data->z[i] + l*meshsize[2];
	}
      }
    }
  }

  maxmaterial = 0;
  if( diffmats ) {
    for(i=1;i<=data->noelements;i++) 
      if(data->material[i] > maxmaterial) maxmaterial = data->material[i];
    if(info ) printf("Material offset for cloning set to: %d\n",maxmaterial);
  }

  for(l=0;l<ncopies[2];l++) {
    for(k=0;k<ncopies[1];k++) {
      for(j=0;j<ncopies[0];j++) {
	for(i=1;i<=data->noelements;i++) {
	  ncopy = j+k*ncopies[0]+l*ncopies[1]*ncopies[0];
	  ind =  i + ncopy*data->noelements;

	  newmaterial[ind] = data->material[i] + diffmats*maxmaterial*ncopy;
	  newelementtypes[ind] = data->elementtypes[i]; 
	  nonodes = newelementtypes[i]%100;
	  for(m=0;m<nonodes;m++)
	    newtopo[ind][m] = data->topology[i][m] + ncopy*data->noknots;
	}
      }
    }
  }

  maxtype = 0;
  if( diffmats ) {
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      for(i=1; i <= bound[j].nosides; i++) 
	if(maxtype < bound[j].types[i]) maxtype = bound[j].types[i];
    }
    if(info ) printf("Boundary offset for cloning set to: %d\n",maxtype);
  }

  for(bndr=0;bndr < MAXBOUNDARIES;bndr++) {

    if(!bound[bndr].created) continue;

    nosides = totcopies * bound[bndr].nosides;

    vparent = Ivector(1, nosides);
    vparent2 = Ivector(1, nosides);
    vside = Ivector(1, nosides);
    vside2 = Ivector(1, nosides);
    vmaterial = Ivector(1, nosides);
    vtypes = Ivector(1, nosides);
    vnormal = Ivector(1, nosides);

    if(bound[bndr].ediscont) { 
      vdiscont = Ivector(1, nosides);
      for(i=1; i <= nosides; i++) 
	vdiscont[i] = 0;
    }

    for(l=0;l<ncopies[2];l++) {
      for(k=0;k<ncopies[1];k++) {
	for(j=0;j<ncopies[0];j++) {
	  for(i=1; i <= bound[bndr].nosides; i++) {

	    ncopy = j+k*ncopies[0]+l*ncopies[1]*ncopies[0];
	    ind = i + ncopy * bound[bndr].nosides;
	    
	    vparent[ind] = bound[bndr].parent[i] + ncopy * data->noelements;
	    vside[ind] = bound[bndr].side[i]; 

	    if(bound[bndr].parent2[i]) {
	      vparent2[ind] = bound[bndr].parent2[i] + ncopy * data->noelements;
	      vside2[ind] = bound[bndr].side2[i]; 
	    }
	    else {
	      vparent2[ind] = 0.0;
	      vside2[ind] = 0.0;
	    }

	    vnormal[ind] = bound[bndr].normal[i]; 

	    if(bound[bndr].ediscont) 
	      vdiscont[ind] = bound[bndr].discont[i]; 

	    vtypes[ind] = bound[bndr].types[i] + diffmats * ncopy * maxtype;

	    vmaterial[ind] = bound[bndr].material[i] + ncopy * maxmaterial;
	  }
	}
      }
    }
   
    bound[bndr].nosides = nosides;
    bound[bndr].side = vside;

    bound[bndr].side2 = vside2;
    bound[bndr].parent = vparent;
    bound[bndr].parent2 = vparent2;
    bound[bndr].types = vtypes;
    bound[bndr].material = vmaterial;
    if(bound[bndr].ediscont) 
      bound[bndr].discont = vdiscont;
  }

  free_Imatrix(data->topology,1,data->noelements,0,data->maxnodes-1);
  free_Ivector(data->material,1,data->noelements);
  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);

  data->noelements = noelements;
  data->noknots  = noknots;
  data->topology = newtopo;
  data->material = newmaterial;

  data->elementtypes = newelementtypes; 
  data->x = newx;
  data->y = newy;
  data->z = newz;

  if( data->bodynamesexist || data->boundarynamesexist ) {
    printf("Cloning cannot treat names yet, omitting treatment of names for now!\n");
    data->bodynamesexist = FALSE;
    data->boundarynamesexist = FALSE;
  } 

  if(info) printf("The mesh was copied to several identical meshes\n");

  return(0);
}


int MirrorMeshes(struct FemType *data,struct BoundaryType *bound,
		 int *symmaxis,int diffmats,Real *meshsize,int symmbound,int info)
/* Makes a mirror image of a mesh and unites it with the original mesh */
{
  int i,j,m;
  int noelements,noknots,nonodes,totcopies,ind,maxnodes;
  int **newtopo=NULL,*newmaterial=NULL,*newelementtypes=NULL;
  int maxtype,bndr,nosides;
  Real *newx=NULL,*newy=NULL,*newz=NULL;
  Real maxcoord[3],mincoord[3];
  int ind0,elem0,axis1,axis2,axis3,symmcount;

  int *vparent=NULL,*vparent2=NULL,*vside=NULL,*vside2=NULL;
  int *vtypes=NULL,*vmaterial=NULL,*vnormal=NULL,*vdiscont=NULL;
  
  printf("MirrorMeshes: making a symmetric mapping of the mesh\n");

  if(symmaxis[0]) symmaxis[0] = 1;
  if(symmaxis[1]) symmaxis[1] = 1;
  if(symmaxis[2]) symmaxis[2] = 1;
  if(data->dim < 3) symmaxis[2] = 0;

  maxcoord[0] = mincoord[0] = data->x[1];
  maxcoord[1] = mincoord[1] = data->y[1];
  maxcoord[2] = mincoord[2] = data->z[1];

  for(i=1;i<=data->noknots;i++) {
    if(data->x[i] > maxcoord[0]) maxcoord[0] = data->x[i]; 
    if(data->x[i] < mincoord[0]) mincoord[0] = data->x[i]; 
    if(data->y[i] > maxcoord[1]) maxcoord[1] = data->y[i]; 
    if(data->y[i] < mincoord[1]) mincoord[1] = data->y[i]; 
    if(data->z[i] > maxcoord[2]) maxcoord[2] = data->z[i]; 
    if(data->z[i] < mincoord[2]) mincoord[2] = data->z[i]; 
  }

  for(i=0;i<3;i++) {
    if(maxcoord[i]-mincoord[i] > meshsize[i]) meshsize[i] = maxcoord[i]-mincoord[i];
  }

  if(diffmats) diffmats = 1;
  
  totcopies = 1;
  for(i=0;i<3;i++)
    if(symmaxis[i]) totcopies *= 2;

  noknots = totcopies * data->noknots;
  noelements  = totcopies * data->noelements;
  maxnodes = data->maxnodes;

  printf("Mirroring the mesh to %d symmetrical domains.\n",totcopies);

  data->maxnodes = maxnodes;
  newtopo = Imatrix(1,noelements,0,maxnodes-1);
  newmaterial = Ivector(1,noelements);
  newelementtypes = Ivector(1,noelements);
  newx = Rvector(1,noknots);
  newy = Rvector(1,noknots);
  newz = Rvector(1,noknots);

  ind0 = 0;
  

  for(axis1=0;axis1 <= symmaxis[0];axis1++) {
    for(axis2=0;axis2 <= symmaxis[1];axis2++) {
      for(axis3=0;axis3 <= symmaxis[2];axis3++) {

	for(i=1;i<=data->noknots;i++) {
	  ind = i + ind0;
	  
	  newx[ind] = (1-2*axis1) * data->x[i];
	  newy[ind] = (1-2*axis2) * data->y[i];
	  newz[ind] = (1-2*axis3) * data->z[i];
	  
	  newmaterial[ind] = data->material[i];
	  newelementtypes[ind] = data->elementtypes[i]; 
	}
	ind0 += data->noknots;
      }
    }
  }

  elem0 = 0;
  ind0 = 0;

  for(axis1=0;axis1 <= symmaxis[0];axis1++) {
    for(axis2=0;axis2 <= symmaxis[1];axis2++) {
      for(axis3=0;axis3 <= symmaxis[2];axis3++) {
	
	for(i=1;i<=data->noelements;i++) {
	  ind =  i + elem0;
	  newmaterial[ind] = data->material[i];
	  newelementtypes[ind] = data->elementtypes[i]; 
	  nonodes = newelementtypes[i]%100;
	  for(m=0;m<nonodes;m++)
	    newtopo[ind][m] = data->topology[i][m] + ind0;
	}
	
	elem0 += data->noelements;
	ind0 += data->noknots;
	printf("elem0=%d ind0=%d\n",elem0,ind0);
      }
    }
  }

  maxtype = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    for(i=1; i <= bound[j].nosides; i++) 
      if(maxtype < bound[j].types[i]) maxtype = bound[j].types[i];
  }

  for(bndr=0;bndr < MAXBOUNDARIES;bndr++) {

    if(!bound[bndr].created) continue;
    nosides = totcopies * bound[bndr].nosides;
    ind = 0;

    vparent = Ivector(1, nosides);
    vparent2 = Ivector(1, nosides);
    vside = Ivector(1, nosides);
    vside2 = Ivector(1, nosides);
    vmaterial = Ivector(1, nosides);
    vtypes = Ivector(1, nosides);
    vnormal = Ivector(1, nosides);

    if(bound[bndr].ediscont) { 
      vdiscont = Ivector(1, nosides);
      for(i=1;i<=nosides;i++)
	vdiscont[i] = 0;
    }

    symmcount = 0;
    elem0 = 0;
    ind0 = 0;

    for(axis1=0;axis1 <= symmaxis[0];axis1++) {
      for(axis2=0;axis2 <= symmaxis[1];axis2++) {
	for(axis3=0;axis3 <= symmaxis[2];axis3++) {
	  
	  for(i=1; i <= bound[bndr].nosides; i++) {
	    
	    if(bound[bndr].types[i] == symmbound) continue;
	    ind++;
	    
	    vparent[ind] = bound[bndr].parent[i] + elem0;
	    vparent2[ind] = bound[bndr].parent2[i] + elem0;
	    vside[ind] = bound[bndr].side[i]; 
	    vside2[ind] = bound[bndr].side2[i]; 
	    
	    vnormal[ind] = bound[bndr].normal[i]; 

	    if(bound[bndr].ediscont) 
	      vdiscont[ind] = bound[bndr].discont[i]; 

	    vtypes[ind] = bound[bndr].types[i] + diffmats * symmcount * maxtype;
	    
	    vmaterial[ind] = bound[bndr].material[i];
	  }

	  symmcount++;
	  elem0 += data->noelements;
	}
      }
    }

    nosides = ind;
    bound[bndr].nosides = nosides;
    bound[bndr].side = vside;
    bound[bndr].side2 = vside2;
    bound[bndr].parent = vparent;
    bound[bndr].parent2 = vparent2;
    bound[bndr].types = vtypes;
    bound[bndr].material = vmaterial;
    if(bound[bndr].ediscont) 
      bound[bndr].discont = vdiscont;
  }

  free_Imatrix(data->topology,1,data->noelements,0,data->maxnodes-1);
  free_Ivector(data->material,1,data->noelements);
  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);

  data->noelements = noelements;
  data->noknots  = noknots;
  data->topology = newtopo;
  data->material = newmaterial;
  data->elementtypes = newelementtypes; 
  data->x = newx;
  data->y = newy;
  data->z = newz;

  if( data->bodynamesexist || data->boundarynamesexist ) {
    printf("Mirroring cannot treat names yet, omitting treatment of names for now!\n");
    data->bodynamesexist = FALSE;
    data->boundarynamesexist = FALSE;
  } 

  if(info) printf("The mesh was copied to several identical meshes\n");

  return(0);
}



static void ReorderAutomatic(struct FemType *data,int iterations,
			     int *origindx,Real corder[],int info)
{
  int i,j,k,l,nonodes,maxnodes,noelements,noknots,minwidth,indexwidth;
  int **neighbours=NULL,*newrank=NULL,*newindx=NULL,*oldrank=NULL,*oldindx=NULL;
  int nocands,*cands=NULL,ind,ind2,cantdo;
  int elemtype,indready,iter,*localorder=NULL,*localtmp=NULL,nolocal;
  Real *localdist=NULL,dx,dy,dz;

  iterations = 3;
  iter = 0;
  maxnodes = 8;

  noelements = data->noelements;
  noknots = data->noknots;

  cantdo = FALSE;
  for(j=1;j<=noelements;j++) {
    elemtype = data->elementtypes[j];
    if(elemtype != 404 && elemtype != 303 && elemtype != 808) cantdo = elemtype;
  }
  if(cantdo) {
    printf("Automatic reordering not specified for elementtype %d\n",cantdo);
    return;
  }

  printf("Allocating...\n");

  cands = Ivector(1,maxnodes);
  localorder = Ivector(1,maxnodes);
  localtmp = Ivector(1,maxnodes);
  localdist = Rvector(1,maxnodes);

  neighbours = Imatrix(1,noknots,1,maxnodes);
  newrank = Ivector(1,noknots);
  oldrank = Ivector(1,noknots);
  newindx = Ivector(1,noknots);
  oldindx = Ivector(1,noknots);

  for(i=1;i<=noknots;i++)
    oldindx[i] = origindx[i];

  for(i=1;i<=noknots;i++)
    oldrank[origindx[i]] = i;

  minwidth = CalculateIndexwidth(data,TRUE,oldrank);
  if(info) printf("Indexwidth of the initial node order is %d.\n",minwidth);

  for(j=1;j<=noknots;j++) 
    for(i=1;i<=maxnodes;i++) 
      neighbours[j][i] = 0;


  if(info) printf("Initializing neighbours\n");

  for(j=1;j<=noelements;j++) {
    elemtype = data->elementtypes[j];
    nonodes = elemtype%100;
    nocands = 0;

    for(i=0;i<nonodes;i++) {
      ind = data->topology[j][i];

      if(elemtype == 404 || elemtype == 303) {
	nocands = 2;
	cands[1] = (i+1)%nonodes;
	cands[2] = (i+nonodes-1)%nonodes;
      }
      else if(elemtype == 808) {
	nocands = 3;
	if(i<4) {
	  cands[1] = (i+1)%4;
	  cands[2] = (i+3)%4;
	  cands[3] = i+4;
	}
	else {
	  cands[1] = (i-4+1)%4+4;
	  cands[2] = (i-4+3)%4+4;
	  cands[3] = i-4;
	}
      }

      for(k=1;k<=nocands;k++) { 
	ind2 = data->topology[j][cands[k]];
	for(l=1;l<=maxnodes;l++) {
	  if(neighbours[ind][l] == 0) break;
	  if(neighbours[ind][l] == ind2) ind2 = 0;
	}
	if(ind2) neighbours[ind][l] = ind2;
      }
    }
  }
  
#if 0
  for(j=1;j<=noknots;j++) {
    printf("neighbours[%d]= ",j);
    for(l=1;l<=maxnodes;l++) 
      printf("%d ",neighbours[j][l]);
    printf("\n");
  }
#endif

  if(info) printf("Reordering neighbours table\n");

  for(j=1;j<=noknots;j++) {

    nolocal = 0;
    dz = 0.0;

    for(l=1;l<=maxnodes;l++){
      if(ind = neighbours[j][l]) {
	nolocal++;
	localtmp[l] = ind;
	dx = data->x[l] - data->x[ind];
	dy = data->y[l] - data->y[ind];
	dz = data->z[l] - data->z[ind];
	localdist[l] = corder[0]*fabs(dx) + corder[1]*fabs(dy) + corder[2]*fabs(dz);
      }
    }    

    SortIndex(nolocal,localdist,localorder);
    
    for(l=1;l<=nolocal;l++) 
      neighbours[j][l] = localtmp[localorder[l]];

#if 0
    for(l=1;l<=nolocal;l++) 
      printf("j=%d l=%d dist=%.3le order=%d  %d\n",
	     j,l,localdist[l],localorder[l],neighbours[j][l]);
#endif
  } 



  for(iter=1;iter<=iterations;iter++) {

    if(info) printf("Optimal topology testing %d\n",iter);

    for(i=1;i<=noknots;i++) 
      newrank[i] = 0;
    
    ind = 0;
    indready = 1;
    
    do {
       if(indready > ind) {
	for(l=noknots;l>=1;l--)
	  if(j = oldindx[l]) break;
	if(info) printf("Starting over from node %d when ind=%d indready=%d\n",j,ind,indready);
      }
      else {
	j = newindx[indready] ;
      }

      for(l=1;ind2 = neighbours[j][l];l++) {
	if(ind2) { 
	  if(!newrank[ind2]) {
	    ind++;
	    newrank[ind2] = ind; 
	    newindx[ind]  = ind2;
	    oldindx[oldrank[ind2]] = 0;
	    oldrank[ind2] = 0;
	  }
	}
      }
      indready++;

    } while(ind < noknots);

    indexwidth = CalculateIndexwidth(data,TRUE,newrank);
    if(info) printf("Indexwidth of the suggested node order is %d.\n",indexwidth);
   
    for(i=1;i<=noknots;i++) 
      oldrank[i] = newrank[i];

    for(i=1;i<=noknots;i++) 
      oldindx[i] = newindx[i];

    if(indexwidth < minwidth) {
      for(i=1;i<=noknots;i++) 
	origindx[i] = newindx[i];
      minwidth = indexwidth;
    }
  }

  free_Ivector(cands,1,maxnodes);
  free_Ivector(localorder,1,maxnodes);
  free_Ivector(localtmp,1,maxnodes);
  free_Rvector(localdist,1,maxnodes);

  free_Imatrix(neighbours,1,noknots,1,maxnodes);
  free_Ivector(newrank,1,noknots);
  free_Ivector(oldrank,1,noknots);
  free_Ivector(newindx,1,noknots);
  free_Ivector(oldindx,1,noknots);
}

	


void ReorderElements(struct FemType *data,struct BoundaryType *bound,
		     int manual,Real corder[],int info)
{
  int i,j,k;
  int noelements,noknots,nonodes,length;
  int **newtopology=NULL,*newmaterial=NULL,*newelementtypes=NULL;
  int *indx=NULL,*revindx=NULL,*elemindx=NULL,*revelemindx=NULL;
  int oldnoknots, oldnoelements;
  Real *newx=NULL,*newy=NULL,*newz=NULL,*arrange=NULL;
  Real dx,dy,dz,cx,cy,cz,cbase;
  
  noelements = oldnoelements = data->noelements;
  noknots = oldnoknots = data->noknots;

  if(info) printf("Reordering %d knots and %d elements in %d-dimensions.\n",
		  noknots,noelements,data->dim);

  if(noelements > noknots)
    length = noelements;
  else
    length = noknots;

  arrange = Rvector(1,length);
  indx = Ivector(1,noknots);
  revindx = Ivector(1,noknots);
  elemindx = Ivector(1,noelements);
  revelemindx = Ivector(1,noelements);
  
  if(manual == 1) {
    cx = corder[0];
    cy = corder[1];
    cz = corder[2];
  }
  else {
    Real xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = xmax = data->x[1];
    ymin = ymax = data->y[1];
    zmin = zmax = data->z[1];

    for(i=1;i<=data->noknots;i++) {
      if(xmin > data->x[i]) xmin = data->x[i];
      if(xmax < data->x[i]) xmax = data->x[i];
      if(ymin > data->y[i]) ymin = data->y[i];
      if(ymax < data->y[i]) ymax = data->y[i];
      if(zmin > data->z[i]) zmin = data->z[i];
      if(zmax < data->z[i]) zmax = data->z[i];
    }
    dx = xmax-xmin;
    dy = ymax-ymin;
    dz = zmax-zmin;

    /* The second strategy seems to be better in many cases */
    cbase = 100.0;
    cx = pow(cbase,1.0*(dx>dy)+1.0*(dx>dz));
    cy = pow(cbase,1.0*(dy>dx)+1.0*(dy>dz));
    cz = pow(cbase,1.0*(dz>dx)+1.0*(dz>dx));

    corder[0] = cx;
    corder[1] = cy;
    corder[2] = cz;
  }

  if(info) printf("Ordering with (%.3lg*x + %.3lg*y + %.3lg*z)\n",cx,cy,cz);
  for(i=1;i<=noknots;i++) {
    arrange[i] = cx*data->x[i] + cy*data->y[i] + cz*data->z[i];  
  }
  SortIndex(noknots,arrange,indx);

  if(manual == 2) ReorderAutomatic(data,0,indx,corder,TRUE);

  for(i=1;i<=noknots;i++) 
    revindx[indx[i]] = i;


  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j]%100;
    arrange[j] = 0.0;
    for(i=0;i<nonodes;i++) {
      k = data->topology[j][i];
      arrange[j] += cx*data->x[k] + cy*data->y[k] + cz*data->z[k];
    }
  }

  SortIndex(noelements,arrange,elemindx);
  for(i=1;i<=noelements;i++) 
    revelemindx[elemindx[i]] = i;


#if 0
  for(i=1;i<=noknots;i++) 
    printf("i=%d  indx=%d  revi=%d  f=%.2lg\n",
	   i,indx[i],revindx[i],arrange[indx[i]]);
#endif  

  if(info) printf("Moving knots to new positions\n");
  newx = Rvector(1,data->noknots);
  newy = Rvector(1,data->noknots);
  newz = Rvector(1,data->noknots);

  for(i=1;i<=data->noknots;i++) {
    newx[i] = data->x[indx[i]];
    newy[i] = data->y[indx[i]];
    newz[i] = data->z[indx[i]];
  }

  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);

  data->x = newx;
  data->y = newy;
  data->z = newz;

  if(info) printf("Moving the elements to new positions\n");

  newtopology = Imatrix(1,noelements,0,data->maxnodes-1);
  newmaterial = Ivector(1,noelements);
  newelementtypes = Ivector(1,noelements);

  for(j=1;j<=noelements;j++) {
    newmaterial[j] = data->material[elemindx[j]];
    newelementtypes[j] = data->elementtypes[elemindx[j]];
    nonodes = newelementtypes[j]%100;
    for(i=0;i<nonodes;i++) {
      k = data->topology[elemindx[j]][i];
      newtopology[j][i] = revindx[k];
    }
  }

  data->material = newmaterial;
  data->elementtypes = newelementtypes;
  data->topology = newtopology;


  printf("Moving the parents of the boundary nodes.\n");
  for(j=0;j < MAXBOUNDARIES;j++) {

    if(!bound[j].created) continue;
    
    for(i=1; i <= bound[j].nosides; i++) {
      
      bound[j].parent[i] = revelemindx[bound[j].parent[i]];

      if(bound[j].parent2[i]) 
	bound[j].parent2[i] = revelemindx[bound[j].parent2[i]];
    }
  }

  i = CalculateIndexwidth(data,FALSE,indx);
  printf("Indexwidth of the new node order is %d.\n",i);

  free_Rvector(arrange,1,length);
  free_Ivector(indx,1,oldnoknots);
  free_Ivector(revindx,1,oldnoknots);
  free_Ivector(elemindx,1,oldnoelements);
  free_Ivector(revelemindx,1,oldnoelements);
}



int RemoveUnusedNodes(struct FemType *data,int info)
{
  int i,j;
  int noelements,noknots,nonodes,activeknots;
  int *indx;
  
  noelements = data->noelements;
  noknots = data->noknots;
  
  indx = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) indx[i] = 0;
  
  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j] % 100;
    for(i=0;i<nonodes;i++) 
      indx[ data->topology[j][i] ] = 1;
  }
  
  activeknots = 0;
  for(i=1;i<=noknots;i++) {
    if(indx[i]) {
      activeknots += 1;
      indx[i] = activeknots;
    }  
  }
  
  if( noknots == activeknots) {
    if(info) printf("All %d nodes were used by the mesh elements\n",noknots);
    return(1);
  }

  if(info) printf("Removing %d unused nodes (out of %d) from the mesh\n",noknots-activeknots,noknots);

  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j] % 100;
    for(i=0;i<nonodes;i++) 
      data->topology[j][i] = indx[ data->topology[j][i] ];
  }

  for(i=1;i<=noknots;i++) {
    j = indx[i];
    if(!j) continue;
    data->x[j] = data->x[i];
    data->y[j] = data->y[i];
    data->z[j] = data->z[i];
  }
  data->noknots = activeknots;
  
  free_Ivector(indx,1,noknots);

  return(0);
}




void RenumberBoundaryTypes(struct FemType *data,struct BoundaryType *bound,
			   int renumber, int bcoffset, int info)
{
  int i,j,k,doinit,isordered;
  int minbc=0,maxbc=0,**mapbc;
  int elemdim=0,elemtype=0,sideind[MAXNODESD1];
  int bctype,havename,isname;

  if(renumber) {
    if(0) printf("Renumbering boundary types\n");
    
    doinit = TRUE;
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;

      for(i=1;i<=bound[j].nosides;i++) {
	if(doinit) {
	  maxbc = minbc = bound[j].types[i];
	  doinit = FALSE;
	}
	maxbc = MAX(maxbc,bound[j].types[i]);
	minbc = MIN(minbc,bound[j].types[i]);     
      }
    }
    if(doinit) return;
    
    if(info) printf("Initial boundary interval [%d,%d]\n",minbc,maxbc);

    mapbc = Imatrix(minbc,maxbc,0,2);
    for(i=minbc;i<=maxbc;i++) 
      for(j=0;j<=2;j++)
	mapbc[i][j] = 0;
      
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      for(i=1;i<=bound[j].nosides;i++) {
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],data,sideind,&elemtype);
	if(!elemtype) printf("could not find boundary element: %d %d %d\n",i,j,bound[j].parent[i]);
	elemdim = GetElementDimension(elemtype);
	bctype = bound[j].types[i];
	
	if(0) printf("type and dim: %d %d %d\n",elemtype,elemdim,bctype);
       	
	mapbc[bctype][elemdim] += 1;
      }
    }

    if(0) {
      for(i=minbc;i<=maxbc;i++) 
	for(j=0;j<=2;j++)
	  if(mapbc[i][j]) printf("bc map1: %d %d\n",i,mapbc[i][j]);
    }
      
    j = 0;
    /* Give the larger dimension always a smaller BC type */
    isordered = TRUE;
    for(isname=1;isname>=0;isname--) {      
      for(elemdim=2;elemdim>=0;elemdim--) {
	for(i=minbc;i<=maxbc;i++) {
	  if(!mapbc[i][elemdim]) continue;

	  /* Give index 1st to the named entities, then to unnamed. */
	  havename = FALSE;
	  if(data->boundarynamesexist) {
	    if(i<MAXBCS) { 
	      if(data->boundaryname[i]) havename = TRUE;
	    }
	  }
	  if(havename != isname) break;

	  j++;
	  if(i == j) {
	    if( isname ) {
	      printf("BC index unaltered %d in %d %dD elements of %s\n",i,mapbc[i][elemdim],elemdim,data->boundaryname[i]); 
	    }
	    else {	      
	      printf("BC index unaltered %d in %d %dD elements\n",i,mapbc[i][elemdim],elemdim);
	    }
	  }
	  else {
	    isordered = FALSE;
	    if( isname ) {
	      printf("BC index changed %d -> %d in %d %dD elements of %s\n",i,j,mapbc[i][elemdim],elemdim,data->boundaryname[i]);
	    }
	    else {
	      printf("BC index changed %d -> %d in %d %dD elements\n",i,j,mapbc[i][elemdim],elemdim);
	    }
	  }
	  mapbc[i][elemdim] = j;
	}
      }
    }
      
    if(0) {
      for(i=minbc;i<=maxbc;i++) 
	for(j=0;j<=2;j++)
	  if(mapbc[i][j]) printf("bc map2: %d %d\n",i,mapbc[i][j]);
    }

    if(isordered) {
      if(info) printf("Numbering of boundary types is already ok\n");
    }
    else {
      if(info) printf("Mapping boundary types from [%d %d] to [%d %d]\n",minbc,maxbc,1,j);    
    
      for(j=0;j < MAXBOUNDARIES;j++) {
	if(!bound[j].created) continue;
	for(i=1;i<=bound[j].nosides;i++) {
	  GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],data,sideind,&elemtype);
	  elemdim = GetElementDimension(elemtype);
	  bound[j].types[i] = mapbc[bound[j].types[i]][elemdim];
	}
      }
      
      if(data->boundarynamesexist) {
	char *boundaryname0[MAXBCS];

	for(j=0;j<MAXBODIES;j++)
	  boundaryname0[j] = NULL;
	
	/* We need some temporal place is name mapping might not be unique */
	for(j=minbc;j<=MIN(maxbc,MAXBCS-1);j++) {
	  k = 0;
	  for(elemdim=2;elemdim>=0;elemdim--) {	    
	    k = mapbc[j][elemdim];
	    if(k) break;
	  }
	  if(k) {
	    if(data->boundaryname[j]) {
	      boundaryname0[j] = Cvector(0,MAXNAMESIZE);	    
	      strcpy(boundaryname0[j],data->boundaryname[j]);
	      free_Cvector(data->boundaryname[j],0,MAXNAMESIZE);
	      data->boundaryname[j] = NULL;
	    }
	  }
	}
	
	for(j=minbc;j<=MIN(maxbc,MAXBCS-1);j++) {
	  k = 0;
	  for(elemdim=2;elemdim>=0;elemdim--) {	    	   	    
	    k = mapbc[j][elemdim];
	    if(k) break;
	  }
	  if(k) {
	    if(boundaryname0[j]) {
	      if(!data->boundaryname[k]) 
		data->boundaryname[k] = Cvector(0,MAXNAMESIZE);
	      strcpy(data->boundaryname[k],boundaryname0[j]);
	    }
	  }
	}
      }
    }
    free_Imatrix(mapbc,minbc,maxbc,0,2);
  }

  if(bcoffset) {
    if(info) printf("Adding offset of %d to the BCs\n",bcoffset);
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      for(i=1;i<=bound[j].nosides;i++)
	bound[j].types[i] += bcoffset;
    }
    if(data->boundarynamesexist) {
      for(j=MAXBCS-bcoffset-1;j>=0;j--) {
	k = j+bcoffset;
	if(!data->boundaryname[k]) data->boundaryname[k] = Cvector(0,MAXNAMESIZE);
	strcpy(data->boundaryname[k],data->boundaryname[j]);
      }
    }
  }
}
  


void RenumberMaterialTypes(struct FemType *data,struct BoundaryType *bound,int info)
{     
  int i,j,k,noelements,doinit;
  int minmat=0,maxmat=0,*mapmat;
  
  if(0) printf("Setting new material types\n");
  
  noelements = data->noelements;
  if(noelements < 1) {
    printf("There are no elements to set!\n");
    return;
  }

  doinit = TRUE;
  for(j=1;j<=noelements;j++) {
    if(doinit) {
      maxmat = minmat = data->material[j];
      doinit = FALSE;
    }
    maxmat = MAX(maxmat,data->material[j]);
    minmat = MIN(minmat,data->material[j]);    
  }

  if(info) printf("Initial body interval [%d,%d]\n",minmat,maxmat);

  mapmat = Ivector(minmat,maxmat);
  for(i=minmat;i<=maxmat;i++) mapmat[i] = 0;
  
  for(j=1;j<=noelements;j++) 
    mapmat[data->material[j]] += 1;
  
  j = 0;
  for(i=minmat;i<=maxmat;i++) {
    if(mapmat[i]) {
      j++;
      if(i != j) printf("body index changed %d -> %d in %d elements\n",i,j,mapmat[i]); 
      mapmat[i] = j;
    }
  }  

  if(maxmat - minmat >= j || minmat != 1) { 
    if(info) printf("Mapping material types from [%d %d] to [%d %d]\n",
		    minmat,maxmat,1,j);    
    for(j=1;j<=noelements;j++) 
      data->material[j] = mapmat[data->material[j]];

    if(data->bodynamesexist) {
      if(info) printf("Mapping entity names to follow material indexes\n");
      for(j=minmat;j<=MIN(maxmat,MAXBODIES-1);j++) {
	k = mapmat[j];
	if(k) {
	  if(data->bodyname[j]) {
	    if(!data->bodyname[k]) data->bodyname[k] = Cvector(0,MAXNAMESIZE);
	    strcpy(data->bodyname[k],data->bodyname[j]);
	  }
	}
      }
    }
  }
  else {
    if(info) printf("Numbering of bodies is already ok\n");
  }
  
  free_Ivector(mapmat,minmat,maxmat);

  if(info) printf("Renumbering of material types completed!\n");
}



int RemoveLowerDimensionalBoundaries(struct FemType *data,struct BoundaryType *bound,int info)
{
  int i,j,noelements;
  int elemtype,maxelemdim,minelemdim,elemdim;
  int parent, side, sideind[MAXNODESD1],sideelemtype;
  int nosides, oldnosides,newnosides;
    
  if(info) printf("Removing lower dimensional boundaries\n");
  
  noelements = data->noelements;
  if(noelements < 1) return(1);

  elemtype = GetMaxElementType(data);

  maxelemdim = GetElementDimension(elemtype);

  if(info) printf("Maximum elementtype is %d and dimension %d\n",elemtype,maxelemdim);
  
  elemtype = GetMinElementType(data);
  minelemdim = GetElementDimension(elemtype);  
  if(info) printf("Minimum elementtype is %d and dimension %d\n",elemtype,minelemdim);

  /* Nothing to remove if the bulk mesh has 1D elements */
  if(minelemdim < 2) return(2);
  
  oldnosides = 0;
  newnosides = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    nosides = 0;
    if(!bound[j].created) continue;
    for(i=1;i<=bound[j].nosides;i++) {

      oldnosides++;
      parent =  bound[j].parent[i];

      side = bound[j].side[i];

      GetBoundaryElement(i,&bound[j],data,sideind,&sideelemtype);
      /* Old: GetElementSide(parent,side,1,data,sideind,&sideelemtype); */
      
      elemdim = GetElementDimension(sideelemtype);

      /* if(maxelemdim - elemdim > 1) continue; */
      /* This was changed as we want to maintain 1D BCs of a hybrid 2D/3D mesh. */
      if(minelemdim - elemdim > 1) continue;
      
      nosides++;      
      if(nosides == i) continue;

      bound[j].parent[nosides] = bound[j].parent[i];
      bound[j].parent2[nosides] = bound[j].parent2[i];
      bound[j].side[nosides] = bound[j].side[i];
      bound[j].side2[nosides] = bound[j].side2[i];
      bound[j].types[nosides] = bound[j].types[i];
    }
    bound[j].nosides = nosides;
    newnosides += nosides;
  }
  
  if(info) printf("Removed %d (out of %d) less than %dD boundary elements\n",
		  oldnosides-newnosides,oldnosides,minelemdim-1);
  return(0);
}  
  
int RemoveInternalBoundaries(struct FemType *data,struct BoundaryType *bound,int info)
{
  int i,j;
  int parent,parent2;
  int nosides,oldnosides,newnosides;
    
  if(info) printf("Removing internal boundaries\n");
  
  if( data->noelements < 1 ) return(1);

  oldnosides = 0;
  newnosides = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    nosides = 0;
    if(!bound[j].created) continue;
    for(i=1;i<=bound[j].nosides;i++) {

      oldnosides++;
      parent =  bound[j].parent[i];
      parent2 = bound[j].parent2[i];

      if( parent > 0 && parent2 > 0 ) continue;
      
      nosides++;      
      if(nosides == i) continue;

      bound[j].parent[nosides] = bound[j].parent[i];
      bound[j].parent2[nosides] = bound[j].parent2[i];
      bound[j].side[nosides] = bound[j].side[i];
      bound[j].side2[nosides] = bound[j].side2[i];
      bound[j].types[nosides] = bound[j].types[i];
    }
    bound[j].nosides = nosides;
    newnosides += nosides;
  }

  if(info) printf("Removed %d (out of %d) internal boundary elements\n",
		  oldnosides-newnosides,oldnosides);
  return(0);
}  


#if 0
static void FindEdges(struct FemType *data,struct BoundaryType *bound,
		      int material,int sidetype,int info)
{
  int i,j,side,identical,noelements,element;
  int noknots,nomaterials,nosides,newbound;
  int maxelementtype,maxedgenodes,elemedges,maxelemedges,edge,dosides;
  int **edgetable,sideind[MAXNODESD1],sideelemtype,allocated;
  int *indx;
  Real *arrange;
  
 
  newbound = 0;
  nomaterials = 0;
  maxelementtype = 0;
  noelements = data->noelements;

  printf("FindEdges: Finding edges of bulk elements of type %d\n",material);
  maxelementtype = GetMaxElementType(data);

  if(maxelementtype/100 > 4) {
    printf("FindEdges: Implemented only for 2D elements!\n");
    dosides = 0;
    return;
  } 

  if(maxelementtype/100 <= 2) maxedgenodes = 1;
  else if(maxelementtype/100 <= 4) maxedgenodes = 2;
  maxelemedges = maxelementtype/100;

  edgetable = Imatrix(1,maxelemedges*nomaterials,0,maxedgenodes+1);
  for(i=1;i<=maxelemedges*nomaterials;i++) 
    for(j=0;j<=maxedgenodes+1;j++) 
      edgetable[i][j] = 0;
  
  edge = 0;
  for(element=1;element<=noelements;element++) {
    if(data->material[element] != material) continue;
    
    elemedges = data->elementtypes[element]/100;
    
    for(side=0;side<elemedges;side++) {
      edge++;
      
      GetElementSide(element,side,1,data,sideind,&sideelemtype);
      edgetable[edge][maxedgenodes] = element;
      edgetable[edge][maxedgenodes+1] = side;
      
      if(maxedgenodes == 1) 
	edgetable[edge][0] = sideind[0];
      else if(maxedgenodes == 2) {
	if(sideind[0] > sideind[1]) {
	    edgetable[edge][0] = sideind[0];
	    edgetable[edge][1] = sideind[1];
	}
	else {
	  edgetable[edge][1] = sideind[0];
	  edgetable[edge][0] = sideind[1];
	}
      }
    }
  }

  noknots = edge;
  arrange = Rvector(1,noknots);
  for(i=1;i<=noknots;i++) 
    arrange[i] = 0.0;
  for(i=1;i<=noknots;i++) 
    arrange[i] = edgetable[i][0];
  indx = Ivector(1,noknots);

  SortIndex(noknots,arrange,indx);

  allocated = FALSE;

omstart:
  nosides = 0;

  for(i=1;i<=noknots;i++) {
    identical = FALSE;
    if(maxedgenodes == 1) {
      for(j=i+1;j<=noknots && edgetable[indx[i]][0] == edgetable[indx[j]][0];j++) 
	identical = TRUE;
      for(j=i-1;j>=1 && edgetable[indx[i]][0] == edgetable[indx[j]][0];j--) 
	identical = TRUE;
    }
    else if(maxedgenodes == 2) {
      for(j=i+1;j<=noknots && edgetable[indx[i]][0] == edgetable[indx[j]][0];j++) 
	if(edgetable[indx[i]][1] == edgetable[indx[j]][1]) 
	  identical = TRUE;
      for(j=i-1;j>=1 && edgetable[indx[i]][0] == edgetable[indx[j]][0];j--) 
	if(edgetable[indx[i]][1] == edgetable[indx[j]][1]) 
	  identical = TRUE;
    }

    if(identical) continue;
    nosides++;
    if(allocated) {
      bound[newbound].parent[nosides] = edgetable[indx[i]][maxedgenodes];
      bound[newbound].parent2[nosides] = 0;
      bound[newbound].side[nosides] = edgetable[indx[i]][maxedgenodes+1];
      bound[newbound].side2[nosides] = 0;
      bound[newbound].types[nosides] = sidetype;
    }
  }
  
  if(!allocated) {
    for(j=0;j < MAXBOUNDARIES && bound[j].created;j++); 
    newbound = j;
    AllocateBoundary(&bound[newbound],nosides);
    allocated = TRUE;  
    if(info) printf("Created boundary %d of type %d and size %d for material %d\n",
		    newbound,sidetype,nosides,material);
    goto omstart;
  }

  free_Ivector(indx,1,noknots);
  free_Imatrix(edgetable,1,maxelemedges*nomaterials,0,maxedgenodes+1);
}
#endif

static int CompareIndexes(int elemtype,int *ind1,int *ind2)
{
  int i,j,same,nosides,hits;
  
  hits = 0;
  nosides = elemtype / 100;
  for(i=0;i<nosides;i++) 
    for(j=0;j<nosides;j++) 
      if(ind1[i] == ind2[j]) hits++;

  same = (hits == nosides);
  return(same);
}



int FindNewBoundaries(struct FemType *data,struct BoundaryType *bound,
		      int *boundnodes,int suggesttype,int dimred,int info)
{
  int i,j,side,identical,element,lowerdim,dim,minedge,maxedge;
  int noelements,noknots,nonodes,nosides,newbound;
  int sideind[MAXNODESD1],sideind0[MAXNODESD1],sideelemtype,sideelemtype0,allocated;
  int noboundnodes,sameside,newtype,elemtype;

  newtype = 0;
  allocated = FALSE;
  dim = data->dim;
  if(dimred) 
    lowerdim = dim - dimred;
  else 
    lowerdim = dim-1;

  noknots = data->noknots;
  noelements = data->noelements;
  noboundnodes = 0;
  newbound = 0;
  maxedge = 0;
  minedge = 0;

  for(i=1;i<=noknots;i++) 
    if(boundnodes[i]) noboundnodes++;
  if(!noboundnodes) {
    printf("FindNewBoundaries: no nonzero entries in boundnodes vector!\n");
    return(1);
  }
  else {
    if(info) printf("There are %d nonzero entries in boundnodes vector!\n",noboundnodes);
  }

 omstart:

  nosides = 0;
  for(element=1;element<=noelements;element++) {
    
    elemtype = data->elementtypes[element];
    if(dim == 1) {
      minedge = 0;
      maxedge = elemtype/100 -1;
    }
    else if(dim == 2) {
      if(lowerdim == 1) {
	minedge = 0;
	maxedge = elemtype/100 -1;
      }
      else if(lowerdim == 0) {
	minedge = elemtype/100;
	maxedge = minedge + 1;
      }
    }    
    else if(dim == 3) {
      if(lowerdim == 2) {
	minedge = 0;
	if(elemtype/100 == 5) maxedge = 3;
	else if(elemtype/100 == 6 || elemtype/100 == 7) maxedge = 4;
	else if(elemtype/100 == 8) maxedge = 5;
      }
      else if(lowerdim == 1) {
	if(elemtype/100 == 8) {
	  minedge = 6;
	  maxedge = 17;
	}
	else if(elemtype/100 == 5) {
	  minedge = 4;
	  maxedge = 9;
	}
	else 
	  printf("FindNewBoundaries: not implemented for all 3d boundaries\n");
      }
      else if(lowerdim == 0) {
	if(elemtype/100 == 8) {
	  minedge = 18;
	  maxedge = 25;
	}
      }
    }

    for(side=minedge;side<=maxedge;side++) {

      GetElementSide(element,side,1,data,sideind,&sideelemtype);

      nonodes = sideelemtype % 100;
      identical = TRUE;
      for(i=0;i<nonodes;i++)
	if(!boundnodes[sideind[i]]) identical = FALSE;

      if(!identical) continue;
      nosides++;

      if(allocated) {
	for(i=1;i<nosides;i++) {
	  if(bound[newbound].parent2[i]) continue;

	  GetElementSide(bound[newbound].parent[i],bound[newbound].side[i],
			 1,data,sideind0,&sideelemtype0);
	  if(sideelemtype0 != sideelemtype) continue;
	  sameside = CompareIndexes(sideelemtype,sideind0,sideind);
	  if(sameside) {
	    bound[newbound].parent2[i] = element;
	    bound[newbound].side2[i] = side;
	    nosides--;
	    goto foundsameside;
	  }
	}

	bound[newbound].types[nosides] = newtype;
	bound[newbound].parent[nosides] = element;
	bound[newbound].side[nosides] = side;
	bound[newbound].types[nosides] = newtype;

      foundsameside:
	continue;
      }
    }
  }

  if(nosides) {
    if(!allocated) {
      newtype = suggesttype;
      for(j=0;j < MAXBOUNDARIES && bound[j].created;j++) {
	newbound = j;
	if(suggesttype) continue;
	for(i=1;i<=bound[j].nosides;i++) 
	  if(bound[j].types[i] > newtype) newtype = bound[j].types[i];
      }    
      newbound++;      
      if(!suggesttype) newtype++;

      AllocateBoundary(&bound[newbound],nosides);
      allocated = TRUE;  
      if(info) printf("Allocating for %d sides of boundary %d\n",nosides,newtype);
      goto omstart;
    }

    bound[newbound].nosides = nosides;
    if(info) printf("Found %d sides of dim %d to define boundary %d\n",nosides,lowerdim,newtype);
    
    for(i=1;i<=nosides;i++) {
      if(j = bound[newbound].parent2[i]) {
	if(bound[newbound].parent[i] > bound[newbound].parent2[i]) {
	  bound[newbound].parent2[i] = bound[newbound].parent[i];
	  bound[newbound].parent[i] = j;
	  j = bound[newbound].side2[i];
	  bound[newbound].side2[i] = bound[newbound].side[i];
	  bound[newbound].side[i] = j;
	}
      }
    }
  }
  else {
    if(lowerdim == 0) {
      printf("The nodes do not form a boundary!\n");
      return(2);
    }
    else {
      lowerdim--;
      printf("The nodes do not form a boundary, trying with %d-dimensional elements.\n",lowerdim);
      goto omstart;
    }
  }

  return(0);
}



int FindBulkBoundary(struct FemType *data,int mat1,int mat2,
		     int *boundnodes,int *noboundnodes,int info)
{
  int i,j,k;
  int nonodes,maxnodes,minnodes,material;
  Real ds,xmin=0.0,xmax=0.0,ymin=0.0,ymax=0.0,zmin=0.0,zmax=0.0,eps;
  int *visited,elemdim,*ind;
  Real *anglesum,dx1,dx2,dy1,dy2,dz1,dz2,ds1,ds2,dotprod;
  
  eps = 1.0e-4;
  *noboundnodes = 0;

  if(mat1 < 1 && mat2 < 1) {
    printf("FindBulkBoundary: Either of the materials must be positive\n");
    return(1);
  }
  else if(mat1 < 1) {
    i = mat1;
    mat1 = mat2;
    mat2 = i;
  }
  if(info) printf("Finding nodes between bulk elements of material %d and %d\n",mat1,mat2);

  visited = Ivector(1,data->noknots);
  for(i=1;i<=data->noknots;i++)
    visited[i] = 0;

  for(i=1;i<=data->noknots;i++)
    boundnodes[i] = 0;

  elemdim = 0;
  for(i=1;i<=data->noelements;i++) {
    material = data->material[i];
    if(material == mat1) {
      nonodes = data->elementtypes[i] % 100;
      k = data->elementtypes[i]/100;
      if(k > elemdim) elemdim = k;

      for(j=0;j<nonodes;j++) {
	k = data->topology[i][j];
	visited[k] += 1;
      }
    }
  }
  maxnodes = minnodes = visited[1];
  for(i=1;i<=data->noknots;i++) {
    if(visited[i] > maxnodes) maxnodes = visited[i];
    if(visited[i] < minnodes) minnodes = visited[i];
  }

  if(elemdim == 3 || elemdim == 4) {
    anglesum = Rvector(1, data->noknots);
    for(i=1;i<=data->noknots;i++)
      anglesum[i] = 0.0;

    for(i=1;i<=data->noelements;i++) {
      material = data->material[i];
      if(material == mat1) {
	nonodes = data->elementtypes[i]/100;
	ind = data->topology[i];

	if(nonodes == 3 || nonodes == 4) {
	  for(k=0;k<nonodes;k++) {
	    dx1 = data->x[ind[(k+1)%nonodes]] - data->x[ind[k]];
	    dy1 = data->y[ind[(k+1)%nonodes]] - data->y[ind[k]];
	    dz1 = data->z[ind[(k+1)%nonodes]] - data->z[ind[k]];
	    dx2 = data->x[ind[(k+nonodes-1)%nonodes]] - data->x[ind[k]];
	    dy2 = data->y[ind[(k+nonodes-1)%nonodes]] - data->y[ind[k]];
	    dz2 = data->z[ind[(k+nonodes-1)%nonodes]] - data->z[ind[k]];
	    ds1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
	    ds2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
	    dotprod = dx1*dx2 + dy1*dy2 + dz1*dz2;

	    anglesum[ind[k]] += acos(dotprod / (ds1*ds2));
	  }
	}
 
      }
    }
    j = 0;
    for(i=1;i<=data->noknots;i++) {    
      anglesum[i] /= 2.0 * FM_PI;
      if(anglesum[i] > 0.99) visited[i] = 0;
      if(anglesum[i] > 1.01) printf("FindBulkBoundary: surprisingly large angle %.3e in node %d\n",anglesum[i],i);
      if(visited[i]) j++;
    }
    if(0) printf("There are %d boundary node candidates\n",j);
    free_Rvector(anglesum,1,data->noknots);
  }

  else {
    for(i=1;i<=data->noknots;i++) 
      if(visited[i] == maxnodes) visited[i] = 0;
    
    if(maxnodes < 2) {
      printf("FindBulkBoundary: Nodes must belong to more than %d elements.\n",maxnodes);
      return(2);
    }
  }  

  if(mat2 == 0) {
    for(i=1;i<=data->noelements;i++) {
      material = data->material[i];
      if(material == mat1) continue;

      nonodes = data->elementtypes[i] % 100;
      for(j=0;j<nonodes;j++) {
	k = data->topology[i][j];
	boundnodes[k] += 1;
      }
    }
    for(k=1;k<=data->noknots;k++) {
      if(!visited[k]) 
	boundnodes[k] = 0;
      else if(visited[k] < boundnodes[k])
	boundnodes[k] = 0;
      else if(visited[k] + boundnodes[k] < maxnodes) 
	boundnodes[k] = 1;
      else 
	boundnodes[k] = 0;
    }
  }
  else if(mat2 == -10) {
    for(i=1;i<=data->noknots;i++) 
      if(visited[i]) boundnodes[i] = 1;
  }
  else if(mat2 == -11 || mat2 == -12 || mat2 > 0) {
    for(i=1;i<=data->noelements;i++) {
      material = data->material[i];

      if(material == mat1) continue;
      if(mat2 > 0 && material != mat2) continue;
      if(mat2 == -11 && material < mat1) continue;
      if(mat2 == -12 && material > mat1) continue;

      nonodes = data->elementtypes[i]%100;
      for(j=0;j<nonodes;j++) {
	k = data->topology[i][j];
	if(visited[k]) boundnodes[k] = 1;
      }
    }
  }
  else if(mat2 >= -2*data->dim && mat2 <= -1) {

    j = TRUE;
    for(i=1;i<=data->noknots;i++) 
      if(visited[i]) {
	if(j) {
	  xmax = xmin = data->x[i];
	  ymax = ymin = data->y[i];
	  zmax = zmin = data->z[i];
	  j = FALSE;
	}
	else {
	  if(data->x[i] > xmax) xmax = data->x[i];
	  if(data->x[i] < xmin) xmin = data->x[i];
	  if(data->y[i] > ymax) ymax = data->y[i];
	  if(data->y[i] < ymin) ymin = data->y[i];
	  if(data->z[i] > zmax) zmax = data->z[i];
	  if(data->z[i] < zmin) zmin = data->z[i];
	}
      }

    ds = (xmax-xmin)*(xmax-xmin) +
      (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin);
    
    ds = sqrt(ds);
    eps = 1.0e-5 * ds;

  
    for(i=1;i<=data->noknots;i++) 
      if(visited[i] < maxnodes && visited[i]) {
	
	if(data->dim == 1) {
	  if(mat2 == -1 && fabs(data->x[i]-xmin) < eps) boundnodes[i] = 1;
	  else if(mat2 == -2 && fabs(data->x[i]-xmax) < eps) boundnodes[i] = 1;
	}
	if(data->dim >= 2) {
	  if(mat2 == -1 && (fabs(data->y[i]-ymin) < eps)) boundnodes[i] = 1;
	  else if(mat2 == -3 && (fabs(data->y[i]-ymax) < eps)) boundnodes[i] = 1;
	  else if(mat2 == -4 && (fabs(data->x[i]-xmin) < eps)) boundnodes[i] = 1;
	  else if(mat2 == -2 && (fabs(data->x[i]-xmax) < eps)) boundnodes[i] = 1;
	}	  
	if(data->dim >= 3) {
	  if(mat2 == -5 && fabs(data->z[i]-zmin) < eps) boundnodes[i] = 1;
	  else if(mat2 == -6 && fabs(data->z[i]-zmax) < eps) boundnodes[i] = 1;
	}
      }
      
  }
  else {
    printf("FindBulkBoundary: unknown option %d for finding a side\n",mat2);
    return(2);
  }


  *noboundnodes = 0;
  for(i=1;i<=data->noknots;i++)
    if(boundnodes[i]) *noboundnodes += 1;

  if(info) printf("Located %d nodes at the interval between materials %d and %d\n",
		  *noboundnodes,mat1,mat2);

  free_Ivector(visited,1,data->noknots);
  return(0);
}



int FindBoundaryBoundary(struct FemType *data,struct BoundaryType *bound,int mat1,int mat2,
			 int *boundnodes,int *noboundnodes,int info)
{
  int i,j,k,l;
  int hits,nonodes,nocorners,maxnodes,minnodes,elemtype,material,bounddim;
  Real ds,xmin=0.0,xmax=0.0,ymin=0.0,ymax=0.0,zmin=0.0,zmax=0.0;
  Real eps,dx1,dx2,dy1,dy2,dz1,dz2,ds1,ds2,dotprod;
  Real *anglesum=NULL;
  int *visited,sideind[MAXNODESD2],elemind[MAXNODESD2];
  
  eps = 1.0e-4;
  *noboundnodes = 0;

  if(mat1 < 1 && mat2 < 1) {
    printf("FindBoundaryBoundary: Either of the boundaries must be positive\n");
    return(1);
  }
  else if(mat1 < 1) {
    i = mat1;
    mat1 = mat2;
    mat2 = i;
  }
  if(info) printf("Finding nodes between boundary elements of type %d and %d\n",mat1,mat2);

  visited = Ivector(1,data->noknots);
  for(i=1;i<=data->noknots;i++)
    visited[i] = 0;

  for(i=1;i<=data->noknots;i++)
    boundnodes[i] = 0;
  
  bounddim = 0;
  /* Set a tag to all nodes that are part of the other boundary */
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    for(i=1; i <= bound[j].nosides; i++) {
      
      if(bound[j].types[i] == mat1) {
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,sideind,&elemtype);

	nonodes = elemtype % 100;
	nocorners = elemtype / 100;
	
	for(k=0;k<nocorners;k++) 
	  visited[sideind[k]] += 1;
	for(k=nocorners;k<nonodes;k++)
	  visited[sideind[k]] -= 1;
	
	if(nocorners == 3 || nocorners == 4) {
	  if(bounddim < 2) {
	    anglesum = Rvector(1, data->noknots);
	    for(k=1;k<=data->noknots;k++)
	      anglesum[k] = 0.0;	    
	    bounddim = 2;
	  }
	  nonodes = nocorners;
	  for(k=0;k<nonodes;k++) {
	    dx1 = data->x[sideind[(k+1)%nonodes]] - data->x[sideind[k]];
	    dy1 = data->y[sideind[(k+1)%nonodes]] - data->y[sideind[k]];
	    dz1 = data->z[sideind[(k+1)%nonodes]] - data->z[sideind[k]];
	    dx2 = data->x[sideind[(k+nonodes-1)%nonodes]] - data->x[sideind[k]];
	    dy2 = data->y[sideind[(k+nonodes-1)%nonodes]] - data->y[sideind[k]];
	    dz2 = data->z[sideind[(k+nonodes-1)%nonodes]] - data->z[sideind[k]];
	    ds1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
	    ds2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
	    dotprod = dx1*dx2 + dy1*dy2 + dz1*dz2;
	    
	    anglesum[sideind[k]] += acos(dotprod / (ds1*ds2));
	  }
	}

      } 
    }
  }
  
  maxnodes = minnodes = abs(visited[1]);
  for(i=1;i<=data->noknots;i++) {
    j = abs( visited[i] );
    maxnodes = MAX( j, maxnodes );
    minnodes = MIN( j, minnodes );
  }
  if(info) printf("There are from %d to %d hits per node\n",minnodes,maxnodes);
  if(maxnodes < 2) {
    printf("FindBulkBoundary: Nodes must belong to more than %d elements.\n",maxnodes);
    return(2);
  }

  if(bounddim == 2) {
    /* For corner nodes eliminate the ones with full angle */
    for(i=1;i<=data->noknots;i++) {
      anglesum[i] /= 2.0 * FM_PI;
      if(anglesum[i] > 0.99) visited[i] = 0;
      if(anglesum[i] > 1.01) printf("FindBulkBoundary: surprisingly large angle %.3e in node %d\n",anglesum[i],i);
    }
    free_Rvector(anglesum,1,data->noknots);

    /* For higher order nodes eliminate the ones with more than one hits */
    k = 0;
    for(i=1;i<=data->noknots;i++) {
      if(visited[i] == -1) 
	visited[i] = 1;
      else if(visited[i] < -1) {
	k++;
	visited[i] = 0;
      }
    }    
    if(k && info) printf("Removed %d potential higher order side nodes from the list.\n",k);
  }

  if(bounddim == 1) {
    if(visited[i] == maxnodes || visited[i] < 0) visited[i] = 0;      
  }

  /* Neighbour to anything */
  if(mat2 == 0) {
    for(k=1;k<=data->noknots;k++) 
      if(visited[k]) 
	boundnodes[k] = 1;
  }
  /* Neighbour to other BCs */
  else if(mat2 == -11 || mat2 == -12 || mat2 == -10 || mat2 > 0) {
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      for(i=1; i <= bound[j].nosides; i++) {
	
	material = bound[j].types[i];
	
	if(material == mat1) continue;
	if(mat2 > 0 && material != mat2) continue;
	if(mat2 == -11 && material < mat1) continue;
	if(mat2 == -12 && material > mat1) continue;

	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,sideind,&elemtype);
	nonodes = elemtype%100;
	for(k=0;k<nonodes;k++) {
	  l = sideind[k];
	  if(visited[l]) boundnodes[l] = 1;
	}
      }
    }
  }

  /* Neighbour to major coordinate directions */
  else if(mat2 >= -2*data->dim && mat2 <= -1) {
    
    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      for(i=1; i <= bound[j].nosides; i++) {
	
	material = bound[j].types[i];
	if(material != mat1) continue;
	
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,sideind,&elemtype);	
	nonodes = elemtype%100;
	
	hits =  0;
	for(k=0;k<nonodes;k++) {
	  l = sideind[k];
	  if(visited[l] < maxnodes) hits++;
	}
	if(!hits) continue;
	
	l = sideind[0];
	xmax = xmin = data->x[l];
	ymax = ymin = data->y[l];
	zmax = zmin = data->z[l];
	
	for(k=1;k<nonodes;k++) {
	  l = sideind[k];
	  if(data->x[l] > xmax) xmax = data->x[l];
	  if(data->x[l] < xmin) xmin = data->x[l];
	  if(data->y[l] > ymax) ymax = data->y[l];
	  if(data->y[l] < ymin) ymin = data->y[l];
	  if(data->z[l] > zmax) zmax = data->z[l];
	  if(data->z[l] < zmin) zmin = data->z[l];
	}

	ds = (xmax-xmin)*(xmax-xmin) +
	  (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin);
	ds = sqrt(ds);
	eps = 1.0e-3 * ds;
	
	for(k=0;k<nonodes;k++) {
	  elemind[k] = 0;
	  l = sideind[k];
	  if(!visited[l]) continue; 
	  
	  if(data->dim == 1) {
	    if(mat2 == -1 && fabs(data->x[l]-xmin) < eps) boundnodes[l] = 1;
	    else if(mat2 == -2 && fabs(data->x[l]-xmax) < eps) boundnodes[l] = 1;
	  }
	  if(data->dim >= 2) {
	    if(mat2 == -1 && (fabs(data->y[l]-ymin) < eps)) elemind[l] = 1;
	    else if(mat2 == -3 && (fabs(data->y[l]-ymax) < eps)) elemind[l] = 1;
	    else if(mat2 == -4 && (fabs(data->x[l]-xmin) < eps)) elemind[l] = 1;
	    else if(mat2 == -2 && (fabs(data->x[l]-xmax) < eps)) elemind[l] = 1;
	  }	  
	  if(data->dim >= 3) {
	    if(mat2 == -5 && fabs(data->z[l]-zmin) < eps) elemind[l] = 1;
	    else if(mat2 == -6 && fabs(data->z[l]-zmax) < eps) elemind[l] = 1;
	  }
	}
      
	if(data->dim > 1) {
	  hits = 0;
	  for(k=0;k<nonodes;k++)
	    hits += elemind[k];
	  
	  if(hits > 1) for(k=0;k<nonodes;k++) 
	    if(elemind[k]) boundnodes[sideind[k]] = 1;
	}
      }
    }
  }
  else {
    printf("FindBoundaryBoundary: unknown option %d for finding a side\n",mat2);
    return(2);
  }
    
  *noboundnodes = 0;
  for(i=1;i<=data->noknots;i++)
    if(boundnodes[i]) *noboundnodes += 1;

  if(info) printf("Located %d nodes at the interval between boundaries %d and %d\n",
		  *noboundnodes,mat1,mat2);

  free_Ivector(visited,1,data->noknots);
  return(0);
}



int IncreaseElementOrder(struct FemType *data,int info)
{
  int i,j,side,element,maxcon,con,newknots,ind,ind2;
  int noelements,noknots,nonodes,maxnodes,maxelemtype,hit,node;
  int elemtype,stat;
  int **newnodetable=NULL,inds[2],**newtopo=NULL;
  Real *newx=NULL,*newy=NULL,*newz=NULL;
  
  if(info) printf("Trying to increase the element order of current elements\n");

  CreateNodalGraph(data,FALSE,info);

  noknots = data->noknots;
  noelements = data->noelements;
  maxcon = data->nodalmaxconnections;
  maxnodes = 0;

  newnodetable = Imatrix(0,maxcon-1,1,noknots);
  for(i=1;i<=noknots;i++) 
    for(j=0;j<maxcon;j++) 
      newnodetable[j][i] = 0;

  newknots = 0;
  for(i=1;i<=noknots;i++) {
    for(j=0;j<maxcon;j++) {
      con = data->nodalgraph[j][i];
      if(con > i) {
	newknots++;
	newnodetable[j][i] = noknots + newknots;
      }
    }
  }

  if(info) printf("There will be %d new nodes in the elements\n",newknots);

  newx = Rvector(1,noknots+newknots);
  newy = Rvector(1,noknots+newknots);
  newz = Rvector(1,noknots+newknots);


  for(i=1;i<=noknots;i++) {
    newx[i] = data->x[i];
    newy[i] = data->y[i];
    newz[i] = data->z[i];
  }
  for(i=1;i<=noknots;i++) {
    for(j=0;j<maxcon;j++) {
      con = data->nodalgraph[j][i];
      ind = newnodetable[j][i];
      if(con && ind) {
	newx[ind] = 0.5*(data->x[i] + data->x[con]);
	newy[ind] = 0.5*(data->y[i] + data->y[con]);
	newz[ind] = 0.5*(data->z[i] + data->z[con]);
      }
    }
  }  
    

  maxelemtype = GetMaxElementType(data);

  if(maxelemtype <= 303) 
    maxnodes = 6;
  else if(maxelemtype == 404) 
    maxnodes = 8;
  else if(maxelemtype == 504) 
    maxnodes = 10;
  else if(maxelemtype == 605) 
    maxnodes = 13;
  else if(maxelemtype == 706) 
    maxnodes = 15;
  else if(maxelemtype == 808) 
    maxnodes = 20;
  else {
    printf("Not implemented for elementtype %d\n",maxelemtype);
    bigerror("IncreaseElementOrder: Cannot continue the subroutine");
  }

  if(info) printf("New leading elementtype is %d\n",100*(maxelemtype/100)+maxnodes);

  newtopo = Imatrix(1,noelements,0,maxnodes-1);
    
  for(element=1;element<=noelements;element++) {
    elemtype = data->elementtypes[element];
    for(i=0;i<elemtype%100;i++)
      newtopo[element][i] = data->topology[element][i];
  }

 
  for(element=1;element<=data->noelements;element++) {
    elemtype = data->elementtypes[element];

    nonodes = data->elementtypes[element] % 100;
    for(side=0;;side++) {
      hit = GetElementGraph(element,side,data,inds);

      if(!hit) break;
      if(inds[0] > inds[1]) {
	ind = inds[1];
	ind2 = inds[0];
      }
      else {
	ind = inds[0];
	ind2 = inds[1];
      }
      for(j=0;j<maxcon;j++) {
	con = data->nodalgraph[j][ind];

	if(con == ind2) {
	  node = newnodetable[j][ind];
	  newtopo[element][nonodes+side] = node;
	}
      }
    }

    elemtype = 100*(elemtype/100)+nonodes+side;
    data->elementtypes[element] = elemtype;
  }

  stat = DestroyNodalGraph(data,info);

  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);
  free_Imatrix(data->topology,1,data->noelements,0,data->maxnodes);
  free_Imatrix(newnodetable,0,maxcon-1,1,noknots);

  data->x = newx;
  data->y = newy;
  data->z = newz;
  data->topology = newtopo;

  data->noknots += newknots;
  data->maxnodes = maxnodes;

  if(info) printf("Increased the element order from 1 to 2\n");

  return(0);
}



static void CylindricalCoordinateTransformation(struct FemType *data,Real r1,Real r2,
						int rectangle)
{
  int i,j,j2,ind1,ind2,nonodes1;
  Real x,y,r,f,z,q,x2,y2,z2,dx,dy,dz,eps,mult;
  int hits,trials,tests;
  int candidates,*candidatelist=NULL,*indx=NULL;

  if(rectangle) {
    printf("Rectangular geometry with r1=%.4lg for %d nodes.\n",
	   r1,data->noknots);
  }
  else {
    printf("Cylindrical geometry with r1=%.4lg r2=%.4lg for %d nodes.\n",
	 r1,r2,data->noknots);
  }
  

  for(i=1;i<=data->noknots;i++) {
    r = data->x[i];
    z = data->y[i];
    f = data->z[i];

    data->z[i] = z;

    if(r >= r2) {
      data->x[i] = cos(f)*r;
      data->y[i] = sin(f)*r;
    }
    else if(r <= r2) {

      mult = r/r1;

      if(r > r1) {
	q = (r-r1)/(r2-r1);
	r = r1;
      }
      else {
	q = -1.0;
      }

      if(f <= 0.25*FM_PI) {
	data->x[i] = r;
	data->y[i] = r1*4*(f-0.00*FM_PI)/FM_PI;
      }
      else if(f <= 0.75*FM_PI) {
	data->y[i] = r;
	data->x[i] = -r1*4*(f-0.5*FM_PI)/FM_PI;
      }
      else if(f <= 1.25*FM_PI) {
	data->x[i] = -r;
	data->y[i] = -r1*4*(f-1.0*FM_PI)/FM_PI;
      }
      else if(f <= 1.75*FM_PI){
	data->y[i] = -r;
	data->x[i] = r1*4*(f-1.5*FM_PI)/FM_PI;
      }
      else {
	data->x[i] = r;
	data->y[i] = r1*4*(f-2.0*FM_PI)/FM_PI;
      }

      if(!rectangle && q > 0.0) {
	data->x[i] = (1-q)*data->x[i] + q*cos(f)*r2;
	data->y[i] = (1-q)*data->y[i] + q*sin(f)*r2;
      }
      else if(rectangle && mult > 1.0) {
	data->y[i] *= mult;
	data->x[i] *= mult;
      }
    } /* r <= r2 */
  }

  eps = 1.0e-3 * data->minsize;

  candidates = 0;
  candidatelist = Ivector(1,data->noknots);
  indx = Ivector(1,data->noknots);

  for(i=1;i<=data->noknots;i++) 
    indx[i] = 0;

  for(j=1;j<=data->noelements;j++) {
    nonodes1 = data->elementtypes[j]%100;
    for(i=0;i<nonodes1;i++) {
      ind2 = data->topology[j][i];
      indx[ind2] = ind2;
    }
  }
      

  for(i=1;i<=data->noknots;i++) {
    if(!indx[i]) continue;

    x = data->x[i];
    y = data->y[i];
    if(fabs(y)  > r1+eps) continue;
    if((fabs(x) > r1+eps)  && (fabs(y) > eps) ) continue;
    if((fabs(x) > eps) && (fabs(y) > eps) &&
       (fabs(fabs(x)-r1) > eps) && (fabs(fabs(y)-r1) > eps)) continue;

    candidates++;
    candidatelist[candidates] = i;
  }
  printf("%d/%d candidates for duplicate nodes.\n",candidates,data->noknots);

  hits = tests = trials = 0;
  for(j=1;j<=candidates;j++) {
    ind1 = indx[candidatelist[j]];
    x = data->x[ind1];
    y = data->y[ind1];
    z = data->z[ind1];
    
    for(j2=j+1;j2<=candidates;j2++) {
      ind2 = indx[candidatelist[j2]];

      x2 = data->x[ind2];
      y2 = data->y[ind2];
      z2 = data->z[ind2];

      dx = x-x2;
      dy = y-y2;
      dz = z-z2;

      tests++;
      if(dx*dx + dy*dy + dz*dz < eps*eps) {
	if(ind2 != ind1) {
	  indx[candidatelist[j2]] = ind1;
	  hits++;
	}
      }
    }
  }
  printf("Found %d double nodes in %d tests.\n",hits,tests);

  for(j=1;j<=data->noelements;j++) {
    nonodes1 = data->elementtypes[j]%100;
    for(i=0;i<nonodes1;i++) {
      ind2 = data->topology[j][i];
      if(ind2 != indx[ind2]) {
	trials++;
	data->topology[j][i] = indx[ind2];
      }
    }
  }
  free_Ivector(indx,1,data->noknots);
  free_Ivector(candidatelist,1,data->noknots);

  printf("Eliminated %d nodes from topology.\n",trials);
}


static void CylindricalCoordinateImprove(struct FemType *data,Real factor,
					 Real r1,Real r2)
{
  int i;
  Real x,y,r,q,q2,c,cmin,cmax,eps;

  printf("Cylindrical coordinate for r1=%.4lg and r2=%.4lg.\n",r1,r2);

  eps = 1.0e-10;
 
  cmin = 1.0/(3.0-sqrt(3.));
  cmax = 1.0;

  if(factor > 1.0) factor=1.0;
  else if(factor < 0.0) factor=0.0;

  c = cmin+(cmax-cmin)*factor;

  if(fabs(c-1.0) < eps) return;

  printf("Improving cylindrical mesh quality r1=%.4lg, r2=%.4lg and c=%.4lg\n",r1,r2,c);

  for(i=1;i<=data->noknots;i++) {
    x = data->x[i];
    y = data->y[i];

    r = sqrt(x*x+y*y);
    if(r >= r2) continue;
    if(r < eps) continue;
    
    if(fabs(x) <= r1+eps && fabs(y) <= r1+eps) {
      if(fabs(x) < fabs(y)) {
	q = fabs(x/y);
	data->x[i] = (c*q+(1.-q))*x;
	data->y[i] = (c*q+(1.-q))*y;
      }
      else {
	q = fabs(y/x);
	data->x[i] = (c*q+(1.-q))*x;
	data->y[i] = (c*q+(1.-q))*y;
      }
    }
    else {
      if(fabs(x) < fabs(y)) {
	q = fabs(x/y);
	q2 = (fabs(y)-r1)/(r2*fabs(y/r)-r1);
	data->x[i] = (c*q+(1.-q)) *x*(1-q2) + q2*x;
	data->y[i] = (c*q+(1.-q)) *(1-q2)*y + q2*y;
      }
      else {
	q = fabs(y/x);
	q2 = (fabs(x)-r1)/(r2*fabs(x/r)-r1);
	data->x[i] = (c*q+(1.-q))*(1-q2)*x + q2*x;
	data->y[i] = (c*q+(1.-q))*y*(1-q2) + q2*y;
      }
    }
  } 
}


void CylindricalCoordinateCurve(struct FemType *data,
				Real zet,Real rad,Real angle)
{
  int i;
  Real x,y,z;
  Real z0,z1,f,f0,z2,x2,r0;
  
  printf("Cylindrical coordinate curve, zet=%.3lg  rad=%.3lg  angle=%.3lg\n",
	 zet,rad,angle);
  
  r0 = rad;
  f0 = FM_PI*(angle/180.);
  z0 = zet;
  z1 = z0+r0*f0;
  
  for(i=1;i<=data->noknots;i++) {

    if(data->dim == 2) {
      z = data->x[i];
      x = data->y[i];
    }
    else {
      x = data->x[i];
      y = data->y[i];
      z = data->z[i];
    }    

    if(z <= z0) continue;
    
    if(z >= z1) {
      z2 = z0 + sin(f0)*(r0+x) + cos(f0)*(z-z1);
      x2 = (cos(f0)-1.0)*r0 + cos(f0)*x - sin(f0)*(z-z1);
    }
    else {
      f = (z-z0)/r0;
      z2 = z0 + sin(f)*(r0+x);
      x2 = (cos(f)-1.0)*r0 + cos(f)*x;
    }          

    if( data->dim == 2) {
      data->x[i] = z2;
      data->y[i] = x2;      
    }
    else {
      data->z[i] = z2;
      data->x[i] = x2;      
    }

  }
}


void SeparateCartesianBoundaries(struct FemType *data,struct BoundaryType *bound,int info)
{
  int i,j,k,l,type,maxtype,addtype,elemsides,totsides,used,hit;
  int sideelemtype,sideind[MAXBOUNDARIES];
  Real x,y,z,sx,sy,sz,sxx,syy,szz,dx,dy,dz;
  Real bclim[MAXBOUNDARIES];
  int bc[MAXBOUNDARIES],bcdim[MAXBOUNDARIES];
  Real eps=1.0e-4;

  maxtype = 0;
  totsides = 0;
  for(j=0;j<MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    if(!bound[j].nosides) continue;

    for(i=1;i<=bound[j].nosides;i++) {    
      totsides++;
      for(k=1;k<=bound[j].nosides;k++)
	if(maxtype < bound[j].types[k]) maxtype = bound[j].types[k];
    }
  }

  if(info) {
    printf("Maximum boundary type is %d\n",maxtype);
    printf("Number of boundaries is %d\n",totsides);
  }
  addtype = maxtype;

  for(type=1;type<=maxtype;type++) {

    for(i=0;i<MAXBOUNDARIES;i++)
      bclim[i] = 0.0;
    for(i=0;i<MAXBOUNDARIES;i++)
      bc[i] = bcdim[i] = 0;
    used = FALSE;

    for(j=0;j<MAXBOUNDARIES;j++) {

      if(!bound[j].created) continue;
      if(!bound[j].nosides) continue;

      for(k=1;k<=bound[j].nosides;k++) {

	if(bound[j].types[k] != type) continue;
	GetElementSide(bound[j].parent[k],bound[j].side[k],bound[j].normal[k],
		       data,sideind,&sideelemtype);
	
	sx = sy = sz = 0.0;
	sxx = syy = szz = 0.0;
	elemsides = sideelemtype%100;
	
	/* Compute the variance within each axis */
	for(l=0;l<elemsides;l++) {
	  x = data->x[sideind[l]];
	  y = data->y[sideind[l]];
	  z = data->z[sideind[l]];
	  sx += x;
	  sy += y;
	  sz += z;
	  sxx += x*x;
	  syy += y*y;
	  szz += z*z;
	}
	sx /= elemsides;
	sy /= elemsides;
	sz /= elemsides;
	sxx /= elemsides;
	syy /= elemsides;
	szz /= elemsides;
	dx = sqrt(sxx-sx*sx);
	dy = sqrt(syy-sy*sy);
	dz = sqrt(szz-sz*sz);

	if(sideelemtype < 300 && dz < eps) {

	  if(dx < eps * dy) {
	    hit = FALSE;
	    for(i=0;i<MAXBOUNDARIES && bcdim[i];i++) { 
	      if(bcdim[i] == 1 && fabs(bclim[i]-sx) < eps*fabs(dy)) {
		bound[j].types[k] = bc[i];
		hit = TRUE;
		break;
	      }
	    }	    
	    
	    if(!hit) {
	      if(used) {
		addtype++;
		printf("Adding new BC %d in Y-direction\n",addtype);
		bc[i] = addtype;
		bound[j].types[k] = bc[i];
	      }
	      else {
		bc[i] = bound[j].types[k];
	      }
	      bcdim[i] = 1;
	      bclim[i] = sx;
	      used = TRUE;
	    }
	  }

	  if(dy < eps * dx) {
	    hit = FALSE;
	    for(i=0;i<MAXBOUNDARIES && bcdim[i];i++) { 
	      if(bcdim[i] == 2 && fabs(bclim[i]-sy) < eps*fabs(dx)) {
		bound[j].types[k] = bc[i];
		hit = TRUE;
		break;
	      }
	    }
	    if(!hit) {
	      if(used) {
		addtype++;
		printf("Adding new BC %d in X-direction\n",addtype);
		bc[i] = addtype;
		bound[j].types[k] = bc[i];
	      }
	      else {
		bc[i] = bound[j].types[k];
	      }
	      bcdim[i] = 2;
	      bclim[i] = sy;
	      used = TRUE;
	    }
	  }
	}
	else {
	  if(dx < eps*dy && dx < eps*dz) {
	  }
	  else if(dy < eps*dx && dy < eps*dz) {
	  }
	}
      }
    }
  }
}




void SeparateMainaxisBoundaries(struct FemType *data,struct BoundaryType *bound)
{
  int i,j,k,l,maxtype,addtype,elemsides;
  int sideelemtype,sideind[MAXNODESD1];
  int axistype[4],axishit[4],axissum,axismax,done;
  Real x,y,z,sx,sy,sz,sxx,syy,szz,dx,dy,dz;
  Real eps=1.0e-6;

  maxtype = 0;
  addtype = 0;

  for(j=0;j<data->noboundaries;j++) {

    if(!bound[j].created) continue;
    if(!bound[j].nosides) continue;

    for(i=1;i<=bound[j].nosides;i++) {    
      for(k=1;k<=bound[j].nosides;k++)
	if(maxtype < bound[j].types[k]) maxtype = bound[j].types[k];
    }
  }
  printf("Maximum boundary type is %d\n",maxtype);

#if 0
  for(j=0;j<data->noboundaries;j++) {
    if(!bound[j].created) continue;
    if(!bound[j].nosides) continue;
    if(bound[j].type) {
      bound[j].types = Ivector(1,bound[j].nosides);
      for(k=1;k<=bound[j].nosides;k++)
	bound[j].types[k] = bound[j].type;
      bound[j].type = 0;
    }
  }
#endif

  for(j=0;j<data->noboundaries;j++) {
    if(!bound[j].created) continue;
    if(!bound[j].nosides) continue;

    for(k=0;k<4;k++) axishit[k] = 0;

    done = 0;
    
  omstart:

    for(k=1;k<=bound[j].nosides;k++) {

      GetElementSide(bound[j].parent[k],bound[j].side[k],bound[j].normal[k],
		     data,sideind,&sideelemtype);

      sx = sy = sz = 0.0;
      sxx = syy = szz = 0.0;
      elemsides = sideelemtype%100;

      /* Compute the variance within each axis */
      for(l=0;l<elemsides;l++) {
	x = data->x[sideind[l]];
	y = data->y[sideind[l]];
	z = data->z[sideind[l]];
	sx += x;
	sy += y;
	sz += z;
	sxx += x*x;
	syy += y*y;
	szz += z*z;
      }
      sx /= elemsides;
      sy /= elemsides;
      sz /= elemsides;
      sxx /= elemsides;
      syy /= elemsides;
      szz /= elemsides;
      dx = sqrt(sxx-sx*sx);
      dy = sqrt(syy-sy*sy);
      dz = sqrt(szz-sz*sz);
	 
      if(dx < eps*dy && dx < eps*dz) {
	if(sx > 0.0) {
	  if(done) {
	    if(axistype[0]) bound[j].types[k] = maxtype + axistype[0];
	  }
	  else
	    axishit[0] += 1;
	}
	if(sx < 0.0) {
	  if(done) {
	    if(axistype[1]) bound[j].types[k] = maxtype + axistype[1];
	  }
	  else
	    axishit[1] += 1;
	}
      }
      else if(dy < eps*dx && dy < eps*dz) {
	if(sy > 0.0) {
	  if(done) {
	    if(axistype[2]) bound[j].types[k] = maxtype + axistype[2];
	  }
	  else
	    axishit[2] += 1;
	}
	if(sy < 0.0) {
	  if(done) {
	    if(axistype[3]) bound[j].types[k] = maxtype + axistype[3];
	  }
	  else
	    axishit[3] += 1;
	}
      }
    }

    /* All this is done to select the sidetype appropriately */
    if(!done) {
      axissum = 0;
      axismax = 0;
      
      for(k=0;k<4;k++) {
	axissum += axishit[k];
	if(axishit[k]) addtype++;
      }

      if(axissum) {
	for(k=0;k<4;k++) {
	  axismax = 0;
	  for(l=0;l<4;l++) {
	    if(axishit[l] > axishit[axismax])
	      axismax = l;
	  }
	  axistype[axismax] = k+1;
	  axishit[axismax] = -(k+1);
	}
	
	if(axissum == bound[j].nosides) {
	  for(k=0;k<4;k++) 
	    axistype[k] -= 1;
	  addtype--;
	}
	
	if(addtype) {
	  printf("Separating %d rectangular boundaries from boundary %d.\n",addtype,j);
	  done = 1;
	  goto omstart;
	}
	else done = 0;
      }
    }
    maxtype += addtype;
  }
}


void CreateKnotsExtruded(struct FemType *dataxy,struct BoundaryType *boundxy,
			 struct GridType *grid,
			 struct FemType *data,struct BoundaryType *bound,
			 int info)
/* Create mesh from 2D mesh either by extrusion or by rotation. 
   Also create the additional boundaries using automated numbering. */
{
#define MAXNEWBC 200
  int i,j,k,l,m,n,knot0,knot1,knot2,elem0,size,kmax,noknots,origtype;
  int nonodes3d,nonodes2d,bclevel,bcset;
  int cellk,element,level,side,parent,parent2,layers,elemtype,material_too_large;
  int material,material2,ind1,ind2;
  int *indx=NULL,*topo=NULL;
  int sideelemtype,sideind[MAXNODESD1],sidetype,minsidetype,maxsidetype,cummaxsidetype,newbounds;
  int refmaterial1[MAXNEWBC],refmaterial2[MAXNEWBC],refsidetype[MAXNEWBC],indxlength;
  Real z,*newx=NULL,*newy=NULL,*newz=NULL,corder[3];
  Real meanx,meany;
  int layerbcoffset;
  int usenames;

  if(grid->rotate)
    SetElementDivisionCylinder(grid,info);
  else if(grid->dimension == 3)
    SetElementDivisionExtruded(grid,info);
  else {
    printf("CreateKnotsExtruded: unknown option!\n");
    return;
  }  

  InitializeKnots(data);

  data->dim = 3;

  origtype = 0;
  for(i=1;i<=dataxy->noelements;i++) 
    origtype = MAX( origtype, dataxy->elementtypes[i]);

  if(origtype == 202)
    elemtype = 404;
  else if(origtype == 303)  
    elemtype = 706;
  else if(origtype == 404)  
    elemtype = 808;
  else if(origtype == 408)
    elemtype = 820;
  else if(origtype == 409)
    elemtype = 827;
  else {
    printf("CreateKnotsExtruded: not implemented for elementtypes %d!\n",origtype);
    return;
  }
  printf("Maximum elementtype %d extruded to type %d.\n",origtype,elemtype);

  nonodes2d = origtype%100;
  data->maxnodes = nonodes3d = elemtype%100;
  if(nonodes3d <= 8) 
    layers = 1;
  else 
    layers = 2;

  /* Initialize the 3D mesh structure */
  data->noknots = noknots = dataxy->noknots*(layers*grid->totzelems+1);
  data->noelements = dataxy->noelements * grid->totzelems;
  data->coordsystem = dataxy->coordsystem;
  data->numbering = dataxy->numbering;
  data->noboundaries = dataxy->noboundaries;
  data->maxsize = dataxy->maxsize;
  data->minsize = dataxy->minsize;
  data->partitionexist = FALSE;
  data->periodicexist = FALSE;
  data->nodeconnectexist = FALSE;
  data->elemconnectexist = FALSE;

  usenames = dataxy->bodynamesexist || dataxy->boundarynamesexist; 
  if( usenames ) {
    if( grid->zmaterialmapexists ) {
      printf("Cannot extrude names when there is a given material mapping!\n");
      usenames = FALSE;
    }
    else {
      if(info) printf("Trying to maintain entity names in extrusion\n");
    }
  }
  

  
  maxsidetype = 0;

  AllocateKnots(data);
  indxlength = MAX(data->noknots,data->noelements);
  indx = Ivector(0,indxlength);
  for(i=0;i<=indxlength;i++)
    indx[i] = 0;

  newbounds = 0;
  if(grid->dimension == 3) 
    newbounds = grid->zcells+1;
  else if(grid->rotate) {
    if(grid->rotateblocks < 4) 
      newbounds = 4;
    if(grid->rotatecartesian)
      newbounds += grid->rotateblocks;
  }

  /* Initialize the boundaries of the 3D mesh */
  for(j=0;j<data->noboundaries+newbounds;j++) {
    if(boundxy[j].created || j>=data->noboundaries) {
      bound[j] = boundxy[j];
      bound[j].created = TRUE;

      size = bound[j].nosides = boundxy[j].nosides * grid->totzelems; 
      if(j >= data->noboundaries) size = dataxy->noelements;

      bound[j].coordsystem = COORD_CART3;
      bound[j].side = Ivector(1,size);
      bound[j].side2 = Ivector(1,size);
      bound[j].material = Ivector(1,size);    
      bound[j].parent = Ivector(1,size);
      bound[j].parent2 = Ivector(1,size);
      bound[j].types = Ivector(1,size);
      bound[j].normal = Ivector(1,size);

      for(i=1;i<=size;i++) {
	bound[j].types[i] = 0;
	bound[j].side[i] = 0;
	bound[j].side2[i] = 0;
	bound[j].parent[i] = 0;
	bound[j].parent2[i] = 0;
	bound[j].material[i] = 0;
	bound[j].normal[i] = 1;
      }
    }
  }
  if(info) printf("Allocated for %d new BC lists\n",j);
  
  knot0 = 0;
  knot1 = layers*dataxy->noknots;
  if(layers == 2) 
    knot2 = dataxy->noknots; 
  else 
    knot2 = 0;
  elem0 = 0;
  level = 0;
  material_too_large = 0;

  /* Set the element topology of the extruded mesh */
  for(cellk=1;cellk <= grid->zcells ;cellk++) {  
    
    kmax = grid->zelems[cellk];

    for(k=1;k<=kmax; k++) {
      
      level++;
      
      for(element=1;element <= dataxy->noelements;element++)  {

	origtype = dataxy->elementtypes[element];
	nonodes2d = origtype % 100;

	if(origtype == 202)
	  elemtype = 404;
	else if(origtype == 303)  
	  elemtype = 706;
	else if(origtype == 404)  
	  elemtype = 808;
	else if(origtype == 408)
	  elemtype = 820;
	else if(origtype == 409)
	  elemtype = 827;

	if( grid->zmaterialmapexists ) {
	  material = dataxy->material[element];
	  if(material > grid->maxmaterial ) {
	    material_too_large += 1;
	    continue;
	  }	  
	  material = grid->zmaterialmap[cellk][material]; 
	  if(material <= 0 ) continue;
	}
	else {
	  if(dataxy->material[element] < grid->zfirstmaterial[cellk]) continue;
	  if(dataxy->material[element] > grid->zlastmaterial[cellk]) continue;

	  if(grid->zmaterial[cellk]) 
	    material = grid->zmaterial[cellk];
	  else 
	    material = dataxy->material[element];
	}

	if(grid->rotate) {
	  meanx = 0.0;
	  for(i=0;i<nonodes2d;i++)
	    meanx += dataxy->x[dataxy->topology[element][i]];
	  meanx = fabs(meanx/nonodes2d);
	}
	if(grid->rotate && meanx < 0.0) continue;
	if(grid->rotate && cellk%2==0 && meanx < grid->rotateradius1) continue;
	
	elem0++;
	/* Vector telling the new element order. */
	indx[(level-1)*dataxy->noelements+element] = elem0;	
	data->elementtypes[elem0] = elemtype;     
	data->material[elem0] = material;

	if(elemtype == 706) {
	  for(i=0;i<3;i++) {
	    data->topology[elem0][i]   = dataxy->topology[element][i]+knot0;
	    data->topology[elem0][i+3] = dataxy->topology[element][i]+knot1;
	  }
	}
	else if(elemtype == 808) {
	  for(i=0;i<4;i++) {
	    data->topology[elem0][i]   = dataxy->topology[element][i]+knot0;
	    data->topology[elem0][i+4] = dataxy->topology[element][i]+knot1;
	  }
	}
	if(elemtype == 820 || elemtype == 827) {
	  for(i=0;i<4;i++) {
	    data->topology[elem0][i]   = dataxy->topology[element][i]+knot0;
	    data->topology[elem0][i+4] = dataxy->topology[element][i]+knot1;
	    data->topology[elem0][i+8] = dataxy->topology[element][i+4]+knot0;
	    data->topology[elem0][i+12] = dataxy->topology[element][i]+knot2;
	    data->topology[elem0][i+16] = dataxy->topology[element][i+4]+knot1;
	  }
	}
	if(elemtype == 827) {
	  for(i=0;i<4;i++) 
	      data->topology[elem0][20+i] = dataxy->topology[element][4+i]+knot2;
	  data->topology[elem0][24] = dataxy->topology[element][8]+knot0;
	  data->topology[elem0][25] = dataxy->topology[element][8]+knot1;
	  data->topology[elem0][26] = dataxy->topology[element][8]+knot2;
	}
	else if(elemtype == 404) {
	  data->topology[elem0][0] = dataxy->topology[element][0]+knot0;
	  data->topology[elem0][1] = dataxy->topology[element][1]+knot0;
	  data->topology[elem0][2] = dataxy->topology[element][1]+knot1;
	  data->topology[elem0][3] = dataxy->topology[element][0]+knot1;
	}	
      }
      knot0 += layers*dataxy->noknots;
      knot1 += layers*dataxy->noknots;
      knot2 += layers*dataxy->noknots;
    }
  }
  data->noelements = elem0;
  printf("Extruded mesh has %d elements in %d levels.\n",elem0,level);
  printf("Simple extrusion would have %d elements\n",level*dataxy->noelements);

  if( material_too_large > 0 ) {
    printf("Material index exceeded %d the size of material permutation table (%d)!\n",
	   material_too_large,grid->maxmaterial);
    printf("Give the max material with > Extruded Max Material < , if needed\n");
  }

  
  if(elem0 == 0) bigerror("No use to continue with zero elements!");

  /* Set the nodal coordinates of the extruded mesh. */
  knot0 = 0;
  for(cellk=1;cellk <= grid->zcells ;cellk++) {  

    if(cellk == 1) k=0;
    else k=1;
    for(;k<=grid->zelems[cellk]; k++) {
      
      if(grid->zlinear[cellk]) { 
	z = grid->z[cellk-1] + k*grid->dz[cellk];
      }
      else if(grid->zexpand[cellk] > 0.0) {
	z = grid->z[cellk-1] + grid->dz[cellk] * 
	  (1.- pow(grid->zratios[cellk],(Real)(k))) / (1.-grid->zratios[cellk]);
      }
      else if(grid->zelems[cellk] <= 2) {
	z = grid->z[cellk-1] + k*grid->dz[cellk];
      }
      else {
	if(k<=grid->zelems[cellk]/2) {
	  z = grid->z[cellk-1] + grid->dz[cellk] *
	    (1.- pow(grid->zratios[cellk],(Real)(k))) / (1.-grid->zratios[cellk]);
	}
	else {
	  z = grid->z[cellk] - grid->dz[cellk] * 
	    (1.- pow(grid->zratios[cellk],(Real)(grid->zelems[cellk]-k))) / (1.-grid->zratios[cellk]);
	}
      }

      for(i=1;i <= dataxy->noknots;i++)  {
	data->x[i+knot0] = dataxy->x[i];
	data->y[i+knot0] = dataxy->y[i];
	data->z[i+knot0] = z;
      }
      knot0 += layers * dataxy->noknots;
    }
  }

  
  /* Set the coordinates for the middle nodes in case 
     of quadratic elements. */  
  if(elemtype == 820 || elemtype == 827) {
    for(element=1;element <= data->noelements;element++)  {
      topo = data->topology[element];
      for(i=0;i<4;i++) {
	data->x[topo[i+12]] = 0.5*(data->x[topo[i]]+data->x[topo[i+4]]);
	data->y[topo[i+12]] = 0.5*(data->y[topo[i]]+data->y[topo[i+4]]);
	data->z[topo[i+12]] = 0.5*(data->z[topo[i]]+data->z[topo[i+4]]);
      }
      if(elemtype == 827) {
	for(i=0;i<4;i++) {
	  data->x[topo[i+20]] = 0.5*(data->x[topo[12+i]]+data->x[topo[12+(i+1)%4]]);
	  data->y[topo[i+20]] = 0.5*(data->y[topo[12+i]]+data->y[topo[12+(i+1)%4]]);
	  data->z[topo[i+20]] = 0.5*(data->z[topo[12+i]]+data->z[topo[12+(i+1)%4]]);
	}	
	data->x[topo[26]] = 0.5*(data->x[topo[0]]+data->x[topo[6]]);
	data->y[topo[26]] = 0.5*(data->y[topo[0]]+data->y[topo[6]]);
	data->z[topo[26]] = 0.5*(data->z[topo[0]]+data->z[topo[6]]);
      }
    }
  }

  /* Perform cylindrical coordinate transformation */
  if(grid->rotate) 
    CylindricalCoordinateTransformation(data,grid->rotateradius1,
					grid->rotateradius2,grid->rotatecartesian);

  cummaxsidetype = 0;
  sidetype = 0;

  /* Extrude the 2D boundary conditions. Initially BCs typically have parents with 
     different material. If due to selective extrusion they become the same then
     the extruded BC does not have that component. */
  for(j=0;j<data->noboundaries;j++) {
    if(!bound[j].created) continue;

    maxsidetype = 0;
    minsidetype = INT_MAX; 
    side  = 0;
    level = 0;

    for(cellk=1;cellk <= grid->zcells ;cellk++) {  
      for(k=1;k<=grid->zelems[cellk]; k++) {
	level++;
	
	for(i=1;i<=boundxy[j].nosides;i++){

	  /* Find the parent element indexes and the corresponding material indexes */
	  ind1 = (level-1)*dataxy->noelements + boundxy[j].parent[i];
	  parent = indx[ind1];

	  if(parent) material = data->material[parent];
	  else material  = 0;

	  if(boundxy[j].parent2[i]) {
	    ind2 = (level-1)*dataxy->noelements + boundxy[j].parent2[i];
	    parent2 = indx[ind2];
	  }
	  else 
	    parent2 = 0;

	  if(parent2) material2 = data->material[parent2];
	  else material2 = 0;

	  if((parent || parent2) && (material != material2)) {
	    side++;

	    if(!parent & !parent2) printf("no parent = %d %d %d %d %d\n",parent,parent2,ind1,ind2,level);

	    sidetype = boundxy[j].types[i];
	    bound[j].types[side] = sidetype;

	    maxsidetype = MAX( maxsidetype, sidetype );
	    minsidetype = MIN( minsidetype, sidetype );	      

	    if(parent) {
	      bound[j].parent[side] = parent;
	      bound[j].parent2[side] = parent2;
	      bound[j].side[side] = boundxy[j].side[i];
	      bound[j].side2[side] = boundxy[j].side2[i];
	      bound[j].material[side] = material;
	    }
	    else {
	      bound[j].parent[side] = parent2;
	      bound[j].parent2[side] = parent;
	      bound[j].side[side] = boundxy[j].side2[i];
	      bound[j].side2[side] = boundxy[j].side[i];
	      bound[j].material[side] = material2;	      
	    }

	    /* The sides have different convention for 1D initial elements */
	    if(elemtype == 404) {
	      if(bound[j].side[side] == 0) bound[j].side[side] = 3; 
	      if(bound[j].side2[side] == 0) bound[j].side2[side] = 3; 
	    }
	  }
	}
      }
    }
    bound[j].nosides = side;
    cummaxsidetype = MAX( maxsidetype, cummaxsidetype );
    
    if(info) {
      if(side) 
	printf("Extruded BCs list %d of types [%d,%d] has %d elements.\n",
	       j,minsidetype,maxsidetype,side);
      else
	printf("Extruded BCs list %d has no elements!\n",j);
    }

  }

  bcset = dataxy->noboundaries-1;


  if( usenames ) {
    for(i=1;i< MAXBODIES;i++) {
      if(dataxy->bodyname[i]) {
	if(!data->bodyname[i]) data->bodyname[i] = Cvector(0,MAXNAMESIZE);
	strcpy(data->bodyname[i],dataxy->bodyname[i]);
      }
    }
    for(i=1;i< MAXBCS;i++) { 
      if(dataxy->boundaryname[i]) {
	if(!data->boundaryname[i]) data->boundaryname[i] = Cvector(0,MAXNAMESIZE);
	strcpy(data->boundaryname[i],dataxy->boundaryname[i]);
      }
    }
    data->bodynamesexist = TRUE;
    data->boundarynamesexist = TRUE;
  }

  
  /* Find the BCs that are created for constant z-levels. 
     Here number all parent combinations so that each pair gets 
     a new BC index. They are numbered by their order of appearance. */
  layerbcoffset = grid->layerbcoffset;
  
  if(grid->layeredbc) {

    if( !layerbcoffset ) sidetype = maxsidetype;

    /* Find the BCs between layers. */
    if(grid->dimension == 3 || grid->rotatecartesian) {
      side = 0;
      level = 0;
      bclevel = 0;


      /* Go through extruded cells */
      for(cellk=1;cellk <= grid->zcells ;cellk++) {  
	int swap,redo;
	redo = FALSE;
	
      redolayer:
	maxsidetype = 0;
	minsidetype = INT_MAX; 

	/* Go through element layers within cells */
	for(k=1;k<=grid->zelems[cellk]; k++) {
	  level++;
	  if(!(k == 1) && !(cellk == grid->zcells && k==grid->zelems[cellk])) continue;
	  
	  /* Last cell in case of last just one element layer gives rise to two BCs */
	  if(cellk == grid->zcells && k == grid->zelems[cellk]) {
	    if(grid->zelems[cellk] == 1) 
	      redo = TRUE;
	    else {
	      level++;
	    }
	  }

	  if(grid->rotatecartesian && cellk % 2 == 1) continue; 
	  if(grid->rotatecartesian && k != 1) continue; 

	  /* If layred bc offset is defined then the BCs are numbered deterministically 
	     otherwise there is a complicated method of defining the BC index so that 
	     indexes would be used in order. */
	  if(!layerbcoffset) {
	    for(i=0;i<MAXNEWBC;i++) {
	      refmaterial1[i] = 0;
	      refmaterial2[i] = 0;
	      refsidetype[i] = 0;
	    }
	  }
	  side = 0;
	  bcset++;
	  bclevel++;
	  maxsidetype = 0;
	  minsidetype = INT_MAX;
	  
	  for(i=1;i<=dataxy->noelements;i++){
	    origtype = dataxy->elementtypes[i];
	    nonodes2d = origtype % 100;

	    if(origtype == 202)
	      elemtype = 404;
	    else if(origtype == 303)  
	      elemtype = 706;
	    else if(origtype == 404)  
	      elemtype = 808;
	    else if(origtype == 408)
	      elemtype = 820;
	    else if(origtype == 409)
	      elemtype = 827;
	    
	    /* Check the parent elements of the layers. Only create a BC if the parents are 
	       different. */
	    ind1 = (level-2)*dataxy->noelements + i;
	    if(ind1 < 1) 
	      parent = 0;
	    else
	      parent = indx[ind1];
	    
	    ind2 = (level-1)*dataxy->noelements + i;
	    if(ind2 > indxlength) 
	      parent2 = 0;
	    else
	      parent2 = indx[ind2];
	    
	    /* If only 2nd parent is given swap the order */
	    if(parent == 0 && parent2 != 0) {
	      parent = parent2;
	      parent2 = 0;
	      swap = 1;
	    } 
	    else {
	      swap = 0;
	    }	    

	    if(!parent) continue;
	    
	    /* Get the materials related to the parents */
	    material = data->material[parent];
	    if(parent2) 
	      material2 = data->material[parent2];
	    else 
	      material2 = 0;
	    
	    if(grid->rotatecartesian && !material2) {
	      if(origtype == 303) GetElementSide(parent,4-swap,1,data,sideind,&sideelemtype);
	      else GetElementSide(parent,5-swap,1,data,sideind,&sideelemtype);
	      meanx = meany = 0.0;
	      if(cellk%4 == 2) {
		for(l=0;l<sideelemtype%100;l++) {
		  meanx += data->y[sideind[l]];
		  meany += data->z[sideind[l]];
		}
	      }
	      else { 
		for(l=0;l<sideelemtype%100;l++) {
		  meanx += data->x[sideind[l]];
		  meany += data->z[sideind[l]];
		}
	      }
	      meanx = fabs(meanx)/(sideelemtype%100);
	      meany = fabs(meany)/(sideelemtype%100);
	      
	      if(fabs(meanx - grid->rotateradius1) > 1.0e-12) {
		material2 = material;
	      }
	      else {
		for(m=0;m<grid->xcells && grid->x[m]+1.0e-12 < meanx;m++);
		for(n=0;n<grid->ycells && grid->y[n]+1.0e-12 < meany;n++);
		material2 = grid->structure[n][m+1];
	      }	    
	    }
	    
	    /* Create bc index only if the materials are different */
	    if(material != material2) {	     	      
	      side++;

	      bound[bcset].nosides = side;
	      bound[bcset].parent[side] = parent;
	      bound[bcset].parent2[side] = parent2;
	      bound[bcset].material[side] = material;

	      if(origtype == 303) {
		bound[bcset].side[side] = 4-swap;
		bound[bcset].side2[side] = 3+swap;
	      }
	      else {
		bound[bcset].side[side] = 5-swap;
		bound[bcset].side2[side] = 4+swap;
	      }	      

	      /* Simple and deterministic, and complex and continuous numbering */
	      if(layerbcoffset) {
		sidetype = bclevel * layerbcoffset + dataxy->material[i];
		bound[bcset].types[side] = sidetype;
		maxsidetype = MAX( sidetype, maxsidetype );
		minsidetype = MIN( sidetype, minsidetype );
	      }
	      else {
		for(m=0;m<MAXNEWBC;m++) {
		  if(refmaterial1[m] == material && refmaterial2[m] == material2) {
		    break;
		  }
		  else if(refmaterial1[m] == 0 && refmaterial2[m] == 0) {
		    refmaterial1[m] = material;
		    refmaterial2[m] = material2;		  
		    sidetype++;
		    maxsidetype = MAX( sidetype, maxsidetype );
		    minsidetype = MIN( sidetype, minsidetype );
		    refsidetype[m] = sidetype;
		    break;
		  }
		  else if(m == MAXNEWBC-1) {
		    printf("Layer includes more than %d new BCs!\n",MAXNEWBC);
		  }
		}
		l = refsidetype[m];
		bound[bcset].types[side] = l;

		
		if( usenames ) {
		  if(!data->boundaryname[l]) data->boundaryname[l] = Cvector(0,MAXNAMESIZE);		  
		  if( bclevel == 1 ) 
		    sprintf(data->boundaryname[l],"%s%s",
			    dataxy->bodyname[dataxy->material[i]],"_Start");		  
		  else if( cellk == grid->zcells )
		    sprintf(data->boundaryname[l],"%s%s",
			    dataxy->bodyname[dataxy->material[i]],"_End");		  
		  else
		    sprintf(data->boundaryname[l],"%s%s%d",
			    dataxy->bodyname[dataxy->material[i]],"_Level",bclevel);		  
		}


	      }

	    }
	  }

	  if(info) {
	    if(side) 
	      printf("Layer BCs list %d of types [%d,%d] has %d elements.\n",
		     bcset,minsidetype,maxsidetype,side);    
	    else	      
	      printf("Layer BCs list %d has no elements!\n",bcset);
	  }

	  if(redo == TRUE) {
	    goto redolayer;
	  }
	}
      }
    } 
  }


  /* Create four additional boundaries that may be used to force
     symmetry constraints. These are only created if the object
     is only partially rotated. */

  bcset++;
  if(grid->rotate && grid->rotateblocks < 4) {
    int o,p;
    int blocks, maxradi,addtype;
    Real eps,fii,rad,meanrad,maxrad,xc,yc,dfii,fii0,rads[4],fiis[4];

    o = p = 0;
    eps = 1.0e-3;
    blocks = grid->rotateblocks;

    for(element=1;element<=data->noelements;element++) {

      for(side=0;side<6;side++) {
	GetElementSide(element,side,1,data,&sideind[0],&sideelemtype);	
	
	meanrad = 0.0;
	maxrad = 0.0;
	maxradi = 0;

	j = sideelemtype/100;
	if(j==1) break;
	
	for(i=0;i<j;i++) {
	  xc = data->x[sideind[i]];
	  yc = data->y[sideind[i]];

	  rad = sqrt(yc*yc+xc*xc);
	  fii = 2*atan2(yc,xc)/FM_PI;  /* Map fii to [0 4] */

	  rads[i] = rad;
	  fiis[i] = fii;

	  if(rad > maxrad) {
	    maxrad = rad;
	    maxradi = i;
	  }
	  meanrad += rad / j;
	}

       	fii0 = fiis[maxradi];
	dfii = 0.0;
	for(i=0;i<4;i++) {
	  if(rads[i] > eps * maxrad) {
	    if( fabs(fiis[i]-fii0) >  dfii) dfii = fabs(fiis[i]-fii0);
	  }
	}

	if(dfii > eps) continue;

	addtype = -1;

	/* BCs for zero angle */
	if(fabs(fii0) < eps) {
	  o++;
	  if(meanrad < grid->rotateradius2) 
	    addtype = 0;
	  else
	    addtype = 2;
	}
	/* BCs for angles 90, 180 or 270. */
	else if(fabs(fii0-blocks) < eps) {
	  p++;
	  if(meanrad < grid->rotateradius2) 
	    addtype = 1;
	  else
	    addtype = 3;
	}	  

	if( addtype >= 0) {	
	  bound[bcset+addtype].nosides++;
	  k = bound[bcset+addtype].nosides;
	  bound[bcset+addtype].side[k] = side;
	  bound[bcset+addtype].parent[k] = element;
	  bound[bcset+addtype].types[k] = sidetype+addtype+1;
	}
      }
    }
    
    for(addtype=0;addtype<4;addtype++) {
      l = bcset+addtype;
      if(bound[l].nosides == 0) {
	bound[l].created = FALSE;
      }
      else {
	bound[l].created = TRUE;
	if(info) {
	  if(bound[l].nosides) 
	    printf("Symmetry BCs list %d of type %d has %d elements.\n",
		   l,sidetype+addtype+1,bound[l].nosides);    
	  else
	    printf("Symmetry BCs list %d has no elements!\n",l);
	}	    
      }
    }
    bcset += 4;
  }
  data->noboundaries = bcset+1;


  /* Renumber the element nodes so that all integers are used.
     Allocate new space for the new nodes and their coordinates. */

  for(i=1;i<=data->noknots;i++)
    indx[i] = 0;
  
  for(element=1;element<=data->noelements;element++) {
    nonodes3d = data->elementtypes[element] % 100;
    for(i=0;i<nonodes3d;i++)
      indx[data->topology[element][i]] = 1;
  }  

  j = 0;
  for(i=1;i<=data->noknots;i++)
    if(indx[i])
      indx[i] = ++j;
  
  if(j < data->noknots) {
    printf("%d original nodes moved to %d new ones.\n",data->noknots,j); 
    newx = Rvector(1,j);
    for(i=1;i<=data->noknots;i++) 
      newx[indx[i]] = data->x[i];

    newy = data->x;
    data->x = newx;
    for(i=1;i<=data->noknots;i++) 
      newy[indx[i]] = data->y[i];

    newz = data->y;
    data->y = newy;
    for(i=1;i<=data->noknots;i++) 
      newz[indx[i]] = data->z[i];

    free_Rvector(data->z,1,data->noknots);
    data->z = newz;
    data->noknots = j;

    for(element=1;element<=data->noelements;element++) {
      nonodes3d = data->elementtypes[element] % 100;
      for(i=0;i<nonodes3d;i++)
	data->topology[element][i] = indx[data->topology[element][i]];
    }
  }

  if(grid->rotate) {
    ReorderElements(data,bound,FALSE,corder,info);    

    CylindricalCoordinateImprove(data,grid->rotateimprove,
				 grid->rotateradius1,grid->rotateradius2);

    if(0 && grid->rotatecurve)
      CylindricalCoordinateCurve(data,grid->curvezet,
				 grid->curverad,grid->curveangle);

    if(grid->rotatecartesian) 
      SeparateMainaxisBoundaries(data,bound);

    printf("Created %d elements and %d nodes by rotation of %d degrees.\n",
	   data->noelements,data->noknots,90*grid->rotateblocks);
  }
  else if(grid->dimension == 3)
    if(info) printf("Created %d elements and %d nodes by extruding the 2D geometry\n",
	   data->noelements,data->noknots);

  free_Ivector(indx,0,indxlength);


  /* Enforce constant helicity for the mesh if requested */
  if( grid->zhelicityexists ) {
    Real helicity,fii,x,y,z,minz,maxz;

    helicity = (FM_PI/180.0)*grid->zhelicity;

    minz = maxz = data->z[1];
    for(i=1;i<=data->noknots;i++) {
      minz = MIN(minz,data->z[i]);
      maxz = MAX(maxz,data->z[i]);
    }
    for(i=1;i<=data->noknots;i++) {
      x = data->x[i];
      y = data->y[i];
      z = data->z[i];
      fii = helicity*(z-minz)/(maxz-minz);

      data->x[i] = cos(fii)*x - sin(fii)*y;
      data->y[i] = sin(fii)*x + cos(fii)*y;
    }
    if(info) printf("Applied helicity of %12.5le degrees\n",grid->zhelicity);
  }

}



void ReduceElementOrder(struct FemType *data,int matmin,int matmax)
/* Reduces the element order at material interval [matmin,matmax] */
{
  int i,j,element,material,elemcode1,elemcode2,maxnode,reduced;
  int *indx=NULL;
  Real *newx=NULL,*newy=NULL,*newz=NULL;

  indx = Ivector(0,data->noknots);
  for(i=0;i<=data->noknots;i++)
    indx[i] = 0;
  reduced = 0;

  for(element=1;element<=data->noelements;element++) {
    elemcode1 = data->elementtypes[element];
    material = data->material[element];
    elemcode2 = elemcode1;
    if(material >= matmin && material <= matmax) 
      elemcode2 = 101*(elemcode1/100);
    if(elemcode2 == 505) elemcode2 = 504; /* tetrahedron */
    else if(elemcode2 == 606) elemcode2 = 605; /* pyramid */
    else if(elemcode2 == 707) elemcode2 = 706; /* prism */
#if 0
    printf("element=%d  codes=[%d,%d]\n",element,elemcode1,elemcode2);
    printf("mat=%d  interval=[%d,%d]\n",material,matmin,matmax);
#endif
    if(elemcode2 < elemcode1) 
      reduced++;
    maxnode = elemcode2%100;
    for(i=0;i<maxnode;i++)
      indx[data->topology[element][i]] = 1;
    data->elementtypes[element] = elemcode2;
  }

  printf("The element order is reduced in %d elements at interval [%d,%d]\n",
	 reduced,matmin,matmax);
  
  j = 0;
  for(i=1;i<=data->noknots;i++)
    if(indx[i])
      indx[i] = ++j;

  printf("%d original nodes moved to %d new ones.\n",data->noknots,j); 

  newx = Rvector(1,j);
  newy = Rvector(1,j);
  newz = Rvector(1,j);

  for(i=1;i<=data->noknots;i++) {
    newx[indx[i]] = data->x[i];
    newy[indx[i]] = data->y[i];
    newz[indx[i]] = data->z[i];
  }

  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);
  
  data->x = newx;
  data->y = newy;
  data->z = newz;
  data->noknots = j;

  for(element=1;element<=data->noelements;element++) {
    maxnode = data->elementtypes[element]%100;
    for(i=0;i<maxnode;i++) 
      data->topology[element][i] = indx[data->topology[element][i]];
  }
}




void MergeElements(struct FemType *data,struct BoundaryType *bound,
		   int manual,Real corder[],Real eps,int mergebounds,int info)
{
  int i,j,k,l;
  int noelements,noknots,newnoknots,nonodes;
  int *mergeindx=NULL,*doubles=NULL;
  Real *newx=NULL,*newy=NULL,*newz=NULL;
  Real cx,cy,cz,dx,dy,dz,cdist,dist;
  
  ReorderElements(data,bound,manual,corder,TRUE);

  /* The known ordering by vector corder[] is used to 
     reduce the cost of finding the merged nodes. */

  cx = corder[0];
  cy = corder[1];
  cz = corder[2];

  /* Normalizing for future use */
  cdist = sqrt(cx*cx+cy*cy+cz*cz);
  cx /= cdist;
  cy /= cdist;
  cz /= cdist;

  noelements  = data->noelements;
  noknots = data->noknots;
  newnoknots = noknots;

  mergeindx = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    mergeindx[i] = 0;

  doubles = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    doubles[i] = 0;

  if(info) printf("Merging nodes close (%.3lg) to one another.\n",eps);

  dz = 0.0;
  for(i=1;i<noknots;i++) {
    if(mergeindx[i]) continue;

    for(j=i+1; j<=noknots;j++) {
      if(mergeindx[j]) continue;

      dx = data->x[i] - data->x[j];
      dy = data->y[i] - data->y[j];
      dz = data->z[i] - data->z[j];

      if(fabs(cx*dx+cy*dy+cz*dz) > eps) break;

      dist = dx*dx + dy*dy + dz*dz;

      if(dist < eps*eps) {
	doubles[i] = doubles[j] = TRUE;
	mergeindx[j] = -i;	  
	newnoknots--;
      }
    }
  }

  if(mergebounds) MergeBoundaries(data,bound,doubles,info);


  j = 0;
  for(i=1;i<=noknots;i++) 
    if(mergeindx[i] == 0) 
      mergeindx[i] = ++j;

  for(i=1;i<=noknots;i++) {
    if(mergeindx[i] < 0) 
      mergeindx[i] = mergeindx[-mergeindx[i]];
  }

  printf("%d original nodes merged to %d new nodes.\n",
	 noknots,newnoknots);

  newx = Rvector(1,newnoknots);
  newy = Rvector(1,newnoknots);
  newz = Rvector(1,newnoknots);

  for(i=1;i<=noknots;i++) {
    newx[mergeindx[i]] = data->x[i];
    newy[mergeindx[i]] = data->y[i];
    newz[mergeindx[i]] = data->z[i];
  }

  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);
  
  data->x = newx;
  data->y = newy;
  data->z = newz;

#if 0
  if(info) printf("Merging the topologies.\n");
#endif

  l = 0;
  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j] % 100;
    for(i=0;i<nonodes;i++) {
      k = data->topology[j][i];
      data->topology[j][i] = mergeindx[k];
    }
  }

  data->noknots = newnoknots;
  free_Ivector(mergeindx,1,noknots);  

  if(info) printf("Merging of nodes is complete.\n");
}



void MergeBoundaries(struct FemType *data,struct BoundaryType *bound,int *doubles,int info)
{
  int i,i2,j,k,l,totsides,newsides,sidenodes,sideelemtype,side;
  int parent,sideind[MAXNODESD1];

  totsides = 0;
  newsides = 0;

  if(info) printf("Eliminating boundaries at joined nodes\n");

  for(j=0;j<MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    if(!bound[j].nosides) continue;
    
    i2 = 0;
    for(i=1;i<=bound[j].nosides;i++) {    
      
      parent = bound[j].parent[i];
      side = bound[j].side[i];
      
      GetElementSide(parent,side,1,data,sideind,&sideelemtype);
      sidenodes = sideelemtype % 100;

      l = 0;
      for(k=0;k<sidenodes;k++)
	if(doubles[sideind[k]]) l++;

      if(l < sidenodes) {
	i2++;

	if(i != i2) {
	  bound[j].parent[i2] = bound[j].parent[i];
	  bound[j].parent2[i2] = bound[j].parent2[i];
	  bound[j].side[i2] = bound[j].side[i];
	  bound[j].side2[i2] = bound[j].side2[i];
	  bound[j].types[i2] = bound[j].types[i];
	  bound[j].normal[i2] = bound[j].normal[i];
	}
      }

    }
    totsides += bound[j].nosides;
    newsides += i2;
    bound[j].nosides = i2;
    if(!i2) bound[j].created = FALSE;
  }

  if(info) printf("Eliminated %d boundaries from original set of %d.\n",totsides-newsides,totsides);
		  
}



void IsoparametricElements(struct FemType *data,struct BoundaryType *bound,
			   int bcstoo,int info)
{
  int i,j,k;
  int noelements,noknots;
  int element,side,sideelemtype,sidenodes,elemtype;
  int *bcindx=NULL,*topo=NULL,sideind[MAXNODESD1];
  Real *x=NULL,*y=NULL,*z=NULL;
  
  noelements  = data->noelements;
  noknots = data->noknots;
  x = data->x;
  y = data->y;
  z = data->z;

  bcindx = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    bcindx[i] = FALSE;

  for(j=0;j < MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    
    for(i=1; i <= bound[j].nosides; i++) {
      element = bound[j].parent[i];
      side = bound[j].side[i];
      
      GetElementSide(element,side,1,data,sideind,&sideelemtype);
      
      sidenodes = sideelemtype%100;
      
      for(k=0;k<sidenodes;k++) 
	bcindx[sideind[k]] = TRUE;
    }
  }
  
  for(j=1;j<=noelements;j++) {
    elemtype = data->elementtypes[j];
    topo = data->topology[j];
    
    if(elemtype == 306) {
      for(i=0;i<3;i++) {
	if(!bcindx[topo[i+3]]) {
	  x[topo[i+3]] = 0.5*(x[topo[i]]+x[topo[(i+1)%3]]);
	  y[topo[i+3]] = 0.5*(y[topo[i]]+y[topo[(i+1)%3]]);
	}
      }
      
    }
    else if(elemtype == 310) {
      for(i=0;i<3;i++) {
	if(!bcindx[topo[2*i+3]]) {
	  x[topo[2*i+3]] = (2.0*x[topo[i]]+1.0*x[topo[(i+1)%3]])/3.0;
	  x[topo[2*i+4]] = (1.0*x[topo[i]]+2.0*x[topo[(i+1)%3]])/3.0;
	  y[topo[2*i+3]] = (2.0*y[topo[i]]+1.0*y[topo[(i+1)%3]])/3.0;
	  y[topo[2*i+4]] = (1.0*y[topo[i]]+2.0*y[topo[(i+1)%3]])/3.0;
	}
      }      
      x[topo[9]] = (x[topo[0]]+x[topo[1]]+x[topo[2]])/3.0;
      y[topo[9]] = (y[topo[0]]+y[topo[1]]+y[topo[2]])/3.0;
    }
    else if(elemtype == 408 || elemtype == 409) {
      for(i=0;i<4;i++) {
	if(!bcindx[topo[i+4]]) {
	  x[topo[i+4]] = 0.5*(x[topo[i]]+x[topo[(i+1)%4]]);
	  y[topo[i+4]] = 0.5*(y[topo[i]]+y[topo[(i+1)%4]]);
	}
      }
      if(elemtype == 409) {
	x[topo[8]] = 0.25*(x[topo[0]]+x[topo[1]]+x[topo[2]]+x[topo[3]]);
	y[topo[8]] = 0.25*(y[topo[0]]+y[topo[1]]+y[topo[2]]+y[topo[3]]);
      }
    }    
    else if(elemtype == 412 || elemtype == 416) {
      for(i=0;i<4;i++) {
	if(!bcindx[topo[2*i+4]]) {
	  x[topo[2*i+4]] = (2.0*x[topo[i]]+1.0*x[topo[(i+1)%4]])/3.0;
	  x[topo[2*i+5]] = (1.0*x[topo[i]]+2.0*x[topo[(i+1)%4]])/3.0;
	  y[topo[2*i+4]] = (2.0*y[topo[i]]+1.0*y[topo[(i+1)%4]])/3.0;
	  y[topo[2*i+5]] = (1.0*y[topo[i]]+2.0*y[topo[(i+1)%4]])/3.0;
	}
      }      
      if(elemtype == 416) {
	Real xmean,ymean;
	xmean = (x[topo[0]]+x[topo[1]]+x[topo[2]]+x[topo[3]])/4.0;
	ymean = (y[topo[0]]+y[topo[1]]+y[topo[2]]+y[topo[3]])/4.0;
	for(i=0;i<4;i++) {
	  x[topo[11+i]] = (2.*xmean + 1.0*x[i]) / 3.0;
	  y[topo[11+i]] = (2.*ymean + 1.0*y[i]) / 3.0;
	}
      }
    }
    else {
      printf("IsoparametricElements: Not implemented for elementtype %d\n",elemtype);
    }
  }
  
  if(info) printf("The elements were forced to be isoparametric\n");
}



void ElementsToBoundaryConditions(struct FemType *data,
				  struct BoundaryType *bound,int retainorphans,int info)
{
  int i,j,k,l,sideelemtype,sideelemtype2,elemind,elemind2,sideelem,sameelem;
  int sideind[MAXNODESD1],sideind2[MAXNODESD1],elemsides,side,hit,same,minelemtype;
  int sidenodes,sidenodes2,maxelemtype,elemtype,elemdim,sideelements,material;
  int *moveelement=NULL,*parentorder=NULL,*possible=NULL,**invtopo=NULL;
  int noelements,maxpossible,noknots,maxelemsides,twiceelem,sideelemdim,minelemdim,maxelemdim;
  int debug,unmoved,removed,elemhits,loopdim,lowdimbulk;
  int notfound,*notfounds=NULL,movenames;


  if(info) {
    printf("Moving bulk elements to boundary elements\n");
    if(0) printf("Trying to retain orphans: %d\n",retainorphans);
  }
  
  for(j=0;j < MAXBOUNDARIES;j++) 
    bound[j].created = FALSE;
  for(j=0;j < MAXBOUNDARIES;j++) 
    bound[j].nosides = 0;

  noelements = data->noelements;
  noknots = data->noknots;
  
  movenames = (data->bodynamesexist && !data->boundarynamesexist);
  if(data->bodynamesexist && info) {
    if(movenames)
      printf("Moving boundarynames together with elements\n");
    else
      printf("Assuming that boundaries names are already Ok!\n");
  }

  maxelemtype = GetMaxElementType(data);
  if(info) printf("Leading bulk elementtype is %d\n",maxelemtype);

  minelemtype = GetMinElementType(data);
  if(info) printf("Trailing bulk elementtype is %d\n",minelemtype);

  maxelemdim = GetElementDimension(maxelemtype);
  minelemdim = GetElementDimension(minelemtype);
  if( maxelemdim - minelemdim == 0) {
    if(info) printf("No lower dimensional elements present!\n");
    return;
  }

  if(maxelemdim-minelemdim > 1 ) {
    int **tagcount;
    int mintag,maxtag,tag,overlap;

    mintag=maxtag=tag=overlap=-1;
    
    if(info) printf("Checking that different dimensions have unique boundary tags!\n");

    for(k=0;k<=2;k++) {
      for(i=1;i<=noelements;i++) {
	elemdim = GetElementDimension(data->elementtypes[i]);      
	tag = data->material[i];
	
	/* Get the tag interval for all elements */
	if(k==0) {
	  if(maxtag==-1) {
	    mintag = maxtag = tag;
	  }
	  else {
	    mintag = MIN(mintag,tag);
	    maxtag = MAX(maxtag,tag);
	  }
	}

	/* Count the number of tags for each dimensional */
	else if(k==1) {
	  tagcount[elemdim][tag] += 1;
	}

	/* Set the new tags for lower dimensions */
	else if(k==2) {
	  if(elemdim < maxelemdim ) { 
	    data->material[i] = tagcount[elemdim][tag];
	  }
	}
      }
      
      if(k==0) {
	if(info) printf("Tag interval for boundaries: [%d %d]\n",mintag,maxtag);
	tagcount = Imatrix(0,maxelemdim,mintag,maxtag);
	for(i=0;i<=maxelemdim;i++)
	  for(j=mintag;j<=maxtag;j++)
	    tagcount[i][j] = 0;
      }
      else if(k==1) {
	for(j=mintag;j<=maxtag;j++) {
	  overlap = 0;
	  for(i=0;i<maxelemdim;i++) 
	    if(tagcount[i][j]) overlap++;
	  if(overlap>1) break;
	}
	if(overlap>1) {
	  if(info) printf("We have an overlap, applying offsets!\n");
	  tag = 0;
	  for(i=maxelemdim-1;i>=0;i--)
	    for(j=mintag;j<=maxtag;j++) {
	      if(tagcount[i][j]) {
		if(tag+1 <= j) {
		  tag = j;
		}
		else  {
		  tag++;
		  if(info) printf("Replacing tag in %d-dim %d -> %d\n",i,j,tag);
		}
		tagcount[i][j] = tag;
	      }
	    }
	}	  
	else  {	  
	  if(info) printf("No overlap, no offsets needed!\n");
	  break;
	}	
      }
      else if(k==2) {
	if(info) printf("Renumbered tags for boundary elements!\n");
      }
    }
    free_Imatrix(tagcount,0,maxelemdim,mintag,maxtag);
  }
  
  moveelement = Ivector(1,noelements); 

  sideelements = 0;
  maxelemtype = 0;
  maxelemsides = 0;
  unmoved = 0;
  removed = 0;
  notfound = 0;
  lowdimbulk = 0;

  for(i=1;i<=noelements;i++) {
    moveelement[i] = FALSE;
    sideelemdim = GetElementDimension(data->elementtypes[i]);
    
    /* Lower dimensional elements are candidates to become BC elements */
    moveelement[i] = maxelemdim - sideelemdim;
    if(moveelement[i]) sideelements++;
  }
  if(info) printf("There are %d (out of %d) lower dimensional elements.\n",
		  sideelements,noelements);
  if(sideelements == 0) return;

  AllocateBoundary(bound,sideelements);

  /* Compute maximum number of hits for inverse topology */
  possible = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) possible[i] = 0; 
  for(elemind=1;elemind <= data->noelements;elemind++) { 
    /* if(moveelement[elemind]) continue; */
    elemtype = data->elementtypes[elemind];
    if(elemtype < 200 ) continue; 
    for(i=0;i<data->elementtypes[elemind]%100;i++) {
      j = data->topology[elemind][i];
      possible[j] += 1;
    }
  }

  j = 1;
  maxpossible = possible[1];
  for(i=1;i<=noknots;i++) {
    if(maxpossible < possible[i]) {
      maxpossible = possible[i]; 
      j = i;
    }  
  }
  if(info) printf("Node %d belongs to maximum of %d elements\n",j,maxpossible);

  /* Make a table showing to which elements a node belongs to 
     Include only the potential parents which are not to be moved to BCs. */
  invtopo = Imatrix(1,noknots,1,maxpossible);
  for(i=1;i<=noknots;i++)
    for(j=1;j<=maxpossible;j++) 
      invtopo[i][j] = 0;

  for(elemind=1;elemind <= data->noelements;elemind++) { 
    /* if(moveelement[elemind]) continue; */    
    elemtype = data->elementtypes[elemind];
    if(elemtype < 200 ) continue;
    for(i=0;i<elemtype%100;i++) {
      k = data->topology[elemind][i];
      for(l=1;invtopo[k][l];l++);
      invtopo[k][l] = elemind;
    }
  }

  sideelem = 0;
  sameelem = 0;
  twiceelem = 0;

  debug = FALSE;

  /* Go through boundary element candidates starting from higher dimension */
  for(loopdim=maxelemdim-1;loopdim>=0;loopdim--) {

    if(debug) printf("loopdim = %d\n",loopdim);
    
    for(elemind=1;elemind <= data->noelements;elemind++) {

      if(!moveelement[elemind]) continue;

      same = FALSE;
      sideelemtype = data->elementtypes[elemind];
      
      /* Only check the elements that have right dimension */
      sideelemdim = GetElementDimension(sideelemtype);
      if(sideelemdim != loopdim ) continue;
      
      sidenodes = sideelemtype % 100;
      for(i=0;i<sidenodes;i++) 
	sideind[i] = data->topology[elemind][i];
      elemhits = 0;    

      if(debug) printf("Finding elem: %d %d %d\n",elemind,sideelemtype,sideelemdim);

      
      for(l=1;l<=maxpossible;l++) {
	elemind2 = invtopo[sideind[0]][l];
	
	if(!elemind2) continue;

	/* The parent should be an element that will not become BC element */
	if(moveelement[elemind2]) continue;
	
	elemtype = data->elementtypes[elemind2];
	elemdim = GetElementDimension(elemtype);
	
	/* Owner element should have higher dimension */
	if(elemdim <= sideelemdim ) continue;

	hit = 0;
	for(i=0;i<sidenodes;i++)
	  for(j=0;j<elemtype%100;j++)
	    if(sideind[i] == data->topology[elemind2][j]) hit++;

	if(hit < sidenodes) continue;
	
	if(hit > sidenodes) printf("Strange: elemhits %d vs. elemnodes %d\n",hit,sidenodes);
	if(hit >= sidenodes) elemhits++;
	
	for(side=0;side<=100;side++) {	  
	  GetElementSide(elemind2,side,1,data,&sideind2[0],&sideelemtype2);
	  
	  if(sideelemtype2 == 0 ) break;

	  if(sideelemtype2 < 300 && sideelemtype > 300) break;	
	  if(sideelemtype2 < 200 && sideelemtype > 200) break;		
	  if(sideelemtype2 != sideelemtype ) continue;	  

	  sidenodes2 = sideelemtype2 % 100;	
	  if(sidenodes != sidenodes2) continue;
	  if(sidenodes2 == 1 && sidenodes > 1) break;
	  
	  hit = 0;
	  for(i=0;i<sidenodes;i++) 
	    for(j=0;j<sidenodes2;j++) 
	      if(sideind[i] == sideind2[j]) hit++;
	  
	  if(0) printf("%d hits in element %d\n",hit,sideelemtype2);
	  if(hit == sidenodes) break;
 
	  if(sideelemtype != sideelemtype2) {
	    printf("Hits in element after mismatch: %d vs. %d\n",sideelemtype,sideelemtype2);
	    continue;
	  }
	}
	if( sidenodes == 1 && !hit) {
	  printf("elemind = %d sideind %d vs. ",elemind,sideind[0]);
	  for(j=0;j<sidenodes2;j++) 
	    printf("%d ",sideind2[j]);
	  printf("\n");						 
	}
	
	if(hit < sidenodes || !sideelemtype2) {
	  if(0) printf("Preliminary hit but not really: %d %d\n",hit,sidenodes);
	  continue;
	}
	  
	if(same) {
	  sameelem += 1;
	  bound->parent2[sideelem] = elemind2;
	  bound->side2[sideelem] = side;

	  if(debug) printf("  Found 2nd: %d %d %d\n",elemind,elemind2,side);
	  goto foundtwo;
	}
	else {
	  sideelem += 1;
	  same = TRUE;
	  if(debug) printf("  Found 1st: %d %d %d %d %d\n",elemind,elemind2,side,sideelemtype,sideelemtype2);

	  bound->parent[sideelem] = elemind2;
	  bound->side[sideelem] = side;
	  bound->parent2[sideelem] = 0;
	  bound->side2[sideelem] = 0;	    
	  material = data->material[elemind];
	  bound->types[sideelem] = material;
	    	    
	  if(sidenodes == 2) {
	    if((sideind[0]-sideind[1])*(sideind2[0]-sideind2[1])<0) 	      
	      bound->normal[sideelem] = -1;
	  }		
	  if(movenames) {
	    data->boundarynamesexist = TRUE;
	    if(material < MAXBODIES && material < MAXBCS) {
	      if(!data->boundaryname[material]) {
		data->boundaryname[material] = Cvector(0,MAXNAMESIZE);
		if(data->bodyname[material]) {
		  strcpy(data->boundaryname[material],data->bodyname[material]);
		  free_Cvector(data->bodyname[material],0,MAXNAMESIZE);
		  data->bodyname[material] = NULL;
		}
		else
		  sprintf(data->boundaryname[material],"body%d",material);
	      }
	    }
	    if(!strncmp(data->boundaryname[material],"body",4)) {
	      strncpy(data->boundaryname[material],"bnry",4);
	    }
	  }

	  /* Only try to find two parents if the boundary element is one degree smaller than maximum dimension */
	  if(moveelement[elemind] > 1) goto foundtwo;
	}
      }
      
      if(!same) {			
	/* If the element is of dimension DIM-1 then create a table showing where they are */
	if(retainorphans ) {
	  /* If we have only one degree smaller unfound element then keep them as bulk elements. */
	  if( moveelement[elemind] == 1) {	
	    moveelement[elemind] = 0;
	    lowdimbulk++;
	    if(debug) printf("  Bulk: %d\n",elemind);
	  }
	  else {
	    if(!notfound) {
	      notfounds = Ivector(1,noelements); 
	      for(i=1;i<=noelements;i++)
		notfounds[i] = FALSE;
	    }
	    notfound++;
	    notfounds[elemind] = TRUE;

	    if(0) {
	      printf("element: elemind = %d type = %d nodes = %d elemhits = %d\n",
		     elemind,sideelemtype,sidenodes,elemhits);
	      printf("         inds =");
	      for(i=0;i<sidenodes;i++)
		printf(" %d ",sideind[i]);
	      printf("\n");
	    }
	    
	    if(debug) printf("  Unfound: %d\n",elemind);	      
	  }
	}
	else {
	  if(debug) printf("  Removed: %d\n",elemind);

	  moveelement[elemind] = -1;
	  removed += 1;
	}	
      }

    foundtwo:
      continue;
      
    }

    if(0) printf("Intermediate results: %d %d %d %d\n",twiceelem,sameelem,sideelem,removed);
  }

  if(twiceelem) printf("Found %d sides that were multiply given\n",twiceelem);
  if(sameelem) printf("Found %d side elements that have two parents.\n",sameelem);


  if(sideelem == sideelements) {
    printf("Found correctly %d side elements.\n",sideelem);
  }
  else {
    printf("Studied %d lower dimensional elements\n",sideelements);
    printf("Defined %d side elements\n",sideelem);
    printf("Defined %d lower dimensional bulk elements\n",lowdimbulk);

    bound->nosides = sideelem;
    
    printf("Removing %d lower dimensional elements from the element list\n",removed);
    if(notfound) {
      if(0) printf("************************** WARNING **********************\n");
      if(retainorphans) {
	printf("Adding %d elements to boundary without parent information\n",notfound);
	
	bound->elementtypes = Ivector(sideelem+1,sideelements);
	for(i=sideelem+1;i<=sideelements;i++) bound->elementtypes[i] = 0;
	
	bound->topology = Imatrix(sideelem+1,sideelements,0,MAXNODESD2-1);
	
	for(elemind=1;elemind <= data->noelements;elemind++) {
	  if(!notfounds[elemind]) continue;
	  sideelem++;

	  j = data->elementtypes[elemind];
	  bound->elementtypes[sideelem] = j;	  

	  for(i=0;i<j%100;i++)
	    bound->topology[sideelem][i] = data->topology[elemind][i];
	  
	  /* Adding some constant here could be used for debugging */
	  bound->types[sideelem] = data->material[elemind] + 1*10;
	  bound->parent[sideelem] = 0;
	}
	bound->nosides = sideelem;
      }
      else {
	printf("Removing %d lower dimensional elements without parent information\n",notfound);	
      }
    }
  }

  /* Reorder remaining bulk elements */
  parentorder = Ivector(1,noelements);
  for(i=1;i<=noelements;i++)
    parentorder[i] = 0;
    
  j = 0;
  for(i=1;i<=noelements;i++) {
    if(moveelement[i] == 0) {
      k = data->elementtypes[i];

      j++;
      parentorder[i] = j;

      if(debug) printf("Bulk is: %d %d\n",i,j);

      if( i != j ) {
	data->material[j] = data->material[i];
	data->elementtypes[j] = data->elementtypes[i];	
	for(l=0;l<k%100;l++) 
	  data->topology[j][l] = data->topology[i][l];     
      }
    }
  }
  data->noelements = j;
  if(info) printf("Parent elements were reordered up to index %d.\n",j);


  /* Reorder boundary to point at the new arrangement of master elements */
  for(i=1;i<=bound->nosides;i++) {
    if(!bound->parent[i]) continue;

    if( !parentorder[bound->parent[i]] ) {
      printf("Zero reorder: %d %d %d\n",i,bound->parent[i],bound->side[i]);
      bigerror("Sorry folks!");
    }
    
    if(bound->parent[i])  bound->parent[i] = parentorder[bound->parent[i]];
    if(bound->parent2[i])  bound->parent2[i] = parentorder[bound->parent2[i]];
    
    GetElementSide(bound->parent[i],bound->side[i],1,data,&sideind2[0],&sideelemtype2);

    if(0) GetBoundaryElement(i,&bound[j],data,&sideind2[0],&sideelemtype2);
  }

  if(info) printf("Moved %d elements (out of %d) to new positions\n",j,noelements);

  free_Ivector(parentorder,1,noelements);

  free_Ivector(moveelement,1,noelements); 
  free_Ivector(possible,1,noknots);
  free_Imatrix(invtopo,1,noknots,1,maxpossible);
  if(notfound) free_Ivector(notfounds,1,noelements);

  if(debug) printf("All done\n");
  
  return;
}


int SideAndBulkMappings(struct FemType *data,struct BoundaryType *bound,struct ElmergridType *eg,int info)
{
  int i,j,l,currenttype;
  

  if(eg->sidemappings) {
    if(info) printf("Renumbering boundary types with %d mappings\n",eg->sidemappings);

    for(l=0;l<eg->sidemappings;l++) 
      if(info) printf("Setting boundary types between %d and %d to %d\n",
		      eg->sidemap[3*l],eg->sidemap[3*l+1],eg->sidemap[3*l+2]);

    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
	
	for(i=1; i <= bound[j].nosides; i++) {
	  if(currenttype = bound[j].types[i]) {
	    for(l=0;l<eg->sidemappings;l++) {
	      if(currenttype >= eg->sidemap[3*l] && currenttype <= eg->sidemap[3*l+1]) {
		bound[j].types[i] = eg->sidemap[3*l+2];
		currenttype = -1;
	      }
	    }
	  }
	}
    }
    if(info) printf("Renumbering boundary types finished\n");
  }
  
  if(eg->bulkmappings) {
    if(info) printf("Renumbering bulk types with %d mappings\n",eg->bulkmappings);

    for(l=0;l<eg->bulkmappings;l++) 
      if(info) printf("Setting material types between %d and %d to %d\n",
		      eg->bulkmap[3*l],eg->bulkmap[3*l+1],eg->bulkmap[3*l+2]);
    for(j=1;j<=data->noelements;j++) {
      currenttype = data->material[j];
      for(l=0;l<eg->bulkmappings;l++) {
	if(currenttype >= eg->bulkmap[3*l] && currenttype <= eg->bulkmap[3*l+1]) {
	  data->material[j] = eg->bulkmap[3*l+2];
	  currenttype = -1;
	}
      }
    }
    if(info) printf("Renumbering material indexes finished\n");
  }
  return(0);
}



int SideAndBulkBoundaries(struct FemType *data,struct BoundaryType *bound,struct ElmergridType *eg,int info)
{
  int l;
  int *boundnodes,noboundnodes;
  boundnodes = Ivector(1,data->noknots);
      
  if(eg->bulkbounds) {
    for(l=0;l<eg->bulkbounds;l++) {
      FindBulkBoundary(data,eg->bulkbound[3*l],eg->bulkbound[3*l+1],
		       boundnodes,&noboundnodes,info);
      FindNewBoundaries(data,bound,boundnodes,eg->bulkbound[3*l+2],1,info);
    }
  }
  if(eg->boundbounds) {
    for(l=0;l<eg->boundbounds;l++) {	
      FindBoundaryBoundary(data,bound,eg->boundbound[3*l],eg->boundbound[3*l+1],
			   boundnodes,&noboundnodes,info);
      FindNewBoundaries(data,bound,boundnodes,eg->boundbound[3*l+2],2,info);
    }
  }
  free_Ivector(boundnodes,1,data->noknots);

  return(0);
}


void NodesToBoundaryChain(struct FemType *data,struct BoundaryType *bound,
			  int *bcinds,int *bctags,int nbc,int bccount,
			  int info)
{
  int i,j,k,l,sideelemtype,sideelemtype2,elemind,elemind2,sideelem,sameelem;
  int sideind[MAXNODESD1],sideind2[MAXNODESD1],elemsides,side,hit,same,minelemtype;
  int sidenodes,sidenodes2,elemtype,elemdim,sideelements,material;
  int *possible=NULL,**invtopo=NULL;
  int noelements,maxpossible,noknots,twiceelem,sideelemdim;
  int elemhits,bci;


  if(info) printf("Creating boundary elements from boundary nodes\n");

  for(j=0;j < MAXBOUNDARIES;j++) 
    bound[j].created = FALSE;
  for(j=0;j < MAXBOUNDARIES;j++) 
    bound[j].nosides = 0;
  
  noelements = data->noelements;
  noknots = data->noknots;
  
  sideelements = nbc - bccount;
  printf("Expected number of BC elements: %d\n",sideelements);

  AllocateBoundary(bound,sideelements);
  
  /* Calculate how may times a node appears */
  possible = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) possible[i] = 0; 
  for(elemind=1;elemind <= data->noelements;elemind++) { 
    for(i=0;i<data->elementtypes[elemind]%100;i++) {
      j = data->topology[elemind][i];
      possible[j] += 1;
    }
  }
  
  j = 1;
  maxpossible = possible[1];
  for(i=1;i<=noknots;i++) {
    if(maxpossible < possible[i]) {
      maxpossible = possible[i]; 
      j = i;
    }  
  }
  if(info) printf("Node %d belongs to maximum of %d elements\n",j,maxpossible);
  
  /* Make a table showing to which elements a node belongs to 
     Include only the potential parents which are not to be moved to BCs. */
  invtopo = Imatrix(1,noknots,1,maxpossible);

  for(i=1;i<=noknots;i++)
    for(j=1;j<=maxpossible;j++) 
      invtopo[i][j] = 0;
  
  for(elemind=1;elemind <= data->noelements;elemind++) { 
    elemtype = data->elementtypes[elemind];
    for(i=0;i<elemtype%100;i++) {
      k = data->topology[elemind][i];
      for(l=1;invtopo[k][l];l++); /* Yes, this is really ok. We look for unset entry. */
      invtopo[k][l] = elemind;
    }
  }
    
  sideelem = 0;
  sameelem = 0;
  twiceelem = 0;

  /* These are here by construction because we are looking for a chain of nodes
     and trying to create 202 elements of them! */
  sidenodes = 2;
  sideelemtype = 202;
  
  for(bci=1;bci<nbc;bci++) {

    same = FALSE;

    if( bctags[bci] != bctags[bci+1] ) continue; 

    sideind[0] = bcinds[bci];
    sideind[1] = bcinds[bci+1];
    material = bctags[bci];
    
    elemhits = 0;    
    
    /* Go through potential parents elements using the inverse topology */
    for(l=1;l<=maxpossible;l++) {
      elemind2 = invtopo[sideind[0]][l];
    
      if(!elemind2) continue;      

      elemtype = data->elementtypes[elemind2];
      hit = 0;
      for(i=0;i<sidenodes;i++)
	for(j=0;j<elemtype%100;j++)
	  if(sideind[i] == data->topology[elemind2][j]) hit++;

      /* We must have all hits to have a chance of finding bc */
      if(hit < sidenodes) continue;

      elemhits++;

      /* Now find on which side the bc is */
      for(side=0;side<3;side++) {
	GetElementSide(elemind2,side,1,data,&sideind2[0],&sideelemtype2);
	if( sideelemtype2 != sideelemtype ) printf("This should not happen!\n");

	hit = 0;
	for(i=0;i<sidenodes;i++) 
	  for(j=0;j<sidenodes;j++) 
	    if(sideind[i] == sideind2[j]) hit++;
	
	if(hit < sidenodes) continue;

	if(same) {
	  /* Ok, we found the other parent for this already */
	  sameelem += 1;
	  bound->parent2[sideelem] = elemind2;
	  bound->side2[sideelem] = side;		  
	  goto foundtwo;
	}
	else {
	  /* We haven't found parents for this bc elements yet */
	  sideelem += 1;
	  same = TRUE;
	  bound->parent[sideelem] = elemind2;
	  bound->side[sideelem] = side;
	  bound->parent2[sideelem] = 0;
	  bound->side2[sideelem] = 0;
	  bound->types[sideelem] = material;
	  if(sidenodes == 2) {
	    if((sideind[0]-sideind[1])*(sideind2[0]-sideind2[1])<0) 	      
	      bound->normal[sideelem] = -1;
	  }		
	}
      }
    }
  foundtwo:
    continue;
  }
  
  if(twiceelem) printf("Found %d sides that were multiply given\n",twiceelem);
  if(sameelem) printf("Found %d side elements that have two parents.\n",sameelem);
  
  
  if(sideelem == sideelements) {
    printf("Found correctly %d side elements.\n",sideelem);
  }
  else {
    printf("Found %d side elements, could have found %d\n",sideelem,sideelements);
  }
        
  bound->nosides = sideelem;

  free_Ivector(possible,1,noknots);
  free_Imatrix(invtopo,1,noknots,1,maxpossible);

  return;
}




int FindPeriodicNodes(struct FemType *data,int periodicdim[],int info)
{
  int i,j,i2,j2,dim;
  int noknots,hit,tothits;
  int *topbot=NULL,*indxper=NULL;
  int botn,topn,*revindtop=NULL,*revindbot=NULL;
  Real eps,dist,dx,dy,dz,coordmax,coordmin;
  Real *coord=NULL,*toparr=NULL,*botarr=NULL,epsmin;


  if(data->dim < 3) periodicdim[2] = 0;
  if(!periodicdim[0] && !periodicdim[1] && !periodicdim[2]) return(1);

  if(data->periodicexist) {
    printf("FindPeriodicNodes: Subroutine is called for second time?\n");
    return(2);
  }

  noknots = data->noknots;
  tothits = 0;

  data->periodicexist = TRUE;
  indxper = Ivector(1,noknots);
  data->periodic = indxper;
  topbot = Ivector(1,noknots);
  
  
  for(i=1;i<=noknots;i++) 
    indxper[i] = i;
  
  for(dim=1;dim<=3;dim++) {
    if(!periodicdim[dim-1]) continue;

    if(info) printf("Finding periodic nodes in direction %d\n",dim);
    
    if(dim==1) coord = data->x;
    else if(dim==2) coord = data->y;
    else coord = data->z;
    
    coordmax = coordmin = coord[1];
    
    for(i=1;i<=data->noknots;i++) {
      if(coordmax < coord[i]) coordmax = coord[i];
      if(coordmin > coord[i]) coordmin = coord[i];
    }

    if(info) printf("Coordinate in dimension %d is at the interval [%.3lg, %.3lg]\n",
		    dim,coordmin,coordmax);

    if(coordmax-coordmin < 1.0e-10) continue;
    eps = 1.0e-5 * (coordmax-coordmin);

    topn = botn = 0;
    for(i=1;i<=data->noknots;i++) {
      if(fabs(coord[i]-coordmax) < eps) {
	topn++;
	topbot[i] = topn;
      }
      else if(fabs(coord[i] - coordmin) < eps) {
	botn++;
	topbot[i] = -botn;
      }
      else {
	topbot[i] = 0;
      }
    }
    
    if(topn != botn) {
      printf("There should be equal number of top and bottom nodes (%d vs. %d)!\n",topn,botn);
      return(3);
    }
    else {
      if(info) printf("Looking for %d periodic nodes\n",topn);
    }

    toparr = Rvector(1,topn);
    botarr = Rvector(1,botn);
    revindtop = Ivector(1,topn);
    revindbot = Ivector(1,botn);
    
    topn = botn = 0;
    for(i=1;i<=noknots;i++) {
      j = topbot[i];
      if(j > 0) {
	topn++;
	revindtop[topn] = i;  
      }
      else if(j < 0) {
	j = abs(j);
	botn++;
	revindbot[botn] = i;  
      }
    }
    
    if(data->dim == 2) {
      for(i=1;i<=botn;i++) {
	j = revindbot[i];
	hit = FALSE;
	for(i2=1;i2<=topn;i2++) {
	  j2 = revindtop[i2];
	  if(dim == 1) 
	    dist = fabs(data->y[j] - data->y[j2]);
	  else 
	    dist = fabs(data->x[j] - data->x[j2]);
	  if(dist < eps) {
	    hit = TRUE;
	    goto hit2d;
	  }
	}

      hit2d:
	if(hit) {
	  tothits++;
	  if(indxper[j] == j) indxper[j2] = j;
	  else if(indxper[indxper[j]]==indxper[j]) {
	    indxper[j2] = indxper[j];
	  }
	  else {
	    printf("unknown 2d case!\n");
	  }
	}
	else {
	  printf("Couldn't find a periodic counterpart for node %d at [%.3lg %.3lg]]\n",
		 j,data->x[j],data->y[j]);
	}
      }
    }
    else if(data->dim == 3) {
      dx = dy = dz = 0.0;
      for(i=1;i<=botn;i++) {
	j = revindbot[i];
	hit = FALSE;
	epsmin = coordmax - coordmin;
	
	for(i2=1;i2<=topn;i2++) {
	  j2 = revindtop[i2];
	  if(dim == 1) {
	    dy = data->y[j] - data->y[j2];
	    dz = data->z[j] - data->z[j2];
	  }
	  else if(dim == 2) {
	    dx = data->x[j] - data->x[j2];
	    dz = data->z[j] - data->z[j2];
	  }
	  else {
	    dx = data->x[j] - data->x[j2];
	    dy = data->y[j] - data->y[j2];
	  }
	  if(dx*dx+dy*dy+dz*dz < eps*eps) {
	    hit = TRUE;
	    goto hit3d;
	  }
	}

      hit3d:
	if(hit) {
	  tothits++;	  
          indxper[j2] = indxper[j];
	}
	else {
	  printf("The periodic counterpart for node %d was not found!\n",j);
	}
      }
    }

    free_Rvector(toparr,1,topn);
    free_Rvector(botarr,1,botn);
    free_Ivector(revindtop,1,topn);
    free_Ivector(revindbot,1,botn);
  }

  if(info) printf("Found all in all %d periodic nodes.\n",tothits);

  free_Ivector(topbot,1,noknots);

  return(0);
}




int FindPeriodicParents(struct FemType *data,struct BoundaryType *bound,int info)
{
  int i,j,k,k2,l,l2,totsides,newsides,sidenodes,sideelemtype,side;
  int noknots,maxhits,nodes,hits,hits2,targets,mappings,targetnode;
  int parent,parent2,sideind[MAXNODESD1],sideind2[MAXNODESD1];
  int **periodicparents=NULL, *periodichits=NULL,*periodictarget=NULL,*indexper=NULL;
  
  totsides = 0;
  newsides = 0;
  targets = 0;
  parent2 = 0;
  
  if(info) printf("Finding secondary periodic parents for boundary elements\n");
  
  if(!data->periodicexist) {
    printf("FindPeriodicParents: Periodic nodes are not defined\n");
    return(2);
  }
  
  indexper = data->periodic;
  
  /* Set pointers that point to the periodic nodes */
  noknots = data->noknots;
  periodictarget = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    periodictarget[i] = 0;

  mappings = 0;
  for(i=1;i<=noknots;i++) {
    j = indexper[i];
    if( j != i) {
      mappings++;
      periodictarget[j] = i;      
    }
  } 

  if(0) for(i=1;i<=noknots;i++)
    printf("indexes(%d) : %d %d\n",i,indexper[i],periodictarget[i]);


  if(info) printf("Number of potential periodic mappings is %d\n",mappings);
  for(i=1;i<=noknots;i++) 
    if(periodictarget[i]) targets++;
  if(info) printf("Number of potential periodic targets is %d\n",targets);
  

  /* Vector telling how many elements are associated with the periodic nodes */
  maxhits = 0;
  periodichits = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    periodichits[i] = 0;

  /* Create the matrix telling which elements are associated with the periodic nodes */
 setparents:
  for(j=1;j <= data->noelements;j++) {
    nodes = data->elementtypes[j] % 100;    
    for(i=0;i<nodes;i++) {
      k = data->topology[j][i];
      if( k != indexper[k] ) {
	periodichits[k] += 1;
	if( maxhits > 0 ) {
	  periodicparents[k][periodichits[k]] = j;
	}
      }
    }
  }

  if( maxhits == 0 )   {
    for(i=1;i<=noknots;i++) 
      maxhits = MAX( maxhits, periodichits[i] );

    printf("Maximum number of elements associated with periodic nodes is %d\n",maxhits);
    periodicparents = Imatrix(1,noknots,1,maxhits);
    for(i=1;i<=noknots;i++) {
      periodichits[i] = 0;
      for(j=1;j<=maxhits;j++) 
        periodicparents[i][j] = 0;
    }
    goto setparents;
  }      

  for(j=0;j<MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    if(!bound[j].nosides) continue;
    
    for(i=1;i<=bound[j].nosides;i++) {    
      
      /* If secondary parent already set skip */
      if(bound[j].parent2[i]) continue;

      parent = bound[j].parent[i];
      if(0) printf("1st parent %d\n",parent);

      side = bound[j].side[i];
      
      GetElementSide(parent,side,1,data,sideind,&sideelemtype);
      sidenodes = sideelemtype % 100;

      /* Some node must be periodic target and others either target or mapped nodes */
      hits = hits2 = 0;      
      for(k=0;k<sidenodes;k++) {
        l = sideind[k];
        if( periodictarget[l] ) 
          hits++;
        else if(indexper[l] != l) 
          hits2++;
      }
      if(!hits || hits + hits2 < sidenodes) continue;

      if(0) printf("Trying to find other parent for boundary %d and parent %d\n",
		   bound[j].types[i],parent);
      if(0) printf("hits = %d %d %d\n",hits,hits2,sidenodes);

      totsides++;

      /* The parent is the one element that has exactly the same set of periodic nodes */
      for(l=0;l<sidenodes;l++) {
        targetnode = periodictarget[sideind[l]];
	if(!targetnode) continue;
    
	for(l2=1;l2<=periodichits[targetnode];l2++) {
	  int side,elemtype,elemsides,sideelemtype2;
	  
	  parent2 = periodicparents[targetnode][l2];
	  if(parent == parent2) continue;
	  
	  elemtype = data->elementtypes[parent2];
	  elemsides = GetElementFaces(elemtype);
	  
	  for(side=0;side<elemsides;side++) { 
	    GetElementSide(parent2,side,1,data,sideind2,&sideelemtype2);
	    if( sideelemtype != sideelemtype2 ) continue;
	    	    
	    hits = 0;
	    for(k=0;k<sidenodes;k++) {
	      for(k2=0;k2<sidenodes;k2++) {
		if( indexper[sideind[k]] == indexper[sideind2[k2]]) {
		  hits++;
		  break;
		}
	      }
	    }
	    if(hits == sidenodes) goto found;
	  }
	}
      }

    found:     
      if(hits == sidenodes) {
	newsides++;
	if(0) printf("Parents joined by boundary element: %d %d\n",parent,parent2);
	bound[j].parent2[i] = -parent2;
      }
      else {
	printf("Could not find a periodic counterpart: %d/%d/%d\n",j,i,parent);
	printf("ind = %d ",sideind[0]);
	for(k=1;k<sidenodes;k++)
	  printf("%d ",sideind[k]);
	printf("\n");
      }
    }
  }

  free_Ivector(periodictarget,1,noknots);
  free_Ivector(periodichits,1,noknots);
  free_Imatrix(periodicparents,1,noknots,1,maxhits);

  if(info) printf("Found %d secondary parents for %d potential sides.\n",newsides,totsides); 
  return(0);
}




int CreateBoundaryLayer(struct FemType *data,struct BoundaryType *bound,
			int nolayers, int *layerbounds, int *layernumber,
			Real *layerratios, Real *layerthickness, int *layerparents,
			int maxfilters, Real layereps, int info)
/* Create Boundary layers that may be used to solve accurately fluid
   flow problems and similar equations. */
{
  int i,j,k,l,m,n,i2,i3,nonodes,maxbc,newbc;
  int noknots,noelements,elemindx,nodeindx,elemtype;
  int oldnoknots,oldnoelements,maxelemtype,oldmaxnodes;
  int nonewnodes,nonewelements,dolayer,dim,order,midpoints;
  int checkmaterials,parent,parent2,use2,second;
  Real dx,dy,ds,ratio,q,p,rectfactor;
  Real *newx=NULL,*newy=NULL,*newz=NULL,*oldx=NULL,*oldy=NULL,*elemwidth=NULL;
  Real e1x,e1y,e2x,e2y;
  int sideelemtype,ind[MAXNODESD2],sidebc[MAXNODESD1];
  int *layernode=NULL,*newelementtypes=NULL,**newtopo=NULL,**oldtopo=NULL;
  int *topomap=NULL,*newmaterial=NULL,*herit=NULL,*inside=NULL,*nonlin=NULL;
  int endbcs, *endparents=NULL, *endtypes=NULL, *endnodes=NULL, *endnodes2=NULL, *endneighbours=NULL;

  if(0) printf("maxfilters=%d layereps=%.3e\n",maxfilters,layereps);

  if(!maxfilters) maxfilters = 1000;
  if(layereps < 1.0e-20) layereps = 1.0e-3;
  rectfactor = 1.0e2;
  midpoints = FALSE;
  order = 1;
  dim = data->dim;
  
  maxelemtype = GetMaxElementType(data);
  if(maxelemtype > 409) {
    printf("Subroutine implemented only up to 2nd degree in 2D!\n");
    bigerror("Cannot continue");
  }

  if(info) printf("Largest elementtype is %d\n",maxelemtype);
  
  second = FALSE;
  checkmaterials = FALSE;
  for(k=0;k<nolayers;k++) 
    if(layerparents[k]) checkmaterials = TRUE;


omstart:

  oldnoelements = noelements = data->noelements;
  oldnoknots = noknots = data->noknots;
  oldmaxnodes = data->maxnodes;

  layernode = Ivector(1,oldnoknots);
  for(i=1;i<=oldnoknots;i++) layernode[i] = 0;


  /* Go through all the boundaries with boundary layer definitions and compute 
     the number of new nodes and new elements. */
  nonewnodes = 0;
  nonewelements = 0;
  maxbc = 0;

  /* Go through the layers and check which ones are active */
  for(j=0;j<MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;

    for(i=1;i<=bound[j].nosides;i++) {
      dolayer = FALSE;
      parent = bound[j].parent[i];
      use2 = FALSE;
      if(bound[j].types[i] > maxbc) maxbc = bound[j].types[i];

      for(k=0;k<nolayers;k++) {

	if(bound[j].types[i] == layerbounds[k]) {
	  if(checkmaterials) {
	    if(layerparents[k] < 0) continue;

	    if(data->material[parent] == layerparents[k])
	      dolayer = k + 1;
	    else if(parent = bound[j].parent2[i]) {
	      if(data->material[parent] == layerparents[k]) {
		use2 = TRUE;
		dolayer = k + 1;
	      }
	    }
	  }
	  else {
	    dolayer = k + 1;
	  }
	}
      }

      if(!dolayer) continue;


      /* We have found an boundary element to extrude */
      if(use2) 
	GetElementSide(bound[j].parent2[i],bound[j].side2[i],bound[j].normal[i],
		       data,ind,&sideelemtype);
      else 
	GetElementSide(parent,bound[j].side[i],bound[j].normal[i],
		       data,ind,&sideelemtype);

      nonewelements += layernumber[dolayer-1];
      
      midpoints = FALSE;
      if(sideelemtype == 202) {
	order = 1;
      }
      else if(sideelemtype == 203) {
	order = 2;
	if(maxelemtype > 408) midpoints = TRUE;
      }
      
      for(l=0;l<sideelemtype%100;l++) {

	/* No layer has yet been created for this node */
	if(!layernode[ind[l]]) {

	  layernode[ind[l]] = -(noknots + nonewnodes);	  

	  if(l < sideelemtype/100 || midpoints) 
	    nonewnodes += order * layernumber[dolayer-1];
	  else 
	    nonewnodes += layernumber[dolayer-1];	      
	}
	else {
	  layernode[ind[l]] = abs(layernode[ind[l]]);
	}
      }
    }
  }
  
  if(!nonewelements) {
    if(info) printf("Found no active boundary layers!\n");
    return(0);
  }

  /* For higher order elements remove the middlenodes from the list of cornernodes */
  if(maxelemtype%100 > 4) {
    for(j=1;j<=noelements;j++) {      
      elemtype = data->elementtypes[j];	  
      for(i=elemtype/100;i<elemtype%100;i++) {
	k = data->topology[j][i];
	layernode[k] = abs(layernode[k]);
      }
    }
  }


  /* Negative indexed means that the node is an end node of the newly created boundary */
  endbcs = 0;
  for(i=1;i<=noknots;i++) 
    if(layernode[i] < 0) endbcs++;

  if(endbcs) {
    endparents = Ivector(1,endbcs);
    endtypes = Ivector(1,endbcs);
    endnodes = Ivector(1,endbcs);
    endnodes2 = Ivector(1,endbcs);

    endneighbours = Ivector(1,2*endbcs);
    for(i=1;i<=endbcs;i++)
      endparents[i] = endtypes[i] = endnodes[i] = endnodes2[i] = 0;

    endbcs = 0;
    for(i=1;i<=noknots;i++) {
      if(layernode[i] < 0) {
	endbcs++;
	endparents[endbcs] = i;
      }
    }
  }


  /* Check if the new boundary is already connected to some one, 
     however it must be different from the extruded boundary */
  for(i2=1;i2<=endbcs;i2++) {
    for(j=0;j<MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      
      for(i=1;i<=bound[j].nosides;i++) {
	
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,ind,&sideelemtype);
	
	/* Check that the node is one of the single nodes */
	dolayer = FALSE;
	for(i3=0;i3<sideelemtype%100;i3++)
	  if(ind[i3] == endparents[i2]) {
	    dolayer = TRUE;	
	    break;
	  }
	if(!dolayer) continue;
		
	/* First check that the found boundary has a correct parent material */
	dolayer = FALSE;
	if(checkmaterials) {
	  for(k=0;k<nolayers;k++) {
	    if(layerparents[k] < 0) continue;
	    parent = bound[j].parent[i];
	    if(data->material[parent] == layerparents[k])
	      dolayer = TRUE;
	    else if(parent = bound[j].parent2[i]) {
	      if(data->material[parent] == layerparents[k]) {
		dolayer = TRUE;
	      }
	    }  
	  }
	}
	if(!dolayer) continue;
	
	/* Finally check that this is not one of the extruded boundaries */
	dolayer = FALSE;
	for(k=0;k<nolayers;k++) {
	  if(layerparents[k] < 0) continue;
	  if(bound[j].types[i] == layerbounds[k]) dolayer = TRUE;
	}
	if(dolayer) {
	  endneighbours[2*i2-1] = ind[1-i3]; 	  
	  continue;
	}	

	endtypes[i2] = bound[j].types[i];
	dx = fabs(data->x[ind[0]] - data->x[ind[1]]);
	dy = fabs(data->y[ind[0]] - data->y[ind[1]]);

	if(dx < rectfactor * dy && dy < rectfactor * dx) {
	  endnodes[i2] = ind[i3];
	  if(sideelemtype%100 > 2) endnodes2[i2] = ind[2];
	  endneighbours[2*i2] = ind[1-i3]; 	    
	}

	if(info) printf("Found an existing boundary %d for the single node %d %d\n",
			bound[j].types[i],endparents[i2],endnodes[i2]); 	

	goto foundbc;
      }
    }

  foundbc:
    
    if(!endtypes[i2]) {
      maxbc++;
      endtypes[i2] = maxbc;
    }
  }


  /* Find the first unused bc */
  for(j=0;j<MAXBOUNDARIES;j++) 
    if(!bound[j].created) {
      newbc = j;
      bound[newbc].nosides = 0;
      break;
    }

  /* Find the maximum of layers */
  i = 0;
  for(k=0;k<nolayers;k++) 
    if(layernumber[k] > i) i = layernumber[k];

  if(endbcs) {
    if(info) {
      printf("Allocating for additional %d boundary elements into bc %d.\n",
	     bound[newbc].nosides,newbc);
    }  
    AllocateBoundary(&bound[newbc],i*endbcs);
    bound[newbc].created = FALSE;
    bound[newbc].nosides = 0;
  }


  /* The size of new mesh */
  noknots = data->noknots + nonewnodes;
  noelements = data->noelements + nonewelements;

  oldnoelements = data->noelements;
  oldnoknots = data->noknots;

  if(info) {
    printf("Creating additional %d elements and %d nodes.\n",nonewelements,nonewnodes);
    printf("Boundary layer mesh has %d elements and %d nodes.\n",noelements,noknots);
  }

  /* there will be more nodes if the original mesh consists of triangles */
  if(maxelemtype <= 303) 
    data->maxnodes = 4;
  else if(maxelemtype == 306)
    data->maxnodes = 8;

  /* Allocate more space for the enlarged data set */
  newtopo = Imatrix(1,noelements,0,data->maxnodes-1);
  newmaterial = Ivector(1,noelements);
  newelementtypes = Ivector(1,noelements);
  newx = Rvector(1,noknots);
  newy = Rvector(1,noknots);
  newz = Rvector(1,noknots);
  for(i=1;i<=noknots;i++) newz[i] = 0.0;

  elemwidth = Rvector(1,nonewelements);
  for(i=1;i<=nonewelements;i++) elemwidth[i] = 0.0;

  herit = Ivector(1,noknots);
  for(i=1;i<=oldnoknots;i++) herit[i] = i;
  for(i=oldnoknots+1;i<=noknots;i++) herit[i] = 0;


  /* Set the old topology */
  for(j=1;j<=data->noelements;j++) {
    newmaterial[j] = data->material[j];
    newelementtypes[j] = data->elementtypes[j];
    for(i=0;i<data->elementtypes[j]%100;i++) 
      newtopo[j][i] = data->topology[j][i];
  }

  /* Set the old nodes */
  for(i=1;i<=data->noknots;i++) {
    newx[i] = data->x[i];
    newy[i] = data->y[i];
  }

  topomap = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) topomap[i] = i;

  inside = Ivector(1,noelements);
  for(i=1;i<=noelements;i++) inside[i] = FALSE;

  /* Set the new node topology and nodes */
  elemindx = data->noelements;
  for(j=0;j<MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    
    for(i=1;i<=bound[j].nosides;i++) {

      dolayer = FALSE;
      parent = bound[j].parent[i];
      parent2 = bound[j].parent2[i];
      use2 = FALSE;

      for(k=0;k<nolayers;k++) {
	if(bound[j].types[i] == layerbounds[k]) {
	  if(checkmaterials) {
	    if(layerparents[k] < 0) continue;

	    if(data->material[parent] == layerparents[k]) {
	      dolayer = k + 1;
	    }
	    else if(parent2) {
	      l = parent;
	      parent = parent2;
	      parent2 = l;
	      if(data->material[parent] == layerparents[k]) {
		use2 = TRUE;
		dolayer = k + 1;
	      }
	    }
	  }
	  else dolayer = k + 1;
	}
      }


      if(!dolayer) continue;

      if(use2) 
	GetElementSide(bound[j].parent2[i],bound[j].side2[i],bound[j].normal[i],
		       data,ind,&sideelemtype);
      else
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,ind,&sideelemtype);

      inside[parent] = 1;

      if(sideelemtype == 202) 
	order = 1;
      else if(sideelemtype == 203) 
	order = 2;

      /* Check if some node should result into additional BC */
      for(i2=0;i2<sideelemtype%100;i2++) {
	sidebc[i2] = FALSE;
	if(i2 < 2 && layernode[ind[i2]] < 0) {
	  layernode[ind[i2]] = abs(layernode[ind[i2]]);
	  sidebc[i2] = TRUE;
	}
      }

      /* Define the normal of the surface */
      dy = -(data->x[ind[1]] - data->x[ind[0]]);
      dx = data->y[ind[1]] - data->y[ind[0]];
      ds = sqrt(dx*dx+dy*dy);
      dx /= ds;
      dy /= ds;

      n = layernumber[dolayer-1];
      ds = -layerthickness[dolayer-1];

      for(l=0;l < n;l++) {
	elemindx++;

	newmaterial[elemindx] = data->material[parent];
	inside[elemindx] = 1;

	if(n <= 1 || fabs(layerratios[dolayer-1]-1.0) < 0.001) {
	  q = (1.0*(l+1))/n;
	  elemwidth[elemindx-oldnoelements] = ds / n;	  
	}
	else {
	  ratio = pow(layerratios[dolayer-1],-1./(n-1.));
	  q = (1.- pow(ratio,(Real)(l+1))) / (1.-pow(ratio,(Real)(n)));
	  p = (1.- pow(ratio,(Real)(l))) / (1.-pow(ratio,(Real)(n)));
	  elemwidth[elemindx-oldnoelements] = (q-p) * ds;
	}


	for(m=0;m<sideelemtype%100;m++) {
	  
	  /* Make the possible additional BC appearing at side of the BL */
	  if(sidebc[m]) {
	    
	    bound[newbc].nosides += 1;
	    i2 = bound[newbc].nosides;
	    bound[newbc].parent[i2] = elemindx;
	    bound[newbc].parent2[i2] = 0;
	    bound[newbc].side[i2] = 3 - 2*m;
	    bound[newbc].side2[i2] = 0;	      
	    
	    for(i3=1;i3<=endbcs;i3++) 
	      if(ind[m] == endparents[i3]) {
		bound[newbc].types[i2] = endtypes[i3];		  
		endneighbours[2*i3-1] = layernode[ind[m]] + 1;
		break;
	      }
	  }

	  /* Set the node coordinates */
	  if(m < 2) {
	    nodeindx = layernode[ind[m]] + order*(l+1);	 
	  }
	  else {
	    nodeindx = layernode[ind[m]] + (1+midpoints)*(l+1);
	  }
	  e1x = dx * q * ds;
	  e1y = dy * q * ds;

	  /* Compute the normal of a joined node */
	  if(herit[nodeindx] != 0) {
	    
	    e2x = newx[nodeindx] - data->x[ind[m]];
	    e2y = newy[nodeindx] - data->y[ind[m]];
	      
	    p = (e1x*e2x + e1y*e2y)/(sqrt(e1x*e1x+e1y*e1y)*sqrt(e2x*e2x+e2y*e2y));

	    newx[nodeindx] += e1x - p * e2x;
	    newy[nodeindx] += e1y - p * e2y;
	  }
	  else {
	    herit[nodeindx] = ind[m];	    
	    newx[nodeindx] = data->x[ind[m]] + e1x;
	    newy[nodeindx] = data->y[ind[m]] + e1y;
	  }
	}

	/* Create the bulk elements */
	if(l==0) {
	  newtopo[elemindx][3] = ind[0];
	  newtopo[elemindx][2] = ind[1];
	  if(order == 2) newtopo[elemindx][6] = ind[2];
	}    
	else {
	  newtopo[elemindx][3] = layernode[ind[0]] + order*l;
	  newtopo[elemindx][2] = layernode[ind[1]] + order*l;
	  if(order == 2) newtopo[elemindx][6] = layernode[ind[2]] + (midpoints+1)*l;
	}
	newtopo[elemindx][0] = layernode[ind[0]] + order*(l+1);
	newtopo[elemindx][1] = layernode[ind[1]] + order*(l+1);

	if(order == 2) {	  
	  newtopo[elemindx][7] = layernode[ind[0]] + order*l+1;
	  newtopo[elemindx][5] = layernode[ind[1]] + order*l+1;	  
	  newtopo[elemindx][4] = layernode[ind[2]] + (midpoints+1)*(l+1);	  
	  if(midpoints) newtopo[elemindx][8] = layernode[ind[2]] + 2*l+1;
	}

	if(order == 1) {
	  newelementtypes[elemindx] = 404;	  
	}
	else if(midpoints) {
	  newelementtypes[elemindx] = 409;	    
	}
	else {
	  newelementtypes[elemindx] = 408;	    
	}


	if(l == n-1 && parent2) {

	  elemtype = data->elementtypes[parent2];	  
	  inside[parent2] = 2;
	  
	  for(i2=0;i2<elemtype%100;i2++) {
	    for(i3=0;i3<sideelemtype%100;i3++) {
	      if(data->topology[parent2][i2] == ind[i3]) {
		if(i3 < 2) {
		  topomap[ind[i3]] = layernode[ind[i3]] + order * n;
		}
		else {
		  topomap[ind[i3]] = layernode[ind[i3]] + (midpoints+1) * n;		  
		}
	      }
	    }
	  }
	}
      }

      /* Finally set the BC to point to the new boundary */
      if(use2) {
	bound[j].side2[i] = 0;
	bound[j].parent2[i] = elemindx;
      }
      else {
	bound[j].side[i] = 0;
	bound[j].parent[i] = elemindx;
      }
    }
  }


  {
    int *inside2;
    inside2 = Ivector(1,noknots);
    for(i=1;i<=noknots;i++) inside2[i] = 0;
    
    /* Put a marker to all nodes that belong to elements that are on the outside */
    for(j=1;j<=noelements;j++) {
      if(inside[j] == 2) {
	elemtype = data->elementtypes[j];
	for(i=0;i<elemtype/100;i++) {
	  inside2[newtopo[j][i]] = TRUE;
	}
      }
    }
    
    /* Now check other outside elements that have at least 2 nodes that are also on outside */
    for(j=1;j<=noelements;j++) {
      if(!inside[j]) {
	elemtype = data->elementtypes[j];
	k = 0;
	for(i=0;i<elemtype/100;i++) 
	  if(inside2[newtopo[j][i]]) k++;
	if(k > 1) inside[j] = 2;       
      }
    }
    free_Ivector(inside2,1,noknots);

     /* Still, go through all elements and if they are not on the list of
	active materials assume them outside */
    if(checkmaterials) {
      for(j=1;j<=oldnoelements;j++) {
	dolayer = FALSE;
	for(k=0;k<nolayers;k++) 
	  if(data->material[j] == layerparents[k]) dolayer = TRUE;

	if(!dolayer) {
	  if(inside[j] == 1) printf("Element %d of material %d should be in the inside\n",
				    j,data->material[j]);
	  inside[j] = 2;
	}
      }
    }
    
    /* And finally remap the nodes that are on the outside */
    for(j=1;j<=noelements;j++) {
      if(inside[j] == 2) {
	elemtype = data->elementtypes[j];
	for(i=0;i<elemtype%100;i++) 
	  newtopo[j][i] = topomap[data->topology[j][i]];
      }
    }
  }
  

  /* Put the pointers to the enlarged data set and destroy the old data */
  oldx = data->x;
  oldy = data->y;
  oldtopo = data->topology;

  data->noelements = noelements;
  data->noknots = noknots;
  data->x = newx;
  data->y = newy;
  data->z = newz;

  free_Ivector(data->elementtypes,1,oldnoelements);
  data->elementtypes = newelementtypes;

  free_Ivector(data->material,1,oldnoelements);
  data->material = newmaterial;
  data->topology = newtopo;


  /* In case one wants to fit the mesh inside the original mesh 
     the mesh nodes may be put to new positions using an appropriate filter. */


  /* For higher order elements remove the middlenodes from the list of cornernodes */
  if(maxelemtype%100 > 4) {
    if(info) printf("Marking the higher order nodes\n");

    nonlin = Ivector(1,noknots);
    for(i=1;i<=noknots;i++) nonlin[i] = FALSE;

    for(j=1;j<=noelements;j++) {      
      elemtype = data->elementtypes[j];	  
      for(i=elemtype/100;i<elemtype%100;i++) {
	k = data->topology[j][i];
	nonlin[k] = TRUE;
      }
    }
  }


  if(maxfilters) {
    int method,iter;
    int ind1,ind2,ind3,*fixedx=NULL,*fixedy=NULL;
    Real *aidx=NULL,*aidy=NULL,*weights=NULL;
    Real maxerror=0.0,minds,dx2,dy2,ds2,fii;
    
    /* There are three methods how to put the weight in the filter,
       1) 1/s, 2) fii/s, 3) sin(fii)/s, the second option seems to be best. */
    method = 2;
    
    if(info) printf("Filtering the mesh to meet the original geometry\n");
    
    fixedx = Ivector(1,noknots);
    fixedy = Ivector(1,noknots);
    weights = Rvector(1,noknots);
    aidx = Rvector(1,noknots);
    aidy = Rvector(1,noknots);
    
    /* Set all the fixed boundaries */
    for(i=1;i<=noknots;i++) fixedx[i] = fixedy[i] = 0;
    
    /* First, make all other materials except the ones with BL to be fixed */
    if(checkmaterials) {
      for(j=1;j<=noelements;j++) {
	
	elemtype = data->elementtypes[j];	  
	dolayer = FALSE;
	for(k=0;k<nolayers;k++) 
	  if(data->material[j] == layerparents[k]) dolayer = TRUE;
	
	for(i=0;i<elemtype/100;i++) {
	  ind1 = data->topology[j][i];
	  if(dolayer && fixedx[ind1]%2 == 0) 
	    fixedx[ind1] += 1;
	  if(!dolayer && fixedx[ind1] < 2) 
	    fixedx[ind1] += 2;
	}
      }
      for(i=1;i<=noknots;i++) {
	if(fixedx[i] == 2) 
	  fixedy[i] = 2;
	else
	  fixedx[i] = 0;
      }
    }
    
    /* Then set all BC:s fixed except the tangential ones */
    for(j=0;j<MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      
      for(i=1;i<=bound[j].nosides;i++) {
	
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,ind,&sideelemtype);
	
	dx = fabs(newx[ind[0]] - newx[ind[1]]);
	dy = fabs(newy[ind[0]] - newy[ind[1]]);
	if(dx > rectfactor * dy) {
	  for(l=0;l<sideelemtype%100;l++) {
	    fixedy[ind[l]] = TRUE;
	  }
	}
	else if(dy > rectfactor * dx) {
	  for(l=0;l<sideelemtype%100;l++) {
	    fixedx[ind[l]] = TRUE;
	  }
	}
	else {
	  for(l=0;l<sideelemtype%100;l++) {
	    fixedy[ind[l]] = TRUE;
	    fixedx[ind[l]] = TRUE;
	  }
	}	
      }
    }
    
    
    /* Then set possibly all remaining active boundaries to be fixed */
    for(j=0;j<MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      
      for(i=1;i<=bound[j].nosides;i++) {
	
	dolayer = FALSE;
	parent = bound[j].parent[i];
	parent2 = bound[j].parent2[i];
	use2 = FALSE;
	
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],
		       data,ind,&sideelemtype);
	
	for(k=0;k<nolayers;k++) {
	  if(bound[j].types[i] == layerbounds[k]) {
	    if(checkmaterials) {
	      if(layerparents[k] < 0) continue;
	      
	      if(data->material[parent] == layerparents[k]) {
		dolayer = k + 1;
	      }
	      else if(parent2) {
		l = parent;
		parent = parent2;
		parent2 = l;
		if(data->material[parent] == layerparents[k]) {
		  use2 = TRUE;
		  dolayer = k + 1;
		}
	      }
	    }
	    else dolayer = k + 1;
	  }
	}
	
	if(dolayer) {
	  for(l=0;l<sideelemtype%100;l++) {
	    fixedy[ind[l]] = TRUE;
	    fixedx[ind[l]] = TRUE;
	  }
	}
      }
    }

    /* Finally loose the problematic triple nodes */
    for(j=1;j<=endbcs;j++) {
      k = endnodes[j];
      if(k) {
	fixedx[k] = FALSE;
	fixedy[k] = FALSE;
      }
      
      /* for second order elements */
      k = endnodes2[j];
      if(k) {
	fixedx[k] = FALSE;
	fixedy[k] = FALSE;
      }
    }
    
    
    j = 0;
    for(i=1;i<=noknots;i++) if(fixedx[i]) j += 1;
    if(info) printf("Number of fixed nodes in x-direction is %d\n",j);
    
    j = 0;
    for(i=1;i<=noknots;i++) if(fixedy[i]) j += 1;
    if(info) printf("Number of fixed nodes in y-direction is %d\n",j);
    
    for(j=1;j<=noknots;j++) {	  

      if(fixedx[j]) {
	if(j <= oldnoknots) 
	  newx[j] = aidx[j] = oldx[j]; 
	else 
	  newx[j] = aidx[j] = oldx[herit[j]]; 
      }
      if(fixedy[j]) {
	if(j <= oldnoknots)
	  newy[j] = aidy[j] = oldy[j];
	else 
	  newy[j] = aidy[j] = oldy[herit[j]];
      }
    }	    
    
    
    for(iter=1;iter<=maxfilters;iter++) {
      maxerror = 0.0;
      minds = 1.0e10;
      
      for(j=1;j<=noknots;j++) {	  
	
	weights[j] = 0.0;
	
	if(!fixedx[j]) {	  
	  aidx[j] = newx[j];
	  newx[j] = 0.0;
	}
	if(!fixedy[j]) {
	  aidy[j] = newy[j];	  
	  newy[j] = 0.0;
	}
      }
      
      for(j=1;j<=noelements;j++) {
	elemtype = data->elementtypes[j];
	nonodes = elemtype / 100;
	
	for(i=0;i<nonodes;i++) {
	  
	  i2 = (i+1)%nonodes;
	  i3 = (i+2)%nonodes;
	  
	  ind1 = data->topology[j][i];
	  ind2 = data->topology[j][i2];
	  ind3 = data->topology[j][i3];
	  
	  if(j<=oldnoelements) {
	    dx = oldx[oldtopo[j][i2]] - oldx[oldtopo[j][i]];
	    dy = oldy[oldtopo[j][i2]] - oldy[oldtopo[j][i]];
	    ds = sqrt(dx*dx+dy*dy);
	  }
	  else {
	    ds = fabs(elemwidth[j-oldnoelements]);
	  }
	  if(ds < minds) minds = ds;

	  
	  if(j<=oldnoelements) {
	    dx2 = oldx[oldtopo[j][i2]] - oldx[oldtopo[j][i3]];
	    dy2 = oldy[oldtopo[j][i2]] - oldy[oldtopo[j][i3]];
	    ds2 = sqrt(dx2*dx2+dy2*dy2);
	  }
	  else {
	    ds2 = fabs(elemwidth[j-oldnoelements]);
	  }
	  
	  if(j <= oldnoelements && ds * ds2 < 1.0e-50) {
	    printf("problem elem %d and nodes %d (%d %d)\n",j,i2,i,i3);
	    printf("dist ds=%.3e ds2=%.3e\n",ds,ds2);
	    printf("coord: %.3e %.3e\n",oldx[oldtopo[j][i2]], oldy[oldtopo[j][i2]]);
	    continue;
	  }
	  
	  if(abs(method) == 2 && j<=oldnoelements) {
	    fii = acos((dx*dx2+dy*dy2)/(ds*ds2)) / (FM_PI/2.0);
	  }
	  else if(abs(method) == 3 && j<=oldnoelements) {
	    fii = acos((dx*dx2+dy*dy2)/(ds*ds2));
	    fii = sin(fii);
	  }
	  else {
	    fii = 1.0;
	  }
	  
	  
	  /* Eliminate the very difficult triple nodes */
	  dolayer = FALSE;
	  for(k=1;k<=endbcs;k++)
	    if(ind2 == endnodes[k]) dolayer = k;
	  
	  if(dolayer) {
	    for(k=1;k<=2;k++) {
	      if(endneighbours[2*(dolayer-1)+k] == ind1) {
		weights[ind2] += fii / ds;		    
		if(!fixedx[ind2]) newx[ind2] += aidx[ind1] * fii / ds;	      	      
		if(!fixedy[ind2]) newy[ind2] += aidy[ind1] * fii / ds;	      
	      }
	    }	
	    for(k=1;k<=2;k++) {
	      if(endneighbours[2*(dolayer-1)+k] == ind3) {
		weights[ind2] += fii / ds2;
		if(!fixedx[ind2]) newx[ind2] += aidx[ind3] * fii / ds2;	      
		if(!fixedy[ind2]) newy[ind2] += aidy[ind3] * fii / ds2;
	      }
	    }	
	  }		
	  else {
	    if(ind2 <= oldnoknots || herit[ind1] == herit[ind2]) {
	      weights[ind2] += fii / ds;		    
	      if(!fixedx[ind2]) newx[ind2] += aidx[ind1] * fii / ds;	      	      
	      if(!fixedy[ind2]) newy[ind2] += aidy[ind1] * fii / ds;	      
	    }
	    
	    if(ind2 <= oldnoknots || herit[ind3] == herit[ind2]) {
	      weights[ind2] += fii / ds2;
	      if(!fixedx[ind2]) newx[ind2] += aidx[ind3] * fii / ds2;	      
	      if(!fixedy[ind2]) newy[ind2] += aidy[ind3] * fii / ds2;
	    }
	  }
	}
      }
      
      if(maxelemtype%100 > 4) {  
	for(j=1;j<=noknots;j++) {	  
	  if(nonlin[j]) continue;

	  if(weights[j] > 1.0e-50) {	  
	    if(!fixedx[j]) newx[j] /= weights[j];
	    if(!fixedy[j]) newy[j] /= weights[j];
	  }
	  else if(iter==1) {
	    printf("no weight for index %d\n",j);
	  }
	  
	  dx = newx[j] - aidx[j];
	  dy = newy[j] - aidy[j];

	  ds = dx*dx + dy*dy;
	  if(ds > maxerror) maxerror = ds;
	}
      } 
      else {
	for(j=1;j<=noknots;j++) {	  
	  if(!fixedx[j]) newx[j] /= weights[j];
	  if(!fixedy[j]) newy[j] /= weights[j];
	  
	  dx = newx[j]-aidx[j];
	  dy = newy[j]-aidy[j];

	  ds = dx*dx + dy*dy;
	  if(ds > maxerror) maxerror = ds;
	}
      }

      maxerror = sqrt(maxerror) / minds;
      if(maxerror < layereps) break;
    }

    if(info) {
      printf("Filtered the new node coordinates %d times with final error %.3e.\n",
	     iter-1,maxerror);
    }
    
    /* In higher order elements map the middle nodes so that they lie in between
       the corner nodes */

    
    if(maxelemtype%100 > 4) {
      for(j=1;j<=noelements;j++) {
	elemtype = data->elementtypes[j];
	if(elemtype%100 <= elemtype/100) continue;
	
	if(elemtype == 306) {
	  for(k=0;k<3;k++) {
	    if(!fixedx[newtopo[j][k+3]]) {
	      newx[newtopo[j][k+3]] = 0.5 * (newx[newtopo[j][k]] + newx[newtopo[j][(k+1)%3]]);
	    }
	    if(!fixedy[newtopo[j][k+3]]) {
	      newy[newtopo[j][k+3]] = 0.5 * (newy[newtopo[j][k]] + newy[newtopo[j][(k+1)%3]]);
	    }
	  }
	}
	
	else if(elemtype == 408 || elemtype == 409) {
	  
	  if(elemtype == 409) {
	    newx[newtopo[j][8]] = 0.0;
	    newy[newtopo[j][8]] = 0.0;
	  }
	  
	  for(k=0;k<4;k++) {
	    if(!fixedx[newtopo[j][k+4]]) {
	      newx[newtopo[j][k+4]] = 0.5 * (newx[newtopo[j][k]] + newx[newtopo[j][(k+1)%4]]);
	    }
	    if(!fixedy[newtopo[j][k+4]]) {
	      newy[newtopo[j][k+4]] = 0.5 * (newy[newtopo[j][k]] + newy[newtopo[j][(k+1)%4]]);
	    }
	    if(elemtype == 409) {
	      newx[newtopo[j][8]] += 0.25 * newx[newtopo[j][k]];
	      newy[newtopo[j][8]] += 0.25 * newy[newtopo[j][k]];
	    }
	  }
	}
	else {
	  printf("Unknown elementtype %d\n",elemtype);
	}
      }
    }
  
    free_Ivector(fixedx,1,noknots);
    free_Ivector(fixedy,1,noknots);
    
    free_Rvector(aidx,1,noknots);
    free_Rvector(aidy,1,noknots);   
    free_Rvector(weights,1,noknots);
  }

  if(bound[newbc].nosides > 0) 
    bound[newbc].created = TRUE;
  
  
  /* In higher order elements map the middle nodes so that they lie in between
     the corner nodes. Elemtypes must be 408 or 409 since they are created in this
     subroutine */
  
  if(!maxfilters && maxelemtype%100 > 4) {
    if(info) printf("Making the higher order nodes to lie in between\n");
    
    for(j=oldnoelements+1;j<=noelements;j++) {
      
      elemtype = data->elementtypes[j];
      if(elemtype%100 <= elemtype/100) continue;
      
      if(elemtype == 408 || elemtype == 409) {
	
	if(elemtype == 409) {
	  newx[newtopo[j][8]] = 0.0;
	  newy[newtopo[j][8]] = 0.0;
	}
	
	for(k=0;k<4;k++) {
	  newx[newtopo[j][k+4]] = 0.5 * (newx[newtopo[j][k]] + newx[newtopo[j][(k+1)%4]]);
	  newy[newtopo[j][k+4]] = 0.5 * (newy[newtopo[j][k]] + newy[newtopo[j][(k+1)%4]]);
	  
	  if(elemtype == 409) {
	    newx[newtopo[j][8]] += 0.25 * newx[newtopo[j][k]];
	    newy[newtopo[j][8]] += 0.25 * newy[newtopo[j][k]];
	  }
	}
      }
    }
  }


#if 0
  ReorderElements(data,bound,FALSE,corder,info);
#endif

  free_Imatrix(oldtopo,1,oldnoelements,0,oldmaxnodes-1);
  free_Ivector(layernode,1,oldnoknots);
  free_Rvector(oldx,1,oldnoknots);
  free_Rvector(oldy,1,oldnoknots);

  if(info) printf("Boundary layers created successfully.\n");

  if(checkmaterials && !second) {    
    for(k=0;k<nolayers;k++) {
      if(layerparents[k] < 0) second = TRUE;
      layerparents[k] = -layerparents[k];
    }
    if(second) {
      if(info) printf("\nPerforming boundary layer generation again for negative materials\n");
      goto omstart;
    }
  }

  return(0);
}





int CreateBoundaryLayerDivide(struct FemType *data,struct BoundaryType *bound,
			      int nolayers, int *layerbounds, int *layernumber,
			      Real *layerratios, Real *layerthickness, int *layerparents,
			      int info)
/* Create Boundary layers that may be used to solve accurately fluid
   flow problems and similar equations. In this subroutine the boundary layer
   is created by dividing the elements close to boundary. */
{
  int i,j,k,l,dim,maxbc,maxelemtype,dolayer,parent,nlayer,sideelemtype,elemind,side;
  int noelements,noknots,oldnoknots,oldnoelements,oldmaxnodes,nonewnodes,nonewelements;
  int maxcon,elemsides,elemdone,midpoints,order,bcnodes,elemhits,elemtype,goforit;
  int ind[MAXNODESD2],baseind[2],topnode[2],basenode[2];
  int *layernode=NULL,*newelementtypes=NULL,**newtopo=NULL,**oldtopo=NULL;
  int *newmaterial=NULL,**edgepairs=NULL,*sharednode=NULL;
  Real dx[2],dy[2],x0[2],y0[2];
  Real *newx=NULL,*newy=NULL,*newz=NULL,*oldx=NULL,*oldy=NULL,*oldz=NULL;
  Real slayer,qlayer,ratio,q;


  dim = data->dim;
  maxelemtype = GetMaxElementType(data);

  if(maxelemtype > 409) {
    printf("Subroutine implemented only up to 2nd degree!\n");
    return(2);
  }

  if(info) printf("Largest elementtype is %d\n",maxelemtype);
  

  oldnoelements = noelements = data->noelements;
  oldnoknots = noknots = data->noknots;
  oldmaxnodes = data->maxnodes;

  layernode = Ivector(1,oldnoknots);
  for(i=1;i<=oldnoknots;i++) 
    layernode[i] = 0;

  sharednode = Ivector(1,oldnoknots);
  for(i=1;i<=oldnoknots;i++) 
    sharednode[i] = 0;


  /* Go through all the boundaries with boundary layer definitions and compute 
     the number of nodes at the surface. */

  maxbc = 0;
  qlayer = 0.0;
  slayer = 0.0;
  nlayer = 0;

  /* Go through the layers and check which ones are active */
  for(j=0;j<MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;

    for(i=1;i<=bound[j].nosides;i++) {
      dolayer = FALSE;
      parent = bound[j].parent[i];
      if(bound[j].types[i] > maxbc) maxbc = bound[j].types[i];

      for(k=0;k<nolayers;k++) {
	if(bound[j].types[i] == layerbounds[k]) {
	  nlayer = layernumber[k];
	  slayer = layerthickness[k];
	  qlayer = layerratios[k];	  
	  dolayer = TRUE;
	}
      }
      if(!dolayer) continue;

      /* We have found an active boundary layer */
      GetElementSide(parent,bound[j].side[i],bound[j].normal[i],
		     data,ind,&sideelemtype);

      midpoints = FALSE;
      if(sideelemtype == 202) {
	order = 1;
      }
      else if(sideelemtype == 203) {
	order = 2;
	if(maxelemtype > 408) midpoints = TRUE;
      }
      
      for(l=0;l<sideelemtype%100;l++)  
	layernode[ind[l]] += 1;
    }
  }

  if(slayer > 1.0 || slayer < 1.0e-20) 
    slayer = 1.0;

  bcnodes = 0;
  maxcon = 0;
  for(i=1;i<=data->noknots;i++) {
    if(layernode[i]) bcnodes++;
    maxcon = MAX(maxcon, layernode[i]);
  }  

  if(info) printf("Found %d new nodes in the boundary layers!\n",bcnodes);
  if(!bcnodes) return(0);  
  if(info) printf("There are %d connections at maximum\n",maxcon);

  /* there will be more nodes if the original mesh consists of triangles */
  if(maxelemtype <= 303) 
    data->maxnodes = 4;
  else if(maxelemtype == 306)
    data->maxnodes = 8;

  /* Compute the number of new elements */
  nonewelements = 0;
  for(j=1;j<=data->noelements;j++) {
    elemhits = 0;
    elemtype = data->elementtypes[j];
    for(i=0;i<elemtype%100;i++) {
      k = data->topology[j][i];
      if( layernode[k]) {
	sharednode[k] += 1;
	elemhits++;
      }
    }
    if(elemhits) {
      nonewelements += nlayer ;
      if(elemhits != 2) nonewelements += nlayer + 1;
    }
  }
  printf("There will %d new elements\n",nonewelements);

  /* This is a conservative estimate */
  nonewnodes = 2*nonewelements;

  edgepairs = Imatrix(1,nonewnodes,1,3);
  for(j=1;j<=nonewnodes;j++) 
    edgepairs[j][1] = edgepairs[j][2] = edgepairs[j][3] = 0;


  /* The size of new mesh */
  oldnoelements = data->noelements;
  oldnoknots = data->noknots;
  oldtopo = data->topology;
  oldx = data->x;
  oldy = data->y;
  oldz = data->z;

  noknots = oldnoknots + nonewnodes;
  noelements = oldnoelements + nonewelements;

  if(info) {
    printf("Creating additional %d elements and %d nodes.\n",nonewelements,nonewnodes);
    printf("Boundary layer mesh has at maximum %d elements and %d nodes.\n",noelements,noknots);
  }

  /* Allocate more space for the enlarged data set */
  newtopo = Imatrix(1,noelements,0,data->maxnodes-1);
  newmaterial = Ivector(1,noelements);
  newelementtypes = Ivector(1,noelements);
  newx = Rvector(1,noknots);
  newy = Rvector(1,noknots);
  newz = Rvector(1,noknots);

  /* Set the old topology */
  for(j=1;j<=data->noelements;j++) {
    newmaterial[j] = data->material[j];
    newelementtypes[j] = data->elementtypes[j];
    for(i=0;i<data->elementtypes[j]%100;i++) 
      newtopo[j][i] = data->topology[j][i];
  }

  /* Set the old nodes */
  for(i=1;i<=data->noknots;i++) {
    newx[i] = data->x[i];
    newy[i] = data->y[i];
    newz[i] = data->z[i];
  }

  noelements = data->noelements;
  elemind = noelements;
  noknots = data->noknots;

  /* Go through elements and make the new elements and nodes */
  for(j=1;j<=data->noelements;j++) {
    elemhits = 0;
    elemtype = data->elementtypes[j];
    elemsides = elemtype % 100;
    for(i=0;i<elemsides;i++) 
      if( layernode[ data->topology[j][i] ]) elemhits++;
    if(!elemhits) continue;
    
    if(elemtype == 404) {      
      elemdone = FALSE;

      for(side=0;side<elemsides;side++) {
	
	goforit = FALSE;	
	if(elemhits == 2 || elemhits == 3) 
	  if(layernode[oldtopo[j][side]] && layernode[oldtopo[j][(side+1)%elemsides]]) goforit = TRUE;
	if(elemhits == 1) 
	  if(layernode[oldtopo[j][side]]) goforit = TRUE;
	if(!goforit) continue;

	/* Treat the special case of three hits 
	   In case of corners find the single node that is not on the boundary */
	if(elemhits == 3) {
	  for(k=0;k<4;k++)
	    if(!layernode[oldtopo[j][k]]) break; 
	  if(0) printf("Special node %d in corner %d\n",oldtopo[j][k],k);
	  
	  basenode[0] = oldtopo[j][side];
	  basenode[1] = oldtopo[j][(side+1)%elemsides];
	  topnode[0] = oldtopo[j][k];
	  topnode[1] = oldtopo[j][k];
	}
	else if(elemhits == 2) {
	  basenode[0] = oldtopo[j][side];
	  basenode[1] = oldtopo[j][(side+1)%elemsides];
	  topnode[0] = oldtopo[j][(side+3)%elemsides];	 	  
	  topnode[1] = oldtopo[j][(side+2)%elemsides];	 	  
	}	    
	else if(elemhits == 1) {
	  basenode[0] = oldtopo[j][side];
	  basenode[1] = basenode[0];
	  topnode[0] = oldtopo[j][(side+3)%elemsides];	 	  
	  topnode[1] = oldtopo[j][(side+1)%elemsides];	 	  
	}
	
	for(k=0;k<=1;k++) {
	  for(i=1;i<=nonewnodes;i++) { 
	    if(!edgepairs[i][1]) break;	    
	    if(basenode[k] == edgepairs[i][1] && topnode[k] == edgepairs[i][2]) break;
	  }
	  if(!edgepairs[i][1]) {
	    edgepairs[i][1] = basenode[k];
	    edgepairs[i][2] = topnode[k];
	    baseind[k] = noknots;
	    edgepairs[i][3] = baseind[k];
	    noknots += nlayer;
	  }
	  else {
	    if(0) printf("Using existing nodes\n");
	    baseind[k] = edgepairs[i][3];
	  }
	  x0[k] = oldx[basenode[k]];
	  y0[k] = oldy[basenode[k]];
	  dx[k] = oldx[topnode[k]] - x0[k];
	  dy[k] = oldy[topnode[k]] - y0[k];
	
	  for(i=1;i<=nlayer;i++) {
	    if(nlayer <= 1 || fabs(qlayer-1.0) < 0.001) {
	      q = (1.0*i) / (nlayer+1);
	    }
	    else {
	      ratio = pow(qlayer,1.0/nlayer);
	      q = (1.- pow(ratio,1.0*i)) /  (1.- pow(ratio,1.0+nlayer));
	    }	      
	    q *= slayer; 
	    newx[baseind[k]+i] = x0[k] + q * dx[k];
	    newy[baseind[k]+i] = y0[k] + q * dy[k];
	  }
	}	

	/* 0:th element */
	if(elemhits == 1) {
	  newelementtypes[j] = 303;	  
	  newtopo[j][0] = basenode[0];
	  newtopo[j][1] = baseind[1] + 1;
	  newtopo[j][2] = baseind[0] + 1;	  	  	  
	}
	else if(elemhits == 3 && elemdone) {
	  elemind++;
	  newelementtypes[elemind] = 404;
	  newmaterial[elemind] = newmaterial[j];	    
	  newtopo[elemind][side] = basenode[0];
	  newtopo[elemind][(side+1)%elemsides] = basenode[1];
	  newtopo[elemind][(side+2)%elemsides] = baseind[1] + 1;
	  newtopo[elemind][(side+3)%elemsides] = baseind[0] + 1;
	}
	else {
	  newtopo[j][(side+2)%elemsides] = baseind[1] + 1;
	  newtopo[j][(side+3)%elemsides] = baseind[0] + 1;
	}
	
	for(i=1;i<nlayer;i++) {
	  elemind++;
	  newelementtypes[elemind] = 404;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = baseind[0] + i;
	  newtopo[elemind][1] = baseind[1] + i;
	  newtopo[elemind][2] = baseind[1] + i+1;
	  newtopo[elemind][3] = baseind[0] + i+1;
	}
	
	/* n:th element */
	if(elemhits == 3) {
	  elemind++;
	  newelementtypes[elemind] = 303;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = baseind[0] + nlayer;
	  newtopo[elemind][1] = baseind[1] + nlayer;
	  newtopo[elemind][2] = topnode[0];
	}
	else if(elemhits == 2 || elemhits == 1) {
	  elemind++;
	  newelementtypes[elemind] = 404;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = baseind[0] + nlayer;
	  newtopo[elemind][1] = baseind[1] + nlayer;
	  newtopo[elemind][2] = topnode[1];
	  newtopo[elemind][3] = topnode[0];
	}	    
	/* n+1:th element */
	if(elemhits == 1) {
	  elemind++;
	  newelementtypes[elemind] = 303;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = topnode[1];
	  newtopo[elemind][1] = oldtopo[j][(side+2)%elemsides];
	  newtopo[elemind][2] = topnode[0];
	}
	
	elemdone = TRUE;
      }
      if(!elemdone) 
	printf("cannot handle quadrilaterals with %d hits\n",elemhits);
    }


    else if(elemtype == 303) {      
      elemdone = FALSE;

      for(side=0;side<elemsides;side++) {	

	goforit = FALSE;	
	if(elemhits == 2) {
	  if(layernode[oldtopo[j][side]] && layernode[oldtopo[j][(side+1)%elemsides]]) goforit = TRUE;
	}
	else if(elemhits == 1) {
	  if(layernode[oldtopo[j][side]]) goforit = TRUE;
	}
	else if(elemhits == 3) {
	  if(sharednode[oldtopo[j][side]] == 1) goforit = TRUE;

	  printf("The boundary layer creation for certain corner triangles is omitted\n");
	  goforit = FALSE;
	}
	if(!goforit) continue;

	if(elemhits == 3) {
	  if(1) printf("Special node %d in corner %d\n",oldtopo[j][side],side);	  
	  basenode[0] = oldtopo[j][side];
	  basenode[1] = basenode[0];
	  topnode[0] = oldtopo[j][(side+2)%elemsides];
	  topnode[1] = oldtopo[j][(side+1)%elemsides];
	}
	else if(elemhits == 2) {
	  basenode[0] = oldtopo[j][side];
	  basenode[1] = oldtopo[j][(side+1)%elemsides];
	  topnode[0] = oldtopo[j][(side+2)%elemsides];	 	  
	  topnode[1] = topnode[0];
	}	    
	else if(elemhits == 1) {
	  basenode[0] = oldtopo[j][side];
	  basenode[1] = basenode[0];
	  topnode[0] = oldtopo[j][(side+2)%elemsides];	 	  
	  topnode[1] = oldtopo[j][(side+1)%elemsides];	 	  
	}
	
	for(k=0;k<=1;k++) {
	  for(i=1;i<=nonewnodes;i++) { 
	    if(!edgepairs[i][1]) break;	    
	    if(basenode[k] == edgepairs[i][1] && topnode[k] == edgepairs[i][2]) break;
	  }
	  if(!edgepairs[i][1]) {
	    edgepairs[i][1] = basenode[k];
	    edgepairs[i][2] = topnode[k];
	    baseind[k] = noknots;
	    edgepairs[i][3] = baseind[k];
	    noknots += nlayer;
	  }
	  else {
	    if(0) printf("Using existing nodes\n");
	    baseind[k] = edgepairs[i][3];
	  }

	  x0[k] = oldx[basenode[k]];
	  y0[k] = oldy[basenode[k]];
	  dx[k] = oldx[topnode[k]] - x0[k];
	  dy[k] = oldy[topnode[k]] - y0[k];
	
	  for(i=1;i<=nlayer;i++) {
	    if(nlayer <= 1 || fabs(qlayer-1.0) < 0.001) {
	      q = (1.0*i) / (nlayer+1);
	    }
	    else {
	      ratio = pow(qlayer,1.0/nlayer);
	      q = (1.- pow(ratio,1.0*i)) /  (1.- pow(ratio,1.0*nlayer));
	    }	      
	    q *= slayer;
	    newx[baseind[k]+i] = x0[k] + q * dx[k];
	    newy[baseind[k]+i] = y0[k] + q * dy[k];
	  }
	}	

	/* 0:th element */
	if(elemhits == 1 || elemhits == 3) {
	  newelementtypes[j] = 303;	  
	  newtopo[j][0] = basenode[0];
	  newtopo[j][1] = baseind[1] + 1;
	  newtopo[j][2] = baseind[0] + 1;	  	  	  
	}
	else if(elemhits == 2) {
	  newelementtypes[j] = 404;	  
	  newtopo[j][side] = basenode[0];
	  newtopo[j][(side+1)%4] = basenode[1];
	  newtopo[j][(side+2)%4] = baseind[1] + 1;
	  newtopo[j][(side+3)%4] = baseind[0] + 1;	  	  	  
	}

	for(i=1;i<nlayer;i++) {
	  elemind++;
	  newelementtypes[elemind] = 404;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = baseind[0] + i;
	  newtopo[elemind][1] = baseind[1] + i;
	  newtopo[elemind][2] = baseind[1] + i+1;
	  newtopo[elemind][3] = baseind[0] + i+1;
	}
	
	/* n:th element */
	if(elemhits == 1 || elemhits == 3) {
	  elemind++;
	  newelementtypes[elemind] = 404;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = baseind[0] + nlayer;
	  newtopo[elemind][1] = baseind[1] + nlayer;
	  newtopo[elemind][2] = topnode[1];
	  newtopo[elemind][3] = topnode[0];
	}	    
	else if(elemhits == 2) {
	  elemind++;
	  newelementtypes[elemind] = 303;
	  newmaterial[elemind] = newmaterial[j];
	  newtopo[elemind][0] = baseind[0] + nlayer;
	  newtopo[elemind][1] = baseind[1] + nlayer;
	  newtopo[elemind][2] = topnode[1];
	}
	elemdone = TRUE;
      }
      if(!elemdone) 
	printf("cannot handle triangles with %d hits\n",elemhits);
    }

    else {
      printf("Not implemented for element %d\n",elemtype);
    }
  }
  noelements = elemind;
 
  data->x = newx;
  data->y = newy;
  data->topology = newtopo;
  data->material = newmaterial;
  data->elementtypes = newelementtypes;
  data->noknots = noknots;
  data->noelements = elemind;

  printf("The created boundary layer mesh has at %d elements and %d nodes.\n",noelements,noknots);

  return(0);
}



int RotateTranslateScale(struct FemType *data,struct ElmergridType *eg,int info)
{
  int i;
  Real x,y,z,xz,yz,yx,zx,zy,xy,cx,cy,cz;
  Real xmin, xmax, ymin, ymax, zmin, zmax;

  if(eg->scale) {
    if(info) printf("Scaling mesh with vector [%.3lg %.3lg %.3lg]\n",
	   eg->cscale[0],eg->cscale[1],eg->cscale[2]);
    for(i=1;i<=data->noknots;i++) {
      data->x[i] *= eg->cscale[0]; 
      data->y[i] *= eg->cscale[1]; 
      data->z[i] *= eg->cscale[2]; 
    }
    if(0) printf("Scaling of mesh finished.\n");
  }
  
  if(eg->rotate) {
    if(info) printf("Rotating mesh with degrees [%.3lg %.3lg %.3lg]\n",
		    eg->crotate[0],eg->crotate[1],eg->crotate[2]);
    cx = FM_PI * eg->crotate[0]/180.0;
    cy = FM_PI * eg->crotate[1]/180.0;
    cz = FM_PI * eg->crotate[2]/180.0;

    for(i=1;i<=data->noknots;i++) {

      x = data->x[i];
      y = data->y[i];
      z = data->z[i];

      xz = x*cos(cz) + y*sin(cz);
      yz = -x*sin(cz) + y*cos(cz);
      
      if( fabs(cx) > 1.0e-8 || fabs(cy) > 1.0e-8 ) {
	yx = yz*cos(cx) + z*sin(cx);
	zx = -yz*sin(cx) + z*cos(cx);
	
	zy = zx*cos(cy) + xz*sin(cy);
	xy = -zx*sin(cy) + xz*cos(cy);
	
	data->x[i] = xy;
	data->y[i] = yx;
	data->z[i] = zy;
      }	
      else {
	data->x[i] = xz;
	data->y[i] = yz;  
      }
    }
    if(0) printf("Rotation of mesh finished.\n");
  }

  if(eg->translate) {
    if(info) printf("Translating the mesh with vector [%.3lg %.3lg %.3lg]\n",
		    eg->ctranslate[0],eg->ctranslate[1],eg->ctranslate[2]);
    for(i=1;i<=data->noknots;i++) {
      data->x[i] += eg->ctranslate[0];
      data->y[i] += eg->ctranslate[1];
      data->z[i] += eg->ctranslate[2];
    }
    if(0) printf("Translation of mesh finished.\n");
  }

  if(eg->center) {
    xmin = xmax = data->x[1];
    ymin = ymax = data->y[1];
    zmin = zmax = data->z[1];

    for(i=1;i<=data->noknots;i++) {
      xmax = MAX( xmax, data->x[i] );
      xmin = MIN( xmin, data->x[i] );
      ymax = MAX( ymax, data->y[i] );
      ymin = MIN( ymin, data->y[i] );
      zmax = MAX( zmax, data->z[i] );
      zmin = MIN( zmin, data->z[i] );
    }
    cx = 0.5 * (xmin + xmax);
    cy = 0.5 * (ymin + ymax);
    cz = 0.5 * (zmin + zmax);

    if(info) printf("Setting new center to %.3e %.3e %.3e\n",cx,cy,cz);

    for(i=1;i<=data->noknots;i++) {
      data->x[i] -= cx;
      data->y[i] -= cy;
      data->z[i] -= cz;
    }    
  }

  return(0);
}



int CreateNodalGraph(struct FemType *data,int full,int info)
{
  int i,j,k,l,m,totcon,noelements, noknots,elemtype,nonodes,hit,ind,ind2;
  int maxcon,percon,edge;

  printf("Creating a nodal graph of the finite element mesh\n");  

  if(data->nodalexists) {
    printf("The nodal graph already exists!\n");
    smallerror("Nodal graph not done");
    return(1);
  }

  maxcon = 0;
  totcon = 0;
  percon = 0;
  noelements = data->noelements;
  noknots = data->noknots;

  for(i=1;i<=noelements;i++) {
    elemtype = data->elementtypes[i];

    /* This sets only the connections resulting from element edges */
    if(!full) {
      int inds[2];
      for(edge=0;;edge++) {
	if( !GetElementGraph(i,edge,data,&inds[0]) ) break;

	ind = inds[0];
	ind2 = inds[1];

	hit = FALSE;
	for(l=0;l<maxcon;l++) { 
	  if(data->nodalgraph[l][ind] == ind2) hit = TRUE;
	  if(data->nodalgraph[l][ind] == 0) break;
	}
	if(!hit) {
	  if(l >= maxcon) {
	    data->nodalgraph[maxcon] = Ivector(1,noknots);
	    for(m=1;m<=noknots;m++)
	      data->nodalgraph[maxcon][m] = 0;
	    maxcon++;
	  }
	  data->nodalgraph[l][ind] = ind2;
	  totcon++;
	}

	/* Make also so symmetric connection */
	for(l=0;l<maxcon;l++) { 
	  if(data->nodalgraph[l][ind2] == ind) hit = TRUE;
	  if(data->nodalgraph[l][ind2] == 0) break;
	}
	if(!hit) {
	  if(l >= maxcon) {
	    data->nodalgraph[maxcon] = Ivector(1,noknots);
	    for(m=1;m<=noknots;m++)
	      data->nodalgraph[maxcon][m] = 0;
	    maxcon++;
	  }
	  data->nodalgraph[l][ind2] = ind;
	  totcon++;
	}
      }
    }

    /* This sets all elemental connections */
    else {
      nonodes = data->elementtypes[i] % 100;
      for(j=0;j<nonodes;j++) {
	ind = data->topology[i][j];
	for(k=0;k<nonodes;k++) {
	  ind2 = data->topology[i][k];
	  if(ind == ind2) continue;
	  
	  hit = FALSE;
	  for(l=0;l<maxcon;l++) { 
	    if(data->nodalgraph[l][ind] == ind2) hit = TRUE;
	    if(data->nodalgraph[l][ind] == 0) break;
	  }
	  if(!hit) {
	    if(l >= maxcon) {
	      data->nodalgraph[maxcon] = Ivector(1,noknots);
	      for(m=1;m<=noknots;m++)
		data->nodalgraph[maxcon][m] = 0;
	      maxcon++;
	    }
	    data->nodalgraph[l][ind] = ind2;
	    totcon++;
	  }
	}
      }
    }

  }

  /* This adds the periodic connections */
  if( data->periodicexist ) {
    for(ind=1;ind<=noknots;ind++) {
      ind2 = data->periodic[ind];      
      if(ind == ind2) continue;

      hit = FALSE;
      for(l=0;l<maxcon;l++) { 
	if(data->nodalgraph[l][ind] == ind2) hit = TRUE;
	if(data->nodalgraph[l][ind] == 0) break;
      }
      if(!hit) {
	if(l >= maxcon) {
	  data->nodalgraph[maxcon] = Ivector(1,noknots);
	  for(m=1;m<=noknots;m++)
	    data->nodalgraph[maxcon][m] = 0;
	  maxcon++;
	}
	data->nodalgraph[l][ind] = ind2;
	totcon++;
	percon++;
      }
    }
  }

  data->nodalmaxconnections = maxcon;
  data->nodalexists = TRUE;
  
  if(info) {
    printf("There are at maximum %d connections in nodal graph.\n",maxcon);
    printf("There are at all in all %d connections in nodal graph.\n",totcon);
    if(percon) printf("There are %d periodic connections in nodal graph.\n",percon);
  }

  return(0);
}


int DestroyNodalGraph(struct FemType *data,int info)
{
  int i,maxcon, noknots;
  
  if(!data->nodalexists) {
    printf("You tried to destroy a non-existing nodal graph\n");
    return(1);
  }

  maxcon = data->nodalmaxconnections;
  noknots = data->noknots;

  for(i=0;i<maxcon;i++)    
    free_Ivector(data->nodalgraph[i],1,noknots);

  data->nodalmaxconnections = 0;
  data->nodalexists = FALSE; 

  if(info) printf("The nodal graph was destroyed\n");
  return(0);
}



int CreateInverseTopology(struct FemType *data,int info)
{
  int i,j,k,l,m,noelements,noknots,elemtype,nonodes,ind;
  int *neededby,minneeded,maxneeded;
  int step,totcon;
  int *rows,*cols;
  struct CRSType *invtopo;

  invtopo = &data->invtopo;
  if(invtopo->created) {
    if(0) printf("The inverse topology already exists!\n");
    return(0);
  }

  printf("Creating an inverse topology of the finite element mesh\n");  

  noelements = data->noelements;
  noknots = data->noknots;

  neededby = Ivector(1,noknots);
  totcon = 0;

  for(step=1;step<=2;step++) {

    for(i=1;i<=noknots;i++)
      neededby[i] = 0;

    for(i=1;i<=noelements;i++) {
      elemtype = data->elementtypes[i];
      nonodes = data->elementtypes[i] % 100;
      
      for(j=0;j<nonodes;j++) {
	ind = data->topology[i][j];

	if( step == 1 ) {
	  neededby[ind] += 1;
	  totcon += 1;
	}
	else {
	  k = rows[ind-1] + neededby[ind];
	  cols[k] = i-1;
	  neededby[ind] += 1;
	}
      }
    }

    if( step == 1 ) {
      rows = Ivector( 0, noknots );
      rows[0] = 0;
      for(i=1;i<=noknots;i++) 
	rows[i] = rows[i-1] + neededby[i];

      cols = Ivector( 0, totcon-1 );
      for(i=0;i<totcon;i++)
	cols[i] = 0;

      invtopo->cols = cols;
      invtopo->rows = rows;
      invtopo->colsize = totcon;
      invtopo->rowsize = noknots;
      invtopo->created = TRUE;
    }
  }    

  minneeded = maxneeded = neededby[1];
  for(i=1;i<=noknots;i++) {
    minneeded = MIN( minneeded, neededby[i]);
    maxneeded = MAX( maxneeded, neededby[i]);
  }

  free_Ivector(neededby,1,noknots);

  if(info) {
    printf("There are from %d to %d connections in the inverse topology.\n",minneeded,maxneeded);
    printf("Each node is in average in %.3f elements\n",1.0*totcon/noknots);
  }

  return(0);
}



int CreateDualGraph(struct FemType *data,int unconnected,int info)
{
  int totcon,dcon,noelements,noknots,elemtype,nonodes,i,j,k,l,i2,m,ind,hit,ci,ci2;
  int dualmaxcon,invmaxcon,showgraph,freeelements,step,orphanelements,stat;
  int *elemconnect,*neededby;
  int *dualrow,*dualcol,dualsize,dualmaxelem,allocated;
  int *invrow,*invcol;
  struct CRSType *dualgraph;

  printf("Creating a dual graph of the finite element mesh\n");  

  dualgraph = &data->dualgraph;
  if(dualgraph->created) {
    printf("The dual graph already exists!\n");
    return(1);
  }

  CreateInverseTopology(data,info);

  noelements = data->noelements;
  noknots = data->noknots;
  freeelements = noelements;
  orphanelements = 0;
  
  /* If a dual graph only for the unconnected nodes is requested do that.
     Basically the connected nodes are omitted in the graph. */
  if( unconnected ) {
    printf("Removing connected nodes from the dual graph\n");
    if( data->nodeconnectexist ) {
      if(info) printf("Creating connected elements list from the connected nodes\n");
      SetConnectedElements(data,info);
    }
    if( data->elemconnectexist ) {
      elemconnect = data->elemconnect;
      freeelements -= data->elemconnectexist;
    }
    else {
      unconnected = FALSE;
    }
    if(info) printf("List of unconnected elements created\n");
  }

  showgraph = FALSE;
  if(showgraph) printf("elemental graph ij pairs\n");

  data->dualexists = TRUE;
  dualmaxcon = 0;
  dualmaxelem = 0;
 
  invrow = data->invtopo.rows;
  invcol = data->invtopo.cols;


  /* This marker is used to identify the connections already accounted for */  
  neededby = Ivector(1,freeelements);
  for(i=1;i<=freeelements;i++)
    neededby[i] = 0;

  allocated = FALSE;
 omstart: 

  totcon = 0;

  for(i=1;i<=noelements;i++) {
    if(showgraph) printf("%d :: ",i);

    dcon = 0;
    elemtype = data->elementtypes[i];    
    nonodes = data->elementtypes[i] % 100;

    if( unconnected ) {
      ci = elemconnect[i];
      if( ci < 0 ) continue;
    }
    else {
      ci = i;
    }      
    if(allocated) dualrow[ci-1] = totcon;

    if(0) printf("i=%d %d\n",i,elemtype);

    for(step=1;step<=2;step++) {
      for(j=0;j<nonodes;j++) {
	ind = data->topology[i][j];
	
	if(0) printf("ind=%d\n",ind);


	for(k=invrow[ind-1];k<invrow[ind];k++) {
	  i2 = invcol[k]+1;
	  
	  if( i2 == i ) continue;
	  
	  if( unconnected ) {
	    ci2 = elemconnect[i2];
	    if( ci2 < 0 ) continue;
	  }
	  else {
	    ci2 = i2;
	  }

	  /* In the first cycle mark the needed connections,
	     and in the second cycle set the marker to zero for next round. */
	  if( step == 1 ) {
	    if( neededby[ci2] ) continue;
	    neededby[ci2] = TRUE;

	    if(0) printf("ci ci2 = %d %d\n",ci,ci2);

	    /* If the dual graph has been allocated populate it */
	    if(allocated) {	      
	      dualcol[totcon] = ci2-1;
	    }

	    dcon += 1;
	    totcon += 1;
	    if( dcon > dualmaxcon ) {
	      dualmaxcon = dcon;
	      dualmaxelem = i;
	    }
	  }
	  else {
	    neededby[ci2] = FALSE;
	  }
	}
      }
    }
    if( dcon == 0 && allocated ) {
      orphanelements += 1;
    }
  }	

  if(allocated) {
    dualrow[dualsize] = totcon;
  }
  else {
    dualsize = freeelements;
    dualrow = Ivector(0,dualsize);
    for(i=1;i<=dualsize;i++) 
      dualrow[i] = 0;

    dualcol = Ivector(0,totcon-1);
    for(i=0;i<totcon;i++) 
      dualcol[i] = 0;

    dualgraph->cols = dualcol;
    dualgraph->rows = dualrow;
    dualgraph->rowsize = dualsize;
    dualgraph->colsize = totcon;
    dualgraph->created = TRUE;

    allocated = TRUE;

    goto omstart;
  } 

  if( orphanelements ) {
    printf("There are %d elements in the dual mesh that are not connected!\n",orphanelements);
    if(unconnected) printf("The orphan elements are likely caused by the hybrid partitioning\n");
  }


#if 0
  j = totcon; k = 0;
  for(i=1;i<=dualsize;i++){
    l = dualrow[i]-dualrow[i-1];
    if(l <= 0 ) printf("row(%d) = %d %d %d\n",i,l,dualrow[i],dualrow[i-1]);
    j = MIN(j,l);
    k = MAX(k,l);
  }
  printf("range dualrow: %d %d\n",j,k);

  j = totcon; k = 0;
  for(i=0;i<totcon;i++) {
    j = MIN(j,dualcol[i]);
    k = MAX(k,dualcol[i]);
  }
  printf("range dualcol: %d %d\n",j,k);
#endif


  if(info) {
    printf("There are at maximum %d connections in dual graph (in element %d).\n",dualmaxcon,dualmaxelem);
    printf("There are at all in all %d connections in dual graph.\n",totcon);
    printf("Average connection per active element in dual graph is %.3f\n",1.0*totcon/freeelements);
  }  

  free_Ivector( neededby,1,freeelements);

  /* Inverse topology is created for partitioning only and then the direct
     topology is needed elsewhere as well. Do do not destroy it. */ 
  if(0) stat = DestroyInverseTopology(data,info);

  return(0);
}


int DestroyCRSMatrix(struct CRSType *sp) {

  if(sp->created) {
    if(0) printf("You tried to destroy a non-existing sparse matrix\n");
    return(1);
  }

  free_Ivector( sp->rows, 0, sp->rowsize );
  free_Ivector( sp->cols, 0, sp->colsize-1);
  sp->rowsize = 0;
  sp->colsize = 0;
  sp->created = FALSE;

  return(0);
}


int DestroyInverseTopology(struct FemType *data,int info)
{
  int stat;
  stat = DestroyCRSMatrix( &data->invtopo );
  return(stat);
}

int DestroyDualGraph(struct FemType *data,int info)
{
  int stat;
  stat = DestroyCRSMatrix( &data->dualgraph );
  return(stat);
}




int MeshTypeStatistics(struct FemType *data,int info)
{
  int i,elemtype,maxelemtype,minelemtype;
  int *elemtypes=NULL;

  maxelemtype = minelemtype = data->elementtypes[1];

  for(i=1;i<=data->noelements;i++) {
    elemtype = data->elementtypes[i];
    maxelemtype = MAX( maxelemtype, elemtype );
    minelemtype = MIN( minelemtype, elemtype );
  }

  elemtypes = Ivector(minelemtype,maxelemtype);
  for(i=minelemtype;i<=maxelemtype;i++)
    elemtypes[i] = 0;

  for(i=1;i<=data->noelements;i++) {
    elemtype = data->elementtypes[i];
    elemtypes[elemtype] += 1;
  }

  if(info) {
    printf("Number of different elementtypes\n");
    for(i=minelemtype;i<=maxelemtype;i++)
      if(elemtypes[i]) printf("\t%d\t%d\n",i,elemtypes[i]);
  }

  free_Ivector(elemtypes,minelemtype,maxelemtype);
  return(0);
}

