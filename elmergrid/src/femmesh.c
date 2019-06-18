/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter R�back
   Email: Peter.Raback@csc.fi
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


 /* -----------------------:  femmesh.c  :----------------------

   These subroutines are used to create a mesh. They operate on 
   variables that are of type MeshType and GridType. The routines 
   are used to initialize the grid used for the calculation. 
   */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "nrutil.h"
#include "common.h"
#include "femdef.h"
#include "femtypes.h"
#include "femmesh.h"

#define DEBUG 0


void InitGrid(struct GridType *grid)
/* Initializes the grid of a specific mesh. A grid can differ 
   between various differential equations. 
   */
{
  int i,j;

  grid->layered = FALSE;
  grid->layeredbc = TRUE;
  grid->layerbcoffset = 0;
  grid->triangles = FALSE;
  grid->triangleangle = 0.0;
  grid->partitions = FALSE;
  grid->wantedelems = 0;
  grid->wantedelems3d = 0;
  grid->wantednodes3d = 0;
  grid->limitdx = 0.0;
  grid->limitdxverify = FALSE;

  grid->dimension = 2;
  grid->xcells = grid->ycells = grid->zcells = 0;
  grid->nocells = 0;
  grid->noknots = grid->noelements = 0;
  grid->totxelems = grid->totyelems = grid->totzelems = 0;
  grid->minxelems = grid->minyelems = 1;
  grid->minzelems = 2;
  grid->firstmaterial = 1;
  grid->lastmaterial = INT_MAX;
  grid->autoratio = 1;
  grid->xyratio = 1.0;
  grid->xzratio = 1.0;
  grid->nonodes = 4;
  grid->numbering = NUMBER_XY;
  grid->rotate = FALSE;
  grid->rotateradius1 = 0.0;
  grid->rotateradius2 = 0.0;
  grid->rotateimprove = 1.0;
  grid->rotateblocks = 4;
  grid->rotatecartesian = FALSE;
  grid->reduceordermatmin = 0;
  grid->reduceordermatmax = 0;
  grid->polarradius = 1.0;
  grid->elemorder = 1;
  grid->elemmidpoints = FALSE;

  grid->rotatecurve = FALSE;
  grid->curverad = 0.5;
  grid->curveangle = 90.0;
  grid->curvezet = 1.0;

  for(i=0;i<=MAXCELLS;i++) {
    grid->xelems[i] = 0;
    grid->xlinear[i] = TRUE;
    grid->xratios[i] = 1.;
    grid->xexpand[i] = 1.;
    grid->xdens[i] = 1.;
    grid->x[i] = 0.; 
  }
  for(i=0;i<=MAXCELLS;i++) {
    grid->yelems[i] = 0;
    grid->ylinear[i] = TRUE;
    grid->yratios[i] = 1.0;
    grid->yexpand[i] = 1.0;
    grid->ydens[i] = 1.0;
    grid->y[i] = 0.;
  }
  for(i=0;i<=MAXCELLS;i++) {
    grid->zelems[i] = 0;
    grid->zlinear[i] = TRUE;
    grid->zratios[i] = 1.0;
    grid->zexpand[i] = 1.0;
    grid->zdens[i] = 1.0;
    grid->z[i] = 0.;
    grid->zfirstmaterial[i] = 1;
    grid->zlastmaterial[i] = INT_MAX;
    grid->zmaterial[i] = 0; 
  }

  grid->zmaterialmapexists = FALSE;
  grid->zhelicityexists = FALSE;
  grid->zhelicity = 0.0;

  /* Initilizes the numbering of the cells. */
  for(j=0;j<=MAXCELLS+1;j++)
    for(i=0;i<=MAXCELLS+1;i++) 
      grid->structure[j][i] = MAT_NOTHING; 

  for(j=0;j<=MAXCELLS+1;j++)
    for(i=0;i<=MAXCELLS+1;i++) 
      grid->numbered[j][i] = 0; 

  grid->noboundaries = 0;
  for(i=0;i<MAXBOUNDARIES;i++) {
    grid->boundint[i] = 0;
    grid->boundext[i] = 0;
    grid->boundsolid[i] = 0;
    grid->boundtype[i] = 0;
  }

  grid->mappings = 0;
  for(i=0;i<MAXMAPPINGS;i++) {
    grid->mappingtype[i] = 0;
    grid->mappingline[i] = 0;    
    grid->mappinglimits[2*i] = 0;    
    grid->mappinglimits[2*i+1] = 0;    
    grid->mappingpoints[i] = 0;
  }
}


void ExampleGrid1D(struct GridType **grids,int *nogrids,int info) 
/* Creates an example grid that might be used to analyze 
   flow through a step. */
{
  int j;
  struct GridType *grid;

  (*nogrids) = 1;
  (*grids) = (struct GridType*)
    malloc((size_t) (*nogrids)*sizeof(struct GridType)); 

  grid = grids[0];

  InitGrid(grid);

  grid->nonodes = 2;
  grid->xcells = 3;
  grid->ycells = 1;
  grid->wantedelems = 40;
  grid->firstmaterial = 1;
  grid->lastmaterial = 1;
  grid->autoratio = 1;
  grid->coordsystem = COORD_CART1;
  grid->numbering = NUMBER_1D;

  grid->x[0] = 0.0;
  grid->x[1] = 1.0;
  grid->x[2] = 3.0;
  grid->x[3] = 4.0;

  grid->xexpand[1] = 2.0;
  grid->xexpand[3] = 0.5;
  grid->xdens[2] = 0.5;

  for(j=1;j<=3;j++)
    grid->structure[1][j] = 1;

  grid->noboundaries = 2;
  grid->boundint[0] = 1;
  grid->boundint[1] = 1;

  grid->boundext[0] = -1;
  grid->boundext[1] = -2;

  grid->boundtype[0] = 1;
  grid->boundtype[1] = 2;

  grid->boundsolid[0] = 1;
  grid->boundsolid[1] = 1;

  grid->mappings = 0;
  
  if(info) printf("A simple 1D example mesh was created\n");
}


void ExampleGrid2D(struct GridType **grids,int *nogrids,int info) 
/* Creates an example grid that might be used to analyze 
   flow through a step. */
{
  int j;
  struct GridType *grid;

  (*nogrids) = 1;
  (*grids) = (struct GridType*)
    malloc((size_t) (*nogrids)*sizeof(struct GridType)); 

  grid = grids[0];

  InitGrid(grid);

  grid->xcells = 4;
  grid->ycells = 3;
  grid->wantedelems = 100;
  grid->totxelems = grid->totyelems = grid->totzelems = 10;
  grid->firstmaterial = 1;
  grid->lastmaterial = 1;
  grid->autoratio = 1;
  grid->coordsystem = COORD_CART2;
  grid->numbering = NUMBER_XY;

  grid->x[0] = -1.0;
  grid->x[1] = 0.0;
  grid->x[2] = 2.0;
  grid->x[3] = 6.0;
  grid->x[4] = 7.0;

  grid->y[0] = 0.0;
  grid->y[1] = 0.5;
  grid->y[2] = 0.75;
  grid->y[3] = 1.0;

  grid->xexpand[2] = 0.4;
  grid->xexpand[3] = 5.0;

  grid->yexpand[2] = 2.0;
  grid->yexpand[3] = 0.5;

  for(j=1;j<=3;j++)
    grid->structure[j][1] = 2;
  for(j=1;j<=3;j++)
    grid->structure[j][2] = 1;
  grid->structure[1][2] = 0;
  for(j=1;j<=3;j++)
    grid->structure[j][3] = 1;
  for(j=1;j<=3;j++)
    grid->structure[j][4] = 3;

  grid->noboundaries = 3;
  grid->boundint[0] = 1;
  grid->boundint[1] = 1;
  grid->boundint[2] = 1;

  grid->boundext[0] = 2;
  grid->boundext[1] = 3;
  grid->boundext[2] = 0;

  grid->boundtype[0] = 1;
  grid->boundtype[1] = 2;
  grid->boundtype[2] = 3;

  grid->boundsolid[0] = 1;
  grid->boundsolid[1] = 1;
  grid->boundsolid[2] = 1;

  grid->mappings = 2;
  grid->mappingtype[0] = 1;
  grid->mappingtype[1] = 1;
  grid->mappingline[0] = 0;
  grid->mappingline[1] = 3;
  grid->mappinglimits[0] = grid->mappinglimits[1] = 0.5;
  grid->mappinglimits[2] = grid->mappinglimits[3] = 0.5;
  grid->mappingpoints[0] = grid->mappingpoints[1] = 4;
  grid->mappingparams[0] = Rvector(0,grid->mappingpoints[0]);
  grid->mappingparams[1] = Rvector(0,grid->mappingpoints[1]);
  grid->mappingparams[0][0] = 2.0;
  grid->mappingparams[0][1] = 0.0;
  grid->mappingparams[0][2] = 6.0;
  grid->mappingparams[0][3] = -0.2;
  grid->mappingparams[1][0] = 0.0;
  grid->mappingparams[1][1] = 0.0;
  grid->mappingparams[1][2] = 6.0;
  grid->mappingparams[1][3] = 0.3;
  
  if(info) printf("A simple 2D example mesh was created\n");
}

  
void ExampleGrid3D(struct GridType **grids,int *nogrids,int info) 
/* Creates an example grid that might be used to analyze 
   a simple accelaration sensor. */
{
  int i,j;
  struct GridType *grid;

  (*nogrids) = 1;
  (*grids) = (struct GridType*)
    malloc((size_t) (*nogrids)*sizeof(struct GridType)); 

  grid = grids[0];

  InitGrid(grid);

  grid->dimension = 3;
  grid->coordsystem = COORD_CART3;
  grid->layered = TRUE;

  grid->xcells = 4;
  grid->ycells = 5;
  grid->zcells = 5;

  grid->wantedelems = 1000;
  grid->firstmaterial = 1;
  grid->lastmaterial = 3;
  grid->autoratio = 1;
  grid->coordsystem = COORD_CART3;
  grid->numbering = NUMBER_XY;

  grid->x[0] = 0.0;
  grid->x[1] = 0.1;
  grid->x[2] = 0.4;
  grid->x[3] = 0.5;
  grid->x[4] = 0.6;

  grid->y[0] = 0.0;
  grid->y[1] = 0.1;
  grid->y[2] = 0.2;
  grid->y[3] = 0.4;
  grid->y[4] = 0.5;
  grid->y[5] = 0.6;

  grid->z[0] = 0.0;
  grid->z[1] = 0.1;
  grid->z[2] = 0.2;
  grid->z[3] = 0.4;
  grid->z[4] = 0.5;
  grid->z[5] = 0.6;

  for(j=1;j<=5;j++)
    for(i=1;i<=4;i++)
      grid->structure[j][i] = 1;

  grid->structure[2][2] = 3;
  grid->structure[2][3] = 3;
  grid->structure[3][3] = 3;
  grid->structure[4][3] = 3;
  grid->structure[4][2] = 3;

  grid->structure[3][2] = 2;

  grid->zlastmaterial[5] = 3;
  grid->zlastmaterial[2] = 1;
  grid->zlastmaterial[3] = 2;
  grid->zlastmaterial[4] = 1;
  grid->zlastmaterial[5] = 3;

  grid->zmaterial[1] = 1;
  grid->zmaterial[2] = 2;
  grid->zmaterial[3] = 3;
  grid->zmaterial[4] = 4;
  grid->zmaterial[5] = 5;

  grid->noboundaries = 4;
  grid->boundint[0] = 1;
  grid->boundint[1] = 1;
  grid->boundint[2] = 1;
  grid->boundint[2] = 2;

  grid->boundext[0] = 0;
  grid->boundext[1] = 2;
  grid->boundext[2] = 3;
  grid->boundext[3] = 3;

  grid->boundtype[0] = 1;
  grid->boundtype[1] = 2;
  grid->boundtype[2] = 3;
  grid->boundtype[3] = 4;

  grid->boundsolid[0] = 1;
  grid->boundsolid[1] = 1;
  grid->boundsolid[2] = 1;
  grid->boundsolid[3] = 1;

  if(info) printf("A simple 3D example mesh was created\n");
}


void SetElementDivision(struct GridType *grid,Real relh,int info)
/* Given the densities and ratios in each cell finds the 
   optimum way to devide the mesh into elements. 
   The procedure is the following:
   For each subcell set the minimum number of elements 
   then add one element at a time till the number of 
   elements is the desired one. The added element
   is always set to the subcell having the sparsest mesh.
   Also numbers the cells taking into consideration only 
   materials that have indeces in interval [firstmat,lastmat]. 
   */
{
  int i,j,nx,ny,nxmax = 0,nymax = 0;
  int sumxelems,sumyelems,sumxyelems;
  int wantedelems,wantedelemsx,wantedelemsy;
  Real ratio,linearlimit;
  Real dxmax = 0,dymax = 0,dx = 0,dy = 0,dxlimit = 0;

  if(0) printf("SetElementDivision\n");

  if(grid->dimension == 1) 
    grid->numbering = NUMBER_1D;

  nx = grid->xcells;
  ny = grid->ycells;
  linearlimit = 0.001;

  /* Lets number the cells from left to right and up to down. */
  grid->nocells = 0;
  for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++) {
      if (grid->structure[j][i] >= grid->firstmaterial  &&  
	  grid->structure[j][i] <= grid->lastmaterial)   
        grid->numbered[j][i] = ++grid->nocells; 
    }
  if(0) printf("The mesh is devided into %d separate subcells.\n",grid->nocells);
  
  /* Put the linearity flags. */
  for(i=1; i<= nx ;i++) {
    if (fabs(1.-grid->xexpand[i]) < linearlimit) 
      grid->xlinear[i] = TRUE; 
    else if (fabs(1.+grid->xexpand[i]) < linearlimit) 
      grid->xlinear[i] = TRUE; 
    else 
      grid->xlinear[i] = FALSE;
  }

  for(i=1; i <= ny ;i++) {
    if(fabs(1.-grid->yexpand[i]) < linearlimit) 
      grid->ylinear[i] = TRUE;
    else if(fabs(1.+grid->yexpand[i]) < linearlimit) 
      grid->ylinear[i] = TRUE;
    else 
      grid->ylinear[i] = FALSE;
  }


 /* If there are no materials no elements either.
     Otherwise at least grid->minelems elments 
     If you want to number elements even if they are 
     not there change this initialization to minelems. */
  if(grid->autoratio) {
    for(i=1;i<=nx;i++)  
      grid->xelems[i] = grid->minxelems;
    for(i=1;i<=ny;i++)  
      grid->yelems[i] = grid->minyelems;
  }    
  else {
    for(i=1;i<=nx;i++)  
      if(grid->xelems[i] < grid->minxelems) grid->xelems[i] = grid->minxelems;
    for(i=1;i<=ny;i++)  
      if(grid->yelems[i] < grid->minyelems) grid->yelems[i] = grid->minyelems;
  }

  sumxelems = 0;
  for(i=1;i<=nx;i++) {
    if(grid->dimension == 1 && !grid->numbered[1][i]) continue;
    sumxelems += grid->xelems[i];
  }
  sumyelems = 0;
  for(i=1;i<=ny;i++)  
    sumyelems += grid->yelems[i];

  if(grid->dimension == 1) {
    grid->autoratio = 2;
    grid->totxelems = grid->wantedelems;
    grid->totyelems = 1;
  }

  /* Allocate elements for both axis separately */
  if(grid->autoratio == 2) {

    wantedelemsx = (int) (1.0*grid->totxelems / relh);
    wantedelemsy = (int) (1.0*grid->totyelems / relh);

    if(sumxelems > wantedelemsx) {
      printf("SetElementDivision: %d is too few elements in x-direction\n",grid->totxelems);
      wantedelemsx = sumxelems+1;
    }      
    if(sumyelems > wantedelemsy ) {
      printf("SetElementDivision: %d is too few elements in y-direction\n",grid->totyelems);
      wantedelemsy = sumyelems+1;
    }   
    
    for(;sumxelems < wantedelemsx;) {
      dxmax = 0.0;
      for(i=1;i<=nx;i++) {
	if(grid->xelems[i] == 0) continue;
	
	/* Don't put elements to passive subcells */
	if(grid->dimension == 1 && !grid->numbered[1][i]) continue; 
	
	if(grid->xlinear[i] == TRUE || grid->xelems[i]==1) 
	  dx = (grid->x[i] - grid->x[i-1])/grid->xelems[i];
	else {
	  if(grid->xexpand[i] > 0.0) {
	    ratio = pow(grid->xexpand[i],1./(grid->xelems[i]-1.));
	    dx = (grid->x[i] - grid->x[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i])));
	    if(ratio < 1.)   
	      dx *= grid->xexpand[i];
	  }
	  else {	    
	    if(grid->xelems[i]==2) {
	      dx = (grid->x[i] - grid->x[i-1])/grid->xelems[i];
	    } 
	    else if(grid->xelems[i]%2 == 0) {
	      ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2-1.));
	      dx = 0.5 * (grid->x[i] - grid->x[i-1]) * 
		(1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i]/2)));
	    }
	    else if(grid->xelems[i]%2 == 1) {
	      ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2));
	      dx = (grid->x[i] - grid->x[i-1]) / 
		(2.0*(1.-pow(ratio,(Real)(grid->xelems[i]/2)))/
		 (1-ratio) + pow(ratio,(Real)(grid->xelems[i]/2+0.5)));
	    }
	  }
	}
	dx *= grid->xdens[i];
	if(dx > dxmax) {
	  dxmax = dx;
	  nxmax = i;
	}
      }      
      grid->xelems[nxmax] += 1;
      sumxelems++;
    }
  

    if(grid->dimension > 1) {
      for(;sumyelems < wantedelemsy;) {
	dymax = 0.0;
	for(i=1;i<=ny;i++) {
	  if(grid->yelems[i] == 0) continue;
	  if(grid->ylinear[i] || grid->yelems[i]==1) 
	    dy = (grid->y[i] - grid->y[i-1])/grid->yelems[i];
	  else {
	    ratio = pow(grid->yexpand[i],1./(grid->yelems[i]-1));
	    dy = (grid->y[i] - grid->y[i-1]) * 
	      (1.-ratio)/(1.-pow(ratio,(Real)(grid->yelems[i])));
	    if(ratio < 1.)   
	      dy *= grid->yexpand[i];
	  }
	  dy *= grid->ydens[i];
	  if(dy > dymax) {
	    dymax = dy;
	    nymax = i;
	  }
	}      
	grid->yelems[nymax] += 1;
	sumyelems++;
      }
    }
  }   


  /* Both axis dependently */
  if(grid->autoratio == 1)  {
    
    wantedelems = (int) (1.0*grid->wantedelems / (relh*relh));

    sumxyelems = 0;
    for(;;) {
      dxmax = 0.0;

      for(i=1;i<=nx;i++) {
	
	if(grid->xelems[i] == 0) continue;
	if(grid->xlinear[i] || grid->xelems[i]==1) 
	  dx = (grid->x[i] - grid->x[i-1])/grid->xelems[i];
	else {
	  if(grid->xexpand[i] > 0.0) {
	    ratio = pow(grid->xexpand[i],1./(grid->xelems[i]-1.));
	    dx = (grid->x[i] - grid->x[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i])));
	    if(ratio < 1.)   
	      dx *= grid->xexpand[i];
	  }
	  else if(grid->xelems[i]==2) {
	    dx = (grid->x[i] - grid->x[i-1])/grid->xelems[i];
	  } 
	  else if(grid->xelems[i]%2 == 0) {
	    ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2-1.));
	    dx = 0.5 * (grid->x[i] - grid->x[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i]/2)));
	  }
	  else if(grid->xelems[i]%2 == 1) {
	    ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2));
	    dx = (grid->x[i] - grid->x[i-1]) / 
	      (2.0*(1.-pow(ratio,(Real)(grid->xelems[i]/2)))/
	       (1-ratio) + pow(ratio,(Real)(grid->xelems[i]/2+0.5)));
	  }
	}

	dx *= grid->xdens[i];
	if(dx > dxmax) {
	  dxmax = dx;
	  nxmax = i;
	}
      }      


      dymax = 0.0;
      for(i=1;i<=ny;i++) {

	if(grid->yelems[i] == 0) continue;
	if(grid->ylinear[i] || grid->yelems[i]==1) 
	  dy = (grid->y[i] - grid->y[i-1])/grid->yelems[i];
	else {
	  if(grid->yexpand[i] > 0.0) {
	    ratio = pow(grid->yexpand[i],1./(grid->yelems[i]-1));
	    dy = (grid->y[i] - grid->y[i-1]) * 
	      (1.-ratio)/(1.-pow(ratio,(Real)(grid->yelems[i])));
	    if(ratio < 1.)   
	      dy *= grid->yexpand[i];
	  }
	  else if(grid->yelems[i]==2) {
	    dy = (grid->y[i] - grid->y[i-1])/grid->yelems[i];
	  } 
	  else if(grid->yelems[i]%2 == 0) {
	    ratio = pow(-grid->yexpand[i],1./(grid->yelems[i]/2-1.));
	    dy = 0.5 * (grid->y[i] - grid->y[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->yelems[i]/2)));
	  }
	  else if(grid->yelems[i]%2 == 1) {
	    ratio = pow(-grid->yexpand[i],1./(grid->yelems[i]/2));
	    dy = (grid->y[i] - grid->y[i-1]) / 
	      (2.0*(1.-pow(ratio,(Real)(grid->yelems[i]/2)))/
	       (1-ratio) + pow(ratio,(Real)(grid->yelems[i]/2+0.5)));
	  }
	}
	dy *= grid->ydens[i];
	if(dy > dymax) {
	  dymax = dy;
	  nymax = i;
	}
      }      
      dymax /= grid->xyratio;

      if(dxmax > dymax) {
	grid->xelems[nxmax] += 1;
	sumxelems++;
	dxlimit = dxmax;
      }
      else {
	grid->yelems[nymax] += 1;
	sumyelems++;
	dxlimit = dymax;
      }
      
      sumxyelems = 0;
      for(j=1;j<=grid->ycells;j++)
	for(i=1;i<=grid->xcells;i++) 
	  if(grid->numbered[j][i]) sumxyelems += grid->xelems[i] * grid->yelems[j];

      if(wantedelems <= sumxyelems) break;
    }
  }    



  if(grid->autoratio == 3 || grid->limitdxverify)  {
    
    if(grid->autoratio == 3) {
      dxlimit = relh * grid->limitdx;
      dxmax = dymax = dxlimit;
    }

    for(i=1;i<=nx;i++) {

      for(;;) {       
	if(grid->xlinear[i] || grid->xelems[i]==0) 
	  dx = (grid->x[i] - grid->x[i-1])/(grid->xelems[i]+1);
	else {
	  if(grid->xexpand[i] > 0.0) {
	    ratio = pow(grid->xexpand[i],1./grid->xelems[i]);
	    dx = (grid->x[i] - grid->x[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i]+1)));
	    if(ratio < 1.)   
	      dx *= grid->xexpand[i];
	  }
	  else if(grid->xelems[i]==1) {
	    dx = (grid->x[i] - grid->x[i-1])/(grid->xelems[i]+1);
	  } 
	  else if((grid->xelems[i]+1)%2 == 0) {
	    ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2));
	    dx = 0.5 * (grid->x[i] - grid->x[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i]/2+1)));
	  }
	  else if((grid->xelems[i]+1)%2 == 1) {
	    ratio = pow(-grid->xexpand[i],1./((grid->xelems[i]+1)/2));
	    dx = (grid->x[i] - grid->x[i-1]) / 
	      (2.0*(1.-pow(ratio,(Real)((grid->xelems[i]+1)/2)))/
	       (1-ratio) + pow(ratio,(Real)((grid->xelems[i]+1)/2+0.5)));
	  }
	}
	dx *= grid->xdens[i];
	
	/* choose the best fit for desired density */
	if(fabs(dx-dxlimit) < fabs( dx*(1+1.0/grid->xelems[i]) -dxlimit) )
	  grid->xelems[i] += 1;
	else
	  break;
      }
      sumxelems += grid->xelems[i];
    }      


    for(i=1;i<=ny;i++) {
      
      for(;;) {
	if(grid->ylinear[i] || grid->yelems[i]==0) 
	  dy = (grid->y[i] - grid->y[i-1])/(grid->yelems[i]+1);
	else {
	  if(grid->yexpand[i] > 0.0) {
	    ratio = pow(grid->yexpand[i],1./grid->yelems[i]);
	    dy = (grid->y[i] - grid->y[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->yelems[i]+1)));
	    if(ratio < 1.)   
	      dy *= grid->yexpand[i];
	  }
	  else if(grid->yelems[i]==1) {
	    dy = (grid->y[i] - grid->y[i-1])/(grid->yelems[i]+1);
	  } 
	  else if((grid->yelems[i]+1)%2 == 0) {
	    ratio = pow(-grid->yexpand[i],1./(grid->yelems[i]/2));
	    dy = 0.5 * (grid->y[i] - grid->y[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->yelems[i]/2+1)));
	  }
	  else if((grid->yelems[i]+1)%2 == 1) {
	    ratio = pow(-grid->yexpand[i],1./((grid->yelems[i]+1)/2));
	    dy = (grid->y[i] - grid->y[i-1]) / 
	      (2.0*(1.-pow(ratio,(Real)((grid->yelems[i]+1)/2)))/
	       (1-ratio) + pow(ratio,(Real)((grid->yelems[i]+1)/2+0.5)));
	  }
	}
	
	dy *= grid->ydens[i] / grid->xyratio;
	/* choose the best fit for desired density */
	if(fabs(dy-dxlimit) < fabs( dy*(1+1.0/grid->yelems[i]) -dxlimit) )
	  grid->yelems[i] += 1;
	else
	  break; 
      }
      sumyelems += grid->yelems[i];
    }      
    
    sumxyelems = 0;
    for(j=1;j<=grid->ycells;j++)
      for(i=1;i<=grid->xcells;i++) 
	if(grid->numbered[j][i]) sumxyelems += grid->xelems[i] * grid->yelems[j];
  }


  /* Put the linearity flags if there is only one division. */
  for(i=1; i<= nx ;i++) 
    if (grid->xelems[i] == 1) 
      grid->xlinear[i] = TRUE; 
  for(i=1; i <= ny ;i++) 
    if(grid->yelems[i] == 1) 
      grid->ylinear[i] = TRUE;
 

  /* Set the size of the first element within each subcell */
  for(i=1;i<=nx;i++) {

    if(grid->xlinear[i] == TRUE) 
      grid->dx[i] = (grid->x[i] - grid->x[i-1])/grid->xelems[i];
    else {
      if(grid->xexpand[i] > 0.0) {
	ratio = pow(grid->xexpand[i],1./(grid->xelems[i]-1.));
	grid->xratios[i] = ratio;
	grid->dx[i] = (grid->x[i] - grid->x[i-1]) * 
	  (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i])));
      }
      else if(grid->xelems[i] == 2) {
	grid->dx[i] = (grid->x[i] - grid->x[i-1])/grid->xelems[i];
      }
      else if(grid->xelems[i]%2 == 0) {
	ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2-1.));
	grid->xratios[i] = ratio;
	grid->dx[i] = 0.5 * (grid->x[i] - grid->x[i-1]) * 
	  (1.-ratio) / (1.-pow(ratio,(Real)(grid->xelems[i]/2)));
      }
      else if(grid->xelems[i]%2 == 1) {
	ratio = pow(-grid->xexpand[i],1./(grid->xelems[i]/2));
	grid->xratios[i] = ratio;
	grid->dx[i] = (grid->x[i] - grid->x[i-1]) / 
	  (2.0*(1.-pow(ratio,(Real)(grid->xelems[i]/2)))/
	   (1-ratio) + pow(ratio,(Real)(grid->xelems[i]/2+0.5)));
      }
    }	
  }

  for(i=1;i<=ny;i++) {

    if(grid->ylinear[i] == TRUE) 
      grid->dy[i] = (grid->y[i] - grid->y[i-1])/grid->yelems[i];
    else {
      if(grid->yexpand[i] > 0.0) {
      ratio = pow(grid->yexpand[i],1./(grid->yelems[i]-1.));
      grid->yratios[i] = ratio;
      grid->dy[i] = (grid->y[i] - grid->y[i-1]) * 
	(1.-ratio)/(1.-pow(ratio,(Real)(grid->yelems[i])));
      }
      else if(grid->yelems[i] == 2) {
	grid->dy[i] = (grid->y[i] - grid->y[i-1])/grid->yelems[i];
      }
      else if(grid->yelems[i]%2 == 0) {
	ratio = pow(-grid->yexpand[i],1./(grid->yelems[i]/2-1.));
	grid->yratios[i] = ratio;
	grid->dy[i] = 0.5 * (grid->y[i] - grid->y[i-1]) * 
	  (1.-ratio) / (1.-pow(ratio,(Real)(grid->yelems[i]/2)));
      }
      else if(grid->yelems[i]%2 == 1) {
	ratio = pow(-grid->yexpand[i],1./(grid->yelems[i]/2));
	grid->yratios[i] = ratio;
	grid->dy[i] = (grid->y[i] - grid->y[i-1]) / 
	  (2.0*(1.-pow(ratio,(Real)(grid->yelems[i]/2)))/
	   (1-ratio) + pow(ratio,(Real)(grid->yelems[i]/2+0.5)));
      }
    }
  }

  /* Calculate the total number of elements */
  sumxyelems = 0;
  for(j=1;j<=grid->ycells;j++)
    for(i=1;i<=grid->xcells;i++) 
      if(grid->numbered[j][i]) sumxyelems += grid->xelems[i] * grid->yelems[j];

  grid->noelements = sumxyelems;

  if(0) printf("Created a total of %d elements\n",sumxyelems);

  grid->dx0 = dxmax;
  grid->dy0 = dymax;
}



void SetCellData(struct GridType *grid,struct CellType *cell,int info)
/* Sets the data that can directly be derived from type GridType. 
 */
{
  int i,j,cnew=0;

  for(j=1;j<= grid->ycells ;j++)                   /* cells direction up    */
    for(i=1;i<= grid->xcells; i++)                 /* cells direction right */      
      if( cnew = grid->numbered[j][i] ) {          /* if cell is occupied   */

        /* Initialize the numbering to zero  */
        cell[cnew].left1st = 0;
        cell[cnew].left2nd = 0;
        cell[cnew].leftlast = 0;
        cell[cnew].levelwidth = 0;
        cell[cnew].levelwidthcenter = 0;
	cell[cnew].leftcenter = 0;
	cell[cnew].elemwidth = 0;
	cell[cnew].elem1st = 0;
	
	/* Set the element type */
	cell[cnew].nonodes = grid->nonodes;
	cell[cnew].numbering = grid->numbering;
	
	/* Set the cell coordinates */
        cell[cnew].xcorner[TOPLEFT] = cell[cnew].xcorner[BOTLEFT] = grid->x[i-1];
	cell[cnew].xcorner[TOPRIGHT] = cell[cnew].xcorner[BOTRIGHT] = grid->x[i];
        cell[cnew].ycorner[BOTRIGHT] = cell[cnew].ycorner[BOTLEFT] = grid->y[j-1];
        cell[cnew].ycorner[TOPRIGHT] = cell[cnew].ycorner[TOPLEFT] = grid->y[j];
	
	/* Set the cell width in both directions */	
	cell[cnew].xwidth = grid->x[i] - grid->x[i-1];
	cell[cnew].ywidth = grid->y[j] - grid->y[j-1];
	
	/* Set the number of elements in a subcell */
	cell[cnew].xelem = grid->xelems[i];
	cell[cnew].yelem = grid->yelems[j];
	
	/* Set the the ratio of elements when going left to right and down to up */
	cell[cnew].xratio = grid->xratios[i];
	cell[cnew].yratio = grid->yratios[j];
	cell[cnew].xlinear = grid->xlinear[i];
	cell[cnew].ylinear = grid->ylinear[j];
	cell[cnew].xlinear = grid->xlinear[i];
	cell[cnew].ylinear = grid->ylinear[j];
	if(grid->xexpand[i] < 0.0 && grid->xlinear[i] == FALSE) {
	  if(grid->xelems[i] == 2) cell[cnew].xlinear = 1;
	  else cell[cnew].xlinear = 2;
	}
	if(grid->yexpand[j] < 0.0 && grid->ylinear[j] == FALSE) {
	  if(grid->yelems[j] == 2) cell[cnew].ylinear = 1;
	  else cell[cnew].ylinear = 2;
	}
	
	/* Set the length of the first element in both directions  */
	cell[cnew].dx1 = grid->dx[i];
	cell[cnew].dy1 = grid->dy[j];
	
	/* Set the material in question */
	cell[cnew].material = grid->structure[j][i];
	
	/* Set the boundary types for sides and corners */
	cell[cnew].boundary[LEFT] = grid->structure[j][i-1];
	cell[cnew].boundary[DOWN] = grid->structure[j-1][i];
	cell[cnew].boundary[RIGHT]= grid->structure[j][i+1];  
	cell[cnew].boundary[UP]   = grid->structure[j+1][i];
	cell[cnew].boundary[4] = grid->structure[j-1][i-1];
	cell[cnew].boundary[5] = grid->structure[j-1][i+1];
	cell[cnew].boundary[6] = grid->structure[j+1][i+1];
	cell[cnew].boundary[7] = grid->structure[j+1][i-1];

	/* Set the neighbour cell indices for sides and corners */ 
	cell[cnew].neighbour[LEFT] = grid->numbered[j][i-1];
	cell[cnew].neighbour[DOWN] = grid->numbered[j-1][i];
	cell[cnew].neighbour[RIGHT]= grid->numbered[j][i+1];  
	cell[cnew].neighbour[UP]   = grid->numbered[j+1][i];
	cell[cnew].neighbour[4] = grid->numbered[j-1][i-1];
	cell[cnew].neighbour[5] = grid->numbered[j-1][i+1];
	cell[cnew].neighbour[6] = grid->numbered[j+1][i+1];
	cell[cnew].neighbour[7] = grid->numbered[j+1][i-1];
      }

  if(info) printf("%d cells were created.\n",grid->nocells);
}



int SetCellKnots(struct GridType *grid, struct CellType *cell,int info)
/* Uses given mesh to number the knots present in the cells.
   The knots are numbered independently of the cells from left 
   to right and up to down. Only the numbers of four knots at the 
   left side of each cell are saved for later use. The number of 
   each knot can later be easily recalled using simple functions. 
   The subroutine was initially written for ordering the knots
   from left to right and down to up. However, this does not always
   produce reasonably small bandwidth and therefore a numbering 
   scheme for the other order was later added. Therefore the algorithm 
   may not be that clear for this other scheme.
   Numbers the knots in rectangular elements. There may be 
   4, 5, 8, 9, 12 or 16 nodes in each elements.
   */
{
  int i,j,level,center;
  int degree,centernodes,sidenodes,nonodes;
  int cnew=0,cup=0,cleft=0,cleftup=0;
  int elemno,knotno;
  int maxwidth,width,numbering;
  int xcells,ycells,*yelems=NULL,*xelems=NULL;

  numbering = grid->numbering;
  nonodes = grid->nonodes;
  knotno  = 0;
  elemno  = 0;

  if(numbering == NUMBER_XY) {
    xcells = grid->xcells;
    ycells = grid->ycells;
    xelems = grid->xelems;
    yelems = grid->yelems;
  }
  else if(numbering == NUMBER_YX) {
    xcells = grid->ycells;
    ycells = grid->xcells;
    xelems = grid->yelems;
    yelems = grid->xelems;
  }
  else { 
    printf("No %d numbering scheme exists!\n",numbering); 
    return(1);
  }

  switch (nonodes) {
  case 4:
    centernodes = FALSE;
    sidenodes = FALSE;
    degree = 1;
    break;
  case 5:
    centernodes = TRUE;
    sidenodes = FALSE;
    degree = 2;
    break;
  case 8:
    centernodes = FALSE;
    sidenodes = TRUE;
    degree = 2;
    break;
  case 9:
    centernodes = TRUE;
    sidenodes = TRUE;
    degree = 2;
    break;
  case 12:
    centernodes = FALSE;
    sidenodes = TRUE;
    degree = 3;
    break;
  case 16:
    centernodes = TRUE;
    sidenodes = TRUE;
    degree = 3;
    break;
  default:
    printf("CreateCells: No numbering scheme for %d-node elements.\n",
	   grid->nonodes);
    return(2);
  }

  for(j=1;j<= ycells ;j++)                     /* cells direction up     */
    for(level=0; level <= yelems[j] ;level++)  /* lines inside cells     */
      for(center=0; center < degree; center++) 
      for(i=1;i<= xcells; i++)                 /* cells direction right  */
	{        
	  if(numbering == NUMBER_XY) {	  
	    cnew = grid->numbered[j][i];
	    cleft= grid->numbered[j][i-1];
	    cup  = grid->numbered[j+1][i];
	    cleftup = grid->numbered[j+1][i-1];
	    if(cnew) {
	      cell[cnew].xind = i;
	      cell[cnew].yind = j;	      
	    }
	  }	  
	  else if(numbering == NUMBER_YX) {
	    cnew = grid->numbered[i][j];
	    cleft= grid->numbered[i-1][j];
	    cup  = grid->numbered[i][j+1];
	    cleftup = grid->numbered[i-1][j+1];
	    if(cnew) {
	      cell[cnew].xind = j;
	      cell[cnew].yind = i;	      
	    }
	  }
	  
	  if(center == 0) {
	    /* the first line of an occupied cell is unnumbered, 
	       this happens only in the bottom */
	    if (cnew && level==0 && !cell[cnew].left1st) { 
	      if(!cleft) knotno++;
	      cell[cnew].left1st = knotno;
	      if(sidenodes) 
		knotno += degree * xelems[i];
	      else
		knotno += xelems[i];
	    } /* level=0 */
	    
	    /* the n:th line of an occupied cell */
	    else if (cnew && level > 0 && level < yelems[j]) {   
	      if(!cleft) knotno++;
	      if(level==1) 
		cell[cnew].left2nd = knotno;
	      if(level==2) 
		cell[cnew].levelwidth = knotno - cell[cnew].left2nd; 	      
	      if(sidenodes) 
		knotno += degree * xelems[i];
	      else 
		knotno += xelems[i];
	    } /* 1 <= level < n */	  
	  
	    /* last line of this cell or first line of upper cell */
	    else if ((cnew || cup) && level == yelems[j]) { 
	      if(!cleft && !cleftup) 
		knotno++;
	      if(cnew)
		cell[cnew].leftlast = knotno;
	      if(level==2) 
		cell[cnew].levelwidth = knotno - cell[cnew].left2nd; 	      
	      if(yelems[j] == 1) 
		cell[cnew].left2nd = cell[cnew].leftlast;
	      if(cup)
		cell[cup].left1st = knotno;
	      if(sidenodes) 
		knotno += degree * xelems[i];
	      else 
		knotno += xelems[i];
	    } /* level=n */

	    /* Number the elements */
	    if(cnew && level > 0) {
	      if(level==1)  
		cell[cnew].elem1st = elemno+1;	      
	      if(level==2) 
		cell[cnew].elemwidth  = (elemno+1) - cell[cnew].elem1st;
	      elemno += xelems[i];	    
	    }
	  }

	  if(center && level < yelems[j]) {
	    if (cnew) { 
	      if(!cleft && sidenodes) 
		knotno++;
	      if(level==0 && center==1) 
		cell[cnew].leftcenter = knotno;
	      if(level==0 && center==2)
		cell[cnew].left2center = knotno;		
	      if(level==1 && center==1)
		cell[cnew].levelwidthcenter = knotno - cell[cnew].leftcenter;
	      if(centernodes && sidenodes)
		knotno += degree * xelems[i];
	      else 
		knotno += xelems[i];
	    } 
	  }
	
	} /* x-cell loop */
  

  maxwidth = 0;
  for(j=1;j<= ycells ;j++) 
    for(i=1;i<= xcells; i++) {
      if(numbering == NUMBER_XY)
	cnew = grid->numbered[j][i];
      else if(numbering == NUMBER_YX)
	cnew = grid->numbered[i][j];      
      if(cnew) {
	width = cell[cnew].left2nd - cell[cnew].left1st;
	if(width > maxwidth) 
	  maxwidth = width;
	width = cell[cnew].levelwidth;
	if(width > maxwidth) 
	  maxwidth = width;
	width = cell[cnew].leftlast-cell[cnew].left2nd
	  -(yelems[j]-2)*cell[cnew].levelwidth;
	if(width > maxwidth) 
	  maxwidth = width;
      }  
    }
  maxwidth += 2;

  grid->maxwidth = maxwidth;
  grid->noknots = knotno;
  grid->noelements = elemno;

  if(info) {
    printf("Numbered %d knots in %d %d-node elements.\n",knotno,elemno,nonodes);
    if(numbering == NUMBER_XY) 
      printf("Numbering order was <x><y> and max levelwidth %d.\n",
	     maxwidth);
    else if(numbering == NUMBER_YX) 
      printf("Numbering order was <y><x> and max levelwidth %d.\n",
	     maxwidth);
  }

  return(0);
}



int SetCellKnots1D(struct GridType *grid, struct CellType *cell,int info)
{
  int i;
  int degree,nonodes;
  int cnew,cleft;
  int elemno,knotno;
  int maxwidth;
  int xcells,*xelems;

  nonodes = grid->nonodes;
  knotno  = 0;
  elemno  = 0;

  xcells = grid->xcells;
  xelems = grid->xelems;

  switch (nonodes) {
  case 2:
    degree = 1;
    break;
  case 3:
    degree = 2;
    break;
  case 4:
    degree = 3;
    break;

  default:
    printf("CreateCells: No numbering scheme for %d-node elements.\n",
	   grid->nonodes);
    return(2);
  }

  for(i=1;i<= xcells; i++) {        
    
    cnew = grid->numbered[1][i];
    cleft= grid->numbered[1][i-1];
    if(cnew) cell[cnew].xind = i;
    
    if (cnew) { 
      /* Number the nodes */
      if(!cleft) knotno++;
      cell[cnew].left1st = knotno;
      knotno += degree * xelems[i];
      
      /* Number the elements */
      cell[cnew].elem1st = elemno+1;	      
      elemno += xelems[i];	    
    }
  } /* x-cell loop */
  
  maxwidth = degree;
  
  grid->maxwidth = maxwidth;
  grid->noknots = knotno;
  grid->noelements = elemno;
  
  if(info) {
    printf("Numbered %d knots in %d %d-node elements.\n",knotno,elemno,nonodes);
  }

  return(0);
}



void CreateCells(struct GridType *grid,struct CellType **cell,int info)
{
  (*cell) = (struct CellType*)
    malloc((size_t) (grid->nocells+1)*sizeof(struct CellType)); 

  SetCellData(grid,*cell,info);

  if(grid->dimension == 1) 
    SetCellKnots1D(grid,*cell,info);
  else
    SetCellKnots(grid,*cell,info);
}


void DestroyCells(struct CellType **cell)
{
  free(cell);
}



int GetKnotIndex(struct CellType *cell,int i,int j)
/* Given the cell and knot indices gives the corresponding 
   global knot index. The indices are expected to be in the 
   range [0..n] and [0..m]. Requires only the structure CellType. 
   */
{
  int ind,aid,maxj;

  if(cell->numbering == NUMBER_1D) {
    ind = cell->left1st;
    if(cell->nonodes == 2) 
      ind += i;
    else if(cell->nonodes == 3) 
      ind += 2*i;
    else if(cell->nonodes == 4) 
      ind += 3*i;
    return(ind);
  }

  if(cell->numbering == NUMBER_XY) 
    maxj = cell->yelem;
  else if(cell->numbering == NUMBER_YX) {
    aid = j; j = i; i = aid;
    maxj = cell->xelem;
  }
  else {
    maxj = 0;
    bigerror("GetKnotIndex: Unknown numbering scheme!");
  }

  if(j == 0) 
    ind = cell->left1st;
  else if (j == maxj)
    ind = cell->leftlast;
  else
    ind = cell->left2nd + (j-1) * cell->levelwidth;
  
  if(cell->nonodes == 4) 
    ind += i;
  else if(cell->nonodes == 5)
    ind += i;
  else if(cell->nonodes == 8 || cell->nonodes == 9) 
    ind += 2*i;
  else if(cell->nonodes == 12 || cell->nonodes == 16) 
    ind += 3*i;

  return(ind);
}


int GetKnotCoordinate(struct CellType *cell,int i,int j,Real *x,Real *y)
/* Given the cell and element indices inside the cell gives the 
   corresponding global knot numbers. The indices are expected 
   to be in the range [0..n] and [0..m]. Requires only the 
   structure CellType.   
   */
{
  int ind;

  if(cell->xlinear == 1)
    (*x) = cell->xcorner[BOTLEFT] + cell->dx1 * i;
  else if(cell->xlinear == 0)
    (*x) = cell->xcorner[BOTLEFT] + cell->dx1 * 
    (1.- pow(cell->xratio,(Real)(i))) / (1.-cell->xratio);
  else if(cell->xlinear == 2) {
    if(i<=cell->xelem/2) {
      (*x) = cell->xcorner[BOTLEFT] + cell->dx1 * 
	(1.- pow(cell->xratio,(Real)(i))) / (1.-cell->xratio);
    }
    else {
      (*x) = cell->xcorner[BOTRIGHT] - cell->dx1 * 
	(1.- pow(cell->xratio,(Real)(cell->xelem-i))) / (1.-cell->xratio);
    }
  }

  if(cell->ylinear == 1)
    (*y) = cell->ycorner[BOTLEFT] + cell->dy1 * j;
  else if(cell->ylinear == 0) 
    (*y) = cell->ycorner[BOTLEFT] + cell->dy1 * 
           (1.- pow(cell->yratio,(Real)(j))) / (1.-cell->yratio);
  else if(cell->ylinear == 2) {
    if(j<=cell->yelem/2) {
      (*y) = cell->ycorner[BOTLEFT] + cell->dy1 * 
	(1.- pow(cell->yratio,(Real)(j))) / (1.-cell->yratio);
    }
    else {
      (*y) = cell->ycorner[TOPLEFT] - cell->dy1 * 
	(1.- pow(cell->yratio,(Real)(cell->yelem-j))) / (1.-cell->yratio);
    }
  }

  ind = GetKnotIndex(cell,i,j);

  return(ind);
}



int GetElementIndices(struct CellType *cell,int i,int j,int *ind)
/* For given element gives the coordinates and index for each knot.
   The indices i and j are expected to range from [1..n] and [1..m]. 
   requires only the structure CellType.
   */
{
  int nonodes,numbering,elemind=0;
  
  nonodes = cell->nonodes;
  numbering = cell->numbering;

  if(numbering != NUMBER_1D) {
    ind[TOPRIGHT] = GetKnotIndex(cell,i,j);
    ind[BOTRIGHT] = GetKnotIndex(cell,i,j-1);
    ind[TOPLEFT]  = GetKnotIndex(cell,i-1,j);
    ind[BOTLEFT]  = GetKnotIndex(cell,i-1,j-1);
  }    

  if(numbering == NUMBER_XY) {
    elemind = cell->elem1st+(i-1) + (j-1)*cell->elemwidth;
    if(nonodes == 4) return(elemind);
        
    if(nonodes == 5) {
      ind[4] = cell->leftcenter + (j-1) * cell->levelwidthcenter + i;
    }
    else if(nonodes == 8) {
      ind[4] = ind[0]+1;
      ind[6] = ind[3]+1;
      ind[5] = cell->leftcenter + (j-1) * cell->levelwidthcenter + i;
      ind[7] = ind[5]-1;
    }
    else if(nonodes == 9) {
      ind[4] = ind[0]+1;
      ind[6] = ind[3]+1;
      ind[5] = cell->leftcenter + (j-1) * cell->levelwidthcenter + 2*i;
      ind[7] = ind[5]-2;
      ind[8] = ind[5]-1;
    }
    else if(nonodes == 12) {
      ind[4]  = ind[0]+1;
      ind[5]  = ind[0]+2;
      ind[9]  = ind[3]+1;
      ind[8]  = ind[3]+2;
      ind[6]  = cell->leftcenter + (j-1) * cell->levelwidthcenter + i;
      ind[11] = ind[6]-1;
      ind[7]  = cell->left2center + (j-1) * cell->levelwidthcenter + i;
      ind[10] = ind[7]-1;
    }
    else if(nonodes == 16) {
      ind[4]  = ind[0]+1;
      ind[5]  = ind[0]+2;
      ind[9]  = ind[3]+1;
      ind[8]  = ind[3]+2;
      ind[6]  = cell->leftcenter + (j-1) * cell->levelwidthcenter + 3*i;
      ind[11] = ind[6]-3;
      ind[7]  = cell->left2center + (j-1) * cell->levelwidthcenter + 3*i;
      ind[10] = ind[7]-3;
      ind[13] = ind[6]-1;
      ind[12] = ind[6]-2;
      ind[14] = ind[7]-1;
      ind[15] = ind[7]-2;
    }
    else 
      printf("GetElementIndices: not implemented for %d nodes.\n",nonodes);    
  }

  else if(numbering == NUMBER_YX) { 
    elemind = cell->elem1st+(j-1) + (i-1)*cell->elemwidth;
    
    if(nonodes == 4) return(elemind);

    if(nonodes == 5) {
      ind[4] = cell->leftcenter + (i-1) * cell->levelwidthcenter + j;
    }
    else if (nonodes==8) {
      ind[7] = ind[0]+1;
      ind[5] = ind[1]+1;
      ind[6] = cell->leftcenter + (i-1) * cell->levelwidthcenter + j;
      ind[4] = ind[6]-1;
    }
    else if(nonodes == 9) {
      ind[7] = ind[0]+1;
      ind[5] = ind[1]+1;
      ind[6] = cell->leftcenter + (i-1) * cell->levelwidthcenter + 2*j;
      ind[4] = ind[6]-2;
      ind[8] = ind[6]-1;
    }
    else if(nonodes == 12) {
      ind[11]= ind[0]+1;
      ind[10]= ind[0]+2;
      ind[6] = ind[1]+1;
      ind[7] = ind[1]+2;
      ind[9] = cell->leftcenter + (i-1) * cell->levelwidthcenter + j;
      ind[4] = ind[9]-1;
      ind[8] = cell->left2center + (i-1) * cell->levelwidthcenter + j;
      ind[5] = ind[8]-1;
    }
    else if(nonodes == 16) {
      ind[11]= ind[0]+1;
      ind[10]= ind[0]+2;
      ind[6] = ind[1]+1;
      ind[7] = ind[1]+2;
      ind[9] = cell->leftcenter + (i-1) * cell->levelwidthcenter + 3*j;
      ind[4] = ind[9]-3;
      ind[15]= ind[9]-1;
      ind[12]= ind[9]-2;
      ind[8] = cell->left2center + (i-1) * cell->levelwidthcenter + 3*j;
      ind[5] = ind[8]-3;
      ind[14]= ind[8]-1;
      ind[13]= ind[8]-2;
    }
    else 
      printf("GetElementIndices: not implemented for %d nodes.\n",nonodes);
  }

  else if(numbering == NUMBER_1D) {
    elemind = cell->elem1st+(i-1);
    ind[0] = GetKnotIndex(cell,i-1,1);
    if(nonodes == 2) {
      ind[1] = ind[0] + 1;
    }
    else if(nonodes == 3) {
      ind[2] = ind[0] + 1;
      ind[1] = ind[0] + 2;
    }
    else if(nonodes == 4) {
      ind[2] = ind[0] + 1;
      ind[3] = ind[0] + 2;
      ind[1] = ind[0] + 3;
    }
  }  
  return(elemind);
}


int GetElementIndex(struct CellType *cell,int i,int j)
/* For given element gives the element index. 
   The indices i and j are expected to range from [1..n] and [1..m]. 
   requires only the structure CellType.
   */
{
  int elemind=0;
 
  if(cell->numbering == NUMBER_XY) 
    elemind = cell->elem1st+(i-1) + (j-1)*cell->elemwidth;
  else if(cell->numbering == NUMBER_YX) 
    elemind = cell->elem1st+(j-1) + (i-1)*cell->elemwidth;

  return(elemind);
}



int GetElementCoordinates(struct CellType *cell,int i,int j,
			  Real *globalcoord,int *ind)
/* For given element gives the coordinates and index for each knot.
   The indices i and j are expected to range from [1..n] and [1..m]. 
   requires only the structure CellType
   Note that this subroutine assumes that the elements are
   rectangular.
   */
{
  int k,nonodes,numbering,elemind=0;
  Real xrat,yrat;

  k = nonodes = cell->nonodes;
  numbering = cell->numbering;

  if(numbering == NUMBER_1D) {
    elemind = cell->elem1st+(i-1);    
    ind[0] = GetKnotCoordinate(cell,i-1,j-1,&globalcoord[0],
			       &globalcoord[nonodes]);
    ind[1] = GetKnotCoordinate(cell,i,j-1,&globalcoord[1],
			       &globalcoord[nonodes+1]);    
    
    if(nonodes == 3) {
      globalcoord[2] = (globalcoord[0]+globalcoord[1])/2.0;
      globalcoord[5] = (globalcoord[3]+globalcoord[4])/2.0;
      ind[2] = ind[0] + 1;
    }
    else if(nonodes == 4) {
      globalcoord[2] = (2.0*globalcoord[0]+globalcoord[1])/3.0;
      globalcoord[6] = (2.0*globalcoord[4]+globalcoord[5])/3.0;
      ind[2] = ind[0] + 1;
      globalcoord[3] = (globalcoord[0]+2.0*globalcoord[1])/3.0;
      globalcoord[7] = (globalcoord[4]+2.0*globalcoord[5])/3.0;
      ind[3] = ind[0] + 2;
    }

    return(elemind);
  }
  else if(numbering == NUMBER_XY) 
    elemind = cell->elem1st+(i-1) + (j-1)*cell->elemwidth;
  else if(numbering == NUMBER_YX)
    elemind = cell->elem1st+(j-1) + (i-1)*cell->elemwidth;


  ind[TOPRIGHT] = GetKnotCoordinate(cell,i,j,&globalcoord[TOPRIGHT],
				    &globalcoord[TOPRIGHT+nonodes]);
  ind[BOTRIGHT] = GetKnotCoordinate(cell,i,j-1,&globalcoord[BOTRIGHT],
				    &globalcoord[BOTRIGHT+nonodes]);
  ind[TOPLEFT] = GetKnotCoordinate(cell,i-1,j,&globalcoord[TOPLEFT],
				    &globalcoord[TOPLEFT+nonodes]);
  ind[BOTLEFT] = GetKnotCoordinate(cell,i-1,j-1,&globalcoord[BOTLEFT],
				    &globalcoord[BOTLEFT+nonodes]);
  if(nonodes == 4) return(elemind);

  GetElementIndices(cell,i,j,ind);

#if 0  
  if(cell->xlinear) 
    xrat = 1.0;
  else
    xrat = sqrt(cell->xratio);

  if(cell->ylinear) 
    yrat = 1.0;
  else 
    yrat = sqrt(cell->yratio);
#endif
  xrat = yrat = 1.0;

  if(nonodes == 5) {
    globalcoord[4] = 0.5*(globalcoord[0]+globalcoord[2]);
    globalcoord[4+k] = 0.5*(globalcoord[k]+globalcoord[2+k]);
  }
  else if(nonodes == 8 || nonodes == 9) {
    globalcoord[4] = (xrat*globalcoord[0]+globalcoord[1])/(1+xrat);
    globalcoord[4+k] = globalcoord[k];    
    globalcoord[5] = globalcoord[1];
    globalcoord[5+k] = 
      (yrat*globalcoord[1+k]+globalcoord[2+k])/(1+yrat);    
    globalcoord[6] = globalcoord[4];
    globalcoord[6+k] = globalcoord[2+k];
    globalcoord[7] = globalcoord[3];
    globalcoord[7+k] = globalcoord[5+k];
    if(nonodes == 9) {
      globalcoord[8] = globalcoord[4];
      globalcoord[8+k] = globalcoord[5+k];
    }
  }
  else if(nonodes == 12 || nonodes == 16) {
    globalcoord[4] = (2.*globalcoord[0]+globalcoord[1])/3.0;
    globalcoord[4+k] = (2.*globalcoord[k]+globalcoord[1+k])/3.0;
    globalcoord[5] = (2.*globalcoord[1]+globalcoord[0])/3.0;
    globalcoord[5+k] = (2.*globalcoord[1+k]+globalcoord[k])/3.0;
    globalcoord[6] = (2.*globalcoord[1]+globalcoord[2])/3.0;
    globalcoord[6+k] = (2.*globalcoord[1+k]+globalcoord[2+k])/3.0;
    globalcoord[7] = (2.*globalcoord[2]+globalcoord[1])/3.0;
    globalcoord[7+k] = (2.*globalcoord[2+k]+globalcoord[1+k])/3.0;
    globalcoord[8] = (2.*globalcoord[2]+globalcoord[3])/3.0;
    globalcoord[8+k] = (2.*globalcoord[2+k]+globalcoord[3+k])/3.0;
    globalcoord[9] = (2.*globalcoord[3]+globalcoord[2])/3.0;
    globalcoord[9+k] = (2.*globalcoord[3+k]+globalcoord[2+k])/3.0;
    globalcoord[10] = (2.*globalcoord[3]+globalcoord[0])/3.0;
    globalcoord[10+k] = (2.*globalcoord[3+k]+globalcoord[k])/3.0;
    globalcoord[11] = (2.*globalcoord[0]+globalcoord[3])/3.0;
    globalcoord[11+k] = (2.*globalcoord[k]+globalcoord[3+k])/3.0;
    if(nonodes == 16) {
      globalcoord[12] = (2.*globalcoord[11]+globalcoord[6])/3.0;
      globalcoord[12+k] = (2.*globalcoord[11+k]+globalcoord[6+k])/3.0;
      globalcoord[13] = (2.*globalcoord[6]+globalcoord[11])/3.0;
      globalcoord[13+k] = (2.*globalcoord[6+k]+globalcoord[11+k])/3.0;
      globalcoord[14] = (2.*globalcoord[7]+globalcoord[10])/3.0;
      globalcoord[14+k] = (2.*globalcoord[7+k]+globalcoord[10+k])/3.0;
      globalcoord[15] = (2.*globalcoord[10]+globalcoord[7])/3.0;
      globalcoord[15+k] = (2.*globalcoord[10+k]+globalcoord[7+k])/3.0;      
    }
  }

#if 0
  for(i=0;i<k;i++) 
    printf("ind[%d]=%d  x=%.4lg  y=%.4lg\n",i,ind[i],
	   globalcoord[i],globalcoord[i+k]);
#endif

  return(elemind);
}


int GetSideInfo(struct CellType *cell,int cellno,int side,int element,
		int *elemind)
/* Given the side and element numbers returns the indices of 
   the side element. When the end of the side is reached the function 
   returns FALSE, otherwise TRUE.
   */
{
  int more,sideno;                 
  int ind[MAXNODESD2];

  more = TRUE;

  if(cell[cellno].numbering == NUMBER_1D) {

    switch(side) {
    case 0:
      more = FALSE;
      elemind[0] = GetElementIndices(&(cell[cellno]),1,element,&(ind[0]));
      if(sideno = cell[cellno].neighbour[LEFT])
	elemind[1] = GetElementIndices(&(cell[sideno]),cell[sideno].xelem,element,&(ind[0]));
      else 
	elemind[1] = 0;
      break;
     
      
    case 1:
      more = FALSE; 
      elemind[0] = GetElementIndices(&(cell)[cellno],cell[cellno].xelem,element,&(ind[0]));
      if(sideno = cell[cellno].neighbour[RIGHT])
	elemind[1] = GetElementIndices(&(cell[sideno]),1,element,&(ind[0]));
      else 
	elemind[1] = 0;
      break;
    }

    return(more);
  }


  switch(side) {

  case LEFT:
    if(element == cell[cellno].yelem) more = FALSE;
    elemind[0] = GetElementIndices(&(cell[cellno]),1,element,&(ind[0]));
    if(sideno = cell[cellno].neighbour[LEFT])
      elemind[1] = GetElementIndices(&(cell[sideno]),cell[sideno].xelem,element,&(ind[0]));
    else 
      elemind[1] = 0;
    break;

  case RIGHT:
    if(element == cell[cellno].yelem) more = FALSE; 
    elemind[0] = GetElementIndices(&(cell)[cellno],cell[cellno].xelem,element,&(ind[0]));
    if(sideno = cell[cellno].neighbour[RIGHT])
      elemind[1] = GetElementIndices(&(cell[sideno]),1,element,&(ind[0]));
    else 
      elemind[1] = 0;
    break;

  case DOWN:
    if(element == cell[cellno].xelem) more = FALSE;
    elemind[0] = GetElementIndices(&(cell)[cellno],element,1,&(ind[0]));
    if(sideno = cell[cellno].neighbour[DOWN])
      elemind[1] = GetElementIndices(&(cell[sideno]),element,cell[sideno].yelem,&(ind[0]));
    else 
      elemind[1] = 0;
    break;

  case UP:
    if(element == cell[cellno].xelem) more = FALSE; 
    elemind[0] = GetElementIndices(&(cell)[cellno],element,cell[cellno].yelem,&(ind[0]));
    if(sideno = cell[cellno].neighbour[UP])
      elemind[1] = GetElementIndices(&(cell[sideno]),element,1,&(ind[0]));
    else 
      elemind[1] = 0;
    break;
    
  default:
    printf("Impossible option in GetSideInfo.\n");
  }

  return(more);
}



void SetElementDivisionExtruded(struct GridType *grid,int info)
{
  int i,nzmax = 0,sumzelems;
  Real ratio,linearlimit;
  Real dzmax = 0,dz = 0;
  
  linearlimit = 0.001;

  if(grid->autoratio) {
    for(i=1;i<=grid->zcells;i++)  
      grid->zelems[i] = grid->minzelems;
  }
  else {
    for(i=1;i<=grid->zcells;i++)  
      if(grid->zelems[i] < grid->minzelems) grid->zelems[i] = grid->minzelems;
  }

  sumzelems = grid->zcells * grid->minzelems;
  if(sumzelems > grid->totzelems) {
#if 0
    printf("SetElementDivision: %d is too few elements in z-direction (min. %d)\n",
	   grid->totzelems,sumzelems);
#endif
    grid->totzelems = sumzelems;
  }      

  /* Put the linearity flags. */
  for(i=1; i<= grid->zcells ;i++) {
    if (fabs(1.-grid->zexpand[i]) < linearlimit) 
      grid->zlinear[i] = TRUE; 
    else 
      grid->zlinear[i] = FALSE;
  }

  if(grid->autoratio) {
    int active;
    for(;;) {
      dzmax = 0.0;
      active = FALSE;

      for(i=1;i<=grid->zcells;i++) {
	if(grid->zelems[i] == 0) continue;
	if(grid->zlinear[i] == TRUE)
	  dz = (grid->z[i] - grid->z[i-1])/(grid->zelems[i]+1);
	else {
	  if(grid->zexpand[i] > 0.0) {
	    ratio = pow(grid->zexpand[i],1./(1.*grid->zelems[i]));
	    dz = (grid->z[i] - grid->z[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)(grid->zelems[i]+1)));
	    if(ratio < 1.)   
	      dz *= grid->zexpand[i];
	  }
	  else if(grid->zelems[i]==1) {
	    dz = (grid->z[i] - grid->z[i-1])/(grid->zelems[i]+1);
	  } 
	  else if((grid->zelems[i]+1)%2 == 0) {
	    ratio = pow(-grid->zexpand[i],1./((grid->zelems[i]+1)/2-1.));
	    dz = 0.5 * (grid->z[i] - grid->z[i-1]) * 
	      (1.-ratio) / (1.-pow(ratio,(Real)((grid->zelems[i]+1)/2)));
	  }
	  else if((grid->zelems[i]+1)%2 == 1) {
	    ratio = pow(-grid->zexpand[i],1./((grid->zelems[i]+1)/2));
	    dz = (grid->z[i] - grid->z[i-1]) / 
	      (2.0*(1.-pow(ratio,(Real)((grid->zelems[i]+1)/2)))/
	       (1-ratio) + pow(ratio,(Real)((grid->zelems[i]+1)/2+0.5)));
	  }
	}
	dz *= grid->zdens[i] / grid->xzratio;

	if(dz > dzmax) {
	  dzmax = dz;
	  nzmax = i;
	}

	if(grid->autoratio) {
	  if(fabs(dz - grid->dx0) < fabs( dz*(1.0/(1.0-1.0/grid->zelems[i])) - grid->dx0) ) {
	    grid->zelems[i] += 1;
	    sumzelems++;
	    active = TRUE;
	  }
	}
      }

      if(grid->autoratio) {
	if(!active) break;
      }
      else {
	grid->zelems[nzmax] += 1;
	sumzelems++;
	if(sumzelems >= grid->totzelems) break;
      }
    }
  }

  /* Set the size of the first element within each subcell */
  grid->totzelems = 0;
  for(i=1;i<=grid->zcells;i++) {
    grid->totzelems += grid->zelems[i];
    if(grid->zlinear[i] == TRUE || grid->zelems[i]==1) 
      grid->dz[i] = (grid->z[i] - grid->z[i-1])/grid->zelems[i];
    else if(grid->zexpand[i] > 0.0) {
      ratio = pow(grid->zexpand[i],1./(grid->zelems[i]-1.));
      grid->zratios[i] = ratio;
      grid->dz[i] = (grid->z[i] - grid->z[i-1]) * 
	(1.-ratio) / (1.-pow(ratio,(Real)(grid->zelems[i])));
    }	
    else if(grid->zelems[i] == 2) {
      grid->dz[i] = (grid->z[i] - grid->z[i-1])/grid->zelems[i];
    }
    else if(grid->zelems[i]%2 == 0) {
      ratio = pow(-grid->zexpand[i],1./(grid->zelems[i]/2-1.));
      grid->zratios[i] = ratio;
	grid->dz[i] = 0.5 * (grid->z[i] - grid->z[i-1]) * 
	  (1.-ratio) / (1.-pow(ratio,(Real)(grid->zelems[i]/2)));
    }
    else if(grid->zelems[i]%2 == 1) {
      ratio = pow(-grid->zexpand[i],1./(grid->zelems[i]/2));
      grid->zratios[i] = ratio;
      grid->dz[i] = (grid->z[i] - grid->z[i-1]) / 
	(2.0*(1.-pow(ratio,(Real)(grid->zelems[i]/2)))/
	 (1-ratio) + pow(ratio,(Real)(grid->zelems[i]/2+0.5)));
    }
  }

  if(info) printf("Created %d extruded divisions.\n",
		  grid->totzelems);
  grid->dz0 = dzmax;
}



void SetElementDivisionCylinder(struct GridType *grid,int info)
{
  int i,k;
  Real ratio,eps;
  Real dzmax;
  
  eps = 1.0e-8;

  grid->zcells = 2*grid->rotateblocks;
  grid->z[0] = 0.0;
  for(i=0;grid->x[i]<eps;i++); k=i;  
  grid->rotateradius1 = grid->x[k];

  if(grid->rotateradius2 < 0.0) {
    grid->rotatecartesian = TRUE;
    grid->rotateradius2 = 1.0e6;
  }
  if(grid->rotateradius2 <= sqrt(2.0)*grid->rotateradius1) {
    if(k+1 <= grid->xcells) 
      grid->rotateradius2 = grid->x[k+1];
    grid->rotateradius2 = MAX(sqrt(2.0)*grid->rotateradius1,
			      grid->rotateradius2);
  }

  if(!grid->xlinear[k]) 
    printf("SetElementDivisionCylinder: The division must be linear!\n");

  for(i=1;i<=grid->zcells;i++) {
    grid->z[i] = i*(2.0*FM_PI)/8.0; 
    grid->zelems[i]  = grid->xelems[k];
    grid->zlinear[i] = grid->xlinear[k];
    grid->zratios[i] = grid->xratios[k];
    grid->zexpand[i] = grid->xexpand[k];
    grid->zmaterial[i] = 0;
  }

  grid->totzelems =  grid->zcells * grid->xelems[k];

  dzmax = 0.0;
  for(i=1;i<=grid->zcells;i++) {
    if(grid->zlinear[i] == TRUE || grid->zelems[i]==1) 
      grid->dz[i] = (grid->z[i] - grid->z[i-1])/grid->zelems[i];
    else if(grid->zexpand[i] > 0.0) {
      ratio = pow(grid->zexpand[i],1./(grid->zelems[i]-1.));
      grid->zratios[i] = ratio;
      grid->dz[i] = (grid->z[i] - grid->z[i-1]) * 
	(1.-ratio) / (1.-pow(ratio,(Real)(grid->zelems[i])));
    }	
    else if(grid->zelems[i] == 2) {
      grid->dz[i] = (grid->z[i] - grid->z[i-1])/grid->zelems[i];
    }
    else if(grid->zelems[i]%2 == 0) {
      ratio = pow(-grid->zexpand[i],1./(grid->zelems[i]/2-1.));
      grid->zratios[i] = ratio;
	grid->dz[i] = 0.5 * (grid->z[i] - grid->z[i-1]) * 
	  (1.-ratio) / (1.-pow(ratio,(Real)(grid->zelems[i]/2)));
    }
    else if(grid->zelems[i]%2 == 1) {
      ratio = pow(-grid->zexpand[i],1./(grid->zelems[i]/2));
      grid->zratios[i] = ratio;
      grid->dz[i] = (grid->z[i] - grid->z[i-1]) / 
	(2.0*(1.-pow(ratio,(Real)(grid->zelems[i]/2)))/
	 (1-ratio) + pow(ratio,(Real)(grid->zelems[i]/2+0.5)));
    }

    if(dzmax < grid->dz[i]) dzmax = grid->dz[i];
  }

  if(info) printf("Created %d divisions in %d cells for rotation [%.2g  %.2g].\n",
		  grid->totzelems,grid->zcells,
		  grid->rotateradius1,grid->rotateradius2);
  grid->dz0 = dzmax;
}


int InspectElement(struct FemType *data,int idx)
{
  int elemtype,nonodes,i,j,k;
  Real *x,*y,*z;
  Real ds,minds,xc,yc,zc;
  Real x0,y0,z0;
  int ii,jj,mini,minj,kk;
  
  elemtype = data->elementtypes[idx];
  printf("Inspecting element %d of type %d\n",idx,elemtype);

  nonodes = elemtype % 100;
  x = data->x;
  y = data->y;
  z = data->z;

  for(k=0;k<nonodes;k++) {
    kk = data->topology[idx][k];
    x0 = x[kk];
    y0 = y[kk];
    z0 = z[kk];
    mini = -1;
    minj = -1;
    minds = 1.0e20;

    for(i=0;i<nonodes;i++) {
      ii = data->topology[idx][i];
      for(j=i+1;j<nonodes;j++) {
	jj = data->topology[idx][j];
	xc = (x[ii]+x[jj])/2;
	yc = (y[ii]+y[jj])/2;
	zc = (z[ii]+z[jj])/2;
	ds = sqrt( pow(x0-xc,2) + pow(y0-yc,2) + pow(z0-zc,2) );
	if( mini == -1 || ds < minds) {
	  mini = i;
	  minj = j;
	  minds = ds;
	}
      }
    }      

    printf("k=%d mini=%d minj=%d ds=%.3e\n",k+1,mini+1,minj+1,minds);
  }

  return(0);
}
