/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author:  Peter Raback
   Email:   elmeradm@csc.fi
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

/* --------------------:  egextra.c  :--------------------------

   These functions provide the user of the program information 
   about how the mesh was created and what the resulting sparse 
   matrix will be alike and also present the results of the 
   calculations in various ways. These subroutines don't affect 
   the operation and results of the program and can thus be 
   omitted if not needed. 
   */
   
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>

#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egmesh.h"
#include "egnative.h"
#include "egextra.h"
#include "egparallel.h" 



int SaveCellInfo(struct GridType *grid,struct CellType *cell,
		 char *prefix,int info)
/* Save some (int) values for each cell in structure CellType. 
   The resulting matrix may be used to check that the division to 
   subcells is as wanted.
   */
{
  int i; 
  FILE *out;
  char filename[MAXFILESIZE];

  AddExtension(prefix,filename,"cell");
  out = fopen(filename,"w");

  fprintf(out,"%-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n",
	  "1st","2nd","last","level","center","mat","xlin","ylin","num");
  for(i=1;i<=grid->nocells;i++) {
    fprintf(out,"%-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d\n",
	    cell[i].left1st,cell[i].left2nd,cell[i].leftlast,
	    cell[i].levelwidth,cell[i].leftcenter,
	    cell[i].material,cell[i].xlinear,cell[i].ylinear,cell[i].numbering);
  }
  fclose(out);

  if(info) printf("The cell information was saved to file %s.\n",filename);

  return(0);
}



int SaveBoundary(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info)
/* This function saves the given boundary to an ascii-file.
   The data is saved in format [x1 x2 y1 y2 v1 v2].
   */
{
  int i,k,sideelemtype;
  FILE *out;
  char filename[MAXFILESIZE];
  int sideind[MAXNODESD1]; 

  if(!bound->created) {
    printf("SaveBoundary: You tried to save a nonexisting boundary.\n");
    return(1);
  }
  if(bound->nosides == 0) return(0);

  AddExtension(prefix,filename,"bound");
  out = fopen(filename,"w");

  for(i=1; i <= bound->nosides; i++) {

    GetElementSide(bound->parent[i],bound->side[i],bound->normal[i],data,sideind,&sideelemtype);

    fprintf(out,"%-12.4le %-12.4le %-12.4le %-12.4le ",
	    data->x[sideind[0]],data->x[sideind[1]],
	    data->y[sideind[0]],data->y[sideind[1]]);
    for(k=0;k<MAXVARS;k++) 
      if(bound->evars[k]) {
	if(bound->points[k] == 1) 
	  fprintf(out,"%-10.4le ",bound->vars[k][i]);
      }		
    fprintf(out,"\n");
  }

  fclose(out);

  if(info) printf("Boundary information was saved to file %s.\n",filename);

  return(0);
}


int CreateBoundaryChain(struct FemType *data,struct BoundaryType *bound,int info)
{
  int i,sideind[MAXNODESD1],sideind2[MAXNODESD1]; 
  int size,indfirst,indsecond,side,setchain,sideelemtype;

  if(bound->created == FALSE) {
    if(info) printf("CreateBoundaryChain: boundary not created!\n");
    return(1);
  }
  if(bound->nosides < 2) {
    if(info) printf("CreateBoundaryChain: there must be at least 2 sides!\n");
    return(2);
  }
  if(bound->echain == TRUE) {
    if(0) printf("CreateBoundaryChain: chain already exists!\n");
    return(3);
  }

  setchain = FALSE;
restart:
  /* First determine the way that the chain starts. */
  GetElementSide(bound->parent[1],bound->side[1],bound->normal[1],data,sideind,&sideelemtype);
  GetElementSide(bound->parent[2],bound->side[2],bound->normal[2],data,sideind2,&sideelemtype);
  if(sideind[1] == sideind2[0]) {
    indfirst = sideind[0];
    indsecond = sideind[1];
  }
  else if(sideind[0] == sideind2[1]) {
    indfirst = sideind[1];
    indsecond = sideind[0];
  }
  else {
    printf("CreateBoundaryChain: Impossible case!\n");
    indfirst = sideind[0];
    indsecond = sideind[1];
  }

  if(setchain) {
    bound->chain[0] = indfirst;
    bound->chain[1] = indsecond;
  }
  size = 1;
  side = 1;

  /* The chain is followed in positive direction while successful,
     if not the direction is changed and so on... */
  for(;;) {	
    for(i=side+1;i<=bound->nosides;i++) {
      GetElementSide(bound->parent[i],bound->side[i],bound->normal[i],data,sideind,&sideelemtype);
      if(sideind[0] == indsecond) {
	indsecond = sideind[1];
	goto jump;
      }
      else if(sideind[1] == indsecond) {
	indsecond = sideind[0];
	goto jump;
      }
    }
	
    for(i=side-1;i>=1;i--) {
      GetElementSide(bound->parent[i],bound->side[i],bound->normal[i],data,sideind,&sideelemtype);
      if(sideind[0] == indsecond) {
	indsecond = sideind[1];
	goto jump;
      }
      else if(sideind[1] == indsecond) {
	indsecond = sideind[0];
	goto jump;
      }
    } 
    goto theend;
    
  jump:
    size++;
    if(setchain) 
      bound->chain[size] = indsecond;
    side = i;
    if(indsecond == indfirst) goto theend;
  }
  
theend:  
  if(setchain == FALSE) {
    if(size != bound->nosides) 
      printf("CreateBoundaryChain: the boundary is not continuous (%d vs. %d)\n",
	     size,bound->nosides);
    else 
      printf("CreateBoundaryChain: the boundary is continuous with %d nodes.\n",size+1);
    bound->chain = Ivector(0,size);
    bound->echain = TRUE;
    bound->chainsize = size;
    setchain = TRUE;
    goto restart;
  }
  return(0);
}


int SaveBoundariesChain(struct FemType *data,struct BoundaryType *bound,
			char *prefix,int info)
/* This function saves the given boundary to an ascii-file.
   The data is saved in format [x y v].
   This may be used particularly for boundaries that 
   are continues i.e. the element sides constitute a full chain. 
   */
{
  int i,j,k,ind,length,col;
  FILE *out;
  char filename[MAXFILESIZE];
  char filename2[MAXFILESIZE];

  for(j=0;j<MAXBOUNDARIES;j++) 
    if(bound[j].created && bound[j].nosides >0) {

      CreateBoundaryChain(data,&bound[j],info);
      length = bound[j].chainsize;
      if(length < 2) continue;

      sprintf(filename,"%s%d%s",prefix,j,".side");
      out = fopen(filename,"w");

      for(i=0;i<=length;i++) {
	ind = bound[j].chain[i];
	fprintf(out,"%-10.4le %-10.4le %-6d ",
		data->x[ind],data->y[ind],ind);
	for(k=0;k<MAXVARS;k++) 
	  if(bound[j].evars[k]) {
	    if(bound[j].points[k] == 0)
	      fprintf(out,"%-10.4le ",bound[j].vars[k][i]);
	    else if(bound[j].points[k] == 1) {
	      if(i==0)
		fprintf(out,"%-10.4le ",bound[j].vars[k][1]);		
	      else if(i==length)
		fprintf(out,"%-10.4le ",bound[j].vars[k][length]);		
	      else
		fprintf(out,"%-10.4le ",
			0.5*(bound[j].vars[k][i]+bound[j].vars[k][i+1]));	
	    }
	  }

	for(k=0;k<MAXDOFS;k++) {
	  if(data->edofs[k] == 1) 
	    fprintf(out,"%-10.4le  ",data->dofs[k][ind]);
	  if(data->edofs[k] == 2) 
	    fprintf(out,"%-10.4le  %-10.4le  ",
		    data->dofs[k][2*ind-1],data->dofs[k][2*ind]);
	}
	
	fprintf(out,"\n");
      }
      fclose(out);

      sprintf(filename2,"%s%d%s",prefix,j,".sidetxt");
      out = fopen(filename2,"w");
      fprintf(out,"Degrees of freedom in file %s are as follows:\n",filename);
      if(bound->coordsystem == COORD_CART2) {
	fprintf(out,"col1: X coordinate\n");
	fprintf(out,"col2: Y coordinate\n");
      }
      else if(bound->coordsystem == COORD_AXIS) {
	fprintf(out,"col1: R coordinate\n");
	fprintf(out,"col2: Z coordinate\n");
      }
      else if(bound->coordsystem == COORD_POLAR) {
	fprintf(out,"col1: R coordinate\n");
	fprintf(out,"col2: F coordinate\n");
      }
      fprintf(out,"col3: node indices\n");
      col = 3;
      for(k=0;k<MAXVARS;k++) 
	if(bound[j].evars[k] && bound[j].points[k] <= 1) 
	  fprintf(out,"col%d: %s\n",++col,bound[j].varname[k]);	  
      for(k=0;k<MAXDOFS;k++) {
	if(data->edofs[k] == 1) 
	  fprintf(out,"col%d: %s\n",++col,data->dofname[k]);	  
	if(data->edofs[k] == 2) {
	  fprintf(out,"col%d: %s1\n",++col,data->dofname[k]);	  
	  fprintf(out,"col%d: %s2\n",++col,data->dofname[k]);	  
	}
      }
      fclose(out);

      if(info) printf("Boundary info was saved to files %s and %s.\n",
		      filename,filename2);
    }

  return(0);
}



int SaveBoundaryForm(struct FemType *data,struct CellType *cell,
		     char* prefix,int info)
/* Saves the form of the boundary as given by function GetSideInfo 
   in form [x,y,index]. 
   */
{
  int sideknots,elemno,side,more,elemind[2],sideelemtype;
  int no,sideind[MAXNODESD1];
  FILE *out;
  char filename[MAXFILESIZE];

  if(data->created == FALSE) {
    printf("SaveBoundaryForm: structure FemType not created\n");
    return(1);
  }

  sideknots = 0;
  more = FALSE;

  AddExtension(prefix,filename,"boundary");
  out = fopen(filename,"w");

  /* Go through all pairs of points and save them into a matrix. */
  for(no=1; no <= data->nocells; no++)
    for(side=0; side < 4; side++) 
      if(cell[no].material != cell[no].boundary[side]) {
        elemno = 0; 
        do { 
          elemno++;
	  sideknots++;
          more = GetSideInfo(cell,no,side,elemno,elemind);
	  GetElementSide(elemind[0],side,1,data,sideind,&sideelemtype);
	  
	  fprintf(out,"%-12.4le %-12.4le %-12.4le %-12.4le\n",
		  data->x[sideind[0]],data->x[sideind[1]],
		  data->y[sideind[0]],data->y[sideind[1]]);
        } while(more);
      }

  fclose(out);

  if(info) printf("%d boundaries between materials were saved to file %s.\n",
		  sideknots,filename);
  return(0);
}


int SaveBoundaryLine(struct FemType *data,int direction,
		     Real c0,char* prefix,int info)
/* Saves the nodes forming a vertical or a horizontal line. 
   The format is [x y v1 v2 v3 ...].
   */
{
  int k,no,points;
  FILE *out;
  char filename[MAXFILESIZE];
  Real c,c1,eps;

  if(data->created == FALSE) {
    printf("SaveBoundaryLine: structure FemType not created\n");
    return(1);
  }

  eps = 1.0e-6;
  points = 0;

  AddExtension(prefix,filename,"line");
  out = fopen(filename,"w");

  /* Go through all pairs of points and save them into amatrix. */
  
  c1 = data->x[1];
  for(no=1; no <= data->noknots; no++) {
    if(direction > 0) 
      c = data->x[no];
    else
      c = data->y[no];
    if(fabs(c-c0) < fabs(c1-c0))
      c1 = c;
  }
  for(no=1; no <= data->noknots; no++) {
    if(direction > 0) 
      c = data->x[no];
    else
      c = data->y[no];
    if(fabs(c-c1) < eps) {
      if(direction > 0) 
	fprintf(out,"%-12.7le %-12.7le ",c,data->y[no]);
      else 
	fprintf(out,"%-12.7le %-12.7le ",data->x[no],c);	
      for(k=0;k<MAXDOFS;k++) {
	if(data->edofs[k] == 1) 
	  fprintf(out,"%-12.7le ",data->dofs[k][no]);
	if(data->edofs[k] == 2) 
	  fprintf(out,"%-12.7le %-12.7le ",
		  data->dofs[k][2*no-1],data->dofs[k][2*no]);
      }	
      fprintf(out,"\n");
      points++;
    }
  }
  if(info) printf("Line (c=%.3lg) with %d nodes was saved to file %s.\n",
		  c1,points,filename);
  
  fclose(out);
  
  return(0);
}


int SaveSubcellForm(struct FemType *data,struct CellType *cell, 
		    char* prefix,int info)
/* Saves the form of the boundary as given by function GetSideInfo 
   in form [x,y,index]. 
   */
{
  int more,sideknots,elemno,elemind[2],side,i;
  int no,sideind[MAXNODESD1],nosidenodes,sidelemtype;
  FILE *out;
  char filename[MAXFILESIZE];

  sideknots = 0;
  more = FALSE;

  if(data->created == FALSE) {
    printf("SaveSubcellForm: structure FemType not created\n");
    return(1);
  }

  AddExtension(prefix,filename,"subcell");
  out = fopen(filename,"w");

  /* Go through all pairs of points and save them into amatrix. */
  for(no=1; no <= data->nocells; no++)
    for(side=0; side < 4; side++) 
      if(cell[no].boundary[side]) {
        elemno = 0; 
        do { 
          elemno++;
	  sideknots++;
          more = GetSideInfo(cell,no,side,elemno,elemind);
	  GetElementSide(elemind[0],side,1,data,sideind,&sidelemtype);
	  nosidenodes = sidelemtype%100;
	  for(i=0;i<nosidenodes;i++) 
	    fprintf(out,"%-12.4le %-12.4le %-8d\n",
		    data->x[sideind[i]],data->y[sideind[i]],sideind[i]);
        } while(more);
      }
  fclose(out);

  if(info) printf("There are %d sideknots in the elements.\n",sideknots);
  if(info) printf("The positions of the sideknots were saved in file %s.\n",filename);

  return(0);
}



int ShowCorners(struct FemType *data,int variable,Real offset)
{
  int i,ind,unknowns;
  Real *solution;

  if(data->nocorners < 1) 
    return(1);

  unknowns = data->edofs[variable];
  if(unknowns == 0) return(2);

  solution = data->dofs[variable];

  printf("Variable %s at free corners:\n",data->dofname[variable]);
  for(i=1;i<=data->nocorners;i++) {
    ind = data->topology[data->corners[2*i-1]][data->corners[2*i]];
    if(data->order[variable][(ind-1)*unknowns+1])
       printf("\t%-2d: %-7.1f  at (%.1f , %.1f )\n",
	      i,solution[(ind-1)*unknowns+1]-offset,
	      data->x[ind],data->y[ind]);
  }
  return(0);
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

    printf("k=%d mini=%d minj=%d ds=%.3le\n",k+1,mini+1,minj+1,minds);
  }

  return(0);
}


int LoadSolutionElmer(struct FemType *data,int results,char *prefix,int info)
/* This procedure reads the solution in a form that is understood 
   by the programs ElmerPost, created
   by Juha Ruokolainen at CSC - IT Center for Science Ltd.. 
   This procedure is not by far general.
   */
{
  int noknots,noelements,novctrs,open;
  int timesteps,i,j,k,grp;
  Real r;
  FILE *in;
  char line[MAXLINESIZE],filename[MAXFILESIZE],text[MAXNAMESIZE];
  int iostat;
  
  AddExtension(prefix,filename,"ep");
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadSolutionElmer: The opening of the Elmer-file %s wasn't successful!\n",
	   filename);
    return(1);
  }
  else 
    printf("Loading Elmer data from %s\n",filename);

  InitializeKnots(data);

  getline;
  sscanf(line,"%d %d %d %d",&noknots,&noelements,&novctrs,&timesteps);

  data->dim = 3;
  data->maxnodes = MAXNODESD2;
  data->noknots = noknots;
  data->noelements = noelements;
  data->timesteps = timesteps;
  
  if(timesteps > 1) 
    printf("LoadSolutionElmer: The subroutine may crash with %d timesteps\n",
	   timesteps);
  if(timesteps < 1) timesteps = 1;
    
  if(info) printf("Allocating for %d knots and %d elements.\n",
		  noknots,noelements);
  AllocateKnots(data);

  if(results) {
    if(timesteps > 1) 
      data->times = Rvector(0,timesteps-1);
    for(i=1;i<=novctrs;i++) {
      sprintf(text,"var%d",i);
      CreateVariable(data,i,timesteps,0.0,text,FALSE);    
    }
  }

  if(info) printf("Reading %d coordinates.\n",noknots);
  for(i=1; i <= noknots; i++) {
    getline;
    sscanf(line,"%le %le %le",
	   &(data->x[i]),&(data->y[i]),&(data->z[i]));
  }

  if(info) printf("Reading %d element topologies.\n",noelements);

  grp = 0;
  open = FALSE;
  for(i=1; i <= noelements; i++) {
    iostat = fscanf(in,"%s",text);
    if(strstr(text,"#group")) {
      grp++;
      printf("Starting a new element group\n");
      iostat = fscanf(in,"%s",text);      
      iostat = fscanf(in,"%s",text);
      open = TRUE;
    }
    if(strstr(text,"#end")) {
      printf("Ending an element group\n");
      iostat = fscanf(in,"%s",text);      
      open = FALSE;
    }
    iostat = fscanf(in,"%d",&(data->elementtypes[i]));
    data->material[i] = grp;
    for(j=0;j< data->elementtypes[i]%100 ;j++) {
      iostat = fscanf(in,"%d",&(data->topology[i][j]));
      data->topology[i][j] += 1;
    }
  }
  if(open) {    
    do {
      iostat = fscanf(in,"%s",text);
    } while (!strstr(text,"#end"));
    iostat = fscanf(in,"%s",text);
    printf("Ending an element group\n");   
    open = FALSE;
  }

  if(results == 0) 
    return(0);

  if(info) printf("Reading %d degrees of freedom for %d knots.\n",
		  novctrs,noknots);
  if (timesteps <= 1) {
    for(i=1; i <= noknots; i++) 
      for(j=1;j <= novctrs;j++) 
	iostat = fscanf(in,"%le",&(data->dofs[j][i]));
  }
  else for(k=0;k<timesteps;k++) {
    iostat = fscanf(in,"%s",text);
    if(iostat < 0) goto end;
    iostat = fscanf(in,"%d",&i);
    iostat = fscanf(in,"%d",&j);
    iostat = fscanf(in,"%le",&r);

    if(0) printf("Loading steps i=%d  j=%d  k=%d  r=%.3lg\n",i,j,k,r);

    for(i=1; i <= noknots; i++) 
      for(j=1;j <= novctrs;j++) 
	iostat = fscanf(in,"%le",&(data->dofs[j][k*noknots+i]));
  }

end:
  data->timesteps = timesteps;

  fclose(in);

  return(0);
}



int SaveSolutionElmer(struct FemType *data,struct BoundaryType *bound,
		      int nobound,char *prefix,int decimals,int info)
/* This procedure saves the solution in a form that is understood 
   by the programs Funcs and ElmerPost, created
   by Juha Ruokolainen at CSC - IT Center for Science Ltd.. 
   */
{
  int material,noknots,noelements,bulkelems,novctrs,sideelems,sideelemtype,elemtype,boundtype;
  char filename[MAXFILESIZE],outstyle[MAXFILESIZE];
  int i,j,k,l,nodesd1,timesteps,nodesd2;
  int ind[MAXNODESD1];
  Real *rpart;
  FILE *out;

  if(!data->created) {
    printf("SaveSolutionElmer: You tried to save points that were never created.\n");
    return(1);
  }

  /* Make a variable showing the owner partition */
  if(data->partitionexist) {
    l = 0;
    do l++; while (data->edofs[l]);
    CreateVariable(data,l,1,0.0,"Partition",FALSE);      
    rpart = data->dofs[l];
    for(i=1;i<=data->noknots;i++) 
      rpart[i] = 1.0 * data->nodepart[i];
  }

  if(data->variables == 0) {
    printf("SaveSolutionElmer: there are no dofs to save!\n");
    return(2);
  }
  
  sideelems = 0;
  if(nobound) {
    for(i=0;i<nobound;i++) {
      if(bound[i].created) sideelems += bound[i].nosides; 
    }
  }

  noknots = data->noknots;
  bulkelems = data->noelements;
  if(nobound)
    noelements = bulkelems + sideelems;
  else
    noelements = bulkelems;
  timesteps = data->timesteps;
  if(timesteps < 1) timesteps = 1;

  novctrs = 0;
  for(i=0;i<MAXDOFS;i++) {
    if(data->edofs[i] == 1) novctrs += 1; 
    if(data->edofs[i] == 2) novctrs += 3; 
    if(data->edofs[i] == 3) novctrs += 3; 
  }

  AddExtension(prefix,filename,"ep");
  if(info) printf("Saving ElmerPost data to %s.\n",filename);  

  out = fopen(filename,"w");
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(3);
  }

  fprintf(out,"%d %d %d %d",noknots,noelements,novctrs,timesteps);

  for(i=0; i<MAXDOFS; i++) {
    if(data->edofs[i] == 1) 
      fprintf(out," scalar: %s",data->dofname[i]);
    else if(data->edofs[i] > 1) 
      fprintf(out," vector: %s",data->dofname[i]);
  }
  fprintf(out,"\n");

  if(info) printf("Saving %d node coordinates.\n",noknots);
  
  sprintf(outstyle,"%%.%dg %%.%dg %%.%dg\n",decimals,decimals,decimals);
  for(i=1; i <= noknots; i++) 
    fprintf(out,outstyle,data->x[i],data->y[i],data->z[i]);      

  printf("Saving %d bulk element topologies.\n",bulkelems);

  for(i=1;i<=bulkelems;i++) {
    elemtype = data->elementtypes[i];
    material = data->material[i];

    if(data->bodynamesexist) 
      fprintf(out,"body_%d_%s %d ",material,data->bodyname[material],elemtype);
    else if(elemtype/100 > 4) 
      fprintf(out,"vol%d %d ",material,elemtype);
    else if(elemtype/100 > 2) 
      fprintf(out,"surf%d %d ",material,elemtype);
    else if(elemtype/100 > 1) 
      fprintf(out,"line%d %d ",material,elemtype);
    else 
      fprintf(out,"pnt%d %d ",material,elemtype);

    nodesd2 = data->elementtypes[i]%100;
    for(j=0;j<nodesd2;j++) 
      fprintf(out,"%d ",data->topology[i][j]-1);
    fprintf(out,"\n");    
  }

  if(nobound) {
    printf("Saving %d side element topologies.\n",sideelems);
    for(j=0;j<nobound;j++) {
      if(bound[j].created == FALSE) continue;
      
      for(i=1;i<=bound[j].nosides;i++) {

	GetBoundaryElement(i,&bound[j],data,ind,&sideelemtype); 

	boundtype = bound[j].types[i];

	if(data->boundarynamesexist) 
	  fprintf(out,"bc_%d_%s %d ",boundtype,data->boundaryname[boundtype],sideelemtype);	  
	else if(sideelemtype/100 > 2) 
	  fprintf(out,"bcside%d %d ",boundtype,sideelemtype);
	else if(sideelemtype/100 > 1) 
	  fprintf(out,"bcline%d %d ",boundtype,sideelemtype);
	else 
	fprintf(out,"bcpnt%d %d ",boundtype,sideelemtype);

	nodesd1 = sideelemtype%100;
	for(k=0;k<nodesd1;k++)
	  fprintf(out,"%d ",ind[k]-1);
	fprintf(out,"\n");
      }
    }
  }

  printf("Saving %d degrees of freedom for each knot.\n",novctrs);
  for(k=0;k<timesteps;k++) {
    for(i=1;i<=noknots;i++){
      for(j=0;j<MAXDOFS;j++) {
	if(data->edofs[j] == 1) 
	  fprintf(out,"%.6lg ",data->dofs[j][k*noknots+i]);
	if(data->edofs[j] == 2) 
	  fprintf(out,"%.6lg %.6lg 0.0 ",
		  data->dofs[j][2*(k*noknots+i)-1],data->dofs[j][2*(k*noknots+i)]);
	if(data->edofs[j] == 3) 
	  fprintf(out,"%.6lg %.6lg %.6lg ",
		  data->dofs[j][3*(k*noknots+i)-2],
		  data->dofs[j][3*(k*noknots+i)-1],
		  data->dofs[j][3*(k*noknots+i)]);
      }
      fprintf(out,"\n");
    }
  }
  fclose(out);

  return(0);
}

int SaveSizeInfo(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info)
{
  int nosides,j;
  FILE *out;
  char filename[MAXFILESIZE];

  if(!data->created) {
    printf("You tried to save points that were never created.\n");
    return(1);
  }

  nosides = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound[j].created == FALSE) continue;
    nosides += bound[j].nosides;
  }

  AddExtension(prefix,filename,"size");
  if(info) printf("Saving size info into file %s.\n",filename);

  out = fopen(filename,"w");
  fprintf(out,"%d\n",data->noknots);
  fprintf(out,"%d\n",data->noelements);
  fprintf(out,"%d\n",nosides);
  fprintf(out,"%d\n",data->nopartitions);

  fclose(out);

  return(0);
}



int SaveElmerInputFemBem(struct FemType *data,struct BoundaryType *bound,
			 char *prefix,int decimals,int info)
/* Saves the mesh in a form that may be used as input 
   in Elmer calculations. Tailored to work with FEM/BEM coupling,
   or other problems with mixed dimension of bulk elements.
   */
{
  int noknots,noelements,material,sumsides,elemtype,fail,nobulkelements,bctype;
  int sideelemtype,nodesd1,nodesd2,newtype,elemdim,maxelemdim,cdstat;
  int i,j,k,l,bulktypes[MAXELEMENTTYPE+1],sidetypes[MAXELEMENTTYPE+1],tottypes;
  int ind[MAXNODESD1],bodyperm[MAXBODIES],bcperm[MAXBCS];
  FILE *out;
  char filename[MAXFILESIZE], outstyle[MAXFILESIZE];
  char directoryname[MAXFILESIZE];

  if(!data->created) {
    printf("You tried to save points that were never created.\n");
    return(1);
  }

  if(data->nodeconnectexist) smallerror("Connectivity data is not saved in the FEM/BEM version");
  if(data->nopartitions > 1) smallerror("Partitioning data is not saved in the FEM/BEM version");

  noelements = data->noelements;
  noknots = data->noknots;
  sumsides = 0;

  for(i=0;i<=MAXELEMENTTYPE;i++)
    bulktypes[i] = sidetypes[i] = 0;

  sprintf(directoryname,"%s",prefix);

  if(info) printf("Saving mesh in ElmerSolver format to directory %s.\n",
		  directoryname);

  cdstat = chdir(directoryname);
  if(cdstat) {
#ifdef MINGW32
    fail = mkdir(directoryname);
#else
    fail = mkdir(directoryname,0750);
#endif
    if(fail) {
      printf("Could not create a result directory!\n");
      return(1);
    }
    else {
      cdstat = chdir(directoryname);
    }
  }
  else {
    printf("Reusing an existing directory\n");
  }

  sprintf(filename,"%s","mesh.nodes");
  out = fopen(filename,"w");

  if(info) printf("Saving %d coordinates to %s.\n",noknots,filename);  
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(2);
  }

  sprintf(outstyle,"%%d %%d %%.%dg %%.%dg %%.%dg\n",decimals,decimals,decimals);
  for(i=1; i <= noknots; i++) 
    fprintf(out,outstyle,i,-1,data->x[i],data->y[i],data->z[i]);    
  fclose(out);


  sprintf(filename,"%s","mesh.elements");
  out = fopen(filename,"w");
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(3);
  }

  maxelemdim = GetMaxElementDimension(data);
  nobulkelements = 0;

  for(i=0;i<MAXBODIES;i++) bodyperm[i] = FALSE;
  for(i=1;i<=noelements;i++) {
    elemtype = data->elementtypes[i];
    elemdim = GetElementDimension(elemtype);
    material = data->material[i];

    if(elemdim == maxelemdim) 
      bodyperm[material] = 1;
    else 
      bodyperm[material] = -1;      
  }
  j = 0;
  k = 0;
  for(i=0;i<MAXBODIES;i++) {
    if(bodyperm[i] > 0) bodyperm[i] = ++j;
    if(bodyperm[i] < 0) bodyperm[i] = --k;
  }


  for(i=1;i<=noelements;i++) {
    elemtype = data->elementtypes[i];

    elemdim = GetElementDimension(elemtype);
    if(elemdim < maxelemdim) continue;

    nobulkelements++;
    material = data->material[i];
    material = bodyperm[material];
    fprintf(out,"%d %d %d",i,material,elemtype);
    
    bulktypes[elemtype] += 1;
    nodesd2 = elemtype%100;
    for(j=0;j < nodesd2;j++) 
      fprintf(out," %d",data->topology[i][j]);
    fprintf(out,"\n");          
  }
  fclose(out);

  if(info) printf("Saving %d (of %d) body elements to mesh.boundary\n",nobulkelements,noelements);


  sprintf(filename,"%s","mesh.boundary");
  out = fopen(filename,"w");
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(4);
  }

  sumsides = 0;
  newtype = 0;
  for(i=0;i<MAXBCS;i++) bcperm[i] = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound[j].created == FALSE) continue;
    for(i=1; i <= bound[j].nosides; i++) 
      bcperm[bound[j].types[i]] = TRUE;
  }
  j = 0;
  for(i=0;i<MAXBCS;i++) 
    if(bcperm[i]) bcperm[i] = ++j;
  newtype = j;


  /* Save normal boundaries */
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound[j].created == FALSE) continue;
    if(bound[j].nosides == 0) continue;
    
    for(i=1; i <= bound[j].nosides; i++) {
      GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],data,ind,&sideelemtype); 
      sumsides++;
      material = bound[j].types[i];
      material = bcperm[material];
      
      fprintf(out,"%d %d %d %d %d",
	      sumsides,material,bound[j].parent[i],bound[j].parent2[i],sideelemtype);
      
      sidetypes[sideelemtype] += 1;

      nodesd1 = sideelemtype % 100;
      for(l=0;l<nodesd1;l++)
	fprintf(out," %d",ind[l]);
      fprintf(out,"\n");
    }
  }

  
  for(i=1;i<=noelements;i++) {
    elemtype = data->elementtypes[i];

    elemdim = GetElementDimension(elemtype);
    if(elemdim == maxelemdim) continue;

    sumsides++;

    material = data->material[i];
    bctype = abs(bodyperm[material]) + newtype;
    
    fprintf(out,"%d %d 0 0 %d",sumsides,bctype,elemtype);
    sidetypes[elemtype] += 1;
    nodesd1 = elemtype%100;
    for(l=0;l<nodesd1;l++)
      fprintf(out," %d",data->topology[i][l]);
    fprintf(out,"\n");
  }

  fclose(out);
  if(info) printf("Saving %d boundary elements to mesh.boundary\n",sumsides);


  tottypes = 0;
  for(i=0;i<=MAXELEMENTTYPE;i++) 
    if(bulktypes[i] || sidetypes[i]) tottypes++;

  sprintf(filename,"%s","mesh.header");
  out = fopen(filename,"w");
  if(info) printf("Saving header info with %d elementtypes to %s.\n",tottypes,filename);  
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(4);
  }
  fprintf(out,"%-6d %-6d %-6d\n",
	  noknots,nobulkelements,sumsides);
  fprintf(out,"%-6d\n",tottypes);
  for(i=0;i<=MAXELEMENTTYPE;i++) 
    if(bulktypes[i] || sidetypes[i]) 
      fprintf(out,"%-6d %-6d\n",i,bulktypes[i]+sidetypes[i]);
  fclose(out);

  if(info) printf("Body and boundary numbers were permutated\n");

  if(data->boundarynamesexist || data->bodynamesexist) {
    sprintf(filename,"%s","mesh.names");
    out = fopen(filename,"w");
    if(info) printf("Saving names info to %s.\n",filename);  
    if(out == NULL) {
      printf("opening of file was not successful\n");
      return(5);
    }
    
    if(data->bodynamesexist) {
      fprintf(out,"! ----- names for bodies -----\n");
      for(i=1;i<MAXBODIES;i++) 
	if(bodyperm[i] > 0) fprintf(out,"$ %s = %d\n",data->bodyname[i],bodyperm[i]);
    }     
    if(data->boundarynamesexist) {
      fprintf(out,"! ----- names for boundaries -----\n");
      /* The reordered numbering of the original BCs */
      for(i=1;i<MAXBCS;i++) {
	if(bcperm[i]) 
	  fprintf(out,"$ %s = %d\n",data->boundaryname[i],bcperm[i]);
      }
      /* Numbering of bodies that are actually saved as boundaries */
      for(i=1;i<MAXBODIES;i++) 
	if(bodyperm[i] < 0) 
	  fprintf(out,"$ %s = %d\n",data->bodyname[i],abs(bodyperm[i])+newtype);
    }
    fclose(out);
  }
  
  cdstat = chdir("..");
  
  return(0);
}


int SolutionFromMeshToMesh(struct CellType *cell1, struct GridType *grid1, 
			   struct FemType *data1,
			   struct CellType *cell2, struct GridType *grid2, 
			   struct FemType *data2,
			   int mapgeo,int variable,int info)
/* Copies variable values from data1 to data2 for two meshes that have similar
   geometry, but different number of elements. Its assumed that all the 
   elements are rectangular in shape. The subroutine may, however, be 
   applied to nonrectangular geometries also, but then the accuracy is 
   difficult to estimate. Note that this subroutine holds only for 
   linear elements.
   */
{
  int xcell,ycell,i1,j1,i2,j2,no1,no2;
  int ind1[MAXNODESD2],ind2,k;
  int nonodes1,nonodes2,unknowns;
  int elem1,elem2,mapres,fast;
  Real coord1[DIM*MAXNODESD2],x2,y2,rx,ry;
  Real epsilon=1.0e-20;
  Real *vector1=NULL,*vector2=NULL;
  
  if(!data1->edofs[variable] || !data2->edofs[variable])
    return(1);

  unknowns = data1->edofs[variable];
  vector1 = data1->dofs[variable];
  vector2 = data2->dofs[variable];

  nonodes1 = grid1->nonodes;
  nonodes2 = grid2->nonodes;

  if(nonodes1 != 4 || nonodes2 != 4) {
    if(info) printf("SolutionFromMeshToMesh: algorithm defined only for 4-node elements\n");
    return(2);
  }

  if(mapgeo && data2->mapgeo < data1->mapgeo) 
    data2->mapgeo = data1->mapgeo;
  else 
    mapgeo = FALSE;

  if(data2->iterdofs[variable] < data1->iterdofs[variable]) {
    mapres = TRUE;
    data2->iterdofs[variable] = data1->iterdofs[variable];
  }
  else
    mapres = FALSE;
      
  if(!mapres && !mapgeo) {
    if(info) 
      printf("SolutionFromMeshToMesh: no mapping for geometry or variable %s (%d vs. %d)!\n",
	     data1->dofname[variable],data1->iterdofs[variable],data2->iterdofs[variable]);
    return(0);
  }

  if(grid1->noknots == data1->noknots && grid2->noknots == data2->noknots) 
    fast = TRUE;
  else
    fast = FALSE;

  for(xcell=1;xcell<=MAXCELLS;xcell++)
    for(ycell=1;ycell<=MAXCELLS;ycell++) 

      /* Go through cells that are common to both grids. */
      if( (no1= grid1->numbered[ycell][xcell]) && (no2= grid2->numbered[ycell][xcell]) ) {

	if(0) printf("xcell=%d  ycell=%d  no1=%d  no=%d\n",xcell,ycell,no1,no2); 

        j1 = 1;
        for(j2=0; j2 <= cell2[no2].yelem; j2++) {
	  i1 = 1;
          for(i2=0; i2 <= cell2[no2].xelem; i2++) {
	    
	    /* ind2 is the original node number for the node */
	    ind2 = GetKnotCoordinate(&(cell2)[no2],i2,j2,&x2,&y2);
	    GetElementCoordinates(&(cell1)[no1],i1,j1,coord1,ind1);

	    /* Find j1 and i1 in the rectangular mesh so that 
	       the node ind2 lies in the element */
	    if(coord1[TOPRIGHT+nonodes1]+epsilon < y2 && j1 < cell1[no1].yelem)
	      do {
	        j1++;
	        GetElementCoordinates(&(cell1)[no1],i1,j1,coord1,ind1);
	      } while(j1 < cell1[no1].yelem  &&  coord1[TOPRIGHT+nonodes1] < y2);
	    
 	    if(coord1[TOPRIGHT]+epsilon < x2 && i1 < cell1[no1].xelem)
	      do {
	        i1++;
	        GetElementCoordinates(&(cell1)[no1],i1,j1,coord1,ind1);
	      } while(i1 < cell1[no1].xelem  &&  coord1[TOPRIGHT] < x2);
	    
	    if(0) printf("j2=%d  i2=%d  j1=%d  i1=%d\n",j2,i2,j1,i1);

	    rx = (coord1[BOTRIGHT]-x2) 
	      / (coord1[BOTRIGHT]-coord1[BOTLEFT]);

	    ry = (coord1[TOPLEFT+nonodes1]-y2) 
	      / (coord1[TOPLEFT+nonodes1]-coord1[BOTLEFT+nonodes1]);

	    rx = MIN(rx,1.0);
	    rx = MAX(rx,0.0);
	    ry = MIN(ry,1.0);
	    ry = MAX(ry,0.0);

	    /* Find the new indices. If there are no discontinuous 
	       boundaries they are the same as the old ones. */
	    if(!fast) {
	      elem1=GetElementIndex(&(cell1)[no1],i1,j1);
	      for(k=0;k<4;k++)
		ind1[k] = data1->topology[elem1][k];

	      if(0) printf("elem1=%d  ind1=[%d %d %d %d]\n",elem1,ind1[0],ind1[1],ind1[2],ind1[3]);
	      
	      if(i2>0 && j2>0) {
		elem2 = GetElementIndex(&(cell2)[no2],i2,j2);
		ind2 = data2->topology[elem2][TOPRIGHT];
	      }
	      else if(i2==0 && j2>0) {
		elem2 = GetElementIndex(&(cell2)[no2],i2+1,j2);
		ind2 = data2->topology[elem2][TOPLEFT];
	      }
	      else if(i2>0 && j2==0) {
		elem2 = GetElementIndex(&(cell2)[no2],i2,j2+1);
		ind2 = data2->topology[elem2][BOTRIGHT];
	      }
	      else {
		elem2 = GetElementIndex(&(cell2)[no2],i2+1,j2+1);
		ind2 = data2->topology[elem2][BOTLEFT];
	      }
	    }

	    if(mapgeo) {
	      data2->x[ind2] 
		= data1->x[ind1[BOTLEFT]] * rx * ry
		+ data1->x[ind1[BOTRIGHT]] * (1.-rx) * ry
		+ data1->x[ind1[TOPLEFT]] * rx * (1.-ry)
		+ data1->x[ind1[TOPRIGHT]] * (1.-rx) * (1.-ry);
	      data2->y[ind2] 
		= data1->y[ind1[BOTLEFT]] * rx * ry
		+ data1->y[ind1[BOTRIGHT]] * (1.-rx) * ry
		+ data1->y[ind1[TOPLEFT]] * rx * (1.-ry)
		+ data1->y[ind1[TOPRIGHT]] * (1.-rx) * (1.-ry);
	    }
	    
	    if(mapres) 
	      for(k=1;k<=unknowns;k++) {
		vector2[unknowns*(ind2-1)+k] 
		  = vector1[unknowns*(ind1[BOTLEFT]-1)+k] * rx * ry
		  + vector1[unknowns*(ind1[BOTRIGHT]-1)+k]* (1.-rx) * ry
		  + vector1[unknowns*(ind1[TOPLEFT]-1)+k] * rx * (1.-ry)
		  + vector1[unknowns*(ind1[TOPRIGHT]-1)+k]* (1.-rx) * (1.-ry);
	      }
	  }
        }
      }
  
  if(info) {
    if(mapgeo) 
      printf("Geometry was mapped from one mesh to another!\n");
    if(mapres) 
      printf("Results of %s was mapped from one mesh to another!\n",
	     data1->dofname[variable]);
  }  

  return(0);
}

void InspectVector(Real *vector,int first,int last,Real *min,
		    Real *max,int *mini,int *maxi)
/* Returns the value and position of the smallest and largest element
   in vector[first,last]. 
   */
{
  int i;
  *min  = vector[first];
  *mini = first;
  *max  = vector[first];
  *maxi = first;

  for(i=first+1;i<=last;i++) {
    if (vector[i]>*max) {
      *max=vector[i];
      *maxi=i;
    }
    if (vector[i]<*min) {
      *min=vector[i];
      *mini=i;
    }
  }
}


int Steepest(Real *vector,int first,int last)
/* Finds the position where vector is at its steepest */
{
  int i,steep;
  Real aid,sub=0.0;

  steep = -1;
  for(i=first;i<last;i++) {
    aid=fabs(vector[i+1]-vector[i]);
    if ( aid > sub) {
      sub=aid;
      steep=i;
    }
  }
  return(steep);
}


Real MeanVector(Real *vector,int first,int last)
/* Calculates the mean of vector[first,last] */
{
  Real sum=0.0;
  int i;

  for(i=first;i<=last;i++)
    sum += vector[i];

  return(sum/(last-first+1));
}



Real AbsMeanVector(Real *vector,int first,int last)
/* Calculates the absolute mean of vector[first,last] */
{
  Real sum=0.0;
  int i;

  for(i=first;i<=last;i++)
    sum += fabs(vector[i]);

  return(sum/(last-first+1));
}


Real DifferVector(Real *vector1,Real *vector2,int first,int last)
/* Calcultes the mean of the relative difference of two vectors */
{
  Real sum=0.0, eps=1.0E-50;
  int i,n;

  for (i=first;i<=last;i++) {
      if ( fabs(vector1[i]+vector2[i]) > eps)
      sum += fabs(2*(vector1[i]-vector2[i])/(vector1[i]+vector2[i]));
    }
  n=last-first+1;

  return ( sum/(Real)(n) );
}


void ReformVector(Real *vector1,int n1,Real *vector2,int n2)
/* Adjusts the values of a vector to another vector with a different number
   of elements 
   */
{
  int i1,i2;
  Real x1,d1;

  for(i2=0;i2<n2;i2++) {
    x1=(n1)*((Real)(i2)/(Real)(n2));
    i1=(int)(x1);
    d1=x1-(Real)(i1);
    vector2[i2]=d1*vector1[i1+1]+(1.0-d1)*vector1[i1];
  }
  vector2[n2]=vector1[n1];
}



void AdjustVector(Real max,Real min,Real *vector,int first,int last)
/* Scales the values of a vector to range [min,max] using a linear model. */
{
  int i;
  Real oldmax,oldmin;

  oldmax=vector[first];
  oldmin=oldmax;
  for(i=first+1;i<=last;i++) {
    if (oldmax < vector[i]) 
      oldmax=vector[i];
    if (oldmin > vector[i])
      oldmin=vector[i];
  }
  for(i=first;i<=last;i++)
    vector[i]=(max-min)*(vector[i]-oldmin)/(oldmax-oldmin)+min;
}



int ReadRealVector(Real *vector,int first,int last,char *filename)
/* Reads a Real vector from an ascii-file with a given name. */
{
  int i;
  FILE *in;
  Real num;

  if ((in = fopen(filename,"r")) == NULL) {
    printf("The opening of the real vector file '%s' wasn't succesfull !\n",filename);
    return(1);
  }
  for(i=first;i<=last;i++) {
    fscanf(in,"%le\n",&num);
    vector[i]=num;
    }
  fclose(in);

  return(0);
}



void SaveRealVector(Real *vector,int first,int last,char *filename)
/* Saves an Real vector to an ascii-file with a given name. */
{
  int i;
  FILE *out;

  out = fopen(filename,"w");
  for (i=first;i<=last;i++) {
    fprintf(out,"%.6le",vector[i]);
    fprintf(out,"\n");
    }
  fclose(out);
}



int ReadIntegerVector(int *vector,int first,int last,char *filename)
{
  int i;
  FILE *in;
  int num;

  if ((in = fopen(filename,"r")) == NULL) {
    printf("The opening of the int vector file '%s' wasn't succesfull !\n",filename);
    return(1);
  }
  for(i=first;i<=last;i++) {
    fscanf(in,"%d\n",&num);
    vector[i]=num;
    }
  fclose(in);

  return(0);
}



void SaveIntegerVector(int *vector,int first,int last,char *filename)
{
  int i;
  FILE *out;

  out = fopen(filename,"w");
  for (i=first;i<=last;i++) {
    fprintf(out,"%d",vector[i]);
    fprintf(out,"\n");
    }
  fclose(out);
}



int ReadRealMatrix(Real **matrix,int row_first,int row_last,
		int col_first,int col_last,char *filename)
{
  int i,j;
  FILE *in;
  Real num;

  if ((in = fopen(filename,"r")) == NULL) {
    printf("The opening of the real matrix file '%s' wasn't succesfull!\n",filename);
    return(1);
  }

  for(j=row_first;j<=row_last;j++) {
    for(i=col_first;i<=col_last;i++) {
      fscanf(in,"%le\n",&num);
      matrix[j][i]=num;
    }
  }
  fclose(in);

  return(0);
}



void SaveRealMatrix(Real **matrix,int row_first,int row_last,
		    int col_first,int col_last,char *filename)
{
  int i,j;
  FILE *out;

  out = fopen(filename,"w");
  for (j=row_first;j<=row_last;j++) {
    for (i=col_first;i<=col_last;i++) {
      fprintf(out,"%-14.6lg",matrix[j][i]);
      fprintf(out,"\t");
    }
    fprintf(out,"\n");
  }
  fclose(out);
}



int ReadIntegerMatrix(int **matrix,int row_first,int row_last,
		int col_first,int col_last,char *filename)
{
  int i,j;
  FILE *in;
  int num;

  if ((in = fopen(filename,"r")) == NULL) {
    printf("The opening of the int matrix file '%s' wasn't succesfull!\n",filename);
    return(FALSE);
  }

  for(j=row_first;j<=row_last;j++) {
    for(i=col_first;i<=col_last;i++) {
      fscanf(in,"%d\n",&num);
      matrix[j][i]=num;
    }
  }
  fclose(in);

  return(TRUE);
}



void SaveIntegerMatrix(int **matrix,int row_first,int row_last,
		int col_first,int col_last,char *filename)
{
  int i,j;
  FILE *out;

  out = fopen(filename,"w");
  for (j=row_first;j<=row_last;j++) {
    for (i=col_first;i<=col_last;i++) {
      fprintf(out,"%-8d",matrix[j][i]);
      fprintf(out,"\t");
    }
    fprintf(out,"\n");
  }
  fclose(out);
}


void SaveNonZeros(Real **matrix,int row_first,int row_last,
		  int col_first,int col_last,char *filename)
/* Saves the nonzero elements in an file. */
{
  int i,j;
  FILE *out;
  Real nearzero=1.0e-20;

  out = fopen(filename,"w");
  for (j=row_first;j<=row_last;j++) 
    for (i=col_first;i<=col_last;i++) 
      if (fabs(matrix[j][i]) > nearzero) 
        fprintf(out,"%d\t %d\t %-12.6le\n",j,i,matrix[j][i]);
  
  fclose(out);
}
