/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter Råback
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

/* --------------------:  femfilein.c  :-------------------------- */

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
/*#include <unistd.h>*/

#include "nrutil.h"
#include "common.h"
#include "femdef.h"
#include "femtypes.h"
#include "femknot.h"
#include "femfileout.h"


int SaveAbaqusInput(struct FemType *data,char *prefix,int info)
/* Saves the grid in a format that can be read by ABAQUS 
   program designed for sructural mechanics. 
   The elementtype is set to be that of thermal conduction. 
   */
{
  int noknots,noelements,nonodes;
  char filename[MAXFILESIZE];
  int i,j;
  FILE *out;

  noknots = data->noknots;
  noelements = data->noelements;
  nonodes = data->maxnodes;

  if(nonodes != 4 && nonodes != 8) {
    printf("SaveAbaqusInput: not designed for %d-node elements\n",nonodes);
    return(1);
  }

  AddExtension(prefix,filename,"inp");
  out = fopen(filename,"w");

  if(info) printf("Saving ABAQUS data to %s.\n",filename);  

  fprintf(out,"*HEADING\n");
  fprintf(out,"Abaqus input file creator by Peter.Raback@csc.fi\n");

  fprintf(out,"*NODE, SYSTEM=");
  if(data->coordsystem == COORD_CART2) fprintf(out,"R\n");
  else if(data->coordsystem == COORD_AXIS) fprintf(out,"C\n");
  else if(data->coordsystem == COORD_POLAR) fprintf(out,"P\n");

  for(i=1; i <= noknots; i++) 
    fprintf(out,"%8d, %12.4e, %12.4e, 0.0\n",i,data->x[i],data->y[i]);

  fprintf(out,"*ELEMENT,TYPE=");
  if(nonodes == 4) fprintf(out,"DC2D4 ");
  else if(nonodes == 8) fprintf(out,"DC2D8 ");
  else printf("SaveAbaqusInput: Not defined for %d-node elements\n",nonodes);

  fprintf(out,",ELSET=SETTI1\n");
  for(i=1;i<=noelements;i++) {
    fprintf(out,"%8d, ",i);
    for(j=1;j<nonodes;j++) 
      fprintf(out,"%6d, ",data->topology[i][j-1]);
    fprintf(out,"%6d\n",data->topology[i][nonodes-1]);
  }

  fprintf(out,"*SOLID SECTION, ELSET=SETTI1, MATERIAL=MAT1\n");

  fprintf(out,"*MATERIAL, NAME=MAT1\n");
  fprintf(out,"*CONDUCTIVITY\n");
  fprintf(out,"1.0\n");

  fprintf(out,"*RESTART, WRITE, FREQUENCY=1\n");
  fprintf(out,"*FILE FORMAT, ASCII\n");

  fprintf(out,"*STEP\n");
  fprintf(out,"*HEAT TRANSFER, STEADY STATE\n");
  fprintf(out,"*END STEP\n");

  fclose(out);

  if(info) printf("Wrote the mesh in ABAQUS restart format.\n");

  return(0);
}


int SaveFidapOutput(struct FemType *data,char *prefix,int info,
		    int vctrs,Real *vect1, ...)
/* This procedure saves the solution in a form that is understood by 
   the CFD program Fidap. The routine that reads this file seems
   to be stupidly place-dependent. The best manual for this 
   subroutine is provided in the appendix E of FIDAP users manual.
   */
{
  int noknots,noelements,dim,nonodes;
  int i,j,no,nogroup,material,cellelem;
  int elemcodes[MAT_MAXNUMBER+1],mat[MAT_MAXNUMBER+1];
  FILE *out;
  Real *dofs[MAXDOFS];
  char filename[MAXFILESIZE];
  va_list ap;

  if(!data->created) {
    printf("You tried to save points that were never created.\n");
    return(1);
  }

  printf("Saving results in FIDAP neutral format.\n");

  noknots = data->noknots;
  noelements = data->noelements;
  dim = data->dim;

  if(vctrs > 3) {
    printf("SaveSolutionFidap: Maximum of 3 d.o.f.\n");
    vctrs = 3;
  }

  /* Read the list of pointers to vectors to be saved. */
  if(vctrs > 0) {
    va_start(ap,vect1);
    dofs[1] = vect1;
    for(i=2;i<=vctrs;i++)
      dofs[i] = va_arg(ap,Real*);
    va_end(ap);
  }

  /* Create groups by placing the same materials to same groups. */
  nogroup = 0;
  for(i=1;i<=MAT_MAXNUMBER;i++) 
    elemcodes[i] = mat[i] = 0;
  for(j=1;j <= noelements;j++)  {
    material = data->material[j];
    if(material > 0 && material<=MAT_MAXNUMBER) mat[material] += 1;
    elemcodes[material] = data->elementtypes[j];
  } 
  for(i=1;i<=MAT_MAXNUMBER;i++) 
    if(mat[i] > 0) nogroup++;

  AddExtension(prefix,filename,"fidap");
  out = fopen(filename,"w");

  /* Control information */
  fprintf(out,"** FIDAP NEUTRAL FILE\n");
  fprintf(out,"Fidap input file creator by Peter.Raback@csc.fi\n");
  fprintf(out,"VERSION %7.2f\n",7.52); /* Fidap version */
  fprintf(out," 1 Dec 96    12:00:00\n");
  fprintf(out,"%15s%15s%15s%15s%15s\n","NO. OF NODES",
	  "NO. ELEMENTS","NO. ELT GROUPS","NDFCD","NDFVL");
  fprintf(out,"%15d%15d%15d%15d%15d\n",noknots,noelements,nogroup,dim,dim);
  fprintf(out,"%15s%15s%15s%15s\n",
	  "STEADY/TRANS","TURB. FLAG","FREE SURF FLAG","COMPR. FLAG");
  fprintf(out,"%15d%15d%15d%15d\n",0,0,0,0);
  fprintf(out,"%s\n","TEMPERATURE/SPECIES FLAGS");
  fprintf(out,"%s\n"," 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
  fprintf(out,"%s\n","PRESSURE FLAGS - IDCTS, IPENY MPDF");
  fprintf(out,"%10d%10d%10d\n",0,0,1);

  fprintf(out,"NODAL COORDINATES\n");
  for(i=1; i <= noknots; i++) 
    fprintf(out,"%10d%20.10e%20.10e\n",i,data->y[i],data->x[i]);

  /* Boundary and initial conditions */
  fprintf(out,"BOUNDARY CONDITIONS\n");
  fprintf(out,"%10d%10d%10d%20.10e\n",0,0,5,0.0);

  fprintf(out,"ELEMENT GROUPS\n");
  nogroup = 0;
  for(no=1;no<=MAT_MAXNUMBER;no++)
    if(mat[no] > 0) {
      nonodes = elemcodes[no]%100;
      nogroup++;
      cellelem = mat[no];
      fprintf(out,"GROUP:    %5d ELEMENTS:%10d NODES:   %10d GEOMETRY:%5d TYPE:%4d\n",
	      nogroup,cellelem,nonodes,1,1);
      fprintf(out,"ENTITY NAME:   %s%d\n","material",nogroup);

      for(j=1;j <= noelements;j++)  {
	if(data->material[j] == no) {    
	  fprintf(out,"%8d\n",j);
	  for(i=0;i < nonodes;i++)
	    fprintf(out,"%8d",data->topology[j][i]);
	  fprintf(out,"\n");
	}
      } 
    }    

  fprintf(out,"TIMESTEP: %5d TIME:     %15.7e INCRMNT: %15.7e\n",1,1.0,1.0);

  fprintf(out,"VELOCITY\n");            
  if(vctrs < 2) 
    for(i=1;i<=2*noknots;i++) {
      fprintf(out,"%16.9e",0.0);
      if(i%5==0) fprintf(out,"\n");
    }  
  else 
    for(i=1;i<=2*noknots;i++) {
      fprintf(out,"%16.9e",dofs[1][i]);
      if((2*i-1)%5 == 0) fprintf(out,"\n");
      fprintf(out,"%16.9e",dofs[2][i]);
      if((2*i)%5 == 0) fprintf(out,"\n");   
  }
  if((2*noknots)%5 != 0) fprintf(out,"\n");

  fprintf(out,"TEMPERATURE\n");

  if(vctrs == 2) 
    for(i=1;i<=noknots;i++) {
      fprintf(out,"%16.9e",0.0);
      if(i%5==0) fprintf(out,"\n");
    }
  else 
    for(i=1;i<=noknots;i++) {
      fprintf(out,"%16.9e",dofs[vctrs][i]);
      if(i%5==0) fprintf(out,"\n");
    }
  if(noknots%5 != 0) fprintf(out,"\n");
  fprintf(out,"ENDOFTIMESTEP\n");       

  fclose(out);
  if(info) printf("Results were saved in FIDAP neutral format to file %s.\n",filename);

  return(0);
}



static int ElmerToGmshType(int elmertype)
{
  int gmshtype = 0;

  switch (elmertype) {
      
  case 202:        
    gmshtype = 1;
    break;
  case 303:        
    gmshtype = 2;
    break;
  case 404:        
    gmshtype = 3;
    break;
  case 504:        
    gmshtype = 4;
    break;
  case 808:        
    gmshtype = 5;
    break;
  case 706:        
    gmshtype = 6;
    break;
  case 605:        
    gmshtype = 7;
    break;
  case 203:        
    gmshtype = 8;
    break;
  case 306:        
    gmshtype = 9;
    break;
  case 409:        
    gmshtype = 10;
    break;
  case 510:        
    gmshtype = 11;
    break;
  case 827:        
    gmshtype = 12;
    break;
  case 101:        
    gmshtype = 15;
    break;
  case 408:
    gmshtype = 16;
    break;
  case 820:
    gmshtype = 17;
    break;
  case 715:
    gmshtype = 18;
    break;
  case 613:
    gmshtype = 19;
    break;
  case 310:
    gmshtype = 21;
    break;

  default:
    printf("Elmer element %d does not have an Gmsh counterpart!\n",elmertype);
  }

  return(gmshtype);
}




static void ElmerToGmshIndx(int elemtype,int *topology)
{
  int i=0,nodes=0,oldtopology[MAXNODESD2];
  int reorder, *porder;

  int order510[]={0,1,2,3,4,5,6,7,9,8};
  int order613[]={0,1,2,3,4,5,8,10,6,7,9,11,12};
  int order715[]={0,1,2,3,4,5,6,9,7,8,10,11,12,14,13};
  int order820[]={0,1,2,3,4,5,6,7,8,11,12,9,10,12,14,15,16,18,19,17};


  reorder = FALSE;

  switch (elemtype) {
      
  case 510:        
    reorder = TRUE;
    porder = &order510[0];
    break;

  case 613:        
    reorder = TRUE;
    porder = &order613[0];
    break;

  case 715:        
    reorder = TRUE;
    porder = &order715[0];
    break;

  case 820:        
    reorder = TRUE;
    porder = &order820[0];
    break;

  }

  if( reorder ) {
    nodes = elemtype % 100;
    for(i=0;i<nodes;i++) 
      oldtopology[i] = topology[i];
    for(i=0;i<nodes;i++) 
      topology[porder[i]] = oldtopology[i];
  }
}



int SaveMeshGmsh(struct FemType *data,struct BoundaryType *bound,
		 int nobound,char *prefix,int decimals,int info)
/* This procedure saves the mesh in a format understood by Gmsh */
{
  int material,noknots,noelements,bulkelems,sideelems,gmshtype,elemtype,boundtype;
  char filename[MAXFILESIZE],outstyle[MAXFILESIZE];
  int i,j,k,nodesd2,elemind;
  int ind[MAXNODESD2];
  FILE *out;

  if(!data->created) {
    printf("SaveMeshGmsh: You tried to save points that were never created.\n");
    return(1);
  }

  /* Compute the number of boundary elements and register the minimum BC type */
  sideelems = 0;
  if(nobound) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created) 
	sideelems += bound[j].nosides; 
    }
    if(0) printf("number of boundary elements: %d\n",sideelems);
  }
 
  noknots = data->noknots;
  bulkelems = data->noelements;
  if(nobound)
    noelements = bulkelems + sideelems;
  else
    noelements = bulkelems;

  AddExtension(prefix,filename,"msh");
  if(info) printf("Saving Gmsh mesh to %s.\n",filename);  

  out = fopen(filename,"w");
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(3);
  }

  fprintf(out,"$MeshFormat\n");
  fprintf(out,"2.2 0 %ld\n",sizeof(double));
  fprintf(out,"$EndMeshFormat\n");

  if(info) printf("Saving %d node coordinates.\n",noknots);
  fprintf(out,"$Nodes\n");
  fprintf(out,"%d\n",noknots);
  if(data->dim == 1) {
    sprintf(outstyle,"%%d %%.%dg 0 0\n",decimals);
    for(i=1; i <= noknots; i++) 
      fprintf(out,outstyle,i,data->x[i]);
  }
  else if(data->dim == 2) {
    sprintf(outstyle,"%%d %%.%dg %%.%dg 0\n",decimals,decimals);
    for(i=1; i <= noknots; i++) 
      fprintf(out,outstyle,i,data->x[i],data->y[i]);
  }
  else if(data->dim == 3) {
    sprintf(outstyle,"%%d %%.%dg %%.%dg %%.%dg\n",decimals,decimals,decimals);
    for(i=1; i <= noknots; i++) 
      fprintf(out,outstyle,i,data->x[i],data->y[i],data->z[i]);      
  }
  fprintf(out,"$EndNodes\n");
  

  printf("Saving %d element topologies.\n",bulkelems);

  fprintf(out,"$Elements\n");
  fprintf(out,"%d\n",bulkelems+sideelems);
  for(i=1;i<=bulkelems;i++) {
    elemtype = data->elementtypes[i];
    material = data->material[i];

    gmshtype = ElmerToGmshType( elemtype );
    
    fprintf(out,"%d %d %d %d %d",i,gmshtype,2,0,material);

    nodesd2 = data->elementtypes[i]%100;

    for(j=0;j<nodesd2;j++) 
      ind[j] = data->topology[i][j];

    ElmerToGmshIndx(elemtype,ind);

    for(j=0;j<nodesd2;j++) 
      fprintf(out," %d",ind[j]);
    fprintf(out,"\n");    
  }

  elemind = bulkelems;
  if(nobound) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created == FALSE) continue;
      
      for(i=1;i<=bound[j].nosides;i++) {

	GetBoundaryElement(i,&bound[j],data,ind,&elemtype); 
	boundtype = bound[j].types[i];

	gmshtype = ElmerToGmshType( elemtype );
	elemind += 1;

	fprintf(out,"%d %d %d %d %d",elemind,gmshtype,2,0,boundtype);
	nodesd2 = elemtype%100;
	
	ElmerToGmshIndx(elemtype,ind);

	for(k=0;k<nodesd2;k++) 
	  fprintf(out," %d",ind[k]);
	fprintf(out,"\n");	
      }
    }
  }
  fprintf(out,"$EndElements\n");


#if 0
  if(data->bodynamesexist || data->boundarynamesexist)  {
    fprintf(out,"$PhysicalNames\n");
    

    /* this is not really coded yet, just some bits to continue from */

    $PhysicalNames
     number-of-names
     physical-dimension physical-number "physical-name"
     ...

      if(data->boundarynamesexist) 
	fprintf(out,"bc_%d_%s %d ",boundtype,data->boundaryname[boundtype],sideelemtype);	  
    
    if(data->bodynamesexist) 
      fprintf(out,"body_%d_%s %d ",material,data->bodyname[material],elemtype);
      fprintf(out,"$EndPhysicalNames\n");
  }

  timesteps = data->timesteps;
  if(timesteps < 1) timesteps = 1;

  if(data->variables == 0) {
    printf("SaveMeshGmsh: there are no dofs to save!\n");
    return(2);
  }
 

  novctrs = 0;
  for(i=0;i<MAXDOFS;i++) {
    if(data->edofs[i] == 1) novctrs += 1; 
    if(data->edofs[i] == 2) novctrs += 3; 
    if(data->edofs[i] == 3) novctrs += 3; 
  }

  if(data->partitionexist) {
    l = 0;
    do l++; while (data->edofs[l]);
    CreateVariable(data,l,1,0.0,"Partition",FALSE);      
    rpart = data->dofs[l];
    for(i=1;i<=data->noknots;i++) 
      rpart[i] = 1.0 * data->nodepart[i];
  }



  for(i=0; i<MAXDOFS; i++) {
    if(data->edofs[i] == 1) 
      fprintf(out," scalar: %s",data->dofname[i]);
    else if(data->edofs[i] > 1) 
      fprintf(out," vector: %s",data->dofname[i]);
  }
  fprintf(out,"\n");

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
#endif

  fclose(out);

  printf("SaveMeshGmsh: All done\n");

  return(0);
}


static int ElmerToVtkType(int elmertype)
{
  int vtktype = 0;

  switch (elmertype) {
      
  case 101:
    vtktype = 1;
  case 202:        
    vtktype = 3;
    break;
  case 203:        
    vtktype = 21;
    break;
  case 303:        
    vtktype = 5;
    break;
  case 306:        
    vtktype = 22;
    break;
  case 404:        
    vtktype = 9;
    break;
  case 408:        
    vtktype = 23;
    break;
  case 409:        
    vtktype = 28;
    break;
  case 504:        
    vtktype = 10;
    break;
  case 510:        
    vtktype = 24;
    break;
  case 605:        
    vtktype = 14;
    break;
  case 613:        
    vtktype = 27;
    break;
  case 706:        
    vtktype = 13;
    break;
  case 715:        
    vtktype = 26;
    break;
  case 808:        
    vtktype = 12;
    break;
  case 820:
    vtktype = 25;
    break;
  case 827:
    vtktype = 29;
    break;

  default:
    printf("Elmer element %d does not have an Vtk counterpart!\n",elmertype);
  }

  return(vtktype);
}



static void ElmerToVtkIndx(int elemtype,int *topology)
{
  int i=0,nodes=0,oldtopology[MAXNODESD2];
  int reorder, *porder;

  int order820[]={0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
  int order827[]={0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,23,21,20,22,24,25,26};

  reorder = FALSE;

  switch (elemtype) {
    
  case 820:        
    reorder = TRUE;
    porder = &order820[0];
    break;
    
  case 827:        
    reorder = TRUE;
    porder = &order827[0];
    break;
  }

  if( reorder ) {
    nodes = elemtype % 100;
    for(i=0;i<nodes;i++) 
      oldtopology[i] = topology[i];
    for(i=0;i<nodes;i++) 
      topology[i] = oldtopology[porder[i]];
  }
}



int SaveMeshVtu(struct FemType *data,struct BoundaryType *bound,
		 int nobound,char *prefix,int dummyzero,int info)
/* This procedure saves the mesh in the VTU format understood by Paraview and Visit, for example. */
{
  int material,noknots,noelements,bulkelems,sideelems,vtktype,elemtype,boundtype;
  char filename[MAXFILESIZE],outstyle[MAXFILESIZE];
  int i,j,k,nodesd2,elemind,idoffset,di;
  int ind[MAXNODESD2];
  int LittleEnd,PrecBits,elemoffset;
  FILE *out;

  /* If we create dummy zero node the real indexes start from one */
  if( dummyzero ) 
    di = 0;
  else
    di = 1;
  
  
  if(!data->created) {
    printf("SaveMeshVtk: You tried to save points that were never created.\n");
    return(1);
  }

  /* Compute the number of boundary elements and register the minimum BC type */
  sideelems = 0;
  if(nobound) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created) 
	sideelems += bound[j].nosides; 
    }
    if(0) printf("number of boundary elements: %d\n",sideelems);
  }
 
  noknots = data->noknots;
  bulkelems = data->noelements;
  if(nobound)
    noelements = bulkelems + sideelems;
  else
    noelements = bulkelems;

  AddExtension(prefix,filename,"vtu");
  if(info) printf("Saving VTU mesh to %s.\n",filename);  

  out = fopen(filename,"w");
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(3);
  }

  idoffset = 100;
  LittleEnd = FALSE;
  PrecBits = 64; /* 32 for single precision */
  

  if(info) printf("Saving Elmer mesh in ascii VTU format\n");
  fprintf(out,"<?xml version=\"1.0\"?>\n");
  if( LittleEnd ) 
    fprintf(out,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  else
    fprintf(out,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");

  fprintf(out,"  <UnstructuredGrid>\n");
  fprintf(out,"    <Piece NumberOfPoints=\"%d\"  NumberOfCells=\"%d\">\n",noknots+dummyzero,noelements);

  /* Write out the nodal indexes, this is mainly just on example */
  fprintf(out,"      <PointData>\n");
  if(0) { /* Here we as as floats - just for testing */
    fprintf(out,"        <DataArray type=\"Float%d\" Name=\"NodeNumber\" NumberOfComponents=\"1\" format=\"ascii\">\n",PrecBits);
    if( dummyzero ) fprintf(out,"%12.6le ",0.0);
    for(i=1;i<=noknots;i++) fprintf(out,"%12.6le ",1.0*i);
  }
  else { /* And here as integers - what they really are */
    fprintf(out,"        <DataArray type=\"Int32\" Name=\"NodeNumber\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    if( dummyzero ) fprintf(out,"%d ",0);
    for(i=1;i<=noknots;i++) fprintf(out,"%d ",i);
  }
  fprintf(out,"\n");
  fprintf(out,"        </DataArray>\n");
  fprintf(out,"      </PointData>\n");


  printf("Saving cell data (Element numbers and Geometry Ids).\n");
  fprintf(out,"      <CellData>\n");
  /* Write out the element indexes, this is mainly just on example */
  if(0) { /* This as floats - just for testing */
    fprintf(out,"        <DataArray type=\"Float%d\" Name=\"ElementNumber\" NumberOfComponents=\"1\" format=\"ascii\">\n",PrecBits);
    for(i=1;i<=noelements;i++) 
      fprintf(out,"%12.6le ",1.0*i);
    fprintf(out,"\n");
    fprintf(out,"        </DataArray>\n");
  }
  else { /* This as integers - what they are */
    fprintf(out,"        <DataArray type=\"Int32\" Name=\"ElementNumber\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(i=1;i<=noelements;i++) 
      fprintf(out,"%d ",i);
    fprintf(out,"\n");
    fprintf(out,"        </DataArray>\n");
  }
    
  /* Write out the geometry Ids */
  fprintf(out,"        <DataArray type=\"Int32\" Name=\"GeometryIds\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(i=1;i<=bulkelems;i++) 
    fprintf(out,"%d ",data->material[i]);

  if(nobound ) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created == FALSE) continue;      
      for(i=1;i<=bound[j].nosides;i++) {
	boundtype = bound[j].types[i];
	fprintf(out,"%d ",idoffset + bound[j].types[i]);
      }
    }
  }
  fprintf(out,"\n");
  fprintf(out,"        </DataArray>\n");


  /* Write out the geometry Ids */
  if(data->partitionexist) {
    fprintf(out,"        <DataArray type=\"Int32\" Name=\"Partition\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(i=1;i<=bulkelems;i++) 
      fprintf(out,"%d ",data->elempart[i]);
    
    if(nobound ) {
      for(j=0;j<nobound;j++) {
	if(bound[j].created == FALSE) continue;      
	for(i=1;i<=bound[j].nosides;i++) {
	  k = bound[j].parent[i];
	  if(k) 
	    fprintf(out,"%d ",data->elempart[k]);
	  else
	    fprintf(out,"%d ",0);
	}
      }
    }
    fprintf(out,"\n");
    fprintf(out,"        </DataArray>\n");
  }

  fprintf(out,"      </CellData>\n");


  if(info) printf("Saving %d nodal coordinates\n",noknots);
  fprintf(out,"      <Points>\n");  
  fprintf(out,"        <DataArray type=\"Float%d\" Name=\"ElementNumber\" NumberOfComponents=\"3\" format=\"ascii\">\n",PrecBits);
  if( dummyzero ) {    
    if(info) printf("Creating dummy duplicate for first index to allow numbering from one\n");
    fprintf(out,"%12.6le %12.6le %12.6le ",data->x[1],data->y[1],data->z[1]);  
  }
  for(i=1;i<=noknots;i++) 
    fprintf(out,"%12.6le %12.6le %12.6le ",data->x[i],data->y[i],data->z[i]);
  fprintf(out,"\n");
  fprintf(out,"        </DataArray>\n");
  fprintf(out,"      </Points>\n");


  printf("Saving %d element topologies.\n",noelements);
  fprintf(out,"      <Cells>\n");  


  fprintf(out,"        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(i=1;i<=bulkelems;i++) {
    elemtype = data->elementtypes[i];
    for(j=0;j<elemtype%100;j++)
      ind[j] = data->topology[i][j];
    ElmerToVtkIndx( elemtype, ind );
    for(j=0;j<elemtype%100;j++)
      fprintf(out,"%d ",ind[j]-di);
  }
  if(nobound ) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created == FALSE) continue;      
      for(i=1;i<=bound[j].nosides;i++) {
	GetBoundaryElement(i,&bound[j],data,ind,&elemtype); 
	for(k=0;k<elemtype%100;k++)
	  fprintf(out,"%d ",ind[k]-di);
      }
    }
  }
  fprintf(out,"\n");
  fprintf(out,"        </DataArray>\n");


  fprintf(out,"        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  elemoffset = 0;
  for(i=1;i<=bulkelems;i++) {
    elemtype = data->elementtypes[i];
    elemoffset += elemtype % 100;
    fprintf(out,"%d ",elemoffset );
  }
  if(nobound ) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created == FALSE) continue;      
      for(i=1;i<=bound[j].nosides;i++) {
	GetBoundaryElement(i,&bound[j],data,ind,&elemtype); 
	elemoffset += elemtype % 100;
	fprintf(out,"%d ",elemoffset );
      }
    }
  }
  fprintf(out,"\n");
  fprintf(out,"        </DataArray>\n");


  fprintf(out,"        <DataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(i=1;i<=bulkelems;i++) {
    elemtype = data->elementtypes[i];
    vtktype = ElmerToVtkType( elemtype );
    fprintf(out,"%d ",vtktype );
  }
  if(nobound ) {
    for(j=0;j<nobound;j++) {
      if(bound[j].created == FALSE) continue;      
      for(i=1;i<=bound[j].nosides;i++) {
	GetBoundaryElement(i,&bound[j],data,ind,&elemtype); 
	vtktype = ElmerToVtkType( elemtype );
	fprintf(out,"%d ",vtktype );
      }
    }
  }
  fprintf(out,"\n");
  fprintf(out,"        </DataArray>\n");

 
  fprintf(out,"      </Cells>\n");
  fprintf(out,"    </Piece>\n");
  fprintf(out,"  </UnstructuredGrid>\n");
  fprintf(out,"</VTKFile>\n");

  fclose(out);

  if(info) printf("Saving of Elmer mesh to VTK format finished\n");

  return(0);
}

