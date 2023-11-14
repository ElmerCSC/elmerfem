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

/****************************************************************************
*                                                                           *  
*                              Elmergrid                                    *
*                                                                           *
*  This program creates very easily a structured 2D meshes with BCs.        *
*  The element types include 4, 5, 8, 9, 12 and 16-node rectangles and      *
*  3, 6 and 10-node triangles. There is also limited 3D functionality       *
*  with 8, 20 and 27-node cubes and 6-node prisms.                          *
*                                                                           *
*  The program may also be used as a mesh import and export utility. It     *
*  is able to read several different formats and writes mainly Elmer input  *
*  and output formats. The meshes may also be given some simple operations. *
*                                                                           *
*  Note: this software was initially part of my first fem implementation    *
*  the Pirfem code, then later called Quickmesh, and finally renamed to     *
*  Elmergrid. The code has never been designed and with new features the    *
*  code has eventually become very dirty and does not present my view of    *
*  good programming.                                                        *
*                                                                           *
****************************************************************************/


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <locale.h>

#define EXE_MODE 0
#define LIB_MODE 1


#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egmesh.h"
#include "egnative.h"
#include "egconvert.h"


#if EXE_MODE
#include "egoutput.h"
#else
#include "src/meshtype.h"
#endif


static struct GridType *grids;
static struct FemType data[MAXCASES];
static struct BoundaryType *boundaries[MAXCASES];
static struct ElmergridType eg;
static char Filename[MAXFILESIZE];
static int Inmethod;


int info=TRUE,nogrids=0,nomeshes=0,activemesh=0;

#if EXE_MODE
static int PartitionMesh(int nofile) 
{
  /* Partititioning related stuff */

  int noopt = 0;

  if(eg.partitions) {
    if(eg.partopt % 2 == 0) 
      PartitionSimpleElements(&data[nofile],eg.partdim,eg.periodicdim,eg.partorder,eg.partcorder,info);	
    else 
      PartitionSimpleNodes(&data[nofile],eg.partdim,eg.periodicdim,eg.partorder,eg.partcorder,info);	
    noopt = eg.partopt / 2;      
  }
#if HAVE_METIS
  if(eg.metis) {
    if(eg.partopt % 5 <= 1) 
      PartitionMetisElements(&data[nofile],eg.metis,eg.partopt % 5,info);
    else
      PartitionMetisNodes(&data[nofile],eg.metis,eg.partopt % 5,info);      
    noopt = eg.partopt / 5;      
  }
#endif
  if(eg.partitions || eg.metis ) 
    OptimizePartitioning(&data[nofile],boundaries[nofile],noopt,info);
}



static int ExportMeshDefinition(int inmethod,int outmethod,int nofile,char *filename)
{
  int i;

  switch (outmethod) {
  case 1:
    SaveElmergrid(grids,nogrids,filename,info);
    break; 
    
  case 2:
    if(data[nofile].nopartitions > 1) 
      SaveElmerInputPartitioned(&data[nofile],boundaries[nofile],filename,eg.decimals,
				eg.partitionhalo,eg.partitionindirect,info);
    else
      SaveElmerInput(&data[nofile],boundaries[nofile],filename,eg.decimals,info);
    break;
    
  case 22:  
    SaveElmerInputFemBem(&data[nofile],boundaries[nofile],filename,eg.decimals,info);
    break;
  
  case 3:
    /* Create a variable so that when saving data in ElmerPost format there is something to visualize */
    if(data[nofile].variables == 0) {
      CreateVariable(&data[nofile],1,1,0.0,"Number",FALSE);
      for(i=1;i<=data[nofile].alldofs[1];i++)
	data[nofile].dofs[1][i] = (Real)(i);	
    }
    SaveSolutionElmer(&data[nofile],boundaries[nofile],eg.saveboundaries ? MAXBOUNDARIES:0,
		      filename,eg.decimals,info=TRUE);
    break;
    
  default:
    Instructions();
    break;
  }    
  Goodbye();
}


#else


int ConvertEgTypeToMeshType(struct FemType *dat,struct BoundaryType *bound,mesh_t *mesh)
{
  int i,j,k,allocated,surfaces,elemdim;
  int sideelemtype,ind[MAXNODESD1];
 
  node_t *n;
  surface_t *b;
  element_t *e;
  edge_t *s;
  
  if(!dat->created) {
    printf("Data is not created!\n");
    return(1);
  }

  printf("Converting ElmerGrid data to ElmerGUI data\n");
  
  elemdim =  GetMaxElementDimension(dat);
  printf("Setting elements of %ddim\n",elemdim); 

  /* for mapped surfaces elemdim and space dimension may differ! */
  mesh->setDim(MAX(data->dim, elemdim));
  mesh->setNodes(dat->noknots);
  mesh->newNodeArray(mesh->getNodes());

  printf("Setting number of nodes %d\n",mesh->getNodes());   
  for(i=0; i < mesh->getNodes(); i++) {
    n = mesh->getNode(i);
    n->setX(0, dat->x[i+1]);
    n->setX(1, dat->y[i+1]);
    n->setX(2, dat->z[i+1]);
    n->setIndex(-1);
  }
  
  /* For 3D save bulk elements & boundary elements */
  if(elemdim == 3) {

    mesh->setElements(dat->noelements);
    mesh->newElementArray(mesh->getElements());
    
    printf("Setting number of 3d elements %d\n",mesh->getElements());   
    for(i = 0; i < mesh->getElements(); i++) {
      e = mesh->getElement(i);
      e->setCode(dat->elementtypes[i+1]);
      e->setNodes(e->getCode() % 100);
      e->newNodeIndexes(e->getNodes());
      e->setNature(PDE_BULK);

      for(j = 0; j < e->getNodes(); j++)
	e->setNodeIndex(j, dat->topology[i+1][j]-1);
      e->setIndex(dat->material[i+1]);
    }


    if(eg.saveboundaries) {

      allocated = FALSE;
    do_b:    surfaces = 0;
      
      
      for(j=0;j<MAXBOUNDARIES;j++) {
	if(bound[j].created == FALSE) continue;
	
	for(i=1;i<=bound[j].nosides;i++) {
	  GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],dat,ind,&sideelemtype); 
	  
	  if(sideelemtype / 100 < 3 || sideelemtype / 100 > 4) continue;
	  surfaces += 1;
	  
	  if(allocated) {
	    b = mesh->getSurface(surfaces-1);
	    b->setElements(0);
	    b->newElementIndexes(2);
	    
	    if(bound[j].parent[i]) {
	      b->setElements(b->getElements() + 1);
	      b->setElementIndex(0, bound[j].parent[i]-1);
	    }
	    else {
	      b->setElementIndex(0, -1);
	    }
	    
	    if(bound[j].parent2[i]) {
	      b->setElements(b->getElements() + 1);
	      b->setElementIndex(1, bound[j].parent2[i]-1);
	    } 
	    else {
	      b->setElementIndex(1, -1);
	    }
	    
	    b->setNormal(0, 0.0);
	    b->setNormal(1, 0.0);
	    b->setNormal(2, -1.0);
	    
	    b->setNature(PDE_BOUNDARY);
	    b->setCode(sideelemtype);
	    b->setNodes(b->getCode() % 100);
	    b->newNodeIndexes(b->getNodes());
	    for(k = 0; k < b->getNodes(); k++) 
	      b->setNodeIndex(k, ind[k]-1);
	    b->setIndex(bound[j].types[i]);
	    
	    b->setEdges(b->getNodes());
	    b->newEdgeIndexes(b->getEdges());
	    for(k=0; k < b->getEdges(); k++)
	      b->setEdgeIndex(k, -1);
	  }
	}
      }
      
      if(!allocated) {
	mesh->setSurfaces(surfaces);
	mesh->newSurfaceArray(mesh->getSurfaces());
	allocated = TRUE;
	goto do_b;
      }
    }
  }

  else if(elemdim == 2) {
    mesh->setElements(0);
    
    mesh->setSurfaces(dat->noelements);
    mesh->newSurfaceArray(mesh->getSurfaces());

    printf("Setting number of 2d elements %d\n",mesh->getSurfaces());       
    for(i = 0; i < mesh->getSurfaces(); i++) {
      b = mesh->getSurface(i);
      
      b->setElements(0);
      b->newElementIndexes(2);
      b->setElementIndex(0, -1);
      b->setElementIndex(1, -1);
      
      b->setNormal(0, 0.0);
      b->setNormal(1, 0.0);
      b->setNormal(2, -1.0);
      
      b->setCode(dat->elementtypes[i+1]);
      b->setNodes(b->getCode() % 100);
      b->newNodeIndexes(b->getNodes());
      b->setNature(PDE_BULK);

      for(j = 0; j < b->getNodes(); j++) 
	b->setNodeIndex(j, dat->topology[i+1][j]-1);
      b->setIndex(dat->material[i+1]);

      b->setEdges(b->getNodes());
      b->newEdgeIndexes(b->getEdges());
      for(k=0; k < b->getEdges(); k++)
	b->setEdgeIndex(k, -1);
    }

    allocated = FALSE;
  do_s:    surfaces = 0;
    
    for(j=0;j<MAXBOUNDARIES;j++) {
      if(bound[j].created == FALSE) continue;
      
      for(i=1;i<=bound[j].nosides;i++) {
	GetElementSide(bound[j].parent[i],bound[j].side[i],bound[j].normal[i],dat,ind,&sideelemtype); 
	
	if(sideelemtype / 100 != 2) continue;
	surfaces += 1;
	
	if(allocated) {
	  s = mesh->getEdge(surfaces-1);
	  s->setSurfaces(0);
	  s->newSurfaceIndexes(2);
	  
	  if(bound[j].parent[i]) {
	    s->setSurfaces(s->getSurfaces() + 1);
	    s->setSurfaceIndex(0, bound[j].parent[i]-1);
	  }
	  else {
	    s->setSurfaceIndex(0, -1);
	  }
	  if(bound[j].parent2[i]) {
	    s->setSurfaces(s->getSurfaces() + 1);
	    s->setSurfaceIndex(1, bound[j].parent2[i]-1);
	  }
	  else {
	    s->setSurfaceIndex(1, -1);
	  }
	  s->setCode(sideelemtype);
	  s->setNodes(s->getCode() % 100);
	  s->newNodeIndexes(s->getNodes());
	  s->setNature(PDE_BOUNDARY);

	  for(k = 0; k < s->getNodes(); k++) 
	    s->setNodeIndex(k, ind[k]-1);
	  s->setIndex(bound[j].types[i]);
	}
      }
   }

    if(!allocated) {
      mesh->setEdges(surfaces);
      mesh->newEdgeArray(mesh->getEdges());
      allocated = TRUE;
      goto do_s;
    }
  }
  else {
    printf("Implemented only for element dimensions 2 and 3 (not %d)\n",elemdim);
  }

  printf("Done converting mesh\n");
  
  return(0);
}
#endif


static int DetermineFileType(const char *filename,int info)
{
  int mode;
  mode = -1;

  if(!strstr(filename,".")) {
    if(info) printf("There cannot be a filetype suffix without a dot: %s\n",filename);
    return(mode);
  }
  else {
    if(strstr(filename,".eg")) mode = 0;
    else if(strstr(filename,".grd")) mode = 1;
    else if(strstr(filename,".ep")) mode = 3;
    else if(strstr(filename,".inp")) mode = 5;
    else if(strstr(filename,".nas")) mode = 6;
    else if(strstr(filename,".fdneut") || strstr(filename,".FDNEUT")) mode = 7;
    else if(strstr(filename,".unv")) mode = 8;
    else if(strstr(filename,".mphtxt")) mode = 9;
    else if(strstr(filename,".msh")) mode = 14;
    else if(strstr(filename,".plt")) mode = 16;
  }
  if(mode == -1) {
    if(info) printf("Could not determine the filetype based on the suffix\n");
  }

  if(info) printf("Filetype determined by suffix: %d\n",mode);

  return(mode);
}



static int ImportMeshDefinition(int inmethod,int nofile,char *filename,int *nogrids)
{
  int i,k,errorstat = 0,dim;
  static int visited = FALSE;


  *nogrids = 0;
  if(!visited) {
    printf("Initializing structures for max. %d meshes\n",MAXCASES);
    for(k=0;k<MAXCASES;k++) {
      boundaries[k] = (struct BoundaryType*)
	malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
      for(i=0;i<MAXBOUNDARIES;i++) {
	boundaries[k][i].created = FALSE; 
	boundaries[k][i].nosides = 0;
      }    
    }
    for(k=0;k<MAXCASES;k++) data[k].created=FALSE;
       
    visited = TRUE;
  }

  /* Native format of ElmerGrid gets specieal treatment */
  switch (inmethod) {
    
  case 1: 
    printf("Loading ElmerGrid format file\n");
    info = TRUE;
    errorstat = LoadElmergrid(&grids,nogrids,eg.filesin[nofile],info);
    if(errorstat == 1) {
      dim = eg.dim;
      CreateExampleGrid(dim,&grids,nogrids,info);
      SaveElmergrid(grids,*nogrids,eg.filesin[nofile],info); 
      printf("Because file %s didn't exist, it was created for you.\n",eg.filesin[nofile]);
      return(errorstat);
    }
    if(*nogrids) LoadCommands(eg.filesin[nofile],&eg,grids,2,info); 
    break;
    
#if EXE_MODE
  case 2: 
    errorstat = LoadElmerInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 3: 
    errorstat = LoadSolutionElmer(&(data[nofile]),TRUE,eg.filesin[nofile],info);
    break;
#endif

  case 4:
    errorstat = LoadAnsysInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 5: 
    errorstat = LoadAbaqusInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 6:
    errorstat = LoadNastranInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 7:
    errorstat = LoadFidapInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    eg.bulkorder = TRUE;
    eg.boundorder = TRUE;
    break;

  case 8:
    errorstat = LoadUniversalMesh(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 9:
    errorstat = LoadComsolMesh(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 10:
    errorstat = LoadFieldviewInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 11:
    errorstat = LoadTriangleInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 12:
    errorstat = LoadMeditInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 13:
    errorstat = LoadGidInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
    break;

  case 14:
    data[nofile].dim = 3; /* default dim 3 with gmsh*/
    eg.dim = 3;
    errorstat = LoadGmshInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],
			      eg.multidim,info);
    break;

   
  case 16:
    errorstat = LoadCGsimMesh(&(data[nofile]),eg.filesin[nofile],info);
    break;

#if EXE_MODE
  case 15: 
    if(info) printf("Partitioned solution is fused on-the-fly therefore no other operations may be performed.\n");
    FuseSolutionElmerPartitioned(eg.filesin[nofile],eg.filesout[nofile],eg.decimals,
				 eg.saveinterval[0],eg.saveinterval[1],eg.saveinterval[2],info);
    if(info) printf("Finishing with the fusion of partitioned Elmer solutions\n");
    Goodbye();
    break;

  default:
    errorstat = 1;
    Instructions();
#endif

  }  

#if LIB_MODE
  /* ElmerGrid cannot deal with three different dimensions or orphan nodes */
  eg.removelowdim = TRUE;
  eg.removeunused = TRUE;
#endif

  if(!errorstat) *nogrids = MAX(*nogrids,1);

  return(errorstat);
}





static int ManipulateMeshDefinition(int inmethod,int outmethod,Real relh)
{
  static int visited = FALSE;
  int i,j,k;
  Real mergeeps;

  printf("Manipulate mesh definition with formats: %d %d\n",inmethod,outmethod);
  
  if(inmethod == 1 && outmethod != 1) {
    if(visited) {
      printf("Deallocating structures from previous call\n");
      for(k=0;k<MAXCASES;k++) {	    
	if(data[k].created) {
	  DestroyKnots(&data[k]);
	  for(i=0;i<MAXBOUNDARIES;i++) {
	    DestroyBoundary(&boundaries[k][i]);
	  }
	}
      }
    }
    printf("Starting to create ElmerGrid meshes\n");
    for(k=0;k<nogrids;k++) 
      CreateElmerGridMesh(&(grids[k]),&(data[k]),boundaries[k],relh,info);
    nomeshes = nogrids;
    printf("Created %d ElmerGrid meshes\n",nomeshes);

    /* Cancel the automatic effect of -autoclean flag for this input format */
    eg.bulkorder = FALSE;
    eg.boundorder = FALSE;
  }

  visited = TRUE;

  /* At first instance perform operations that should rather be done before extrusion 
     or mesh union. */
  for(k=0;k<nomeshes;k++) {
    
    /* Make the discontinuous boundary needed, for example, in poor thermal conduction */
    if(!eg.discont) {
      for(j=0;j<grids[k].noboundaries;j++) 
	if(grids[k].boundsolid[j] == 2) {
	  eg.discontbounds[eg.discont] = grids[k].boundtype[j];
	  eg.discont++;	  
	}
    }
    if(eg.discont) {
      for(i=1;i<=eg.discont;i++) 
	SetDiscontinuousBoundary(&(data[k]),boundaries[k],eg.discontbounds[i-1],2,info);
    }
    
    /* Make a connected boundary (specific to Elmer format) needed in linear constraints */
    /* No longer needed in serial at least!! */
    /* for(i=1;i<=eg.connect;i++) 
       SetConnectedBoundary(&(data[k]),boundaries[k],eg.connectbounds[i-1],i,info); */
  
    /* Divide quadrilateral meshes into triangular meshes */
    if(eg.triangles || grids[k].triangles == TRUE) {
      Real criticalangle;
      criticalangle = MAX(eg.triangleangle, grids[k].triangleangle);
      ElementsToTriangles(&data[k],boundaries[k],criticalangle,info);
    }

    /* Make a boundary layer with two different methods */
    if(eg.layers > 0) 
      CreateBoundaryLayer(&data[k],boundaries[k],eg.layers,
			  eg.layerbounds, eg.layernumber, eg.layerratios, eg.layerthickness,
			  eg.layerparents, eg.layermove, eg.layereps, info);
    else if(eg.layers < 0) 
      CreateBoundaryLayerDivide(&data[k],boundaries[k],abs(eg.layers),
				eg.layerbounds, eg.layernumber, eg.layerratios, eg.layerthickness,
				eg.layerparents, info);
  }

  /* Take up the infor on rotation */
  for(k=0;k<nogrids;k++) 
    if( grids[k].rotatecurve ) {
      eg.rotatecurve = TRUE;
      eg.curvezet = grids[k].curvezet;
      eg.curverad = grids[k].curverad;
      eg.curveangle = grids[k].curveangle;
    }
  
  if(outmethod != 1 && eg.dim != 2) { 
    j = MAX(1,nogrids);
    for(k=0;k<j;k++) {
      if(grids[k].dimension == 3 || grids[k].rotate) {
	CreateKnotsExtruded(&(data[k]),boundaries[k],&(grids[k]),
			    &(data[j]),boundaries[j],info);
#if LIB_MODE
	activemesh = j;
	nomeshes = j+1;
#endif
#if EXE_MODE
	data[k] = data[j];
	boundaries[k] = boundaries[j];
#endif
      }
    }
  }

  
  /* Unite meshes if there are several of them */
  if(eg.unitemeshes) {
    for(k=1;k<nomeshes;k++)
      UniteMeshes(&data[0],&data[k],boundaries[0],boundaries[k],eg.unitenooverlap,info);
    nomeshes = 1;
  }

  
  for(k=0;k<nomeshes;k++) {

    /* If the original mesh was given in polar coordinates make the transformation into cartesian ones */
    if(eg.polar || data[k].coordsystem == COORD_POLAR) {
      if(!eg.polar) eg.polarradius = grids[k].polarradius;
      PolarCoordinates(&data[k],eg.polarradius,info);
    }
    /* If the original mesh was given in cylindrical coordinates make the transformation into cartesian ones */
    if(eg.cylinder || data[k].coordsystem == COORD_CYL) {
      CylinderCoordinates(&data[k],info);
    }
    if(eg.clone[0] || eg.clone[1] || eg.clone[2]) {
      CloneMeshes(&data[k],boundaries[k],eg.clone,eg.clonesize,FALSE,info);
      mergeeps = fabs(eg.clonesize[0]+eg.clonesize[1]+eg.clonesize[2]) * 1.0e-8;
      MergeElements(&data[k],boundaries[k],eg.order,eg.corder,mergeeps,TRUE,TRUE);
    }
    
    /* Reduce element order if requested */
    if(nogrids && grids[k].reduceordermatmax) {
      eg.reduce = TRUE;
      eg.reducemat1 = grids[k].reduceordermatmin;
      eg.reducemat2 = grids[k].reduceordermatmax;
    }
    if(eg.reduce) 
      ReduceElementOrder(&data[k],eg.reducemat1,eg.reducemat2);

    /* Increase element order */
    if(eg.increase) 
      IncreaseElementOrder(&data[k],TRUE);

    if(eg.merge) 
      MergeElements(&data[k],boundaries[k],eg.order,eg.corder,eg.cmerge,FALSE,TRUE);
#if HAVE_METIS
    else if(eg.order == 3) 
      ReorderElementsMetis(&data[k],TRUE);
#endif
    else if(eg.order) 
      ReorderElements(&data[k],boundaries[k],eg.order,eg.corder,TRUE);
    
    if(eg.bulkbounds || eg.boundbounds) 
      SideAndBulkBoundaries(&data[k],boundaries[k],&eg,info);

    RotateTranslateScale(&data[k],&eg,info);
    
    if(eg.rotatecurve) 
      CylindricalCoordinateCurve(&data[k],eg.curvezet,eg.curverad,eg.curveangle);
    
    if(eg.removelowdim) 
      RemoveLowerDimensionalBoundaries(&data[k],boundaries[k],info);

    if(eg.removeunused) 
      RemoveUnusedNodes(&data[k],info);

    if(eg.sidemappings || eg.bulkmappings)  
      SideAndBulkMappings(&data[k],boundaries[k],&eg,info);

    if(eg.boundorder || eg.bcoffset) 
      RenumberBoundaryTypes(&data[k],boundaries[k],eg.boundorder,eg.bcoffset,info);

    if(eg.bulkorder) 
      RenumberMaterialTypes(&data[k],boundaries[k],info);

    if(eg.periodicdim[0] || eg.periodicdim[1] || eg.periodicdim[2]) 
      FindPeriodicNodes(&data[k],eg.periodicdim,info);
  }
  return 0;
}


  

#if LIB_MODE
int eg_loadmesh(const char *filename)
{
  static int inmethod,errorstat,info;

  setlocale(LC_ALL, "C");  
    
  strcpy(Filename,filename);
  info = TRUE;
  if(info) printf("\nElmerGrid checking filename suffix for file: %s\n",filename);

  inmethod = DetermineFileType(filename,info);
  Inmethod = inmethod;

  if(inmethod < 0) 
    errorstat = 1;
  else
    errorstat = 0;
  
  if(info) printf("Initialized the filetype\n");
  return(errorstat);
}



int eg_transfermesh(mesh_t *mesh,const char *str)
{
  int i,k,inmethod,outmethod,errorstat,nofile;
  info = TRUE;
  static char arguments[10][15],**argv;
  int argc;
  static int visited = FALSE;
  char filename[MAXFILESIZE];

  setlocale(LC_ALL, "C");  
  
  activemesh = 0;
  nofile = 0;
  nomeshes = 0;
  nogrids = 0;
  info = TRUE;

  if(!visited) {
    grids = (struct GridType*)malloc((size_t) (MAXCASES)*sizeof(struct GridType));     
    argv = (char**) malloc((size_t) 10*sizeof(char*));
  }
  visited++;

  InitParameters(&eg);
  InitGrid(grids);

  strcpy(filename,Filename);
  if(info) printf("\nElmerGrid loading data from file: %s\n",filename);

  inmethod = Inmethod;
  if(inmethod < 0) return(1);

  if(inmethod == 0) {
    errorstat = LoadCommands(filename,&eg,grids,1,info);
    inmethod = eg.inmethod;
    info = !eg.silent;  
  } 
  else {
    eg.inmethod = inmethod;
  }
  strcpy(eg.filesin[0],filename);

  printf("Import mesh definition\n");
  errorstat = ImportMeshDefinition(inmethod,nofile,filename,&nogrids);

  if(errorstat) return(errorstat);
  nomeshes += nogrids;
  
  if(info) printf("\nElmerGrid manipulating and importing data\n");

  mesh->setNodes(0);
  mesh->setPoints(0);
  mesh->setEdges(0);
  mesh->setSurfaces(0);
  mesh->setElements(0);

  if(nomeshes == 0) {
    printf("No mesh to work with!\n");
    return(1);
  }
  
  /* Checking in-line parameters */
  argc = StringToStrings(str,arguments,10,' ');
  for(i=0;i<argc;i++) argv[i] = &arguments[i][0];

  printf("Number of inline arguments: %d\n",argc);
  
  errorstat = InlineParameters(&eg,argc,argv,0,info);
  if(errorstat) printf("Errorstat for inline parameters: %d\n",errorstat);
  
  inmethod = eg.inmethod;
  outmethod = 0;

  printf("Input method: %d\n",inmethod);
  
  ManipulateMeshDefinition(inmethod,outmethod,eg.relh);
  
  errorstat = ConvertEgTypeToMeshType(&data[activemesh],boundaries[activemesh],mesh);
  return(errorstat);

  for(k=0;k<MAXCASES;k++) {
    DestroyKnots(&data[k]);
    for(i=0;i<MAXBOUNDARIES;i++) 
      DestroyBoundary(&boundaries[k][i]);
  }
  if(info) printf("Done destroying structures\n");
}


#if 0
int main(int argc, char *argv[])
{
  char prefix[MAXFILESIZE];
  class mesh_t mesh;

  eg_loadmesh("apu.grd");
  eg_transfermesh("test",&mesh); 
}
#endif
#endif


#if EXE_MODE
int main(int argc, char *argv[])
{
  char prefix[MAXFILESIZE],filename[MAXFILESIZE];
  Real relh;  
  static int i,j,k,l,inmethod,outmethod,errorstat;
  static int nofile,dim;
  static Real mergeeps;
  long ii;

  if(info) printf("\nStarting program Elmergrid, compiled on %s\n", __DATE__ );
  
  InitParameters(&eg);
  grids = (struct GridType*)malloc((size_t) (MAXCASES)*sizeof(struct GridType));     
  InitGrid(grids);
  info = TRUE;

  if(argc <= 2) {
    errorstat = LoadCommands(argv[1],&eg,grids,argc-1,info);     
    if(errorstat) {
      if(argc <= 1) Instructions();
      Goodbye();
    }
  }
  else if(argc == 3) {
    Instructions();
    Goodbye();
  } 
  else {
    errorstat = InlineParameters(&eg,argc,argv,4,info);
    if(errorstat) Goodbye();
  }
  inmethod = eg.inmethod;
  outmethod = eg.outmethod;  
  info = !eg.silent;  
  dim = eg.dim;
  relh = eg.relh;  
  if(!outmethod || !inmethod) {
    printf("Please define the input and output formats\n");
    Goodbye();
  }

  /**********************************/
  if(info) printf("\nElmergrid loading data:\n");
   
  nofile = 0;
  nomeshes = 0;
  nogrids = 0;

  for(nofile=0;nofile<eg.nofilesin;nofile++) {
    errorstat = ImportMeshDefinition(inmethod,nofile,eg.filesin[nofile],&nogrids);
    if(errorstat) Goodbye();
    nomeshes += nogrids;
  }

  /***********************************/
  if(info) printf("\nElmergrid creating and manipulating meshes:\n");
  ManipulateMeshDefinition(inmethod,outmethod,relh);

  /* Partititioning related stuff */
  for(k=0;k<nomeshes;k++) 
    PartitionMesh(nofile); 

  /********************************/
  if(info) printf("\nElmergrid saving data:\n");
  sprintf(prefix,"%s",eg.filesout[0]);
  for(nofile=0;nofile<nomeshes;nofile++) {
    if(nomeshes == 1) 
      sprintf(filename,"%s",prefix);
    else 
      sprintf(filename,"%s%d",prefix,nofile+1);     
    ExportMeshDefinition(inmethod,outmethod,nofile,filename);
  }
  Goodbye();
}
#endif



