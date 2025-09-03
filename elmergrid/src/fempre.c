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

#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egmesh.h"
#include "egextra.h"
#include "egnative.h"
#include "egparallel.h"
#include "egconvert.h"
#include "egexport.h"


int main(int argc, char *argv[])
{
  int i,j,k,l,inmethod,outmethod,info,errorstat;
  int nogrids,nogrids0,nomeshes,nofile,dim,elementsredone=0;
  int nodes3d,elements3d,showmem;
  Real mergeeps;
  char prefix[MAXFILESIZE];
  struct GridType *grids;
  struct CellType *cell[MAXCASES];
  struct FemType data[MAXCASES];
  struct BoundaryType *boundaries[MAXCASES];
  struct ElmergridType eg;

  showmem = TRUE;

  printf("\nStarting program Elmergrid, compiled on %s\n", __DATE__ );

  InitParameters(&eg);

  grids = (struct GridType*)malloc((size_t) (MAXCASES)*sizeof(struct GridType));     
  InitGrid(grids);
  info = TRUE;

  if(argc <= 1) {
    errorstat = LoadCommands(argv[1],&eg,grids,argc-1,info);     
    if(errorstat) {
      Instructions();
      Goodbye();
    }
  }
  else if(argc == 2) {
    errorstat = LoadCommands(argv[1],&eg,grids,argc-1,info);     
    if(errorstat) Goodbye();
  }
  else if(argc < 4) {
    Instructions();
    Goodbye();
  } 
  else {
    errorstat = InlineParameters(&eg,argc,argv,4,info);
    if(errorstat) Goodbye();
  }


  if(!eg.outmethod || !eg.inmethod) {
    printf("Please define the input and output formats\n");
  }
  if(eg.inmethod != 1) {
    if(eg.outmethod == 1 || eg.outmethod == 8 || eg.outmethod == 9 || eg.outmethod == 10) {
      printf("input of type %d can't create output of type %d\n",
	     eg.inmethod,eg.outmethod);
      errorstat++;
      Goodbye();
    }
  }

  if(eg.timeron) timer_activate(eg.infofile);

  /**********************************/
  printf("\nElmergrid loading data:\n");
  printf(  "-----------------------\n");

  dim = eg.dim;
  nofile = 0;
  nomeshes = 0;
  nogrids = 0;
  inmethod = eg.inmethod;
  outmethod = eg.outmethod;


 read_another_file:    

  timer_show();
  
  switch (inmethod) {

  case 1:        
    nogrids0 = nogrids;
    if(LoadElmergrid(&grids,&nogrids,eg.filesin[nofile],info) == 1) {   
      CreateExampleGrid(eg.dim,&grids,&nogrids,info);
      SaveElmergrid(grids,nogrids,eg.filesin[nofile],info); 
      printf("Because file %s didn't exist, it was created for you.\n",eg.filesin[nofile]);
      Goodbye();
    }
    LoadCommands(eg.filesin[nofile],&eg,grids,2,info); 
    for(k=nogrids0;k < nogrids;k++) {
      SetElementDivision(&(grids[k]),eg.relh,info);
    } 

    break;

  case 2: 
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if(LoadElmerInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],
		      !eg.usenames,info))
      Goodbye();


    nomeshes++;
    break;

  case 3: 
    if(LoadSolutionElmer(&(data[nofile]),TRUE,eg.filesin[nofile],info)) 
      Goodbye();
    nomeshes++;
    break;

  case 4:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    if(LoadAnsysInput(&(data[0]),boundaries[0],eg.filesin[nofile],info)) 
      Goodbye();
    nomeshes++;
    break;

  case 5: 
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if(LoadAbaqusInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE)) 
      Goodbye();
    nomeshes++;
    break;

  case 6:
    if(LoadAbaqusOutput(&(data[nofile]),eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 7:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if(LoadFidapInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    
    eg.bulkorder = TRUE;
    eg.boundorder = TRUE;

    if(!eg.usenames) data[nofile].boundarynamesexist = data[nofile].bodynamesexist = FALSE;
  
    nomeshes++;
    break;

  case 8:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadUniversalMesh(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

 case 9:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }   
    if(LoadComsolMesh(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info)) 
      Goodbye();
    nomeshes++;
    break;

  case 10:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	    
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if(LoadFieldviewInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 11:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadTriangleInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 12:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadMeditInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 13:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadGidInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 14:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    data[nofile].dim = (eg.dim >= 1 && eg.dim <= 3) ? eg.dim : 3; /* default dim 3 with gmsh*/
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }

    if (LoadGmshInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],
		      eg.multidim,TRUE))
      Goodbye();
    nomeshes++;    
    break;

  case 15: 
    if(info) printf("Partitioned solution is fused on-the-fly therefore no other operations may be performed.\n");
    FuseSolutionElmerPartitioned(eg.filesin[nofile],eg.filesout[nofile],eg.decimals,eg.partjoin,
				 eg.saveinterval[0],eg.saveinterval[1],eg.saveinterval[2],info);
    if(info) printf("Finishing with the fusion of partitioned Elmer solutions\n");
    Goodbye();
    break;

  case 16:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadFvcomMesh(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;
    
  case 17:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadNastranInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 18:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
   
    if(LoadCGsimMesh(&(data[nofile]),eg.filesin[nofile],info))
       Goodbye();
    nomeshes++;
    break;

  case 19:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadGeoInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 20:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadFluxMesh(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  case 21:
    boundaries[nofile] = (struct BoundaryType*)
      malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
    for(i=0;i<MAXBOUNDARIES;i++) {
      boundaries[nofile][i].created = FALSE; 
      boundaries[nofile][i].nosides = 0;
    }
    if (LoadFluxMesh3D(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],TRUE))
      Goodbye();
    nomeshes++;
    break;

  default:
    Instructions();
    Goodbye();
  }  

  nofile++;
  if(nofile < eg.nofilesin) {
    printf("\nElmergrid loading data from another file:\n");
    goto read_another_file;
  }

  /***********************************/


 redoelements:

  printf("\nElmergrid creating and manipulating meshes:\n");
  printf(  "-------------------------------------------\n");
  timer_show();


  if(nogrids > nomeshes && outmethod != 1) { 

    nomeshes = nogrids;
    for(k=0;k<nogrids;k++) {

      CreateCells(&(grids[k]),&(cell[k]),info);  
      CreateKnots(&(grids[k]),cell[k],&(data[k]),0,0);

      boundaries[k] = (struct BoundaryType*)
	malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 

      for(j=0;j<MAXBOUNDARIES;j++) {
	boundaries[k][j].created = FALSE;
	boundaries[k][j].nosides = FALSE;
      }

      if(grids[k].noboundaries > 0) {
	for(j=0;j<grids[k].noboundaries;j++) {
	  if(grids[k].boundsolid[j] < 4) {
	    CreateBoundary(cell[k],&(data[k]),&(boundaries[k][j]),
			   grids[k].boundext[j],grids[k].boundint[j],
			   1,grids[k].boundtype[j],info);  
	  } 
	  else { 
	    CreatePoints(cell[k],&(data[k]),&(boundaries[k][j]),
			 grids[k].boundext[j],grids[k].boundint[j],
			 grids[k].boundsolid[j],grids[k].boundtype[j],info); 	    
	  }
	}
      }
    }
  }

  /* In some formats the dimension for curved 2D meshes seems to be set to 2.
     This should fix the problem for all input types. */
  if( data->dim < 3 ) {
    data->dim = GetCoordinateDimension(data,info);
  }


  /* Make the discontinuous boundary needed, for example, in poor thermal conduction */
  for(k=0;k<nomeshes;k++) {
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
  }


  /* Divide quadrilateral meshes into triangular meshes */
  for(k=0;k<nomeshes;k++) 
    if(nogrids && (eg.triangles || grids[k].triangles == TRUE)) {
      Real criticalangle;
      criticalangle = MAX(eg.triangleangle , grids[k].triangleangle);
      ElementsToTriangles(&data[k],boundaries[k],criticalangle,info);
    }


  /* Make a boundary layer with two different methods */
  if(eg.layers > 0) 
    for(k=0;k<nomeshes;k++) 
      CreateBoundaryLayer(&data[k],boundaries[k],eg.layers,
			  eg.layerbounds, eg.layernumber, eg.layerratios, eg.layerthickness,
			  eg.layerparents, eg.layermove, eg.layereps, info);

  else if(eg.layers < 0) 
    for(k=0;k<nomeshes;k++) 
      CreateBoundaryLayerDivide(&data[k],boundaries[k],abs(eg.layers),
				eg.layerbounds, eg.layernumber, eg.layerratios, eg.layerthickness,
				eg.layerparents, info);

  /* Take up the infor on rotation */
  for(k=0;k<nogrids;k++) 
    if( grids[k].rotatecurve ) {
      eg.rotatecurve = TRUE;
      eg.curvezet = grids[k].curvezet;
      eg.curverad = grids[k].curverad;
      eg.curveangle = grids[k].curveangle;
    }


  if(outmethod != 1 && dim != 2 && eg.dim != 2) { 
    j = MAX(nogrids,1);


    for(k=0;k<j;k++) {
      if(grids[k].dimension == 3 || grids[k].rotate) {

	boundaries[j] = (struct BoundaryType*)
	  malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
	
	for(i=0;i<MAXBOUNDARIES;i++) 
	  boundaries[j][i].created = FALSE;

	CreateKnotsExtruded(&(data[k]),boundaries[k],&(grids[k]),
			    &(data[j]),boundaries[j],info);

	if(nogrids) {
	  elements3d = MAX(eg.elements3d, grids[k].wantedelems3d);
	  nodes3d = MAX(eg.nodes3d, grids[k].wantednodes3d);

	  if(elements3d) {
	    if( abs(data[j].noelements - elements3d) / (1.0*elements3d) > 0.01 && elementsredone < 5 ) {
	      grids[k].wantedelems *= pow(1.0*elements3d / data[j].noelements, (2.0/3.0));
	      elementsredone++;
	    }
	    else elementsredone = 0;
	  }
	  else if(nodes3d) {
	    if( abs(data[j].noknots - nodes3d) / (1.0*nodes3d) > 0.01 && elementsredone < 5 ) {
	      grids[k].wantedelems *= pow(1.0*nodes3d / data[j].noknots, (2.0/3.0));
	      elementsredone++;
	    }
	    else elementsredone = 0;
	  }

	  if(elementsredone) {
	    nomeshes = 0;
	    for(i=0;i < nogrids;i++) SetElementDivision(&(grids[i]),eg.relh,info);
	    
	    DestroyKnots(&data[j]);
	    DestroyKnots(&data[k]);
	    free(cell[k]);
	    
	    if(info) printf("Iteration %d of elements number targeting %d in 2D\n",
			    elementsredone,grids[k].wantedelems);
	    goto redoelements;
	  }
	}	

	data[k] = data[j];
	boundaries[k] = boundaries[j];
      }
    }
  }

  /* If the original mesh was given in polar coordinates make the transformation into cartesian ones */
  for(k=0;k<nomeshes;k++) {
    if(eg.polar || data[k].coordsystem == COORD_POLAR) {
      if(!eg.polar) eg.polarradius = grids[k].polarradius;
      PolarCoordinates(&data[k],eg.polarradius,info);
    }
  }

  /* If the original mesh was given in cylindrical coordinates make the transformation into cartesian ones */
  for(k=0;k<nomeshes;k++) {
    if(eg.cylinder || data[k].coordsystem == COORD_CYL) {
      CylinderCoordinates(&data[k],info);
    }
  }

  if(1) for(k=0;k<nomeshes;k++) 
    RotateTranslateScale(&data[k],&eg,info);


  /* Rotate may apply to 2d and 3d geometries as well */
  for(k=0;k<nomeshes;k++) 
    if(eg.rotatecurve) 
      CylindricalCoordinateCurve(&data[k],eg.curvezet,eg.curverad,eg.curveangle);

  /* Unite meshes if there are several of them */
  if(eg.unitemeshes) {
    for(k=1;k<nomeshes;k++)
      UniteMeshes(&data[0],&data[k],boundaries[0],boundaries[k],eg.unitenooverlap,info);
    nomeshes = nogrids = 1;
  }
  
  if(eg.clone[0] || eg.clone[1] || eg.clone[2]) {
    for(k=0;k<nomeshes;k++) {
      CloneMeshes(&data[k],boundaries[k],eg.clone,eg.clonesize,eg.cloneinds,info);
    }
  }

  if(eg.mirror[0] || eg.mirror[1] || eg.mirror[2]) {
    for(k=0;k<nomeshes;k++) {
      MirrorMeshes(&data[k],boundaries[k],eg.mirror,FALSE,eg.clonesize,eg.mirrorbc,info);
      mergeeps = fabs(eg.clonesize[0]+eg.clonesize[1]+eg.clonesize[2]) * 1.0e-8;
      MergeElements(&data[k],boundaries[k],eg.order,eg.corder,mergeeps,FALSE,TRUE);
    }
  }

  /* Naming convection for the case of several meshes */
  if(nomeshes > 1) {
    strcpy(prefix,eg.filesout[0]);
    for(k=0;k<nomeshes;k++)
      sprintf(eg.filesout[k],"%s%d",prefix,k+1);
  }

  for(k=0;k<nomeshes;k++) {
    if(nogrids && grids[k].reduceordermatmax) {
      eg.reduce = TRUE;
      eg.reducemat1 = grids[k].reduceordermatmin;
      eg.reducemat2 = grids[k].reduceordermatmax;
    }
    if(eg.reduce) 
      ReduceElementOrder(&data[k],eg.reducemat1,eg.reducemat2);
  }

  for(k=0;k<nomeshes;k++) 
    if(eg.increase) IncreaseElementOrder(&data[k],TRUE);
 
  for(k=0;k<nomeshes;k++) {
    if(eg.merge) 
      MergeElements(&data[k],boundaries[k],eg.order,eg.corder,eg.cmerge,FALSE,TRUE);
    else if(eg.order == 3) 
#if USE_METIS
      ReorderElementsMetis(&data[k],TRUE);
#else
      printf("Cannot order nodes by Metis as it is not even compiled!\n");
#endif    
    else if(eg.order) 
      ReorderElements(&data[k],boundaries[k],eg.order,eg.corder,TRUE);
    
    if(eg.isoparam) 
      IsoparametricElements(&data[k],boundaries[k],TRUE,info);
  }  

  for(k=0;k<nomeshes;k++) {
    if(eg.bulkbounds || eg.boundbounds) 
      SideAndBulkBoundaries(&data[k],boundaries[k],&eg,info);
#if 0
    if(eg.bulkbounds || eg.boundbounds)
      SeparateCartesianBoundaries(&data[k],boundaries[k],info);
#endif
  }

  if(0) for(k=0;k<nomeshes;k++) 
    RotateTranslateScale(&data[k],&eg,info);

  if(eg.removelowdim) 
    for(k=0;k<nomeshes;k++)
      RemoveLowerDimensionalBoundaries(&data[k],boundaries[k],info);

  if(eg.removeintbcs) 
    for(k=0;k<nomeshes;k++)
      RemoveInternalBoundaries(&data[k],boundaries[k],info);

  if(eg.removeunused) 
    for(k=0;k<nomeshes;k++)
      RemoveUnusedNodes(&data[k],info);

  if(eg.boundorder || eg.bcoffset) 
    for(k=0;k<nomeshes;k++) 
      RenumberBoundaryTypes(&data[k],boundaries[k],eg.boundorder,eg.bcoffset,info);

  if(eg.bulkorder) 
    for(k=0;k<nomeshes;k++) 
      RenumberMaterialTypes(&data[k],boundaries[k],info);

  if(eg.sidemappings || eg.bulkmappings)  
    for(k=0;k<nomeshes;k++) 
      SideAndBulkMappings(&data[k],boundaries[k],&eg,info);
  
  if(eg.coordinatemap[0] && eg.coordinatemap[1] && eg.coordinatemap[2] ) {
    Real *tmpcoord[3];

    if(info) printf("Mapping coordinates with [%d %d %d]\n",
		    eg.coordinatemap[0],eg.coordinatemap[1],eg.coordinatemap[2]);
    for(k=0;k<nomeshes;k++) {
      tmpcoord[0] = data[k].x;
      tmpcoord[1] = data[k].y;
      tmpcoord[2] = data[k].z;    
      
      data[k].x = tmpcoord[eg.coordinatemap[0]-1];
      data[k].y = tmpcoord[eg.coordinatemap[1]-1];
      data[k].z = tmpcoord[eg.coordinatemap[2]-1];      

      if(eg.coordinatemap[2] != 3) data[k].dim = 3;
    }
  }


  for(k=0;k<nomeshes;k++) {
    int partoptim, partbcoptim, partopt, fail, partdual;

    if( eg.metis == 1 ) {
      if(info) printf("One Metis partition requested, enforcing serial mode\n");
      eg.metis = 0;
    }

    if( eg.partitions == 1 ) {
      if(!eg.connect) {
	if(info) printf("One geometric partition requested, enforcing serial mode\n");
	eg.partitions = 0;
      }
    }


    partoptim = eg.partoptim;
    partbcoptim = eg.partbcoptim;
    partdual = eg.partdual;

    if(eg.partitions || eg.metis) {
      printf("\nElmergrid partitioning meshes:\n");
      printf(  "------------------------------\n");
      timer_show();

      partopt = eg.partopt;

      /* Make a connected boundary needed to enforce nodes to same partitioning */
      for(i=1;i<=eg.connect;i++) 
	SetConnectedNodes(&(data[k]),boundaries[k],eg.connectbounds[i-1],eg.connectboundsset[i-1],info);      

      if(eg.periodicdim[0] || eg.periodicdim[1] || eg.periodicdim[2]) 
	FindPeriodicNodes(&data[k],eg.periodicdim,info);

      if(eg.partitions) {
	if( partopt == -1 ) partopt = partdual;
	if(partopt == 0) 
	  PartitionSimpleElements(&data[k],&eg,boundaries[k],eg.partdim,eg.periodicdim,
				  eg.partorder,eg.partcorder,eg.parttol,info);	
	else if(partopt == 2) 
	  PartitionSimpleElementsNonRecursive(&data[k],eg.partdim,eg.periodicdim,info);	
	else if(partopt == 3) 
	  PartitionSimpleElementsRotational(&data[k],eg.partdim,eg.periodicdim,info);	
	else
	  PartitionSimpleNodes(&data[k],eg.partdim,eg.periodicdim,eg.partorder,
			       eg.partcorder,eg.parttol,info);	
      }
#if USE_METIS
      if(eg.metis) {
	if( partopt < 0 || partopt > 4 ) {
	  printf("Metis optional parameter should be in range [0,4], not %d\n",partopt);
	  bigerror("Cannot perform partitioning");
	}
	if( partopt <= 1 && eg.connect ) {
	  printf("Elemental Metis partition cannot deal with constraints!\n");
	  printf("Using Metis algorithms based on the dual graph\n");
	  partopt = 2;	  
	}
       	
	if(partopt <= 1) {
	  if(!partdual) partdual = partopt;
	  fail = PartitionMetisMesh(&data[k],&eg,eg.metis,partdual,info);
	}
	else {
	  PartitionMetisGraph(&data[k],boundaries[k],&eg,eg.metis,partopt,partdual,info);
	} 
      }
#endif
      if(data[k].periodicexist) 
	FindPeriodicParents(&data[k],boundaries[k],info);	

      timer_show();
      /* This is the only routine that affects the ownership of elements */
      if( partbcoptim ) {
	OptimizePartitioningAtBoundary(&data[k],boundaries[k],info);
      }

      OptimizePartitioning(&data[k],boundaries[k],!partoptim,eg.partbw,info);

      if(data[k].periodicexist) {
	free_Ivector(data[k].periodic,1,data[k].noknots);
	data[k].periodicexist = FALSE;
      }
      if(info) printf("Partitioning routines finished!\n");	
    }
    timer_show();
  }

  if(eg.timeron) {
    if(info) printf("Saving size info for timing purposes\n");
    for(k=0;k<nomeshes;k++) 
      SaveSizeInfo(&data[k],boundaries[k],eg.infofile,info);
  }

  if(info) for(k=0;k<nomeshes;k++)
    BoundingBox(&data[k],k,nomeshes,info);

  if(info) for(k=0;k<nomeshes;k++)
    MeshPieces(&data[k],k,nomeshes,info);

  if(eg.nosave) {
    Goodbye();
    return(0);
  }


  /********************************/
  printf("\nElmergrid saving data with method %d:\n",outmethod);
  printf(  "-------------------------------------\n");



  switch (outmethod) {
  case 1:
    SaveElmergrid(grids,nogrids,eg.filesout[0],info);
    break; 

  case 2:
    for(k=0;k<nomeshes;k++) {
      if(data[k].nopartitions > 1) 
	SaveElmerInputPartitioned(&data[k],boundaries[k],eg.filesout[k],eg.decimals,
				  eg.parthalo,eg.partitionindirect,eg.parthypre,
				  MAX(eg.partbcz,eg.partbcr),eg.nooverwrite,info);
      else
	SaveElmerInput(&data[k],boundaries[k],eg.filesout[k],eg.decimals,eg.nooverwrite,info);
    }
    break;

  case 22:

    for(k=0;k<nomeshes;k++) {
      SaveElmerInputFemBem(&data[k],boundaries[k],eg.filesout[k],eg.decimals,info);
    }
    break;


  case 3: 
      /* Create a variable so that when saving data in ElmerPost format there is something to visualize */
    for(k=0;k<nomeshes;k++) {
      if(data[k].variables == 0) {
	CreateVariable(&data[k],1,1,0.0,"Number",FALSE);
	for(i=1;i<=data[k].alldofs[1];i++)
	  data[k].dofs[1][i] = (Real)(i);	
      }
      SaveSolutionElmer(&data[k],boundaries[k],eg.saveboundaries ? MAXBOUNDARIES:0,
			eg.filesout[k],eg.decimals,info);
    }
    break;

  case 4:
    for(k=0;k<nomeshes;k++) {
      SaveMeshGmsh(&data[k],boundaries[k],eg.saveboundaries ? MAXBOUNDARIES:0,
		   eg.filesout[k],eg.decimals,info);
    }
    break;

  case 5:
    for(k=0;k<nomeshes;k++) {
      SaveMeshVtu(&data[k],boundaries[k],eg.saveboundaries ? MAXBOUNDARIES:0,
		   eg.filesout[k],eg.vtuone,info);
    }
    break;

  case 6:
    for(k=0;k<nomeshes;k++)
      SaveAbaqusInput(&data[k],eg.filesout[k],info); 
    break;
    
  case 7:
    for(k=0;k<nomeshes;k++)
      SaveFidapOutput(&data[k],eg.filesout[k],info,1,data[k].dofs[1]);
    break;


    /* Some obsolete special formats related to mapping, view factors etc. */

  case 101:
    for(k=0;k<nogrids;k++) {
      for(i=0;i<grids[k].noboundaries;i++)
	if(boundaries[k][i].created == TRUE) {
	  sprintf(prefix,"%s%d",eg.filesout[k],i+1);
	  SaveBoundary(&data[k],&boundaries[k][i],prefix,info);
	}
    }
    break;
    
  default:
    Instructions();
    break;
  }    

  timer_show();

  Goodbye();
  return(0);
}
