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

const char *IOmethods[] = {
  /*0*/ "EG",
  /*1*/ "ELMERGRID",
  /*2*/ "ELMERSOLVER",
  /*3*/ "ELMERPOST",
  /*4*/ "ANSYS",
  /*5*/ "IDEAS",
  /*6*/ "NASTRAN",
  /*7*/ "FIDAP",
  /*8*/ "UNV",
  /*9*/ "COMSOL",
  /*10*/ "FIELDVIEW",
  /*11*/ "TRIANGLE",
  /*12*/ "MEDIT",
  /*13*/ "GID",
  /*14*/ "GMSH",
  /*15*/ "PARTITIONED",
  /*16*/ "CGSIM",
};



int InlineParameters(struct ElmergridType *eg,int argc,char *argv[],const char *IOmethods[],int first,int info)
{
  int arg,i,dim;
  char command[MAXLINESIZE];
  
  dim = eg->dim;

  /* Type of input file */
  if(first > 3) {
    for(i=0;i<MAXLINESIZE;i++) command[i] = ' ';

    strcpy(command,argv[1]);
    for(i=0;i<MAXLINESIZE;i++) command[i] = toupper(command[i]);
    for(i=0;i<=MAXFORMATS;i++) {
      if(strstr(command,IOmethods[i])) {
	eg->inmethod = i;
	break;
      }
    }
    if(i>MAXFORMATS) eg->inmethod = atoi(argv[1]);
    

    /* Type of output file (fewer options) */
    strcpy(command,argv[2]);
    for(i=0;i<MAXLINESIZE;i++) command[i] = toupper(command[i]);
    for(i=1;i<=MAXFORMATS;i++) {
      if(strstr(command,IOmethods[i])) {
	eg->outmethod = i;
	break;
      }
    }
    if(i>MAXFORMATS) eg->outmethod = atoi(argv[2]);
 
    /* Name of output file */
    strcpy(eg->filesin[0],argv[3]);
    strcpy(eg->filesout[0],eg->filesin[0]);
    strcpy(eg->mapfile,eg->filesin[0]);
  }


  /* The optional inline parameters */

  for(arg=first;arg <argc; arg++) {

    if(strcmp(argv[arg],"-silent") == 0) {
      eg->silent = TRUE;
      info = FALSE;
    }

    if(strcmp(argv[arg],"-in") ==0 ) {
      if(arg+1 >= argc) {
	printf("The secondary input file name is required as a parameter\n");
	return(1);
      }
      else {
	strcpy(eg->filesin[eg->nofilesin],argv[arg+1]);
	printf("A secondary input file %s will be loaded.\n",eg->filesin[eg->nofilesin]);
	eg->nofilesin++;
      }
    }

    if(strcmp(argv[arg],"-out") == 0) {
      if(arg+1 >= argc) {
	printf("The output name is required as a parameter\n");
	return(2);
      }
      else {
	strcpy(eg->filesout[0],argv[arg+1]);
      }
    }

    if(strcmp(argv[arg],"-decimals") == 0) {
      eg->decimals = atoi(argv[arg+1]);
    }

    if(strcmp(argv[arg],"-triangles") ==0) {
      eg->triangles = TRUE;
      printf("The rectangles will be split to triangles.\n");
      if(arg+1 < argc) {
	if(strcmp(argv[arg+1],"-")) {
	  eg->triangleangle = atof(argv[arg+1]);
	}
      }
    }

    if(strcmp(argv[arg],"-merge") == 0) {
      if(arg+1 >= argc) {
	printf("Give a parameter for critical distance.\n");
	return(3);
      }
      else {
	eg->merge = TRUE;
	eg->cmerge = atof(argv[arg+1]);
      }
    }

    if(strcmp(argv[arg],"-relh") == 0) {
      if(arg+1 >= argc) {
	printf("Give a relative mesh density related to the specifications\n");
	return(3);
      }
      else {
	eg->relh = atof(argv[arg+1]);
      }
    }

    if(strcmp(argv[arg],"-order") == 0) {
      if(arg+dim >= argc) {
	printf("Give %d parameters for the order vector.\n",dim);
 	return(4);
      }
      else {
	eg->order = TRUE;
	eg->corder[0] = atof(argv[arg+1]);
	eg->corder[1] = atof(argv[arg+2]);
	if(dim==3) eg->corder[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-autoorder") == 0) {
      eg->order = 2;
    }

    if(strcmp(argv[arg],"-halo") == 0) {
      eg->partitionhalo = TRUE;
    }
    if(strcmp(argv[arg],"-indirect") == 0) {
      eg->partitionindirect = TRUE;
    }
    if(strcmp(argv[arg],"-metisorder") == 0) {
      eg->order = 3;
    }
    if(strcmp(argv[arg],"-centralize") == 0) {
      eg->center = TRUE;
    }
    if(strcmp(argv[arg],"-scale") == 0) {
      if(arg+dim >= argc) {
	printf("Give %d parameters for the scaling.\n",dim);
 	return(5);
     }
      else {
	eg->scale = TRUE;
	eg->cscale[0] = atof(argv[arg+1]);
	eg->cscale[1] = atof(argv[arg+2]);
	if(dim==3) eg->cscale[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-translate") == 0) {
      if(arg+dim >= argc) {
	printf("Give %d parameters for the translate vector.\n",dim);
	return(6);
      }
      else {
	eg->translate = TRUE;
	eg->ctranslate[0] = atof(argv[arg+1]);
	eg->ctranslate[1] = atof(argv[arg+2]);
	if(dim == 3) eg->ctranslate[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-saveinterval") == 0) {
      if(arg+dim >= argc) {
	printf("Give min, max and step for the interval.\n");
	return(7);
      }
      else {
	eg->saveinterval[0] = atoi(argv[arg+1]);
	eg->saveinterval[1] = atoi(argv[arg+2]);
	eg->saveinterval[2] = atoi(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-rotate") == 0 || strcmp(argv[arg],"-rotate") == 0) {
      if(arg+dim >= argc) {
	printf("Give three parameters for the rotation angles.\n");
	return(8);
      }
      else {
	eg->rotate = TRUE;
	eg->crotate[0] = atof(argv[arg+1]);
	eg->crotate[1] = atof(argv[arg+2]);
	eg->crotate[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-clone") == 0) {
      if(arg+dim >= argc) {
	printf("Give the number of clones in each %d directions.\n",dim);
 	return(9);
     }
      else {
	eg->clone[0] = atoi(argv[arg+1]);
	eg->clone[1] = atoi(argv[arg+2]);
	if(dim == 3) eg->clone[2] = atoi(argv[arg+3]);
      }
    }
    if(strcmp(argv[arg],"-clonesize") == 0) {
      if(arg+dim >= argc) {
	printf("Give the clone size in each %d directions.\n",dim);
 	return(10);
      }
      else {
	eg->clonesize[0] = atof(argv[arg+1]);
	eg->clonesize[1] = atof(argv[arg+2]);
	if(dim == 3) eg->clonesize[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-unite") == 0) {
      eg->unitemeshes = TRUE;
      printf("The meshes will be united.\n");
    }   

    if(strcmp(argv[arg],"-names") == 0) {
      eg->usenames = TRUE;
      printf("Names will be conserved when possible\n");
    }   

    if(strcmp(argv[arg],"-removelowdim") == 0) {
      eg->removelowdim = TRUE;
      printf("Lower dimensional boundaries will be removed\n");
    }   

    if(strcmp(argv[arg],"-removeunused") == 0) {
      eg->removeunused = TRUE;
      printf("Nodes that do not appear in any element will be removed\n");
    }   

    if(strcmp(argv[arg],"-autoclean") == 0) {
      eg->removelowdim = TRUE;
      eg->bulkorder = TRUE;
      eg->boundorder = TRUE;
      eg->removeunused = TRUE;
      printf("Lower dimensional boundaries will be removed\n");
      printf("Materials and boundaries will be renumbered\n");
      printf("Nodes that do not appear in any element will be removed\n");
    }   

    if(strcmp(argv[arg],"-polar") == 0) {
      eg->polar = TRUE;
      printf("Making transformation to polar coordinates.\n");
      if(arg+1 >= argc) {
	printf("The preferred radius is required as a parameter\n");
	eg->polarradius = 1.0;
      }
      else {
	eg->polarradius = atoi(argv[arg+1]);
      }
    }

    if(strcmp(argv[arg],"-cylinder") == 0) {
      eg->cylinder = TRUE;
      printf("Making transformation from cylindrical to cartesian coordinates.\n");
    }

    if(strcmp(argv[arg],"-reduce") == 0) {
      if(arg+2 >= argc) {
	printf("Give two material for the interval.\n");
 	return(12);
      }
      else {
	eg->reduce = TRUE;      
	eg->reducemat1 = atoi(argv[arg+1]);
	eg->reducemat2 = atoi(argv[arg+2]);
      }
    }
    if(strcmp(argv[arg],"-increase") == 0) {
      eg->increase = TRUE;
    }
    if(strcmp(argv[arg],"-bulkorder") == 0) {
      eg->bulkorder = TRUE;
    }
    if(strcmp(argv[arg],"-boundorder") == 0) {
      eg->boundorder = TRUE;
    }
    if(strcmp(argv[arg],"-pelem") == 0) {
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) 
	eg->pelemmap[3*eg->pelems+i-1-arg] = atoi(argv[i]);
      eg->pelems++;
    } 
    if(strcmp(argv[arg],"-belem") == 0) {
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) 
	eg->belemmap[3*eg->belems+i-1-arg] = atoi(argv[i]);
      eg->belems++;
    } 
    if(strcmp(argv[arg],"-partition") == 0) {
      if(arg+dim >= argc) {
	printf("The number of partitions in %d dims is required as parameters.\n",dim);
	return(13);
      }
      else {
	eg->partitions = 1;
	eg->partdim[0] = atoi(argv[arg+1]);
	eg->partdim[1] = atoi(argv[arg+2]);
	if(dim == 3) eg->partdim[2] = atoi(argv[arg+3]);
	eg->partitions = 1;
	for(i=0;i<3;i++) {
	  if(eg->partdim[i] == 0) eg->partdim[i] = 1;
	  eg->partitions *= eg->partdim[i];
	}
	eg->partopt = 0;
	if(arg+4 < argc) 
	  if(argv[arg+4][0] != '-') eg->partopt = atoi(argv[arg+4]);

	printf("The mesh will be partitioned with simple division to %d partitions.\n",
	       eg->partitions);
      }
    }
    if(strcmp(argv[arg],"-partorder") == 0) {
      if(arg+dim >= argc) {
	printf("Give %d parameters for the order vector.\n",dim);
 	return(14);
      }
      else {
	eg->partorder = 1;
	eg->partcorder[0] = atof(argv[arg+1]);
	eg->partcorder[1] = atof(argv[arg+2]);
	if(dim==3) eg->partcorder[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-metis") == 0) {
#if HAVE_METIS
      if(arg+1 >= argc) {
	printf("The number of partitions is required as a parameter\n");
	return(15);
      }
      else {
	eg->metis = atoi(argv[arg+1]);
	printf("The mesh will be partitioned with Metis to %d partitions.\n",eg->metis);
	eg->partopt = 0;
	if(arg+2 < argc) 
	  if(argv[arg+2][0] != '-') eg->partopt = atoi(argv[arg+2]);
      }
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif     
    }

    if(strcmp(argv[arg],"-periodic") == 0) {
      if(arg+dim >= argc) {
	printf("Give the periodic coordinate directions (e.g. 1 1 0)\n");
 	return(16);
      }
      else {
	eg->periodicdim[0] = atoi(argv[arg+1]);
	eg->periodicdim[1] = atoi(argv[arg+2]);
	if(dim == 3) eg->periodicdim[2] = atoi(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-discont") == 0) {
      if(arg+1 >= argc) {
	printf("Give the discontinuous boundary conditions.\n");
 	return(17);
      }
      else {
	eg->discontbounds[eg->discont] = atoi(argv[arg+1]);
	eg->discont++;
      }
    }

    if(strcmp(argv[arg],"-connect") == 0) {
      if(arg+1 >= argc) {
	printf("Give the connected boundary conditions.\n");
 	return(10);
      }
      else {
	eg->connectbounds[eg->connect] = atoi(argv[arg+1]);
	eg->connect++;
      }
    } 
 
    if(strcmp(argv[arg],"-boundbound") == 0) {
      for(i=arg+1;i<=arg+3 && i<argc; i++) {
	eg->boundbound[3*eg->boundbounds+i-(1+arg)] = atoi(argv[i]);
	if((i-arg)%3 == 0) eg->boundbounds++;
      }
    } 
    if(strcmp(argv[arg],"-bulkbound") == 0) {
      for(i=arg+1;i<=arg+3 && i<argc; i++) {
	eg->bulkbound[3*eg->bulkbounds+i-(1+arg)] = atoi(argv[i]);
	if((i-arg)%3 == 0) eg->bulkbounds++;
      }
    } 
    if(strcmp(argv[arg],"-boundtype") == 0) {
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) 
	eg->sidemap[3*eg->sidemappings+i-1-arg] = atoi(argv[i]);
      eg->sidemappings++;
    } 
    if(strcmp(argv[arg],"-bulktype") == 0) {
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) 
	eg->bulkmap[3*eg->bulkmappings+i-1-arg] = atoi(argv[i]);
      eg->bulkmappings++;
    } 

    if(strcmp(argv[arg],"-layer") == 0) {
      if(arg+4 >= argc) {
	printf("Give four parameters for the layer: boundary, elements, thickness, ratio.\n");
	return(18);
      }
      else if(eg->layers == MAXBOUNDARIES) {
	printf("There can only be %d layers, sorry.\n",MAXBOUNDARIES);
	return(19);
      }
      else {
	eg->layerbounds[eg->layers] = atoi(argv[arg+1]);
	eg->layernumber[eg->layers] = atoi(argv[arg+2]);
	eg->layerthickness[eg->layers] = atof(argv[arg+3]);
	eg->layerratios[eg->layers] = atof(argv[arg+4]);
	eg->layerparents[eg->layers] = 0;
	eg->layers++;
      }
    }
    
    if(strcmp(argv[arg],"-layermove") == 0) {
      if(arg+1 >= argc) {
	printf("Give maximum number of Jacobi filters.\n");
 	return(20);
      }
      else {
	eg->layermove = atoi(argv[arg+1]);
      }
    }

    /* This uses a very dirty trick where the variables related to argument -layer are used 
       with a negative indexing */ 
    if(strcmp(argv[arg],"-divlayer") == 0) {
      if(arg+4 >= argc) {
	printf("Give four parameters for the layer: boundary, elements, relative thickness, ratio.\n");
	return(21);
      }
      else if(abs(eg->layers) == MAXBOUNDARIES) {
	printf("There can only be %d layers, sorry.\n",MAXBOUNDARIES);
	return(22);
      }
      else {
	eg->layerbounds[abs(eg->layers)] = atoi(argv[arg+1]);
	eg->layernumber[abs(eg->layers)] = atoi(argv[arg+2]);
	eg->layerthickness[abs(eg->layers)] = atof(argv[arg+3]);
	eg->layerratios[abs(eg->layers)] = atof(argv[arg+4]);
	eg->layerparents[abs(eg->layers)] = 0;
	eg->layers--;
      }
    }

    if(strcmp(argv[arg],"-3d") == 0) {
      eg->dim = dim = 3;
    }
    if(strcmp(argv[arg],"-2d") == 0) {
      eg->dim = dim = 2;
    }
    if(strcmp(argv[arg],"-1d") == 0) {
      eg->dim = dim = 1;
    }
    
    if(strcmp(argv[arg],"-isoparam") == 0) {
      eg->isoparam = TRUE;
    }
    if(strcmp(argv[arg],"-nobound") == 0) {
      eg->saveboundaries = FALSE;
    }
    
    /* The following keywords are not actively used */
    
    if(strcmp(argv[arg],"-map") ==0) {
      if(arg+1 >= argc) {
	printf("Give the name of the mapping file\n");
	return(23);
      }
      else {
	strcpy(eg->mapfile,argv[arg+1]);
	printf("Mapping file is %s\n",eg->mapfile);
      }
    }
    if(strcmp(argv[arg],"-bcoffset") == 0) {
      eg->bcoffset = atoi(argv[arg+1]);
    }
    if(strcmp(argv[arg],"-noelements") == 0) {
      eg->elements3d = atoi(argv[arg+1]);
    }
    if(strcmp(argv[arg],"-nonodes") == 0) {
      eg->nodes3d = atoi(argv[arg+1]);
    }
    
    if(strcmp(argv[arg],"-sidefind") == 0) {
      eg->findsides = 0;
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) {
	eg->sidebulk[i-1-arg] = atoi(argv[i]);
	eg->findsides++;
      }
    } 
    if(strcmp(argv[arg],"-findbound") == 0) {
      eg->findsides = 0;
      for(i=arg+1;i+1<argc && strncmp(argv[i],"-",1); i += 2) {
	eg->sidebulk[i-1-arg] = atoi(argv[i]);
	eg->sidebulk[i-arg] = atoi(argv[i+1]);
	eg->findsides++;
      }
    } 
  }

  {
    char *ptr1;
    ptr1 = strchr(eg->filesout[0], '.');
    if (ptr1) *ptr1 = '\0';
    ptr1 = strchr(eg->mapfile, '.');
    if (ptr1) *ptr1 = '\0';
  }

  return(0);
}


#if EXE_MODE
static void Goodbye()
{
  printf("\nThank you for using Elmergrid!\n");
  printf("Send bug reports and feature wishes to peter.raback@csc.fi\n");
  exit(0);
}

static void Instructions()
{
  printf("****************** Elmergrid ************************\n");
  printf("This program can create simple 2D structured meshes consisting of\n");
  printf("linear, quadratic or cubic rectangles or triangles. The meshes may\n");
  printf("also be extruded and revolved to create 3D forms. In addition many\n");
  printf("mesh formats may be imported into Elmer software. Some options have\n");
  printf("not been properly tested. Contact the author if you face problems.\n\n");

  printf("The program has two operation modes\n");
  printf("A) Command file mode which has the command file as the only argument\n");
  printf("   'ElmerGrid commandfile.eg'\n\n");

  printf("B) Inline mode which expects at least three input parameters\n");
  printf("   'ElmerGrid 1 3 test'\n\n");
  printf("The first parameter defines the input file format:\n");
  printf("1)  .grd      : Elmergrid file format\n");
  printf("2)  .mesh.*   : Elmer input format\n");
  printf("3)  .ep       : Elmer output format\n");
  printf("4)  .ansys    : Ansys input format\n");
  printf("5)  .inp      : Abaqus input format by Ideas\n");
  printf("6)  .msh      : Nastran format\n");
  printf("7)  .FDNEUT   : Gambit (Fidap) neutral file\n");
  printf("8)  .unv      : Universal mesh file format\n");
  printf("9)  .mphtxt   : Comsol Multiphysics mesh format\n");
  printf("10) .dat      : Fieldview format\n");
  printf("11) .node,.ele: Triangle 2D mesh format\n");
  printf("12) .mesh     : Medit mesh format\n");
  printf("13) .msh      : GID mesh format\n");
  printf("14) .msh      : Gmsh mesh format\n");
  printf("15) .ep.i     : Partitioned ElmerPost format\n");

  printf("\nThe second parameter defines the output file format:\n");
  printf("1)  .grd      : ElmerGrid file format\n");
  printf("2)  .mesh.*   : ElmerSolver format (also partitioned .part format)\n");
  printf("3)  .ep       : ElmerPost format\n");

  printf("\nThe third parameter is the name of the input file.\n");
  printf("If the file does not exist, an example with the same name is created.\n");
  printf("The default output file name is the same with a different suffix.\n\n");

  printf("There are several additional in-line parameters that are\n");
  printf("taken into account only when applicable to the given format.\n");

  printf("-out str             : name of the output file\n");
  printf("-in str              : name of a secondary input file\n");
  printf("-silent              : do not echo run time information\n");
  printf("-decimals            : number of decimals in the saved mesh (eg. 8)\n");
  printf("-triangles           : rectangles will be divided to triangles\n");
  printf("-relh real           : give relative mesh density parameter for ElmerGrid meshing\n");
  printf("-merge real          : merges nodes that are close to each other\n");
  printf("-order real[3]       : reorder elements and nodes using c1*x+c2*y+c3*z\n");
  printf("-centralize          : set the center of the mesh to origin\n");
  printf("-scale real[3]       : scale the coordinates with vector real[3]\n");
  printf("-translate real[3]   : translate the nodes with vector real[3]\n");
  printf("-rotate real[3]      : rotate around the main axis with angles real[3]\n");
  printf("-clone int[3]        : make ideantilcal copies of the mesh\n");
  printf("-clonesize real[3]   : the size of the mesh to be cloned if larger to the original\n");
  printf("-unite               : the meshes will be united\n");
  printf("-polar real          : map 2D mesh to a cylindrical shell with given radius\n");
  printf("-cylinder            : map 2D/3D cylindrical mesh to a cartesian mesh\n");
  printf("-reduce int[2]       : reduce element order at material interval [int1 int2]\n");
  printf("-increase            : increase element order from linear to quadratic\n");
  printf("-bcoffset int        : add an offset to the boundary conditions\n");
  printf("-discont int         : make the boundary to have secondary nodes\n");
  printf("-connect int         : make the boundary to have internal connection among its elements\n");
  printf("-removelowdim        : remove boundaries that are two ranks lower than highest dim\n");
  printf("-removeunused        : remove nodes that are not used in any element\n");
  printf("-bulkorder           : renumber materials types from 1 so that every number is used\n");
  printf("-boundorder          : renumber boundary types from 1 so that every number is used\n");
  printf("-autoclean           : this performs the united action of the three above\n");
  printf("-bulkbound int[3]    : set the union of materials [int1 int2] to be boundary int3\n");
  printf("-boundbound int[3]   : set the union of boundaries [int1 int2] to be boundary int3\n");
  printf("-bulktype int[3]     : set material types in interval [int1 int2] to type int3\n");
  printf("-boundtype int[3]    : set sidetypes in interval [int1 int2] to type int3\n");
  printf("-layer int[2] real[2]: make a boundary layer for given boundary\n");
  printf("-layermove int       : apply Jacobi filter int times to move the layered mesh\n");
  printf("-divlayer int[2] real[2]: make a boundary layer for given boundary\n");
  printf("-3d / -2d / -1d      : mesh is 3, 2 or 1-dimensional (applies to examples)\n");
  printf("-isoparam            : ensure that higher order elements are convex\n");
  printf("-nobound             : disable saving of boundary elements in ElmerPost format\n");

  printf("\nThe following keywords are related only to the parallel Elmer computations.\n");
  printf("-partition int[4]    : the mesh will be partitioned in main directions\n");
  printf("-partorder real[3]   : in the above method, the direction of the ordering\n");
#if HAVE_METIS
  printf("-metis int[2]        : the mesh will be partitioned with Metis\n");
#endif
  printf("-halo                : create halo for the partitioning\n");
  printf("-indirect            : create indirect connections in the partitioning\n");
  printf("-periodic int[3]     : decleare the periodic coordinate directions for parallel meshes\n");
  printf("-saveinterval int[3] : the first, last and step for fusing parallel data\n");

  if(0) printf("-names               : conserve name information where applicable\n");
}







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

  printf("Done converting\n");
  
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
    for(k=0;k<MAXCASES;k++) {
      boundaries[k] = (struct BoundaryType*)
	malloc((size_t) (MAXBOUNDARIES)*sizeof(struct BoundaryType)); 	
      for(i=0;i<MAXBOUNDARIES;i++) {
	boundaries[k][i].created = FALSE; 
	boundaries[k][i].nosides = 0;
      }    
    }
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
    if(*nogrids) LoadCommands(eg.filesin[nofile],&eg,grids,2,IOmethods,info); 
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
    printf("Loading Gmsh format file\n");
    errorstat = LoadGmshInput(&(data[nofile]),boundaries[nofile],eg.filesin[nofile],info);
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

  if(!errorstat) *nogrids = MAX(*nogrids,1);

  return(errorstat);
}





static int ManipulateMeshDefinition(int inmethod,int outmethod,Real relh)
{
  static int visited = FALSE;
  int i,j,k;
  Real mergeeps;

  printf("Manipulate mesh definition: %d %d %12.3le\n",inmethod,outmethod,relh);
  
  if(inmethod == 1 && outmethod != 1) {
    if(visited) {
      printf("Deallocating structures from previous call\n");
      for(k=0;k<MAXCASES;k++) {	    
	if(data[k].created) {
	  DestroyKnots(&data[k]);
	  for(i=0;i<MAXBOUNDARIES;i++) 
	    DestroyBoundary(&boundaries[k][i]);
	}
      }
    }
    printf("Starting to create ElmerGrid meshes\n");
    for(k=0;k<nogrids;k++) 
      CreateElmerGridMesh(&(grids[k]),&(data[k]),boundaries[k],relh,info);
    nomeshes = nogrids;
    printf("Created %d ElmerGrid meshes\n",nomeshes);
  }

  visited = TRUE;

  /* At first instance perform operations that should rather be done before extrusion 
     or mesh union. */
  for(k=0;k<nomeshes;k++) {
    
    /* Make the discontinous boundary needed, for example, in poor thermal conduction */
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
    for(i=1;i<=eg.connect;i++) 
      SetConnectedBoundary(&(data[k]),boundaries[k],eg.connectbounds[i-1],i,info);
  
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
      UniteMeshes(&data[0],&data[k],boundaries[0],boundaries[k],info);
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
  static char arguments[10][10],**argv;
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
    errorstat = LoadCommands(filename,&eg,grids,1,IOmethods,info);
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

  errorstat = InlineParameters(&eg,argc,argv,IOmethods,0,info);
  if(errorstat) printf("Errorstat for inline parameters: %d\n",errorstat);
  
  inmethod = eg.inmethod;
  outmethod = 0;

  printf("Input method: %d\n",inmethod);
  
  ManipulateMeshDefinition(inmethod,outmethod,eg.relh);

  errorstat = ConvertEgTypeToMeshType(&data[activemesh],boundaries[activemesh],mesh);

  if(info) printf("Done converting mesh\n");
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

  if(info) printf("\nStarting program Elmergrid\n");
  
  InitParameters(&eg);
  grids = (struct GridType*)malloc((size_t) (MAXCASES)*sizeof(struct GridType));     
  InitGrid(grids);
  info = TRUE;

  if(argc <= 2) {
    errorstat = LoadCommands(argv[1],&eg,grids,argc-1,IOmethods,info);     
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
    errorstat = InlineParameters(&eg,argc,argv,IOmethods,4,info);
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



