/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author:  Peter Råback
   Email:   Peter.Raback@csc.fi
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

/* --------------------:  feminfo.c  :--------------------------

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

#include "common.h"
#include "nrutil.h"
#include "femdef.h"
#include "femtypes.h"
#include "femmesh.h"
#include "femknot.h"
#include "feminfo.h"
#include "femelmer.h"

int matcactive=FALSE, iodebug=FALSE;

int MAXINMETHODS = 15;
char *InMethods[] = {
  /*0*/ "EG",
  /*1*/ "ELMERGRID",
  /*2*/ "ELMERSOLVER",
  /*3*/ "ELMERPOST",
  /*4*/ "ANSYS",
  /*5*/ "IDEAS",
  /*6*/ "ABAQUS",
  /*7*/ "FIDAP",
  /*8*/ "UNV",
  /*9*/ "COMSOL",
  /*10*/ "FIELDVIEW",
  /*11*/ "TRIANGLE",
  /*12*/ "MEDIT",
  /*13*/ "GID",
  /*14*/ "GMSH",
  /*15*/ "PARTITIONED",
};

int MAXOUTMETHODS = 5;
char *OutMethods[] = {
  /*0*/ "EG",
  /*1*/ "ELMERGRID",
  /*2*/ "ELMERSOLVER",
  /*3*/ "ELMERPOST",
  /*4*/ "GMSH",
  /*5*/ "VTU",
};



#ifndef DISABLE_MATC
char *mtc_domath(const char *);
void mtc_init(FILE * input, FILE * output, FILE * error);
#endif

static int Getline(char *line1,FILE *io) 
{
  int i,isend;
  char line0[MAXLINESIZE],*charend,*matcpntr,*matcpntr0;

  for(i=0;i<MAXLINESIZE;i++) 
    line0[i] = ' ';

 newline:

  charend = fgets(line0,MAXLINESIZE,io);
  isend = (charend == NULL);

  if(isend) return(1);

  if(line0[0] == '#' || line0[0] == '%' || line0[0] == '!') goto newline;
  if(!matcactive && line0[0] == '*') goto newline;

#ifndef DISABLE_MATC
  if(matcactive) {
    matcpntr0 = strchr(line0,'$');
    if(matcpntr0) {
      matcpntr = mtc_domath(&matcpntr0[1]);
      if(matcpntr) {
	strcpy(matcpntr0, matcpntr);
	if(0) printf("A: %s\n%s\n",matcpntr0,matcpntr);
      }
    }
  }
#endif 

  if(strstr(line0,"subcell boundaries")) goto newline;
  if(strstr(line0,"material structure")) goto newline;
  if(strstr(line0,"mode")) goto newline;
  if(strstr(line0,"type")) goto newline;

  for(i=0;i<MAXLINESIZE;i++) 
    line1[i] = toupper(line0[i]);

  if(iodebug) {
    printf("line: ");
    for(i=0;i<40;i++) printf("%c",line1[i]);
    printf("\n");
  }

  return(0);
}


static int GetCommand(char *line1,char *line2,FILE *io) 
/* Line1 for commands and line2 for arguments. */
{
  int i,j,isend;
  char line0[MAXLINESIZE],*charend,*matcpntr0,*matcpntr;
  int gotlinefeed;

 newline:

  for(i=0;i<MAXLINESIZE;i++) 
    line2[i] = line1[i] = line0[i] = ' ';

  charend = fgets(line0,MAXLINESIZE,io);
  isend = (charend == NULL);

  if(isend) return(1);

  if(strlen(line0)<=strspn(line0," \t\r\n")) goto newline;

  if(line0[0] == '#' || line0[0] == '%' || line0[0] == '!' || line0[0] == '\n') goto newline;
  if(!matcactive && line0[0] == '*') goto newline;

#ifndef DISABLE_MATC
  if(matcactive) {
    matcpntr0 = strchr(line0,'$');
    if(matcpntr0) {
      matcpntr = mtc_domath(&matcpntr0[1]);
      if(matcpntr) {
	strcpy(matcpntr0, matcpntr);
	if(0) printf("B: %s\n%s\n",matcpntr0,matcpntr);
      }
      else {
	if(0) printf("B0: %s\n",matcpntr0);
	goto newline;
      }
    }
  }
#endif 
  
  gotlinefeed = FALSE;
  j = 0;
  for(i=0;i<MAXLINESIZE;i++) {
    if(line0[i] == '\n' ) {
      gotlinefeed = TRUE;
      break;
    }
    if(line0[i] == '=') {
      j = i;
      break;
    }
    line1[i] = toupper(line0[i]);
  }

  /* After these commands there will be no nextline even though there is no equality sign */
  if(strstr(line1,"END")) return(0);
  if(strstr(line1,"NEW MESH")) return(0);

  if(j) { /* Arguments are actually on the same line after '=' */
    for(i=j+1;i<MAXLINESIZE;i++) {
      line2[i-j-1] = line0[i];    
      if( line0[i] == '\n' ) {
	gotlinefeed = TRUE;
	break;
      }
    }  
    if(!gotlinefeed) {
      printf("There is a risk that somethings was missed in line:\n");
      printf("%s\n",line0);
      smallerror("Check your output line length!\n");
    }
  }
  else { /* rguments are on the next line */
  newline2:
    charend = fgets(line2,MAXLINESIZE,io);
    isend = (charend == NULL);
    if(isend) return(2);
    if(line2[0] == '#' || line2[0] == '%' || line2[0] == '!') goto newline2;
    if(!matcactive && line2[0] == '*') goto newline2;

    gotlinefeed = FALSE;
    for(i=0;i<MAXLINESIZE;i++) {
      if(line2[i] == '\n' ) {
	gotlinefeed = TRUE;
	break;
      }
    }
    if(!gotlinefeed) {
      printf("There is a risk that somethings was missed in line:\n");
      printf("%s\n",line2);
      smallerror("Check your output line length!\n");
    }

#ifndef DISABLE_MATC
    if(matcactive) {
      matcpntr0 = strchr(line2,'$');
      if(matcpntr0) {
	matcpntr = mtc_domath(&matcpntr0[1]);
	if(matcpntr) {
	  strcpy(matcpntr0, matcpntr);
	  if(0) printf("C: %s\n%s\n",matcpntr0,matcpntr);
	}
      }
    }
#endif 
  }

  if(iodebug) {
    printf("command: ");
    for(i=0;i<40;i++) printf("%c",line1[i]);
    printf("\nparams: ");
    for(i=0;i<40;i++) printf("%c",line2[i]);
    printf("\n");
  }
  
  return(0);
}



void InitParameters(struct ElmergridType *eg)
{
  int i;

  eg->relh = 1.0;
  eg->inmethod = 0;
  eg->outmethod = 0;
  eg->nofilesin = 1;
  eg->unitemeshes = FALSE;
  eg->triangles = FALSE;
  eg->triangleangle = 0.0;
  eg->rotate = FALSE;
  eg->polar = FALSE;
  eg->cylinder = FALSE;
  eg->usenames = TRUE;
  eg->layers = 0;
  eg->layereps = 0.0;
  eg->layermove = 0;
  eg->partitions = 0;
  eg->elements3d = 0;
  eg->nodes3d = 0;
  eg->metis = 0;
  eg->metiscontig = FALSE;
  eg->partopt = 0;
  eg->partoptim = FALSE;
  eg->partbcoptim = TRUE;
  eg->partjoin = 0;
  eg->partitionhalo = 0;
  eg->partitionindirect = FALSE;
  eg->reduce = FALSE;
  eg->increase = FALSE;
  eg->translate = FALSE;
  eg->isoparam = FALSE;
  eg->removelowdim = FALSE;
  eg->removeintbcs = FALSE;
  eg->removeunused = FALSE;
  eg->dim = 3;
  eg->center = FALSE;
  eg->scale = FALSE;
  eg->order = FALSE;
  eg->boundbounds = 0;
  eg->saveinterval[0] = eg->saveinterval[1] = eg->saveinterval[2] = 0;
  eg->bulkbounds = 0;
  eg->partorder = FALSE;
  eg->findsides = FALSE;
  eg->parthypre = FALSE;
  eg->partdual = FALSE;
  eg->partbcz = 0;
  eg->partbcr = 0;
  eg->partbclayers = 1;
  eg->partbcmetis = 0;
  eg->partbw = FALSE;
  eg->saveboundaries = TRUE;
  eg->vtuone = FALSE;
  eg->timeron = FALSE;
  eg->nosave = FALSE;
  eg->nooverwrite = FALSE;
  eg->merge = FALSE;
  eg->bcoffset = FALSE;
  eg->periodic = 0;
  eg->periodicdim[0] = 0;
  eg->periodicdim[1] = 0;
  eg->periodicdim[2] = 0;
  eg->bulkorder = FALSE;
  eg->boundorder = FALSE;
  eg->sidemappings = 0;
  eg->bulkmappings = 0;
  eg->coordinatemap[0] = eg->coordinatemap[1] = eg->coordinatemap[2] = 0;
  eg->clone[0] = eg->clone[1] = eg->clone[2] = 0;
  eg->mirror[0] = eg->mirror[1] = eg->mirror[2] = 0;
  eg->cloneinds = FALSE;
  eg->mirrorbc = 0;
  eg->decimals = 12;
  eg->discont = 0;
  eg->connect = 0;
  eg->connectboundsnosets = 0;

  eg->rotatecurve = FALSE;
  eg->curverad = 0.5;
  eg->curveangle = 90.0;
  eg->curvezet = 0.0;
  
  for(i=0;i<MAXSIDEBULK;i++) 
    eg->sidebulk[i] = 0;
}




int InlineParameters(struct ElmergridType *eg,int argc,char *argv[])
{
  int arg,i,dim;
  char command[MAXLINESIZE];
  
  dim = eg->dim;

  printf("Elmergrid reading in-line arguments\n");

  /* Type of input file */
  strcpy(command,argv[1]);
  for(i=0;i<MAXLINESIZE;i++) command[i] = toupper(command[i]);
  for(i=0;i<=MAXINMETHODS;i++) {
    if(strstr(command,InMethods[i])) {
      eg->inmethod = i;
      break;
    }
  }
  if(i>MAXINMETHODS) eg->inmethod = atoi(argv[1]);


  /* Type of output file (fewer options) */
  strcpy(command,argv[2]);
  for(i=0;i<MAXLINESIZE;i++) command[i] = toupper(command[i]);
  for(i=1;i<=MAXOUTMETHODS;i++) {
    if(strstr(command,OutMethods[i])) {
      eg->outmethod = i;
      break;
    }
  }
  if(i>MAXOUTMETHODS) eg->outmethod = atoi(argv[2]);
 

  /* Default names of output file are derived from input file name */
  strcpy(eg->filesin[0],argv[3]);
  strcpy(eg->filesout[0],eg->filesin[0]);
  strcpy(eg->infofile,eg->filesin[0]);


  /* The optional inline parameters */

  for(arg=4;arg <argc; arg++) {



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
	eg->corder[2] = atof(argv[arg+3]);
      }
    }

    if(strcmp(argv[arg],"-autoorder") == 0) {
      eg->order = 2;
    }

    if(strcmp(argv[arg],"-halo") == 0) {
      eg->partitionhalo = 1;
    }
    if(strcmp(argv[arg],"-halobc") == 0) {
      eg->partitionhalo = 2;
    }
    if(strcmp(argv[arg],"-haloz") == 0) {
      eg->partitionhalo = 3;
    }
    if(strcmp(argv[arg],"-halor") == 0) {
      eg->partitionhalo = 3;
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
	eg->cscale[2] = atof(argv[arg+3]);
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
	eg->ctranslate[2] = atof(argv[arg+3]);
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
    if(strcmp(argv[arg],"-cloneinds") == 0) {
      eg->cloneinds = TRUE;
    }
    if(strcmp(argv[arg],"-mirror") == 0) {
      if(arg+dim >= argc) {
	printf("Give the symmetry of the coordinate directions, eg. 1 1 0\n");
      }
      else {
	eg->mirror[0] = atoi(argv[arg+1]);
	eg->mirror[1] = atoi(argv[arg+2]);
	if(dim == 3) eg->mirror[2] = atoi(argv[arg+3]);
      }
    }
    if(strcmp(argv[arg],"-mirrorbc") == 0) {
      if(arg+1 >= argc) {
	printf("Give the number of symmetry BC.\n");
 	return(11);
      }
      else {
	eg->mirrorbc = atoi(argv[arg+1]);
      }
    }

    if(strcmp(argv[arg],"-unite") == 0) {
      eg->unitemeshes = TRUE;
      printf("The meshes will be united.\n");
    }   

    if(strcmp(argv[arg],"-nonames") == 0) {
      eg->usenames = FALSE;
      printf("Names will be omitted even if they would exist\n");
    }   

    if(strcmp(argv[arg],"-removelowdim") == 0) {
      eg->removelowdim = TRUE;
      printf("Lower dimensional boundaries will be removed\n");
    }   

    if(strcmp(argv[arg],"-removeintbcs") == 0) {
      eg->removeintbcs = TRUE;
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
    if(strcmp(argv[arg],"-partition") == 0  ||
       strcmp(argv[arg],"-partcell") == 0  || 
       strcmp(argv[arg],"-partcyl") == 0 ) {
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
	eg->partopt = -1;
	if( strcmp(argv[arg],"-partition") == 0 ) {
	  if(arg+4 < argc) 
	    if(argv[arg+4][0] != '-') eg->partopt = atoi(argv[arg+4]);
	}
	else if( strcmp( argv[arg],"-partcell") == 0 )  {
	  eg->partopt = 2;
	} else if( strcmp( argv[arg],"-partcyl") == 0 ) {
	  eg->partopt = 3;
	}
	  
	printf("The mesh will be partitioned geometrically to %d partitions.\n",
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
    if(strcmp(argv[arg],"-partoptim") == 0) {
      eg->partoptim = TRUE;
      printf("Aggressive optimization will be applied to node sharing.\n");
    }
    if(strcmp(argv[arg],"-partnobcoptim") == 0) {
      eg->partbcoptim = FALSE;
      printf("Aggressive optimization will not be applied to parent element sharing.\n");
    }
    if(strcmp(argv[arg],"-partbw") == 0) {
      eg->partbw = TRUE;
      printf("Bandwidth will be optimized for partitions.\n");
    }
    if(strcmp(argv[arg],"-parthypre") == 0) {
      eg->parthypre = TRUE;
      printf("Numbering of partitions will be made continous.\n");
    }
    if(strcmp(argv[arg],"-partdual") == 0) {
      eg->partdual = TRUE;
      printf("Using dual (elemental) graph in partitioning.\n");
    }

    if(strcmp(argv[arg],"-metis") == 0 ||
       strcmp(argv[arg],"-metisrec") == 0 ||
       strcmp(argv[arg],"-metiskway") == 0 ) {
#if PARTMETIS
      if(arg+1 >= argc) {
	printf("The number of partitions is required as a parameter\n");
	return(15);
      }
      else {
	eg->metis = atoi(argv[arg+1]);
	printf("The mesh will be partitioned with Metis to %d partitions.\n",eg->metis);
	eg->partopt = 0;
	if(strcmp(argv[arg],"-metisrec") == 0)
	  eg->partopt = 2;
	else if(strcmp(argv[arg],"-metiskway") == 0 )
	  eg->partopt = 3;
	else if(arg+2 < argc) 
	  if(argv[arg+2][0] != '-') eg->partopt = atoi(argv[arg+2]);
      }    
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif     
    }

    if(strcmp(argv[arg],"-partjoin") == 0) {
      if(arg+1 >= argc) {
	printf("The number of partitions is required as a parameter!\n");
	return(15);
      }
      else {
	eg->partjoin = atoi(argv[arg+1]);
	printf("The results will joined using %d partitions.\n",eg->partjoin);
      }
    }

    if(strcmp(argv[arg],"-partconnect") == 0 || strcmp(argv[arg],"-partzbc") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of 1D partitions is required as a parameter!\n");
	return(15);
      }
      else {
	eg->partbcz = atoi(argv[arg+1]);
	printf("The connected BCs will be partitioned to %d partitions in Z.\n",eg->partbcz);
      }
    }

    if(strcmp(argv[arg],"-partrbc") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of 1D partitions is required as a parameter!\n");
	return(15);
      }
      else {
	eg->partbcr = atoi(argv[arg+1]);
	printf("The connected BCs will be partitioned to %d partitions in R.\n",eg->partbcr);
      }
    }

    if(strcmp(argv[arg],"-partlayers") == 0) {
      if(arg+1 >= argc) {
	printf("The number of layers to be extended is required as a parameter\n");
	return(15);
      }
      else {
	eg->partbclayers = atoi(argv[arg+1]);
	printf("The boundary partitioning will be extended by %d layers.\n",eg->partbclayers);
      }
    }

    if(strcmp(argv[arg],"-metiscontig") == 0 ) {
      eg->metiscontig = TRUE;
    }
    
    if(strcmp(argv[arg],"-metisconnect") == 0 || strcmp(argv[arg],"-metisbc") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of Metis partitions is required as a parameter\n");
	return(15);
      }
      else {
	eg->partbcmetis = atoi(argv[arg+1]);
	printf("The connected BCs will be partitioned to %d partitions by Metis.\n",eg->partbcmetis);
      }
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
	eg->connectboundsnosets += 1;
	for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) { 
	  eg->connectbounds[eg->connect] = atoi(argv[i]);
	  eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
	  eg->connect++;
	}
      }
    } 

    if(strcmp(argv[arg],"-connectall") == 0) {
      eg->connectboundsnosets += 1;
      eg->connectbounds[eg->connect] = -1;
      eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
      eg->connect++;
    }

    if(strcmp(argv[arg],"-connectint") == 0) {
      eg->connectboundsnosets += 1;
      eg->connectbounds[eg->connect] = -2;
      eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
      eg->connect++;
    }

    if(strcmp(argv[arg],"-connectfree") == 0) {
      eg->connectboundsnosets += 1;
      eg->connectbounds[eg->connect] = -3;
      eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
      eg->connect++;
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
    if(strcmp(argv[arg],"-coordinatemap") == 0) {
      if( arg+3 >= argc ) {
	printf("Give three parameters for the index permutation\n");
	return(18);
      }
      else {
	for(i=0;i<3;i++) 
	  eg->coordinatemap[i] = atoi(argv[arg+1+i]);
      }
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
    if(strcmp(argv[arg],"-vtuone") == 0) {
      eg->vtuone = TRUE;
    }
    if(strcmp(argv[arg],"-nosave") == 0) {
      eg->nosave = TRUE;
    }
    if(strcmp(argv[arg],"-nooverwrite") == 0) {
      eg->nooverwrite = TRUE;
    }
    if(strcmp(argv[arg],"-timer") == 0) {
      eg->timeron = TRUE;
    }

    if(strcmp(argv[arg],"-infofile") == 0) {
      eg->timeron = TRUE;
      if(arg+1 >= argc) {
	printf("The output name is required as a parameter\n");
	return(2);
      }
      else {
	strcpy(eg->infofile,argv[arg+1]);
      }
    }
    

    /* The following keywords are not actively used */

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
    int badpoint;
    char *ptr1,*ptr2;
    ptr1 = strrchr(eg->filesout[0], '.');
    if(ptr1) {
      badpoint=FALSE;
      ptr2 = strrchr(eg->filesout[0], '/');
      if(ptr2 && ptr2 > ptr1) badpoint = TRUE;
      if(!badpoint) *ptr1 = '\0';
    }
  }

  printf("Output will be saved to file %s.\n",eg->filesout[0]);

  return(0);
}




int LoadCommands(char *prefix,struct ElmergridType *eg,
		 struct GridType *grid, int mode,int info) 
{
  char filename[MAXFILESIZE];
  char command[MAXLINESIZE],params[MAXLINESIZE],*cp;

  FILE *in=NULL;
  int i,j;

  iodebug = FALSE;

  if( mode == 0) {  
    if (in = fopen("ELMERGRID_STARTINFO","r")) {
      fscanf(in,"%s",filename);
      fclose(in);
      printf("Using the file %s defined in ELMERGRID_STARTINFO\n",filename);
      if ((in = fopen(filename,"r")) == NULL) {
	printf("LoadCommands: opening of the file '%s' wasn't succesfull !\n",filename);
	return(1);
      }    
      else printf("Loading ElmerGrid commands from file '%s'.\n",filename);    
    }    
    else 
      return(2);
  }
  else if(mode == 1) { 
    AddExtension(prefix,filename,"eg");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadCommands: opening of the file '%s' wasn't succesfull !\n",filename);
      return(3);
    }    
    if(info) printf("Loading ElmerGrid commands from file '%s'.\n",filename);    
  }
  else if(mode == 2) {
    AddExtension(prefix,filename,"grd");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadCommands: opening of the file '%s' wasn't succesfull !\n",filename);
      return(4);
    }    
    if(info) printf("\nLoading ElmerGrid commands from file '%s'.\n",filename);
  }



  for(;;) {

    if(GetCommand(command,params,in)) {
      printf("Reached the end of command file\n");
      goto end;
    }    

    /* If the mode is the command file mode read also the file information from the command file. */

    if(mode <= 1) {
      if(strstr(command,"INPUT FILE")) {
	sscanf(params,"%s",eg->filesin[0]);
      }

      else if(strstr(command,"OUTPUT FILE")) {
	sscanf(params,"%s",eg->filesout[0]);
      }

      else if(strstr(command,"INPUT MODE")) {
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	
	for(i=0;i<=MAXINMETHODS;i++) {
	  if(strstr(params,InMethods[i])) {
	    eg->inmethod = i;
	    break;
	  }
	}
	if(i>MAXINMETHODS) sscanf(params,"%d",&eg->inmethod);
      }

      else if(strstr(command,"OUTPUT MODE")) {
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	
	/* Type of output file (fewer options) */
	for(i=1;i<=MAXOUTMETHODS;i++) {
	  if(strstr(params,OutMethods[i])) {
	    eg->outmethod = i;
	    break;
	  }
	}
	if(i>MAXOUTMETHODS) sscanf(params,"%d",&eg->outmethod);	
      }
    }    
    /* End of command file specific part */


    if(strstr(command,"DECIMALS")) {
      sscanf(params,"%d",&eg->decimals);
    }
    else if(strstr(command,"BOUNDARY OFFSET")) {
      sscanf(params,"%d",&eg->bcoffset);
    }
    else if(strstr(command,"TRIANGLES CRITICAL ANGLE")) {
      sscanf(params,"%le",&eg->triangleangle);
    }
    else if(strstr(command,"TRIANGLES")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->triangles = TRUE;      
    }
    else if(strstr(command,"MERGE NODES")) {
      eg->merge = TRUE;
      sscanf(params,"%le",&eg->cmerge);
    }
    else if(strstr(command,"UNITE")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->unitemeshes = TRUE;      
    }
    else if(strstr(command,"ORDER NODES")) {
      eg->order = TRUE;
      if(eg->dim == 1) 
	sscanf(params,"%le",&eg->corder[0]);
      else if(eg->dim == 2) 
	sscanf(params,"%le%le",&eg->corder[0],&eg->corder[1]);
      else if(eg->dim == 3) 
	sscanf(params,"%le%le%le",&eg->corder[0],&eg->corder[1],&eg->corder[2]);
    }
    else if(strstr(command,"SCALE")) {
      eg->scale = TRUE;
      if(eg->dim == 1) 
	sscanf(params,"%le",&eg->cscale[0]);
      else if(eg->dim == 2) 
	sscanf(params,"%le%le",&eg->cscale[0],&eg->cscale[1]);
      else if(eg->dim == 3) 
	sscanf(params,"%le%le%le",&eg->cscale[0],&eg->cscale[1],&eg->cscale[2]);
    }
    else if(strstr(command,"CENTRALIZE")) {
      eg->center = TRUE;
    }
    else if(strstr(command,"TRANSLATE")) {
      eg->translate = TRUE;
      if(eg->dim == 1) 
	sscanf(params,"%le",&eg->ctranslate[0]);
      else if(eg->dim == 2) 
	sscanf(params,"%le%le",&eg->ctranslate[0],&eg->ctranslate[1]);
      else if(eg->dim == 3) 
	sscanf(params,"%le%le%le",&eg->ctranslate[0],&eg->ctranslate[1],&eg->ctranslate[2]);
    }
    else if(strstr(command,"ROTATE MESH")) {
      eg->rotate = TRUE;
      sscanf(params,"%le%le%le",&eg->crotate[0],&eg->crotate[1],&eg->crotate[2]);
    }
    else if(strstr(command,"CLONE")) {
      if(strstr(command,"CLONE SIZE")) {
	if(eg->dim == 1) 
	  sscanf(params,"%le",&eg->clonesize[0]);
	else if(eg->dim == 2) 
	  sscanf(params,"%le%le",&eg->clonesize[0],&eg->clonesize[1]);
	else if(eg->dim == 3) 
	  sscanf(params,"%le%le%le",&eg->clonesize[0],&eg->clonesize[1],&eg->clonesize[2]);	
      }
      else {
	if(eg->dim == 1) 
	  sscanf(params,"%d",&eg->clone[0]);
	else if(eg->dim == 2) 
	  sscanf(params,"%d%d",&eg->clone[0],&eg->clone[1]);
	else if(eg->dim == 3) 
	  sscanf(params,"%d%d%d",&eg->clone[0],&eg->clone[1],&eg->clone[2]);
      }
    }
    else if(strstr(command,"MERGE")) {
      eg->merge = TRUE;
      sscanf(params,"%le",&eg->cmerge);
    }
    else if(strstr(command,"MIRROR")) {
      if(eg->dim == 1) 
	sscanf(params,"%d",&eg->mirror[0]);
      else if(eg->dim == 2) 
	sscanf(params,"%d%d",&eg->mirror[0],&eg->mirror[1]);
      else if(eg->dim == 3) 
	sscanf(params,"%d%d%d",&eg->mirror[0],&eg->mirror[1],&eg->mirror[2]);
    }
    else if(strstr(command,"POLAR RADIUS")) {
      eg->polar = TRUE;
      sscanf(params,"%le",&eg->polarradius);
    }
    else if(strstr(command,"CYLINDER")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->cylinder = TRUE;      
    }
    else if(strstr(command,"REDUCE DEGREE")) {
      eg->reduce = TRUE;
      sscanf(params,"%d%d",&eg->reducemat1,&eg->reducemat2);
    }
    else if(strstr(command,"INCREASE DEGREE")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->increase = TRUE;      
    }
    else if(strstr(command,"METIS OPTION")) {
#if PARTMETIS
      sscanf(params,"%d",&eg->partopt);
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"METIS")) {
#if PARTMETIS
      sscanf(params,"%d",&eg->metis);
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"METISKWAY")) {
#if PARTMETIS
      sscanf(params,"%d",&eg->metis);
      eg->partopt = 3;
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"METISREC")) {
#if PARTMETIS
      sscanf(params,"%d",&eg->metis);
      eg->partopt = 2;
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"PARTITION DUAL")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->partdual = TRUE;      
    }
    else if(strstr(command,"PARTJOIN")) {
      sscanf(params,"%d",&eg->partjoin);
    }
    else if(strstr(command,"PARTITION ORDER")) {
      eg->partorder = 1;
      if(eg->dim == 2) sscanf(params,"%le%le",&eg->partcorder[0],&eg->partcorder[1]);
      if(eg->dim == 3) sscanf(params,"%le%le%le",&eg->partcorder[0],
			      &eg->partcorder[1],&eg->partcorder[2]);      
    }
    else if(strstr(command,"PARTITION") || strstr(command,"PARTCYL") || strstr(command,"PARTCELL")) {
      if(eg->dim == 2) sscanf(params,"%d%d",&eg->partdim[0],&eg->partdim[1]);
      if(eg->dim == 3) sscanf(params,"%d%d%d",&eg->partdim[0],&eg->partdim[1],&eg->partdim[2]);
      eg->partitions = 1;
      for(i=0;i<eg->dim;i++) {
	if(eg->partdim[i] < 1) eg->partdim[i] = 1;
	eg->partitions *= eg->partdim[i];
      }
      if( strstr(command,"PARTCYL") ) eg->partopt = 3;
      if( strstr(command,"PARTCCELL") ) eg->partopt = 2;
    }
    else if(strstr(command,"PERIODIC")) {
      if(eg->dim == 2) sscanf(params,"%d%d",&eg->periodicdim[0],&eg->periodicdim[1]);
      if(eg->dim == 3) sscanf(params,"%d%d%d",&eg->periodicdim[0],
			      &eg->periodicdim[1],&eg->periodicdim[2]);
    }
    else if(strstr(command,"HALO")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->partitionhalo = 1;      
    }
    else if(strstr(command,"BOUNDARY HALO")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->partitionhalo = 2;
    }
    else if(strstr(command,"EXTRUDED HALO")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->partitionhalo = 3;
    }
    else if(strstr(command,"PARTBW")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->partbw = TRUE;      
    }
    else if(strstr(command,"PARTHYPRE")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->parthypre = TRUE;      
    }
    else if(strstr(command,"INDIRECT")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->partitionindirect = TRUE;      
    }
    else if(strstr(command,"BOUNDARY TYPE MAPPINGS")) {
      for(i=0;i<MAXMAPPINGS;i++) {
	if(i>0) Getline(params,in);
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	cp = params;      
	sscanf(params,"%d%d%d",&eg->sidemap[3*i],&eg->sidemap[3*i+1],&eg->sidemap[3*i+2]);
      }
      printf("Found %d boundary type mappings\n",i);
      eg->sidemappings = i;
    }
    else if(strstr(command,"BULK TYPE MAPPINGS")) {
      for(i=0;i<MAXMAPPINGS;i++) {
	if(i>0) Getline(params,in);
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	cp = params;      
	sscanf(params,"%d%d%d",&eg->bulkmap[3*i],&eg->bulkmap[3*i+1],&eg->bulkmap[3*i+2]);
      }
      printf("Found %d bulk type mappings\n",i);
      eg->bulkmappings = i;
    }
    else if(strstr(command,"COORDINATE MAPPING")) {
      sscanf(params,"%d%d%d",&eg->coordinatemap[0],&eg->coordinatemap[1],&eg->coordinatemap[2]);
    }	
    else if(strstr(command,"BOUNDARY BOUNDARY")) {
      for(i=0;i<MAXBOUNDARIES;i++) {
	if(i>0) Getline(params,in);
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	cp = params;      
	sscanf(params,"%d%d%d",&eg->boundbound[3*i+2],&eg->boundbound[3*i],&eg->boundbound[3*i+1]);
      }
      printf("Found %d boundary boundary definitions\n",i);
      eg->boundbounds = i;
    }
    else if(strstr(command,"MATERIAL BOUNDARY")) {
      for(i=0;i<MAXBOUNDARIES;i++) {
	if(i>0) Getline(params,in);
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	cp = params;      
	sscanf(params,"%d%d%d",&eg->bulkbound[3*i+2],&eg->bulkbound[3*i],&eg->bulkbound[3*i+1]);
      }
      printf("Found %d material boundary definitions\n",i);
      eg->bulkbounds = i;
    }

    else if(strstr(command,"RENUMBER BOUNDARY")) {
      for(i=0;i<MAXBOUNDARIES;i++) {
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	cp = params;      
	sscanf(params,"%d%d%d",&eg->sidemap[3*i],&eg->sidemap[3*i+1],&eg->sidemap[3*i+2]);
      }
      printf("Found %d boundary mappings\n",i);
      eg->sidemappings = i;
    }
    else if(strstr(command,"RENUMBER MATERIAL")) {
      for(i=0;i<MAXBOUNDARIES;i++) {
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	cp = params;      
	sscanf(params,"%d%d%d",&eg->bulkmap[3*i],&eg->bulkmap[3*i+1],&eg->bulkmap[3*i+2]);
      }
      printf("Found %d material mappings\n",i);
      eg->bulkmappings = i;
    }

    else if(strstr(command,"BOUNDARY LAYER")) {
      if(strstr(command,"BOUNDARY LAYER MOVE")) {
	sscanf(params,"%d",&eg->layermove);
      }
      else if(strstr(command,"BOUNDARY LAYER EPSILON")) {
	sscanf(params,"%le",&eg->layereps);
      }
      else {
	for(i=0;i<MAXBOUNDARIES;i++) {
	  if(i>0) Getline(params,in);
	  for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	  cp = params;      

	  if(strstr(params,"END") || strstr(params,"End") ) break;
	  eg->layerbounds[i] = next_int(&cp);
	  eg->layernumber[i] = next_int(&cp);
	  eg->layerthickness[i] = next_real(&cp);
	  eg->layerratios[i] = next_real(&cp);
	  eg->layerparents[i] = next_int(&cp);	  
	}
	printf("Found %d boundary layers\n",i);
	eg->layers = i;
      }
    }
    else if(strstr(command,"REMOVE LOWER DIMENSIONAL BOUNDARIES")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->removelowdim = TRUE; 
    }
    else if(strstr(command,"REMOVE INTERNAL BOUNDARIES")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->removeintbcs = TRUE; 
    }
    else if(strstr(command,"REMOVE UNUSED NODES")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->removeunused = TRUE; 
    }
    else if(strstr(command,"NO MESH NAMES")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->usenames = FALSE; 
    }
    else if(strstr(command,"REORDER MATERIAL")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->bulkorder = TRUE; 
    }
    else if(strstr(command,"REORDER BOUNDARY")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->boundorder = TRUE; 
    }
    else if(strstr(command,"DIMENSION")) {
      sscanf(params,"%d",&eg->dim);
    }
    else if(strstr(command,"ISOPARAMETRIC")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->isoparam = TRUE;
    }
    else if(strstr(command,"NO BOUNDARY")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->saveboundaries = FALSE;
    }
    else if(strstr(command,"LAYERED BOUNDARIES")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"TRUE")) grid->layeredbc = 1;
      if(strstr(params,"FALSE")) grid->layeredbc = 0;
    }
    else if(strstr(command,"EXTRUDED")) {
      grid->dimension = 3;

      if(strstr(command,"EXTRUDED DIVISIONS")) {
	sscanf(params,"%d",&grid->zcells);		
      }
      if(strstr(command,"EXTRUDED BC OFFSET")) {
	sscanf(params,"%d",&grid->layerbcoffset);		
      }
      else if(strstr(command,"EXTRUDED LIMITS")) {
	cp = params;
	for(i=0;i<=grid->zcells;i++ ) { 
	  grid->z[i] = next_real(&cp);
	  if(i > 0 && grid->z[i] < grid->z[i-1]) {
	    printf("Extruded limits %d: %12.6le %12.6le\n",i,grid->z[i],grid->z[i-1]);
	    bigerror("Values for limits should be a growing series, existing\n");
	  }
	}
      }
      else if(strstr(command,"EXTRUDED SIZES")) {
	cp = params;
	for(i=1;i<=grid->zcells;i++) grid->z[i] = next_real(&cp);
	for(i=1;i<=grid->zcells;i++) grid->z[i] = grid->z[i-1] + grid->z[i];
      }     
      else if(strstr(command,"EXTRUDED ELEMENTS")) {
	cp = params;
	for(i=1;i<=grid->zcells;i++) grid->zelems[i] = next_int(&cp);
	grid->autoratio = FALSE;    
      } 
      else if(strstr(command,"EXTRUDED RATIOS")) {
	cp = params;
	for(i=1;i<=grid->zcells;i++) grid->zexpand[i] = next_real(&cp);
      }
      else if(strstr(command,"EXTRUDED DENSITIES")) {
	cp = params;
	for(i=1;i<=grid->zcells;i++) grid->zdens[i] = next_real(&cp);
      }
      else if(strstr(command,"EXTRUDED STRUCTURE")) {
	for(i=1;i<= grid->zcells;i++) {
	  if(i>1) Getline(params,in);
	  sscanf(params,"%d %d %d\n",
		 &grid->zfirstmaterial[i],&grid->zlastmaterial[i],&grid->zmaterial[i]); 
	}
      }
      else if(strstr(command,"EXTRUDED MAX MATERIAL")) {
	sscanf(params,"%d",&grid->maxmaterial);		
      }
      else if(strstr(command,"EXTRUDED MATERIAL MAPPINGS")) {
	grid->zmaterialmap = Imatrix(1,grid->zcells,1,grid->maxmaterial);
	for(i=1;i<=grid->zcells;i++) {
	  if(i>1) Getline(params,in); 
	  cp = params;
	  for(j=1;j<=grid->maxmaterial;j++)
	    grid->zmaterialmap[i][j] = next_int(&cp);
	}
	grid->zmaterialmapexists = TRUE;
      }
      else if(strstr(command,"EXTRUDED HELICITY")) {
	sscanf(params,"%le",&grid->zhelicity);			
	grid->zhelicityexists = TRUE;
      }

    }
  }

end:
  printf("Read commands from a file\n");

  return(0);
}



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

    fprintf(out,"%-12.4e %-12.4e %-12.4e %-12.4e ",
	    data->x[sideind[0]],data->x[sideind[1]],
	    data->y[sideind[0]],data->y[sideind[1]]);
    for(k=0;k<MAXVARS;k++) 
      if(bound->evars[k]) {
	if(bound->points[k] == 1) 
	  fprintf(out,"%-10.4e ",bound->vars[k][i]);
      }		
    fprintf(out,"\n");
  }

  fclose(out);

  if(info) printf("Boundary information was saved to file %s.\n",filename);

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
	fprintf(out,"%-10.4e %-10.4e %-6d ",
		data->x[ind],data->y[ind],ind);
	for(k=0;k<MAXVARS;k++) 
	  if(bound[j].evars[k]) {
	    if(bound[j].points[k] == 0)
	      fprintf(out,"%-10.4e ",bound[j].vars[k][i]);
	    else if(bound[j].points[k] == 1) {
	      if(i==0)
		fprintf(out,"%-10.4e ",bound[j].vars[k][1]);		
	      else if(i==length)
		fprintf(out,"%-10.4e ",bound[j].vars[k][length]);		
	      else
		fprintf(out,"%-10.4e ",
			0.5*(bound[j].vars[k][i]+bound[j].vars[k][i+1]));	
	    }
	  }

	for(k=0;k<MAXDOFS;k++) {
	  if(data->edofs[k] == 1) 
	    fprintf(out,"%-10.4e  ",data->dofs[k][ind]);
	  if(data->edofs[k] == 2) 
	    fprintf(out,"%-10.4e  %-10.4e  ",
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
    printf("SaveBoundaryForm: stucture FemType not created\n");
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
	  
	  fprintf(out,"%-12.4e %-12.4e %-12.4e %-12.4e\n",
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
    printf("SaveBoundaryLine: stucture FemType not created\n");
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
	fprintf(out,"%-12.7e %-12.7e ",c,data->y[no]);
      else 
	fprintf(out,"%-12.7e %-12.7e ",data->x[no],c);	
      for(k=0;k<MAXDOFS;k++) {
	if(data->edofs[k] == 1) 
	  fprintf(out,"%-12.7e ",data->dofs[k][no]);
	if(data->edofs[k] == 2) 
	  fprintf(out,"%-12.7e %-12.7e ",
		  data->dofs[k][2*no-1],data->dofs[k][2*no]);
      }	
      fprintf(out,"\n");
      points++;
    }
  }
  if(info) printf("Line (c=%.3g) with %d nodes was saved to file %s.\n",
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
    printf("SaveSubcellForm: stucture FemType not created\n");
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
	    fprintf(out,"%-12.4e %-12.4e %-8d\n",
		    data->x[sideind[i]],data->y[sideind[i]],sideind[i]);
        } while(more);
      }
  fclose(out);

  if(info) printf("There are %d sideknots in the elements.\n",sideknots);
  if(info) printf("The positions of the sideknots were saved in file %s.\n",filename);

  return(0);
}



int SaveElmergrid(struct GridType *grid,int nogrids,char *prefix,int info)
{
  int sameline,maxsameline;
  int i,j,dim;
  FILE *out;
  char filename[MAXFILESIZE];

  AddExtension(prefix,filename,"grd");
  out = fopen(filename,"w");
  dim = grid->dimension;
  if(grid->coordsystem == COORD_CART1) dim = 1;

  j = 0;
  sameline = TRUE;
  maxsameline = 6;
  if(grid->xcells > maxsameline) sameline = FALSE;
  if(dim >= 2 && grid->ycells > maxsameline) sameline = FALSE;
  if(dim >= 3 && grid->zcells > maxsameline) sameline = FALSE;
  
  fprintf(out,"#####  ElmerGrid input file for structured grid generation  ######\n");
  fprintf(out,"Version = 210903\n");

  fprintf(out,"Coordinate System = ");
  if(grid->coordsystem == COORD_AXIS)
    fprintf(out,"2D Axisymmetric\n");
  else if(grid->coordsystem == COORD_POLAR)
    fprintf(out,"2D Polar\n");
  else 
    fprintf(out,"Cartesian %dD\n",dim);
 
  fprintf(out,"Subcell Divisions in %dD = ",dim);
  if(dim >= 1) fprintf(out,"%d ",grid->xcells);
  if(dim >= 2) fprintf(out,"%d ",grid->ycells);
  if(dim >= 3) fprintf(out,"%d ",grid->zcells);
  fprintf(out,"\n");

  fprintf(out,"Subcell Limits 1 %s",sameline ? "= ":"\n  ");
  for(i=0;i <= grid->xcells;i++) 
    fprintf(out,"%.5g ",grid->x[i]); 
  fprintf(out,"\n");
    
  if(dim >= 2) {
    fprintf(out,"Subcell Limits 2 %s",sameline ? "= ":"\n  ");
    for(i=0;i <= grid->ycells;i++) 
      fprintf(out,"%.5g ",grid->y[i]); 
    fprintf(out,"\n");
  }
  
  if(dim >= 3) {
    fprintf(out,"Subcell Limits 3 %s",sameline ? "= ":"\n  ");
    for(i=0;i <= grid->zcells;i++) 
      fprintf(out,"%.5g ",grid->z[i]); 
    fprintf(out,"\n");
  }  

  fprintf(out,"Material Structure in %dD\n",dim==1 ? 1:2);  
  for(j=grid->ycells;j>=1;j--) {
    fprintf(out,"  ");
    for(i=1;i<=grid->xcells;i++) 
      fprintf(out,"%-5d",grid->structure[j][i]);
    fprintf(out,"\n");
  }
  fprintf(out,"End\n");

  if(grid->mappings > 0) {
    fprintf(out,"Geometry Mappings\n");
    fprintf(out,"# mode  line  limits(2)   Np  params(Np)\n");
    for(i=0;i<grid->mappings;i++) {
      fprintf(out,"  %-5d %-5d %-7.5g %-7.5g %-3d ",
	      grid->mappingtype[i],grid->mappingline[i],
	      grid->mappinglimits[2*i],grid->mappinglimits[2*i+1],
	      grid->mappingpoints[i]);
      for(j=0;j<grid->mappingpoints[i];j++) 
	fprintf(out,"%.4g ",grid->mappingparams[i][j]);
      fprintf(out,"\n");
    }
    fprintf(out,"End\n");
  }

  j = 0;
  if(grid[j].rotate) {
    fprintf(out,"Revolve Blocks = %d\n",grid[j].rotateblocks);
    fprintf(out,"Revolve Radius = %-8.3g\n",grid[j].rotateradius2);
    if(fabs(grid[j].rotateimprove-1.0) > 1.0e-10)
      fprintf(out,"Revolve Improve = %-8.3g\n",grid[j].rotateimprove);
    
  }
  if(grid[j].rotatecurve) {
    fprintf(out,"Revolve Curve Direct = %-8.3g\n",grid[j].curvezet);
    fprintf(out,"Revolve Curve Radius = %-8.3g\n",grid[j].curverad);
    fprintf(out,"Revolve Curve Angle = %-8.3g\n",grid[j].curveangle);
  }

  if(grid[j].coordsystem == COORD_POLAR) {
    fprintf(out,"Polar Radius = %.3g\n",grid[j].polarradius);
  } 

  for(j=0;j<nogrids;j++) {
    
    if(j>0) fprintf(out,"\nStart New Mesh\n");
  
    fprintf(out,"Materials Interval = %d %d\n",
	    grid[j].firstmaterial,grid[j].lastmaterial);
  
    if(dim == 3) {
      fprintf(out,"Extruded Structure\n");
      fprintf(out,"# %-8s %-8s %-8s\n","1stmat", "lastmat","newmat");
      for(i=1;i<=grid[j].zcells;i++) 
	fprintf(out,"  %-8d %-8d %-8d\n",
		grid[j].zfirstmaterial[i],grid[j].zlastmaterial[i],
		grid[j].zmaterial[i]); 
      fprintf(out,"End\n");    
    }

    if(grid[j].noboundaries > 0) {
      fprintf(out,"Boundary Definitions\n");
      fprintf(out,"# %-8s %-8s %-8s\n","type","out","int"); 
      for(i=0;i<grid[j].noboundaries;i++)
	fprintf(out,"  %-8d %-8d %-8d %-8d\n",
		grid[j].boundtype[i],grid[j].boundext[i],
		grid[j].boundint[i], grid[j].boundsolid[i]);
      fprintf(out,"End\n");
    }

    if(grid->numbering == NUMBER_XY)
      fprintf(out,"Numbering = Horizontal\n");
    if(grid->numbering == NUMBER_YX)
      fprintf(out,"Numbering = Vertical\n");
        
    fprintf(out,"Element Degree = %d\n",grid[j].elemorder);
    fprintf(out,"Element Innernodes = %s\n",grid[j].elemmidpoints ? "True" : "False");
    fprintf(out,"Triangles = %s\n",grid[j].triangles ? "True" : "False");
    if(grid[j].autoratio) 
      fprintf(out,"Surface Elements = %d\n",grid[j].wantedelems);
    if(dim == 3 && grid[j].wantedelems3d) 
      fprintf(out,"Volume Elements = %d\n",grid[j].wantedelems3d);
    if(dim == 3 && grid[j].wantednodes3d) 
      fprintf(out,"Volume Nodes = %d\n",grid[j].wantednodes3d);

    if(dim == 2)
      fprintf(out,"Coordinate Ratios = %-8.3g\n",grid[j].xyratio);
    if(dim == 3)
      fprintf(out,"Coordinate Ratios = %-8.3g %-8.3g\n",
	      grid[j].xyratio,grid[j].xzratio);
 
    fprintf(out,"Minimum Element Divisions = %d",grid[j].minxelems);
    if(dim >= 2) fprintf(out," %d",grid[j].minyelems);
    if(dim >= 3) fprintf(out," %d",grid[j].minzelems);
    fprintf(out,"\n");

    fprintf(out,"Element Ratios 1 %s",sameline ? "= ":"\n  ");
    for(i=1;i<=grid[j].xcells;i++) 
      fprintf(out,"%.3g ",grid[j].xexpand[i]); 
    fprintf(out,"\n");
    if(dim >= 2) {
      fprintf(out,"Element Ratios 2 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].ycells;i++) 
	fprintf(out,"%.3g ",grid[j].yexpand[i]); 
      fprintf(out,"\n");
    }
    if(dim >= 3) {
      fprintf(out,"Element Ratios 3 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].zcells;i++) 
	fprintf(out,"%.3g ",grid[j].zexpand[i]); 
      fprintf(out,"\n");
    }

    if(grid[j].autoratio) {
      fprintf(out,"Element Densities 1 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].xcells;i++) 
	fprintf(out,"%.3g ",grid[j].xdens[i]); 
      fprintf(out,"\n");
      if(dim >= 2) {
	fprintf(out,"Element Densities 2 %s",sameline ? "= ":"\n  ");
	for(i=1;i<=grid[j].ycells;i++) 
	  fprintf(out,"%.3g ",grid[j].ydens[i]); 
	fprintf(out,"\n");
      }
      if(dim >= 3) {       
	fprintf(out,"Element Densities 3 %s",sameline ? "= ":"\n  ");
	for(i=1;i<=grid[j].zcells;i++) 
	  fprintf(out,"%.3g ",grid[j].zdens[i]); 
	fprintf(out,"\n");
      }
    }
    else {
      fprintf(out,"Element Divisions 1 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].xcells;i++) 
	fprintf(out,"%d ",grid[j].xelems[i]); 
      fprintf(out,"\n");
      if(dim >= 2) {
	fprintf(out,"Element Divisions 2 %s",sameline ? "= ":"\n  ");
	for(i=1;i<=grid[j].ycells;i++) 
	  fprintf(out,"%d ",grid[j].yelems[i]); 
	fprintf(out,"\n");
      }
      if(dim >= 3) {       
	fprintf(out,"Element Divisions 3 %s",sameline ? "= ":"\n  ");
	for(i=1;i<=grid[j].zcells;i++) 
	  fprintf(out,"%d ",grid[j].zelems[i]); 
	fprintf(out,"\n");
      }
    }
    

  }

  if(info) printf("The Elmergrid input was saved to file %s.\n",filename);
  fclose(out);

  return(0);
}




int LoadElmergridOld(struct GridType **grid,int *nogrids,char *prefix,int info) 
{
  char filename[MAXFILESIZE];
  FILE *in;
  int i,j,k,l,error=0;
  Real scaling;
  char *cp;
  int mode,noknots,noelements,dim,axisymmetric;
  int elemcode,maxnodes,totelems,nogrids0;
  int minmat,maxmat;
  char line[MAXLINESIZE];


  AddExtension(prefix,filename,"grd");
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadElmergrid: opening of the file '%s' wasn't succesfull !\n",filename);
    return(1);
  }

  if(info) printf("Loading the geometry from file '%s'.\n",filename);

  InitGrid(grid[*nogrids]);
  k = *nogrids;
  nogrids0 = *nogrids;

  mode = 0;
  noknots = 0;
  noelements = 0;
  dim = 0;
  axisymmetric = FALSE;
  elemcode = 0;
  maxnodes = 4;
  totelems = 0;
  scaling = 1.0;



  Getline(line,in);
  for(;;) {
    if(Getline(line,in)) goto end;
    if(line[0]=='\0') goto end;

    if(strstr(line,"END")) goto end;
    if(strstr(line,"RESULTS")) goto end;

    /* Control information */
    if(strstr(line,"VERSION")) mode = 1;
    else if(strstr(line,"GEOMETRY")) mode = 2;
    else if(strstr(line,"MAPPINGS IN")) mode = 31;
    else if(strstr(line,"MAPPINGS OUT")) mode = 32;
    else if(strstr(line,"MAPPINGS")) mode = 3;
    else if(strstr(line,"NUMBERING")) mode = 4;
    else if(strstr(line,"MESHING")) mode = 5;
    else if(strstr(line,"ELEMENTS")) mode = 6;
    else if(strstr(line,"ELEMENT NUMBER")) mode = 29;
    else if(strstr(line,"NODES")) mode = 7;
    else if(strstr(line,"TRIANGLE")) mode = 8;
    else if(strstr(line,"SQUARE")) mode = 17;
    else if(strstr(line,"COORDINATE RATIO"))  mode = 10;
    else if(strstr(line,"MATERIALS")) mode = 11;
    else if(strstr(line,"LAYERED ST")) mode = 12;
    else if(strstr(line,"ELEMENT RAT")) mode = 13;
    else if(strstr(line,"ELEMENT DENS")) mode = 14;
    else if(strstr(line,"ELEMENT MINIMUM")) mode = 27;
    else if(strstr(line,"BOUNDARY COND")) mode = 15;
    else if(strstr(line,"ELEMENTTYPE") || strstr(line,"ELEMENTCODE")) mode = 16;
    else if(strstr(line,"ROTATE")) mode = 20;
    else if(strstr(line,"ROTRAD")) mode = 21;
    else if(strstr(line,"ROTBLOCK")) mode = 22;
    else if(strstr(line,"ROTIMP")) mode = 24;
    else if(strstr(line,"ROTCURVE")) mode = 25;
    else if(strstr(line,"REDUCE ELEMENT")) mode = 26;
    else if(strstr(line,"SCALING")) mode = 23;
    else if(strstr(line,"LAYERED BO")) mode = 28;
    else if(strstr(line,"POLAR RADIUS")) mode = 30;


    switch (mode) {
    case 1: 
      printf("Loading Elmergrid file: %s\n",line);
      mode = 0;
      break;
      
    case 2:
      grid[k]->dimension = 2;
      if(strstr(line,"CARTES") && strstr(line,"1D")) {
	grid[k]->coordsystem = COORD_CART1;
	grid[k]->dimension = 1;
      }
      else if(strstr(line,"CARTES") && strstr(line,"2D")) 
	grid[k]->coordsystem = COORD_CART2;
      else if(strstr(line,"AXIS") && strstr(line,"2D")) 
	grid[k]->coordsystem = COORD_AXIS;
      else if(strstr(line,"POLAR") && strstr(line,"2D")) 
	grid[k]->coordsystem = COORD_POLAR;
      else if(strstr(line,"CARTES") && strstr(line,"3D")) {
	grid[k]->coordsystem = COORD_CART3;
	grid[k]->dimension = 3;
      }
      else if(strstr(line,"CYLINDRICAL")) {
	grid[k]->coordsystem = COORD_CYL;
	grid[k]->dimension = 3;
      }
      else printf("Unknown coordinate system: %s\n",line);
      printf("Defining the coordinate system (%d-DIM).\n",grid[k]->dimension);

      Getline(line,in);

      if(grid[k]->dimension == 1) {
	sscanf(line,"%d",&(*grid)[k].xcells);
	grid[k]->ycells = 1;	
      }
      if(grid[k]->dimension == 2) 
	sscanf(line,"%d %d",&(*grid)[k].xcells,&(*grid)[k].ycells);
      if(grid[k]->dimension == 3) 
	sscanf(line,"%d %d %d",&(*grid)[k].xcells,&(*grid)[k].ycells,&(*grid)[k].zcells);      
      if(grid[k]->xcells >= MAXCELLS || grid[k]->ycells >= MAXCELLS || 
	 grid[k]->zcells >= MAXCELLS) {
	printf("LoadGrid: Too many subcells [%d %d %d] vs. %d:\n",
	       grid[k]->xcells,grid[k]->ycells,grid[k]->zcells,MAXCELLS);
      }
      
      if(grid[k]->dimension == 1) {
	printf("Loading [%d] subcell intervals in 1D\n",
	       grid[k]->xcells);
      }
      else if(grid[k]->dimension == 2) {
	printf("Loading [%d %d] subcell intervals in 2D\n",
	       grid[k]->xcells,grid[k]->ycells);   
      } else {
	printf("Loading [%d %d %d] subcell intervals in 3D\n",
	       grid[k]->xcells,grid[k]->ycells,grid[k]->zcells);   
      }


      for(j=1;j<=grid[k]->dimension;j++) {
	Getline(line,in);
	cp=line;

	if(j==1) for(i=0;i<=grid[k]->xcells;i++) grid[k]->x[i] = next_real(&cp);
	if(j==2) for(i=0;i<=grid[k]->ycells;i++) grid[k]->y[i] = next_real(&cp);
	if(j==3) for(i=0;i<=grid[k]->zcells;i++) grid[k]->z[i] = next_real(&cp);
      }

      printf("Loading material structure\n");

      for(j=grid[k]->ycells;j>=1;j--) {
	
	Getline(line,in);
	cp=line;
	
	for(i=1;i<=grid[k]->xcells;i++) 
	  grid[k]->structure[j][i] = next_int(&cp);
      }

      minmat = maxmat = grid[k]->structure[1][1];
      for(j=grid[k]->ycells;j>=1;j--) 
	for(i=1;i<=grid[k]->xcells;i++) {
	  if(minmat > grid[k]->structure[j][i])
	    minmat = grid[k]->structure[j][i];
	  if(maxmat < grid[k]->structure[j][i])
	    maxmat = grid[k]->structure[j][i];
	}      
      if(minmat < 0) 
	printf("LoadElmergrid: please use positive material indices.\n");

      mode = 0;
      break;

    case 3:
    case 31:
    case 32:

      /* I dont know how to set this, luckily this piece of code should be obsolite */
      l = 1;
      for(i=grid[k]->mappings;i<grid[k]->mappings+l;i++) {
	Getline(line,in);
	cp=line; 

	grid[k]->mappingtype[i] = next_int(&cp);
	if(mode == 32) grid[k]->mappingtype[i] += 50*SGN(grid[k]->mappingtype[i]);

	grid[k]->mappingline[i] = next_int(&cp);
	grid[k]->mappinglimits[2*i] = next_real(&cp);
	grid[k]->mappinglimits[2*i+1] = next_real(&cp);
	grid[k]->mappingpoints[i] = next_int(&cp);
	grid[k]->mappingparams[i] = Rvector(0,grid[k]->mappingpoints[i]);
	for(j=0;j<grid[k]->mappingpoints[i];j++) 
	  grid[k]->mappingparams[i][j] = next_real(&cp);
      }
      
      printf("Loaded %d geometry mappings\n",l);
      grid[k]->mappings += l;

      mode = 0;
      break;
      
    case 4: /* NUMBERING */
      if(strstr(line,"HORIZ")) grid[k]->numbering = NUMBER_XY;
      if(strstr(line,"VERTI")) grid[k]->numbering = NUMBER_YX;
      mode = 0;
      break;

    case 5: /* MESHING */
      if((*nogrids) >= MAXCASES) {
	printf("There are more grids than was allocated for!\n"); 
	printf("Ignoring meshes starting from %d\n.",(*nogrids)+1);
	goto end;
      }
      (*nogrids)++;
      printf("Loading element meshing no %d\n",*nogrids);
      k = *nogrids - 1;	           
      if(k > nogrids0) (*grid)[k] = (*grid)[k-1];	 
      mode = 0;
      break;

    case 6: /* ELEMENTS */
      sscanf(line,"%d",&(*grid)[k].wantedelems);
      mode = 0;
      break;

    case 7: /* NODES */
      sscanf(line,"%d",&(*grid)[k].nonodes);      
      
      (*grid)[k].elemmidpoints = FALSE;
      if((*grid)[k].nonodes == 4) 
	(*grid)[k].elemorder = 1;
      if((*grid)[k].nonodes == 8) 
	(*grid)[k].elemorder = 2;
      if((*grid)[k].nonodes == 16) 
	(*grid)[k].elemorder = 3;

      if((*grid)[k].nonodes == 9) { 
	(*grid)[k].elemorder = 2;
	(*grid)[k].elemmidpoints = TRUE;
      }
      if((*grid)[k].nonodes == 12) { 
	(*grid)[k].elemorder = 3;
	(*grid)[k].elemmidpoints = TRUE;
      }


      mode = 0;
      break;

    case 8: /* TRIANGLES */
      (*grid)[k].triangles = TRUE;
      mode = 0;
      break;

    case 17: /* SQUARES */
      (*grid)[k].triangles = FALSE;
      mode = 0;
      break;

    case 16: /* ELEMENTTYPE and ELEMENTCODE */
      sscanf(line,"%d",&elemcode);
      if(elemcode/100 == 2) {
	(*grid)[k].triangles = FALSE;      
	(*grid)[k].nonodes = elemcode%100;
      }
      else if(elemcode/100 == 4) {
	(*grid)[k].triangles = FALSE;      
	(*grid)[k].nonodes = elemcode%100;
      }
      else if(elemcode/100 == 3) {  
	(*grid)[k].triangles = TRUE;      
	if(elemcode%100 == 3)       (*grid)[k].nonodes = 4;
	else if(elemcode%100 == 6)  (*grid)[k].nonodes = 9;
	else if(elemcode%100 == 10) (*grid)[k].nonodes = 16;	
      }

      (*grid)[k].elemmidpoints = FALSE;
      if((*grid)[k].nonodes == 4) 
	(*grid)[k].elemorder = 1;
      if((*grid)[k].nonodes == 8) 
	(*grid)[k].elemorder = 2;
      if((*grid)[k].nonodes == 16) 
	(*grid)[k].elemorder = 3;

      if((*grid)[k].nonodes == 9) { 
	(*grid)[k].elemorder = 2;
	(*grid)[k].elemmidpoints = TRUE;
      }
      if((*grid)[k].nonodes == 12) { 
	(*grid)[k].elemorder = 3;
	(*grid)[k].elemmidpoints = TRUE;
      }

      mode = 0;
      break;

    case 10: /* COORDINATE RATIO */
      if((*grid)[k].dimension == 2) 
	sscanf(line,"%le",&(*grid)[k].xyratio);
      if((*grid)[k].dimension == 3) 
	sscanf(line,"%le %le",&(*grid)[k].xyratio,&(*grid)[k].xzratio);      
      mode = 0;
      break;

    case 11: /* MATERIALS */
      sscanf(line,"%d %d",&(*grid)[k].firstmaterial,&(*grid)[k].lastmaterial);      
      mode = 0;
      break;

    case 12: /* LAYERES */
      for(i=1;i<=(*grid)[k].zcells;i++) {
	Getline(line,in);
	sscanf(line,"%d %d %d\n",
		&(*grid)[k].zfirstmaterial[i],&(*grid)[k].zlastmaterial[i],&(*grid)[k].zmaterial[i]); 
      }
      mode = 0;
      break;

    case 13: /* ELEMENT RATIOS */
      printf("Loading element ratios\n");

      for (j=1;j<=(*grid)[k].dimension;j++) {
	Getline(line,in);
	cp = line;

	if(j==1) for(i=1;i<=(*grid)[k].xcells;i++) (*grid)[k].xexpand[i] = next_real(&cp);
	if(j==2) for(i=1;i<=(*grid)[k].ycells;i++) (*grid)[k].yexpand[i] = next_real(&cp);
	if(j==3) for(i=1;i<=(*grid)[k].zcells;i++) (*grid)[k].zexpand[i] = next_real(&cp);
      }
      mode = 0;
      break;

    case 29: /* ELEMENT NUMBER */
      printf("Loading element numbers\n");

      for (j=1;j<=(*grid)[k].dimension;j++) {
	Getline(line,in);
	cp = line;
	if(j==1) for(i=1;i<=(*grid)[k].xcells;i++) (*grid)[k].xelems[i] = next_int(&cp);
	if(j==2) for(i=1;i<=(*grid)[k].ycells;i++) (*grid)[k].yelems[i] = next_int(&cp);
	if(j==3) for(i=1;i<=(*grid)[k].zcells;i++) (*grid)[k].zelems[i] = next_int(&cp);
      }
      (*grid)[k].autoratio = 0;
      mode = 0;
      break;

    case 27: /* ELEMENT MINIMUM */
      printf("Loading minimum number of elements\n");
      if((*grid)[k].dimension == 1) 
	sscanf(line,"%d",&(*grid)[k].minxelems);
      if((*grid)[k].dimension == 2) 
	sscanf(line,"%d %d",&(*grid)[k].minxelems,&(*grid)[k].minyelems);
      if((*grid)[k].dimension == 3) 
	sscanf(line,"%d %d %d",&(*grid)[k].minxelems,&(*grid)[k].minyelems,&(*grid)[k].minzelems);
      mode = 0;
      break;

    case 14: /* ELEMENT DENSITIES */
      printf("Loading element densities\n");
      for (j=1;j<=(*grid)[k].dimension;j++) {
	Getline(line,in);
	cp = line;

	if(j==1) for(i=1;i<=(*grid)[k].xcells;i++) (*grid)[k].xdens[i] = next_real(&cp);
	if(j==2) for(i=1;i<=(*grid)[k].ycells;i++) (*grid)[k].ydens[i] = next_real(&cp);
	if(j==3) for(i=1;i<=(*grid)[k].zcells;i++) (*grid)[k].zdens[i] = next_real(&cp);
      }
      mode = 0;
      break;

    case 15: /* BOUNDARY CONDITIONS */
      sscanf(line,"%d",&(*grid)[k].noboundaries);
      printf("Loading %d boundary conditions\n",(*grid)[k].noboundaries);

      for(i=0;i<(*grid)[k].noboundaries;i++) {
	Getline(line,in);
	sscanf(line,"%d %d %d %d",
	       &(*grid)[k].boundtype[i],&(*grid)[k].boundext[i],
	       &(*grid)[k].boundint[i],&(*grid)[k].boundsolid[i]);
      }  
      mode = 0;
      break;

    case 20: /* ROTATE */
      (*grid)[k].rotate = TRUE;
      mode = 0;
      break;

    case 21: /* ROTRAD */
      sscanf(line,"%le",&(*grid)[k].rotateradius2);
      mode = 0;
      break;

    case 22: /* ROTBLOCK */
      sscanf(line,"%d",&(*grid)[k].rotateblocks);
      if(0) printf("Reading blocks %d\n",(*grid)[k].rotateblocks);
      mode = 0;
      break;

    case 24: /* ROTIMP */
      sscanf(line,"%le",&(*grid)[k].rotateimprove);
      mode = 0;
      break;

    case 30: /* POLAR RADIUS */
      sscanf(line,"%le",&(*grid)[k].polarradius);
      mode = 0;
      break;

    case 25: /* ROTCURVE */
      (*grid)[k].rotatecurve = TRUE;
      sscanf(line,"%le%le%le",&(*grid)[k].curvezet,
	     &(*grid)[k].curverad,&(*grid)[k].curveangle);
      mode = 0;
      break;

    case 26: /* REDUCE ELEMENT */
      sscanf(line,"%d%d",&(*grid)[k].reduceordermatmin,
	     &(*grid)[k].reduceordermatmax);
      mode = 0;
      break;

    case 28: /* LAYERED BO */
      sscanf(line,"%d",&(*grid)[k].layeredbc);
      mode = 0;
      break;

    case 23: /* SCALING */
      sscanf(line,"%le",&scaling);
      for(i=0;i<=(*grid)[k].xcells;i++) (*grid)[k].x[i] *= scaling;
      if((*grid)[k].dimension > 1) 
	for(i=0;i<=(*grid)[k].ycells;i++) (*grid)[k].y[i] *= scaling;
      if((*grid)[k].dimension == 3) 
	for(i=0;i<=(*grid)[k].ycells;i++) (*grid)[k].z[i] *= scaling;

      (*grid)[k].rotateradius2 *= scaling;
      (*grid)[k].curverad *= scaling;
      (*grid)[k].curvezet *= scaling;
      mode = 0;
      break;

    default:
      if(0) printf("Unknown case: %s",line);
    }

  }

end:

  if(info) printf("Found %d divisions for grid\n",*nogrids);

  for(k=nogrids0;k < (*nogrids) && k<MAXCASES;k++) {
    SetElementDivision(&(*grid)[k],1.0,info);
  }


  fclose(in);
  return(error);
}



int LoadElmergrid(struct GridType **grid,int *nogrids,char *prefix,Real relh,int info) 
{
  char filename[MAXFILESIZE];
  char command[MAXLINESIZE],params[MAXLINESIZE];
  FILE *in;
  int i,j,k;
  char *cp;
  int noknots,noelements,dim,axisymmetric;
  int elemcode,maxnodes,totelems,nogrids0,minmat,maxmat;
  long code;
  Real raid;

  AddExtension(prefix,filename,"grd");
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadElmergrid: opening of the file '%s' wasn't succesfull !\n",filename);
    return(1);
  }

  if(info) printf("Loading the geometry from file '%s'.\n",filename);

  InitGrid(grid[*nogrids]);
  k = *nogrids;
  nogrids0 = *nogrids;

  iodebug = FALSE;

  noknots = 0;
  noelements = 0;
  dim = 0;
  axisymmetric = FALSE;
  elemcode = 0;
  maxnodes = 4;
  totelems = 0;

  matcactive = FALSE;

  for(;;) {
    if(GetCommand(command,params,in)) {
      printf("Reached the end of command file\n");
      goto end;
    }    

    /* Control information */
    if(strstr(command,"VERSION")) {
      if(strstr(command,"080500")) {
	printf("Loading old version of Elmergrid file.\n");
	i = LoadElmergridOld(grid,nogrids,prefix,info);
	return(i);
      }
      else {
	sscanf(params,"%ld",&code);
	if(code == 210903) 
	  printf("Loading ElmerGrid file version: %ld\n",code);
	else {
	  printf("Unknown ElmerGrid file version: %ld\n",code);
	  return(2);
	}
      }
      *nogrids += 1;
    }      

    else if(strstr(command,"DEBUG IO")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"FALSE")) 
	iodebug = FALSE;
      else {
	iodebug = TRUE;
	printf("IO debugging activated\n");
      }
    }
 
    else if(strstr(command,"MATC")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"FALSE")) 
	matcactive = FALSE;
      else {
#ifndef DISABLE_MATC
	matcactive = TRUE;
	mtc_init(NULL, stdout, stderr);
	strcpy(command, "format( 12 )");	
	mtc_domath(command);	 
	printf("MATC language activated with 12 digit accuracy.\n");	
#else
        matcactive = FALSE;
        printf("Unable to activate matc as it is not even compiled.\n");
#endif 
      }
    }
    
    else if(strstr(command,"COORDINATE SYSTEM")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      grid[k]->dimension = 2;
      if(strstr(params,"CARTESIAN 1D")) {
	grid[k]->coordsystem = COORD_CART1;
	grid[k]->dimension = 1;
      }
      else if(strstr(params,"CARTESIAN 2D")) 
	grid[k]->coordsystem = COORD_CART2;
      else if(strstr(params,"AXISYMMETRIC")) 
	grid[k]->coordsystem = COORD_AXIS;
      else if(strstr(params,"POLAR"))
	grid[k]->coordsystem = COORD_POLAR;
      else if(strstr(params,"CARTESIAN 3D")) {
	grid[k]->coordsystem = COORD_CART3;
	grid[k]->dimension = 3;
      }
      else printf("Unknown coordinate system: %s\n",params);
      printf("Defining the coordinate system (%d-DIM).\n",grid[k]->dimension);
    }
    
    else if(strstr(command,"SUBCELL DIVISIONS")) {
      if(grid[k]->dimension == 1) {
	sscanf(params,"%d",&(*grid)[k].xcells);
	grid[k]->ycells = 1;	
      }
      else if(grid[k]->dimension == 2) 
	sscanf(params,"%d %d",&(*grid)[k].xcells,&(*grid)[k].ycells);
      else if(grid[k]->dimension == 3) 
	sscanf(params,"%d %d %d",&(*grid)[k].xcells,&(*grid)[k].ycells,&(*grid)[k].zcells);      

      if(grid[k]->xcells >= MAXCELLS || grid[k]->ycells >= MAXCELLS || grid[k]->zcells >= MAXCELLS) {
	printf("LoadElmergrid: Too many subcells [%d %d %d] vs. %d:\n",
	       grid[k]->xcells,grid[k]->ycells,grid[k]->zcells,MAXCELLS);
     }

      /* Initialize the default stucture with ones */
      for(j=grid[k]->ycells;j>=1;j--) 
	for(i=1;i<=grid[k]->xcells;i++) 
	  grid[k]->structure[j][i] = 1;
    }
    
    else if(strstr(command,"MINIMUM ELEMENT DIVISION")) {
      printf("Loading minimum number of elements\n");

      if((*grid)[k].dimension == 1) 
	sscanf(params,"%d",&(*grid)[k].minxelems); 

      if((*grid)[k].dimension == 2)  
	sscanf(params,"%d %d",&(*grid)[k].minxelems,&(*grid)[k].minyelems);

      if((*grid)[k].dimension == 3) 
	sscanf(params,"%d %d %d",&(*grid)[k].minxelems,
	       &(*grid)[k].minyelems,&(*grid)[k].minzelems);
    }      
    
    else if(strstr(command,"SUBCELL LIMITS 1")) {
      printf("Loading %d subcell limits in X-direction\n",grid[k]->xcells+1);
      cp = params;
      for(i=0;i<=grid[k]->xcells;i++) {
	grid[k]->x[i] = next_real(&cp);
	if(i > 0 && grid[k]->x[i] < grid[k]->x[i-1]) {
	  printf("Subcell limits 1(%d): %12.6le %12.6le\n",i,grid[k]->x[i],grid[k]->x[i-1]);
	  bigerror("Values for limits 1 should be a growing series, existing\n");
	}
      }
    }    
    else if(strstr(command,"SUBCELL LIMITS 2")) {
      printf("Loading %d subcell limits in Y-direction\n",grid[k]->ycells+1);
      cp = params;
      for(i=0;i<=grid[k]->ycells;i++) {
	grid[k]->y[i] = next_real(&cp);
	if(i > 0 && grid[k]->y[i] < grid[k]->y[i-1]) {
	  printf("Subcell limits 2(%d): %12.6le %12.6le\n",i,grid[k]->y[i],grid[k]->y[i-1]);
	  bigerror("Values for limits should be a growing series, existing\n");
	}
      }
    }      
    else if(strstr(command,"SUBCELL LIMITS 3")) {
      printf("Loading %d subcell limits in Z-direction\n",grid[k]->zcells+1);
      cp = params;
      for(i=0;i<=grid[k]->zcells;i++) {
	grid[k]->z[i] = next_real(&cp);
	if(i > 0 && grid[k]->z[i] < grid[k]->z[i-1]) {
	  printf("Subcell limits 3(%d): %12.6le %12.6le\n",i,grid[k]->z[i],grid[k]->z[i-1]);
	  bigerror("Values for limits should be a growing series, existing\n");
	}
      }
    }

    else if(strstr(command,"SUBCELL SIZES 1")) {
      printf("Loading %d subcell sizes in X-direction\n",grid[k]->xcells);
      cp = params;
      for(i=1;i<=grid[k]->xcells;i++) grid[k]->x[i] = next_real(&cp);
      for(i=1;i<=grid[k]->xcells;i++) grid[k]->x[i] = grid[k]->x[i-1] + grid[k]->x[i];
    }      
    else if(strstr(command,"SUBCELL SIZES 2")) {
      printf("Loading %d subcell sizes in Y-direction\n",grid[k]->ycells);
      cp = params;
      for(i=1;i<=grid[k]->ycells;i++) grid[k]->y[i] = next_real(&cp);
      for(i=1;i<=grid[k]->ycells;i++) grid[k]->y[i] = grid[k]->y[i-1] + grid[k]->y[i];
    }      
    else if(strstr(command,"SUBCELL SIZES 3")) {
      printf("Loading %d subcell sizes in Z-direction\n",grid[k]->zcells);
      cp = params;
      for(i=1;i<=grid[k]->zcells;i++) grid[k]->z[i] = next_real(&cp);
      for(i=1;i<=grid[k]->zcells;i++) grid[k]->z[i] = grid[k]->z[i-1] + grid[k]->z[i];
    }

    else if(strstr(command,"SUBCELL ORIGIN 1")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"CENTER")) {
	raid = 0.5 * (grid[k]->x[0] + grid[k]->x[grid[k]->xcells]);
      }
      else if(strstr(params,"LEFT") || strstr(params,"MIN") ) {
	raid = grid[k]->x[0];
      }
      else if(strstr(params,"RIGHT") || strstr(params,"MAX") ) {
	raid = grid[k]->x[grid[k]->xcells];
      }
      else {
	cp = params;
	raid = next_real(&cp);
      }
      for(i=0;i<=grid[k]->xcells;i++) grid[k]->x[i] -= raid;
    }
    else if(strstr(command,"SUBCELL ORIGIN 2")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"CENTER")) {
	raid = 0.5 * (grid[k]->y[0] + grid[k]->y[grid[k]->ycells]);
      }
      else if(strstr(params,"LEFT")) {
	raid = grid[k]->y[0];
      }
      else if(strstr(params,"RIGHT")) {
	raid = grid[k]->y[grid[k]->ycells];
      }
      else {
	cp = params;
	raid = next_real(&cp);
      }      
      for(i=0;i<=grid[k]->ycells;i++) grid[k]->y[i] -= raid;
    }
    else if(strstr(command,"SUBCELL ORIGIN 3")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"CENTER")) {
	raid = 0.5 * (grid[k]->z[0] + grid[k]->z[grid[k]->zcells]);
      }
      else if(strstr(params,"LEFT")) {
	raid = grid[k]->z[0];
      }
      else if(strstr(params,"RIGHT")) {
	raid = grid[k]->z[grid[k]->zcells];
      }
      else {
	cp = params;
	raid = next_real(&cp);
      }
      for(i=0;i<=grid[k]->zcells;i++) grid[k]->z[i] -= raid;      
    }

    else if(strstr(command,"MATERIAL STRUCTURE")) {
      printf("Loading material structure\n");

      /* Initialize the default stucture with zeros */
      for(j=grid[k]->ycells;j>=1;j--) 
	for(i=1;i<=grid[k]->xcells;i++) 
	  grid[k]->structure[j][i] = 0;
     
      for(j=grid[k]->ycells;j>=1;j--) {
	if(j < grid[k]->ycells) Getline(params,in);
	cp=params;
	for(i=1;i<=grid[k]->xcells;i++) 
	  grid[k]->structure[j][i] = next_int(&cp);
      }      
      minmat = maxmat = grid[k]->structure[1][1];
      for(j=grid[k]->ycells;j>=1;j--) 
	for(i=1;i<=grid[k]->xcells;i++) {
	  if(minmat > grid[k]->structure[j][i])
	    minmat = grid[k]->structure[j][i];
	  if(maxmat < grid[k]->structure[j][i])
	    maxmat = grid[k]->structure[j][i];
	}      
      if(minmat < 0) 
	printf("LoadElmergrid: please use positive material indices.\n");
      grid[k]->maxmaterial = maxmat;
    }
    else if(strstr(command,"MATERIALS INTERVAL")) {
      sscanf(params,"%d %d",&(*grid)[k].firstmaterial,&(*grid)[k].lastmaterial);      
    }
     
    else if(strstr(command,"REVOLVE")) {
      if(0) printf("revolve: %s %s\n",command,params);

      if(strstr(command,"REVOLVE RADIUS")) {
	(*grid)[k].rotate = TRUE;
	sscanf(params,"%le",&(*grid)[k].rotateradius2);
      }
      else if(strstr(command,"REVOLVE BLOCKS")) {
	(*grid)[k].rotate = TRUE;
	sscanf(params,"%d",&(*grid)[k].rotateblocks);
      }
      else if(strstr(command,"REVOLVE IMPROVE")) {
	(*grid)[k].rotate = TRUE;
	sscanf(params,"%le",&(*grid)[k].rotateimprove);
      }
      else if(strstr(command,"REVOLVE RADIUS")) {
	sscanf(params,"%le",&(*grid)[k].polarradius);
      }
      else if(strstr(command,"REVOLVE CURVE DIRECT")) {
	(*grid)[k].rotatecurve = TRUE;
	sscanf(params,"%le",&(*grid)[k].curvezet);
      }
      else if(strstr(command,"REVOLVE CURVE RADIUS")) {
	(*grid)[k].rotatecurve = TRUE;
	sscanf(params,"%le",&(*grid)[k].curverad);
      }
      else if(strstr(command,"REVOLVE CURVE ANGLE")) {
	(*grid)[k].rotatecurve = TRUE;
	sscanf(params,"%le",&(*grid)[k].curveangle);
      }
    }

    else if(strstr(command,"REDUCE ORDER INTERVAL")) {
      sscanf(params,"%d%d",&(*grid)[k].reduceordermatmin,
	     &(*grid)[k].reduceordermatmax);
    }
    
    else if(strstr(command,"BOUNDARY DEFINITION")) {
      printf("Loading boundary conditions\n");
      
      for(i=0;i<MAXBOUNDARIES;i++) {
	if(i>0) Getline(params,in);
	for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
	if(strstr(params,"END")) break;
	sscanf(params,"%d %d %d %d",
	       &(*grid)[k].boundtype[i],&(*grid)[k].boundext[i],
	       &(*grid)[k].boundint[i],&(*grid)[k].boundsolid[i]);
      }  
      printf("Found %d boundaries\n",i);
      (*grid)[k].noboundaries = i;
    }
    
    /* else if(strstr(command,"LAYERED BOUNDARIES")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"TRUE")) (*grid)[k].layeredbc = 1;
      if(strstr(params,"FALSE")) (*grid)[k].layeredbc = 0;
    } */
    
    else if(strstr(command,"NUMBERING")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"HORIZONATAL")) (*grid)[k].numbering = NUMBER_XY;
      if(strstr(params,"VERTICAL")) (*grid)[k].numbering = NUMBER_YX;
    }
    
    else if(strstr(command,"ELEMENT DEGREE")) {
      sscanf(params,"%d",&(*grid)[k].elemorder);
    }
    
    else if(strstr(command,"ELEMENT INNERNODES")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"TRUE")) (*grid)[k].elemmidpoints = TRUE;
      if(strstr(params,"FALSE")) (*grid)[k].elemmidpoints = FALSE;
    }
    else if(strstr(command,"ELEMENTTYPE") || strstr(command,"ELEMENTCODE")) {
      sscanf(params,"%d",&elemcode);
    }
    
    else if(strstr(command,"TRIANGLES CRITICAL ANGLE")) {
      sscanf(params,"%le",&(*grid)[k].triangleangle);      
    }
    else if(strstr(command,"TRIANGLES")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"TRUE")) (*grid)[k].triangles = TRUE;
      if(strstr(params,"FALSE")) (*grid)[k].triangles = FALSE;
    }
    
    else if(strstr(command,"PLANE ELEMENTS")) {
      sscanf(params,"%d",&(*grid)[k].wantedelems);
    }
    else if(strstr(command,"SURFACE ELEMENTS")) {
      sscanf(params,"%d",&(*grid)[k].wantedelems);
    }
    else if(strstr(command,"REFERENCE DENSITY")) {
      sscanf(params,"%le",&(*grid)[k].limitdx);
      (*grid)[k].autoratio = 3;     
    }
    else if(strstr(command,"VERIFY DENSITY")) {
      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      if(strstr(params,"TRUE")) (*grid)[k].limitdxverify = TRUE;
      if(strstr(params,"FALSE")) (*grid)[k].limitdxverify = FALSE;
    }
    else if(strstr(command,"VOLUME ELEMENTS")) {
      sscanf(params,"%d",&(*grid)[k].wantedelems3d);
    }
    else if(strstr(command,"VOLUME NODES")) {
      sscanf(params,"%d",&(*grid)[k].wantednodes3d);
    }

    else if(strstr(command,"COORDINATE RATIO")) {
      if((*grid)[k].dimension == 2) 
	sscanf(params,"%le",&(*grid)[k].xyratio);
      if((*grid)[k].dimension == 3) 
	sscanf(params,"%le %le",&(*grid)[k].xyratio,&(*grid)[k].xzratio);      
    }
    
    else if(strstr(command,"ELEMENT RATIOS 1")) {
      cp = params;
      for(i=1;i<=(*grid)[k].xcells;i++) (*grid)[k].xexpand[i] = next_real(&cp);
    }
    else if(strstr(command,"ELEMENT RATIOS 2")) {
      cp = params;
      for(i=1;i<=(*grid)[k].ycells;i++) (*grid)[k].yexpand[i] = next_real(&cp);
    }
    else if(strstr(command,"ELEMENT RATIOS 3")) {
      cp = params;
      for(i=1;i<=(*grid)[k].zcells;i++) (*grid)[k].zexpand[i] = next_real(&cp);
    }
    
    else if(strstr(command,"ELEMENT DENSITIES 1")) {
      cp = params;
      for(i=1;i<=(*grid)[k].xcells;i++) (*grid)[k].xdens[i] = next_real(&cp);
    }
    else if(strstr(command,"ELEMENT DENSITIES 2")) {
      cp = params;
      for(i=1;i<=(*grid)[k].ycells;i++) (*grid)[k].ydens[i] = next_real(&cp);
    }
    else if(strstr(command,"ELEMENT DENSITIES 3")) {
      cp = params;
      for(i=1;i<=(*grid)[k].zcells;i++) (*grid)[k].zdens[i] = next_real(&cp);
    }
    
    else if(strstr(command,"ELEMENT DIVISIONS 1")) {
      cp = params;
      for(i=1;i<=(*grid)[k].xcells;i++) (*grid)[k].xelems[i] = next_int(&cp);
      (*grid)[k].autoratio = 0;
    }
    else if(strstr(command,"ELEMENT DIVISIONS 2")) {
      cp = params;
      for(i=1;i<=(*grid)[k].ycells;i++) (*grid)[k].yelems[i] = next_int(&cp);
      (*grid)[k].autoratio = 0;
    }
    else if(strstr(command,"ELEMENT DIVISIONS 3")) {
      cp = params;
      for(i=1;i<=(*grid)[k].zcells;i++) (*grid)[k].zelems[i] = next_int(&cp);
      (*grid)[k].autoratio = 0;
    }
    
    /* else if(strstr(command,"EXTRUDED STRUCTURE")) {
      for(i=1;i<=(*grid)[k].zcells;i++) {
	if(i>1) Getline(params,in);
	sscanf(params,"%d %d %d\n",
	       &(*grid)[k].zfirstmaterial[i],&(*grid)[k].zlastmaterial[i],&(*grid)[k].zmaterial[i]); 
      }
    } */
    
    else if(strstr(command,"GEOMETRY MAPPINGS")) {     
      if(k > 0) (*grid)[k].mappings = 0;

      for(i=0;i<MAXLINESIZE;i++) params[i] = toupper(params[i]);
      for(i=(*grid)[k].mappings;i<MAXMAPPINGS;i++) {
	if(i>(*grid)[k].mappings) Getline(params,in);

	if(strstr(params,"END")) break;
	cp=params; 
	(*grid)[k].mappingtype[i] = next_int(&cp);	
#if 0
	(*grid)[k].mappingtype[i] += 50*SGN((*grid)[k].mappingtype[i]);
#endif
	(*grid)[k].mappingline[i] = next_int(&cp);
	(*grid)[k].mappinglimits[2*i] = next_real(&cp);
	(*grid)[k].mappinglimits[2*i+1] = next_real(&cp);
	(*grid)[k].mappingpoints[i] = next_int(&cp);
	(*grid)[k].mappingparams[i] = Rvector(0,(*grid)[k].mappingpoints[i]);
	for(j=0;j<(*grid)[k].mappingpoints[i];j++) 
	  (*grid)[k].mappingparams[i][j] = next_real(&cp);
      }      
      printf("Loaded %d geometry mappings\n",i);
      (*grid)[k].mappings = i;      
    }

    else if(strstr(command,"END") ) {      
      if(0) printf("End of field\n");
    }
      
    else if(strstr(command,"START NEW MESH")) {
      if((*nogrids) >= MAXCASES) {
	printf("There are more grids than was allocated for!\n"); 
	printf("Ignoring meshes starting from %d\n.",(*nogrids)+1);
	goto end;
      }
      (*nogrids)++;
      printf("\nLoading element meshing no %d\n",*nogrids);
      k = *nogrids - 1;	           
      if(k > nogrids0) (*grid)[k] = (*grid)[k-1];	 
    }

    else {
      if(0) printf("Unknown command: %s",command);
    }
  }

end:

  if(info) printf("Found %d divisions for grid\n",*nogrids);
  
  for(k=nogrids0;k < (*nogrids) && k<MAXCASES;k++) {

    if(elemcode == 0) {
      if((*grid)[k].dimension == 1) {
	(*grid)[k].nonodes = (*grid)[k].elemorder + 1;
      }
      else if((*grid)[k].elemmidpoints == FALSE) {
	(*grid)[k].nonodes = 4 * (*grid)[k].elemorder;
      }					
      else {
	if((*grid)[k].elemorder == 2) (*grid)[k].nonodes = 9;
	if((*grid)[k].elemorder == 3) (*grid)[k].nonodes = 16;	
      }
    }
    else if(elemcode/100 == 2) {
      (*grid)[k].triangles = FALSE;      
      (*grid)[k].nonodes = elemcode%100;
    }
    else if(elemcode/100 == 4) {
      (*grid)[k].triangles = FALSE;      
      (*grid)[k].nonodes = elemcode%100;
    }
    else if(elemcode/100 == 3) {  
      (*grid)[k].triangles = TRUE;      
      if(elemcode%100 == 3)       (*grid)[k].nonodes = 4;
      else if(elemcode%100 == 6)  (*grid)[k].nonodes = 9;
      else if(elemcode%100 == 10) (*grid)[k].nonodes = 16;	
    }    
  }

  
  for(k=nogrids0;k < (*nogrids) && k<MAXCASES;k++) {
    SetElementDivision(&(*grid)[k],relh,info);
  }

  fclose(in);
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





