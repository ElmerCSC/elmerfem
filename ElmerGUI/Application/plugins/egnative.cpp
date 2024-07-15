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

/* -------------------------------:  egnative.c  :----------------------------
   This module includes routines for I/O of native formats of Elmer. 
*/

#include <stdio.h>
#include <unistd.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <limits.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egmesh.h"
/* #include "egparallel.h" */
#include "egnative.h"
/*#include "../config.h"*/

#define GETLINE ioptr=fgets(line,MAXLINESIZE,in) 
static char *ioptr;


#define DEBUG 0


int matcactive=FALSE, iodebug=FALSE;

#define MAXINMETHODS 21
const char *InMethods[] = {
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
  /*16*/ "FVCOM",
  /*17*/ "NASTRAN",
  /*18*/ "CGSIM",
  /*19*/ "GEO",
  /*20*/ "FLUX2D",
  /*21*/ "FLUX3D",
};


#define MAXOUTMETHODS 5
const char *OutMethods[] = {
  /*0*/ "EG",
  /*1*/ "ELMERGRID",
  /*2*/ "ELMERSOLVER",
  /*3*/ "ELMERPOST",
  /*4*/ "GMSH",
  /*5*/ "VTU",
};


void Instructions()
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
  printf("1)  .grd      : ElmerGrid file format\n");
  printf("2)  .mesh.*   : Elmer input format\n");
  printf("3)  .ep       : Elmer output format\n");
  printf("4)  .ansys    : Ansys input format\n");
  printf("5)  .inp      : Abaqus input format by Ideas\n");
  printf("6)  .fil      : Abaqus output format\n");
  printf("7)  .FDNEUT   : Gambit (Fidap) neutral file\n");
  printf("8)  .unv      : Universal mesh file format\n");
  printf("9)  .mphtxt   : Comsol Multiphysics mesh format\n");
  printf("10) .dat      : Fieldview format\n");
  printf("11) .node,.ele: Triangle 2D mesh format\n");
  printf("12) .mesh     : Medit mesh format\n");
  printf("13) .msh      : GID mesh format\n");
  printf("14) .msh      : Gmsh mesh format\n");
  printf("15) .ep.i     : Partitioned ElmerPost format\n");
  printf("16) .2dm      : 2D triangular FVCOM format\n");
#if 0
  printf("17) .msh      : Nastran format\n");
  printf("18) .msh      : CGsim format\n");
  printf("19) .geo      : Geo format\n");
  printf("20) .tra      : Cedrat Flux 2D format\n");
  printf("21) .pf3      : Cedrat Flux 3D format\n");
#endif 

  printf("\nThe second parameter defines the output file format:\n");
  printf("1)  .grd      : ElmerGrid file format\n");
  printf("2)  .mesh.*   : ElmerSolver format (also partitioned .part format)\n");
  printf("3)  .ep       : ElmerPost format\n");
  printf("4)  .msh      : Gmsh mesh format\n");
  printf("5)  .vtu      : VTK ascii XML format\n");
#if 0
  printf("5)  .inp      : Abaqus input format\n");
  printf("7)  .fidap    : Fidap format\n");
  printf("18) .ep       : Fastcap input format.\n");
#endif

  printf("\nThe third parameter is the name of the input file.\n");
  printf("If the file does not exist, an example with the same name is created.\n");
  printf("The default output file name is the same with a different suffix.\n\n");

  printf("There are several additional in-line parameters that are\n");
  printf("taken into account only when applicable to the given format.\n");

  printf("-out str             : name of the output file\n");
  printf("-in str              : name of a secondary input file\n");
  printf("-decimals            : number of decimals in the saved mesh (eg. 8)\n");
  printf("-relh real           : give relative mesh density parameter for ElmerGrid meshing\n");
  printf("-triangles           : rectangles will be divided to triangles\n");
  printf("-merge real          : merges nodes that are close to each other\n");
  printf("-order real[3]       : reorder elements and nodes using c1*x+c2*y+c3*z\n");
  printf("-centralize          : set the center of the mesh to origin\n");
  printf("-scale real[3]       : scale the coordinates with vector real[3]\n");
  printf("-translate real[3]   : translate the nodes with vector real[3]\n");
  printf("-rotate real[3]      : rotate around the main axis with angles real[3]\n");
  printf("-clone int[3]        : make identical copies of the mesh\n");
  printf("-clonesize real[3]   : the size of the mesh to be cloned if larger to the original\n");
  printf("-mirror int[3]       : copy the mesh around the origin in coordinate directions\n");
  printf("-cloneinds           : when performing cloning should cloned entities be given new indexes\n");
  printf("-unite               : the meshes will be united\n");
  printf("-unitenooverlap      : the meshes will be united without overlap in entity numbering\n");
  printf("-polar real          : map 2D mesh to a cylindrical shell with given radius\n");
  printf("-cylinder            : map 2D/3D cylindrical mesh to a cartesian mesh\n");
  printf("-reduce int[2]       : reduce element order at material interval [int1 int2]\n");
  printf("-increase            : increase element order from linear to quadratic\n");
  printf("-bcoffset int        : add an offset to the boundary conditions\n");
  printf("-discont int         : make the boundary to have secondary nodes\n");
  printf("-connect int         : make the boundary to have internal connection among its elements\n");
  printf("-removeintbcs        : remove internal boundaries if they are not needed\n");
  printf("-removelowdim        : remove boundaries that are two ranks lower than highest dim\n");
  printf("-removeunused        : remove nodes that are not used in any element\n");
  printf("-bulkorder           : renumber materials types from 1 so that every number is used\n");
  printf("-boundorder          : renumber boundary types from 1 so that every number is used\n");
  printf("-autoclean           : this performs the united action of the four above\n");
  printf("-multidim            : keep lower order entities even if they are not boundaries\n");
  printf("-bulkbound int[3]    : set the intersection of materials [int1 int2] to be boundary int3\n");
  printf("-boundbound int[3]   : set the intersection of boundaries [int1 int2] to be boundary int3\n");
  printf("-bulktype int[3]     : set material types in interval [int1 int2] to type int3\n");
  printf("-boundtype int[3]    : set sidetypes in interval [int1 int2] to type int3\n");
  printf("-layer int[2] real[2]: make a boundary layer for given boundary\n");
  printf("-layermove int       : apply Jacobi filter int times to move the layered mesh\n");
  printf("-divlayer int[2] real[2]: make a boundary layer for given boundary\n");
  printf("-3d / -2d / -1d      : mesh is 3, 2 or 1-dimensional (applies to examples)\n");
  printf("-isoparam            : ensure that higher order elements are convex\n");
  printf("-nonames             : disable use of mesh.names even if it would be supported by the format\n");
  printf("-nosave              : disable saving part altogether\n");
  printf("-nooverwrite         : if mesh already exists don't overwrite it\n");
  printf("-vtuone              : start real node indexes in vtu file from one\n");
  printf("-timer               : show timer information\n");
  printf("-infofile str        : file for saving the timer and size information\n");

  printf("\nKeywords are related to mesh partitioning for parallel ElmerSolver runs:\n");
  printf("-partition int[3]    : the mesh will be partitioned in cartesian main directions\n");
  printf("-partorder real[3]   : in the 'partition' method set the direction of the ordering\n");
  printf("-parttol real        : in the 'partition' method set the tolerance for ordering\n");
  printf("-partcell int[3]     : the mesh will be partitioned in cells of fixed sizes\n");
  printf("-partcyl int[3]      : the mesh will be partitioned in cylindrical main directions\n");
#if USE_METIS
  printf("-metiskway int       : mesh will be partitioned with Metis using graph Kway routine\n");
  printf("-metisrec int        : mesh will be partitioned with Metis using graph Recursive routine\n");
  printf("-metiscontig         : enforce that the metis partitions are contiguous\n");
  printf("-metisvol            : minimize total communication volume in Metis\n");
  printf("-metisminconn        : minimize the maximum connectivity count in Metis\n");
  printf("-metisseed int       : random number generator seed for Metis algorithms\n");
  printf("-metisncuts int      : number of competing partitions to generate\n");
#endif
  printf("-partdual            : use the dual graph in partition method (when available)\n");
  printf("-halo                : create halo for the partitioning for DG\n");
  printf("-halobc              : create halo for the partitioning at boundaries only\n");
  printf("-haloz / -halor      : create halo for the the special z- or r-partitioning\n");
  printf("-halogreedy          : create halo being greedy over the partition interfaces\n");
  printf("-indirect            : create indirect connections (102 elements) in the partitioning\n");
  printf("-periodic int[3]     : periodic coordinate directions for parallel & conforming meshes\n");
  printf("-partoptim           : apply aggressive optimization to node sharing\n");
  printf("-partnobcoptim       : do not apply optimization to bc ownership sharing\n");
  printf("-partbw              : minimize the bandwidth of partition-partition couplings\n");
  printf("-parthypre           : number the nodes continuously partitionwise\n");
  printf("-partzbc             : partition connected BCs separately to partitions in Z-direction\n");
  printf("-partrbc             : partition connected BCs separately to partitions in R-direction\n");
#if USE_METIS
  printf("-metisbc             : partition connected BCs separately to partitions by Metis\n");
#endif
  printf("-partlayers int      : extend boundary partitioning by element layers\n");

  printf("\nKeywords are related to (nearly obsolete) ElmerPost format:\n");
  printf("-partjoin int        : number of ElmerPost partitions in the data to be joined\n");
  printf("-saveinterval int[3] : the first, last and step for fusing parallel data\n");
  printf("-nobound             : disable saving of boundary elements in ElmerPost format\n");

  if(0) printf("-names               : conserve name information where applicable\n");
}
 

void Goodbye()
{
  printf("\nThank you for using Elmergrid!\n");
  printf("Send bug reports and feature wishes to elmeradm@csc.fi\n");
  exit(0);
}



#if USE_MATC
char *mtc_domath(const char *);
void mtc_init(FILE * input, FILE * output, FILE * error);
#endif


static int Getline(char *line1,FILE *io) 
{
  int i,isend;
  char line0[MAXLINESIZE],*charend,*matcpntr,*matcpntr0;

  for(i=0;i<MAXLINESIZE;i++)
    line0[i] = 0x00;

 newline:

  charend = fgets(line0,MAXLINESIZE,io);
  isend = (charend == NULL);

  if(isend) return(1);

  if(line0[0] == '#' || line0[0] == '%' || line0[0] == '!') goto newline;
  if(!matcactive && line0[0] == '*') goto newline;

  if(!matcactive && strchr(line0,'$') ) {
    #if USE_MATC
      printf("\n a $ found and MATC has been compiled but not activated,\n");
      printf("either remove the $ or add 'MATC = True' to grd input file.\n");
    #else
      printf("\n a $ found and MATC has not been compiled into ElmerGrid,\n");
      printf("either remove the $ or compile ElmerGrid with MATC\n");
    #endif // USE_MATC
      bigerror("Error: $ found in command and MATC is not active");
  }

#if USE_MATC
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

  if(iodebug) printf("line: %s\n",line1);

  return(0);
}



void InitGrid(struct GridType *grid)
/* Initializes the grid of a specific mesh. A grid can differ 
   between various differential equations. 
   */
#define MAXBODYID 1000
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
  grid->lastmaterial = MAXBODYID;
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
    grid->zlastmaterial[i] = MAXBODYID;
    grid->zmaterial[i] = 0; 
  }

  grid->zmaterialmapexists = FALSE;
  grid->zhelicityexists = FALSE;
  grid->zhelicity = 0.0;

  /* Initializes the numbering of the cells. */
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


static void ExampleGrid1D(struct GridType **grids,int *nogrids,int info) 
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


static void ExampleGrid2D(struct GridType **grids,int *nogrids,int info) 
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

  
static void ExampleGrid3D(struct GridType **grids,int *nogrids,int info) 
/* Creates an example grid that might be used to analyze 
   a simple acceleration sensor. */
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


void CreateExampleGrid(int dim,struct GridType **grids,int *nogrids,int info) 
{
  switch (dim) {
    
  case 1: 
    ExampleGrid1D(grids,nogrids,info);
    break;

  case 2: 
    ExampleGrid2D(grids,nogrids,info);
    break;

  case 3: 
    ExampleGrid3D(grids,nogrids,info);
    break;
  }
}




void SetElementDivision(struct GridType *grid,Real relh,int info)
/* Given the densities and ratios in each cell finds the 
   optimum way to divide the mesh into elements. 
   The procedure is the following:
   For each subcell set the minimum number of elements
   then add one element at a time till the number of 
   elements is the desired one. The added element
   is always set to the subcell having the sparsest mesh.
   Also numbers the cells taking into consideration only 
   materials that have indices in interval [firstmat,lastmat].
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
  if(0) printf("The mesh is divided into %d separate subcells.\n",grid->nocells);
  
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
     Otherwise at least grid->minelems elements 
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



static int SetCellKnots(struct GridType *grid, struct CellType *cell,int info)
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
    printf("SetCellKnots: No numbering scheme for %d-node elements.\n",
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



static int SetCellKnots1D(struct GridType *grid, struct CellType *cell,int info)
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
    printf("SetCellKnots1D: No numbering scheme for %d-node elements.\n",
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
  int ind=0,aid=0,maxj=0;

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

  if(info) printf("Created %d divisions in %d cells for rotation [%.2lg  %.2lg].\n",
		  grid->totzelems,grid->zcells,
		  grid->rotateradius1,grid->rotateradius2);
  grid->dz0 = dzmax;
}



static int GetCommand(char *line1,char *line2,FILE *io) 
/* Line1 for commands and line2 for arguments. */
{
  int i,j,isend;
  char line0[MAXLINESIZE],*charend,*matcpntr0,*matcpntr;
  int gotlinefeed;

 newline:

  for(i=0;i<MAXLINESIZE;i++)
    line2[i] = line1[i] = line0[i] = 0x00;

  charend = fgets(line0,MAXLINESIZE,io);
  isend = (charend == NULL);

  if(isend) return(1);

  if(strlen(line0)<=strspn(line0," \t\r\n")) goto newline;

  if(line0[0] == '#' || line0[0] == '%' || line0[0] == '!' || line0[0] == '\n') goto newline;
  if(!matcactive && line0[0] == '*') goto newline;

  if(!matcactive && strchr(line0,'$') ) {
    #if USE_MATC
      printf("\n a $ found and MATC has been compiled but not activated,\n");
      printf("either remove the $ or add 'MATC = True' to grd input file.\n");
    #else
      printf("\n a $ found and MATC has not been compiled into ElmerGrid,\n");
      printf("either remove the $ or compile ElmerGrid with MATC\n");
    #endif // USE_MATC
      bigerror("Error: $ found in command and MATC is not active");
  }

#if USE_MATC
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
    if(line0[i] == '\n' || line0[i] == '\0' ) {
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
      if( line0[i] == '\n' || line0[i] == '\0' ) {
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
  else { /* Arguments are on the next line */
  newline2:
    charend = fgets(line2,MAXLINESIZE,io);
    isend = (charend == NULL);
    if(isend) return(2);
    if(line2[0] == '#' || line2[0] == '%' || line2[0] == '!') goto newline2;
    if(!matcactive && line2[0] == '*') goto newline2;

    gotlinefeed = FALSE;
    for(i=0;i<MAXLINESIZE;i++) {
      if(line2[i] == '\n' || line2[i] == '\0' ) {
	gotlinefeed = TRUE;
	break;
      }
    }
    if(!gotlinefeed) {
      printf("There is a risk that somethings was missed in line:\n");
      printf("%s\n",line2);
      smallerror("Check your output line length!\n");
    }

    if(!matcactive && strchr(line0,'$') ) {
      #if USE_MATC
        printf("\n a $ found and MATC has been compiled but not activated,\n");
        printf("either remove the $ or add 'MATC = True' to grd input file.\n");
      #else
        printf("\n a $ found and MATC has not been compiled into ElmerGrid,\n");
        printf("either remove the $ or compile ElmerGrid with MATC\n");
      #endif // USE_MATC
        bigerror("Error: $ found in command and MATC is not active");
    }

#if USE_MATC
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
    printf("command: %s\n",line1);
    printf("params:  %s\n",line2);
  }
  
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
    fprintf(out,"%.5lg ",grid->x[i]); 
  fprintf(out,"\n");
    
  if(dim >= 2) {
    fprintf(out,"Subcell Limits 2 %s",sameline ? "= ":"\n  ");
    for(i=0;i <= grid->ycells;i++) 
      fprintf(out,"%.5lg ",grid->y[i]); 
    fprintf(out,"\n");
  }
  
  if(dim >= 3) {
    fprintf(out,"Subcell Limits 3 %s",sameline ? "= ":"\n  ");
    for(i=0;i <= grid->zcells;i++) 
      fprintf(out,"%.5lg ",grid->z[i]); 
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
      fprintf(out,"  %-5d %-5d %-7.5lg %-7.5lg %-3d ",
	      grid->mappingtype[i],grid->mappingline[i],
	      grid->mappinglimits[2*i],grid->mappinglimits[2*i+1],
	      grid->mappingpoints[i]);
      for(j=0;j<grid->mappingpoints[i];j++) 
	fprintf(out,"%.4lg ",grid->mappingparams[i][j]);
      fprintf(out,"\n");
    }
    fprintf(out,"End\n");
  }

  j = 0;
  if(grid[j].rotate) {
    fprintf(out,"Revolve Blocks = %d\n",grid[j].rotateblocks);
    fprintf(out,"Revolve Radius = %-8.3lg\n",grid[j].rotateradius2);
    if(fabs(grid[j].rotateimprove-1.0) > 1.0e-10)
      fprintf(out,"Revolve Improve = %-8.3lg\n",grid[j].rotateimprove);
    
  }
  if(grid[j].rotatecurve) {
    fprintf(out,"Revolve Curve Direct = %-8.3lg\n",grid[j].curvezet);
    fprintf(out,"Revolve Curve Radius = %-8.3lg\n",grid[j].curverad);
    fprintf(out,"Revolve Curve Angle = %-8.3lg\n",grid[j].curveangle);
  }

  if(grid[j].coordsystem == COORD_POLAR) {
    fprintf(out,"Polar Radius = %.3lg\n",grid[j].polarradius);
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
      fprintf(out,"Coordinate Ratios = %-8.3lg\n",grid[j].xyratio);
    if(dim == 3)
      fprintf(out,"Coordinate Ratios = %-8.3lg %-8.3lg\n",
	      grid[j].xyratio,grid[j].xzratio);
 
    fprintf(out,"Minimum Element Divisions = %d",grid[j].minxelems);
    if(dim >= 2) fprintf(out," %d",grid[j].minyelems);
    if(dim >= 3) fprintf(out," %d",grid[j].minzelems);
    fprintf(out,"\n");

    fprintf(out,"Element Ratios 1 %s",sameline ? "= ":"\n  ");
    for(i=1;i<=grid[j].xcells;i++) 
      fprintf(out,"%.3lg ",grid[j].xexpand[i]); 
    fprintf(out,"\n");
    if(dim >= 2) {
      fprintf(out,"Element Ratios 2 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].ycells;i++) 
	fprintf(out,"%.3lg ",grid[j].yexpand[i]); 
      fprintf(out,"\n");
    }
    if(dim >= 3) {
      fprintf(out,"Element Ratios 3 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].zcells;i++) 
	fprintf(out,"%.3lg ",grid[j].zexpand[i]); 
      fprintf(out,"\n");
    }

    if(grid[j].autoratio) {
      fprintf(out,"Element Densities 1 %s",sameline ? "= ":"\n  ");
      for(i=1;i<=grid[j].xcells;i++) 
	fprintf(out,"%.3lg ",grid[j].xdens[i]); 
      fprintf(out,"\n");
      if(dim >= 2) {
	fprintf(out,"Element Densities 2 %s",sameline ? "= ":"\n  ");
	for(i=1;i<=grid[j].ycells;i++) 
	  fprintf(out,"%.3lg ",grid[j].ydens[i]); 
	fprintf(out,"\n");
      }
      if(dim >= 3) {       
	fprintf(out,"Element Densities 3 %s",sameline ? "= ":"\n  ");
	for(i=1;i<=grid[j].zcells;i++) 
	  fprintf(out,"%.3lg ",grid[j].zdens[i]); 
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





static int LoadElmergridOld(struct GridType **grid,int *nogrids,char *prefix,int info) 
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
    printf("LoadElmergrid: opening of the file '%s' wasn't successful !\n",filename);
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

      /* I don't know how to set this, luckily this piece of code should be obsolete */
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



int LoadElmergrid(struct GridType **grid,int *nogrids,char *prefix,int info) 
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
    printf("LoadElmergrid: opening of the file '%s' wasn't successful !\n",filename);
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
      if(info) printf("Reached the end of command file\n");
      goto end;
    }    

    /* Control information */
    if(strstr(command,"VERSION")) {
      if(strstr(command,"080500")) {
	if(info) printf("Loading old version of Elmergrid file.\n");
	i = LoadElmergridOld(grid,nogrids,prefix,info);
	return(i);
      }
      else {
	sscanf(params,"%ld",&code);
	if(code == 210903) 
	  if(info) printf("Loading ElmerGrid file version: %ld\n",code);
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
#if USE_MATC
	matcactive = TRUE;
	mtc_init(NULL, stdout, stderr);
	strcpy(command, "format( 12 )");	
	mtc_domath(command);	 
	printf("MATC language activated with 12 digit accuracy.\n");	
#else
        matcactive = FALSE;
        printf("Unable to enable MATC, as it is not even compiled.\n");
        bigerror("Unable to enable MATC, as it is not even compiled.");
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
      if(info) printf("Defining the coordinate system (%d-DIM).\n",grid[k]->dimension);
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

      /* Initialize the default structure with ones */
      for(j=grid[k]->ycells;j>=1;j--) 
	for(i=1;i<=grid[k]->xcells;i++) 
	  grid[k]->structure[j][i] = 1;
    }
    
    else if(strstr(command,"MINIMUM ELEMENT DIVISION")) {
      if(info) printf("Loading minimum number of elements\n");

      if((*grid)[k].dimension == 1) 
	sscanf(params,"%d",&(*grid)[k].minxelems); 

      if((*grid)[k].dimension == 2)  
	sscanf(params,"%d %d",&(*grid)[k].minxelems,&(*grid)[k].minyelems);

      if((*grid)[k].dimension == 3) 
	sscanf(params,"%d %d %d",&(*grid)[k].minxelems,
	       &(*grid)[k].minyelems,&(*grid)[k].minzelems);
    }      
    
    else if(strstr(command,"SUBCELL LIMITS 1")) {
      if(info) printf("Loading %d subcell limits in X-direction\n",grid[k]->xcells+1);
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
      if(info) printf("Loading %d subcell limits in Y-direction\n",grid[k]->ycells+1);
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
      if(info) printf("Loading %d subcell limits in Z-direction\n",grid[k]->zcells+1);
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
      if(info) printf("Loading %d subcell sizes in X-direction\n",grid[k]->xcells);
      cp = params;
      for(i=1;i<=grid[k]->xcells;i++) grid[k]->x[i] = next_real(&cp);
      for(i=1;i<=grid[k]->xcells;i++) grid[k]->x[i] = grid[k]->x[i-1] + grid[k]->x[i];
    }      
    else if(strstr(command,"SUBCELL SIZES 2")) {
      if(info) printf("Loading %d subcell sizes in Y-direction\n",grid[k]->ycells);
      cp = params;
      for(i=1;i<=grid[k]->ycells;i++) grid[k]->y[i] = next_real(&cp);
      for(i=1;i<=grid[k]->ycells;i++) grid[k]->y[i] = grid[k]->y[i-1] + grid[k]->y[i];
    }      
    else if(strstr(command,"SUBCELL SIZES 3")) {
      if(info) printf("Loading %d subcell sizes in Z-direction\n",grid[k]->zcells);
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
      if(info) printf("Loading material structure\n");

      /* Initialize the default structure with zeros */
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
      if(maxmat > MAXBODYID) 
	printf("LoadElmergrid: material indices larger to %d may create problems.\n",
	       MAXBODYID);
      printf("LoadElmergrid: materials interval is [%d,%d]\n",minmat,maxmat);

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
      if(strstr(params,"HORIZONTAL")) (*grid)[k].numbering = NUMBER_XY;
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

  
  /* for(k=nogrids0;k < (*nogrids) && k<MAXCASES;k++) {
    SetElementDivision(&(*grid)[k],relh,info);
    } */

  fclose(in);
  return(0);
}



void InitParameters(struct ElmergridType *eg)
{
  int i;

  eg->relh = 1.0;
  eg->inmethod = 0;
  eg->outmethod = 0;
  eg->silent = FALSE;
  eg->nofilesin = 1;
  eg->unitemeshes = FALSE;
  eg->unitenooverlap = FALSE;
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
  eg->metis_contig = FALSE;
  eg->metis_volcut = FALSE;
  eg->metis_seed = 0;
  eg->metis_ncuts = 1;
  eg->partopt = 0;
  eg->partoptim = FALSE;
  eg->partbcoptim = TRUE;
  eg->partjoin = 0;
  for(i=0;i<MAXHALOMODES;i++) eg->parthalo[i] = FALSE;
  eg->partitionindirect = FALSE;
  eg->reduce = FALSE;
  eg->increase = FALSE;
  eg->translate = FALSE;
  eg->isoparam = FALSE;
  eg->multidim = FALSE;
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
  eg->parttol = 0.0;
  
  for(i=0;i<MAXSIDEBULK;i++) 
    eg->sidebulk[i] = 0;
}




int InlineParameters(struct ElmergridType *eg,int argc,char *argv[],int first,int info)
{
  int arg,i,dim;
  char command[MAXLINESIZE];
  
  dim = eg->dim;

  printf("Elmergrid reading in-line arguments\n");

  /* Type of input file */
  if(first > 3) {
    for(i=0;i<MAXLINESIZE;i++) command[i] = ' ';

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
 
    /* Name of output file */
    strcpy(eg->filesin[0],argv[3]);
    strcpy(eg->filesout[0],eg->filesin[0]);
    strcpy(eg->infofile,eg->filesin[0]);
  }


  /* The optional inline parameters */

  for(arg=first;arg <argc; arg++) {

    if(strcmp(argv[arg],"-silent") == 0) {
      eg->silent = TRUE;
      info = FALSE;
    }
    else if(strcmp(argv[arg],"-verbose") == 0) {
      eg->silent = FALSE;
      info = TRUE;
    }   
    else if(strcmp(argv[arg],"-in") ==0 ) {
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
    else if(strcmp(argv[arg],"-out") == 0) {
      if(arg+1 >= argc) {
	printf("The output name is required as a parameter\n");
	return(2);
      }
      else {
	strcpy(eg->filesout[0],argv[arg+1]);
      }
    }
    else if(strcmp(argv[arg],"-decimals") == 0) {
      eg->decimals = atoi(argv[arg+1]);
    }
    else if(strcmp(argv[arg],"-triangles") ==0) {
      eg->triangles = TRUE;
      printf("The rectangles will be split to triangles.\n");
      if(arg+1 < argc) {
	if(strcmp(argv[arg+1],"-")) {
	  eg->triangleangle = atof(argv[arg+1]);
	}
      }
    }
    else if(strcmp(argv[arg],"-merge") == 0) {
      if(arg+1 >= argc) {
	printf("Give a parameter for critical distance.\n");
	return(3);
      }
      else {
	eg->merge = TRUE;
	eg->cmerge = atof(argv[arg+1]);
      }
    }
    else if(strcmp(argv[arg],"-relh") == 0) {
      if(arg+1 >= argc) {
	printf("Give a relative mesh density related to the specifications\n");
	return(3);
      }
      else {
	eg->relh = atof(argv[arg+1]);
      }
    }
    else if(strcmp(argv[arg],"-order") == 0) {
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
    else if(strcmp(argv[arg],"-parttol") == 0) {
      if(arg+1 >= argc) {
	printf("Give a tolerance for geometric partition algorithms\n");
	return(3);
      }
      else {
	eg->parttol = atof(argv[arg+1]);
      }
    }
    else if(strcmp(argv[arg],"-autoorder") == 0) {
      eg->order = 2;
    }
    else if(strcmp(argv[arg],"-halo") == 0) {
      eg->parthalo[1] = TRUE;
    }
    else if(strcmp(argv[arg],"-halobc") == 0) {
      eg->parthalo[2] = TRUE;
    }
    else if(strcmp(argv[arg],"-halodb") == 0) {
      eg->parthalo[1] = TRUE;
      eg->parthalo[2] = TRUE;            
    }   
    else if(strcmp(argv[arg],"-haloz") == 0) {
      eg->parthalo[3] = TRUE;
    }
    else if(strcmp(argv[arg],"-halor") == 0) {
      eg->parthalo[3] = TRUE;
    }
    else if(strcmp(argv[arg],"-halogreedy") == 0) {
      eg->parthalo[4] = TRUE;
    }    
    else if(strcmp(argv[arg],"-indirect") == 0) {
      eg->partitionindirect = TRUE;
    }
    else if(strcmp(argv[arg],"-metisorder") == 0) {
      eg->order = 3;
    }    
    else if(strcmp(argv[arg],"-centralize") == 0) {
      eg->center = TRUE;
    }
    else if(strcmp(argv[arg],"-scale") == 0) {
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
    else if(strcmp(argv[arg],"-translate") == 0) {
      if(arg+dim >= argc) {
	printf("Give %d parameters for the translate vector.\n",dim);
	return(6);
      }
      else {
	eg->translate = TRUE;
	eg->ctranslate[0] = atof(argv[arg+1]);
	eg->ctranslate[1] = atof(argv[arg+2]);
	if(dim==3) eg->ctranslate[2] = atof(argv[arg+3]);
      }
    }
    else if(strcmp(argv[arg],"-saveinterval") == 0) {
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
    else if(strcmp(argv[arg],"-rotate") == 0 || strcmp(argv[arg],"-rotate") == 0) {
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
    else if(strcmp(argv[arg],"-clone") == 0) {
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
    else if(strcmp(argv[arg],"-clonesize") == 0) {
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
    else if(strcmp(argv[arg],"-cloneinds") == 0) {
      eg->cloneinds = TRUE;
    }
    else if(strcmp(argv[arg],"-mirror") == 0) {
      if(arg+dim >= argc) {
	printf("Give the symmetry of the coordinate directions, eg. 1 1 0\n");
      }
      else {
	eg->mirror[0] = atoi(argv[arg+1]);
	eg->mirror[1] = atoi(argv[arg+2]);
	if(dim == 3) eg->mirror[2] = atoi(argv[arg+3]);
      }
    }
    else if(strcmp(argv[arg],"-mirrorbc") == 0) {
      if(arg+1 >= argc) {
	printf("Give the number of symmetry BC.\n");
 	return(11);
      }
      else {
	eg->mirrorbc = atoi(argv[arg+1]);
      }
    }
    else if(strcmp(argv[arg],"-unite") == 0) {
      eg->unitemeshes = TRUE;
      printf("The meshes will be united.\n");
    }   
    else if(strcmp(argv[arg],"-unitenooverlap") == 0) {
      eg->unitemeshes = TRUE;
      eg->unitenooverlap = TRUE;
      printf("The meshes will be united without overlap in BCs or bodies.\n");
    }   
    else if(strcmp(argv[arg],"-nonames") == 0) {
      eg->usenames = FALSE;
      printf("Names will be omitted even if they would exist\n");
    }   
    else if(strcmp(argv[arg],"-multidim") == 0) {
      eg->multidim = TRUE;
      printf("Lower dimensional entities may be bulk too!\n");
    }   
    else if(strcmp(argv[arg],"-removelowdim") == 0) {
      eg->removelowdim = TRUE;
      printf("Lower dimensional boundaries will be removed\n");
    }   
    else if(strcmp(argv[arg],"-removeintbcs") == 0) {
      eg->removeintbcs = TRUE;
      printf("Lower dimensional boundaries will be removed\n");
    }   
    else if(strcmp(argv[arg],"-removeunused") == 0) {
      eg->removeunused = TRUE;
      printf("Nodes that do not appear in any element will be removed\n");
    }   
    else if(strcmp(argv[arg],"-autoclean") == 0) {
      eg->removelowdim = TRUE;
      eg->bulkorder = TRUE;
      eg->boundorder = TRUE;
      eg->removeunused = TRUE;
      printf("Lower dimensional boundaries will be removed\n");
      printf("Materials and boundaries will be renumbered\n");
      printf("Nodes that do not appear in any element will be removed\n");
    }   
    else if(strcmp(argv[arg],"-polar") == 0) {
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
    else if(strcmp(argv[arg],"-cylinder") == 0) {
      eg->cylinder = TRUE;
      printf("Making transformation from cylindrical to cartesian coordinates.\n");
    }
    else if(strcmp(argv[arg],"-reduce") == 0) {
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
    else if(strcmp(argv[arg],"-increase") == 0) {
      eg->increase = TRUE;
    }
    else if(strcmp(argv[arg],"-bulkorder") == 0) {
      eg->bulkorder = TRUE;
    }
    else if(strcmp(argv[arg],"-boundorder") == 0) {
      eg->boundorder = TRUE;
    }
    else if(strcmp(argv[arg],"-partition") == 0  ||
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
    else if(strcmp(argv[arg],"-partorder") == 0) {
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
    else if(strcmp(argv[arg],"-partoptim") == 0) {
      eg->partoptim = TRUE;
      printf("Aggressive optimization will be applied to node sharing.\n");
    }
    else if(strcmp(argv[arg],"-partnobcoptim") == 0) {
      eg->partbcoptim = FALSE;
      printf("Aggressive optimization will not be applied to parent element sharing.\n");
    }
    else if(strcmp(argv[arg],"-partbw") == 0) {
      eg->partbw = TRUE;
      printf("Bandwidth will be optimized for partitions.\n");
    }
    else if(strcmp(argv[arg],"-parthypre") == 0) {
      eg->parthypre = TRUE;
      printf("Numbering of partitions will be made continuous.\n");
    }
    else if(strcmp(argv[arg],"-partdual") == 0) {
      eg->partdual = TRUE;
      printf("Using dual (elemental) graph in partitioning.\n");
    }
    else if(strcmp(argv[arg],"-metis") == 0 ||
       strcmp(argv[arg],"-metisrec") == 0 ||
       strcmp(argv[arg],"-metiskway") == 0 ) {
#if USE_METIS
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
    else if(strcmp(argv[arg],"-metisseed") == 0 ) {
      if(arg+1 >= argc) {
	printf("The random number seed is required as parameter for -metisseed!\n");
	return(15);
      }
      else {
	eg->metis_seed = atoi(argv[arg+1]);
	printf("Seed for Metis partitioning routines: %d\n",eg->metis_seed);
      }
    }
    else if(strcmp(argv[arg],"-metisncuts") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of parameters is required as parameter for -metisncuts!\n");
	return(15);
      }
      else {
	eg->metis_ncuts = atoi(argv[arg+1]);
	printf("Number of competing partitions to generate : %d\n",eg->metis_ncuts);
      }
    }   
    else if(strcmp(argv[arg],"-partjoin") == 0) {
      if(arg+1 >= argc) {
	printf("The number of partitions is required as a parameter!\n");
	return(15);
      }
      else {
	eg->partjoin = atoi(argv[arg+1]);
	printf("The results will joined using %d partitions.\n",eg->partjoin);
      }
    }
    else if(strcmp(argv[arg],"-partconnect") == 0 || strcmp(argv[arg],"-partzbc") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of 1D partitions is required as a parameter!\n");
	return(15);
      }
      else {
	eg->partbcz = atoi(argv[arg+1]);
	printf("The connected BCs will be partitioned to %d partitions in Z.\n",eg->partbcz);
      }
    }
    else if(strcmp(argv[arg],"-partrbc") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of 1D partitions is required as a parameter!\n");
	return(15);
      }
      else {
	eg->partbcr = atoi(argv[arg+1]);
	printf("The connected BCs will be partitioned to %d partitions in R.\n",eg->partbcr);
      }
    }
    else if(strcmp(argv[arg],"-partlayers") == 0) {
      if(arg+1 >= argc) {
	printf("The number of layers to be extended is required as a parameter\n");
	return(15);
      }
      else {
	eg->partbclayers = atoi(argv[arg+1]);
	printf("The boundary partitioning will be extended by %d layers.\n",eg->partbclayers);
      }
    }
    else if(strcmp(argv[arg],"-metiscontig") == 0 ) {
      eg->metis_contig = TRUE;
    }
    else if(strcmp(argv[arg],"-metisvol") == 0 ) {
      eg->metis_volcut = TRUE;
    }
    else if(strcmp(argv[arg],"-metisminconn") == 0 ) {
      eg->metis_minconn = TRUE;
    }    
    else if(strcmp(argv[arg],"-metisconnect") == 0 || strcmp(argv[arg],"-metisbc") == 0 ) {
      if(arg+1 >= argc) {
	printf("The number of Metis partitions is required as a parameter\n");
	return(15);
      }
      else {
	eg->partbcmetis = atoi(argv[arg+1]);
	printf("The connected BCs will be partitioned to %d partitions by Metis.\n",eg->partbcmetis);
      }
    }
    else if(strcmp(argv[arg],"-periodic") == 0) {
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
    else if(strcmp(argv[arg],"-discont") == 0) {
      if(arg+1 >= argc) {
	printf("Give the discontinuous boundary conditions.\n");
 	return(17);
      }
      else {
	eg->discontbounds[eg->discont] = atoi(argv[arg+1]);
	eg->discont++;
      }
    }
    else if(strcmp(argv[arg],"-connect") == 0) {
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
    else if(strcmp(argv[arg],"-connectall") == 0) {
      eg->connectboundsnosets += 1;
      eg->connectbounds[eg->connect] = -1;
      eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
      eg->connect++;
    }
    else if(strcmp(argv[arg],"-connectint") == 0) {
      eg->connectboundsnosets += 1;
      eg->connectbounds[eg->connect] = -2;
      eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
      eg->connect++;
    }
    else if(strcmp(argv[arg],"-connectfree") == 0) {
      eg->connectboundsnosets += 1;
      eg->connectbounds[eg->connect] = -3;
      eg->connectboundsset[eg->connect] = eg->connectboundsnosets;
      eg->connect++;
    }
    else if(strcmp(argv[arg],"-boundbound") == 0) {
      for(i=arg+1;i<=arg+3 && i<argc; i++) {
	eg->boundbound[3*eg->boundbounds+i-(1+arg)] = atoi(argv[i]);
	if((i-arg)%3 == 0) eg->boundbounds++;
      }
    } 
    else if(strcmp(argv[arg],"-bulkbound") == 0) {
      for(i=arg+1;i<=arg+3 && i<argc; i++) {
	eg->bulkbound[3*eg->bulkbounds+i-(1+arg)] = atoi(argv[i]);
	if((i-arg)%3 == 0) eg->bulkbounds++;
      }
    } 
    else if(strcmp(argv[arg],"-boundtype") == 0) {
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) 
	eg->sidemap[3*eg->sidemappings+i-1-arg] = atoi(argv[i]);
      eg->sidemappings++;
    } 
    else if(strcmp(argv[arg],"-bulktype") == 0) {
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) 
	eg->bulkmap[3*eg->bulkmappings+i-1-arg] = atoi(argv[i]);
      eg->bulkmappings++;
    } 
    else if(strcmp(argv[arg],"-coordinatemap") == 0) {
      if( arg+3 >= argc ) {
	printf("Give three parameters for the index permutation\n");
	return(18);
      }
      else {
	for(i=0;i<3;i++) 
	  eg->coordinatemap[i] = atoi(argv[arg+1+i]);
      }
    } 
    else if(strcmp(argv[arg],"-layer") == 0) {
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
    else if(strcmp(argv[arg],"-layermove") == 0) {
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
    else if(strcmp(argv[arg],"-divlayer") == 0) {
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
    else if(strcmp(argv[arg],"-3d") == 0) {
      eg->dim = dim = 3;
    }
    else if(strcmp(argv[arg],"-2d") == 0) {
      eg->dim = dim = 2;
    }
    else if(strcmp(argv[arg],"-1d") == 0) {
      eg->dim = dim = 1;
    }
    else if(strcmp(argv[arg],"-isoparam") == 0) {
      eg->isoparam = TRUE;
    }
    else if(strcmp(argv[arg],"-nobound") == 0) {
      eg->saveboundaries = FALSE;
    }
    else if(strcmp(argv[arg],"-vtuone") == 0) {
      eg->vtuone = TRUE;
    }
    else if(strcmp(argv[arg],"-nosave") == 0) {
      eg->nosave = TRUE;
    }
    else if(strcmp(argv[arg],"-nooverwrite") == 0) {
      eg->nooverwrite = TRUE;
    }
    else if(strcmp(argv[arg],"-timer") == 0) {
      eg->timeron = TRUE;
    }
    else if(strcmp(argv[arg],"-infofile") == 0) {
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
    else if(strcmp(argv[arg],"-bcoffset") == 0) {
      eg->bcoffset = atoi(argv[arg+1]);
    }
    else if(strcmp(argv[arg],"-noelements") == 0) {
      eg->elements3d = atoi(argv[arg+1]);
    }
    else if(strcmp(argv[arg],"-nonodes") == 0) {
      eg->nodes3d = atoi(argv[arg+1]);
    }
    else if(strcmp(argv[arg],"-sidefind") == 0) {
      eg->findsides = 0;
      for(i=arg+1;i<argc && strncmp(argv[i],"-",1); i++) {
	eg->sidebulk[i-1-arg] = atoi(argv[i]);
	eg->findsides++;
      }
    } 
    else if(strcmp(argv[arg],"-findbound") == 0) {
      eg->findsides = 0;
      for(i=arg+1;i+1<argc && strncmp(argv[i],"-",1); i += 2) {
	eg->sidebulk[i-1-arg] = atoi(argv[i]);
	eg->sidebulk[i-arg] = atoi(argv[i+1]);
	eg->findsides++;
      }
    }
    else if(strcmp(argv[arg],"-") == 0 ) {
      printf("Unknown in-line argument: %s\n",argv[arg]);
      bigerror("Cannot deal with argument, aborting!");
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
  int i,j,iostat;

  iodebug = FALSE;

  if( mode == 0) {  
    if (in = fopen("ELMERGRID_STARTINFO","r")) {
      iostat = fscanf(in,"%s",filename);
      fclose(in);
      printf("Using the file %s defined in ELMERGRID_STARTINFO\n",filename);
      if ((in = fopen(filename,"r")) == NULL) {
	printf("LoadCommands: opening of the file '%s' wasn't successful !\n",filename);
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
      printf("LoadCommands: opening of the file '%s' wasn't successful !\n",filename);
      return(3);
    }    
    if(info) printf("Loading ElmerGrid commands from file '%s'.\n",filename);    
  }
  else if(mode == 2) {
    AddExtension(prefix,filename,"grd");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadCommands: opening of the file '%s' wasn't successful !\n",filename);
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
    else if(strstr(command,"UNITENOOVERLAP")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) {
	eg->unitemeshes = TRUE;      
	eg->unitenooverlap = TRUE;      
      }
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
#if USE_METIS
      sscanf(params,"%d",&eg->partopt);
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"METIS")) {
#if USE_METIS
      sscanf(params,"%d",&eg->metis);
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"METISKWAY")) {
#if USE_METIS
      sscanf(params,"%d",&eg->metis);
      eg->partopt = 3;
#else
      printf("This version of ElmerGrid was compiled without Metis library!\n");
#endif
    }
    else if(strstr(command,"METISREC")) {
#if USE_METIS
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
      if(strstr(params,"TRUE")) eg->parthalo[1] = TRUE;      
    }
    else if(strstr(command,"BOUNDARY HALO")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->parthalo[2] = TRUE;
    }
    else if(strstr(command,"EXTRUDED HALO")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->parthalo[3] = TRUE;
    }
    else if(strstr(command,"GREEDY HALO")) {
      for(j=0;j<MAXLINESIZE;j++) params[j] = toupper(params[j]);
      if(strstr(params,"TRUE")) eg->parthalo[4] = TRUE;
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
  fclose(in);
  return(0);
}





static int FindParentSide(struct FemType *data,struct BoundaryType *bound,
			  int sideelem,int sideelemtype,int *sideind)
{
  int i,j,sideelemtype2,elemind,parent,normal,elemtype;
  int elemsides,side,sidenodes,nohits,hit,hit1,hit2;
  int sideind2[MAXNODESD1];
  int debug;
  
  hit = FALSE;
  elemsides = 0;
  elemtype = 0;
  hit1 = FALSE;
  hit2 = FALSE;

  debug = FALSE;
  
  for(parent=1;parent<=2;parent++) {
    if(parent == 1) 
      elemind = bound->parent[sideelem];
    else
      elemind = bound->parent2[sideelem];
    
    if(elemind > 0) {
      elemtype = data->elementtypes[elemind];
      elemsides = elemtype / 100;

      if(elemsides == 8) elemsides = 6;
      else if(elemsides == 7) elemsides = 5;
      else if(elemsides == 6) elemsides = 5;
      else if(elemsides == 5) elemsides = 4;
      
      for(normal=1;normal >= -1;normal -= 2) {

	for(side=0;side<elemsides;side++) {

	  if(debug) printf("elem = %d %d %d %d\n",elemind,elemsides,normal,side);

	  GetElementSide(elemind,side,normal,data,&sideind2[0],&sideelemtype2);
	  
	  if(sideelemtype2 < 300 && sideelemtype > 300) break;	
	  if(sideelemtype2 < 200 && sideelemtype > 200) break;		
	  if(sideelemtype != sideelemtype2) continue;

	  sidenodes = sideelemtype / 100;

	  for(j=0;j<sidenodes;j++) {
	    if(debug) printf("sidenode: %d %d %d\n",j,sideind[j],sideind2[j]);

	    hit = TRUE;
	    for(i=0;i<sidenodes;i++) 
	      if(sideind[(i+j)%sidenodes] != sideind2[i]) hit = FALSE;
	    
	    if(hit == TRUE) {
	      if(parent == 1) {
		hit1 = TRUE;
		bound->side[sideelem] = side;
		bound->normal[sideelem] = normal;
	      }
	      else {
		hit2 = TRUE;
		bound->side2[sideelem] = side;	      
	      }
	      goto skip;
	    }
	  }
	}
      }

      
      /* this finding of sides does not guarantee that normals are oriented correctly */
      normal = 1;
      hit = FALSE;
 
      for(side=0;;side++) {

	if(0) printf("side = %d\n",side);

	GetElementSide(elemind,side,normal,data,&sideind2[0],&sideelemtype2);

	if(sideelemtype2 == 0 ) break;
	if(sideelemtype2 < 300 && sideelemtype > 300) break;	
	if(sideelemtype2 < 200 && sideelemtype > 200) break;		
	if(sideelemtype != sideelemtype2) continue;
		
	sidenodes = sideelemtype % 100;

	nohits = 0;
	for(j=0;j<sidenodes;j++) 
	  for(i=0;i<sidenodes;i++) 
	    if(sideind[i] == sideind2[j]) nohits++;
	if(nohits == sidenodes) {
	  hit = TRUE;
	  if(parent == 1) {
	    hit1 = TRUE;
	    bound->side[sideelem] = side;
	  }
	  else {
	    hit2 = TRUE;
	    bound->side2[sideelem] = side;	      
	  }
	  goto skip;
	}
	
      }

    skip:  
      if(!hit) {
	printf("FindParentSide: cannot locate BC element in parent %d: %d\n",parent,elemind);
	printf("BC elem %d of type %d with corner indexes: ",sideelem,sideelemtype);
	for(i=0;i<sideelemtype/100;i++)
	  printf(" %d ",sideind[i]);
	printf("\n");

	printf("Bulk elem %d of type %d with corner indexes: ",elemind,elemtype);
	j = elemtype/100;
	if( j >= 5 && j<=7 ) j = j-1;
	for(i=0;i<j;i++)
	  printf(" %d ",data->topology[elemind][i]);
	printf("\n");             
      }
    }
  }

  if(hit1 || hit2) 
    return(0);
  else
    return(1);
}


static int Getnamerow(char *line,FILE *io,int upper) 
{
  int i,isend;
  char *charend;


  charend = fgets(line,MAXLINESIZE,io);
  isend = (charend == NULL);

  if(isend) 
    return(1);
  else
    return(0);
}



int LoadElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int nonames, int info)
/* This procedure reads the mesh assuming ElmerSolver format.
   */
{
  int noknots,noelements,nosides,maxelemtype,maxnodes,nonodes;
  int sideind[MAXNODESD1],tottypes,elementtype;
  int i,j,k,l,dummyint,cdstat,fail;
  int falseparents,noparents,bctopocreated;
  int activeperm,activeelemperm,mini,maxi,minelem,maxelem,p1,p2;
  int *nodeperm,*elemperm,*invperm,*invelemperm;
  int iostat,noelements0;
  FILE *in;
  char line[MAXLINESIZE],line2[MAXLINESIZE],filename[MAXFILESIZE],directoryname[MAXFILESIZE];
  char *ptr1,*ptr2;


  sprintf(directoryname,"%s",prefix);
  cdstat = chdir(directoryname);

  if(info) {
    if(cdstat) 
      printf("Loading mesh in ElmerSolver format from root directory.\n");
    else
      printf("Loading mesh in ElmerSolver format from directory %s.\n",directoryname);
  }

  InitializeKnots(data);


  sprintf(filename,"%s","mesh.header");
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadElmerInput: The opening of the header-file %s failed!\n",
	   filename);
    return(1);
  }
  else 
    printf("Loading header from %s\n",filename);

  GETLINE;
  sscanf(line,"%d %d %d",&noknots,&noelements,&nosides);
  GETLINE;
  sscanf(line,"%d",&tottypes);

  maxelemtype = 0;
  maxnodes = 0;
  for(i=1;i<=tottypes;i++) {   
    GETLINE;
    sscanf(line,"%d",&dummyint);
    maxelemtype = MAX( dummyint, maxelemtype );
    j = maxelemtype % 100;
    maxnodes = MAX( j, maxnodes );
  }
  printf("Maximum elementtype index is: %d\n",maxelemtype);
  printf("Maximum number of nodes in element is: %d\n",maxnodes);
  fclose(in);

  data->dim = GetElementDimension(maxelemtype);

  data->maxnodes = maxnodes;
  data->noknots = noknots;
  data->noelements = noelements0 = noelements;


  if(info) printf("Allocating for %d knots and %d elements.\n",
		  noknots,noelements);
  AllocateKnots(data);


  sprintf(filename,"%s","mesh.nodes");
  if ((in = fopen(filename,"r")) == NULL) {
    if(info) printf("LoadElmerInput: The opening of the nodes-file %s failed!\n",
		    filename);
    bigerror("Cannot continue without nodes file!\n");
  }
  else 
    printf("Loading %d Elmer nodes from %s\n",noknots,filename);

  activeperm = FALSE;
  for(i=1; i <= noknots; i++) {
    GETLINE;
    sscanf(line,"%d %d %le %le %le",
	   &j, &dummyint, &(data->x[i]),&(data->y[i]),&(data->z[i]));
    if(j != i && !activeperm) {
      printf("LoadElmerInput: The node number (%d) at node %d is not compact, creating permutation\n",j,i);
      activeperm = TRUE;
      nodeperm = Ivector(1,noknots);
      for(k=1;k<i;k++) nodeperm[k] = k;
    }
    if(activeperm) nodeperm[i] = j;
  }
  fclose(in);


  /* Create inverse permutation for nodes */
  if(activeperm) {
    for(i=1;i<=noknots;i++) {
      if(i==1) {
	mini = nodeperm[i];
	maxi = nodeperm[i];
      }
      else {
	mini = MIN(nodeperm[i],mini);
	maxi = MAX(nodeperm[i],maxi);
      }
    }
    if(info) printf("LoadElmerInput: Node index range is: [%d %d]\n",mini,maxi);
    invperm = Ivector(mini,maxi);
    for(i=mini;i<=maxi;i++)
      invperm[i] = -1;
    for(i=1;i<=noknots;i++) {
      j = nodeperm[i];
      if( invperm[j] > 0 ) 
	printf("LoadElmerInput: Node %d is redundant which may be problematic!\n",j);      
      else
	invperm[j] = i;
    }
  }
  else {
    mini = 1;
    maxi = noknots;
  }
  
  
  activeelemperm = FALSE;
  sprintf(filename,"%s","mesh.elements");
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadElmerInput: The opening of the element-file %s failed!\n",
	   filename);
    bigerror("Cannot continue without element file!\n");
  }
  else 
    if(info) printf("Loading %d bulk elements from %s\n",noelements,filename);
  
  for(i=1; i <= noelements; i++) {
    iostat = fscanf(in,"%d",&j);
    if(iostat <= 0 ) {
      printf("LoadElmerInput: Failed reading element line %d, reducing size of element table to %d!\n",i,i-1);
      data->noelements = noelements = i-1;
      break;
    }
    
    if(i != j && !activeelemperm) {
      printf("LoadElmerInput: The element numbering (%d) at element %d is not compact, creating permutation\n",j,i);
      activeelemperm = TRUE;
      elemperm = Ivector(1,noelements0);
      for(k=1; k < i; k++)
	elemperm[k] = k;
    }
    if( activeelemperm ) elemperm[i] = j;
    iostat = fscanf(in,"%d %d",&(data->material[i]),&elementtype);
    if( iostat < 2 ) {
      printf("LoadElmerInput: Failed reading definitions for bulk element %d\n",j);
      bigerror("Cannot continue without this data!\n");
    }
    if(elementtype > maxelemtype ) {
      printf("Invalid bulk elementtype: %d\n",elementtype);
      bigerror("Cannot continue with invalid elements");
    }
    data->elementtypes[i] = elementtype;
    nonodes = elementtype % 100;
    if( nonodes > maxnodes ) {
      printf("Number of nodes %d in element %d is greater than allocated maximum %d\n",nonodes,j,maxnodes);
      bigerror("Cannot continue with invalid elements");
    }
    for(k=0;k<nonodes;k++) {
      iostat = fscanf(in,"%d",&l);
      if( l < mini || l > maxi ) {
	printf("Node %d in element %d is out of range: %d\n",k+1,j,l);
	bigerror("Cannot continue with this node numbering");
      }
      if( activeperm )
	data->topology[i][k] = invperm[l];
      else
	data->topology[i][k] = l;
    }
  }
  fclose(in);
 

  /* Create inverse permutation for bulk elements */
  if(activeelemperm) {
    for(i=1;i<=noelements;i++) {
      if(i==1) {
	minelem = elemperm[i];
	maxelem = elemperm[i];
      }
      else {
	minelem = MIN(elemperm[i],minelem);
	maxelem = MAX(elemperm[i],maxelem);
      }
    }
    if(info) printf("LoadElmerInput: Element index range is: [%d %d]\n",minelem,maxelem);
    invelemperm = Ivector(minelem,maxelem);
    for(i=minelem;i<=maxelem;i++)
      invelemperm[i] = -1;
    for(i=1;i<=noelements;i++) {
      j = elemperm[i];
      if( invelemperm[j] > 0 )
	printf("LoadElmerInput: Element %d is redundant which may be problematic!\n",j);      
      else	
	invelemperm[j] = i;
    }
  }
  else {
    minelem = 1;
    maxelem = noelements;
  }



  falseparents = 0;
  noparents = 0;
  bctopocreated = FALSE;

  sprintf(filename,"%s","mesh.boundary");
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadElmerInput: The opening of the boundary-file %s failed!\n",
	   filename);
    return(4);
  }
  else {
    if(info) printf("Loading %d boundary elements from %s\n",nosides,filename);
  }

  if( nosides > 0 ) {
    AllocateBoundary(bound,nosides);
    data->noboundaries = 1;
  };

  i = 0;
  for(k=1; k <= nosides; k++) {
    
    iostat = fscanf(in,"%d",&dummyint);
    if( iostat < 1 ) {
      printf("LoadElmerInput: Failed reading boundary element line %d, reducing size of element table to %d!\n",k,i);
      bound->nosides = nosides = i;
      break;
    }      
    i++;

    iostat = fscanf(in,"%d %d %d %d",&(bound->types[i]),&p1,&p2,&elementtype);
    if(iostat < 4 ) {
      printf("LoadElmerInput: Failed reading definitions for boundary element %d\n",k);
      bigerror("Cannot continue without this data!\n"); 
    }    
    if( p1 > 0 && (p1 < minelem || p1 > maxelem ) ) {
      printf("Parent in boundary element %d out of range: %d\n",k,p1);    
      bigerror("Cannot continue with bad parents");
    }
    if( p2 > 0 && (p2 < minelem || p2 > maxelem ) ) {
      printf("Parent in boundary element %d out of range: %d\n",k,p2);
      bigerror("Cannot continue with bad parents");
    }
      
    if(activeelemperm) {
      if( p1 > 0 ) p1 = invelemperm[p1];
      if( p2 > 0 ) p2 = invelemperm[p2];
    }
    
    if(elementtype > maxelemtype ) {
      printf("Invalid boundary elementtype: %d\n",elementtype);
      bigerror("Cannot continue with invalid elements");
    }
    nonodes = elementtype % 100;
    if( nonodes > maxnodes ) {
      printf("Number of nodes %d in side element %d is greater than allocated maximum %d\n",nonodes,dummyint,maxnodes);
      bigerror("Cannot continue with invalid elements");
    }
    
    for(j=0;j< nonodes ;j++) { 
      iostat = fscanf(in,"%d",&l);
      if(activeperm) 
	sideind[j] = invperm[l];
      else
	sideind[j] = l;
    }
          
    if( p1 == 0 && p2 != 0 ) {
      bound->parent[i] = p2;
      bound->parent2[i] = p1;
    }
    else {
      bound->parent[i] = p1;
      bound->parent2[i] = p2;
    }
    
    if(bound->parent[i] > 0) {
      fail = FindParentSide(data,bound,i,elementtype,sideind);
      if(fail) falseparents++;      
    }
    else {
#if 0
      printf("Parents not specified for side %d with inds: ",dummyint);
      for(j=0;j< elementtype%100 ;j++) 
	printf("%d ",sideind[j]);
      printf("and type: %d\n",bound->types[i]);   
#endif
      if( !bctopocreated ) {
	bound->elementtypes = Ivector(1,nosides);
	for(j=1;j<=nosides;j++)
	  bound->elementtypes[j] = 0;
	bound->topology = Imatrix(1,nosides,0,data->maxnodes-1);
	bctopocreated = TRUE;
      }
      for(j=0;j< elementtype%100 ;j++) 
	bound->topology[i][j] = sideind[j];
      bound->elementtypes[i] = elementtype;

      printf("elementtype = %d %d %d\n",i,elementtype,sideind[0]);
      noparents++;
    }
  }
  
  if( falseparents ) {
    printf("There seems to be %d false parents in the mesh\n",falseparents);
  }
  if( noparents ) {
    printf("There seems to be %d bc elements without parents in the mesh\n",noparents);
  }

  bound->nosides = i;
  fclose(in); 
  
  /* Save node permutation for later use */
  data->nodepermexist = activeperm;
  if(activeperm) {
    data->nodeperm = nodeperm;
    free_Ivector(invperm,mini,maxi);
  }
  
  /* Element permutation is irrelevant probably for practical purposes (?) and hence it is forgotten. */
  if(activeelemperm) {
    free_Ivector(invelemperm,minelem,maxelem);
    free_Ivector(elemperm,1,noelements0);
  }
  


  sprintf(filename,"%s","mesh.names");
  if (in = fopen(filename,"r") ) {
    int isbody,started,nameproblem;
    
    isbody = TRUE;
    nameproblem = FALSE;

    if( nonames ) {
      printf("Ignoring > mesh.names < because it was explicitly requested!\n");
      goto namesend;
    }

    if(info) printf("Loading names for mesh parts from file %s\n",filename);

    for(;;) {
      if(Getnamerow(line,in,FALSE)) goto namesend;

      if(strstr(line,"names for boundaries")) {
	if(info) printf("Reading names for mesh boundaries\n");
	isbody = FALSE;
	continue;
      }
      else if(strstr(line,"names for bodies")) {
	if(info) printf("Reading names for mesh bodies\n");	
	isbody = TRUE;
	continue;
      }

      /* get position for entity name */
      ptr1 = strchr( line,'$');
      if(!ptr1) continue;
      ptr1++;

      /* get position for entity index and read it */
      ptr2 = strchr( line,'=');
      if(!ptr2) continue;
      ptr2++;      
      j = next_int(&ptr2);

      /* Initialize the entity name by white spaces */
      for(i=0;i<MAXLINESIZE;i++) 
	line2[i] = ' ';

      started = FALSE;
      k = 0;
      for(i=0;i<MAXLINESIZE;i++) {
	if( ptr1[0] == '=' ) {
	  /* remove possible trailing white space */
	  if(line2[k-1] == ' ') k--;	    
	  line2[k] = '\0';
	  break;
	}
	if(started || ptr1[0] != ' ') {
	  /* remove possible leading white space */
	  line2[k] = ptr1[0];
	  started = TRUE;
	  k++;
	}
	ptr1++;
      }

      /* Copy the entityname to mesh structure */
      if( isbody ) {
	if(j < 0 || j >= MAXBODIES ) {
	  if(!nameproblem) printf("Cannot treat names for bodies beyond %d\n",j);
	  nameproblem++;
	}
	else {
	  if(!data->bodyname[j]) data->bodyname[j] = Cvector(0,MAXNAMESIZE);
	  strcpy(data->bodyname[j],line2);	
	  data->bodynamesexist = TRUE;
	}
      }
      else {
	if(j < 0 || j >= MAXBCS ) {
	  if(!nameproblem) printf("Cannot treat names for boundaries beyond %d\n",j);
	  nameproblem++;
	}
	else {
	  if(!data->boundaryname[j]) data->boundaryname[j] = Cvector(0,MAXNAMESIZE);
	  strcpy(data->boundaryname[j],line2);	
	  data->boundarynamesexist = TRUE;
	}
      }
    }
    namesend:

    if( nameproblem ) {
      data->boundarynamesexist = FALSE;
      data->bodynamesexist = FALSE;
      printf("Number of indexes beyond range: %d\n",nameproblem);
      smallerror("Omitting use of names because some indexes are beyond range, code some more...");
    }
  }


  if(!cdstat) cdstat = chdir("..");

  if(info) printf("Elmer mesh loaded successfully\n");

  return(0);
}


int SaveElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int decimals,int nooverwrite, int info)
/* Saves the mesh in a form that may be used as input 
   in Elmer calculations. 
   */
{
  int noknots,noelements,material,sumsides,elemtype,fail,cdstat,bcdim;
  int sideelemtype,conelemtype,nodesd1,nodesd2,newtype;
  int i,j,k,l,bulktypes[MAXELEMENTTYPE+1],sidetypes[MAXELEMENTTYPE+1];
  int alltypes[MAXELEMENTTYPE+1],tottypes;
  int ind[MAXNODESD1],ind2[MAXNODESD1],usedbody[MAXBODIES],usedbc[MAXBCS];
  FILE *out;
  char filename[MAXFILESIZE], outstyle[MAXFILESIZE];
  char directoryname[MAXFILESIZE];

  if(!data->created) {
    printf("You tried to save points that were never created.\n");
    return(1);
  }

  noelements = data->noelements;
  noknots = data->noknots;
  sumsides = 0;

  for(i=0;i<=MAXELEMENTTYPE;i++)
    alltypes[i] = bulktypes[i] = sidetypes[i] = 0;

  for(i=0;i<MAXBODIES;i++)
    usedbody[i] = 0;
  for(i=0;i<MAXBCS;i++)
    usedbc[i] = 0;

  sprintf(directoryname,"%s",prefix);

  if(info) printf("Saving mesh in ElmerSolver format to directory %s.\n",
		  directoryname);

  fail = chdir(directoryname);
  if(fail) {
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
    if(info) printf("Reusing an existing directory\n");
    if(nooverwrite) {
      if (out = fopen("mesh.header", "r")) {
	printf("Mesh seems to already exist, writing is cancelled!\n"); 
	return(1);
      }
    }
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
  if(info) printf("Saving %d element topologies to %s.\n",noelements,filename);
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(3);
  }

  for(i=1;i<=noelements;i++) {
    elemtype = data->elementtypes[i];
    material = data->material[i];

    if(material < MAXBODIES) usedbody[material] += 1;
    fprintf(out,"%d %d %d",i,material,elemtype);

    bulktypes[elemtype] += 1;
    nodesd2 = elemtype%100;
    for(j=0;j < nodesd2;j++) 
      fprintf(out," %d",data->topology[i][j]);
    fprintf(out,"\n");          
  }
  fclose(out);


  sprintf(filename,"%s","mesh.boundary");
  out = fopen(filename,"w");
  if(info) printf("Saving boundary elements to %s.\n",filename);
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(4);
  }

  sumsides = 0;


  /* Save normal boundaries */
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound[j].created == FALSE) continue;
    if(bound[j].nosides == 0) continue;
    
    for(i=1; i <= bound[j].nosides; i++) {
      GetBoundaryElement(i,&bound[j],data,ind,&sideelemtype); 
      sumsides++;
      
      fprintf(out,"%d %d %d %d ",
	      sumsides,bound[j].types[i],bound[j].parent[i],bound[j].parent2[i]);
      fprintf(out,"%d",sideelemtype);

      k = bound[j].types[i];
      if(k < MAXBCS) {	
	bcdim = GetElementDimension(sideelemtype);
	usedbc[k] = MAX(usedbc[k],bcdim+1);
      }
	
      sidetypes[sideelemtype] += 1;
      nodesd1 = sideelemtype%100;
      for(l=0;l<nodesd1;l++)
	fprintf(out," %d",ind[l]);
      fprintf(out,"\n");
    }
  }

  newtype = 0;
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(bound[j].created == FALSE) continue;
    for(i=1; i <= bound[j].nosides; i++) 
      newtype = MAX(newtype, bound[j].types[i]);
  }  
  fclose(out);

  tottypes = 0;
  for(i=0;i<=MAXELEMENTTYPE;i++) {
    alltypes[i] = bulktypes[i] + sidetypes[i];
    if(alltypes[i]) tottypes++;
  }

  sprintf(filename,"%s","mesh.header");
  out = fopen(filename,"w");
  if(info) printf("Saving header info to %s.\n",filename);  
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(4);
  }
  fprintf(out,"%-6d %-6d %-6d\n",
	  noknots,noelements,sumsides);
  fprintf(out,"%-6d\n",tottypes);
  for(i=0;i<=MAXELEMENTTYPE;i++) {
    if(alltypes[i]) 
      fprintf(out,"%-6d %-6d\n",i,bulktypes[i]+sidetypes[i]);
  }
  fclose(out);


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
	if(usedbody[i]) {
	  if(data->bodyname[i])
	    fprintf(out,"$ %s = %d\n",data->bodyname[i],i);
	  else
	    fprintf(out,"$ body%d = %d\n",i,i);	    
	}
    }     
    if(data->boundarynamesexist) {
      fprintf(out,"! ----- names for boundaries -----\n");
      for(i=1;i<MAXBCS;i++) 	
	if(usedbc[i]) {
	  bcdim = usedbc[i]-1;
	  if(data->boundaryname[i]) 
	    fprintf(out,"$ %s = %d\n",data->boundaryname[i],i);
	  else if(bcdim == 2) { 
	    fprintf(out,"$ surf_bc%d = %d\n",i,i);
	  }
	  else if(bcdim == 1) {
	    fprintf(out,"$ line_bc%d = %d\n",i,i);
	  }
	  else if(bcdim == 0) {
	    fprintf(out,"$ node_bc%d = %d\n",i,i);
	  }
	}
    }
    fclose(out);

    sprintf(filename,"%s","entities.sif");
    out = fopen(filename,"w");
    if(info) printf("Saving entities info to %s.\n",filename);  
    if(out == NULL) {
      printf("opening of file was not successful\n");
      return(5);
    }

    if(data->bodynamesexist) {
      fprintf(out,"!------ Skeleton for body section -----\n");
      j = 0;
      for(i=1;i<MAXBODIES;i++) {
	if(usedbody[i]) {
	  j = j + 1;
	  fprintf(out,"Body %d\n",j);
	  if(data->bodyname[i])
	    fprintf(out,"  Name = %s\n",data->bodyname[i]);
	  else
	    fprintf(out,"  Name = body%d\n",i);	    
	  fprintf(out,"End\n\n");
	}
      }
    }

    if(data->boundarynamesexist) {
      fprintf(out,"!------ Skeleton for boundary section -----\n");
      j = 0;
      for(i=1;i<MAXBCS;i++) {
	if(usedbc[i]) {
	  j = j + 1;
	  fprintf(out,"Boundary Condition %d\n",j);
	  if(data->boundaryname[i])
	    fprintf(out,"  Name = %s\n",data->boundaryname[i]);
	  else
	    fprintf(out,"  Name = bc%d\n",i);	    
	  fprintf(out,"End\n\n");
	}
      }
    }
    fclose(out);
  }
  
  if(data->nodepermexist) {
    sprintf(filename,"%s","mesh.nodeperm");
    out = fopen(filename,"w");

    if(info) printf("Saving initial node permutation to %s.\n",filename);  
    if(out == NULL) {
      printf("opening of file was not successful\n");
      return(3);
    }
    for(i=1; i <= noknots; i++) 
      fprintf(out,"%d %d\n",i,data->nodeperm[i]);
  }


  cdstat = chdir("..");
  
  return(0);
}


int CreateElmerGridMesh(struct GridType *grid,
			struct FemType *data,struct BoundaryType *boundaries,
			Real relh,int info) {
  int j;  
  struct CellType *cell;
  
  for(j=0;j<MAXBOUNDARIES;j++) {
    boundaries[j].created = FALSE;
    boundaries[j].nosides = FALSE;
  }

  SetElementDivision(grid,relh,info);
  
  cell = (struct CellType*)
    malloc((size_t) (grid->nocells+1)*sizeof(struct CellType)); 
  SetCellData(grid,cell,info);

  if(grid->dimension == 1) 
    SetCellKnots1D(grid,cell,info);
  else
    SetCellKnots(grid,cell,info);

  CreateKnots(grid,cell,data,0,0);

  if(grid->noboundaries > 0) {
    for(j=0;j<grid->noboundaries;j++) {
      if(grid->boundsolid[j] < 4) {
	CreateBoundary(cell,data,&(boundaries[j]),grid->boundext[j],grid->boundint[j],
		       1,grid->boundtype[j],info);  
      } 
      else { 
	CreatePoints(cell,data,&(boundaries[j]),grid->boundext[j],grid->boundint[j],
		     grid->boundsolid[j],grid->boundtype[j],info); 	    
      }
    }
  }

  free(cell);

  return 0;
}

