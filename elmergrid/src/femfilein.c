/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter RÅÂback
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
#include <limits.h>
/*#include <strings.h>*/
/*#include <unistd.h>*/

#include "nrutil.h"
#include "common.h"
#include "femdef.h"
#include "femtypes.h"
#include "femknot.h"
#include "femfilein.h"

#define GETLINE getlineptr=fgets(line,MAXLINESIZE,in) 
#define GETLONGLINE getlineptr=fgets(longline,LONGLINESIZE,in)

static int linenumber;
static char *getlineptr;


static int Getrow(char *line1,FILE *io,int upper) 
{
  int i,isend;
  char line0[MAXLINESIZE],*charend;

  for(i=0;i<MAXLINESIZE;i++) 
    line0[i] = ' ';

 newline:
  charend = fgets(line0,MAXLINESIZE,io);
  linenumber += 1;

  isend = (charend == NULL);

  if(isend) return(1);

  if(line0[0] == '#' || line0[0] == '!') goto newline;
  if(strstr(line0,"#")) goto newline;

  if(upper) {
    for(i=0;i<MAXLINESIZE;i++) 
      line1[i] = toupper(line0[i]);
  }
  else {
    for(i=0;i<MAXLINESIZE;i++) 
      line1[i] = line0[i];    
  }

  return(0);
}

static int GetrowDouble(char *line1,FILE *io)
{
  int i,isend;
  char line0[MAXLINESIZE],*charend;

  for(i=0;i<MAXLINESIZE;i++) 
    line0[i] = ' ';

 newline:
  charend = fgets(line0,MAXLINESIZE,io);
  linenumber += 1;

  isend = (charend == NULL);

  if(isend) return(1);

  if(line0[0] == '#' || line0[0] == '!') goto newline;
  if(strstr(line0,"#")) goto newline;

  for(i=0;i<MAXLINESIZE;i++) { 

    /* The fortran double is not recognized by C string operators */
    if( line0[i] == 'd' || line0[i] == 'D' ) {
      line1[i] = 'e';
    } else {
      line1[i] = line0[i];    
    }
  }

  return(0);
}


static int Comsolrow(char *line1,FILE *io) 
{
  int i,isend;
  char line0[MAXLINESIZE],*charend;

  for(i=0;i<MAXLINESIZE;i++) 
    line0[i] = ' ';

  charend = fgets(line0,MAXLINESIZE,io);
  linenumber += 1;

  isend = (charend == NULL);

  if(isend) return(1);

  for(i=0;i<MAXLINESIZE;i++) line1[i] = line0[i];    

  return(0);
}



static void FindPointParents(struct FemType *data,struct BoundaryType *bound,
			    int boundarynodes,int *nodeindx,int *boundindx,int info)
{
  int i,j,k,sideelemtype,elemind,*indx;
  int boundarytype,minboundary,maxboundary,minnode,maxnode,sideelem,elemtype;
  int sideind[MAXNODESD1],elemsides,side,sidenodes,hit,nohits;
  int *elemhits;

  
  info = TRUE;

  sideelem = 0;
  maxboundary = minboundary = boundindx[1];
  minnode = maxnode = nodeindx[1]; 

  for(i=1;i<=boundarynodes;i++) {
    if(maxboundary < boundindx[i]) maxboundary = boundindx[i];
    if(minboundary > boundindx[i]) minboundary = boundindx[i];
    if(maxnode < nodeindx[i]) maxnode = nodeindx[i];
    if(minnode > nodeindx[i]) minnode = nodeindx[i];
  }

  if(info) {
    printf("Boundary types are in interval [%d, %d]\n",minboundary,maxboundary);
    printf("Boundary nodes are in interval [%d, %d]\n",minnode,maxnode);
  }

  indx = Ivector(1,data->noknots);

  printf("Allocating hit table of size: %d\n",data->noknots);
  elemhits = Ivector(1,data->noknots);
  for(i=1;i<=data->noknots;i++) elemhits[i] = 0;

  
  for(elemind=1;elemind<=data->noelements;elemind++) {
    elemtype = data->elementtypes[elemind];
    elemsides = elemtype % 100;
    
    for(i=0;i<elemsides;i++) {
      j = data->topology[elemind][i];
      elemhits[j] += 1;
    }
  }

  for(boundarytype=minboundary;boundarytype <= maxboundary;boundarytype++) {
    int boundfirst,bchits,bcsame,sideelemtype2;
    int sideind2[MAXNODESD1];

    boundfirst = 0;
  
    for(i=1;i<=data->noknots;i++) 
      indx[i] = 0;

    for(i=1;i<=boundarynodes;i++) {
      if(boundindx[i] == boundarytype) 
	indx[nodeindx[i]] = TRUE; 
    }


    for(elemind=1;elemind<=data->noelements;elemind++) {
      elemtype = data->elementtypes[elemind];
      elemsides = elemtype / 100;
      if(elemsides == 8) elemsides = 6;
      else if(elemsides == 5) elemsides = 4;
      else if(elemsides == 6 || elemsides == 7) elemsides = 5;
            
      /* Check whether the bc nodes occupy every node in the selected side */
      for(side=0;side<elemsides;side++) {
	GetElementSide(elemind,side,1,data,&sideind[0],&sideelemtype);
	sidenodes = sideelemtype%100;
		
	hit = TRUE;
	nohits = 0;
	for(i=0;i<sidenodes;i++) {
	  if(sideind[i] <= 0) { 
	    if(0) printf("sideind[%d] = %d\n",i,sideind[i]);
	    hit = FALSE;
	  }
	  else if(sideind[i] > data->noknots) { 
	    if(0) printf("sideind[%d] = %d (noknots=%d)\n",i,sideind[i],data->noknots);
	    hit = FALSE;
	  }
	  else if(!indx[sideind[i]]) {
	    hit = FALSE;
	  }
	  else {
	    nohits++;
	  }
	}
	
	if(0) {	  
	  printf("******\n");
	  printf("hits=%d  ind=%d  type=%d  sides=%d\n",nohits,elemind,elemtype,elemsides);
	  printf("sidenodes=%d  side=%d\n",sidenodes,side);

	  for(i=0;i<sidenodes;i++)
	    printf("%d  ",sideind[i]);

	  printf("\neleminds %d  %d\n",elemtype,elemind);
	  for(i=0;i<elemtype%100;i++)
	    printf("%d  ",data->topology[elemind][i]);
	  printf("\n");

	  hit = TRUE;
	}	  



	if(hit == TRUE) {

	  /* If all the points belong to more than one element there may be 
	     another parent for the side element */
	  bcsame = FALSE;
	  bchits = TRUE;
	  for(i=0;i<sidenodes;i++) 
	    if(elemhits[sideind[i]] < 2) bchits = FALSE;

	  if(bchits && boundfirst) {
	    for(j=boundfirst;j<=sideelem;j++) {
	      if(bound->parent2[j]) continue;
	      GetElementSide(bound->parent[j],bound->side[j],1,data,&sideind2[0],&sideelemtype2);
	      if(sideelemtype != sideelemtype2) continue;
	      bcsame = 0;
	      for(i=0;i<sidenodes;i++) 
		for(k=0;k<sidenodes;k++)
		  if(sideind[i] == sideind2[k]) bcsame++;

	      if(bcsame == sidenodes) {
		if(data->material[bound->parent[j]] > data->material[elemind]) {
		  bound->parent2[j] = bound->parent[j];
		  bound->side2[j] = bound->side[j];
		  bound->parent[j] = elemind;
		  bound->side[j] = side;
		}
		else {		  
		  bound->parent2[j] = elemind;
		  bound->side2[j] = side;
		}
		goto bcset;
	      }
	    }
	  }
		
	  sideelem += 1;

	  if( sideelem > bound->nosides ) {
	    printf("There are more side elements than allocated for (%d vs. %d)\n",sideelem,bound->nosides);
	  }
	  bound->parent[sideelem] = elemind;
	  bound->side[sideelem] = side;
	  bound->parent2[sideelem] = 0;
	  bound->side2[sideelem] = 0;
	  bound->types[sideelem] = boundarytype;

	  if(!boundfirst) boundfirst = sideelem;

	bcset:
	  continue;
	}
      }
    }    
  }

  
  free_Ivector(indx,1,data->noknots);

  if(info) printf("Found %d side elements formed by %d points.\n",
		  sideelem,boundarynodes);

  bound->nosides = MIN( sideelem, bound->nosides );

  return;
}




int LoadAbaqusInput(struct FemType *data,struct BoundaryType *bound,
		    char *prefix,int info)
/* Load the grid from a format that can be read by ABAQUS 
   program designed for sructural mechanics. The commands
   understood are only those that IDEAS creates when saving
   results in ABAQUS format.
   */
{
  int noknots,noelements,elemcode,maxnodes,material,maxelem,nodeoffset;
  int mode,allocated,nvalue,nvalue2,maxknot,nosides,elemnodes,ncum;
  int boundarytype,boundarynodes,elsetactive,elmatactive,cont;
  int *nodeindx=NULL,*boundindx=NULL,*materials=NULL,*elemindx=NULL;
  char *pstr;
  char filename[MAXFILESIZE];
  char line[MAXLINESIZE];
  int i,j,k,*ind=NULL;
  FILE *in;
  Real rvalues[MAXDOFS];
  int ivalues[MAXDOFS],ivalues0[MAXDOFS];
  int setmaterial;
  int debug,firstline;
  char entityname[MAXNAMESIZE];

  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"inp");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadAbaqusInput: opening of the ABAQUS-file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }

  printf("Reading input from ABAQUS input file %s.\n",filename);
  InitializeKnots(data);

  allocated = FALSE;
  maxknot = 0;
  maxelem = 0;
  elsetactive = FALSE;
  elmatactive = FALSE;
  
  /* Because the file format doesn't provide the number of elements
     or nodes the results are read twice but registered only in the
     second time. */

  debug = FALSE;
  
omstart:

  mode = 0;
  maxnodes = 0;
  noknots = 0;
  noelements = 0;
  nodeoffset = 0;
  elemcode = 0;
  boundarytype = 0;
  boundarynodes = 0;
  material = 0;
  ivalues0[0] = ivalues0[1] = 0;


  for(;;) {
    /* GETLINE; */

    if (Getrow(line,in,TRUE)) goto end;

    /* if(!line) goto end; */
    /* if(strstr(line,"END")) goto end; */


    if(strrchr(line,'*')) {
      if( mode == 2 ) {
	printf("Number of nodes so far: %d\n",noknots);
      }
      else if( mode == 3 ) {
	printf("Number of elements so far: %d\n",noelements);	
	if(elmatactive) {
	  nodeoffset = noknots;
	  printf("Node offset is %d\n",nodeoffset);
	}
      }
	
      if( strstr(line,"**") ) {
	if( strstr(line,"**HWCOLOR") )
	  mode = 10;
	else if( strstr(line,"**HWNAME") )
	  mode = 10;
	else if( strstr(line,"**HMASSEM") )
	  mode = 10;
	else if( strstr(line,"**HM_COMP") )
	  mode = 10;
	else if( strstr(line,"**HM_PROP") )
	  mode = 10;
	else
	  if(info && !allocated) printf("comment: %s",line);
      }
      else if(strstr(line,"HEAD")) {
	mode = 1;
      }
      else if(strstr(line,"*NODE")) {
	if(pstr = strstr(line,"NODE OUTPUT")) {
	  mode = 10;
	}
	else  {	  
	  if(strstr(line,"SYSTEM=R")) data->coordsystem = COORD_CART2;
	  if(strstr(line,"SYSTEM=C")) data->coordsystem = COORD_AXIS;      
	  if(strstr(line,"SYSTEM=P")) data->coordsystem = COORD_POLAR;      
	  mode = 2;
	}
      }
      else if(strstr(line,"*ELEMENT")) {
	if(pstr = strstr(line,"ELEMENT OUTPUT")) {
	  mode = 10;
	}
	else {
	  if(!(elsetactive || elmatactive)) material++;
	  if(strstr(line,"S3") || strstr(line,"STRI3") || strstr(line,"M3D3"))
	    elemcode = 303;
	  else if(strstr(line,"2D4") || strstr(line,"SP4") || strstr(line,"AX4") 
		  || strstr(line,"S4") || strstr(line,"CPE4")) 
	    elemcode = 404;
	  else if(strstr(line,"2D8") || strstr(line,"AX8") || strstr(line,"DS8") )
	    elemcode = 408;
	  else if(strstr(line,"3D4"))
	    elemcode = 504;
	  else if(strstr(line,"3D5"))
	    elemcode = 605;
	  else if(strstr(line,"3D6"))
	    elemcode = 706;
	  else if(strstr(line,"3D15"))
	    elemcode = 715;
	  else if(strstr(line,"3D8"))
	    elemcode = 808;
	  else if(strstr(line,"3D20"))
	    elemcode = 820;
	  else 
	    printf("Unknown element code: %s\n",line);

	  if(pstr = strstr(line,"ELSET=")) {
	    if(allocated) {
	      printf("Loading element set %d from %s",material,pstr+6);
	    }	    
	  }
	  
	  elemnodes = elemcode % 100;
	  maxnodes = MAX( maxnodes, elemnodes);
	  mode = 3;
	  if(allocated) {
	    printf("Loading elements of type %d starting from element %d.\n",
		   elemcode,noelements);
	    if(!(elsetactive || elmatactive)) {
	      sscanf(pstr+6,"%s",entityname);
	      strcpy(data->bodyname[material],entityname);
	      data->bodynamesexist = TRUE;
	      data->boundarynamesexist = TRUE;
	    }
	  }
	  
	  firstline = TRUE;
	}
      }
      else if( strstr(line,"BOUNDARY") ) {
	boundarytype++;
	mode = 4;
	if(allocated) {
	  printf("Treating keyword BOUNDARY\n");
	}
      }
      else if( strstr(line,"SOLID SECTION") ) {
	/* Have this here since solid section may include ELSET */
	mode = 10;
      }
      else if( strstr(line,"MEMBRANE SECTION") ) {
	/* Have this here since solid section may include ELSET */
	mode = 10;
      }
      else if( strstr(line,"CLOAD") ) {
	
	boundarytype++;
	mode = 4;
	if(allocated) {
	  printf("Treating keyword CLOAD\n");
	}
      }
       else if(pstr = strstr(line,"NSET=")) {
	 if( strstr(line,"ELSET=") ) {
	   /* Skipping association of ELSET to NSET */
	   mode = 10;
	 }
	 else {
	   boundarytype++;
	   mode = 5;	   
	   if(allocated) {
	     printf("Loading boundary node set %d from: %s",boundarytype,pstr+5);
	   }
	 }
      }
      else if(pstr = strstr(line,"ELSET=")) {
	elsetactive = TRUE;
	material += 1;
	mode = 6;

	if(allocated) {
	  printf("Loading element set %d from %s",material,pstr+6);
	  sscanf(pstr+6,"%s",entityname);
	  strcpy(data->bodyname[material],entityname);
	  data->bodynamesexist = TRUE;
	  data->boundarynamesexist = TRUE;
	}
      }
      else if(pstr = strstr(line,"PART, NAME=")) {
	elmatactive = TRUE;
	material += 1;
	mode = 6;

	if(allocated) {
	  printf("Loading part name %d from %s",material,pstr+11);
	  sscanf(pstr+6,"%s",entityname);
	  strcpy(data->bodyname[material],entityname);
	  data->bodynamesexist = TRUE;
	  data->boundarynamesexist = TRUE;
	}
      }
      else if(pstr = strstr(line,"HWCOLOR")) {
	/* unused command */
	mode = 0;
      }
      else {
	if(!allocated) printf("unknown command: %s",line);
	mode = 0;
      }
    }

    else if(mode) {  
      switch (mode) {
	
      case 1:
	if(info) printf("Loading Abaqus input file:\n%s",line);
	break;
	
      case 2: /* NODE */
	nvalue = StringToReal(line,rvalues,MAXNODESD2+1,',');

	if(nvalue != 4) {
	  printf("line: %s\n",line);
	  printf("Invalid nvalue = %d\n",nvalue);
	}
	  
	i = (int)(rvalues[0]+0.5);	
	noknots++;
	
	if(allocated) {
	  if( debug && (i==1 || i==maxknot) ) {
	    printf("debug node: %i %d %.3le %.3le %.3le\n",i,noknots,rvalues[1],rvalues[2],rvalues[3]);
	  }

	  i = MAX( i, noknots );
	  if(i <= 0 || i > maxknot) {
	    printf("Invalid node index = %d\n",i);
	  }
	  else {
	    ind[i] = noknots;
	    data->x[noknots] = rvalues[1];
	    data->y[noknots] = rvalues[2];
	    data->z[noknots] = rvalues[3];
	  }
	}
	else {
	  if(maxknot < i) maxknot = i;
	}
	break;
	
      case 3: /* ELEMENT */
	noelements++;
	
	nvalue = StringToIntegerNoZero(line,ivalues,elemnodes+1,',');

	if(allocated) {
	  if( debug && firstline ) {	  
	    printf("debug elem: %d %d %d %d\n",noelements,ivalues[0],elemcode,material);
	    printf("      topo:");
	    for(i=0;i<nvalue-1;i++)
	      printf(" %d",ivalues[i+1]);
	    printf("\n");
	    firstline = FALSE;
	  }
	  
	  elemindx[noelements] = ivalues[0];
	  data->elementtypes[noelements] = elemcode;

	  data->material[noelements] = material;
	  for(i=0;i<nvalue-1;i++) 
	    data->topology[noelements][i] = ivalues[i+1];	  

	  if( nodeoffset ) {
	    for(i=0;i<nvalue-1;i++) 
	      data->topology[noelements][i] += nodeoffset;
	  }

	}
	else {
	  if( maxelem < ivalues[0] ) maxelem = ivalues[0];
	}
	  
	ncum = nvalue-1;
	
	/* Read 2nd line if needed */
	if(ncum < elemnodes ) {
	  Getrow(line,in,TRUE);
	  nvalue = StringToIntegerNoZero(line,ivalues,elemnodes-ncum,',');
	  if(allocated) {
	    for(i=0;i<nvalue;i++) 
	      data->topology[noelements][ncum+i] = ivalues[i];	  	    
	  }
	  ncum = ncum + nvalue;
	}

	/* Be prepared for 3rd line as well */
	if(ncum < elemnodes ) {
	  Getrow(line,in,TRUE);
	  nvalue = StringToIntegerNoZero(line,ivalues,elemnodes-ncum,',');
	  if(allocated) {
	    for(i=0;i<nvalue;i++) 
	      data->topology[noelements][ncum+i] = ivalues[i];	  	    
	  }
	  ncum = ncum + nvalue;
	}
	if(ncum != elemnodes) printf("ncum = %d vs. %d\n",ncum,elemnodes);

	if( allocated ) {
	  j = FALSE;
	  for(i=0;i<elemnodes;i++)
	    if(!data->topology[noelements][i]) j = TRUE;

	  if(j) {
	    printf("zero in this element\n");
	    printf("element = %d %d\n",noelements,elemnodes);
	    for(i=0;i<elemnodes;i++)
	      printf("%d ",data->topology[noelements][i]);
	    printf("\n");
	  }
	}
	 	
	break;

      case 4:
	nvalue = StringToInteger(line,ivalues,2,',');

	if(ivalues[0] == ivalues0[0] && ivalues[1] != ivalues0[1]) continue;
	ivalues0[0] = ivalues[0];
	ivalues0[1] = ivalues[1];

	boundarynodes++;
	if(allocated) {
	  nodeindx[boundarynodes] = ivalues[0];
	  boundindx[boundarynodes] = boundarytype;
	}
	break;

      case 5: /* NSET */
	nvalue = StringToIntegerNoZero(line,ivalues,10,',');

	if(allocated) {
	  for(i=0;i<nvalue;i++) {
	    boundarynodes += 1;
	    nodeindx[boundarynodes] = ivalues[i];
	    boundindx[boundarynodes] = boundarytype;
	  }
	}
	else
	  boundarynodes += nvalue;
	break;

      case 6: /* ELSET */
	nvalue = StringToIntegerNoZero(line,ivalues,10,',');

	if(allocated) {
	  for(i=0;i<nvalue;i++) {
	    j = ivalues[i];
	    materials[j] = material;
	  }
	}
	break;

      case 10: 
	/* Doing nothing */
	break;
	

      default:
	printf("Unknown case: %d\n",mode);
      }      
    }

  }    
  end:


  if(allocated == TRUE) {
    int errcount,okcount;

    if(info) printf("The mesh was loaded from file %s.\n",filename);

    
    /* ABAQUS format does not expect that all numbers are used
       when numbering the elements. Therefore the nodes must
       be renumberred from 1 to noknots. */

    if(noknots != maxknot) {
      if(info) printf("There are %d nodes but maximum index is %d.\n",
		      noknots,maxknot);
      if(info) printf("Renumbering %d elements\n",noelements);
      errcount = 0;
      okcount = 0;
      for(j=1;j<=noelements;j++) {
	elemcode = data->elementtypes[j];
	elemnodes = elemcode % 100;
	for(i=0;i < elemnodes;i++)  {
	  k = data->topology[j][i];
	  if(k<=0) {
	    printf("err elem ind: %d %d %d %d\n",j,elemcode,i,k);
	    errcount++;
	  }
	  else {
	    data->topology[j][i] = ind[k];
	    okcount++;
	  }
	}
      }
      printf("There are %d positive and %d non-positive indexes in elements!\n",okcount,errcount);
      
      if(info) printf("Renumbering %d nodes in node sets\n",boundarynodes);
      errcount = 0;
      okcount = 0;
      for(j=1;j<=boundarynodes;j++) {
	k = nodeindx[j];
	if(k<=0 || k > maxknot) {
	  printf("err node set ind: %d %d\n",j,k);
	  errcount++;
	}
	else {
	  nodeindx[j] = ind[k];
	  okcount++;
	}
      }
      printf("There are %d positive and %d non-positive indexes in node sets!\n",okcount,errcount);
    }

    if(elsetactive) {
      for(i=1;i<=noelements;i++) {
	j = elemindx[i];
	data->material[i] = materials[j];
      }
    }
      
    ElementsToBoundaryConditions(data,bound,FALSE,info);
    
    free_ivector(ind,1,maxknot);
    free_ivector(materials,1,maxelem);
    free_Ivector(elemindx,1,noelements);

    if( boundarynodes > 0 ) {
      printf("Number of nodes in boundary sets: %d\n",boundarynodes);      
      free_Ivector(nodeindx,1,boundarynodes);
      free_Ivector(boundindx,1,boundarynodes);
    }
    
    fclose(in);

    return(0);
  }

  rewind(in);
  data->noknots = noknots;
  data->noelements = noelements;
  data->maxnodes = maxnodes;
  data->dim = 3; 
  
  if(info) printf("Allocating for %d knots and %d %d-node elements.\n",
		  noknots,noelements,maxnodes);
  AllocateKnots(data);
  
  elemindx = Ivector(1,noelements);
  for(i=1;i<=noelements;i++)
    elemindx[i] = 0;

  printf("Number of boundary nodes: %d\n",boundarynodes);
  if( boundarynodes > 0 ) {
    nodeindx = Ivector(1,boundarynodes);
    boundindx = Ivector(1,boundarynodes);
  }

  printf("Maximum element index in file: %d\n",maxelem);
  maxelem = MAX( maxelem, noelements );
  materials = ivector(1,maxelem);
  for(i=1;i<=maxelem;i++)
    materials[i] = 0;
 
  printf("Maximum node index in file: %d\n",maxknot);
  maxknot = MAX( maxknot, noknots );
  ind = ivector(1,maxknot);
  for(i=1;i<=maxknot;i++)
    ind[i] = 0;
    
  allocated = TRUE;
  goto omstart;
}


static int ReadAbaqusField(FILE *in,char *buffer,int *argtype,int *argno)
/* This subroutine reads the Abaqus file format and tries to make
   sence out of it. 
   */
{
  int i,val,digits;
  static int maxargno=0,mode=0;

  val = fgetc(in);

  if(val==EOF) return(-1);
  if(val=='\n') val = fgetc(in);

  if(val=='*') {
    if(0) printf("start field\n");
    if((*argno) != maxargno) 
      printf("The previous field was of wrong length, debugging time!\n");
    (*argno) = 0;
    mode = 0;
    val = fgetc(in);
    if(val=='\n') val = fgetc(in);
  }

  if(val=='I') {
    for(i=0;i<2;i++) {
      val = fgetc(in);
      if(val=='\n') val = fgetc(in);	
      buffer[i] = val;
    }
    buffer[2] = '\0';
    digits = atoi(buffer);
    for(i=0;i<digits;i++) {
      val = fgetc(in);
      if(val=='\n') val = fgetc(in);
      buffer[i] = val;
    }
    buffer[digits] = '\0';
    (*argno)++;
    (*argtype) = 1;
    if((*argno) == 1) maxargno = atoi(buffer);
    if((*argno) == 2) mode = atoi(buffer);
  }   
  else if(val=='D') {
    for(i=0;i<22;i++) {
      val = fgetc(in);
      if(val=='\n') val = fgetc(in);
      if(val=='D') val = 'E';
      buffer[i] = val;
    }
    buffer[22] = '\0';
    (*argno)++;
    (*argtype) = 2;
  }
  else if(val=='A') {
    for(i=0;i<8;i++) {
      val = fgetc(in);
      if(val=='\n') val = fgetc(in);
      buffer[i] = val;
    }
    buffer[8] = '\0';
    (*argno)++;
    (*argtype) = 3;
  }
  else {
    buffer[0] = val;
    buffer[1] = '\0';
    (*argtype) = 0;
  }
  return(mode);
}



int LoadAbaqusOutput(struct FemType *data,char *prefix,int info)
/* Load the grid from a format that can be read by ABAQUS 
   program designed for sructural mechanics. 
   */
{
  int knotno,elemno,elset,secno;
  int argtype,argno;
  int mode,allocated;
  char filename[MAXFILESIZE];
  char buffer[MAXLINESIZE];
  int i,j,prevdog,nodogs;
  int ignored;
  FILE *in=NULL;
  int *indx=NULL;


  AddExtension(prefix,filename,"fil");

  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadAbaqusOutput: opening of the Abaqus-file '%s' wasn't succesfull !\n",
	   filename);
    return(1);
  }
  if(info) printf("Reading input from ABAQUS output file %s.\n",filename);
  InitializeKnots(data);

  allocated = FALSE;
  mode = 0;
  knotno = 0;
  elemno = 0;
  nodogs = 0;
  prevdog = 0;
  argno = 0;
  elset = 0;
  secno = 0;
  ignored = 0;
  data->maxnodes = 9;
  data->dim = 3;

  for(;;) {

    mode = ReadAbaqusField(in,buffer,&argtype,&argno);
    if(0) printf("%d %d: buffer: %s\n",argtype,argno,buffer);

    switch (mode) {

    case -1:
      goto jump;
      
    case 0:
      break;

    case 1921:
      /* General info */
      if(argno == 3) printf("Reading output file for Abaqus %s\n",buffer);
      else if(argno == 4) printf("Created on %s",buffer);
      else if(argno == 5) printf("%s",buffer);
      else if(argno == 6) printf("%s\n",buffer);
      else if(argno == 7) data->noelements = atoi(buffer);
      else if(argno == 8 && allocated == FALSE) {
	data->noknots = atoi(buffer);
	allocated = TRUE;
	AllocateKnots(data);
	indx = Ivector(0,2 * data->noknots);
	for(i=1;i<=2*data->noknots;i++)
	  indx[i] = 0;
      }
      break;
      
    case 1900:
      /* Element definition */
      if(argno == 3) elemno = atoi(buffer);
      else if(argno == 4) {
	if(strstr(buffer,"2D4") || strstr(buffer,"SP4") || strstr(buffer,"AX4")) 
	  data->elementtypes[elemno] = 404;
	else if(strstr(buffer,"2D8") || strstr(buffer,"AX8") || strstr(buffer,"S8R5"))
	  data->elementtypes[elemno] = 408;
	else if(strstr(buffer,"3D8"))
	  data->elementtypes[elemno] = 808;
	else printf("Unknown element code: %s\n",buffer);
      } 
      else if(argno >= 5)
	data->topology[elemno][argno-5] = atoi(buffer);
      break;

    case 1901:
      /* Node definition */
      if(argno == 3) {
	knotno++;
	if(atoi(buffer) > 2*data->noknots) 
	  printf("LoadAbaqusOutput: allocate more space for indx.\n");
	else 
	  indx[atoi(buffer)] = knotno;
      }
      if(argno == 4) sscanf(buffer,"%le",&(data->x[knotno]));
      if(argno == 5) sscanf(buffer,"%le",&(data->y[knotno]));
      if(argno == 6) sscanf(buffer,"%le",&(data->z[knotno]));
      break;

    case 1933:
      /* Element set */
      if(argno == 3) { 
	elset++;
	strcpy(data->bodyname[elset],buffer);
      }
    case 1934:
      /* Element set continuation */
      if(argno > 3) {
	elemno = atoi(buffer);
	data->material[elemno] = elset;
      }
      break;

    case 2001:
      /* Just ignore */
      break;

    case 1:
      if(argno == 3) knotno = indx[atoi(buffer)];
      if(argno == 5) secno = atoi(buffer);
      break;

    case 2:
      if(prevdog != mode) {
	prevdog = mode;
	nodogs++;
	CreateVariable(data,nodogs,1,0.0,"Temperature",FALSE);
      }
      break;

    /* Read vectors in nodes in elements */
    case 11:
      if(prevdog != mode) {
	prevdog = mode;
	nodogs++;
	CreateVariable(data,nodogs,3,0.0,"Stress",FALSE);
      }
    case 12:
      if(prevdog != mode) {
	prevdog = mode;
	nodogs++;
	CreateVariable(data,nodogs,3,0.0,"Invariants",FALSE);
      }
      if(secno==1 && argno == 3) sscanf(buffer,"%le",&(data->dofs[nodogs][3*knotno-2]));
      if(secno==1 && argno == 4) sscanf(buffer,"%le",&(data->dofs[nodogs][3*knotno-1]));
      if(secno==1 && argno == 5) sscanf(buffer,"%le",&(data->dofs[nodogs][3*knotno]));
      break;

    /* Read vectors in nodes. */ 
    case 101:
      if(prevdog != mode) {
	prevdog = mode;
	nodogs++;
	CreateVariable(data,nodogs,3,0.0,"Displacement",FALSE);
      }
    case 102:
      if(prevdog != mode) {
	prevdog = mode;
	nodogs++;
	CreateVariable(data,nodogs,3,0.0,"Velocity",FALSE);
      }
      if(argno == 3) knotno = indx[atoi(buffer)]; 
      if(argno == 4) sscanf(buffer,"%le",&(data->dofs[nodogs][3*knotno-2]));
      if(argno == 5) sscanf(buffer,"%le",&(data->dofs[nodogs][3*knotno-1]));
      if(argno == 6) sscanf(buffer,"%le",&(data->dofs[nodogs][3*knotno]));
      break;

    default:
      if(ignored != mode) {
	printf("Record %d was ignored!\n",mode);
	ignored = mode;
      }
      break;
    }
  }

jump:

  if(info) printf("Renumbering elements\n");
  for(j=1;j<=data->noelements;j++) 
    for(i=0;i < data->elementtypes[j]%100;i++) 
      data->topology[j][i] = indx[data->topology[j][i]];

  free_ivector(indx,0,2*data->noknots);

  fclose(in);

  if(info) printf("LoadAbacusInput: results were loaded from file %s.\n",filename);

  return(0);
}




int LoadNastranInput(struct FemType *data,struct BoundaryType *bound,
		    char *prefix,int info)
/* Load the grid from a format that in Nastran format 
   */
{
  int noknots,noelements,maxnodes;
  int allocated,maxknot,minknot,nodes;
  char filename[MAXFILESIZE];
  char line[MAXLINESIZE],*cp;
  int j,k;
  FILE *in;


  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"nas");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadNastranInput: opening of the Nastran file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }

  if(info) printf("Reading mesh from Nastran file %s.\n",filename);
  InitializeKnots(data);

  allocated = FALSE;
  maxknot = 0;
  minknot = 10000;

  /* Because the file format doesn't provide the number of elements
     or nodes the results are read twice but registered only in the
     second time. */
omstart:

  maxnodes = 0;
  noknots = 0;
  noelements = 0;

  for(;;) {
    /* GETLINE; */

    if (Getrow(line,in,TRUE)) goto end;

    if(line[0] == '$') {
      if(!allocated) printf("comment: %s",line);
    }
    else if(strstr(line,"GRID")) {
      noknots++;
      if(0) printf("line=%s\n",line);
      cp = &line[5];
      j = next_int(&cp);
      if(0) printf("j=%d\n",j);

      if(allocated) {
	k = next_int(&cp);
	data->x[noknots] = next_real(&cp);
	data->y[noknots] = next_real(&cp);
      }
      else {
	if(j > maxknot) maxknot = j;
	if(j < minknot) minknot = j;
      }

      if(strstr(line,"*")) 
	Getrow(line,in,TRUE);

      if(allocated) {
	cp = &line[4];
	data->z[noknots] = next_real(&cp);
      }

    }
    else if(strstr(line,"TETRA")) {
      noelements++;
      nodes = 4;
      if(nodes > maxnodes) maxnodes = nodes;
      if(allocated) {
	data->elementtypes[noelements] = 504;
	cp = &line[6];
	k = next_int(&cp);
	data->material[noelements] = next_int(&cp) + 1;	
	for(j=0;j<nodes;j++)
	  data->topology[noelements][j] = next_int(&cp);            
      }      
    }
    else if(strstr(line,"PYRAM")) {
      noelements++;
      nodes = 5;
      if(nodes > maxnodes) maxnodes = nodes;
      if(allocated) {
	data->elementtypes[noelements] = 605;
	cp = &line[6];
	k = next_int(&cp);
	data->material[noelements] = next_int(&cp) + 1;	
	for(j=0;j<nodes;j++)
	  data->topology[noelements][j] = next_int(&cp);            
      }      
    }
    else if(strstr(line,"PENTA")) {
      noelements++;
      nodes = 6;
      if(nodes > maxnodes) maxnodes = nodes;
      if(allocated) {
	data->elementtypes[noelements] = 706;
	cp = &line[6];
	k = next_int(&cp);
	data->material[noelements] = next_int(&cp) + 1;	
	for(j=0;j<nodes;j++)
	  data->topology[noelements][j] = next_int(&cp);            
      }      
    }
     else if(strstr(line,"CHEXA")) {
      noelements++;
      nodes = 8;
      if(nodes > maxnodes) maxnodes = nodes;
      if(allocated) {
	data->elementtypes[noelements] = 808;
	cp = &line[5];
	k = next_int(&cp);
	data->material[noelements] = next_int(&cp) + 1;	
	for(j=0;j<6;j++)
	  data->topology[noelements][j] = next_int(&cp);            
      }           
      Getrow(line,in,TRUE);
      if(allocated) {
	cp = &line[1];       
	for(j=6;j<8;j++) 
	  data->topology[noelements][j] = k;
      }
    }
    else if(strstr(line,"ENDDAT")) {
      goto end;
    }
    else {
      printf("unknown command: %s",line);
    }
  }

  end:

  if(allocated == TRUE) {
    if(info) printf("The mesh was loaded from file %s.\n",filename);
    fclose(in);
    return(0);
  }

  rewind(in);
  data->noknots = noknots;
  data->noelements = noelements;
  data->maxnodes = maxnodes;
  data->dim = 3;
  
  printf("maxknot = %d  minknot = %d\n",maxknot,minknot);
  
  if(info) printf("Allocating for %d knots and %d %d-node elements.\n",
		  noknots,noelements,maxnodes);
  AllocateKnots(data);
    
  allocated = TRUE;
  goto omstart;
}




static void ReorderFidapNodes(struct FemType *data,int element,int nodes,int typeflag) 
{
  int i,oldtopology[MAXNODESD2],*topology,dim;
  int order808[]={1,2,4,3,5,6,8,7};
  int order408[]={1,3,5,7,2,4,6,8};
  int order306[]={1,3,5,2,4,6};
  int order203[]={1,3,2};
  int order605[]={1,2,4,3,5};

  dim = data->dim;
  if(typeflag > 10) dim -= 1;

  data->elementtypes[element] = 101*nodes;
  topology = data->topology[element];

  if(dim == 1) {
    if(nodes == 3) {
      data->elementtypes[element] = 203;
      for(i=0;i<nodes;i++) 
	oldtopology[i] = topology[i];
      for(i=0;i<nodes;i++) 
	topology[i] = oldtopology[order203[i]-1];            
    }
  }
  else if(dim == 2) {
    if(nodes == 6) {
      data->elementtypes[element] = 306;
      for(i=0;i<nodes;i++) 
	oldtopology[i] = topology[i];
      for(i=0;i<nodes;i++) 
	topology[i] = oldtopology[order306[i]-1];      
    }
    else if(nodes == 8) {
      data->elementtypes[element] = 408;
      for(i=0;i<nodes;i++) 
	oldtopology[i] = topology[i];
      for(i=0;i<nodes;i++) 
	topology[i] = oldtopology[order408[i]-1];      
    }
  }
  else if(dim == 3) {
    if(nodes == 4) {
      data->elementtypes[element] = 504;      
    }
    else if(nodes == 5) {
      data->elementtypes[element] = 605;      
      for(i=0;i<nodes;i++) 
	oldtopology[i] = topology[i];
      for(i=0;i<nodes;i++) 
	topology[i] = oldtopology[order605[i]-1];
    }
    else if(nodes == 6) {
      data->elementtypes[element] = 706;      
    }
    else if(nodes == 8) {
      for(i=0;i<nodes;i++) 
	oldtopology[i] = topology[i];
      for(i=0;i<nodes;i++) 
	topology[i] = oldtopology[order808[i]-1];
    }
    else {
      printf("Unknown Fidap elementtype with %d nodes.\n",nodes);
    }

  }
  else printf("ReorderFidapNodes: unknown dimension (%d)\n",data->dim);
}




int LoadFidapInput(struct FemType *data,char *prefix,int info)
/* Load the grid from a format that can be read by FIDAP 
   program designed for fluid mechanics. 

   Still under implementation
   */
{
  int noknots,noelements,dim,novel,maxnodes;
  int mode,maxknot,totelems,entity,maxentity;
  char filename[MAXFILESIZE];
  char line[MAXLINESIZE],entityname[MAXNAMESIZE];
  int i,j,k,*ind,geoflag,typeflag;
  int **topology;
  FILE *in;
  Real *vel,*temp;
  int nogroups,knotno;
  char *isio;

  AddExtension(prefix,filename,"fidap");

  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"FDNEUT");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadFidapInput: opening of the Fidap-file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }
 
  InitializeKnots(data);

  entity = 0;
  mode = 0;
  noknots = 0;
  noelements = 0;
  dim = 0;
  maxnodes = 4;
  totelems = 0;
  maxentity = 0;

  for(;;) {
    isio = GETLINE;

    if(!isio) goto end;
    if(strstr(line,"END")) goto end;

    /* Control information */
    if(strstr(line,"FIDAP NEUTRAL FILE")) mode = 1;
    else if(strstr(line,"NO. OF NODES")) mode = 2;
    else if(strstr(line,"TEMPERATURE/SPECIES FLAGS")) mode = 3;
    else if(strstr(line,"PRESSURE FLAGS")) mode = 4;
    else if(strstr(line,"NODAL COORDINATES")) mode = 5;
    else if(strstr(line,"BOUNDARY CONDITIONS")) mode = 6;
    else if(strstr(line,"ELEMENT GROUPS")) mode = 7;
    else if(strstr(line,"GROUP:")) mode = 8;
    else if(strstr(line,"VELOCITY")) mode = 10;
    else if(strstr(line,"TEMPERATURE")) mode = 11;
    else if(0) printf("unknown: %s",line);
    
    switch (mode) {

    case 1: 
      if(info) printf("Loading FIDAP input file %s\n",filename);
      GETLINE;
      if(info) printf("Name of the case: %s",line);
      mode = 0;
      break;

    case 2:
      GETLINE;   
      if(0) printf("reading the header info\n");
      sscanf(line,"%d%d%d%d%d",&noknots,&noelements,
	     &nogroups,&dim,&novel);
      data->noknots = noknots;
      data->noelements = noelements;
      data->maxnodes = maxnodes;
      data->dim = dim;

      mode = 0;
      break;
      
    case 5:
      if(info) printf("Allocating for %d knots and %d %d-node elements.\n",
		      noknots,noelements,maxnodes);
      AllocateKnots(data);
      if(info) printf("reading the nodes\n");
      for(i=1;i<=noknots;i++) {
	GETLINE;
	if (dim == 2)
	  sscanf(line,"%d%le%le",&knotno,
		 &(data->x[i]),&(data->y[i]));
	else if(dim==3) 
	  sscanf(line,"%d%le%le%le",&knotno,
		 &(data->x[i]),&(data->y[i]),&(data->z[i]));
      }
      break;
      
    case 8: 
      {
	int val,group,elems,nodes;
	char *cp;

	i=0;
	do val=line[i++];
	while(val!=':');i++;
	sscanf(&line[i],"%d",&group);
	
	do val=line[i++];
	while(val!=':');i++;
	sscanf(&line[i],"%d",&elems);
	
	do val=line[i++];
	while(val!=':');i++;
	sscanf(&line[i],"%d",&nodes);

	do val=line[i++];
	while(val!=':');i++;
	sscanf(&line[i],"%d",&geoflag);

	do val=line[i++];
	while(val!=':');i++;
	sscanf(&line[i],"%d",&typeflag);
	
	GETLINE;
	i=0;
	do val=line[i++];
	while(val!=':');i++;
	sscanf(&line[i],"%s",entityname);

	if(nodes > maxnodes) {
	  if(info) printf("Allocating a %d-node topology matrix\n",nodes);
	  topology = Imatrix(1,noelements,0,nodes-1);
	  if(totelems > 0) {
	    for(j=1;j<=totelems;j++) 
	      for(i=0;i<data->elementtypes[j] % 100;i++)
		topology[j][i] = data->topology[j][i];
	  }
	  free_Imatrix(data->topology,1,data->noelements,0,data->maxnodes-1);
	  data->maxnodes = maxnodes = nodes;
	  data->topology = topology;
	}

	if(0) printf("reading %d element topologies with %d nodes for %s\n",
			elems,nodes,entityname);

	for(entity=1;entity<=maxentity;entity++) {
#if 0
	  k = strcmp(entityname,entitylist[entity]);
#else
	  k = strcmp(entityname,data->bodyname[entity]);
#endif
	  if(k == 0) break;
	}

	if(entity > maxentity) {
	  maxentity++;
#if 0
	  strcpy(entitylist[entity],entityname);
#else
	  strcpy(data->bodyname[entity],entityname);
#endif
	  if(info) printf("Found new entity: %s\n",entityname);
	}

	for(i=totelems+1;i<=totelems+elems;i++) {
	  GETLINE;

	  cp = line;
	  j = next_int(&cp);
          for(j=0;j<nodes;j++)
	    data->topology[i][j] = next_int(&cp);

	  ReorderFidapNodes(data,i,nodes,typeflag);

	  if(data->elementtypes[i] == 0) {
	    printf("******** nolla\n");
	  }

	  if(entity) data->material[i] = entity;
	  else data->material[i] = group;
	}
	totelems += elems;
      }
    mode = 0;
    break;
      
    case 10:
      if(info) printf("reading the velocity field\n");
      CreateVariable(data,2,dim,0.0,"Velocity",FALSE);
      vel = data->dofs[2];
      for(j=1;j<=noknots;j++) {
	GETLINE;
	if(dim==2) 
	  sscanf(line,"%le%le",&(vel[2*j-1]),&(vel[2*j]));
	if(dim==3) 
	  sscanf(line,"%le%le%le",&(vel[3*j-2]),&(vel[3*j-1]),&(vel[3*j]));
      }
      mode = 0;
      break;
      
    case 11:

      if(info) printf("reading the temperature field\n");
      CreateVariable(data,1,1,0.0,"Temperature",FALSE);
      temp = data->dofs[1];
      for(j=1;j<=noknots;j++) {
	GETLINE;
	sscanf(line,"%le",&(temp[j]));
      }      
      mode = 0;
      break;
      
    default:      
      break;
    }
  }    

end:

  /* Renumber the nodes */
  maxknot = 0;
  for(i=1;i<=noelements;i++) 
    for(j=0;j < data->elementtypes[i]%100;j++) 
      if(data->topology[i][j] > maxknot) maxknot = data->topology[i][j];

  if(maxknot > noknots) {
    if(info) printf("renumbering the nodes from 1 to %d\n",noknots);

    ind = ivector(1,maxknot);
    for(i=1;i<=maxknot;i++)
      ind[i] = 0;

    for(i=1;i<=noelements;i++) 
      for(j=0;j < data->elementtypes[i]%100;j++) 
	ind[data->topology[i][j]] = data->topology[i][j];
    i=0;
    for(j=1;j<=noknots;j++) {
      i++;
      while(ind[i]==0) i++;
      ind[i] = j;
    } 
    for(i=1;i<=noelements;i++) 
      for(j=0;j < data->elementtypes[i]%100;j++) 
	data->topology[i][j] = ind[data->topology[i][j]];
  }


  if(maxentity > 0) data->bodynamesexist = TRUE;

  fclose(in);
  
  if(info) printf("Finished reading the Fidap neutral file\n");

  return(0);
}



static void ReorderAnsysNodes(struct FemType *data,int *oldtopology,
			      int element,int dim,int nodes) 
{
  int i,*topology,elementtype;
  int order820[]={1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16};
  int order504[]={1,2,3,5};
  int order306[]={1,2,3,5,6,8};
  int order510[]={1,2,3,5,9,10,12,17,18,19};
  int order613[]={1,2,3,4,5,9,10,11,12,17,18,19,20};
  int order706[]={1,2,3,5,6,7};
  int order715[]={1,2,3,5,6,7,9,10,12,17,18,19,13,14,16};

  elementtype = 0;
  if(dim == 3) {
    if(nodes == 20) {
      if(oldtopology[2] == oldtopology[3] &&
         oldtopology[4] == oldtopology[5]) elementtype = 510;
      else if(oldtopology[2] == oldtopology[3] &&
              oldtopology[4] != oldtopology[5]) elementtype = 715;
      else if(oldtopology[4] == oldtopology[5]) elementtype = 613;
      else elementtype = 820;
    }
    if(nodes == 8) {
      if(oldtopology[2] == oldtopology[3] &&
         oldtopology[4] == oldtopology[7] &&
         oldtopology[5] == oldtopology[7] &&
         oldtopology[6] == oldtopology[7]) elementtype = 504;
      else if(oldtopology[2] == oldtopology[3]  && 
              oldtopology[6] == oldtopology[7]) elementtype = 706;
      else if(oldtopology[4] == oldtopology[5]) elementtype = 605;
      else elementtype = 808;
    }
    if(nodes == 4) elementtype = 504;
    if(nodes == 10) elementtype = 510;
  }
  else if(dim == 2) {
    if(nodes == 9) elementtype = 408;
    if(nodes == 8) {
      if(oldtopology[3] == oldtopology[6]) 
	elementtype = 306;
      else
	elementtype = 408;
    }
    if(nodes == 4) elementtype = 404;
    if(nodes == 10) elementtype = 310;
    if(nodes == 6) elementtype = 306;
    if(nodes == 3) elementtype = 303;
  }
  else if(dim == 1) {
    if(nodes == 4) elementtype = 204;
    if(nodes == 3) elementtype = 203;
    if(nodes == 2) elementtype = 202;
  }

  if(!elementtype) {
    printf("Unknown elementtype in element %d (%d nodes, %d dim).\n",
	   element,nodes,dim);
  }

  data->elementtypes[element] = elementtype;
  topology = data->topology[element];

  switch (elementtype) {
 
  case 820:
    for(i=0;i<elementtype%100;i++) 
      topology[i] = oldtopology[order820[i]-1];
    break;

  case 504:
    if(nodes == 4)
      for(i=0;i<elementtype%100;i++)
         topology[i] = oldtopology[i];
    else
    for(i=0;i<elementtype%100;i++) 
      topology[i] = oldtopology[order504[i]-1];
    break;


  case 510:
    for(i=0;i<elementtype%100;i++) {
      if(oldtopology[2] == oldtopology[3]) 
	topology[i] = oldtopology[order510[i]-1];
      else 
	topology[i] = oldtopology[i];
    }	
    break;

  case 605:
    for(i=0;i<elementtype%100;i++) {
      topology[i] = oldtopology[i];
    }	
    break;

  case 613:
    for(i=0;i<elementtype%100;i++) {
      if(oldtopology[4] == oldtopology[5]) 
	topology[i] = oldtopology[order613[i]-1];
      else 
	topology[i] = oldtopology[i];
    }	
    break;

  case 306:
    for(i=0;i<elementtype%100;i++) {
      if(oldtopology[3] == oldtopology[6]) 
	topology[i] = oldtopology[order306[i]-1];    
      else 
	topology[i] = oldtopology[i];
    }
    break;

  case 706:
    for(i=0;i<elementtype%100;i++) {
	topology[i] = oldtopology[order706[i]-1];
    }	
    break;

  case 715:
    for(i=0;i<elementtype%100;i++) {
	topology[i] = oldtopology[order715[i]-1];
    }	
    break;

  default:
    for(i=0;i<elementtype%100;i++) 
      topology[i] = oldtopology[i];

  }  
}


int LoadAnsysInput(struct FemType *data,struct BoundaryType *bound,
		      char *prefix,int info)
/* This procedure reads the FEM mesh as written by Ansys. */
{
  int noknots=0,noelements=0,nosides,sidetype,currenttype;
  int maxindx,*indx,*revindx,topology[100],ind;
  int i,j,k,l,imax,*nodeindx,*boundindx,boundarynodes=0;
  int noansystypes,*ansysdim,*ansysnodes,*ansystypes,boundarytypes=0;
  int namesexist,maxside,sides;
  Real x,y,z;
  FILE *in;
  char *cp,line[MAXLINESIZE],filename[MAXFILESIZE],
    text[MAXNAMESIZE],text2[MAXNAMESIZE];


  /* ExportMesh.header */

  sprintf(filename,"%s.header",prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadAnsysInput: The opening of the header-file %s failed!\n",
	   filename);
    return(1);
  }

  if(info) printf("Calculating Ansys elementtypes from %s\n",filename);
  for(i=0;GETLINE;i++);

  noansystypes = i-1;
  printf("There seems to be %d elementytypes in file %s.\n",noansystypes,filename);
  
  ansysdim = Ivector(1,noansystypes);
  ansysnodes = Ivector(1,noansystypes);
  ansystypes = Ivector(1,noansystypes);

  rewind(in);
  for(i=0;i<=noansystypes;i++) {
    Real dummy1,dummy2,dummy3;
    GETLINE;

    /* Ansys writes decimal points also for integers and therefore these 
       values are read in as real numbers. */

    sscanf(line,"%le %le %le",&dummy1,&dummy2,&dummy3);

    if(i==0) {
      noelements = dummy1+0.5;
      noknots = dummy2+0.5;
      boundarytypes = dummy3+0.5;
    }
    else {
      ansysdim[i] = dummy1+0.5;
      ansysnodes[i] = dummy2+0.5;
      ansystypes[i] = dummy3+0.5;
    }
  }
  fclose(in);

  printf("Ansys file has %d elements, %d nodes and %d boundary types.\n",
	 noelements,noknots,boundarytypes);
  
 /* ExportMesh.names */

  sprintf(filename,"%s.names",prefix);
  in = fopen(filename,"r");
  if(in == NULL) 
    namesexist = FALSE;
  else
    namesexist = TRUE;

  if(namesexist) printf("Using names of bodies and boundaries\n");


  /* ExportMesh.node */

  sprintf(filename,"%s.node",prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadAnsysInput: The opening of the nodes-file %s failed!\n",
	   filename);
    return(2);
  }

  if(info) printf("Calculating Ansys nodes from %s\n",filename);
  for(i=0;GETLINE;i++);

  if(info) printf("There seems to be %d nodes in file %s.\n",i,filename);
  if(i != noknots) printf("Conflicting number of nodes %d vs %d!\n",i,noknots);

  /* Make room and initialize the mesh */
  InitializeKnots(data);

  /* Even 2D elements may form a 3D object! */
  data->dim = 3;
  data->maxnodes = 1;

  for(i=1;i<=noansystypes;i++) {
    if(ansysdim[i] > data->dim) data->dim = ansysdim[i];
    if(ansysnodes[i] > data->maxnodes) data->maxnodes = ansysnodes[i];
  }    
  if(data->maxnodes < 8) data->maxnodes = 8;

  data->noknots = noknots;
  data->noelements = noelements;

  if(info) printf("Allocating for %d nodes and %d elements with max. %d nodes in %d-dim.\n",
		  noknots,noelements,data->maxnodes,data->dim);
  AllocateKnots(data);

  indx = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) indx[i] = 0;

  if(info) printf("Loading %d Ansys nodes from %s\n",noknots,filename);
  rewind(in);
  for(i=1;i<=noknots;i++) {
    GETLINE; cp=line;

    indx[i] = next_int(&cp);
    if(cp[0] == '.') cp++;

    x = next_real(&cp);
    if(!cp) x = y = z = 0.0;
    else {
      y = next_real(&cp);
      if(!cp) y = z = 0.0;
      else if(data->dim == 3) {
	z = next_real(&cp);
	if(!cp) z = 0.0;
      }
    }
    data->x[i] = x;
    data->y[i] = y;
    if(data->dim == 3) data->z[i] = z;
  }
  fclose(in);
  
  /* reorder the indexes */
  maxindx = indx[1];
  for(i=1;i<=noknots;i++) 
    if(indx[i] > maxindx) maxindx = indx[i];
  revindx = Ivector(0,maxindx);

  if(maxindx > noknots) {
    printf("There are %d nodes with indexes up to %d.\n",
	   noknots,maxindx);
    for(i=1;i<=maxindx;i++) 
      revindx[i] = 0;
    for(i=1;i<=noknots;i++) 
      revindx[indx[i]] = i;
  }
  else {
    for(i=1;i<=noknots;i++) 
      revindx[i] = i;
  }

  /* ExportMesh.elem */

  sprintf(filename,"%s.elem",prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadAnsysInput: The opening of the element-file %s failed!\n",
	   filename);
    return(4);
  }
   
  if(info) printf("Loading %d Ansys elements from %s\n",noelements,filename);

  for(j=1;j<=noelements;j++) {

    GETLINE; 
    cp=line;

    for(i=0;i<8;i++) {
      ind = next_int(&cp);
      if(cp[0] == '.') cp++;
      topology[i] = revindx[ind];
    }

    data->material[j] = next_int(&cp);
    currenttype = next_int(&cp);
    
    for(k=1;k<=noansystypes;k++) 
      if(ansystypes[k] == currenttype) break;
    if(ansystypes[k] != currenttype) k=1;

    if(ansysnodes[k] > 8) {
      GETLINE; 
      cp=line;

      if(ansysnodes[k] == 10 && topology[2] != topology[3])
	imax = 10;
      else 
	imax = 20;

      for(i=8;i<imax;i++) {
	ind = next_int(&cp);
	if(cp[0] == '.') cp++;
	topology[i] = revindx[ind];
      }
    }
   
    ReorderAnsysNodes(data,&topology[0],j,ansysdim[k],ansysnodes[k]);
  }      
  fclose(in);


  /* ExportMesh.boundary */

  sprintf(filename,"%s.boundary",prefix);
  printf("Calculating nodes in file %s\n",filename);
  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadAnsysInput: The opening of the boundary-file %s failed!\n",
	   filename);
    return(5);
  }

  j = 0;
  for(i=0;GETLINE;i++)
    if(!strstr(line,"Boundary")) j++;

  boundarynodes = j;
  nosides = 6*boundarynodes;

  if(info) printf("There are %d boundary nodes, allocating %d elements\n",
		  boundarynodes,nosides);
  AllocateBoundary(bound,nosides);
  nodeindx = Ivector(1,boundarynodes);
  boundindx = Ivector(1,boundarynodes);

  if(info) printf("Loading Ansys boundary from %s\n",filename);

  for(i=1;i<=boundarynodes;i++) nodeindx[i] = boundindx[i] = 0;

  rewind(in);
  maxside = 0;
  j = 0;
  for(i=0;GETLINE;i++) {
    if(strstr(line,"Boundary")) {
      sscanf(line,"%d",&sidetype);
      maxside = MAX(sidetype,maxside);
    }
    else {
      j++;
      sscanf(line,"%d",&k);
      nodeindx[j] = revindx[k];
      if(!nodeindx[j]) printf("The boundary %dth node %d is not in index list\n",j,k);
      boundindx[j] = sidetype;
    }
  }
  printf("Found %d boundary nodes with %d as maximum side.\n",j,maxside);
  fclose(in);

  FindPointParents(data,bound,boundarynodes,nodeindx,boundindx,info);

  if(namesexist) {
    int bcind,*bctypes=NULL,*bctypeused=NULL,*bcused=NULL,newsides;

    data->bodynamesexist = TRUE;
    if(bound[0].nosides) {
      newsides = 0;
      bctypes = Ivector(1,maxside);
      bctypeused = Ivector(1,maxside);
      bcused = Ivector(1,bound[0].nosides);      
      for(i=1;i<=bound[0].nosides;i++) bcused[i] = FALSE;      
      data->boundarynamesexist = TRUE;
      for(i=1;i<=maxside;i++) 
	bctypeused[i] = FALSE;      
    }

    sprintf(filename,"%s.names",prefix);
    in = fopen(filename,"r");
            
    for(;;) {
      if(Getrow(line,in,TRUE)) break;
      sscanf(line,"%d%s%s%d",&bcind,&text[0],&text2[0],&sides);
      
      if(strstr(text2,"BODY")) {      
	GETLINE;
	sscanf(line,"%d%d",&j,&bcind);
	strcpy(data->bodyname[bcind],text);	
      }
      else if(strstr(text2,"BOUNDARY")) {      
	/* Read the boundary groups belonging to a particular name */
	for(i=1;i<=maxside;i++) 
	  bctypes[i] = 0;
	for(i=1;i<=sides;i++) {
	  GETLINE;
	  sscanf(line,"%d%d",&j,&bcind);
	  bctypes[bcind] = TRUE;
	}
	
	/* Find 1st unsed boundarytype */
	for(i=1;i<=maxside;i++) 
	  if(bctypes[i] && !bctypeused[i]) break;
	
	bcind = i;
	bctypeused[bcind] = TRUE;
	if(0) printf("First unused boundary is of type %d\n",bcind);
	strcpy(data->boundaryname[bcind],text);
	
	/* Check which of the BCs have already been named */
	k = l = 0;
	for(i=1;i<=bound[0].nosides;i++) {
	  j = bound[0].types[i];
	  
	  /* The bc is not given any name, hence it can't be a duplicate */
	  if(!bctypes[j]) continue;
	  
	  if(!bcused[i]) {
	    k++;
	    bcused[i] = bcind;
	  }
	  else {
	    l++;
	    if(newsides == 0) AllocateBoundary(&bound[1],bound[0].nosides);
	    newsides++;
	    bound[1].types[newsides] = bcind;
	    bound[1].parent[newsides] =  bound[0].parent[i];
	    bound[1].side[newsides] = bound[0].side[i];
	    bound[1].parent2[newsides] = bound[0].parent2[i];
	    bound[1].side2[newsides] = bound[0].side2[i];
	    bound[1].normal[newsides] = bound[0].normal[i];
	  }
	}
	if(info) printf("There are %d boundary elements with name %s.\n",k+l,data->boundaryname[bcind]);
      }
    }

    fclose(in);

    /* Put the indexes of all conditions with the same name to be same */
    if(bound[0].nosides) {
      for(i=1;i<=bound[0].nosides;i++) 
	if(bcused[i]) bound[0].types[i] = bcused[i];
      if(newsides) {
	bound[1].nosides = newsides;
	if(info) printf("Created %d additional boundary elements to achieve unique naming.\n",newsides);
      }
      free_Ivector(bctypes,1,maxside);    
      free_Ivector(bctypeused,1,maxside);    
      free_Ivector(bcused,1,bound[0].nosides);
    }
  }

  free_Ivector(boundindx,1,boundarynodes);
  free_Ivector(nodeindx,1,boundarynodes);

  if(info) printf("Ansys mesh loaded succefully\n");

  return(0);
}


static void ReorderFieldviewNodes(struct FemType *data,int *oldtopology,
				  int element,int dim,int nodes) 
{
  int i,*topology,elementtype;
  int order808[]={1,2,4,3,5,6,8,7};
  int order706[]={1,4,6,2,3,5};
  int order404[]={1,2,3,4};
    
  elementtype = 0;
  if(dim == 3) {
    if(nodes == 8) elementtype = 808;
    if(nodes == 6) elementtype = 706;
    if(nodes == 4) elementtype = 504;
  }
  else if(dim == 2) {
    if(nodes == 4) elementtype = 404;
    if(nodes == 3) elementtype = 303;
  }

  if(!elementtype) {
    printf("Unknown elementtype in element %d (%d nodes, %d dim).\n",
	   element,nodes,dim);
    bigerror("Cannot continue");
  }

  data->elementtypes[element] = elementtype;
  topology = data->topology[element];

  for(i=0;i<elementtype%100;i++) {
    if(elementtype == 808) topology[i] = oldtopology[order808[i]-1];
    else if(elementtype == 706) topology[i] = oldtopology[order706[i]-1];
    else if(elementtype == 404) topology[i] = oldtopology[order404[i]-1];
    else topology[i] = oldtopology[i];
  }
}




int LoadFieldviewInput(struct FemType *data,char *prefix,int info)
/* Load the grid from a format that can be read by FieldView
   program by PointWise. This is a suitable format to read files created
   by GridGen. */
{
  int noknots,noelements,maxnodes,mode;
  char filename[MAXFILESIZE];
  char line[MAXLINESIZE],*cp;
  int i,j,k;
  FILE *in;
  Real x,y,z;
  int maxindx;
  char *isio;
  int nobound=0,nobulk=0,maxsidenodes;
  int *boundtypes=NULL,**boundtopos=NULL,*boundnodes=NULL,*origtopology=NULL;

  if ((in = fopen(prefix,"r")) == NULL) {
    AddExtension(prefix,filename,"dat");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadFieldviewInput: opening of the Fieldview-file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }
 
  InitializeKnots(data);
  data->dim = 3;
  data->created = TRUE;

  mode = 0;
  noknots = 0;
  noelements = 0;
  maxnodes = 8;
  maxsidenodes = 4;
  maxindx = 0;

  data->maxnodes = maxnodes;

  for(;;) {
    
    if(mode == 0) {
      isio = GETLINE;
      
      if(!isio) goto end;
      if(strstr(line,"END")) goto end;
      
      /* Control information */
      if(strstr(line,"FIELDVIEW")) mode = 1;
      else if(strstr(line,"CONSTANTS")) mode = 2;
      else if(strstr(line,"GRIDS")) mode = 3;
      else if(strstr(line,"Boundary Table")) mode = 4;
      else if(strstr(line,"Variable Names")) mode = 5;
      else if(strstr(line,"Nodes")) mode = 6;
      else if(strstr(line,"Boundary Faces")) mode = 7;
      else if(strstr(line,"Elements")) mode = 8;
      else if(strstr(line,"Variables")) mode = 9;
      else if(0) printf("unknown: %s",line);
    }      

    switch (mode) {

    case 1: 
      printf("This is indeed a Fieldview input file.\n");
      mode = 0;
      break;
      
    case 2:
      if(0) printf("Constants block\n");
      mode = 0;
      break;

    case 3:
      if(0) printf("Grids block\n");
      mode = 0;
      break;

    case 4:
      if(0) printf("Boundary Table\n");
      mode = 0;
      break;

    case 5:
      if(0) printf("Variable names\n");
      mode = 0;
      break;

    case 6:

      GETLINE;   
      sscanf(line,"%d",&noknots);
      data->noknots = noknots;

      if(info) printf("Loading %d node coordinates\n",noknots);
      
      data->x = Rvector(1,noknots);
      data->y = Rvector(1,noknots);
      data->z = Rvector(1,noknots);
      
      for(i=1;i<=noknots;i++) {
	GETLINE;
	sscanf(line,"%le%le%le",&x,&y,&z);
	data->x[i] = x;
	data->y[i] = y;
	data->z[i] = z;
      }
      mode = 0;
      break;
      
    case 7:
      
      GETLINE;   
      sscanf(line,"%d",&nobound);

      if(info) printf("Loading %d boundary element definitions\n",nobound);

      boundtypes = Ivector(1,nobound);
      boundtopos = Imatrix(1,nobound,0,maxsidenodes-1);
      boundnodes = Ivector(1,nobound);
      
      for(i=1;i<=nobound;i++) {
	GETLINE; cp=line;
	
	boundtypes[i]= next_int(&cp);
	maxsidenodes = next_int(&cp);

	for(j=0;j<maxsidenodes && cp;j++) 
	  boundtopos[i][j] = next_int(&cp);

	boundnodes[i] = j;
      }
      mode = 0;
      break;
      
      
    case 8: 

      if(info) printf("Loading bulk element definitions\n");

      if(maxsidenodes == 4) noelements = noknots + nobound;
      else noelements = 6*noknots + nobound;

      origtopology = Ivector(0,maxnodes-1);
      data->topology = Imatrix(1,noelements,0,maxnodes-1);
      data->material = Ivector(1,noelements);
      data->elementtypes = Ivector(1,noelements);
      
      for(i=0;;) {
	GETLINE; cp=line;

	if(strstr(line,"Variables")) mode = 9;
	if(mode != 8) break;
	
	i++;

	k = next_int(&cp);

	if(k == 2) maxnodes = 8;
	else if(k == 1) maxnodes = 4;

	data->material[i]= next_int(&cp);

	for(j=0;j<maxnodes && cp;j++) {
	  origtopology[j] = next_int(&cp);
	  if(maxindx < origtopology[j]) maxindx = origtopology[j];
	}	

	ReorderFieldviewNodes(data,origtopology,i,3,j);

	if(nobulk+nobound == noelements+1) 
	  printf("Too few elements (%d) were allocated!!\n",noelements);

      }
      nobulk = i;

      printf("Found %d bulk elements\n",nobulk);
      if(nobulk+nobound > noelements) printf("Too few elements (%d) were allocated!!\n",noelements);
      printf("Allocated %.4g %% too many elements\n",
	     noelements*100.0/(nobulk+nobound)-100.0);


      for(i=1;i<=nobound;i++) {
	ReorderFieldviewNodes(data,boundtopos[i],i+nobulk,2,boundnodes[i]);
	data->material[i+nobulk] = boundtypes[i];
      }

      data->noelements = noelements = nobulk + nobound;

      mode = 0;
      break;


    case 9: 

      printf("Variables\n");
      
      mode = 0;
      break;
      
    default:      
      break;
    }
  }    
  
end:

  if(maxindx != noknots) 
    printf("The maximum index %d differs from the number of nodes %d !\n",maxindx,noknots);
  
  return(0);
}




int LoadTriangleInput(struct FemType *data,struct BoundaryType *bound,
		      char *prefix,int info)
/* This procedure reads the mesh assuming Triangle format
   */
{
  int noknots,noelements,maxnodes,elematts,nodeatts,dim;
  int elementtype,bcmarkers,sideelemtype;
  int i,j,k,*boundnodes;
  FILE *in;
  char *cp,line[MAXLINESIZE],elemfile[MAXFILESIZE],nodefile[MAXFILESIZE], 
    polyfile[MAXLINESIZE];
  int *invrow,*invcol;


  if(info) printf("Loading mesh in Triangle format from file %s.*\n",prefix);

  sprintf(nodefile,"%s.node",prefix);
  if ((in = fopen(nodefile,"r")) == NULL) {
    printf("LoadTriangleInput: The opening of the nodes file %s failed!\n",nodefile);
    return(1);
  }
  else 
    printf("Loading nodes from file %s\n",nodefile);

  GETLINE;
  sscanf(line,"%d %d %d %d",&noknots,&dim,&nodeatts,&bcmarkers);
  fclose(in);

  if(dim != 2) {
    printf("LoadTriangleInput assumes that the space dimension is 2, not %d.\n",dim);
    return(2);
  }

  sprintf(elemfile,"%s.ele",prefix);
  if ((in = fopen(elemfile,"r")) == NULL) {
    printf("LoadTriangleInput: The opening of the element file %s failed!\n",elemfile);
    return(3);
  }
  else 
    printf("Loading elements from file %s\n",elemfile);

  GETLINE;
  sscanf(line,"%d %d %d",&noelements,&maxnodes,&elematts);
  fclose(in);


  InitializeKnots(data);
  data->dim = dim;
  data->maxnodes = maxnodes;
  data->noelements = noelements;
  data->noknots = noknots;
  elementtype = 300 + maxnodes;

  if(info) printf("Allocating for %d knots and %d elements.\n",noknots,noelements);
  AllocateKnots(data);

  boundnodes = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) 
    boundnodes[i] = 0;

  in = fopen(nodefile,"r");
  GETLINE;
  for(i=1; i <= noknots; i++) {
    GETLINE;
    cp = line;
    j = next_int(&cp);
    if(j != i) printf("LoadTriangleInput: nodes i=%d j=%d\n",i,j);
    data->x[i] = next_real(&cp);
    data->y[i] = next_real(&cp);
    for(j=0;j<nodeatts;j++) {
      next_real(&cp);
    }
    if(bcmarkers > 0) 
      boundnodes[i] = next_int(&cp);
  }
  fclose(in);

  in = fopen(elemfile,"r");
  GETLINE;
  for(i=1; i <= noelements; i++) {
    GETLINE;
    cp = line;
    data->elementtypes[i] = elementtype;
    j = next_int(&cp);
    if(j != i) printf("LoadTriangleInput: elem i=%d j=%d\n",i,j);
    for(j=0;j<3;j++)
      data->topology[i][j] = next_int(&cp);
    if(maxnodes == 6) {
      data->topology[i][4] = next_int(&cp);
      data->topology[i][5] = next_int(&cp);
      data->topology[i][3] = next_int(&cp);
    }
    data->material[i] = 1;
  }
  fclose(in);


  sprintf(polyfile,"%s.poly",prefix);
  if ((in = fopen(polyfile,"r")) == NULL) {
    printf("LoadTriangleInput: The opening of the poly file %s failed!\n",polyfile);
    return(1);
  }
  else 
    printf("Loading nodes from file %s\n",polyfile);

  {
    int bcelems,markers,ind1,ind2,bctype,j2,k2,hit;
    int elemsides,sideind[2],side,elemind=0;

    bctype = 1;
    elemsides = 3;
    hit = FALSE;

    GETLINE;
    GETLINE;
    sscanf(line,"%d %d",&bcelems,&markers);

    CreateInverseTopology(data,info);
    invrow = data->invtopo.rows;
    invcol = data->invtopo.cols;

    AllocateBoundary(bound,bcelems);

    for(i=1;i<=bcelems;i++) {
      
      GETLINE;
      if(markers)
	sscanf(line,"%d %d %d %d",&j,&ind1,&ind2,&bctype);
      else 
	sscanf(line,"%d %d %d",&j,&ind1,&ind2);
     
      /* find an element which owns both the nodes */
#if 0
      for(j=1;j<=data->maxinvtopo;j++) {
	hit = FALSE;
	k = data->invtopo[j][ind1];
	if(!k) break;

	for(j2=1;j2<=data->maxinvtopo;j2++) { 
	  k2 = data->invtopo[j2][ind2];
	  if(!k2) break;
	  if(k == k2) {
	    hit = TRUE;
	    elemind = k;
	    break;
	  }
	}
	if(hit) break;
      }
#else
      for(j=invrow[ind1-1];j<invrow[ind1];j++) {
	k = invcol[j]+1;
	hit = FALSE;

	for(j2=invrow[ind2-1];j2<invrow[ind2];j2++) {
	  k2 = invcol[j2]+1;
	  if(k == k2) {
	    hit = TRUE;
	    elemind = k;
	    break;
	  }
	}
	if(hit) break;
      }
#endif


      if(!hit) return(1);


      /* Find the correct side of the triangular element */
      for(side=0;side<elemsides;side++) {
	GetElementSide(elemind,side,1,data,&sideind[0],&sideelemtype);
	
	hit = FALSE;
	if(sideind[0] == ind1 && sideind[1] == ind2) hit = TRUE;
	if(sideind[0] == ind2 && sideind[1] == ind1) hit = TRUE;

	if(hit) {
	  bound->parent[i] = elemind;
	  bound->side[i] = side;
	  bound->parent2[i] = 0;
	  bound->side2[i] = 0;
	  bound->types[i] = bctype;
	}
      }
    }
  } 

  printf("Successfully read the mesh from the Triangle input file.\n");

  return(0);
}




int LoadMeditInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int info)
/* This procedure reads the mesh assuming Medit format
   */
{
  int noknots,noelements,maxnodes,dim=0,elementtype;
  int i,j,allocated;
  FILE *in;
  char *cp,line[MAXLINESIZE],nodefile[MAXFILESIZE];


  sprintf(nodefile,"%s.mesh",prefix);
  if(info) printf("Loading mesh in Medit format from file %s\n",prefix);

  if ((in = fopen(nodefile,"r")) == NULL) {
    printf("LoadMeditInput: The opening of the mesh file %s failed!\n",nodefile);
    return(1);
  }

  allocated = FALSE;
  maxnodes = 0;

allocate:

  if(allocated) {
    InitializeKnots(data);
    data->dim = dim;
    data->maxnodes = maxnodes;
    data->noelements = noelements;
    data->noknots = noknots;
    elementtype = 300 + maxnodes;

    if(info) printf("Allocating for %d knots and %d elements.\n",noknots,noelements);
    AllocateKnots(data);
    in = fopen(nodefile,"r");
  }


  for(;;) {
    if(Getrow(line,in,TRUE)) goto end;
    if(strstr(line,"END")) goto end;
    
    if(strstr(line,"DIMENSION")) {
      if(Getrow(line,in,TRUE)) goto end;
      cp = line;
      dim = next_int(&cp);
      printf("dim = %d %s",dim,line);
    }
    else if(strstr(line,"VERTICES")) {
      printf("verts: %s",line);

      if(Getrow(line,in,TRUE)) goto end;
      cp = line;
      noknots = next_int(&cp);      

      printf("noknots = %d %s",noknots,line);

      for(i=1; i <= noknots; i++) {
	GETLINE;
#if 0
	printf("i=%d line=%s",i,line);
#endif
	if(allocated) {
	  cp = line;
#if 0
	  printf("cp = %s",cp);
#endif

	  data->x[i] = next_real(&cp);
	  data->y[i] = next_real(&cp);
	  if(dim > 2) data->z[i] = next_real(&cp);
	}
      }
    }
    else if(strstr(line,"TRIANGLES")) {
      if(Getrow(line,in,TRUE)) goto end;
      cp = line;
      noelements = next_int(&cp);      

      printf("noelements = %d %s",noelements,line);

      elementtype = 303;
      if(maxnodes < 3) maxnodes = 3;

      for(i=1; i <= noelements; i++) {
	GETLINE;
	if(allocated) {
	  cp = line;
	  data->elementtypes[i] = elementtype;
	  for(j=0;j<3;j++)
	    data->topology[i][j] = next_int(&cp);
	  data->material[i] = next_int(&cp);
	}
      }
    }
#if 0
    else printf("unknown command: %s",line);
#endif
  }    

end:
  fclose(in);

printf("ALLOCATED=%d\n",allocated);

  if(!allocated) {
    allocated = TRUE;
    goto allocate;
  }

  printf("Successfully read the mesh from the Medit input file.\n");

  return(0);
}




int LoadGidInput(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info)
/* Load the grid from GID mesh format */
{
  int noknots,noelements,maxnodes,foundsame;
  int mode,allocated,nosides,sideelemtype;
  int boundarytype,side,parent,elemsides,materialtype=0;
  int dim=0, elemnodes=0, elembasis=0, elemtype=0, bulkdone, usedmax=0,hits;
  int minbulk,maxbulk,minbound,maxbound,label,debug;
  int *usedno=NULL, **usedelem=NULL;  
  char filename[MAXFILESIZE],line[MAXLINESIZE],*cp;
  int i,j,k,n,ind,inds[MAXNODESD2],sideind[MAXNODESD1];
  FILE *in;
  Real x,y,z;

  debug = FALSE;

  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"msh");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadGidInput: opening of the GID-file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }

  printf("Reading mesh from GID file %s.\n",filename);
  InitializeKnots(data);

  allocated = FALSE;

  /* Because the file format doesn't provide the number of elements
     or nodes the results are read twice but registered only in the
     second time. */

  minbulk = minbound = 1000;
  maxbulk = maxbound = 0;

omstart:

  mode = 0;
  maxnodes = 0;
  noknots = 0;
  noelements = 0;
  boundarytype = 0;
  nosides = 0;
  bulkdone = FALSE;

  for(;;) {

    if(Getrow(line,in,FALSE)) goto end;

    if(strstr(line,"MESH")) {
      if(debug) printf("MESH\n");

      if(strstr(line,"dimension 3")) 
	dim = 3;
      else if(strstr(line,"dimension 2")) 
	dim = 2;
      else printf("Unknown dimension\n");
      
      if(strstr(line,"ElemType Tetrahedra"))
	elembasis = 500;
      else if(strstr(line,"ElemType Triangle"))
	elembasis = 300;
      else if(strstr(line,"ElemType Linear"))
	elembasis = 200;
      else printf("Unknown elementtype: %s\n",line);

      if(strstr(line,"Nnode 4"))
	elemnodes = 4;
      else if(strstr(line,"Nnode 3"))
	elemnodes = 3;
      else if(strstr(line,"Nnode 2"))
	elemnodes = 2;
      else printf("Unknown elementnode: %s\n",line);

      if(elemnodes > maxnodes) maxnodes = elemnodes;      
      elemtype = elembasis + elemnodes;
      mode = 0;
      
      if(debug) printf("dim=%d elemtype=%d\n",dim,elemtype);
    }
    else if(strstr(line,"Coordinates")) {
      if(debug) printf("Start coords\n");
      mode = 1;
    }
    else if(strstr(line,"end coordinates")) {
      if(debug) printf("End coords\n");
      mode = 0;
    }
    else if(strstr(line,"Elements")) {
      if(bulkdone) { 
	if(debug) printf("Start boundary elems\n");
	mode = 3;
	boundarytype++;
      }
      else {
	if(debug) printf("Start bulk elems\n");	
	mode = 2;
      }
    }
    else if(strstr(line,"end elements")) {
      if(debug) printf("End elems\n");

      mode = 0;
      if(!bulkdone && allocated) {
	usedno = Ivector(1,data->noknots);

	for(j=1;j<=data->noknots;j++)
	  usedno[j] = 0;

	for(j=1;j<=data->noelements;j++) {
	  n = data->elementtypes[j] % 100; 
	  for(i=0;i<n;i++) {
	    ind = data->topology[j][i];
	    usedno[data->topology[j][i]] += 1;
	  }
	}

	usedmax = 0;
	for(i=1;i<=data->noknots;i++)
	  if(usedno[i] > usedmax) usedmax = usedno[i];

	for(j=1;j<=data->noknots;j++)
	  usedno[j] = 0;

	usedelem = Imatrix(1,data->noknots,1,usedmax);
	for(j=1;j<=data->noknots;j++)
	  for(i=1;i<=usedmax;i++)
	    usedelem[j][i] = 0;

	for(j=1;j<=data->noelements;j++) {
	  n = data->elementtypes[j] % 100;
	  for(i=0;i<n;i++) {
	    ind = data->topology[j][i];
	    usedno[ind] += 1;
	    k = usedno[ind];
	    usedelem[ind][k] = j;
	  }
	}
      }
      bulkdone = TRUE;
    }
    else if(!mode) {
      if(debug) printf("mode: %d %s\n",mode,line);
    }
    else if(mode) {  
      switch (mode) {
	
      case 1:
	cp = line;
	ind = next_int(&cp);
	if(ind > noknots) noknots = ind;

	x = next_real(&cp);
	y = next_real(&cp);
	if(dim == 3) z = next_real(&cp);
	
	if(allocated) {
	  data->x[ind] = x;
	  data->y[ind] = y;
	  if(dim == 3) data->z[ind] = z;
	}
	break;

      case 2:
	cp = line;
	ind = next_int(&cp);
	if(ind > noelements) noelements = ind;
	
	for(i=0;i<elemnodes;i++) {
	  k = next_int(&cp);
	  if(allocated) {
	    data->topology[ind][i] = k;
	    data->elementtypes[ind] = elemtype;
	    data->material[ind] = 1;
	  }
	}
	label = next_int(&cp);
	
	if(allocated) {
	  if(label) { 
	    materialtype = label-minbulk+1;
	  }
	  else if(maxbound) {
	    materialtype = maxbulk-minbulk + 2;
	  }
	  else { 
	    materialtype = 1;
	  }
	  data->material[ind] = materialtype;
	}
	else {
	  if(label > maxbulk) maxbulk = label;
	  if(label < minbulk) minbulk = label;	  
	}

	break;
	
      case 3:
	cp = line;
	ind = next_int(&cp);
	nosides++;

	for(i=0;i<elemnodes;i++) {
	  k = next_int(&cp);
	  inds[i] = k;
	}
	label = next_int(&cp);
	
	if(!allocated) {
	  if(label) {
	    if(label > maxbound) maxbound = label;
	    if(label < minbound) minbound = label;
	  }
	}

	if(allocated) {
	  
	  if(label) { 
	    boundarytype = label-minbound+1;
	  }
	  else if(maxbound) {
	    boundarytype = maxbound-minbound + 2;
	  }
	  else { 
	    boundarytype = 1;
	  }

	  foundsame = FALSE;
	  for(i=1;i<=usedno[inds[0]];i++) {
	    parent = usedelem[inds[0]][i];

	    elemsides = data->elementtypes[parent] % 100;
	    
	    for(side=0;side<elemsides;side++) {

	      GetElementSide(parent,side,1,data,&sideind[0],&sideelemtype);

	      if(elemnodes != sideelemtype%100) printf("LoadGidMesh: bug?\n");
	      
	      hits = 0;
	      for(j=0;j<elemnodes;j++) 
		for(k=0;k<elemnodes;k++) 
		  if(sideind[j] == inds[k]) hits++;

	      if(hits < elemnodes) continue;
	      
	      if(!foundsame) {
		foundsame++;
		bound->parent[nosides] = parent;
		bound->side[nosides] = side;
		bound->parent2[nosides] = 0;
		bound->side2[nosides] = 0;
		bound->types[nosides] = boundarytype;	  
	      }
	      else if(foundsame == 1) {
		if(parent == bound->parent[nosides]) continue;
		foundsame++;
		bound->parent2[nosides] = parent;
		bound->side2[nosides] = side;		  
	      }
	      else if(foundsame > 1) {
		printf("Boundary %d has more than 2 parents\n",nosides);
	      }
	    }
	  }
	  if(!foundsame) {
	    printf("Did not find parent for side %d\n",nosides);
	    nosides--;
	  }
	  else {
	    if(debug) printf("Parent of side %d is %d\n",nosides,bound->parent[nosides]);
	  }
	}
      }
    }
  }

end:

  if(!allocated) {
    rewind(in);
    data->noknots = noknots;
    data->noelements = noelements;
    data->maxnodes = maxnodes;
    data->dim = dim;
    
    if(info) {
      printf("Allocating for %d knots and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);
      printf("Initial material indexes are at interval %d to %d.\n",minbulk,maxbulk);
    }  
    AllocateKnots(data);

    if(info) {
      printf("Allocating %d boundary elements\n",nosides);
      printf("Initial boundary indexes are at interval %d to %d.\n",minbound,maxbound);
    }
    AllocateBoundary(bound,nosides);

    bound->nosides = nosides;
    bound->created = TRUE;
    
    nosides = 0;
    bulkdone = FALSE;
    boundarytype = 0;

    allocated = TRUE;
    goto omstart;
  }

  bound->nosides = nosides;
  free_Ivector(usedno,1,data->noknots);
  free_Imatrix(usedelem,1,data->noknots,1,usedmax);

  if(info) printf("The mesh was loaded from file %s.\n",filename);
  return(0);
}



static void ReorderComsolNodes(int elementtype,int *topo) 
{
  int i,tmptopo[MAXNODESD2];
  int order404[]={1,2,4,3};
  int order808[]={1,2,4,3,5,6,8,7};
  int order605[]={1,2,4,3,5};

   
  switch (elementtype) {
 
  case 404:
    for(i=0;i<elementtype%100;i++) 
      tmptopo[i] = topo[i];
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order404[i]-1];
    break;

  case 808:
    for(i=0;i<elementtype%100;i++) 
      tmptopo[i] = topo[i];
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order808[i]-1];
    break;

  case 605:
    for(i=0;i<elementtype%100;i++) 
      tmptopo[i] = topo[i];
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order605[i]-1];
    break;


  default:
    break;

  }  
}



int LoadComsolMesh(struct FemType *data,char *prefix,int info)
/* Load the grid in Comsol Multiphysics mesh format */
{
  int noknots,noelements,maxnodes,material;
  int allocated,dim=0, elemnodes=0, elembasis=0, elemtype;
  int debug,offset,domains,mindom,minbc,elemdim=0;
  char filename[MAXFILESIZE],line[MAXLINESIZE],*cp;
  int i,j,k;
  FILE *in;

  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"mphtxt");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadComsolMesh: opening of the Comsol mesh file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }

  printf("Reading mesh from Comsol mesh file %s.\n",filename);
  InitializeKnots(data);

  debug = FALSE;
  allocated = FALSE;

  mindom = 1000;
  minbc = 1000;
  offset = 1;

omstart:

  maxnodes = 0;
  noknots = 0;
  noelements = 0;
  material = 0;
  domains = 0;


  for(;;) {

    if(Comsolrow(line,in)) goto end;

    if(strstr(line,"# sdim")) {
      cp = line;
      dim = next_int(&cp);
      if(debug) printf("dim=%d\n",dim);
    }

    else if(strstr(line,"# number of mesh points")) {
      cp = line;
      noknots = next_int(&cp);
      if(debug) printf("noknots=%d\n",noknots);
    }

    else if(strstr(line,"# lowest mesh point index")) {
      cp = line;
      offset = 1 - next_int(&cp);
      if(debug) printf("offset=%d\n",offset);
    }

    else if(strstr(line,"# type name")) {
      if(strstr(line,"vtx")) elembasis = 100;
      else if(strstr(line,"edg")) elembasis = 200;
      else if(strstr(line,"tri")) elembasis = 300;
      else if(strstr(line,"quad")) elembasis = 400;
      else if(strstr(line,"tet")) elembasis = 500;
      else if(strstr(line,"pyr")) elembasis = 600;
      else if(strstr(line,"prism")) elembasis = 700;
      else if(strstr(line,"hex")) elembasis = 800;
      else {
	printf("unknown element type = %s",line);
	bigerror("Cannot continue!\n");
      }
    }

    else if(strstr(line,"# number of nodes per element")) {
      cp = line;
      elemnodes = next_int(&cp);
      if(elemnodes > maxnodes) maxnodes = elemnodes;      
      if(debug) printf("elemnodes=%d\n",elemnodes);           
    }

    else if(strstr(line,"# Mesh point coordinates")) {
      printf("Loading %d coordinates\n",noknots);

      for(i=1;i<=noknots;i++) {
	Comsolrow(line,in);	

	if(allocated) {
	  cp = line;
	  data->x[i] = next_real(&cp);
	  data->y[i] = next_real(&cp);
	  if(dim == 3) data->z[i] = next_real(&cp);
	}
      }
    }

    else if(strstr(line,"# number of elements")) {

      cp = line;
      k = next_int(&cp);

      Comsolrow(line,in);	            
      elemtype = elemnodes + elembasis;
      elemdim = GetElementDimension(elemtype);

      if(debug) printf("Loading %d elements of type %d\n",k,elemtype);
      domains = noelements;

      for(i=1;i<=k;i++) {
	Comsolrow(line,in);	

	if(dim == 3 && elembasis < 300) continue;
	if(dim == 2 && elembasis < 200) continue;

	noelements = noelements + 1;
	if(allocated) {
	  data->elementtypes[noelements] = elemtype;
	  data->material[noelements] = 1;	  
	  cp = line;
	  for(j=0;j<elemnodes;j++)
	    data->topology[noelements][j] = next_int(&cp) + offset;

	  if(1) ReorderComsolNodes(elemtype,data->topology[noelements]);

	}
      }
    }

    else if(strstr(line,"# number of geometric entity indices") || 
	    strstr(line,"# number of domains")) {

      cp = line;
      k = next_int(&cp);

      Comsolrow(line,in);	            
      if(debug) printf("Loading %d domains for the elements\n",k);

      for(i=1;i<=k;i++) {
	Comsolrow(line,in);	

	if(dim == 3 && elembasis < 300) continue;
	if(dim == 2 && elembasis < 200) continue;

	domains = domains + 1;
	cp = line;
	material = next_int(&cp);

	if(allocated) {
	  if(elemdim < dim) 
	    material = material - minbc + 1;
	  else 
	    material = material - mindom + 1;
	  data->material[domains] = material;	  
	}
	else {
	  if(elemdim < dim) {
	    if(minbc > material) minbc = material;
	  }
	  else { 
	    if(mindom > material) mindom = material;	  
	  }
	}

      }
    }

    else if(strstr(line,"#")) {
      if(debug) printf("Unused command:  %s",line);
    }

  }

end:

  if(!allocated) {

    if(noknots == 0 || noelements == 0 || maxnodes == 0) {
       printf("Invalid mesh consits of %d knots and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);     
       fclose(in);
       return(2);
    }

    rewind(in);
    data->noknots = noknots;
    data->noelements = noelements;
    data->maxnodes = maxnodes;
    data->dim = dim;
    
    if(info) {
      printf("Allocating for %d knots and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);
    }  
    AllocateKnots(data);
    allocated = TRUE;

    goto omstart;    
  }
  fclose(in);

  if(info) printf("The Comsol mesh was loaded from file %s.\n\n",filename);
  return(0);
}



static int GmshToElmerType(int gmshtype)
{
  int elmertype = 0;

  switch (gmshtype) {
      
  case 1:        
    elmertype = 202;
    break;
  case 2:        
    elmertype = 303;
    break;
  case 3:        
    elmertype = 404;
    break;
  case 4:        
    elmertype = 504;
    break;
  case 5:        
    elmertype = 808;
    break;
  case 6:        
    elmertype = 706;
    break;
  case 7:        
    elmertype = 605;
    break;
  case 8:        
    elmertype = 203;
    break;
  case 9:        
    elmertype = 306;
    break;
  case 10:        
    elmertype = 409;
    break;
  case 11:        
    elmertype = 510;
    break;
  case 12:        
    elmertype = 827;
    break;
  case 15:        
    elmertype = 101;
    break;
  case 16:
    elmertype = 408;
    break;
  case 17:
    elmertype = 820;
    break;
  case 18:
    elmertype = 715;
    break;
  case 19:
    elmertype = 613;
    break;
  case 21:
    elmertype = 310;
    break;

  /* These are supported by Gmsh but not by ElmerSolver */
  case 13:        
    elmertype = 718;
    break;
  case 14:        
    elmertype = 614;
    break;
  case 20:        
    elmertype = 309;
    break;
  case 22:        
    elmertype = 312;
    break;
  case 24:
    elmertype = 315;
    break;
  case 25:
    elmertype = 320;
    break;
  case 26:
    elmertype = 204;
    break;
  case 27:
    elmertype = 205;
    break;
  case 28:
    elmertype = 206;
    break;
  case 29:
    elmertype = 520;
    break;

  default:
    printf("Gmsh element %d does not have an Elmer counterpart!\n",gmshtype);
  }

  return(elmertype);
}


static void GmshToElmerIndx(int elemtype,int *topology)
{
  int i=0,nodes=0,oldtopology[MAXNODESD2];
  int reorder, *porder;

  int order510[]={0,1,2,3,4,5,6,7,9,8};
  int order614[]={0,1,2,3,4,5,8,10,6,7,9,11,12,13};
  int order718[]={0,1,2,3,4,5,6,9,7,8,10,11,12,14,13,15,17,16};
  int order820[]={0,1,2,3,4,5,6,7,8,11,13,9,10,12,14,15,16,18,19,17};
  int order827[]={0,1,2,3,4,5,6,7,8,11,13,9,10,12,14,15,16,18,19,17,21,23,24,22,20,25,26};
  /*             {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26}; */
  

  reorder = FALSE;

  switch (elemtype) {
      
  case 510:        
    reorder = TRUE;
    porder = &order510[0];
    break;

  case 613:        
  case 614:        
    reorder = TRUE;
    porder = &order614[0];
    break;

  case 715:        
  case 718:
    reorder = TRUE;
    porder = &order718[0];
    break;

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


static int LoadGmshInput1(struct FemType *data,struct BoundaryType *bound,
			  char *filename,int info)
{
  int noknots = 0,noelements = 0,maxnodes,dim;
  int elemind[MAXNODESD2],elementtype;
  int i,j,k,allocated,*revindx=NULL,maxindx;
  int elemno, gmshtype, regphys, regelem, elemnodes,maxelemtype;
  FILE *in;
  char *cp,line[MAXLINESIZE];


  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadGmshInput: The opening of the mesh file %s failed!\n",filename);
    return(1);
  }
  if(info) printf("Loading mesh in Gmsh format 1.0 from file %s\n",filename);

  allocated = FALSE;
  dim = data->dim;
  maxnodes = 0;
  maxindx = 0;
  maxelemtype = 0;

allocate:

  if(allocated) {
    InitializeKnots(data);
    data->dim = dim;
    data->maxnodes = maxnodes;
    data->noelements = noelements;
    data->noknots = noknots;

    if(info) printf("Allocating for %d knots and %d elements.\n",noknots,noelements);
    AllocateKnots(data);

    if(maxindx > noknots) {
      revindx = Ivector(1,maxindx);
      for(i=1;i<=maxindx;i++) revindx[i] = 0;
    }
    in = fopen(filename,"r");
  }


  for(;;) {
    if(Getrow(line,in,TRUE)) goto end;
    if(strstr(line,"END")) goto end;

    if(strstr(line,"$NOD")) {
      
      GETLINE;
      cp = line;
      noknots = next_int(&cp);

      for(i=1; i <= noknots; i++) {
	GETLINE;
	cp = line;

	j = next_int(&cp);
	if(allocated) {
	  if(maxindx > noknots) revindx[j] = i;
	  data->x[i] = next_real(&cp);
	  data->y[i] = next_real(&cp);
	  if(dim > 2) data->z[i] = next_real(&cp);
	}
	else {
	  maxindx = MAX(j,maxindx);
	}
      }
      GETLINE;
      if(!strstr(line,"$ENDNOD")) {
	printf("NOD section should end to string ENDNOD\n");
	printf("%s\n",line);
      }
    }
    
    if(strstr(line,"$ELM")) {
      GETLINE;
      cp = line;
      noelements = next_int(&cp);

      for(i=1; i <= noelements; i++) {
	GETLINE;
	
	cp = line;
	elemno = next_int(&cp);
	gmshtype = next_int(&cp);
	regphys = next_int(&cp);
	regelem = next_int(&cp);
	elemnodes = next_int(&cp);

	if(allocated) {
	  elementtype = GmshToElmerType(gmshtype);
	  maxelemtype = MAX(maxelemtype,elementtype);
	  data->elementtypes[i] = elementtype;
	  data->material[i] = regphys;

	  if(elemnodes != elementtype%100) {
	    printf("Conflict in elementtypes %d and number of nodes %d!\n",
		   elementtype,elemnodes);
	  }	  
	  for(j=0;j<elemnodes;j++)
	    elemind[j] = next_int(&cp);
	  
	  GmshToElmerIndx(elementtype,elemind);	  

	  for(j=0;j<elemnodes;j++)
	    data->topology[i][j] = elemind[j];
	}
	else {
	  maxnodes = MAX(elemnodes,maxnodes);
	}

      }
      GETLINE;
      if(!strstr(line,"$ENDELM")) 
	printf("ELM section should end to string ENDELM\n");
    }

  }

 end:

  fclose(in);

  if(!allocated) {
    allocated = TRUE;
    goto allocate;
  }


  if(maxindx > noknots) {
    printf("Renumbering the Gmsh nodes from %d to %d\n",maxindx,noknots);

    for(i=1; i <= noelements; i++) {
      elementtype = data->elementtypes[i];
      elemnodes = elementtype % 100; 

      for(j=0;j<elemnodes;j++) {
	k = data->topology[i][j];
	if(k <= 0 || k > maxindx) 
	  printf("index out of bounds %d\n",k);
	else if(revindx[k] <= 0) 
	  printf("unkonwn node %d %d in element %d\n",k,revindx[k],i);
	else 
	  data->topology[i][j] = revindx[k];
      }      
    }
    free_Ivector(revindx,1,maxindx);
  }

  ElementsToBoundaryConditions(data,bound,FALSE,info);

  printf("Successfully read the mesh from the Gmsh input file.\n");

  return(0);
}


static int LoadGmshInput2(struct FemType *data,struct BoundaryType *bound,
			  char *filename,int info)
{
  int noknots = 0,noelements = 0,nophysical = 0,maxnodes,dim,notags;
  int elemind[MAXNODESD2],elementtype;
  int i,j,k,allocated,*revindx=NULL,maxindx;
  int elemno, gmshtype, tagphys=0, taggeom=0, tagpart, elemnodes,maxelemtype;
  int usetaggeom,tagmat,verno;
  int physvolexist, physsurfexist;
  FILE *in;
  const char manifoldname[4][10] = {"point", "line", "surface", "volume"};
  char *cp,line[MAXLINESIZE];

  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadGmshInput2: The opening of the mesh file %s failed!\n",filename);
    return(1);
  }
  if(info) printf("Loading mesh in Gmsh format 2.0 from file %s\n",filename);

  allocated = FALSE;
  dim = data->dim;
  maxnodes = 0;
  maxindx = 0;
  maxelemtype = 0;
  usetaggeom = FALSE;
  physvolexist = FALSE;
  physsurfexist = FALSE;

omstart:

  for(;;) {
    if(Getrow(line,in,FALSE)) goto end;
    if(strstr(line,"$End")) continue;

    if(strstr(line,"$MeshFormat")) {
      GETLINE;
      cp = line;
      verno = next_int(&cp);

      if(verno != 2) {
	printf("Version number is not compatible with the parser: %d\n",verno);
      }

      GETLINE;
      if(!strstr(line,"$EndMeshFormat")) {
	printf("$MeshFormat section should end to string $EndMeshFormat:\n%s\n",line);
      }      
    }
      
    else if(strstr(line,"$Nodes")) {
      GETLINE;
      cp = line;
      noknots = next_int(&cp);

      for(i=1; i <= noknots; i++) {
	GETLINE;
	cp = line;

	j = next_int(&cp);
	if(allocated) {
	  if(maxindx > noknots) revindx[j] = i;
	  data->x[i] = next_real(&cp);
	  data->y[i] = next_real(&cp);
	  if(dim > 2) data->z[i] = next_real(&cp);
	}
	else {
	  maxindx = MAX(j,maxindx);
	}
      }
      GETLINE;
      if(!strstr(line,"$EndNodes")) {
	printf("$Nodes section should end to string $EndNodes:\n%s\n",line);
      }           
    }
    
    else if(strstr(line,"$Elements")) {
      GETLINE;
      cp = line;
      noelements = next_int(&cp);

      for(i=1; i <= noelements; i++) {
	GETLINE;
	
	cp = line;
	elemno = next_int(&cp);
	gmshtype = next_int(&cp);
	elementtype = GmshToElmerType(gmshtype);

	if(allocated) {
	  elemnodes = elementtype % 100;
	  data->elementtypes[i] = elementtype;

	  /* Point does not seem to have physical properties */
	  notags = next_int(&cp);
	  if(notags > 0) tagphys = next_int(&cp);
	  if(notags > 1) taggeom = next_int(&cp);
	  if(notags > 2) tagpart = next_int(&cp);
	  for(j=4;j<=notags;j++)
	    next_int(&cp);

	  if(tagphys) {
	    tagmat = tagphys;
	  }
	  else {
	    tagmat = taggeom;
	    usetaggeom = TRUE;
	  }

	  data->material[i] = tagmat;
	  for(j=0;j<elemnodes;j++)
	    elemind[j] = next_int(&cp);

	  GmshToElmerIndx(elementtype,elemind);	  

	  for(j=0;j<elemnodes;j++)
	    data->topology[i][j] = elemind[j];
	}
	else {
	  maxelemtype = MAX(maxelemtype,elementtype);
	}
	
      }
      GETLINE;
      if(!strstr(line,"$EndElements")) {
	printf("$Elements section should end to string $EndElements:\n%s\n",line);
      }   
    }
    else if(strstr(line,"$PhysicalNames")) {
      GETLINE;
      cp = line;
      nophysical = next_int(&cp);
      for(i=0;i<nophysical;i++) {
	GETLINE;
        if(allocated) {
            cp = line;
            gmshtype = next_int(&cp);
            tagphys = next_int(&cp);
            if(gmshtype == dim-1) {
                physsurfexist = TRUE;
                if(tagphys < MAXBCS) sscanf(cp," \"%[^\"]\"",data->boundaryname[tagphys]);
                else printf("Index %d too high: ignoring physical %s %s",tagphys,manifoldname[dim-1],cp+1);
            }
            else if(gmshtype == dim) {
                physvolexist = TRUE;
                if(tagphys < MAXBODIES) sscanf(cp," \"%[^\"]\"",data->bodyname[tagphys]);
                else printf("Index %d too high: ignoring physical %s %s",tagphys,manifoldname[dim],cp+1);
            }
            else printf("Physical groups of dimension %d not supported in %d-dimensional mesh: "
                        "ignoring group %d %s",gmshtype,dim,tagphys,cp+1);
        }
      }

      GETLINE;
      if(!strstr(line,"$EndPhysicalNames")) {
	printf("$PhysicalNames section should end to string $EndPhysicalNames:\n%s\n",line);
      }   
    }
    else {
      if(!allocated) printf("Untreated command: %s",line);
    }

  }

 end:


  if(!allocated) {
    maxnodes = maxelemtype % 100;
    InitializeKnots(data);
    data->dim = dim;
    data->maxnodes = maxnodes;
    data->noelements = noelements;
    data->noknots = noknots;

    if(info) printf("Allocating for %d knots and %d elements.\n",noknots,noelements);
    AllocateKnots(data);

    if(maxindx > noknots) {
      revindx = Ivector(1,maxindx);
      for(i=1;i<=maxindx;i++) revindx[i] = 0;
    }
    rewind(in);
    allocated = TRUE;
    goto omstart;
  }

  if(maxindx > noknots) {
    printf("Renumbering the Gmsh nodes from %d to %d\n",maxindx,noknots);

    for(i=1; i <= noelements; i++) {
      elementtype = data->elementtypes[i];
      elemnodes = elementtype % 100; 

      for(j=0;j<elemnodes;j++) {
	k = data->topology[i][j];
	if(k <= 0 || k > maxindx) 
	  printf("index out of bounds %d\n",k);
	else if(revindx[k] <= 0) 
	  printf("unkonwn node %d %d in element %d\n",k,revindx[k],i);
	else 
	  data->topology[i][j] = revindx[k];
      }      
    }
    free_Ivector(revindx,1,maxindx);
  }

  ElementsToBoundaryConditions(data,bound,FALSE,info);

  /* The geometric entities are rather randomly numbered */
  if( usetaggeom ) {
    RenumberBoundaryTypes(data,bound,TRUE,0,info);
    RenumberMaterialTypes(data,bound,info);
  }
  data->bodynamesexist = physvolexist;
  data->boundarynamesexist = physsurfexist;

  if(info) printf("Successfully read the mesh from the Gmsh input file.\n");

  return(0);
}


static int LoadGmshInput4(struct FemType *data,struct BoundaryType *bound,
			  char *filename,int info)
{
  int noknots = 0,noelements = 0,nophysical = 0,maxnodes,dim,notags;
  int elemind[MAXNODESD2],elementtype;
  int i,j,k,l,allocated,*revindx=NULL,maxindx;
  int elemno, gmshtype, tagphys=0, taggeom=0, tagpart, elemnodes,maxelemtype;
  int tagmat,verno;
  int physvolexist, physsurfexist,**tagmap,tagsize;
  FILE *in;
  const char manifoldname[4][10] = {"point", "line", "surface", "volume"};
  char *cp,line[MAXLINESIZE],longline[LONGLINESIZE];

  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadGmshInput4: The opening of the mesh file %s failed!\n",filename);
    return(1);
  }
  if(info) printf("Loading mesh in Gmsh format 4.0 from file %s\n",filename);

  allocated = FALSE;
  dim = data->dim;
  maxnodes = 0;
  maxindx = 0;
  maxelemtype = 0;
  physvolexist = FALSE;
  physsurfexist = FALSE;

omstart:

  for(;;) {
    if(Getrow(line,in,FALSE)) goto end;
    if(strstr(line,"$End")) continue;
 
    if(strstr(line,"$MeshFormat")) {
      GETLINE;
      cp = line;
      verno = next_int(&cp);

      if(verno != 4) {
	printf("Version number is not compatible with the parser: %d\n",verno);
      }

      GETLINE;
      if(!strstr(line,"$EndMeshFormat")) {
	printf("$MeshFormat section should end to string $EndMeshFormat:\n%s\n",line);
      }      
    }
      
    else if(strstr(line,"$Nodes")) {
      int numEntityBlocks,tagEntity,dimEntity,parEntity,numNodes,ind;

      GETLINE;
      cp = line;

      numEntityBlocks = next_int(&cp);
      noknots = next_int(&cp);

      if(allocated && info) printf("Reading %d nodes in %d blocks.\n",noknots,numEntityBlocks);
      
      k = 0;
      
      for(j=1; j <= numEntityBlocks; j++) {
	GETLINE;
	cp = line;

	tagEntity = next_int(&cp);
	dimEntity = next_int(&cp);
	parEntity = next_int(&cp);
	numNodes = next_int(&cp);
	       	
	for(i=1; i <= numNodes; i++) {
	  GETLINE;
	  cp = line;
	  k += 1;
	  
	  ind = next_int(&cp);
	  if(allocated) {
	    if(maxindx > noknots) revindx[ind] = k;
	    data->x[k] = next_real(&cp);
	    data->y[k] = next_real(&cp);
	    if(dim > 2) data->z[k] = next_real(&cp);
	  }
	  else {
	    maxindx = MAX(ind,maxindx);
	  }
	}
      }
      GETLINE;

      if(!strstr(line,"$EndNodes")) {
	printf("$Nodes section should end to string $EndNodes:\n%s\n",line);
      }           
    }

    else if(strstr(line,"$Entities")) {
      int numPoints, numCurves, numSurfaces, numVolumes, numEnt;
      int tag,tagdim,nophys,phystag,maxtag[4];
      int nobound, idum;
      Real rdum;
      
      GETLINE;
      cp = line;
      numPoints = next_int(&cp);
      numCurves = next_int(&cp);
      numSurfaces = next_int(&cp);
      numVolumes = next_int(&cp);

      
      if(allocated) {
	tagsize = 0;
	for(tagdim=0;tagdim<=3;tagdim++)
	  tagsize = MAX( tagsize, maxtag[tagdim]);
	if( tagsize > 0 ) {
	  tagmap = Imatrix(0,3,1,tagsize);
	  for(i=0;i<=3;i++)
	    for(j=1;j<=tagsize;j++)
	      tagmap[i][j] = 0;
	}
      }
      
      for(tagdim=0;tagdim<=3;tagdim++) {	

	
	if( tagdim == 0 ) 
	  numEnt = numPoints;
	else if( tagdim == 1 )
	  numEnt = numCurves;
	else if( tagdim == 2 )
	  numEnt = numSurfaces;
	else if( tagdim == 3 )
	  numEnt = numVolumes;
	
	if(!allocated)
	  maxtag[tagdim] = 0;
	else if( maxtag[tagdim] > 0 )
	  printf("Maximum original tag for %d %dDIM entities is %d\n",numEnt,tagdim,maxtag[tagdim]);

	if(numEnt > 0 && !allocated) printf("Reading %d entities in %dD\n",numEnt,tagdim);

	
	for(i=1; i <= numEnt; i++) {
	  GETLONGLINE;

	  // if( i==1 ) printf("1st line of dim %d with %d entries: %s\n",tagdim,numEnt,line);
	  
	  if( tagdim == 0 ) continue;
	  
	  cp = longline;
	  tag = next_int(&cp);
	 	  
	  if(!allocated)
	    maxtag[tagdim] = MAX( maxtag[tagdim], tag );
	  
	  for(j=1;j<=6;j++) rdum = next_real(&cp);
	  nophys = next_int(&cp);

	  if( nophys > 0 ) phystag = next_int(&cp);
	  
	  if(allocated) tagmap[tagdim][tag] = phystag;


	  // The lines may be too long. So fill the string buffer until we get a newline. 
	  j = k = 0;
	  for(;;) { 
	    for(l=0;l<LONGLINESIZE;l++) {
	      if( longline[l] == '\n') {
		j = l;
	    	break;	    
	      }
	    }
	    if(j) break;
	    k += LONGLINESIZE;
	    GETLONGLINE;
	  }	   	    	      
	  if( k > 0 && !allocated) printf("Entity line %d has length %d.\n",i,k+j);
	  
	  //for(j=2;j<=nophys;j++)
	  //  idum = next_int(&cp);

	  //// if( tagdim == 0 ) continue;

	  //nobound = next_int(&cp);
	  // for(j=1;j<=nobound;j++)
	  //  idum = next_int(&cp);	  
	}
      }
      
      GETLONGLINE;
      if(!strstr(longline,"$EndEntities")) {
	printf("$Entities section should end to string $EndEntities:\n%s\n",longline);
      }           
    }

    else if(strstr(line,"$Elements")) {
      int numEntityBlocks, numElements, tagEntity, dimEntity, typeEle, NumElements;
      
      GETLINE;
      cp = line;

      k = 0;
      numEntityBlocks = next_int(&cp);
      noelements = next_int(&cp);

      if(allocated) printf("Reading %d elements in %d blocks.\n",noelements,numEntityBlocks);

      
      for(j=1; j<= numEntityBlocks; j++ ) {
	
	GETLINE;	
	cp = line;

	tagEntity = next_int(&cp);
	dimEntity = next_int(&cp);

	typeEle = next_int(&cp);
	numElements = next_int(&cp);
	
	elementtype = GmshToElmerType(typeEle);
	elemnodes = elementtype % 100;
	maxelemtype = MAX(maxelemtype,elementtype);
	
	if( allocated && tagsize > 0 ) {
	  printf("Reading %d elements with tag %d of type %d\n", numElements, tagEntity, elementtype);
	  if( tagsize > 0 ) {
	    if( tagmap[dimEntity][tagEntity] ) {
	      printf("Mapping mesh tag %d to physical tag %d in %dDIM\n",tagEntity,tagmap[dimEntity][tagEntity],dimEntity);	    
	      tagEntity = tagmap[dimEntity][tagEntity];
	    }
	    else {
	      printf("Mesh tag %d is not associated to any physical tag!\n",tagEntity);
	    }
	  }
	}
			     
	for(i=1; i <= numElements; i++) {
	  GETLINE;	
	  cp = line;

	  k += 1;
	  	  
	  elemno = next_int(&cp);
	  
	  if(allocated) {
	    data->elementtypes[k] = elementtype;
	    data->material[k] = tagEntity;
	    for(l=0;l<elemnodes;l++)
	      elemind[l] = next_int(&cp);

	    GmshToElmerIndx(elementtype,elemind);	  

	    for(l=0;l<elemnodes;l++)
	      data->topology[k][l] = elemind[l];
	  }	
	}
      }

      GETLINE;
      if(!strstr(line,"$EndElements")) {
	printf("$Elements section should end to string $EndElements:\n%s\n",line);
      }   
    }

    else if(strstr(line,"$PhysicalNames")) {
      GETLINE;
      cp = line;
      nophysical = next_int(&cp);
      for(i=0;i<nophysical;i++) {
	GETLINE;
        if(allocated) {
	  cp = line;
	  gmshtype = next_int(&cp);
	  tagphys = next_int(&cp);
	  if(gmshtype == dim-1) {
	    physsurfexist = TRUE;
	    if(tagphys < MAXBCS) {
	      sscanf(cp," \"%[^\"]\"",data->boundaryname[tagphys]);
	      printf("Boundary name for physical group %d is: %s\n",tagphys,data->boundaryname[tagphys]);
	    }
	    else
	      printf("Index %d too high: ignoring physical %s %s",tagphys,manifoldname[dim-1],cp+1);
	  }
	  else if(gmshtype == dim) {
	    physvolexist = TRUE;
	    if(tagphys < MAXBODIES) {
	      sscanf(cp," \"%[^\"]\"",data->bodyname[tagphys]);
	      printf("Body name for physical group %d is: %s\n",tagphys,data->bodyname[tagphys]);
	    }
	    else
	      printf("Index %d too high: ignoring physical %s %s",tagphys,manifoldname[dim],cp+1);
	  }
	  else printf("Physical groups of dimension %d not supported in %d-dimensional mesh: "
		      "ignoring group %d %s",gmshtype,dim,tagphys,cp+1);
        }
      }

      GETLINE;
      if(!strstr(line,"$EndPhysicalNames")) {
	printf("$PhysicalNames section should end to string $EndPhysicalNames:\n%s\n",line);
      }   
    }
    else if(strstr(line,"$Periodic")) {
      int numPeriodicLinks;
      if(allocated) printf("Reading periodic links but doing nothing with them!\n");
      
      GETLINE;
      cp = line;
      numPeriodicLinks = next_int(&cp);
      for(i=1; i <= numPeriodicLinks; i++) {
	GETLINE;
      }     
      GETLINE;
      if(!strstr(line,"$EndPeriodic")) {
	printf("$Periodic section should end to string $EndPeriodic:\n%s\n",line);
      }           
    }

    else if(strstr(line,"$PartitionedEntities")) {
      if(allocated) printf("Reading partitioned entities but doing nothing with them!\n");      
      for(;;) {
	GETLINE;
	if(strstr(line,"$EndPartitionedEntities")) break;
      }
    }
    else if(strstr(line,"$NodeData")) {
      if(allocated) printf("Reading node data but doing nothing with them!\n");      
      for(;;) {
	GETLINE;
	if(strstr(line,"$EndNodeData")) break;
      }
    }
    else if(strstr(line,"$ElementData")) {
      if(allocated) printf("Reading element data but doing nothing with them!\n");      
      for(;;) {
	GETLINE;
	if(strstr(line,"$EndElementData")) break;
      }
    }
    else if(strstr(line,"$ElementNodeData")) {
      if(allocated) printf("Reading element node data but doing nothing with them!\n");      
      for(;;) {
	GETLINE;
	if(strstr(line,"$EndElementNodeData")) break;
      }
    }
    else if(strstr(line,"$GhostElements")) {
      if(allocated) printf("Reading ghost elements data but doing nothing with them!\n");      
      for(;;) {
	GETLINE;
	if(strstr(line,"$EndGhostElements")) break;
      }
    }
    else if(strstr(line,"$InterpolationScheme")) {
      if(allocated) printf("Reading interpolation scheme but doing nothing with them!\n");      
      for(;;) {
	GETLINE;
	if(strstr(line,"$EndInterpolationScheme")) break;
      }
    }    
    else {
      if(allocated) printf("Untreated command: %s",line);
    }

  }

 end:


  if(!allocated) {
    if( noelements == 0 ) bigerror("No elements to load in Gmsh file!");
    if( noknots == 0 ) bigerror("No nodes to load in Gmsh file!");

    maxnodes = maxelemtype % 100;
    InitializeKnots(data);
    data->dim = dim;
    data->maxnodes = maxnodes;
    data->noelements = noelements;
    data->noknots = noknots;

    if(info) printf("Allocating for %d knots and %d elements.\n",noknots,noelements);
    AllocateKnots(data);

    if(maxindx > noknots) {
      revindx = Ivector(1,maxindx);
      for(i=1;i<=maxindx;i++) revindx[i] = 0;
    }
    rewind(in);
    allocated = TRUE;
    goto omstart;
  }

  if(maxindx > noknots) {
    printf("Renumbering the Gmsh nodes from %d to %d\n",maxindx,noknots);

    for(i=1; i <= noelements; i++) {
      elementtype = data->elementtypes[i];
      elemnodes = elementtype % 100; 

      for(j=0;j<elemnodes;j++) {
	k = data->topology[i][j];
	if(k <= 0 || k > maxindx) 
	  printf("index out of bounds %d\n",k);
	else if(revindx[k] <= 0) 
	  printf("unkonwn node %d %d in element %d\n",k,revindx[k],i);
	else 
	  data->topology[i][j] = revindx[k];
      }      
    }
    free_Ivector(revindx,1,maxindx);
  }

  ElementsToBoundaryConditions(data,bound,FALSE,info);

  data->bodynamesexist = physvolexist;
  data->boundarynamesexist = physsurfexist;
  
  if( tagsize > 0 ) free_Imatrix(tagmap,0,3,1,tagsize);
  
  if(info) printf("Successfully read the mesh from the Gmsh input file.\n");

  return(0);
}



int LoadGmshInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int info)
{
  FILE *in;
  char line[MAXLINESIZE],filename[MAXFILESIZE];
  int errno;

  sprintf(filename,"%s",prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    sprintf(filename,"%s.msh",prefix);
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadGmshInput: The opening of the mesh file %s failed!\n",filename);
      return(1);
    }
  }

  Getrow(line,in,FALSE);

  if(info) {
    printf("Format chosen using the first line: %s",line);
  }

  if(strstr(line,"$")) {
    int verno;
    char *cp;
    
    Getrow(line,in,FALSE);
    cp = line;    
    verno = next_int(&cp);
    fclose(in);
    
    if( verno == 4 )
      errno = LoadGmshInput4(data,bound,filename,info);
    else      
      errno = LoadGmshInput2(data,bound,filename,info);
    
  } else {
    fclose(in);
    printf("*****************************************************\n");
    printf("The first line did not start with $, assuming Gmsh 1 format\n");
    printf("This version of Gmsh format is no longer supported\n");
    printf("Please use Gsmh 2 or 4 versions for output\n");
    printf("*****************************************************\n");
    
    errno = LoadGmshInput1(data,bound,filename,info);
  }     

  return(errno);
}


int LoadGeoInput(struct FemType *data,struct BoundaryType *bound,
		 char *filename,int info)
{
  int noknots = 0,noelements = 0,maxnodes,dim;
  int elemind[MAXNODESD2],elementtype;
  int i,j,k,allocated,*revindx=NULL,maxindx;
  int elemnodes,maxelemtype,elemtype0;
  int usetaggeom,tagmat;
  FILE *in;
  char *cp,line[MAXLINESIZE];


  if ((in = fopen(filename,"r")) == NULL) {
    printf("LoadGeoInput: The opening of the mesh file %s failed!\n",filename);
    return(1);
  }
  if(info) printf("Loading mesh in geo format from file %s\n",filename);

  allocated = FALSE;
  dim = 3;
  maxnodes = 0;
  maxindx = 0;
  maxelemtype = 0;
  usetaggeom = FALSE;

omstart:


  for(;;) {
    if(Getrow(line,in,FALSE)) goto end;
    if(line[0]=='\0') goto end;
    if(strstr(line,"$End")) continue;

    if(strstr(line,"TYPES")) {
      if(!strstr(line,"ALL=TET04")) {
	printf("Only all tets implemnted at the monment!\n");
	return(1);
      }
      elemtype0 = 504;
      GETLINE;
    }
    else if(strstr(line,"COORDINATES")) {
      i = 0;
      for(;;) {
	GETLINE;
	if( strstr(line,"END_COORDINATES")) break;
	cp = line;
	j = next_int(&cp);
	i = i + 1;
	if(allocated) {
	  if(maxindx > noknots) revindx[j] = i;
	  data->x[i] = next_real(&cp);
	  data->y[i] = next_real(&cp);
	  if(dim > 2) data->z[i] = next_real(&cp);
	}
	else {
	  maxindx = MAX(j,maxindx);
	}
      }
      noknots = i;
    }
    else if(strstr(line,"ELEMENTS")) {
      i = 0;
      elementtype = elemtype0;
      tagmat = 1;

      for(;;) {
	GETLINE;
	if( strstr(line,"END_ELEMENTS")) break;
	cp = line;
	j = next_int(&cp);
	i = i + 1;

	if(allocated) {
	  elemnodes = elementtype % 100;
	  data->elementtypes[i] = elementtype;
	  data->material[i] = tagmat;
	  for(k=0;k<elemnodes;k++)
	    elemind[k] = next_int(&cp);
	  for(k=0;k<elemnodes;k++)
	    data->topology[i][k] = elemind[k];
	}
	else {
	  maxelemtype = MAX(maxelemtype,elementtype);
	}
      }
      noelements = i;
    }
    else if ( strstr(line,"BOUNDARIES")) {
      for(;;) {
	GETLINE;	
	if( strstr(line,"END_BOUNDARIES")) break;

	printf("Implement boundaries!\n");
      }      
    }
  }

 end:


  if(!allocated) {
    maxnodes = maxelemtype % 100;
    InitializeKnots(data);
    data->dim = dim;
    data->maxnodes = maxnodes;
    data->noelements = noelements;
    data->noknots = noknots;

    if(info) printf("Allocating for %d knots and %d elements.\n",noknots,noelements);
    AllocateKnots(data);

    if(maxindx > noknots) {
      revindx = Ivector(1,maxindx);
      for(i=1;i<=maxindx;i++) revindx[i] = 0;
    }
    rewind(in);
    allocated = TRUE;
    goto omstart;
  }

  if(maxindx > noknots) {
    printf("Renumbering the Geo nodes from %d to %d\n",maxindx,noknots);

    for(i=1; i <= noelements; i++) {
      elementtype = data->elementtypes[i];
      elemnodes = elementtype % 100; 

      for(j=0;j<elemnodes;j++) {
	k = data->topology[i][j];
	if(k <= 0 || k > maxindx) 
	  printf("index out of bounds %d\n",k);
	else if(revindx[k] <= 0) 
	  printf("unkonwn node %d %d in element %d\n",k,revindx[k],i);
	else 
	  data->topology[i][j] = revindx[k];
      }      
    }
    free_Ivector(revindx,1,maxindx);
  }

  if(0) ElementsToBoundaryConditions(data,bound,FALSE,info);

  if(info) printf("Successfully read the mesh from the Geo input file.\n");

  return(0);
}


/* Mapping between the element type of Universal file format and 
   ElmerSolver element type. */
static int UnvToElmerType(int unvtype)
{
  int elmertype;

  switch (unvtype) {

  case 11: 
  case 21:
    elmertype = 202;
    break;

  case 22:
  case 23:
    elmertype = 203;
    break;

  case 41:
  case 51:
  case 61:
  case 74:
  case 81:
  case 91:
    elmertype = 303;
    break;

  case 42:
  case 52:
  case 62:
  case 72:
  case 82:
  case 92:
    elmertype = 306;
    break;

  case 43:
  case 53:
  case 63:
  case 73:
  case 93:
    elmertype = 310;
    break;

  case 44:
  case 54:
  case 64:
  case 71:
  case 84:
  case 94:
    elmertype = 404;
    break;

  case 45:
  case 46:
  case 56:
  case 66:
  case 76:
  case 96:
    elmertype = 408;
    break;

  case 111:
    elmertype = 504;
    break;

  case 118:
    elmertype = 510;
    break;

  case 112:
    elmertype = 706;
    break;

  case 113:
    elmertype = 715;
    break;

  case 115:
    elmertype = 808;
    break;

  case 116:
    elmertype = 820;
    break;

  default:
    elmertype = 0;
    if(0) printf("Unknown elementtype in universal mesh format: %d\n",unvtype);
  }

  return(elmertype);
}


/* The Universal format supports something as "degenerated" elements. 
   This means that the same node is given multiple times in the element
   topology */
static int UnvRedundantIndexes(int nonodes,int *ind)
{
  int i,j,redundant;
  
  redundant = FALSE;
  for(i=0;i<nonodes;i++) {
    if( ind[i] == 0 ) redundant = TRUE;
    for(j=i+1;j<nonodes;j++) 
      if(ind[i] == ind[j]) redundant = TRUE;
  }
  if( redundant ) {
    printf("Redundant element %d: ",nonodes);
    for(i=0;i<nonodes;i++)
      printf(" %d ",ind[i]);
    printf("\n");
  }

  return(redundant);
}


/* Mapping between the elemental node order of Universal file format to 
   Elmer file format. */
static void UnvToElmerIndx(int elemtype,int *topology)
{
  int i=0,nodes=0,oldtopology[MAXNODESD2];
  int reorder, *porder;

  int order203[]={1,3,2};
  int order306[]={1,3,5,2,4,6};
  int order510[]={1,3,5,10,2,4,6,7,8,9};
  int order408[]={1,3,5,7,2,4,6,8};
  int order820[]={1,3,5,7,13,15,17,19,2,4,6,8,9,10,11,12,14,16,18,20};


  reorder = FALSE;

  switch (elemtype) {
      
  case 203:        
    reorder = TRUE;
    porder = &order203[0];
    break;

  case 306:        
    reorder = TRUE;
    porder = &order306[0];
    break;

  case 510:        
    reorder = TRUE;
    porder = &order510[0];
    break;

  case 408:        
    reorder = TRUE;
    porder = &order408[0];
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
      topology[i] = oldtopology[porder[i]-1];
  }
}




int LoadUniversalMesh(struct FemType *data,struct BoundaryType *bound,
		      char *prefix,int info)
/* Load the grid in universal file format. This format includes thousands of possible 
   fields and hence the parser can never be exhaustive. Just the mostly common used
   fields in FE community are treated. */
{
  int noknots,totknots,noelements,elemcode,maxnodes;
  int allocated,dim,ind,lines;
  int reordernodes,reorderelements,nogroups,maxnodeind,maxelem,elid,unvtype,elmertype;
  int nonodes,group,grouptype,mode,nopoints,nodeind,matind,physind,colorind;
  int minelemtype,maxelemtype,physoffset=0,doscaling=FALSE;
  int debug,mingroup,maxgroup,minphys,maxphys,nogroup,noentities,dummy;
  int *u2eind=NULL,*u2eelem=NULL;
  int *elementtypes;
  char filename[MAXFILESIZE],line[MAXLINESIZE],*cp;
  int i,j,k;
  char entityname[MAXNAMESIZE];
  Real scaling[4];
  FILE *in;


  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"unv");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadUniversalMesh: opening of the universal mesh file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }
 
  printf("Reading mesh from universal mesh file %s.\n",filename);
  InitializeKnots(data);

  dim = 3;
  allocated = FALSE;
  reordernodes = FALSE;
  reorderelements = FALSE;

  debug = FALSE;
  if( debug ){
    elementtypes = Ivector(0,820);
    for(i=0;i<=820;i++) elementtypes[i] = FALSE;
  }

  maxnodeind = 0;
  maxnodes = 0;
  maxelem = 0;
  mingroup = INT_MAX;
  maxgroup = 0;
  minphys = INT_MAX;
  maxphys = 0;
    
omstart:

  /* this is a global variable in the module */
  linenumber = 0;

  if(info) {
    if(allocated) 
      printf("Second round for reading data\n");
    else 
      printf("First round for allocating data\n");
  }

  noknots = 0;
  noelements = 0;
  nogroups = 0;
  nopoints = 0;
  group = 0;
  mode = 0;


  for(;;) { 

    if(0) printf("line: %d  %s\n",mode,line);

  nextline:
    if( !strncmp(line,"    -1",6)) mode = 0;
    if( Getrow(line,in,FALSE)) {
      goto end;
    }
    if(line[0]=='\0') goto end;

    if( !strncmp(line,"    -1",6)) mode = 0;
    else if( !strncmp(line,"  2411",6)) mode = 2411;
    else if( !strncmp(line,"  2412",6)) mode = 2412;
    else if( !strncmp(line,"  2420",6)) mode = 2420;
    else if( !strncmp(line,"  2467",6)) mode = 2467;
    else if( !strncmp(line,"  2435",6)) mode = 2435;
    else if( !strncmp(line,"   781",6)) mode = 781;
    else if( !strncmp(line,"   780",6)) mode = 780;
    else if( !strncmp(line,"   164",6)) mode = 164;
    else if( allocated && strncmp(line,"      ",6)) printf("Unknown mode line %d: %s",linenumber,line);


    if(debug && mode) printf("Current mode is %d\n",mode);

    /* node definition */
    if( mode == 2411 || mode == 781 ) {
      if(allocated && info) printf("Reading node coordinates\n");
      for(;;) {
	GetrowDouble(line,in);
	if( !strncmp(line,"    -1",6)) {
	  if(!allocated && info) printf("There are %d nodes in the mesh\n",noknots);
	  goto nextline;
	}

	cp = line;
	nodeind = next_int(&cp);
	/* Three other fields omitted: two coordinate systems and color */
	noknots += 1;
	GetrowDouble(line,in);

	if(allocated) {
	  if(reordernodes) {
	    if(u2eind[nodeind]) 
	      printf("Reordering node %d already set (%d vs. %d)\n",
		     nodeind,u2eind[nodeind],noknots);
	    else
	      u2eind[nodeind] = noknots;
	  }

	  cp = line;
	  data->x[noknots] = next_real(&cp);
	  data->y[noknots] = next_real(&cp);
	  data->z[noknots] = next_real(&cp);
	}
	else {
	  if(nodeind != noknots) reordernodes = TRUE;
	  maxnodeind = MAX(maxnodeind,nodeind);
	}
      }
    }

    if( mode == 2412 ) {
      minelemtype = INT_MAX;
      maxelemtype = 0;

      if(allocated && info) printf("Reading element topologies\n");
      for(;;) {
	Getrow(line,in,FALSE);
	if( strstr(line,"-1")) {
	  if(info && !allocated) printf("Element type range in mesh [%d,%d]\n",minelemtype,maxelemtype);
	  goto nextline;
	}	

	noelements += 1;
	cp = line;
	elid = next_int(&cp);
	unvtype = next_int(&cp);
	physind = next_int(&cp);
	matind = next_int(&cp);
	colorind = next_int(&cp);
	nonodes = next_int(&cp);

	if(!allocated ) {
	  if(0) printf("elem = %d %d %d %d\n",noelements,unvtype,physind,matind);
	}	

	if (!allocated) {
	  minphys = MIN( minphys, physind );
	  maxphys = MAX( maxphys, physind );	 
	  maxnodes = MAX(maxnodes, nonodes);
	  if(elid != noelements) reorderelements = TRUE;
	  maxelem = MAX(maxelem, elid);
	}
	
	if(unvtype == 11 || unvtype == 21 || unvtype == 22 ) Getrow(line,in,FALSE);
	Getrow(line,in,FALSE);
	cp = line;

	elmertype = UnvToElmerType(unvtype); 
	if(!elmertype) {
	  printf("Unknown elementtype %d %d %d %d %d %d %d\n",
		 noelements,elid,unvtype,physind,matind,colorind,nonodes);
	  printf("line %d: %s\n",linenumber,line);
	  bigerror("done");
	}

	if(elmertype == 510 ) 	   
	  lines = 1;
	else if(elmertype == 820 ) 
	  lines = 2;
	else
	  lines = 0;

	if(allocated) {
	  if(reorderelements) u2eelem[elid] = noelements;

	  if(debug && !elementtypes[elmertype]) {
	    elementtypes[elmertype] = TRUE;
	    printf("new elementtype in elmer: %d (unv: %d)\n",elmertype,unvtype);
	  }

	  if(elmertype % 100 != nonodes) {
	    printf("nonodes = %d elemtype = %d elid = %d\n",nonodes,elmertype,elid);
	    nonodes = elmertype % 100;
	  }
	  
	  data->elementtypes[noelements] = elmertype;
	  for(i=0;i<nonodes;i++) {
	    if( lines > 0 && i >= 8 ) {
	      if( i%8 == 0 ) {
		Getrow(line,in,FALSE);
		cp = line;
	      }
	    }
	    data->topology[noelements][i] = next_int(&cp);
	  }

	  UnvRedundantIndexes(nonodes,data->topology[noelements]);

	  UnvToElmerIndx(elmertype,data->topology[noelements]);	  

	  /* should this be physical property or material property? */
	  data->material[noelements] = physind + physoffset;
	}
	else {
	  minelemtype = MIN( minelemtype, elmertype );
	  maxelemtype = MAX( maxelemtype, elmertype );
	  for(i=1;i<=lines;i++) {
	    Getrow(line,in,FALSE);
	  }	  
	}
      }
    }

    if( mode == 2420 ) {
      int partuid,coordlabel,coordtype;
      Real coeff;
      if(allocated && info) printf("Reading Coordinate system information\n");

      Getrow(line,in,FALSE);
      if( !allocated ) {
	cp = line;
	partuid = next_int(&cp);
	printf("Part UID = %d\n",partuid);
      }	
      Getrow(line,in,FALSE);
      if(!allocated ) {
	sscanf(line,"%s",entityname);
	printf("Part name = %s\n",entityname);
      }      
      Getrow(line,in,FALSE);
      if( !allocated ) {
	cp = line;
	coordlabel = next_int(&cp);
	coordtype = next_int(&cp);
	if( coordtype != 0 ) {
	  printf("Coordinate system is not cartesian: %d\n",coordtype);
	  printf("Code some more if you want to consider this!\n");
	}
      }      

      Getrow(line,in,FALSE);
      if(!allocated ) {
	sscanf(line,"%s",entityname);
	printf("Coord system name = %s\n",entityname);
      }      
      for(i=1;i<=4;i++) {
	Getrow(line,in,FALSE);
	if( !allocated ) {
	  cp = line;
	  if(!cp) printf("Problem reading line %d for coordinate system\n",i);
	  for(j=1;j<= 3;j++) {
	    coeff = next_real(&cp);
	    if( i == j ) {
	      scaling[i] = coeff;
	      if( fabs(coeff) < 1.0e-20) {
		printf("Scaling for component %d too small %le\n",i,coeff);
	      }
	      else if( fabs(coeff-1.0) ) {
		doscaling = TRUE;
		printf("Scaling component %d by %le\n",i,coeff);
	      }
	    }
	    else {
	      if(fabs(coeff) > 1.0e-20 ) {
		printf("Transformation matrix is not diagonal %d%d: %e\n",i,j,coeff);
		smallerror("Code some more...");
	      }
	    }
	  }	  
	}      
      }
      Getrow(line,in,FALSE);
      if( strncmp(line,"    -1",6)) 
	printf("Field 2420 should already be ending: %s\n",line); 
      goto nextline;
    }

    if( mode == 780 ) {
      int physind2,matind2;
      maxelemtype = 0;
      minelemtype = 1000;

      if(allocated && info) printf("Reading element groups in mode %d\n",mode);
      for(;;) {
	Getrow(line,in,FALSE);
	if( !strncmp(line,"    -1",6)) goto nextline;
	
	noelements += 1;
	cp = line;
	elid = next_int(&cp);
	unvtype = next_int(&cp);

	physind = next_int(&cp);
	physind2 = next_int(&cp);
	matind = next_int(&cp);
	matind2 = next_int(&cp);
	colorind = next_int(&cp);
	nonodes = next_int(&cp);
	
	if (!allocated) {
	  maxnodes = MAX(maxnodes, nonodes);
	  if(elid != noelements) reorderelements = TRUE;
	  maxelem = MAX(maxelem, elid);
	  minphys = MIN( minphys, physind );
	  maxphys = MAX( maxphys, physind );
	}
	
	if(unvtype == 11 || unvtype == 21) Getrow(line,in,FALSE);
	Getrow(line,in,FALSE);
	cp = line;
	if(allocated) {
	  if(reorderelements) u2eelem[elid] = noelements;

	  elmertype = UnvToElmerType(unvtype); 
	  maxelemtype = MAX( maxelemtype, elmertype );
	  minelemtype = MIN( minelemtype, elmertype );

	  if(debug && !elementtypes[elmertype]) {
	    elementtypes[elmertype] = TRUE;
	    printf("new elementtype in elmer: %d (unv: %d)\n",elmertype,unvtype);
	  }

	  if(elmertype % 100 != nonodes) {
	    printf("nonodes = %d elemtype = %d elid = %d\n",nonodes,elmertype,elid);
	    nonodes = elmertype % 100;
	  }

	  data->elementtypes[noelements] = elmertype;
	  for(i=0;i<nonodes;i++)
	    data->topology[noelements][i] = next_int(&cp);

	  UnvRedundantIndexes(nonodes,data->topology[noelements]);

	  UnvToElmerIndx(elmertype,data->topology[noelements]);	  

	  /* should this be physical property or material property? */
	  data->material[noelements] = physind + physoffset;
	}
      }
    }  

    if( mode == 2467 || mode == 2435) {
      if(allocated && info) printf("Reading element groups in mode %d\n",mode);
      
      for(;;) {
	Getrow(line,in,FALSE);
	if( !strncmp(line,"    -1",6)) goto nextline;
	
	cp = line;
	nogroup = next_int(&cp);
	maxelemtype = 0;
	minelemtype = 1000;
	for(i=1;i<=6;i++)
	  dummy = next_int(&cp);
	noentities = next_int(&cp);

	if(!allocated) {
	  mingroup = MIN( mingroup, nogroup );
	  maxgroup = MAX( maxgroup, nogroup );
	}

	Getrow(line,in,FALSE);	
	if( !strncmp(line,"    -1",6)) goto nextline;
	
	/* Used for the empty group created by salome */
	/* if( mode == 2467 && !strncmp(line,"            ",12)) continue; */
	
	group++;
	k = 0;
	if(allocated) {
	  sscanf(line,"%s",entityname);
	  strcpy(data->bodyname[nogroup],entityname);
	  data->bodynamesexist = TRUE;
	  data->boundarynamesexist = TRUE;

	  if(info) printf("Reading %d:th group with index %d with %d entities: %s\n",
			  group,nogroup,noentities,entityname);
	}
	if(noentities == 0) Getrow(line,in,FALSE);

	for(i=0;i<noentities;i++) {
	  if(i%2 == 0) {
	    Getrow(line,in,FALSE);
	    if( !strncmp(line,"    -1",6)) goto nextline;
	    cp = line;
	  }	  

	  grouptype = next_int(&cp);
	  ind = next_int(&cp);
	  dummy = next_int(&cp);
	  dummy = next_int(&cp);

	  if(ind == 0) continue;

	  if( grouptype == 8 ) {

	    if(allocated) {
	      if(reorderelements) ind = u2eelem[ind];
	      elemcode = data->elementtypes[ind];
	      maxelemtype = MAX( maxelemtype, elemcode );
	      minelemtype = MIN( minelemtype, elemcode );
	      data->material[ind] = nogroup;
	    }
	  }
	  else if(grouptype == 7) {
	    nopoints += 1;
	    if(allocated) {
	      elemcode = 101;
	      data->material[noelements+nopoints] = nogroup;
	      maxelemtype = MAX( maxelemtype, elemcode );
	      minelemtype = MIN( minelemtype, elemcode );
	      data->elementtypes[noelements+nopoints] = elemcode;	      
	      data->topology[noelements+nopoints][0] = ind;
	    }
	  }
	  else {
	    printf("unknown group type %d\n",grouptype);
	  }
	}
	if(allocated && info) {
	  printf("Element type range in group is [%d %d]\n",minelemtype,maxelemtype);
	}

      }
    }

    if( mode == 164 ) {
      if(!allocated) printf("Units dataset content is currently omitted!\n");
      for(;;) {
	Getrow(line,in,FALSE);
	if( !strncmp(line,"    -1",6)) 
	goto nextline;
      }
    }

  }

end:

  exit;
  if(0) printf("Done reading mesh\n");


  if(!allocated) {

    if(reordernodes) {
      if(info) printf("Reordering %d nodes with indexes up to %d\n",noknots,maxnodeind);
      u2eind = Ivector(1,maxnodeind);
      for(i=1;i<=maxnodeind;i++) u2eind[i] = 0;
    }
    if(reorderelements) {
      if(info) printf("Reordering %d elements with indexes up to %d\n",noelements,maxelem);
      u2eelem = Ivector(1,maxelem);
      for(i=1;i<=maxelem;i++) u2eelem[i] = 0;
    }

    if(noknots == 0 || noelements == 0 || maxnodes == 0) {
      printf("Invalid mesh consists of %d nodes and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);     
      fclose(in);
      return(2);
    }

    rewind(in);
    totknots = noknots;
    data->noknots = noknots;
    data->noelements = noelements + nopoints;
    data->maxnodes = maxnodes;
    data->dim = dim;
    
    if(info) {
      printf("Allocating mesh with %d nodes and %d %d-node elements in %d dims.\n",
	     noknots,noelements,maxnodes,dim);
    }  
    AllocateKnots(data);
    allocated = TRUE;

    /* Set an offset for physical indexes so that the defined groups and 
       existing physical indexes won't mix confusingly */
    if( maxphys >= mingroup && minphys <= maxgroup ) {
      physoffset = maxgroup - minphys + 1;
    }
    else {
      physoffset = 0;
    }

    if(info) {
      printf("Physical index interval is [%d,%d]\n",minphys,maxphys);
      if( maxgroup ) 
	printf("Group index interval is [%d,%d]\n",mingroup,maxgroup);
      if(physoffset) printf("Using offset %d for physical indexes\n",physoffset);
    }


    goto omstart;    
  }
  fclose(in);

  /* If the physical index may be zero, then we have a risk that there is 
     an unset material index. Elmer does not like material indexes of zeros. 
     This could be made prettier as now the almost same thing is done twice. */
  if( minphys + physoffset == 0 ) {
    mingroup = INT_MAX;
    maxgroup = 0;
    for(i=1;i<=data->noelements;i++) {
      mingroup = MIN( mingroup, data->material[i] );
      maxgroup = MAX( maxgroup, data->material[i] );
    }
    if( mingroup == 0 ) {
      if(info) {
	if(!maxgroup) printf("No material groups were successfully applied\n");
	printf("Unset elements were given material index %d\n",maxgroup+1);    
      }
      for(i=1;i<=data->noelements;i++) 
	if(data->material[i] == 0) data->material[i] = maxgroup + 1;
    }
  }    

  /* Elmer likes that node indexes are given so that no integers are missed.
     If this is not the case we need to do renumbering of nodes. */
  if(reordernodes) {
    printf("Reordering nodes continuously\n");
    for(j=1;j<=noelements;j++)
      for(i=0;i<data->elementtypes[j]%100;i++)
	data->topology[j][i] = u2eind[data->topology[j][i]];
    free_Ivector(u2eind,1,maxnodeind);
  }
  if(reorderelements) {
    free_Ivector(u2eelem,1,maxelem);
  }

  /* Do scaling if requested */
  if( doscaling ) {
    Real *coord;
    for(j=1;j<=3;j++) {
      if( j == 1 ) 
	coord = data->x;
      else if( j == 2 ) 
	coord = data->y;
      else 
	coord = data->z;
      
      if( fabs(scaling[j]-1.0) >= 1.0e-20 ) {
	for(i=1;i<=noknots;i++)
	  coord[i] *= scaling[j];
      }
    }
  }

  /* This is here for debugging of the nodal order */
  if(FALSE) for(j=1;j<=noelements;j++) {
    int elemtype = data->elementtypes[j];
    printf("element = %d\n",j);
    for(i=0;elemtype%100;i++) {
      k = data->topology[j][i];
      printf("node i=%d %.3le %.3le %.3le\n",i,data->x[k],data->z[k],data->y[k]);
    }
  }
      
  
  /* Until this far all elements have been listed as bulk elements. 
     Now separate the lower dimensional elements to be boundary elements. */
  ElementsToBoundaryConditions(data,bound,TRUE,info);
 
  if(info) printf("The Universal mesh was loaded from file %s.\n\n",filename);

  return(0);
}





int LoadCGsimMesh(struct FemType *data,char *prefix,int info)
/* Load the mesh from postprocessing format of CGsim */
{
  int noknots,noelements,maxnodes,material,allocated,dim,debug,thismat,thisknots,thiselems;
  char filename[MAXFILESIZE],line[MAXLINESIZE],*cp;
  int i,j,inds[MAXNODESD2],savedofs;
  Real dummyreal;
  FILE *in;


  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"plt");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadCGsimMesh: opening of the CGsim mesh file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }

  printf("Reading mesh from CGsim mesh file %s.\n",filename);
  InitializeKnots(data);

  debug = FALSE;
  allocated = FALSE;
  savedofs = FALSE;

omstart:

  maxnodes = 4;
  noknots = 0;
  noelements = 0;
  material = 0;
  dim = 2;
  thismat = 0;


  for(;;) {
    
    if(Getrow(line,in,FALSE)) goto end;
    if(line[0]=='\0') goto end;
    
    cp = strstr(line,"ZONE");    
    if(!cp) continue;

    thismat += 1;
    cp = strstr(line," N=");
    cp += 3;
    thisknots = next_int(&cp);

    cp = strstr(line,",E=");
    cp += 3;
    thiselems = next_int(&cp);

    if(debug) {
      printf("%s",line);
      printf("thismat = %d knots = %d elems = %d\n",thismat,thisknots,thiselems);
    }

    for(i=1;i<=thisknots;i++) {
      GETLINE;

      if(allocated) {
	cp = line;
	data->x[noknots+i] = next_real(&cp);
	data->y[noknots+i] = next_real(&cp);
	data->z[noknots+i] = 0.0;

	if(savedofs == 1) {
	  for(j=1;j<=4;j++)
	    dummyreal = next_real(&cp);	    
	  data->dofs[1][noknots+i] = next_real(&cp);
	}
	else if(savedofs == 5) {
	  for(j=1;j<=5;j++)	  
	    data->dofs[j][noknots+i] = next_real(&cp);	  
	}

      }
    }

    for(i=1;i<=thiselems;i++) {
      GETLINE;

      if(allocated) {
	cp = line;
	for(j=0;j<4;j++) 
	  inds[j] = next_int(&cp);
	for(j=0;j<4;j++) 
	  data->topology[noelements+i][j] = inds[j]+noknots;
	if(inds[2] == inds[3]) 
	  data->elementtypes[noelements+i] = 303;
	else
	  data->elementtypes[noelements+i] = 404;	  
	data->material[noelements+i] = thismat;
      }
   }

    noknots += thisknots;
    noelements += thiselems;
  }

 end:

  if(!allocated) {
    if(noknots == 0 || noelements == 0 || maxnodes == 0) {
       printf("Invalid mesh consits of %d knots and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);     
       fclose(in);
       return(2);
    }

    rewind(in);
    data->noknots = noknots;
    data->noelements = noelements;
    data->maxnodes = maxnodes;
    data->dim = dim;

    
    if(info) {
      printf("Allocating for %d knots and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);
    }  
    AllocateKnots(data);

    if(savedofs == 1) {
      CreateVariable(data,1,1,0.0,"Temperature",FALSE);
    }
    else if(savedofs == 5) {
      CreateVariable(data,1,1,0.0,"dTdX",FALSE);
      CreateVariable(data,2,1,0.0,"dTdY",FALSE);
      CreateVariable(data,3,1,0.0,"Qx",FALSE);
      CreateVariable(data,4,1,0.0,"Qy",FALSE);
      CreateVariable(data,5,1,0.0,"Temperature",FALSE);
    }

    allocated = TRUE;
    goto omstart;    
  }
  fclose(in);

  if(info) printf("The CGsim mesh was loaded from file %s.\n\n",filename);
  return(0);
}


int FluxToElmerType(int nonodes, int dim) {
  int elmertype;

  elmertype = 0;
  
  if( dim == 2 ) {
    switch( nonodes ) {
    case 3: 
      elmertype = 203;
      break;
    case 6:
      elmertype = 306;
      break;
    case 8:
      elmertype = 408;
      break;
    }
  }
  
  if( !elmertype ) printf("FluxToElmerType could not deduce element type! (%d %d)\n",nonodes,dim);

  return(elmertype);
}





int LoadFluxMesh(struct FemType *data,struct BoundaryType *bound,
		    char *prefix,int info)
/* Load the mesh from format of Flux Cedrat in TRA format. */
{
  int noknots,noelements,maxnodes,dim,elmertype;
  int nonodes,matind,noregions,mode;
  int debug;
  int *elementtypes;
  char filename[MAXFILESIZE],line[MAXLINESIZE],*cp;
  int i,j,k;
  char entityname[MAXNAMESIZE];
  FILE *in;


  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"TRA");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadFluxMesh: opening of the Flux mesh file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }
 
  printf("Reading 2D mesh from Flux mesh file %s.\n",filename);
  InitializeKnots(data);

  debug = FALSE;
  linenumber = 0;
  dim = 2;
  noknots = 0;
  noelements = 0;
  mode = 0;
  maxnodes = 8;



  for(;;) { 

    if(0) printf("line: %d  %s\n",mode,line);

    if( Getrow(line,in,FALSE)) goto end;
    if(line[0]=='\0') goto end;

    if( strstr(line,"Number of nodes")) mode = 1;
    else if( strstr(line,"Total number of elements")) mode = 2;
    else if( strstr(line,"Total number of regions")) mode = 3;

    else if( strstr(line,"Description of elements")) mode = 10;
    else if( strstr(line,"Coordinates of the nodes")) mode = 11;
    else if( strstr(line,"Names of the regions")) mode = 12;

    else if( strstr(line,"Neighbouring element table")) mode = 13;
    else if( strstr(line,"List of boundary nodes")) mode = 14;
    else if( strstr(line,"Physical properties")) mode = 15;
    else if( strstr(line,"Boundary conditions")) mode = 16;
    else {
      if(debug) printf("Unknown mode line %d: %s",linenumber,line);
      mode = 0;
    }

    if(debug && mode) printf("Current mode is %d\n",mode);

    switch( mode ) {
    case 1: 
      noknots = atoi(line);
      break;

    case 2: 
      noelements = atoi(line);
      break;

    case 3: 
      noregions = atoi(line);
      break;

     
    case 10:
      if(info) {
	printf("Allocating mesh with %d nodes and %d %d-node elements in %d dims.\n",
	       noknots,noelements,maxnodes,dim);
      }  

      data->noknots = noknots;
      data->noelements = noelements;
      data->maxnodes = maxnodes;
      data->dim = dim;
      AllocateKnots(data);

      if(info) printf("Reading %d element topologies\n",noelements);
      for(i=1;i<=noelements;i++) {
	Getrow(line,in,FALSE);
	cp = line;
	j = next_int(&cp);
	if( i != j ) {
	  printf("It seems that reordering of elements should be performed! (%d %d)\n",i,j);
	}
	nonodes = next_int(&cp);
	matind = abs( next_int(&cp) ); 
	
	elmertype = FluxToElmerType( nonodes, dim );
	data->elementtypes[i] = elmertype;
	data->material[i] = matind;

	Getrow(line,in,FALSE);
	cp = line;
	for(k=0;k<nonodes;k++) {
	  data->topology[i][k] = next_int(&cp);
	}
      }
      break;

    case 11:
      if(info) printf("Reading %d element nodes\n",noknots);
      for(i=1;i<=noknots;i++) {
	Getrow(line,in,FALSE);
	cp = line;
	j = next_int(&cp);
	if( i != j ) {
	  printf("It seems that reordering of nodes should be performed! (%d %d)\n",i,j);
	}
	data->x[i] = next_real(&cp);
	data->y[i] = next_real(&cp);
	if(dim == 3) data->z[i] = next_real(&cp);
      }
      break;


    case 12:
      if(info) printf("Reading %d names of regions\n",noregions);
      for(i=1;i<=noregions;i++) {
	Getrow(line,in,FALSE);
	cp = line;
	j = next_int(&cp);
	if( i != j ) {
	  printf("It seems that reordering of regions should be performed! (%d %d)\n",i,j);
	}
	sscanf(cp,"%s",entityname);
	strcpy(data->bodyname[i],entityname);
      }
      data->bodynamesexist = TRUE;
      data->boundarynamesexist = TRUE;
      break;


    default:
      if(debug) printf("unimplemented mode: %d\n",mode );
      mode = 0;
      break;
    }
  }

 end:
  fclose(in);

  /* Until this far all elements have been listed as bulk elements. 
     Now separate the lower dimensional elements to be boundary elements. */
  ElementsToBoundaryConditions(data,bound,TRUE,info);
 
  if(info) printf("The Flux mesh was loaded from file %s.\n\n",filename);

  return(0);
}



/* Mapping between the elemental node order of PF3 file format to 
   Elmer file format. */
static void PF3ToElmerPermuteNodes(int elemtype,int *topology)
{
  int i=0, nodes=0, oldtopology[MAXNODESD2];
  int reorder, *porder;
  int debug;
  
  int order303[] = {3,1,2};                //tri
  int order306[] = {3,1,2,6,4,5};          //tri^2
  int order404[] = {3,4,1,2};             //quad
  int order408[] = {3,4,1,2,7,8,5,6};     //quad^2
  int order504[] = {1,2,3,4};             //tetra
  int order510[] = {1,2,3,4,5,8,6,7,10,9};//tetra^2
  int order605[] = {3,2,1,4,5};           //pyramid
  int order613[] = {3,2,1,4,5,7,6,9,8,12,11,10,13};           //pyramid^2
  int order706[] = {6,4,5,3,1,2};         //wedge (prism)  
  int order715[] = {6,4,5,3,1,2,12,10,11,9,7,8,15,13,14};   //wedge^2 (prism^2)  
  int order808[] = {7,8,5,6,3,4,1,2};     //hexa
  int order820[] = {7,8,5,6,3,4,1,2,15,16,13,14,19,20,17,18,11,12,9,10};  //hexa^2

  debug = TRUE;
  
  reorder = FALSE;

  switch (elemtype) {
  
  case 101:
    //nothing to change here
    break;
    
  case 202:
    //nothing to change here
    break;
    
  case 203:
    //nothing to change here
    break;
    
  case 303:
    reorder = TRUE;
    porder = &order303[0];
    break;
    
  case 306:
    reorder = TRUE;
    porder = &order306[0];
    break;
     
  case 404:        
    reorder = TRUE;
    porder = &order404[0];
    break;

  case 408:   
    reorder = TRUE;
    porder = &order408[0];
    break;
    
  case 504:        
    reorder = TRUE;
    porder = &order504[0];
    break;

  case 510:   
    reorder = TRUE;
    porder = &order510[0];
    break;
    
  case 605:        
    reorder = TRUE;
    porder = &order605[0];
    break;

  case 613:   
    reorder = TRUE;
    porder = &order613[0];
    break;

  case 706:        
    reorder = TRUE;
    porder = &order706[0];
    break;

  case 715:     
    reorder = TRUE;
    porder = &order715[0];
    break;
    
  case 808:        
    reorder = TRUE;
    porder = &order808[0];
    break;

  case 820:     
    reorder = TRUE;
    porder = &order820[0];
    break;    

  default:
      if(debug) printf("Warning : Unknown element type: %d\n",elemtype );
      break;
  }

  if( reorder ) {
    nodes = elemtype % 100;
    for(i=0;i<nodes;i++) 
      oldtopology[i] = topology[i];
    for(i=0;i<nodes;i++) 
      topology[i] = oldtopology[porder[i]-1];
  }
}


int FluxToElmerType3D(int nonodes, int dim) {
  int elmertype;

  elmertype = 0;
  
  if( dim == 2 ) {
    switch( nonodes ) {
    case 3: 
      elmertype = 303;
      break;
    case 4: 
      elmertype = 404;
      break;
    case 6:
      elmertype = 306;
      break;
    case 8:
      elmertype = 408;
      break;
    }
  }
  
  if( dim == 3 ) {
    switch( nonodes ) {
    case 4:
      elmertype = 504;
      break;
    case 5:
      elmertype = 605;
      break;
    case 6:
      elmertype = 706;
      break;
    case 8:
      elmertype = 808;
      break;
    case 10: 
      elmertype = 510;
      break;
    case 13:
      elmertype = 613;
      break;
    case 15: 
      elmertype = 715;
      break;
    case 20:
      elmertype = 820;
      break;
    }
  }
    
  if( !elmertype ) printf("FluxToElmerType3D could not deduce element type! (%d %d)\n",nonodes,dim);

  return(elmertype);
}


int LoadFluxMesh3D(struct FemType *data,struct BoundaryType *bound,
		    char *prefix,int info)
/* Load the mesh from format of Flux Cedrat in PF3 format. */
{
  int noknots,noelements,maxnodes,dim,elmertype;
  int nonodes,matind,noregions,mode;
  int dimplusone, maxlinenodes, nodecnt;
  int debug;
  int *elementtypes;
  char filename[MAXFILESIZE],line[MAXLINESIZE],*cp;
  int i,j,k;
  char entityname[MAXNAMESIZE];
  FILE *in;

  strcpy(filename,prefix);
  if ((in = fopen(filename,"r")) == NULL) {
    AddExtension(prefix,filename,"PF3");
    if ((in = fopen(filename,"r")) == NULL) {
      printf("LoadFluxMesh3D: opening of the Flux mesh file '%s' wasn't succesfull !\n",
	     filename);
      return(1);
    }
  }
 
  printf("Reading 3D mesh from Flux mesh file %s.\n",filename);
  InitializeKnots(data);

  debug = FALSE;
  linenumber = 0;
  dim = 3;
  noknots = 0;
  noelements = 0;
  mode = 0;
  maxnodes = 20;      // 15?
	maxlinenodes = 12; //nodes can be located at several lines

  for(;;) { 

    if(0) printf("line: %d  %s\n",mode,line);

    if( Getrow(line,in,FALSE)) goto end;
    if(line[0]=='\0') goto end;
    if( strstr(line,"==== DECOUPAGE  TERMINE")) goto end;

    if( strstr(line,"NOMBRE DE DIMENSIONS DU DECOUPAGE")) mode = 1;
    else if( strstr(line,"NOMBRE  D'ELEMENTS")) mode = 3;
    else if( strstr(line,"NOMBRE DE POINTS")) mode = 2;
    else if( strstr(line,"NOMBRE DE REGIONS")) mode = 4;

    else if( strstr(line,"DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS")) mode = 10;
    else if( strstr(line,"COORDONNEES DES NOEUDS")) mode = 11;
    else if( strstr(line,"NOMS DES REGIONS")) mode = 12;
    else {
      if(debug) printf("Unknown mode line %d: %s",linenumber,line);
      mode = 0;
    }

    if(debug && mode) printf("Current mode is %d\n",mode);

    switch( mode ) {
    case 1: 
      dim = atoi(line);
      break;

    case 2: 
      if( strstr(line,"NOMBRE DE POINTS D'INTEGRATION")) break;/* We are looking for the total number of nodes */
      noknots = atoi(line);
      break;

    case 3: 
      i = atoi(line);
      noelements = MAX(i,noelements); /* We are looking for the total number of elements */
      break;

    case 4: 
      i = atoi(line);
      noregions = MAX(i,noregions); /* We are looking for the total number of regions */
      break;

     
    case 10:
      if(info) {
	printf("Allocating mesh with %d nodes and %d %d-node elements in %d dims.\n",
	       noknots,noelements,maxnodes,dim);
      }  

      data->noknots = noknots;
      data->noelements = noelements;
      data->maxnodes = maxnodes;
      data->dim = dim;
      AllocateKnots(data);

      if(info) printf("Reading %d element topologies\n",noelements);
      for(i=1;i<=noelements;i++) 
      {
	Getrow(line,in,FALSE);
	cp = line;
	j = next_int(&cp);
	if( i != j ) {
	  printf("It seems that reordering of elements should be performed! (%d %d)\n",i,j);
	}
	next_int(&cp);              //2 internal elemnt type description
	next_int(&cp);              //3 internal elemnt type description
	matind = next_int(&cp);     //4 number of the belonging region
	dimplusone = next_int(&cp); //5 dimensiality 4-3D 3-2D
	next_int(&cp);              //6 zero here always
	next_int(&cp);              //7 internal elemnt type description
	nonodes = next_int(&cp);    //8 number of nodes
		
	elmertype = FluxToElmerType3D( nonodes, dimplusone-1 );
	data->elementtypes[i] = elmertype;
	data->material[i] = matind;

	Getrow(line,in,FALSE);
	cp = line;
	nodecnt = 0;
	for(k=0;k<nonodes;k++) {

	  if(nodecnt >= maxlinenodes) {
	    nodecnt = 0;
	  	Getrow(line,in,FALSE);
	    cp = line;
	  }
	  data->topology[i][k] = next_int(&cp);
	  nodecnt+=1;
	}
	
	PF3ToElmerPermuteNodes(elmertype,data->topology[noelements]);	 
	
      }
      break;

    case 11:
      if(info) printf("Reading %d element nodes\n",noknots);
      for(i=1;i<=noknots;i++) {
	Getrow(line,in,FALSE);
	cp = line;
	j = next_int(&cp);
	if( i != j ) {
	  printf("It seems that reordering of nodes should be performed! (%d %d)\n",i,j);
	}
	data->x[i] = next_real(&cp);
	data->y[i] = next_real(&cp);
	data->z[i] = next_real(&cp);
      }
      break;


    case 12:
      if(info) printf("Reading %d names of regions\n",noregions);
      for(i=1;i<=noregions;i++) {
	Getrow(line,in,FALSE);

	/* currently we just cycle trough this and get a new row */
	if( strstr(line,"REGIONS SURFACIQUES")) Getrow(line,in,FALSE);
	if( strstr(line,"REGIONS VOLUMIQUES")) Getrow(line,in,FALSE);

	sscanf(line,"%s",entityname);
	strcpy(data->bodyname[i],entityname);
      }
      data->bodynamesexist = TRUE;
      data->boundarynamesexist = TRUE;
      break;


    default:
      if(debug) printf("unimplemented mode: %d\n",mode );
      mode = 0;
      break;
    }
  }

 end:
  fclose(in);

  /* Until this far all elements have been listed as bulk elements. 
     Now separate the lower dimensional elements to be boundary elements. */
  ElementsToBoundaryConditions(data,bound,TRUE,info);
 
  if(info) printf("The Flux 3D mesh was loaded from file %s.\n\n",filename);

  return(0);
}
