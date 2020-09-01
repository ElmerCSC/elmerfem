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
#include <sys/stat.h>
#include <sys/types.h>

#include "egutils.h"
#include "egdef.h"
#include "egtypes.h"
#include "egmesh.h"
#include "egnative.h"
#include "../config.h"

#define getline fgets(line,MAXLINESIZE,in) 


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

    if(0) printf("Loading steps i=%d  j=%d  k=%d  r=%.3g\n",i,j,k,r);

    for(i=1; i <= noknots; i++) 
      for(j=1;j <= novctrs;j++) 
	iostat = fscanf(in,"%le",&(data->dofs[j][k*noknots+i]));
  }

end:
  data->timesteps = timesteps;

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

  getline;
  sscanf(line,"%d %d %d",&noknots,&noelements,&nosides);
  getline;
  sscanf(line,"%d",&tottypes);

  maxelemtype = 0;
  maxnodes = 0;
  for(i=1;i<=tottypes;i++) {   
    getline;
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
    getline;
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
	if(j < 0 || j > MAXBODIES ) {
	  printf("Cannot treat names for body %d\n",j);
	  nameproblem = TRUE;
	}
	else {
	  strcpy(data->bodyname[j],line2);	
	  data->bodynamesexist = TRUE;
	}
      }
      else {
	if(j < 0 || j > MAXBOUNDARIES ) {
	  printf("Cannot treat names for boundary %d\n",j);
	  nameproblem = TRUE;
	}
	else {
	  strcpy(data->boundaryname[j],line2);	
	  data->boundarynamesexist = TRUE;
	}
      }
    }
    namesend:

    if( nameproblem ) {
      data->boundarynamesexist = FALSE;
      data->bodynamesexist = FALSE;
      printf("Warning: omitting use of names because the indexes are beyond range, code some more...\n");
    }
  }


  if(!cdstat) cdstat = chdir("..");

  if(info) printf("Elmer mesh loaded successfully\n");

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
	  fprintf(out,"%.6g ",data->dofs[j][k*noknots+i]);
	if(data->edofs[j] == 2) 
	  fprintf(out,"%.6g %.6g 0.0 ",
		  data->dofs[j][2*(k*noknots+i)-1],data->dofs[j][2*(k*noknots+i)]);
	if(data->edofs[j] == 3) 
	  fprintf(out,"%.6g %.6g %.6g ",
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


int SaveElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int decimals,int nooverwrite, int info)
/* Saves the mesh in a form that may be used as input 
   in Elmer calculations. 
   */
#define MAXELEMENTTYPE 827
{
  int noknots,noelements,material,sumsides,elemtype,fail,cdstat;
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
      
      if(bound[j].types[i] < MAXBCS) usedbc[bound[j].types[i]] += 1;

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
	if(usedbody[i]) fprintf(out,"$ %s = %d\n",data->bodyname[i],i);
    }     
    if(data->boundarynamesexist) {
      fprintf(out,"! ----- names for boundaries -----\n");
      for(i=1;i<MAXBCS;i++) 
	if(usedbc[i]) fprintf(out,"$ %s = %d\n",data->boundaryname[i],i);
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
	  fprintf(out,"  Name = %s\n",data->bodyname[i]);
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
	  fprintf(out,"  Name = %s\n",data->boundaryname[i]);
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
