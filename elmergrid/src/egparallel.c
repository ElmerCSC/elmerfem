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

/* -------------------------------:  egparallel.c  :----------------------------
   This module includes routines related to parallel file format of Elmer suite. 
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
#include "egparallel.h"
#include "../config.h"

#if USE_METIS
#include "metis-5.1.0/include/metis.h"
#endif


int FuseSolutionElmerPartitioned(char *prefix,char *outfile,int decimals,int parts,
				 int minstep, int maxstep, int dstep, int info)
{
#define LONGLINE 2048
  int *noknots,*noelements,novctrs,elemcode;
  int totknots,totelements,sumknots,sumelements;
  int timesteps,i,j,k,l,step;
  int ind[MAXNODESD3];
  int nofiles,activestep;
  Real *res, x, y, z;
  FILE *in[MAXPARTITIONS+1],*intest,*out;
  char line[LONGLINE],filename[MAXFILESIZE],text[MAXNAMESIZE],outstyle[MAXFILESIZE];
  char *cp, *charend;

  if(minstep || maxstep || dstep) {
    if(info) printf("Saving results in the interval from %d to %d with step %d\n",minstep,maxstep,dstep);
  }

  for(i=0;;i++) {
    sprintf(filename,"%s.ep.%d",prefix,i);
    if ((intest = fopen(filename,"r")) == NULL) break;
    if(i<=MAXPARTITIONS) in[i] = intest;
  }
  if(i > MAXPARTITIONS) {
    printf("**********************************************************\n");
    printf("Only data for %d partitions is fused (%d)\n",MAXPARTITIONS,i);
    printf("**********************************************************\n");
    i = MAXPARTITIONS;
  } 
  nofiles = i;

  if(nofiles < 2) {
    printf("Opening of partitioned data from file %s wasn't successful!\n",
	   filename);
    return(2);
  } else {
    if(parts > 0) nofiles = MIN(parts,nofiles);
    if(info) printf("Loading Elmer results from %d partitions.\n",nofiles);
  }

  if(minstep || maxstep || dstep) {
    if(info) printf("Saving results in the interval from %d to %d with step %d\n",minstep,maxstep,dstep);
  }

  noknots = Ivector(0,nofiles-1);
  noelements = Ivector(0,nofiles-1);
 
  sumknots = 0;
  sumelements = 0;

  for(i=0;i<nofiles;i++) {
    charend = fgets(line,LONGLINE,in[i]);
    if(i==0) {
      cp = line;
      noknots[i] = next_int(&cp);
      noelements[i] = next_int(&cp);
      novctrs = next_int(&cp);
      timesteps = next_int(&cp);
    }
    else {
      sscanf(line,"%d %d",&noknots[i],&noelements[i]);
    }
    sumknots += noknots[i];
    sumelements += noelements[i];
  }
  totknots = sumknots;
  totelements = sumelements;
  res = Rvector(1,novctrs);

  if(info) printf("There are altogether %d nodes and %d elements.\n",totknots,sumelements);


  AddExtension(outfile,filename,"ep");
  if(info) printf("Saving ElmerPost data to %s.\n",filename);  
  out = fopen(filename,"w");
  if(out == NULL) {
    printf("opening of file was not successful\n");
    return(3);
  }

  i = timesteps;
  if(minstep || maxstep || dstep) {
    if ( dstep>1 ) {
      i = 0;
      for(step = minstep; step <= maxstep; step++)
        if((step-minstep)%dstep==0) i++;
    } else i=maxstep-minstep+1;
  }
  fprintf(out,"%d %d %d %d %s %s",totknots,totelements,novctrs+1,i,"scalar: Partition",cp);
 
  if(info) printf("Reading and writing %d coordinates.\n",totknots);
  sprintf(outstyle,"%%.%dg %%.%dg %%.%dg\n",decimals,decimals,decimals);

  for(j=0; j < nofiles; j++) {
    for(i=1; i <= noknots[j]; i++) {
      do {
	charend = fgets(line,LONGLINE,in[j]);
      } while(line[0] == '#');

      sscanf(line,"%le %le %le",&x,&y,&z);
      fprintf(out,outstyle,x,y,z);
    }
  }

  if(info) printf("Reading and writing %d element topologies.\n",totelements);
  sumknots = 0;

  for(j=0; j < nofiles; j++) {
    for(i=1; i <= noelements[j]; i++) {
      do {
	charend = fgets(line,LONGLINE,in[j]);
      } while (line[0] == '#');

      sscanf(line,"%s",text);
      cp = strstr(line," ");

      elemcode = next_int(&cp);

      for(k=0;k< elemcode%100 ;k++) {
	/* Dirty trick for long lines */
	l = strspn(cp," ");
	if( l == 0) {
	  charend = fgets(line,LONGLINE,in[j]);
	  cp = line;
	}
	ind[k] = next_int(&cp);
      }
      if(elemcode == 102) elemcode = 101;

      fprintf(out,"%s %d",text,elemcode);
      for(k=0;k < elemcode%100 ;k++)       
	fprintf(out," %d",ind[k]+sumknots);
      fprintf(out,"\n");
    }
    sumknots += noknots[j];
  }

  if(info) printf("Reading and writing %d degrees of freedom.\n",novctrs);
  sprintf(outstyle,"%%.%dg ",decimals);

  activestep = FALSE;
  if(maxstep) timesteps = MIN(timesteps, maxstep);

  for(step = 1; step <= timesteps; step++) {
        
    if (step>=minstep) {
      if ( dstep>0 ) {
         activestep=((step-minstep)%dstep==0);
      } else activestep=TRUE;
    }

    for(k=0;k<nofiles;k++) 
      for(i=1; i <= noknots[k]; i++) {
	do {
	  charend = fgets(line,LONGLINE,in[k]);
          if (activestep) {
            if(k==0 && strstr(line,"#time")) {
	      fprintf(out,"%s",line);
	      fprintf(stderr,"%s",line);
            }
          }
	}
	while (line[0] == '#');

	if(activestep) {
	  cp = line;
	  for(j=1;j <= novctrs;j++) 
	    res[j] = next_real(&cp);
	  
	  fprintf(out,"%d ",k+1);
	  for(j=1;j <= novctrs;j++) 
	    fprintf(out,outstyle,res[j]);
	  fprintf(out,"\n");
	}

      }
  }


  for(i=0;i<nofiles;i++) 
    fclose(in[i]);
  fclose(out);

  if(info) printf("Successfully fused partitioned Elmer results\n");

  return(0);
}




static int CreatePartitionTable(struct FemType *data,int info)
{
  int i,j,k,m,noelements,noknots,partitions,nonodes,periodic;
  int maxneededtimes,sharings,part,ind,hit,notinany,debug;
  int *indxper;

  printf("Creating a table showing all parenting partitions of nodes.\n");  

  if(data->maxpartitiontable) {
    printf("The partition table already exists!\n");
    smallerror("Partition table not done");
    return(0);
  }

  partitions = data->nopartitions;
  if( partitions <= 1 ) {
    bigerror("Cannot do partition table for less than two partitions");
  }
  

  noelements = data->noelements;
  periodic = data->periodicexist;
  noknots = data->noknots;
  if(periodic) indxper = data->periodic;

  maxneededtimes = 0;
  sharings = 0;

  /* Make the partition list taking into account periodicity */
  for(i=1;i<=noelements;i++) {
    part = data->elempart[i];
    nonodes =  data->elementtypes[i] % 100;

    for(j=0;j < nonodes;j++) {
      ind = data->topology[i][j];
      if(periodic) ind = indxper[ind];

      debug = FALSE;
      if(debug) printf("ind=%d i=%d j=%d part=%d\n",ind,i,j,part); 

      hit = 0;
      for(k=1;k<=maxneededtimes;k++) { 
	if(data->partitiontable[k][ind] == 0) hit = k;
	if(data->partitiontable[k][ind] == part) hit = -k;
	if(hit) break;
      }
      if( hit > 0) {
	data->partitiontable[hit][ind] = part;
	if(hit == 2) sharings++;
      }
      else if(hit == 0) {
	maxneededtimes++;
	data->partitiontable[maxneededtimes] = Ivector(1,noknots);
	for(m=1;m<=noknots;m++)
	  data->partitiontable[maxneededtimes][m] = 0;
	data->partitiontable[maxneededtimes][ind] = part;
	if(maxneededtimes == 2) sharings++;
      }
      if(debug) printf("hit = %d\n",hit);
    }
  }


  /* Make the partitiontable such that the owner node is the first one in the list */  
  notinany = 0;
  for(i=1;i<=noknots;i++) {
    
    /* Skip the periodic nodes and take care of them later */
    if(periodic) 
      if(i != indxper[i]) continue;

    hit = FALSE;
    for(k=1;k<=maxneededtimes;k++) {
      if(!data->partitiontable[k][i]) break;
      if(data->partitiontable[k][i] == data->nodepart[i]) {
	hit = k;
	break;
      }
    }
    if( hit > 1 ) {
      data->partitiontable[hit][i] = data->partitiontable[1][i];
      data->partitiontable[1][i] = data->nodepart[i];
    }
    else if(!hit) {
      if(0) {
	printf("Node %d in partition %d not in the table!\n",i,data->nodepart[i]);
	if(periodic) printf("indexper: %d\n",indxper[i]);
      }
      
      notinany++;
      data->nodepart[i] = data->partitiontable[1][i];
    }
  }


  /* For periodic counterparts copy the table and ownership */
  if(periodic) {
    for(i=1;i<=noknots;i++) {
      ind = indxper[i];
      if(ind == i) continue;
      for(k=1;k<=maxneededtimes;k++) {
	if( !data->partitiontable[k][ind] ) break;
	data->partitiontable[k][i] = data->partitiontable[k][ind];
      }
      data->nodepart[i] = data->nodepart[ind];
    }
  }


  if(info) {
    printf("Nodes belong to %d partitions in maximum\n",maxneededtimes);
    printf("There are %d shared nodes which is %.2f %% of all nodes.\n",
	   sharings,(100.*sharings)/noknots);
    printf("The initial owner was not any of the elements for %d nodes\n",notinany);
  }


  if(0) for(i=1;i<=noknots;i++) {
    if(data->partitiontable[maxneededtimes][i]) {
      printf("node %d parts: ",i);
      for(k=1;k<=maxneededtimes;k++) 
	printf("%d ",data->partitiontable[k][i]);
      printf("\n");
    }
  }
    
  data->maxpartitiontable = maxneededtimes;
  data->partitiontableexists = TRUE;
  return(0);
}



static int PartitionElementsByNodes(struct FemType *data,int info)
/* Given nodal partitioning determine the elemental partitioning. 
   This is usually suboptimal, it is preferable to do the other way around. */
{
  int i,j,noelements,nopartitions,part,maxpart,maxpart2,minpart;
  int *elempart,*nodepart,*nodesinpart,*cuminpart,**knows,**cumknows,set;

  if(!data->partitionexist) return(1);

  noelements = data->noelements;
  nopartitions = data->nopartitions;
  elempart = data->elempart;
  nodepart = data->nodepart;

  nodesinpart = Ivector(1,nopartitions);
  cuminpart = Ivector(1,nopartitions);
  for(j=1;j<=nopartitions;j++) 
    cuminpart[j] = 0;

  knows = Imatrix(1,nopartitions,1,nopartitions);
  cumknows = Imatrix(1,nopartitions,1,nopartitions);
  for(i=1;i<=nopartitions;i++)
    for(j=1;j<=nopartitions;j++)
      knows[i][j] = cumknows[i][j] = 0;

  set = FALSE;

 omstart:

  /* In the first round count the equally joined elements and 
     on the second round split them equally using cumulative numbering */

  for(i=1;i<=noelements;i++) {
    for(j=1;j<=nopartitions;j++) 
      nodesinpart[j] = 0;
    for(j=0;j<data->elementtypes[i] % 100;j++) {
      part = nodepart[data->topology[i][j]];
      nodesinpart[part] += 1;
    }
    maxpart = maxpart2 = 1;
    for(j=1;j<=nopartitions;j++) 
      if(nodesinpart[j] > nodesinpart[maxpart]) maxpart = j;
    if(maxpart == 1) maxpart2 = 2;
    for(j=1;j<=nopartitions;j++) {
      if(j == maxpart) continue;
      if(nodesinpart[j] > nodesinpart[maxpart2]) maxpart2 = j;
    }
    
    if(nodesinpart[maxpart] > nodesinpart[maxpart2]) {
      if(set) 
	elempart[i] = maxpart;    
      else
	cuminpart[maxpart] += 1;
    }
    else {
      if(set) {
	cumknows[maxpart][maxpart2] += 1;
	if( cumknows[maxpart][maxpart2] > knows[maxpart][maxpart2] / 2) {
	  elempart[i] = maxpart2;
	  cuminpart[maxpart2] += 1;
	}
	else {
	  elempart[i] = maxpart;
	  cuminpart[maxpart] += 1;
	}
      }	
      else
	knows[maxpart][maxpart2] += 1;
    }
  }    

  if(!set) {
    set = TRUE;
    goto omstart;
  }

  minpart = maxpart = cuminpart[1];
  for(j=1;j<=nopartitions;j++) {
    minpart = MIN( minpart, cuminpart[j]);
    maxpart = MAX( maxpart, cuminpart[j]);
  }

  if(info) {
    printf("Set the element partitions by the dominating nodal partition\n");
    printf("There are from %d to %d elements in the %d partitions.\n",minpart,maxpart,nopartitions);
  }  

  free_Ivector(nodesinpart,1,nopartitions);
  free_Ivector(cuminpart,1,nopartitions);
  free_Imatrix(knows,1,nopartitions,1,nopartitions);
  free_Imatrix(cumknows,1,nopartitions,1,nopartitions);

  return(0);
}


static int PartitionNodesByElements(struct FemType *data,int info)
/* Given the elemental partitioning set the nodal ownership.
   This is optimal for Elmer since the elemental partitioning is primary. */
{
  int i,j,k,noknots,nopartitions,part,minpart,maxpart;
  int maxpart2,*cuminpart,**knows,**cumknows,set;
  int *elempart,*nodepart,*nodesinpart;
  int *invrow,*invcol;

  if(!data->partitionexist) return(1);

  CreateInverseTopology(data, info);

  invrow = data->invtopo.rows;
  invcol = data->invtopo.cols;

  noknots = data->noknots;
  nopartitions = data->nopartitions;
  elempart = data->elempart;
  nodepart = data->nodepart;

  if(info) printf("Number of nodal partitions: %d\n",nopartitions);

  nodesinpart = Ivector(1,nopartitions);
  cuminpart = Ivector(1,nopartitions);
  for(j=1;j<=nopartitions;j++) 
    cuminpart[j] = 0;

  knows = Imatrix(1,nopartitions,1,nopartitions);
  cumknows = Imatrix(1,nopartitions,1,nopartitions);
  for(i=1;i<=nopartitions;i++)
    for(j=1;j<=nopartitions;j++)
      knows[i][j] = cumknows[i][j] = 0;

  set = FALSE;

 omstart:

  for(i=1;i<=noknots;i++) {

    for(j=1;j<=nopartitions;j++) 
      nodesinpart[j] = 0;

    /* Tag the number of owner partitions */
    for(j=invrow[i-1];j<invrow[i];j++) {
      k = invcol[j]+1;
      part = elempart[k];
      if( part > 0 ) nodesinpart[part] += 1;
    }

    /* Find the partition with maximum number of hits */
    maxpart = maxpart2 = 0;
    for(j=1;j<=nopartitions;j++) 
      if(nodesinpart[j] > 0 ) {
	if(!maxpart) {
	  maxpart = j;
	} 
	else if(nodesinpart[j] > nodesinpart[maxpart]) {
	  maxpart = j;
	}
      }
  
    /* Find the partition with 2nd largest number of hits */
    for(j=1;j<=nopartitions;j++) 
      if(nodesinpart[j] > 0 ) {
	if(j == maxpart) continue;
	if(!maxpart2) {
	  maxpart2 = j;
	} 
	else if(nodesinpart[j] > nodesinpart[maxpart2]) {
	  maxpart2 = j;
	}
      }

    /* If there is a clear dominator use that */
    if(!maxpart && !maxpart2)
      printf("Node is not included in any partition: %d\n",i);
    else if(!maxpart2) {
      nodepart[i] = maxpart;
      cuminpart[maxpart] += 1;
    }
     else
      if(nodesinpart[maxpart] > nodesinpart[maxpart2]) {
	if(set) 
	  nodepart[i] = maxpart;    
	else
	  cuminpart[maxpart] += 1;
      }

    /* Otherwise make a half and half split between the major owners */
      else {
	if(set) {
	  cumknows[maxpart][maxpart2] += 1;
	  if( cumknows[maxpart][maxpart2] > knows[maxpart][maxpart2] / 2) {
	    nodepart[i] = maxpart2;
	    cuminpart[maxpart2] += 1;
	  }
	  else {
	    nodepart[i] = maxpart;
	    cuminpart[maxpart] += 1;
	  }
	}	
	else {
	  knows[maxpart][maxpart2] += 1;
	}
      }
  }
  
  if(!set) {
    set = TRUE;
    goto omstart;
  }

  minpart = maxpart = cuminpart[1];
  for(j=1;j<=nopartitions;j++) {
    minpart = MIN( minpart, cuminpart[j]);
    maxpart = MAX( maxpart, cuminpart[j]);
  }

  if(info) {
    printf("Set the node partitions by the dominating element partition.\n");
    printf("There are from %d to %d nodes in the %d partitions.\n",minpart,maxpart,nopartitions);
  }  

  free_Ivector(nodesinpart,1,nopartitions);
  free_Ivector(cuminpart,1,nopartitions);
  free_Imatrix(knows,1,nopartitions,1,nopartitions);
  free_Imatrix(cumknows,1,nopartitions,1,nopartitions);

  return(0);
}



int PartitionSimpleElements(struct FemType *data,struct ElmergridType *eg,struct BoundaryType *bound,
			    int dimpart[],int dimper[],int partorder, Real corder[],
			    Real parttol, int info)
/* Partition elements recursively in major directions. 
   This may be the optimal method of partitioning for simple geometries. */ 
{
  int i,j,k,l,kprev,ind,minpart,maxpart;
  int noknots,noelements,nonodes,elemsinpart,periodic;
  int partitions1,partitions2,partitions3,partitions;
  int vpartitions1,vpartitions2,vpartitions3,vpartitions;
  int noelements0,noelements1,noparts0;
  int *indx,*nopart,*inpart,*elemconnect;
  Real *arrange;
  Real x,y,z,cx,cy,cz;

  if(info) printf("PartitionSimpleElements\n");

  noelements = data->noelements;
  noknots = data->noknots;
  
  partitions1 = dimpart[0];
  partitions2 = dimpart[1];
  partitions3 = dimpart[2];
  if(data->dim < 3) partitions3 = 1;
  partitions = partitions1 * partitions2 * partitions3;

  if(partitions1 < 2 && partitions2 < 2 && partitions3 < 2 && !eg->connect) {
    printf("Some of the divisions must be larger than one: %d %d %d\n",
	   partitions1, partitions2, partitions3 );
    bigerror("Partitioning not performed");
  }

  if(partitions >= noelements) {
    printf("There must be fever partitions than elements (%d vs %d)!\n",
	   partitions,noelements);
    bigerror("Partitioning not performed");
  }
    
  if( eg->partbcz > 1 || eg->partbcr ) 
    PartitionConnectedElements1D(data,bound,eg,info);
  else if( eg->partbcmetis > 1 ) 
    PartitionConnectedElementsMetis(data,bound,eg,eg->partbcmetis,3,info); 
  else if( eg->connect ) 
   PartitionConnectedElementsStraight(data,bound,eg,info);
     

  if( data->nodeconnectexist ) 
    ExtendBoundaryPartitioning(data,bound,eg->partbclayers,info);

  if(!data->partitionexist) {
    data->partitionexist = TRUE;
    data->elempart = Ivector(1,noelements);
    data->nodepart = Ivector(1,noknots);
    data->nopartitions = partitions;
  }
  inpart = data->elempart;
  
  vpartitions1 = partitions1;
  vpartitions2 = partitions2;
  vpartitions3 = partitions3;

  if(0) printf("connect: %d %d\n",data->elemconnectexist,data->nodeconnectexist);

  noparts0 = 0;
  noelements0 = 0;
  if( data->elemconnectexist || data->nodeconnectexist) {
    elemconnect = data->elemconnect;  
    noparts0 = 0;
    for(i=1;i<=data->noelements;i++) {
      if( elemconnect[i] ) noelements0 = noelements0 + 1;
      noparts0 = MAX( noparts0, elemconnect[i] );
    }
    if(info) printf("There are %d initial partitions in the connected mesh\n",noparts0);
    if(info) printf("There are %d initial elements in the connected mesh\n",noelements0);
  }
  noelements1 = noelements - noelements0; 

  vpartitions1 = partitions1;
  vpartitions2 = partitions2;
  vpartitions3 = partitions3;

  periodic = dimper[0] || dimper[1] || dimper[2];
  if(periodic) {
    if(dimper[0] && partitions1 > 1) vpartitions1 *= 2;
    if(dimper[1] && partitions2 > 1) vpartitions2 *= 2;
    if(dimper[2] && partitions3 > 1) vpartitions3 *= 2;
  }
  vpartitions = vpartitions1 * vpartitions2 * vpartitions3;
  nopart = Ivector(1,vpartitions+noparts0);
  
  if( vpartitions == 1 && data->elemconnectexist ) {
    if(info) printf("Only one regular partitions requested, skipping 2nd part of hybrid partitioning\n");
    for(j=1;j<=noelements;j++) {
      if( elemconnect[j] > 0 ) 
	inpart[j] = elemconnect[j];
      else
	inpart[j] = noparts0 + 1;
    }
    partitions = noparts0 + 1;
    data->nopartitions = partitions;
    goto skippart;
  }



  if(info) printf("Making a simple partitioning for %d elements in %d-dimensions.\n",
		  noelements1,data->dim);

  arrange = Rvector(1,noelements);
  indx = Ivector(1,noelements);

  if(partorder) {
    cx = corder[0];
    cy = corder[1];
    cz = corder[2];    
  }
  else {
    cx = 1.0;
    cy = 0.0001;
    cz = cy*cy;
  }
  z = 0.0;

  for(i=1;i<=noelements;i++) 
    inpart[i] = 1;
  
  if(vpartitions1 > 1) {

    if(info) printf("Ordering 1st direction with (%.3lg*x + %.3lg*y + %.3lg*z)\n",cx,cy,cz);

    for(j=1;j<=noelements;j++) {
      if( data->elemconnectexist ) {
	if( elemconnect[j] > 0 ) {
	  arrange[j] = 1.0e9;
	  continue;
	}
      }
      nonodes = data->elementtypes[j]%100;
      x = y = z = 0.0;
      for(i=0;i<nonodes;i++) {
	k = data->topology[j][i];
	x += data->x[k];
	y += data->y[k];
	z += data->z[k];
      }
      arrange[j] = (cx*x + cy*y + cz*z) / nonodes;
    }

    SortIndex(noelements,arrange,indx);
    
    for(i=1;i<=noelements1;i++) {
      ind = indx[i];
      k = (i*vpartitions1-1)/noelements1+1;
      inpart[ind] = k;
    }
  } 

 
  /* Partition in the 2nd direction taking into account the 1st direction */
  if(vpartitions2 > 1) {
    if(info) printf("Ordering in the 2nd direction.\n");

    for(j=1;j<=noelements;j++) {
      if( data->elemconnectexist ) {
	if( elemconnect[j] > 0 ) {
	  arrange[j] = 1.0e9;
	  continue;
	}
      }
      nonodes = data->elementtypes[j]%100;
      x = y = z = 0.0;
      for(i=0;i<nonodes;i++) {
	k = data->topology[j][i];
	x += data->x[k];
	y += data->y[k];
	z += data->z[k];
      }
      arrange[j] = (-cy*x + cx*y + cz*z) / nonodes;
    }
    SortIndex(noelements,arrange,indx);
    
    for(i=1;i<=vpartitions;i++)
      nopart[i] = 0;
    
    elemsinpart = noelements1 / (vpartitions1*vpartitions2);
    for(i=1;i<=noelements1;i++) {
      j = 0;
      ind = indx[i];
      do {
	j++;
	k = (inpart[ind]-1) * vpartitions2 + j;
      }
      while(nopart[k] >= elemsinpart && j < vpartitions2);
      
      nopart[k] += 1;
      inpart[ind] = (inpart[ind]-1)*vpartitions2 + j;
    }
  }  

  /* Partition in the 3rd direction taking into account the 1st and 2nd direction */
  if(vpartitions3 > 1) {
    if(info) printf("Ordering in the 3rd direction.\n");

    for(j=1;j<=noelements;j++) {
      if( data->elemconnectexist ) {
	if( elemconnect[j] > 0 ) {
	  arrange[j] = 1.0e9;
	  continue;
	}
      }

      nonodes = data->elementtypes[j]%100;
      x = y = z = 0.0;
      for(i=0;i<nonodes;i++) {
	k = data->topology[j][i];
	x += data->x[k];
	y += data->y[k];
	z += data->z[k];
      }
      arrange[j] = (-cz*x - cy*y + cx*z) / nonodes;
    }

    SortIndex(noelements,arrange,indx);

    for(i=1;i<=vpartitions;i++)
      nopart[i] = 0;
    
    elemsinpart = noelements1 / (vpartitions1*vpartitions2*vpartitions3);
    for(i=1;i<=noelements1;i++) {
      j = 0;
      ind = indx[i];
      do {
	j++;
	k = (inpart[ind]-1)*vpartitions3 + j;
      }
      while(nopart[k] >= elemsinpart && j < vpartitions3);
    
      nopart[k] += 1;
      inpart[ind] = (inpart[ind]-1)*vpartitions3 + j;
    }
  }

  /* This piece of code may be used to assign "strides" back to same partition.
     After performing the geometric partitioning we revisit whether elements 
     are too close to each other and yet in different partition in the chosen measure.
     For example, if chosen measure has the smallest coefficient in z-direction then
     a suitable value will revert elements to being on the same stride. */
  if(parttol > 1.0e-20 )  {
    if(0) {
      printf("Original partition counts:\n");
      for(i=1;i<=partitions;i++)
	printf("nopart[%d] = %d\n",i,nopart[i]);
    }
      
    for(j=1;j<=noelements;j++) {
      if( data->elemconnectexist ) {
	if( elemconnect[j] > 0 ) {
	  arrange[j] = 1.0e9;
	  continue;
	}
      }
      nonodes = data->elementtypes[j]%100;
      x = y = z = 0.0;
      for(i=0;i<nonodes;i++) {
	k = data->topology[j][i];
	x += data->x[k];
	y += data->y[k];
	z += data->z[k];
      }
      arrange[j] = (cx*x + cy*y + cz*z) / nonodes;
    }
    SortIndex(noelements,arrange,indx);

    kprev = inpart[indx[1]];
    l = 0;
    for(i=2;i<=noelements1;i++) {
      ind = indx[i];
      k = inpart[ind];
      if(k != kprev ) {
	if(fabs(arrange[i]-arrange[i-1]) < parttol ) {
	  nopart[k] -= 1;
	  nopart[kprev] += 1;
	  inpart[ind] = kprev;
	  k = kprev;
	  if(0) printf("arrange %d %d %d %d %.6lg %.6lg\n",i,ind,k,kprev,arrange[i],arrange[i-1]);
	  l++;
	}
      }
      kprev = k;
    }
    if(0) {
      printf("Updated partition counts:\n");
      for(i=1;i<=partitions;i++)
	printf("nopart[%d] = %d\n",i,nopart[i]);    
    }      
    if(info) printf("Number of partition changes due to tolerance: %d\n",l);
  }

  
  free_Rvector(arrange,1,noelements);
  free_Ivector(indx,1,noelements);


  if( data->elemconnectexist ) {
    for(j=1;j<=noelements;j++) {
      if( elemconnect[j] > 0 ) 
	inpart[j] = elemconnect[j];
      else
	inpart[j] = inpart[j] + noparts0;
    }
    partitions = partitions + noparts0;
    data->nopartitions = partitions;
  }


  /* For periodic systems the number of virtual partitions is larger. Now map the mesh so that the 
     1st and last partition for each direction will be joined */
  if(periodic) {
    int *partmap;
    int p1,p2,p3,q1,q2,q3;
    int P,Q;

    if(data->elemconnectexist ) {
      bigerror("Cannot use connect flag with periodic systems\n");
    }

    p1=p2=p3=1;
    partmap = Ivector(1,vpartitions);
    for(i=1;i<=vpartitions;i++)
      partmap[i] = 0;
    for(p1=1;p1<=vpartitions1;p1++) {
      q1 = p1;
      if(dimper[0] && vpartitions1 > 1) {
	if(q1==vpartitions1) q1 = 0;
	q1 = q1/2 + 1;
      }
      for(p2=1;p2<=vpartitions2;p2++) {
	q2 = p2;
	if(dimper[1] && vpartitions2 > 1) {
	  if(q2==vpartitions2) q2 = 0;
	  q2 = q2/2 + 1;
	}
	for(p3=1;p3<=vpartitions3;p3++) {
	  q3 = p3;
	  if(dimper[2] && vpartitions3 > 1) {
	    if(q3==vpartitions3) q3 = 0;
	    q3 = q3/2 + 1;
	  }
	  
	  P = vpartitions3 * vpartitions2 * (p1 - 1) + vpartitions3 * (p2-1) + p3;
	  Q = partitions3 * partitions2 * (q1 - 1) + partitions3 * (q2-1) + q3;

	  partmap[P] = Q;
	}
      }
    }
    for(i=1;i<=noelements;i++)
      inpart[i] = partmap[inpart[i]];
    free_Ivector(partmap,1,vpartitions);
  }

 skippart:

  for(i=1;i<=partitions;i++)
    nopart[i] = 0;
  for(i=1;i<=noelements;i++) 
    nopart[inpart[i]] += 1;

  minpart = maxpart = nopart[1];
  for(i=1;i<=partitions;i++) {
    minpart = MIN( nopart[i], minpart );
    maxpart = MAX( nopart[i], maxpart );
  }

  free_Ivector(nopart,1,vpartitions);
  PartitionNodesByElements(data,info);

  if(info) printf("Successfully made a partitioning with %d to %d elements.\n",minpart,maxpart);

  return(0);
}


int PartitionSimpleElementsNonRecursive(struct FemType *data,int dimpart[],int dimper[],
			    int info)
/* This is similar to the previous except the partitioning is not recursive. 
   This kind of partitioning is optimal for uniform grids with blockwise structure.
   An example would be L-type geometry where 3/4 of potential partitions would be active. */
{
  int i,j,k,minpart,maxpart;
  int noknots, noelements,nonodes,periodic;
  int partitions1,partitions2,partitions3,partitions;
  int IndX,IndY,IndZ;
  int *nopart,*inpart;
  Real x,y,z,MaxX,MinX,MaxY,MinY,MaxZ,MinZ;

  partitions1 = dimpart[0];
  partitions2 = dimpart[1];
  partitions3 = dimpart[2];
  if(data->dim < 3) partitions3 = 1;
  partitions = partitions1 * partitions2 * partitions3;

  if(partitions1 < 2 && partitions2 < 2 && partitions3 < 2) {
    printf("Some of the divisions must be larger than one: %d %d %d\n",
	   partitions1, partitions2, partitions3 );
    bigerror("Partitioning not performed");
  }
    
  if(!data->partitionexist) {
    data->partitionexist = TRUE;
    data->elempart = Ivector(1,data->noelements);
    data->nodepart = Ivector(1,data->noknots);
    data->nopartitions = partitions;
  }
  inpart = data->elempart;
  noelements = data->noelements;
  noknots = data->noknots;

  periodic = dimper[0] || dimper[1] || dimper[2];
  if(periodic) {
    printf("Implement periodicity for this partitioning routine!\n");
  }

  if(info) {
    printf("Making a simple partitioning for %d elements in %d-dimensions.\n",
	   noelements,data->dim);
    printf("There can be at maximum %d partitions\n",partitions);
  }

  nopart = Ivector(1,partitions);
  for(i=1;i<=partitions;i++)
    nopart[i] = 0;

  z = 0.0;
  IndZ = 1;
  for(i=1;i<=noelements;i++) 
    inpart[i] = 0;

  MaxX = MinX = data->x[1];
  MaxY = MinY = data->y[1];
  MaxZ = MinZ = data->z[1];

  for(i=1;i<=noknots;i++) {
    x = data->x[i];
    y = data->y[i];
    z = data->z[i];
    
    MaxX = MAX( MaxX, x);
    MinX = MIN( MinX, x);
    MaxY = MAX( MaxY, y);
    MinY = MIN( MinY, y);
    MaxZ = MAX( MaxZ, z);
    MinZ = MIN( MinZ, z);
  }
  if( info ) {
    printf("Range in x-direction: %12.5le %12.5le\n",MinX,MaxX);
    printf("Range in y-direction: %12.5le %12.5le\n",MinY,MaxY);
    printf("Range in z-direction: %12.5le %12.5le\n",MinZ,MaxZ);
  }
   
  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j]%100;
    x = y = z = 0.0;
    for(i=0;i<nonodes;i++) {
      k = data->topology[j][i];
      x += data->x[k];
      y += data->y[k];
      z += data->z[k];
    }
    x = x / nonodes;
    y = y / nonodes;
    z = z / nonodes;

    IndX = ceil( partitions1 * ( x - MinX ) / ( MaxX - MinX ) );
    IndY = ceil( partitions2 * ( y - MinY ) / ( MaxY - MinY ) );
    if( data->dim == 3 ) {
      IndZ = ceil( partitions3 * ( z - MinZ ) / ( MaxZ - MinZ ) );
    }
    k = (IndZ-1) * partitions1 * partitions2 + 
      (IndY-1) * partitions1 + IndX;

    inpart[j] = k;
    nopart[k] += 1;
  } 

  maxpart = 0;
  minpart = noelements;
  j = 0;
  if( info ) printf("Renumbering the partitions (showing only 64):\n");
  for(i=1;i<=partitions;i++) {
    k = nopart[i];
    if( k ) {
      maxpart = MAX( k, maxpart );
      minpart = MIN( k, minpart );
      j += 1;
      nopart[i] = j;
      if( info && j<=64) printf("%d -> %d\n",i,j);      
    }
  }
  if(info) printf("There are %d active partitions out of %d possible ones\n",j,partitions);

  data->nopartitions = j;

  for(i=1;i<=noelements;i++) {
    j = inpart[i];
    inpart[i] = nopart[j];
  }

  free_Ivector(nopart,1,partitions);

  PartitionNodesByElements(data,info);

  if(info) printf("Successfully made a partitioning with %d to %d elements.\n",minpart,maxpart);

  return(0);
}


#define MAXCATEGORY 500
int PartitionSimpleElementsRotational(struct FemType *data,int dimpart[],int dimper[],
				      int info) {
  int i,j,k,minpart,maxpart,dim,hit;
  int noknots, noelements,nonodes,periodic;
  int partitions1,partitions2,partitions3,partitions;
  int IndR,IndF,IndZ,connect;
  int *nopart,*inpart,*cumf,*cumr,*cumz = NULL,*elemconnect = NULL;
  Real x,y,z,r,f,MaxR,MinR,MaxF,MinF,MaxZ,MinZ;

  dim = data->dim;
  partitions1 = dimpart[0];
  partitions2 = dimpart[1];
  partitions3 = dimpart[2];
  if(dim < 3) partitions3 = 1;
  partitions = partitions1 * partitions2 * partitions3;

  /* See if there are connected elements */
  connect = FALSE;
  if(data->nodeconnectexist) {
    SetConnectedElements(data,info);
    connect = TRUE;
    elemconnect = data->elemconnect;
    partitions = partitions + partitions3;
  }
    
  if(partitions1 < 2 && partitions2 < 2 && partitions3 < 2) {
    printf("Some of the divisions must be larger than one: %d %d %d\n",
	   partitions1, partitions2, partitions3 );
    bigerror("Partitioning not performed");
  }

  if(!data->partitionexist) {
    data->partitionexist = TRUE;
    data->elempart = Ivector(1,data->noelements);
    data->nodepart = Ivector(1,data->noknots);
    data->nopartitions = partitions;
  }
  inpart = data->elempart;
  noelements = data->noelements;
  noknots = data->noknots;

  periodic = dimper[0] || dimper[1] || dimper[2];
  if(periodic) {
    printf("Implement periodicity for this partitioning routine!\n");
  }

  if(info) {
    printf("Making a simple rotational partitioning for %d elements in %d-dimensions.\n",
	   noelements,data->dim);
    printf("There can be at maximum %d partitions\n",partitions);
  }

  nopart = Ivector(1,partitions);
  for(i=1;i<=partitions;i++)
    nopart[i] = 0;

  z = 0.0;
  IndZ = 1;
  for(i=1;i<=noelements;i++) 
    inpart[i] = 0;

  x = data->x[1];
  y = data->y[1];
  z = data->z[1];

  r = sqrt(x*x+y*y);
  f = 180 * atan2(y,x)/FM_PI;
  if( f < 0.0 ) f = f + 360.0;    
  MaxR = MinR = r;
  MaxF = MinF = f;
  MaxZ = MinZ = z;

  for(i=1;i<=noknots;i++) {
    x = data->x[i];
    y = data->y[i];
    z = data->z[i];

    r = sqrt(x*x+y*y);
    f = 180 * atan2(y,x)/FM_PI;
    if( f < 0.0 ) f = f + 360.0;    

    MaxR = MAX( MaxR, r);
    MinR = MIN( MinR, r);
    MaxF = MAX( MaxF, f);
    MinF = MIN( MinF, f);
    MaxZ = MAX( MaxZ, z);
    MinZ = MIN( MinZ, z);
  }

  if( info ) {
    printf("Range in r-direction: %12.5le %12.5le\n",MinR,MaxR);
    printf("Range in f-direction: %12.5le %12.5le\n",MinF,MaxF);
    printf("Range in z-direction: %12.5le %12.5le\n",MinZ,MaxZ);
  }
  if( MaxF - MinF > 180.0 ) {
    MaxF = 360.0;
    MinF = 0.0;
  }

  /* Zero is the 1st value so that recursive algos can be used. */ 
  cumr = Ivector(0,MAXCATEGORY);
  cumf = Ivector(0,MAXCATEGORY);
  if(dim==3) cumz = Ivector(0,MAXCATEGORY);

  for(i=0;i<=MAXCATEGORY;i++) {
    cumf[i] = cumr[i] = 0;
    if( dim == 3) cumz[i] = 0;
  }   

  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j]%100;
    x = y = z = 0.0;
    for(i=0;i<nonodes;i++) {
      k = data->topology[j][i];
      x += data->x[k];
      y += data->y[k];
      z += data->z[k];
    }
    x = x / nonodes;
    y = y / nonodes;
    z = z / nonodes;
    
    r = sqrt(x*x+y*y);
    f = 180 * atan2(y,x)/FM_PI;
    if( f < 0.0 ) f = f + 360.0;    
    
    IndR = ceil( MAXCATEGORY * ( r - MinR ) / ( MaxR - MinR ) );
    IndF = ceil( MAXCATEGORY * ( f - MinF ) / ( MaxF - MinF ) );
    if(dim==3) IndZ = ceil( MAXCATEGORY * ( z - MinZ ) / ( MaxZ - MinZ ) );
    
    if( IndR < 1 || IndR > MAXCATEGORY ) {
      printf("IndR out of bounds : %d\n",IndR );
      IndR = MIN( MAX( IndR, 1 ), MAXCATEGORY );
    }
    if( IndF < 1 || IndF > MAXCATEGORY ) {
      printf("IndF out of bounds : %d\n",IndF );
      IndF = MIN( MAX( IndF, 1 ), MAXCATEGORY );
    }
    if( dim == 3 ) {
      if( IndZ < 1 || IndZ > MAXCATEGORY ) {
	printf("IndZ out of bounds : %d\n",IndZ );
	IndZ = MIN( MAX( IndZ, 1 ), MAXCATEGORY );
      }
      cumz[IndZ] += 1;
    }
    if( connect && elemconnect[j] < 0 ) continue;

    cumr[IndR] += 1;
    cumf[IndF] += 1;
  }
   
  /* Count the cumulative numbers */
  for(i=1;i<=MAXCATEGORY;i++) {
    cumr[i] = cumr[i] + cumr[i-1];
    cumf[i] = cumf[i] + cumf[i-1];
    if(dim==3) cumz[i] = cumz[i] + cumz[i-1];
  }

  j = noelements - data->elemconnectexist;
  for(i=1;i<=MAXCATEGORY;i++) {
    cumr[i] = ceil( 1.0 * partitions1 * cumr[i] / noelements );
    cumf[i] = ceil( 1.0 * partitions2 * cumf[i] / j );
    if(dim==3) cumz[i] = ceil( 1.0 * partitions3 * cumz[i] / j );
  }


  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j]%100;
    x = y = z = 0.0;
    for(i=0;i<nonodes;i++) {
      k = data->topology[j][i];
      x += data->x[k];
      y += data->y[k];
      z += data->z[k];
    }
    x = x / nonodes;
    y = y / nonodes;
    z = z / nonodes;
    
    r = sqrt(x*x+y*y);
    f = 180 * atan2(y,x)/FM_PI;
    if( f < 0.0 ) f = f + 360.0;    
    
    IndR = ceil( MAXCATEGORY * ( r - MinR ) / ( MaxR - MinR ) );
    IndF = ceil( MAXCATEGORY * ( f - MinF ) / ( MaxF - MinF ) );
    if(dim==3) IndZ = ceil( MAXCATEGORY * ( z - MinZ ) / ( MaxZ - MinZ ) );
    
    IndR = MIN( MAX( IndR, 1 ), MAXCATEGORY );
    IndF = MIN( MAX( IndF, 1 ), MAXCATEGORY );
    if(dim==3) IndZ = MIN( MAX( IndZ, 1 ), MAXCATEGORY );

    IndR = cumr[IndR];
    IndF = cumf[IndF];
    if(dim==3) {
      IndZ = cumz[IndZ];
    }
    else {
      IndZ = 1;
    }

    k = (IndZ-1) * partitions1 * partitions2 + 
      (IndF-1) * partitions1 + IndR;
    if( connect ) {
      if( elemconnect[j] < 0 ) 
	k = IndZ;
      else 
	k += 1;
    }
    inpart[j] = k;
    nopart[k] += 1;
  } 

  maxpart = 0;
  minpart = noelements;
  j = 0;

  /* Check whether some partition was not used */
  hit = FALSE;
  for(i=1;i<=partitions;i++) 
    if(!nopart[i]) hit = TRUE;

  if( hit ) {
    if( info ) printf("Renumbering the partitions (showing only 64):\n");
    for(i=1;i<=partitions;i++) {
      k = nopart[i];
      if( k ) {
	maxpart = MAX( k, maxpart );
	minpart = MIN( k, minpart );
	j += 1;
	nopart[i] = j;
	if( info && i!=j && j<=64) printf("%d -> %d %d\n",i,j,k);      
      }
    }
    if(info) printf("There are %d active partitions out of %d possible ones\n",j,partitions);
    data->nopartitions = j;

    for(i=1;i<=noelements;i++) {
      j = inpart[i];
      inpart[i] = nopart[j];
    }
  }
  else {
    if(info) printf("All possible %d out of %d partitions are active\n",j,partitions);
    data->nopartitions = partitions;
  }


  free_Ivector(nopart,1,partitions);
  free_Ivector( cumr,0,MAXCATEGORY);
  free_Ivector( cumf,0,MAXCATEGORY);
  if(dim==3) free_Ivector( cumz,0,MAXCATEGORY);

  PartitionNodesByElements(data,info);

  if(info) printf("Successfully made a partitioning with %d to %d elements.\n",minpart,maxpart);

  return(0);
}



int PartitionConnectedElementsStraight(struct FemType *data,struct BoundaryType *bound,
				       struct ElmergridType *eg, int info) {
  int i,j,k,l,dim,allocated,debug,partz,hit,bctype;
  int noknots, noelements,bcelem,bc,maxbcelem;
  int IndZ,noconnect,totpartelems,sideelemtype,sidenodes,sidehits,nohits;
  int *cumz,*elemconnect,*partelems,*nodeconnect;
  int sideind[MAXNODESD2];
  Real z,MaxZ,MinZ; 


  debug = FALSE;

  if(info) {
    printf("Link the connected set directly to the partition\n");
  }

  if(!data->nodeconnectexist) {
    printf("There are no connected boundary nodes?\n");
    return(1);
  }
  nodeconnect = data->nodeconnect;
  noknots = data->noknots;
  noelements = data->noelements;
  totpartelems = 0;

  bcelem = 0;  
  data->elemconnect = Ivector(1,noelements);
  elemconnect = data->elemconnect;
  for(i=1;i<=noelements;i++) elemconnect[i] = 0;

  for(bc=0;bc<MAXBOUNDARIES;bc++) {    
    if(bound[bc].created == FALSE) continue;
    if(bound[bc].nosides == 0) continue;
    
    for(i=1;i<=bound[bc].nosides;i++) {
      
      GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
		     data,sideind,&sideelemtype);

      sidenodes = sideelemtype % 100;

      hit = FALSE;
      
      for(k=1;k<=eg->connect;k++) {
	bctype = eg->connectbounds[k-1];
	hit = (bound[bc].types[i] == bctype);
	if(hit) break;
      }	
      if(!hit) continue;
      bcelem += 1;

      k = bound[bc].parent[i];      
      l = bound[bc].parent2[i];
      if(k) {
	if(!elemconnect[k]) {
	  elemconnect[k] = 1;
	  totpartelems += 1;
	}
      }
      if(l) {
	if(!elemconnect[l]) {
	  elemconnect[l] = 1;
	  totpartelems += 1;
	}
      }
    }	
  }

  data->elemconnectexist = totpartelems;  
  
  if(info) printf("Number of constrained boundary elements = %d\n",bcelem);   
  if(info) printf("Number of constrained bulk elements = %d\n",totpartelems);   
  if(info) printf("Successfully made connected elements to a partition\n");

  return(0);
}



int PartitionConnectedElements1D(struct FemType *data,struct BoundaryType *bound,
				 struct ElmergridType *eg, int info) {
  int i,j,k,l,dim,allocated,debug,partz,partr,parts,hit,bctype;
  int noknots, noelements,bcelem,bc,maxbcelem;
  int IndZ,noconnect,totpartelems,sideelemtype,sidenodes,sidehits,nohits;
  int *cumz,*elemconnect,*partelems,*nodeconnect;
  int sideind[MAXNODESD2];
  Real val,z,MaxZ,MinZ; 


  debug = FALSE;

  partz = eg->partbcz;
  partr = eg->partbcr;
  
  if( partz == 0 && partr == 0) return(0);

  parts = MAX( partz, partr ); 

  if(info) {
    if( partz )
      printf("Making a simple 1D partitioning in z for the connected elements only\n");
    else
      printf("Making a simple 1D partitioning in r for the connected elements only\n");     
  }

  if(!data->nodeconnectexist) {
    printf("There are no connected boundary nodes?\n");
    return(1);
  }
  nodeconnect = data->nodeconnect;

  dim = data->dim;
  if( dim < 3 ) {
    printf("PartitionConnectedElements1D is only applicable in 3D\n");
    return(2);
  }

  noknots = data->noknots;
  noelements = data->noelements;
  totpartelems = 0;

  /* Because we don't want to change the code too much use 'z' for the 
     coordinate also in the radial case. */ 
  if( partr )
    z = sqrt( data->x[1] * data->z[1] + data->y[1] * data->y[1]);
  else
    z = data->z[1];
  MaxZ = MinZ = z;
  
  for(i=1;i<=noknots;i++) {    
    if( partr )
      z = sqrt( data->x[i] * data->x[i] + data->y[i] * data->y[i]);
    else
      z = data->z[i];
    MaxZ = MAX( MaxZ, z);
    MinZ = MIN( MinZ, z);
  }

  if( info ) {
    printf("Range in coordinate extent: %12.5le %12.5le\n",MinZ,MaxZ);
  }

  /* Zero is the 1st value so that recursive algos can be used. */ 
  cumz = Ivector(0,MAXCATEGORY);
  for(i=0;i<=MAXCATEGORY;i++) 
    cumz[i] = 0;

  allocated = FALSE;
  bcelem = 0;

 omstart:

  for(bc=0;bc<MAXBOUNDARIES;bc++) {    
    if(bound[bc].created == FALSE) continue;
    if(bound[bc].nosides == 0) continue;
    
    for(i=1;i<=bound[bc].nosides;i++) {
      
      GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
		     data,sideind,&sideelemtype);

      sidenodes = sideelemtype % 100;
      hit = FALSE;
      
      for(k=1;k<=eg->connect;k++) {
	bctype = eg->connectbounds[k-1];
	if( bctype > 0 ) {
	  if(bound[bc].types[i] == bctype) hit = TRUE;
	} 
	else if( bctype == -1 ) {
	  if( bound[bc].parent[i] ) hit = TRUE;
	}
	else if( bctype == -2 ) {
	  if( bound[bc].parent[i] && bound[bc].parent2[i] ) hit = TRUE;
	}
	else if( bctype == -3 ) {
	  if( bound[bc].parent[i] && !bound[bc].parent2[i] ) hit = TRUE;
	}
	if(hit) break;
      }	
      if(!hit) continue;

      z = 0.0; 
      for(j=0;j<sidenodes;j++) {
	k = sideind[j];
	if( partr )
	  val = sqrt( data->x[k]*data->x[k] + data->y[k]*data->y[k]);
	else
	  val = data->z[k];
	z += val;
      }

      z = z / sidenodes;
      IndZ = ceil( MAXCATEGORY * ( z - MinZ ) / ( MaxZ - MinZ ) );

      /* To be on the safe side */
      IndZ = MIN( MAX( IndZ, 1 ), MAXCATEGORY );

     
      if(allocated) {
	IndZ = cumz[IndZ];	
	if(0) partelems[IndZ] += 1;

	/* Also set the connected elements if available */
	k = bound[bc].parent[i];      
	l = bound[bc].parent2[i];
	
	if(k) {
	  elemconnect[k] = IndZ;
	  totpartelems += 1;
	}
	if(l) {
	  elemconnect[l] = IndZ;
	  totpartelems += 1;
	}
      }
      else {      
	bcelem += 1;
	cumz[IndZ] += 1;
      }
    }	
  }

  if( !allocated )  {
    maxbcelem = bcelem;
    if( debug ) {
      printf("Number of constrained boundary elements = %d\n",bcelem);   
      printf("Differential categories (showing only 20 active ones from %d)\n",MAXCATEGORY);
      k = 0;
      for(j=0;j<=MAXCATEGORY;j++) {
	if( cumz[j] > 0 ) {
	  k++;
	  printf("%d : %d\n",j,cumz[j]);
	  if( k == 20 ) break;
	}
      }   
    }      

    /* Count the cumulative numbers and map them to number of partitions */
    for(i=1;i<=MAXCATEGORY;i++) 
      cumz[i] = cumz[i] + cumz[i-1];
    
    if( debug ) {
      printf("Cumulative categories\n");
      for(j=0;j<=MAXCATEGORY;j++) {
	printf("%d : %d\n",j,cumz[j]);
      }   
    }
    
    noconnect = bcelem;
    for(i=1;i<=MAXCATEGORY;i++) 
      cumz[i] = ceil( 1.0 * parts * cumz[i] / noconnect );
    
    if( debug ) {
      printf("Partition categories\n");
      for(j=0;j<=MAXCATEGORY;j++) {
	printf("%d : %d\n",j,cumz[j]);
      }   
    }

    data->elemconnect = Ivector(1,noelements);
    elemconnect = data->elemconnect;

    for(i=1;i<=noelements;i++) elemconnect[i] = 0;

    allocated = TRUE;
    if(info) printf("Allocated and now setting the entries\n");
    goto omstart;
  }
    
  data->elemconnectexist = totpartelems;  
  free_Ivector( cumz,0,MAXCATEGORY);

  if(info) printf("Successfully made a partitioning 1D with %d elements\n",totpartelems);

  return(0);
}


#if USE_METIS
static int SetMetisOptions(idx_t *options, struct ElmergridType *eg,int info)  {
  
  if(info) printf("Setting default Metis options!\n");
  METIS_SetDefaultOptions(options);  

  if(info) printf("Setting user defined Metis options!\n");
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_DBGLVL] = 3;
  options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
  options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
  options[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;    

  if( eg->metis_contig ) options[METIS_OPTION_CONTIG] = 1;
  if( eg->metis_seed ) options[METIS_OPTION_SEED] = eg->metis_seed;
  if( eg->metis_ncuts > 1) options[METIS_OPTION_NCUTS] = eg->metis_ncuts;
  if( eg->metis_volcut) options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  if( eg->metis_minconn) options[METIS_OPTION_MINCONN] = 1;
}


int PartitionConnectedElementsMetis(struct FemType *data,struct BoundaryType *bound,
				    struct ElmergridType *eg, int nparts,int metisopt,int info) {
  int i,j,k,l,n,m,dim;
  int noknots,noelements,sidenodes,ind,nohits,totpartelems;
  int dualmaxcon,invmaxcon,totcon,step,bc,bcelem,bcelem2,hit,set;
  int noconnect,minpartelems,maxpartelems,sideelemtype,con,maxbcelem;
  int *elemconnect,*partelems,*nodeperm,*neededby;
  int *bcdualgraph[MAXCONNECTIONS],*bcinvtopo[MAXCONNECTIONS];
  int sideind[MAXNODESD1];  
  int nn;
  idx_t options[METIS_NOPTIONS]; 
  int *xadj,*adjncy,*vwgt,*adjwgt,wgtflag,*npart;
  int numflag,edgecut,ncon;
  int *nodepart;
  

  if(info) {
    printf("Making a Metis partitioning for the connected BC elements only\n");
  }

  if(!data->nodeconnectexist) {
    printf("There are no connected elements\n");
    return(1);
  }
  
  SetMetisOptions(options,eg,info);
  
  dim = data->dim;
  noknots = data->noknots;
  noelements = data->noelements;
  
  /* Count the connected nodes and set the permutation */
  noconnect = 0;
  for(i=1;i<=noknots;i++) 
    if( data->nodeconnect[i] ) noconnect++;
  if( noconnect == 0 ) {
    printf("There are really no connected nodes?\n");
    return(2);
  }
  if(info) printf("Number of connected nodes is %d\n",noconnect);
  
  /* Permute the nodes so that they only include the connected nodes. */
  nodeperm = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) nodeperm[i] = 0;
  noconnect = 0;
  for(i=1;i<=noknots;i++) 
    if( data->nodeconnect[i] ) {
      noconnect++;
      nodeperm[i] = noconnect;
    }
  if(info) printf("Permuted the connected nodes\n");
  
  /* Cycle over the boundary elements
     1) First cycle create a table showing the elements sharing a node
     2) Second cycle create a table showing the element-to-element connections */
  
  invmaxcon = 0;
  dualmaxcon = 0;
  totcon = 0;
  
  neededby = Ivector(1,noconnect);
  for(i=1;i<=noconnect;i++) neededby[i] = 0;
  
  for(step=1;step<=2;step++) {

    bcelem = 0;
    for(bc=0;bc<MAXBOUNDARIES;bc++) {    
      if(bound[bc].created == FALSE) continue;
      if(bound[bc].nosides == 0) continue;
      
      for(i=1;i<=bound[bc].nosides;i++) {
	
	GetBoundaryElement(i,&bound[bc],data,sideind,&sideelemtype); 
	/* GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
	   data,sideind,&sideelemtype); */
	sidenodes = sideelemtype % 100;
	nohits = 0;
	
	for(j=0;j<sidenodes;j++) 
	  if( nodeperm[sideind[j]] ) nohits++;
	if( nohits < sidenodes ) continue;
	bcelem++;
	
	for(j=0;j<sidenodes;j++) {
	  ind = nodeperm[sideind[j]];
	  
	  /* Create table that shows all the elements sharing a node. */   
	  if( step == 1 ) {
	    neededby[ind] += 1;
	    k = neededby[ind];
	    
	    if(k > invmaxcon) {
	      invmaxcon++;
	      bcinvtopo[invmaxcon] = Ivector(1,noconnect);
	      if(0) printf("allocating invtopo %d %d\n",invmaxcon,noconnect);
	      for(m=1;m<=noconnect;m++)
		bcinvtopo[invmaxcon][m] = 0;
	    }
	    bcinvtopo[k][ind] = bcelem;
	  } 
	  else if( step == 2 ) {	    
	    /* Now create a table that shows how elements are connected to other elements, 
	       this is the dual graph. */
	    for(k=1;k<=invmaxcon;k++) {
	      bcelem2 = bcinvtopo[k][ind];
	
	      if( bcelem2 == 0 ) break;
	      if( bcelem2 == bcelem ) continue;
	      
	      hit = FALSE;
	      for(l=1;l<=dualmaxcon;l++) { 
		if(bcdualgraph[l][bcelem] == bcelem2) hit = TRUE;
		if(bcdualgraph[l][bcelem] == 0) break;
	      }
	      if(!hit) {
		if(l > dualmaxcon) {
		  dualmaxcon++;
		  if( l >= MAXCONNECTIONS ) {
		    printf("Number of connections %d vs. static limit %d\n",l,MAXCONNECTIONS-1);
		    bigerror("Maximum of connections in dual graph larger than the static limit!");
		  }
		  if(0) printf("allocating dual topo %d %d\n",dualmaxcon,maxbcelem);
		  bcdualgraph[dualmaxcon] = Ivector(1,maxbcelem);
		  for(m=1;m<=maxbcelem;m++)
		    bcdualgraph[dualmaxcon][m] = 0;
		}
		bcdualgraph[l][bcelem] = bcelem2;	    
		totcon++;
	      }
	    }	      	      
	  }
	}
      }
    }
    maxbcelem = bcelem;
  }    
  
  if(info) {
    printf("There are %d connected boundary elements.\n",maxbcelem);
    printf("There are %d connections altogether in the dual graph.\n",totcon);
  }  

  /* Create the sparse matrix format graph for Metis */
  xadj = Ivector(0,maxbcelem);
  adjncy = Ivector(0,totcon-1);
  for(i=0;i<totcon;i++) 
    adjncy[i] = 0;
  
  totcon = 0;
  for(i=1;i<=maxbcelem;i++) {
    xadj[i-1] = totcon;
    for(j=1;j<=dualmaxcon;j++) {
      con = bcdualgraph[j][i];
      if(con) {
	adjncy[totcon] = con-1;
	totcon++;
	if(0) printf("con: %d %d %d\n",i,totcon,con);
      }
    }
  }
  xadj[maxbcelem] = totcon;
  
  /* Deallocate the temporal graphs */
  free_Ivector( neededby, 1, noconnect );
  for(i=1;i<=invmaxcon;i++)
    free_Ivector(bcinvtopo[i],1,noconnect);
  for(i=1;i<=dualmaxcon;i++)
    free_Ivector(bcdualgraph[i],1,maxbcelem);
  
  /* Parameters for Metis */
  numflag = 0;
  nn = maxbcelem;
  npart = Ivector(0,nn-1);
  wgtflag = 0;

  /* Optional weights */
  vwgt = NULL;
  adjwgt = NULL;

  if(metisopt == 2) {
    if(info) printf("Starting graph partitioning METIS_PartGraphRecursive.\n");
    METIS_PartGraphRecursive(&nn,&ncon,xadj,adjncy,vwgt,&wgtflag,adjwgt,
			     &nparts,NULL,NULL,options,&edgecut,npart); 
  }
  else if(metisopt == 3 || metisopt == 4) {
    if(info) printf("Starting graph partitioning METIS_PartGraphKway.\n");      
    ncon = 0;
    wgtflag = 0;    
    METIS_PartGraphKway(&nn,&ncon,xadj,adjncy,vwgt,&wgtflag,adjwgt,
			&nparts,NULL,NULL,options,&edgecut,npart); 
  }
  else {
    printf("Unknown Metis option %d\n",metisopt);
  }
  if(0) printf("Finished Metis routine for boundary partitioning\n");


  free_Ivector(xadj, 0,maxbcelem);
  free_Ivector(adjncy, 0,totcon-1);


  data->elemconnect = Ivector(1,noelements);
  elemconnect = data->elemconnect;
  for(i=1;i<=noelements;i++) elemconnect[i] = 0;

  /* Set the bulk element partitions by following the boundary elements. 
     Also set the nodal partitions to apply majority rule later on. */     
  bcelem = 0;
  totpartelems = 0;
  for(bc=0;bc<MAXBOUNDARIES;bc++) {    
    if(bound[bc].created == FALSE) continue;
    if(bound[bc].nosides == 0) continue;
    
    for(i=1;i<=bound[bc].nosides;i++) {
      
      GetBoundaryElement(i,&bound[bc],data,sideind,&sideelemtype); 
      /* GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
	 data,sideind,&sideelemtype); */
      sidenodes = sideelemtype % 100;
      nohits = 0;
      for(j=0;j<sidenodes;j++) 
	if( nodeperm[sideind[j]] ) nohits++;
      if( nohits < sidenodes ) continue;
    
      /* Always set the connected element to a negative value as this information is used
	 in a dirty way later on. */
      k = bound[bc].parent[i];      
      l = bound[bc].parent2[i];

      if(0) printf("kl = %d %d %d %d %d\n",k,l,bcelem,npart[bcelem]+1,elemconnect[k]);
      
      /* Have the addition here since npart has indexes [0,nn-1]. */
      if(k) {
	elemconnect[k] = npart[bcelem]+1;
	totpartelems += 1;
      }
      if(l) {
	elemconnect[l] = npart[bcelem]+1;
	totpartelems += 1;
      }
      bcelem++;
    }
  }

  free_Ivector( nodeperm, 1, noknots ); 
  free_Ivector( npart, 0, nn-1 );

  data->elemconnectexist = totpartelems;  
  if(info) printf("Successfully created boundary partitioning for %d bulk elements\n",totpartelems);

  return(0);
}
#endif


int ExtendBoundaryPartitioning(struct FemType *data,struct BoundaryType *bound,
			       int elemlayers,int info) 
{
  int i,j,k,l,m,n,nonodes,noknots,noelements,totpartelems,nparts,set;
  int minpartelems,maxpartelems,refcount,refpart,part;
  int *partelems,*elemconnect;
  int *invrow,*invcol;

  printf("Extending boundary partitioning by majority rule\n");

  CreateInverseTopology(data, info);
  invrow = data->invtopo.rows;
  invcol = data->invtopo.cols;
  
  noelements = data->noelements;
  noknots = data->noknots;
  
  if( !data->elemconnectexist ) {
    printf("Initial partitioning does not exist!\n");
    return(1);
  }
  elemconnect = data->elemconnect;
  
  nparts = 0;
  for(i=1;i<=data->noelements;i++)
    nparts = MAX( nparts, elemconnect[i] );
  
  partelems = Ivector(1,nparts);
  for(i=1;i<=nparts;i++)
    partelems[i] = 0;
  
  totpartelems = 0;
  for(i=1;i<=noelements;i++) {
    j = elemconnect[i];
    if(j) {
      partelems[j] += 1;
      totpartelems +=1;
    }
  }
  
  if(info) printf("Initial partitioning with %d elements:\n",totpartelems);
  minpartelems = maxpartelems = partelems[1];
  for(j=1;j<=nparts;j++) {
    minpartelems = MIN( minpartelems, partelems[j] );
    maxpartelems = MAX( maxpartelems, partelems[j] );
  }
  printf("Initial partitioning has %d to %d elements in partition\n",minpartelems,maxpartelems);


  /* Counter to find the dominating element */
  for(i=1;i<=nparts;i++)
    partelems[i] = 0;
  
  /* Apply majority rule to set the remaining bulk elements that 
     share nodes on the connected boundary. */
  for(m=1;m<=elemlayers;m++) {
    n = 0;
    for(j=1;j<=noelements;j++) {

      /* This element is already set */
      if( elemconnect[j] > 0 ) continue;
      
      set = FALSE;
      nonodes = data->elementtypes[j] % 100;
      for(i=0;i<nonodes;i++) {	
	k = data->topology[j][i];
	for(l=invrow[k-1];l<invrow[k];l++) {
	  part = elemconnect[invcol[l]+1];	 
	  if( part > 0 ) {
	    set = TRUE;
	    partelems[part] += 1;
	  }
	}
      }
      /* No hits, cannot set */
      if(!set) continue;
      
      refcount = 0;
      for(i=0;i<nonodes;i++) {
	k = data->topology[j][i];
	for(l=invrow[k-1];l<invrow[k];l++) {
	  part = elemconnect[invcol[l]+1];	 
	  if(part > 0 ) {
	    if( partelems[part] > refcount ) {
	      refcount = partelems[part];
	      refpart = part;
	    }
	    partelems[part] = 0;
	  }
	}
      }
      
      if(0) printf("setting: %d %d %d\n",j,refpart,refcount);

      /* Set negative entry first so that this is not used within this loop */
      elemconnect[j] = -refpart;
      n++;
    }

    /* Revert to positive indexes at the end */
    for(j=1;j<=noelements;j++)
      elemconnect[j] = abs( elemconnect[j] );

    if(info) printf("Setting partition for %d elements in layer %d\n",n,m);
    if( n == 0 ) {
      printf("No elements set in the layer, breaking loop %d!\n",m);
      break; 
    }
  }
  
  
  /* Compute the division to partitions */
  for(i=1;i<=nparts;i++)
    partelems[i] = 0;
  
  totpartelems = 0;
  for(i=1;i<=noelements;i++) {
    j = elemconnect[i];
    if(j) {
      partelems[j] += 1;
      totpartelems += 1;
    }
  }
  
  if(info) printf("Extended partitioning of %d elements:\n",totpartelems);
  minpartelems = maxpartelems = partelems[1];
  for(j=1;j<=nparts;j++) {
    minpartelems = MIN( minpartelems, partelems[j] );
    maxpartelems = MAX( maxpartelems, partelems[j] );
    if(0) printf("Part: %d %d\n",j,partelems[j]);
  }
  if(info) printf("Extended partitioning has %d to %d elements in partition\n",minpartelems,maxpartelems);

  data->elemconnectexist = totpartelems;
  free_Ivector( partelems, 1, nparts );
  
  if(info) printf("Successfully extended the partitioning by %d layers\n",elemlayers);
  
  return(0);
 }





int PartitionSimpleNodes(struct FemType *data,int dimpart[],int dimper[],
			 int partorder, Real corder[],Real parttol,int info)
/* Do simple partitioning going for the nodes. Then using the the majority rule
   define the elemental partitioning. This is suboptimal for Elmer because the 
   elemental ordering is primary in Elmer. */
{
  int i,j,k,k0,kprev,l,ind,minpart,maxpart;
  int noknots, noelements,elemsinpart,periodic;
  int partitions1,partitions2,partitions3,partitions;
  int vpartitions1,vpartitions2,vpartitions3,vpartitions,hit;
  int *indx,*nopart,*nodepart;
  Real *arrange;
  Real x,y,z,cx,cy,cz;
  
  partitions1 = dimpart[0];
  partitions2 = dimpart[1];
  partitions3 = dimpart[2];
  if(data->dim < 3) partitions3 = 1;
  partitions = partitions1 * partitions2 * partitions3;

  if(partitions1 < 2 && partitions2 < 2 && partitions3 < 2) {
    printf("Some of the divisions must be larger than one: %d %d %d\n",
	   partitions1, partitions2, partitions3 );
    bigerror("Partitioning not performed");
  }

  if(partitions >= data->noelements) {
    printf("There must be fever partitions than elements (%d vs %d)!\n",
	   partitions,data->noelements);
    bigerror("Partitioning not performed");
  }
    
  if(!data->partitionexist) {
    data->partitionexist = TRUE;
    data->elempart = Ivector(1,data->noelements);
    data->nodepart = Ivector(1,data->noknots);
    data->nopartitions = partitions;
  }
  nodepart = data->nodepart;

  vpartitions1 = partitions1;
  vpartitions2 = partitions2;
  vpartitions3 = partitions3;
  periodic = dimper[0] || dimper[1] || dimper[2];
  if(periodic) {
    if(dimper[0] && partitions1 > 1) vpartitions1 *= 2;
    if(dimper[1] && partitions2 > 1) vpartitions2 *= 2;
    if(dimper[2] && partitions3 > 1) vpartitions3 *= 2;
  }
  vpartitions = vpartitions1 * vpartitions2 * vpartitions3;

  nopart = Ivector(1,vpartitions);
  noelements = data->noelements;
  noknots = data->noknots;

  if(info) printf("Making a simple partitioning for %d nodes in %d-dimensions.\n",
		  noknots,data->dim);

  arrange = Rvector(1,noknots);
  indx = Ivector(1,noknots);

  if(partorder) {
    cx = corder[0];
    cy = corder[1];
    cz = corder[2];    
  }
  else {
    cx = 1.0;
    cy = 0.0001;
    cz = cy*cy;
  }

  z = 0.0;

  for(i=1;i<=noknots;i++) 
    nodepart[i] = 1;  

  if(vpartitions1 > 1) {
    if(info) printf("Ordering 1st direction with (%.3lg*x + %.3lg*y + %.3lg*z)\n",cx,cy,cz);
    for(j=1;j<=noknots;j++) {
      x = data->x[j];
      y = data->y[j];
      z = data->z[j];
      arrange[j] = cx*x + cy*y + cz*z;
    }
    SortIndex(noknots,arrange,indx);

    for(i=1;i<=noknots;i++) {
      ind = indx[i];
      k = (i*vpartitions1-1)/noknots+1;
      nodepart[ind] = k;
    }
  } 

  /* Partition in the 2nd direction taking into account the 1st direction */
  if(vpartitions2 > 1) {
    if(info) printf("Ordering in the 2nd direction.\n");
    for(j=1;j<=noknots;j++) {
      x = data->x[j];
      y = data->y[j];
      z = data->z[j];
      arrange[j] = -cy*x + cx*y + cz*z;
    }
    SortIndex(noknots,arrange,indx);
    
    for(i=1;i<=vpartitions;i++)
      nopart[i] = 0;
    
    elemsinpart = noknots / (vpartitions1*vpartitions2);
    j = 1;
    for(i=1;i<=noknots;i++) {
      ind = indx[i];
      k0 = (nodepart[ind]-1) * vpartitions2;
      for(l=1;l<=vpartitions2;l++) {
	hit = FALSE;
	if( j < vpartitions ) {
	  if( nopart[k0+j] >= elemsinpart ) {
	    j += 1;
	    hit = TRUE;
	  }
	}
	if( j > 1 ) {
	  if( nopart[k0+j-1] < elemsinpart ) {
	    j -= 1;
	    hit = TRUE;
	  }
	}
	if( !hit ) break;
      }
      k = k0 + j;
      nopart[k] += 1;
      nodepart[ind] = k;
    }
  }  

  /* Partition in the 3rd direction taking into account the 1st and 2nd direction */
  if(vpartitions3 > 1) {
    if(info) printf("Ordering in the 3rd direction.\n");
    for(j=1;j<=noknots;j++) {
      x = data->x[j];
      y = data->y[j];
      z = data->z[j];
      arrange[j] = -cz*x - cy*y + cx*z;
    }
    SortIndex(noknots,arrange,indx);

    for(i=1;i<=vpartitions;i++)
      nopart[i] = 0;
    
    elemsinpart = noknots / (vpartitions1*vpartitions2*vpartitions3);
    j = 1;
    for(i=1;i<=noknots;i++) {
      ind = indx[i];
      k0 = (nodepart[ind]-1)*vpartitions3;

      for(l=1;l<=vpartitions;l++) {
	hit = FALSE;
	if( j < vpartitions3 ) {
	  if( nopart[k0+j] >= elemsinpart ) {
	    j += 1;
	    hit = TRUE;
	  }
	}
	if( j > 1 ) {
	  if( nopart[k0+j-1] < elemsinpart ) {
	    j -= 1;
	    hit = TRUE;
	  }
	}
	if( !hit ) break;
      }
      k = k0 + j;
      nopart[k] += 1;
      nodepart[ind] = k;
    }
  }
  

  /* This piece of code may be used to assign "strides" back to same partition. */
  if(parttol > 1.0e-20 )  {
    if(0) {
      printf("Original partition counts:\n");
      for(i=1;i<=partitions;i++)
	printf("nopart[%d] = %d\n",i,nopart[i]);
    }
      
    for(j=1;j<=noknots;j++) {
      x = data->x[j];
      y = data->y[j];
      z = data->z[j];
      arrange[j] = cx*x + cy*y + cz*z;
    }
    SortIndex(noknots,arrange,indx);

    kprev = nodepart[indx[1]];
    l = 0;
    for(i=2;i<=noknots;i++) {
      ind = indx[i];
      k = nodepart[ind];
      if(k != kprev ) {
	if(fabs(arrange[i]-arrange[i-1]) < parttol ) {
	  nopart[k] -= 1;
	  nopart[kprev] += 1;
	  nodepart[ind] = kprev;
	  k = kprev;
	  if(0) printf("arrange %d %d %d %d %.6lg %.6lg\n",i,ind,k,kprev,arrange[i],arrange[i-1]);
	  l++;
	}
      }
      kprev = k;
    }
    if(0) {
      printf("Updated partition counts:\n");
      for(i=1;i<=partitions;i++)
	printf("nopart[%d] = %d\n",i,nopart[i]);    
    }      
    if(info) printf("Number of partition changes due to tolerance: %d\n",l);
  }

  
  /* For periodic systems the number of virtual partitions is larger. Now map the mesh so that the 
     1st and last partition for each direction will be joined */
  if(periodic) {
    int *partmap;
    int p1,p2,p3,q1,q2,q3;
    int P,Q;
    p1=p2=p3=1;
    partmap = Ivector(1,vpartitions);
    for(i=1;i<=vpartitions;i++)
      partmap[i] = 0;
    for(p1=1;p1<=vpartitions1;p1++) {
      q1 = p1;
      if(dimper[0] && vpartitions1 > 1) {
        if(q1==vpartitions1) q1 = 0;
        q1 = q1/2 + 1;
      }
      for(p2=1;p2<=vpartitions2;p2++) {
        q2 = p2;
        if(dimper[1] && vpartitions2 > 1) {
          if(q2==vpartitions2) q2 = 0;
          q2 = q2/2 + 1;
        }
        for(p3=1;p3<=vpartitions3;p3++) {
          q3 = p3;
          if(dimper[2] && vpartitions3 > 1) {
            if(q3==vpartitions3) q3 = 0;
            q3 = q3/2 + 1;
          }
	  
          P = vpartitions3 * vpartitions2 * (p1 - 1) + vpartitions3 * (p2-1) + p3;
          Q = partitions3 * partitions2 * (q1 - 1) + partitions3 * (q2-1) + q3;

          partmap[P] = Q;
        }
      }
    }
    for(i=1;i<=noknots;i++)
      nodepart[i] = partmap[nodepart[i]];
    free_Ivector(partmap,1,vpartitions);
  }


  for(i=1;i<=partitions;i++)
    nopart[i] = 0;
  for(i=1;i<=noknots;i++) 
    nopart[nodepart[i]] += 1;
  
  minpart = maxpart = nopart[1];
  for(i=1;i<=partitions;i++) {
    minpart = MIN( nopart[i], minpart );
    maxpart = MAX( nopart[i], maxpart );
  }

  free_Rvector(arrange,1,noelements);
  free_Ivector(nopart,1,partitions);
  free_Ivector(indx,1,noelements);

  PartitionElementsByNodes(data,info);

  if(info) printf("Successfully made a partitioning with %d to %d nodes.\n",minpart,maxpart);

  return(0);
}


static int LinearNodes(int elemtype)
{
  int elemfamily,nonodes;
  int fam2nodemap[] = {0, 1, 2, 3, 4, 4, 5, 6, 8 }; 
  
  elemfamily = elemtype / 100;
  nonodes = fam2nodemap[elemfamily];

  return(nonodes);
}


int PartitionMetisMesh(struct FemType *data,struct ElmergridType *eg,
		       int partitions,int dual,int info)
/* Perform partitioning using Metis. This uses the elemental routines of Metis that assume that 
   there exists only one elementtype. If this condition is not met then this routine cannot be 
   used. If the elements are higher order nodal elements then use only the linear basis. */
{
  int i,j,k,periodic, noelements, noknots, sides, highorder;
  int nodesd2, etype, numflag,mintype,maxtype,elemtype,minnodes,ninv;
  int *neededby,*indxper,*inpart;
  idx_t *metistopo,*eptr,*npart,*epart;
  idx_t ne,nn,ncommon,edgecut,nparts;
  idx_t options[METIS_NOPTIONS];


  SetMetisOptions(options,eg,info);

  if(info) printf("Making a Metis partitioning for %d elements in %d-dimensions.\n",
		  data->noelements,data->dim);

  noelements = data->noelements;
  noknots = data->noknots;
 
  for(i=1;i<=noelements;i++) {
    elemtype = data->elementtypes[i];
    nodesd2 = LinearNodes( elemtype );
    if(i == 1 ) {
      mintype = maxtype = elemtype;
      minnodes = nodesd2;
    } else {
      mintype = MIN( mintype, elemtype );
      maxtype = MAX( maxtype, elemtype );
      minnodes = MIN( minnodes, nodesd2 );
    }
  }

  if(info) {
    if(mintype == maxtype ) {    
      printf("All elements are of type %d\n",mintype);
    }
    else {
      printf("Minimum element type is %d\n",mintype);
      printf("Maximum element type is %d\n",maxtype);
    }
  }
   
  if( minnodes >= 8) {
    ncommon = 4;
  }
  else if(minnodes > 4 ) {
    ncommon = 3;
  }
  else {
    ncommon = 2;
  }
  if(info) {
    printf("Minimum number of linear nodes in elements %d\n",minnodes);
    printf("Requiring number of nodes in dual graph %d\n",ncommon);
  }
    
  if(!data->partitionexist) {
    data->partitionexist = TRUE;
    data->elempart = Ivector(1,data->noelements);
    data->nodepart = Ivector(1,data->noknots);
    data->nopartitions = partitions;
  }
  inpart = data->elempart;

  /* Are there periodic boundaries. This information is used to join the boundaries. */
  periodic = data->periodicexist;
  if(periodic) {
    if(info) printf("There seems to be periodic boundaries\n");
    indxper = data->periodic;
  }
  
  ne = noelements;
  nn = noknots;

  neededby = Ivector(1,noknots);
  
  numflag = 0;
  nparts = partitions;

  /* Mark the nodes that are needed */
  k = 0;
  for(i=1;i<=noknots;i++) 
    neededby[i] = 0;
  if(periodic) {
    for(i=1;i<=noelements;i++) {
      nodesd2 = LinearNodes( data->elementtypes[i] );
      for(j=0;j<nodesd2;j++) {
        neededby[indxper[data->topology[i][j]]] = 1;
	k += 1;
      }
    }
  }
  else {
    for(i=1;i<=noelements;i++) {
      nodesd2 = LinearNodes( data->elementtypes[i] );
      for(j=0;j<nodesd2;j++) {
        neededby[data->topology[i][j]] = 1;
	k += 1;
      }
    }
  }

  j = 0;
  for(i=1;i<=noknots;i++) 
    if(neededby[i]) 
      neededby[i] = ++j;

  nn = j;
  if(info) {
    if(nn == noknots)
      printf("Using all %d possible nodes in the Metis graph\n",nn);
    else
      printf("Using %d nodes of %d possible nodes in the Metis graph\n",nn,noknots);
    printf("Allocating mesh topology of size %d\n",k);
  }

  eptr = Ivector(0,noelements);
  metistopo = Ivector(0,k-1);

  k = 0;
  eptr[0] = k;
  if(periodic) {
    for(i=1;i<=noelements;i++) {
      nodesd2 = LinearNodes( data->elementtypes[i] );
      for(j=0;j<nodesd2;j++) {
        metistopo[k] = neededby[indxper[data->topology[i][j]]]-1;
	k += 1;
      }
      eptr[i] = k;
    }
  }    
  else if(nn < noknots) {
    for(i=1;i<=noelements;i++) {
      nodesd2 = LinearNodes( data->elementtypes[i] );
      for(j=0;j<nodesd2;j++) {
        metistopo[k] = neededby[data->topology[i][j]]-1;    
	k += 1;
      }
      eptr[i] = k;
    }
  }
  else {
    for(i=1;i<=noelements;i++) {
      nodesd2 = LinearNodes( data->elementtypes[i] );
      for(j=0;j<nodesd2;j++) {
        metistopo[k] = data->topology[i][j]-1;    
	k += 1;
      }
      eptr[i] = k;
    }
  }

  npart = Ivector(0,nn-1);
  for(i=0;i<nn;i++) npart[i] = 0;
  
  epart = Ivector(0,ne-1);
  for(i=0;i<ne;i++) epart[i] = 0;
  
  if(dual) {
    if(info) printf("Starting graph partitioning METIS_PartMeshDual.\n");  
    METIS_PartMeshDual(&ne,&nn,eptr,metistopo,NULL,NULL,&ncommon,
		       &nparts,NULL,options,&edgecut,epart,npart);
    if(info) printf("Finished graph partitioning METIS_PartMeshDual.\n");  
  }
  else {
    if(info) printf("Starting graph partitioning METIS_PartMeshNodal.\n");  
    METIS_PartMeshNodal(&ne,&nn,eptr,metistopo,NULL,NULL,
			&nparts,NULL,options,&edgecut,epart,npart);
    if(info) printf("Finished graph partitioning METIS_PartMeshNodal.\n");  
  }

  /* Set the partition given by Metis for each element. */
  ninv = 0;
  for(i=1;i<=noelements;i++) {
    inpart[i] = epart[i-1]+1;
    if(inpart[i] < 1 || inpart[i] > partitions) {
      if(ninv < 10 ) printf("Invalid partition %d for element %d\n",inpart[i],i);
      ninv++;
    }
  }
  if(ninv) {
    printf("Number of invalid element partitions by Metis %d\n",ninv);
    bigerror("Cannot continue with invalid partitioning!");
  }
	     
  if( highorder ) {
    PartitionNodesByElements(data,info);
  }
  else {
    if(info) printf("Set the partition given by Metis for each node\n");
    for(i=1;i<=noknots;i++) {
      if(periodic) 
	j = neededby[indxper[i]];
      else if(nn < noknots)
	j = neededby[i];
      else
	j = i;
      if(!j) printf("Cannot set partitioning for node %d\n",i);
      data->nodepart[i] = npart[j-1]+1;
      if(data->nodepart[i] < 1 || data->nodepart[i] > partitions) {
        if(ninv < 10) printf("Invalid partition %d for node %d\n",data->nodepart[i],i);
	ninv++;
      }
    }
    if(ninv) {
      printf("Number of invalid node partitions by Metis %d\n",ninv);
      bigerror("Cannot continue with invalid partitioning!");
    }
  }
  
  free_Ivector(neededby,1,noknots);
  free_Ivector(metistopo,0,k-1);
  free_Ivector(epart,0,noelements-1);
  free_Ivector(npart,0,nn-1);
  free_Ivector(eptr,0,noelements+1);

  if(info) printf("Successfully made a Metis partition using the element mesh.\n");

  return(0);
}



#if USE_METIS
int PartitionMetisGraph(struct FemType *data,struct BoundaryType *bound,
			struct ElmergridType *eg,int partitions,int metisopt,
			int dual,int info)
/* Perform partitioning using Metis. One may use either the direct or the dual
   graph. The direct graph means that the nodes are first partitioned and the 
   elements follow. The dual graph means that the elements are partitioned and 
   the ownership of the nodes will follow. The latter is optimal for Elmer. */
{
  int i,j,k,noelements,noknots,errstat;
  int nn,ncon,con,maxcon,totcon;
  int *xadj,*adjncy,*vwgt,*adjwgt,wgtflag,*npart,**graph;
  int numflag,nparts,edgecut,maxconset;
  struct CRSType *dualgraph;
  idx_t options[METIS_NOPTIONS];

  if(info) printf("Making a Metis partitioning for %d nodes in %d-dimensions.\n",
		  data->noknots,data->dim);
  if(partitions < 2 ) {
    bigerror("There should be at least two partitions for partitioning!");
  }

  noknots = data->noknots;
  noelements = data->noelements;
  maxconset = 0;

  if(data->periodicexist && dual ) {
    printf("Dual graph not implemented for periodic system!\n");
    printf("Enforcing nodal graph for partitioning\n");
    dual = FALSE;
  }

  SetMetisOptions(options,eg,info);
      
  nparts = partitions;
  if( dual ) {
    if( eg->partbcz > 1 || eg->partbcr ) 
      PartitionConnectedElements1D(data,bound,eg,info);
    else if( eg->partbcmetis > 1 ) 
      PartitionConnectedElementsMetis(data,bound,eg,eg->partbcmetis,metisopt,info);
    else if( eg->connect ) {
      PartitionConnectedElementsStraight(data,bound,eg,info);
    }    

    if( data->nodeconnectexist ) {
      errstat = ExtendBoundaryPartitioning(data,bound,eg->partbclayers,info);
      if( errstat ) bigerror("Extend boundary partitioning returned error!");
    }
    
    /* Find the number of partitions already used for connected elements */
    if( data->elemconnectexist ) {
      maxconset = 0;
      for(i=1;i<=noelements;i++)
	maxconset = MAX( maxconset, data->elemconnect[i] );
      nparts -= maxconset;
    }

    if( nparts == 1 ) {
      if(info) printf("Just one partition left, skipping 2nd part of hybrid partitioning!\n");
      goto skippart;
    }

    CreateDualGraph(data,TRUE,info);

    /* Create the sparse matrix format graph for Metis */
    dualgraph = &data->dualgraph; 
    
    nn = dualgraph->rowsize;
    xadj = dualgraph->rows;
    adjncy = dualgraph->cols;
    totcon = dualgraph->colsize;
  }
  else {
    CreateNodalGraph(data,TRUE,info);
    maxcon = data->nodalmaxconnections;
    nn = noknots;
    graph = data->nodalgraph;

    /* Compute total number of connections in graph */
    totcon = 0;
    for(i=1;i<=nn;i++) {
      for(j=0;j<maxcon;j++) {
        con = graph[j][i];
        if(con) totcon++;
      }
    }
    
    /* Create the sparse matrix format graph for Metis */
    xadj = Ivector(0,nn);
    adjncy = Ivector(0,totcon-1);
    for(i=0;i<totcon;i++) 
      adjncy[i] = 0;
    
    totcon = 0;
    for(i=1;i<=nn;i++) {
      xadj[i-1] = totcon;
      for(j=0;j<maxcon;j++) {
        con = graph[j][i];
        if(con) {
          adjncy[totcon] = con-1;
          totcon++;
        }
      }
    }
    xadj[nn] = totcon;
  }

  if(info) printf("There are %d connections altogether in the graph.\n",totcon);

  /* Parameters for Metis */
  numflag = 0;
  npart = Ivector(0,nn-1);
  wgtflag = 0;
  ncon = 1;
  
  /* Optional weights */
  vwgt = NULL;
  adjwgt = NULL;

  if( !dual ) {
    /* Create weights if needed */
    if(data->periodicexist || eg->connect) {
      if(info) printf("Creating weight for %d connections.\n",totcon);
      wgtflag = 1;
      adjwgt = Ivector(0,totcon-1);
      for(i=0;i<totcon;i++)
        adjwgt[i] = 1;
      
      if(metisopt != 3) {
        printf("For weighted partitioning Metis subroutine METIS_PartGraphKway is enforced\n");
        metisopt = 3;
      }
    }
    
    /* Make the periodic connections stronger */
    if(data->periodicexist) {
      if(info) printf("Setting periodic connections to dominate %d\n",totcon);
      for(i=0;i<noknots;i++) {
        j = data->periodic[i+1]-1;
        if(j == i) continue;
        for(k=xadj[i];k<xadj[i+1];k++) 
          if(adjncy[k] == j) adjwgt[k] = maxcon;
      }
    }
    
    /* Make the constraint connections stronger */
    if(eg->connect) {
      int maxweight;
      int con,bc,bctype,sideelemtype,sidenodes;
      int j2,ind,ind2;
      int sideind[MAXNODESD1];
      
      maxweight = noknots+noelements;
      printf("Adding weight of %d for constrained nodes\n",maxweight);
      
      for(con=1;con<=eg->connect;con++) {
        bctype = eg->connectbounds[con-1];
	
        for(bc=0;bc<MAXBOUNDARIES;bc++) {    
          if(bound[bc].created == FALSE) continue;
          if(bound[bc].nosides == 0) continue;
	  
          for(i=1;i<=bound[bc].nosides;i++) {
            if(bound[bc].types[i] != bctype) continue;
	    
	    GetBoundaryElement(i,&bound[bc],data,sideind,&sideelemtype); 
	    /* GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
	       data,sideind,&sideelemtype); */

	    sidenodes = sideelemtype%100;
	    
            for(j=0;j<sidenodes;j++) {
              for(j2=0;j2<sidenodes;j2++) {
                if(j==j2) continue;
		
                ind = sideind[j]-1;
                ind2 = sideind[j2]-1;
		
                for(k=xadj[ind];k<xadj[ind+1];k++) 
                  if(adjncy[k] == ind2) adjwgt[k] = maxweight;
              }
            }
          }
        }
      }
    }
  } /* !dual */    
  
  if(metisopt == 2) {
    if(info) printf("Starting graph partitioning METIS_PartGraphRecursive.\n");  
    METIS_PartGraphRecursive(&nn,&ncon,xadj,adjncy,vwgt,&wgtflag,adjwgt,
			     &nparts,NULL,NULL,options,&edgecut,npart); 
  }
  else if( metisopt == 3 || metisopt == 4 ) {
    ncon = 1;
    wgtflag = 0;
    if(info) printf("Starting graph partitioning METIS_PartGraphKway.\n");      
    METIS_PartGraphKway(&nn,           /* number of vertices in the graph */
			&ncon,         /* number of balancing constraints */
			xadj, adjncy,  /* the adjacency structure of the graph */
			vwgt,          /* weights of the vertices */
			&wgtflag,      /* size of the vertices for computing communication */
			adjwgt,        /* weight of the edges */     
			&nparts,       /* number of partitions */
			NULL,          /* weights for each partition and constraint */
			NULL,          /* allowed load imbalance */
			options,       /* array of options */
			&edgecut,      /* the total communication volume */
			npart);        /* partition vector of the graph */
  }
  else {
    printf("Unknown Metis option %d\n",metisopt);
  }

  if(info) printf("Finished Metis graph partitioning call.\n");


  free_Ivector(adjncy,0,totcon-1);
  free_Ivector(xadj,0,nn);
  if(wgtflag == 1)  free_Ivector(adjwgt,0,totcon-1);

 skippart:

  if(!data->partitionexist) {
    data->partitionexist = TRUE;
    data->elempart = Ivector(1,noelements);
    data->nodepart = Ivector(1,noknots);
    data->nopartitions = partitions;
  }

  /* Set the partition given by Metis for each node. */
  if( dual ) {
    if( data->elemconnectexist ) {
      for(i=1;i<=noelements;i++) {
	j = data->elemconnect[i];
	if( nparts == 1 ) {
	  if(j) 
	    data->elempart[i] = j;
	  else 
	    data->elempart[i] = maxconset+1;
	} else {
	  if(j < 0) 
	    data->elempart[i] = -j;	  
	  else
	    data->elempart[i] = npart[j-1]+1+maxconset;  
	}
      }
    }
    else {
      for(i=1;i<=nn;i++) 
	data->elempart[i] = npart[i-1]+1;
    }
    PartitionNodesByElements(data,info);
  }
  else {
    for(i=1;i<=nn;i++) 
      data->nodepart[i] = npart[i-1]+1;
    PartitionElementsByNodes(data,info);

    /* Finally check that the constraint is really honored */
    if(!eg->connect) {
      int con,bc,bctype,sideelemtype,sidenodes,par;
      int ind,sideind[MAXNODESD1];
      int *sidehits,sidepartitions;

      printf("Checking connection integrity\n");
      sidehits = Ivector(1,partitions);

      for(con=1;con<=eg->connect;con++) {
	bctype = eg->connectbounds[con-1];

	for(i=1;i<=partitions;i++)
	  sidehits[i] = 0;

	for(bc=0;bc<MAXBOUNDARIES;bc++) {    
	  if(bound[bc].created == FALSE) continue;
	  if(bound[bc].nosides == 0) continue;
	
	  for(i=1;i<=bound[bc].nosides;i++) {
	    if(bound[bc].types[i] != bctype) continue;
	  
	    if(1)
	      GetBoundaryElement(i,&bound[bc],data,sideind,&sideelemtype); 
	    else
	      GetElementSide(bound[bc].parent[i],bound[bc].side[i],bound[bc].normal[i],
			     data,sideind,&sideelemtype);
	    sidenodes = sideelemtype%100;
      
	    for(j=0;j<sidenodes;j++) {
	      ind = sideind[j];
	      par = data->nodepart[ind];
	      sidehits[par] += 1;
	    }
	  }   
	}

	sidepartitions = 0;
	for(i=1;i<=partitions;i++)
	  if( sidehits[i] ) sidepartitions += 1;

	if(sidepartitions != 1) {
	  printf("PartitionMetisGraph: side %d belongs to %d partitions\n",bctype,sidepartitions);
	  bigerror("Parallel constraints might not be set!");
	}
      }
    }
  }

  if( nparts > 1 ) {
    free_Ivector(npart,0,nn-1);
  }

  if(info) {
    if( dual ) 
      printf("Successfully made a Metis partition using the dual graph.\n");
    else
      printf("Successfully made a Metis partition using the nodal graph.\n");
  }

  return(0);
}

#endif  


static void CheckPartitioning(struct FemType *data,int info)
{
  int i,j,k,partitions,part,part2,noknots,noelements,mini,maxi,sumi,hit,ind,nodesd2,elemtype;
  int *elempart, *nodepart,*elemsinpart,*nodesinpart,*sharedinpart;

  noknots = data->noknots;
  noelements = data->noelements;
  partitions = data->nopartitions;
  elemsinpart = Ivector(1,partitions);
  nodesinpart = Ivector(1,partitions);
  sharedinpart = Ivector(1,partitions);
  for(i=1;i<=partitions;i++)
    elemsinpart[i] = nodesinpart[i] = sharedinpart[i] = 0;

  if(info) printf("Checking for partitioning\n");

  /* Check that division of elements */
  elempart = data->elempart;
  j=0;
  for(i=1;i<=data->noelements;i++) {
    part = elempart[i];
    if(part < 1 || part > partitions) 
      j++;
    else 
      elemsinpart[part] += 1;
  }      
  if(j) {
    printf("Bad initial partitioning: %d elements do not belong anywhere!\n",j);
    bigerror("Can't continue with broken partitioning");
  }    

  /* Check the division of nodes */
  nodepart = data->nodepart; 
  j=0;
  for(i=1;i<=data->noknots;i++) {
    part = nodepart[i];
    if(part < 1 || part > partitions) 
      j++;
    else 
      nodesinpart[part] += 1;
  }
  
  if(j) {
    printf("Bad initial partitioning: %d nodes do not belong anywhere!\n",j);
    bigerror("Can't continue with broken partitioning");
  }

  if(data->partitiontableexists) {
    for(i=1;i<=noknots;i++) {
      part = nodepart[i];
      for(j=1;j<=data->maxpartitiontable;j++) {
	part2 = data->partitiontable[j][i];
	if(!part2) break;
	if(part != part2) sharedinpart[part2] += 1;
      }
    }
  }

  if(info) {
    printf("Information on partition bandwidth\n");
    if(partitions <= 10) {
      printf("Distribution of elements, nodes and shared nodes\n");
      printf("     %-10s %-10s %-10s %-10s\n","partition","elements","nodes","shared");
      for(i=1;i<=partitions;i++)
	printf("     %-10d %-10d %-10d %-10d\n",i,elemsinpart[i],nodesinpart[i],sharedinpart[i]);
    } 
    else {
      mini = maxi = elemsinpart[1];
      for(i=1;i<=partitions;i++) {
	mini = MIN( elemsinpart[i], mini);
	maxi = MAX( elemsinpart[i], maxi);
      }
      printf("Average %d elements with range %d in partition\n",noelements/partitions,maxi-mini);

      mini = maxi = nodesinpart[1];
      for(i=1;i<=partitions;i++) {
	mini = MIN( nodesinpart[i], mini);
	maxi = MAX( nodesinpart[i], maxi);
      }
      printf("Average %d nodes with range %d in partition\n",noknots/partitions,maxi-mini);

      sumi = 0;
      mini = maxi = sharedinpart[1];
      for(i=1;i<=partitions;i++) {
	mini = MIN( sharedinpart[i], mini);
	maxi = MAX( sharedinpart[i], maxi);
	sumi += sharedinpart[i];
      }
      printf("Average %d shared nodes with range %d in partition\n",sumi/partitions,maxi-mini);
    }
  }

  if(!data->maxpartitiontable) return;

  if(0) printf("Checking that each node in elements belongs to nodes\n");
  for(i=1;i<=data->noelements;i++) {
    part = elempart[i];
    elemtype = data->elementtypes[i];
    nodesd2 = elemtype % 100;

    for(j=0;j < nodesd2;j++) {
      ind = data->topology[i][j];
      
      hit = FALSE;
      part2 = 0;
      for(k=1;k<=data->maxpartitiontable;k++) {
	part2 = data->partitiontable[k][ind];
	if( part == part2 ) hit = TRUE;
	if(hit && !part) break;
      }
      if(!hit) {
	printf("******** Warning *******\n");
	printf("Node %d in element %d does not belong to partition %d (%d)\n",ind,i,part,part2);
	printf("elemtype = %d nodesd2 = %d\n",elemtype,nodesd2);
	for(k=0;k < nodesd2;k++) 
	  printf("ind[%d] = %d\n",k,data->topology[i][k]);
      }
    }
  }

  if(0) printf("Checking that each node in partition is shown in partition list\n");
  for(i=1;i<=data->noknots;i++) {
    part = nodepart[i];
    
    hit = FALSE;
    for(j=1;j<=data->maxpartitiontable;j++) {
      part2 = data->partitiontable[j][i];
      if( part == part2 ) hit = TRUE;
	if(hit && !part) break;
    }
    if(!hit) {
      printf("***** Node %d in partition %d is not in partition list\n",i,part);
    }
  }
}


int OptimizePartitioningAtBoundary(struct FemType *data,struct BoundaryType *bound,int info)
{
  int i,j,k,l,boundaryelems,ind,hit1,hit2,fix1,fix2;
  int dompart,part1,part2,newmam,mam1,mam2,nodesd2;
  int *alteredparent;

  if(!data->partitionexist) {
    printf("OptimizePartitioningAtBoundary: this should be called only after partitioning\n");
    bigerror("Optimization not performed!");
  }

  if(info) printf("Optimizing the partitioning at boundaries.\n");

  /* Memorize the original parent */
  alteredparent = Ivector(1,data->noelements);
  for(i=1;i<=data->noelements;i++)
    alteredparent[i] = data->elempart[i];

  
  /* Set the secondary parent to be a parent also because we want all 
     internal BCs to be within the same partition. 
     Also set the nodes of the altered elements to be in the desired partition. */
  k = 0;
  do {
    k++;
    boundaryelems = 0;

    for(j=0;j < MAXBOUNDARIES;j++) {
      if(!bound[j].created) continue;
      for(i=1; i <= bound[j].nosides; i++) {
	if(bound[j].ediscont)
	  if(bound[j].discont[i]) continue;

	mam1 = bound[j].parent[i];
	mam2 = abs(bound[j].parent2[i]);
	if(!mam1 || !mam2) continue;

	part1 = data->elempart[mam1];
	part2 = data->elempart[mam2];
	if(part1 == part2) continue;
	       
	/* Check if both nodes are enforced to be fixed */
	if( data->elemconnectexist ) {
	  fix1 = ( data->elemconnect[mam1] < 0 );
	  fix2 = ( data->elemconnect[mam2] < 0 );
	  if( fix1 && fix2 ) continue;
	}
	else {
	  fix1 = fix2 = 0;
	}

	/* The first iterations check which parents is ruling 
	   thereafter choose pragmatically the other to overcome
	   oscillating solutions. */
	if(fix1 || fix2 ) {
	}
	if(k < 5) {
	  hit1 = hit2 = 0;
	  nodesd2 = data->elementtypes[mam1] % 100;
	  for(l=0;l < nodesd2;l++) {
	    ind = data->topology[mam1][l];
	    if(data->nodepart[ind] == part1) hit1++;
	    if(data->nodepart[ind] == part2) hit2++;
	  }
	  nodesd2 = data->elementtypes[mam2] % 100;    
	  for(l=0;l < nodesd2;l++) {
	    ind = data->topology[mam2][l];
	    if(data->nodepart[ind] == part1) hit1++;
	    if(data->nodepart[ind] == part2) hit2++;
	  }	  
	  fix1 = ( hit1 >= hit2 );
	  fix2 = !fix1;
	} 
	else {
	  fix1 = TRUE;
	  fix2 = FALSE;
	}   

	/* Make the more ruling parent dominate the whole boundary */
	if(fix1) {
	  dompart = part1;
	  newmam = mam2;
	}
	else {
	  dompart = part2;
	  newmam = mam1;
	}
	
	data->elempart[newmam] = dompart;
	boundaryelems++;	    

	/* Move the ownership of all nodes to the leading partition */
	if(0) {
	  nodesd2 =  data->elementtypes[newmam] % 100;
	  for(l=0;l < nodesd2;l++) {
	    ind = data->topology[newmam][l];
	    data->nodepart[ind] = dompart;
	  }
	}
	else if(0) {
	  nodesd2 = data->elementtypes[mam1] % 100;
	  for(l=0;l < nodesd2;l++) {
	    ind = data->topology[mam1][l];
	    data->nodepart[ind] = dompart;
	  }	    
	  nodesd2 = data->elementtypes[mam2] % 100;
	  for(l=0;l < nodesd2;l++) {
	    ind = data->topology[mam2][l];
	    data->nodepart[ind] = dompart;
	  }	    
	}

      }
    }
    if(info && boundaryelems) 
      printf("Round %d: %d bulk elements with BCs removed from interface.\n",
	     k,boundaryelems);
  } while(boundaryelems && k < 10);


  j = 0;
  for(i=1;i<=data->noelements;i++) {
    if(alteredparent[i] == data->elempart[i]) 
      alteredparent[i] = 0;
    else 
      j++;
  }
  free_Ivector(alteredparent,1,data->noelements);
  if(info) printf("Ownership of %d parents was changed at BCs\n",j); 



  /* Remove the negative secondary parents that were only used to optimize the partitioning */ 
  for(j=0;j < MAXBOUNDARIES;j++) {
    if(!bound[j].created) continue;
    for(i=1; i <= bound[j].nosides; i++) 
      bound[j].parent2[i] = MAX(0, bound[j].parent2[i]);
  }

  if(0) printf("The partitioning at discontinuous gaps was optimized.\n");
  return(0);
}


static void Levelize(int n,int level,int *maxlevel,int *levels,int *rows,int *cols,int *done)
{    
   int j,k;

   levels[n] = level;
   done[n] = TRUE;
   *maxlevel = MAX( *maxlevel,level );

   for( j=rows[n]; j<rows[n+1]; j++ ) {
      k = cols[j];
      if ( !done[k] ) Levelize(k,level+1,maxlevel,levels,rows,cols,done);
   }
}


static int RenumberCuthillMckee( int nrows, int *rows, int *cols, int *iperm )
{
  int i,j,k,n,startn,mindegree,maxlevel,newroot,bw_bef,bw_aft;
  int *level,*degree,*done;
  
  done   = Ivector(0,nrows-1);
  level  = Ivector(0,nrows-1);
  degree = Ivector(0,nrows-1);
  
  bw_bef = 0;
  for(i=0; i<nrows; i++ )
    {
      for( j=rows[i]; j<rows[i+1]; j++ )
	bw_bef = MAX( bw_bef, ABS(cols[j]-i)+1 );
      degree[i] = rows[i+1]-rows[i];
    }
  printf( "RenumberCuthillMckee: Bandwidth before: %d\n", bw_bef );
  
  startn = 0;
  mindegree = degree[startn];
  for( i=0; i<nrows; i++ ) {
    if ( degree[i] < mindegree ) {
      startn = i;
      mindegree = degree[i];
    }
    level[i] = 0;
  }
  
  maxlevel = 0;
  for( i=0; i<nrows; i++ ) done[i]=FALSE;
  
  Levelize( startn,0,&maxlevel,level,rows,cols,done );
  
  newroot = TRUE;
  while(newroot) {
    newroot = FALSE;
    mindegree = degree[startn];
    k = startn;
    
    for( i=0; i<nrows; i++ ) {
      if ( level[i] == maxlevel ) {
	if ( degree[i] < mindegree ) {
	  k = i;
	  mindegree = degree[i];
	}
      }
    }
    
    if ( k != startn ) {
      j = maxlevel;
      maxlevel = 0;
      for(i=0; i<nrows; i++ ) done[i]=FALSE;
      
      Levelize( k,0,&maxlevel,level,rows,cols,done );
      
      if ( j > maxlevel ) {
	newroot = TRUE;
	startn = j;
      }
    }
  }
  
  for(i=0; i<nrows; i++ ) done[i]=-1,iperm[i]=-1;
  
  done[0]=startn;
  iperm[startn]=0;
  i=1;
  
  for( j=0; j<nrows; j++ ) {
    if ( done[j]<0 ) {
      for( k=0; k<nrows; k++ ) {
	if ( iperm[k]<0 ) {
	  done[i]=k;
	  iperm[k]=i;
	  i++;
	  break;
	}
      }
    }
    
    for( k=rows[done[j]]; k<rows[done[j]+1]; k++) { 
      n = cols[k];
      if ( iperm[n]<0 ) {
	done[i] = n;
	iperm[n] = i;
	i++;
      }
    }
  }
  
  for( i=0; i<nrows; i++ )
    iperm[done[i]] = nrows-1-i;
  
  bw_aft = 0;
  for(i=0; i<nrows; i++ )
    for( j=rows[i]; j<rows[i+1]; j++ )
      bw_aft = MAX( bw_aft, ABS(iperm[cols[j]]-iperm[i])+1 );
  
  printf( "RenumberCuthillMckee: Bandwidth after: %d\n", bw_aft );
  
  free_Ivector(level,0,nrows-1);
  free_Ivector(done,0,nrows-1);
  free_Ivector(degree,0,nrows-1);

  return bw_aft < bw_bef;
}



static void RenumberPartitions(struct FemType *data,int info)
/* Minimize bandwidth of partition indexing. This could be favourable if the
   communication between neighbouring partitions is faster than between 
   partitions with larger partition index. Probably there is not much 
   difference for normal hardware configurations. */
{
  int i,j,k,l,m,hit,con,totcon,noelements,noknots,partitions;
  int maxneededtimes,totneededtimes;
  int part,part1,part2,bw_reduced;
  int *nodepart,*elempart;
  int *perm;
  int *xadj,*adjncy;
  int *partparttable[MAXCONNECTIONS];


  if(info) printf("Renumbering partitions to minimize bandwidth.\n");  

  partitions = data->nopartitions;
  noelements = data->noelements;
  noknots = data->noknots;
  maxneededtimes = data->maxpartitiontable;


  /* Make the partition-partition list from the node-partition list */
  totneededtimes = 0;
  totcon = 0;
  for(i=1;i<=noknots;i++) {
    if(data->partitiontable[2][i] == 0) continue;
    for(j=1;j<=maxneededtimes;j++) {
      part1 =  data->partitiontable[j][i];
      if(!part1) break;

      for(k=1;k<=maxneededtimes;k++) {
	if(k==j) continue;
	part2 = data->partitiontable[k][i];
	if(!part2) break;

	hit = 0;
	for(l=1;l<=totneededtimes;l++) { 
	  if(partparttable[l][part1] == part2) {
	    hit = -1;
	    break;
	  }
	  else if(partparttable[l][part1] == 0) {
	    totcon++;
	    partparttable[l][part1] = part2;
	    hit = 1;
	    break;
	  }
	}
	if(!hit) {
	  totneededtimes++;
	  partparttable[totneededtimes] = Ivector(1,partitions);
	  for(m=1;m<=partitions;m++)
	    partparttable[totneededtimes][m] = 0;
	  partparttable[totneededtimes][part1] = part2;
	  totcon++;
	}
      }
    }
  }

  if(info) {
    printf("There are %d connections altogether\n",totcon);
    printf("There are %.3f connections between partitions in average\n",1.0*totcon/partitions);
  }


  xadj = Ivector(0,partitions);
  adjncy = Ivector(0,totcon-1);
  for(i=0;i<totcon;i++) 
    adjncy[i] = 0;

  totcon = 0;
  for(i=1;i<=partitions;i++) {
    xadj[i-1] = totcon;
    for(j=1;j<=totneededtimes;j++) {
      con = partparttable[j][i];
      if(!con) continue;
      adjncy[totcon] = con-1;
      totcon++;
    }
  }    
  xadj[partitions] = totcon;

  perm = Ivector(0,partitions-1);
  bw_reduced = RenumberCuthillMckee( partitions, xadj, adjncy, perm );


  /* Print the new order of partitions */
  if(0 && info) {
    printf( "Partition order after Cuthill-McKee bandwidth optimization: \n" );
    for(i=0;i<partitions;i++)
      printf("old=%d new=%d\n",i,perm[i] );
  }

  /* Use the renumbering or not */
  if(bw_reduced) {
    if(info) printf("Successful ordering: moving partitions to new positions\n");
    nodepart = data->nodepart;
    elempart = data->elempart;
    for(i=1;i<=noelements;i++) 
      elempart[i] = perm[elempart[i]-1]+1;
    for(i=1;i<=noknots;i++)
      nodepart[i] = perm[nodepart[i]-1]+1;
    for(i=1;i<=noknots;i++) {
      for(j=1;j<=maxneededtimes;j++) {
	part = data->partitiontable[j][i];
	if(!part) break;
	data->partitiontable[j][i] = perm[part-1]+1;
      }
    }
  }

  for(i=1;i<=totneededtimes;i++)
    free_Ivector(partparttable[i],1,partitions);
  free_Ivector(xadj,0,partitions);
  free_Ivector(adjncy,0,totcon-1);


}


static int CheckSharedDeviation(int *neededvector,int partitions,int info)
{
  int i,minshared,maxshared,dshared;
  Real sumshared, sumshared2, varshared, needed, det, ave;
    
  sumshared = sumshared2 = 0.0;
  minshared = maxshared = neededvector[1];
  for(i=1;i<=partitions;i++) {
    needed = 1.0 * neededvector[i];
    sumshared += needed;
    sumshared2 += (needed * needed);
    maxshared = MAX(maxshared, neededvector[i]);
    minshared = MIN(minshared, neededvector[i]);
  }
  
  dshared = maxshared - minshared;  

  if(info) {
    ave = sumshared / partitions;
    det = 1.0*sumshared2 / partitions - 1.0*(sumshared/partitions)*(sumshared / partitions);
    varshared = sqrt( det );

    printf("Average number of elements in partition %.3le\n",ave);
    printf("Maximum deviation in ownership %d\n",dshared);
    printf("Average deviation in ownership %.3le\n",varshared);      
    printf("Average relative deviation %.2lf %%\n",100.0*varshared/ave);
  }

  return(dshared);
}  




int OptimizePartitioning(struct FemType *data,struct BoundaryType *bound,int noopt,
			 int partbw, int info)
/* Optimize partitioning of elements so that each partition has as closely as possible the 
   desired amount of elements. Also, ir requested, check that there are no odd couplings within 
   elements that are not directly present at the given partition. It is a bit unclear 
   to which extent these checks are needed in the current Elmer version. */
{
  int i,j,k,l,n,m,noelements,partitions,ind,hit;
  int noknots,dshared,dshared0;
  int *elempart,*nodepart,sharings,maxshared;
  int nodesd2,maxneededtimes,*probnodes=NULL,optimize;
  int *neededvector;
  int somethingdone = 0;
  int *elemparts = NULL,*invelemparts = NULL;
  int **knows = NULL;

  if(!data->partitionexist) {
    printf("OptimizePartitioning: this should be called only after partitioning\n");
    bigerror("Optimization not performed!");
  }

  noknots = data->noknots;
  noelements = data->noelements;
  partitions = data->nopartitions;
  if(info) printf("Optimizing for %d partitions\n", partitions);

  elempart = data->elempart;
  nodepart = data->nodepart; 
  maxshared = 0;

  if( partitions < 2 ) {
    printf("OptimizePartitioning: does not make sense for %d partitions\n",partitions);
    bigerror("Optimization not performed!"); 
  }

  /* Create a table showing to which partitions nodes belong to */
  CreatePartitionTable(data,info);
  maxneededtimes = data->maxpartitiontable;

  /* Renumber the bandwidth of partition-partition connections */
  if(partbw) RenumberPartitions(data,info);

  /* Check partitioning after table is created for the first time */
  printf("Checking partitioning before optimization\n");
  CheckPartitioning(data,info);

  /* Calculate how many nodes is owned by each partition */
  neededvector = Ivector(1,partitions);  
  for(i=1;i<=partitions;i++) 
    neededvector[i] = 0;    
  for(i=1;i<=noknots;i++) 
    neededvector[nodepart[i]] += 1;
   
  if(!noopt) {
    optimize = 1;
    probnodes = Ivector(1,noknots);
    for(i=1;i<=noknots;i++)
      probnodes[i] = 0;
    printf("Applying aggressive optimization for load balancing\n");
  }
  else {
    optimize = 1;
  }

 optimizeownership:

  dshared = CheckSharedDeviation(neededvector,partitions,info);

  if(!noopt) {    
    
    /* Distribute the shared nodes as evenly as possible. */

    int nochanges, maxrounds, dtarget;

    maxrounds = 5;
    dtarget = 3;
    nochanges = 0;

    for(n=1;n<=maxrounds;n++) {
      nochanges = 0;

      for(i=1;i<=noknots;i++) {      
	ind = i;

	/* owner partition may only be changed it there a are two of them */
	k = data->partitiontable[2][ind];
	if(!k) continue;
	
	/* only apply the switch to cases with exactly two partitions 
	   to avoid the nasty multiply coupled nodes. */
	if(maxneededtimes > 2) 
	  if(data->partitiontable[3][ind]) continue;
	
	/* Do not change the ownership of nodes that cause topological problems */
	if(probnodes[ind]) continue;	  
	
	j = data->partitiontable[1][ind];

	/* Switch the owner to the smaller owner group if possible */
	if(abs(neededvector[j] - neededvector[k]) < dtarget ) continue;

	if(neededvector[j] < neededvector[k] && nodepart[ind] == k) {
	  nochanges++;
	  neededvector[j] += 1;
	  neededvector[k] -= 1;
	  nodepart[ind] = j;
	}
	else if(neededvector[k] < neededvector[j] && nodepart[ind] == j) {
	  nochanges++;
	  neededvector[k] += 1;
	  neededvector[j] -= 1;
	  nodepart[ind] = k;
	}
      }
      if(info && nochanges) printf("Changed the ownership of %d nodes\n",nochanges);
      
      dshared0 = dshared;
      dshared = CheckSharedDeviation(neededvector,partitions,info);
      
      if(dshared >= dshared0) break;
      somethingdone = somethingdone + nochanges;
    }
  }


  /* optimizesharing: */
  
  if(info) printf("Checking for problematic sharings\n"); 
  m = 0;

  if(partitions > 2) {  
    if(info) printf("Optimizing sharing for %d partitions\n", partitions);

    do {
      
      int i1,i2,e1,e2,owners;
      
      m++;
      sharings = 0;
      i1 = i2 = 0;
      e1 = e2 = 0;

      /* There cannot be more partitions sharing a node than there are
	 node in the elements. Elementtype 827 has 27 nodes. */
      maxshared = 27;
      if(m == 1 && optimize == 1) {
	elemparts = Ivector(1,partitions);
	invelemparts = Ivector(1,maxshared);
	knows = Imatrix(1,maxshared,1,maxshared);
      }
      
      for(j=1;j<=partitions;j++) 
	elemparts[j] = 0;
      
      for(j=1;j<=maxshared;j++) 
	invelemparts[j] = 0;
      
      for(j=1;j<=maxshared;j++) 
	for(k=1;k<=maxshared;k++)
	  knows[j][k] = 0;
      
      for(i=1;i<=noelements;i++) {
	int ownpart;
	
	owners = 0;
	nodesd2 = data->elementtypes[i] % 100;
	ownpart = FALSE;
	
	/* Check the number of owners in an element */
	for(j=0;j < nodesd2;j++) {
	  ind = data->topology[i][j];
	  k = nodepart[ind];
	  
	  /* Mark if the element partition is one of the owners */
	  if( k == elempart[i]) ownpart = TRUE;
	  if(!elemparts[k]) {
	    owners++;
	    elemparts[k] = owners;
	    invelemparts[owners] = k;
	  }
	}
	
	/* One strange owner is still ok. */
	if(owners - ownpart <= 1) {
	  /* Nullify the elemparts vector */
	  for(j=1;j<=owners;j++) {
	    k = invelemparts[j];
	    elemparts[k] = 0;
	  }
	  continue;
	}
	
	/* Check which of the partitions are related by a common node */
	for(j=0;j < nodesd2;j++) {
	  ind = data->topology[i][j];
	  for(l=1;l<=maxneededtimes;l++) {
	    e1 = data->partitiontable[l][ind];
	    if(!e1) break;
	    e1 = elemparts[e1];
	    if(!e1) continue;
	    for(k=l+1;k<=maxneededtimes;k++) {
	      e2 = data->partitiontable[k][ind];
	      if(!e2) break;
	      e2 = elemparts[e2];
	      if(!e2) continue;
	      knows[e1][e2] = knows[e2][e1] = TRUE;
	    }
	  }
	}    
	
	/* Check if there are more complex relations:
	   i.e. two partitions are joined at an element but not at the same node. */
	hit = FALSE;
	for(j=1;j<=owners;j++)
	  for(k=j+1;k<=owners;k++) {
	    if(!knows[j][k]) {
	      hit += 1;
	      i1 = invelemparts[j];
	      i2 = invelemparts[k];
	    }
	  }
	/* Nullify the elemparts vector */
	for(j=1;j<=owners;j++) {
	  k = invelemparts[j];
	  elemparts[k] = 0;
	}
	
	/* Nullify the knows matrix */
	for(j=1;j <= owners;j++) 
	  for(k=1;k <= owners;k++) 
	    knows[j][k] = FALSE;
	
	if(hit) {
	  e1 = e2 = 0;
	  
	  if(info) {
	    if( hit + m > 2 ) printf("Partitions %d and %d in element %d (%d owners) oddly related %d times\n",
				     i1,i2,i,owners,hit);
	  }
	  
	  /* Count the number of nodes with wrong parents */
	  for(j=0;j < nodesd2;j++) {
	    ind = data->topology[i][j];
	    for(l=1;l<=maxneededtimes;l++) {
	      k = data->partitiontable[l][ind];
	      if(k == 0) break;
	      if(k == i1) e1++;
	      if(k == i2) e2++;
	    }
	  }
	  
	  /* Change the owner of those with less sharings */
	  for(j=0;j < nodesd2;j++) {
	    ind = data->topology[i][j];
	    k = nodepart[ind];
	    if((k == i1 && e1 < e2) || (k == i2 && e1 >= e2)) {
	      if(!noopt) probnodes[ind] += 1;
	      nodepart[ind] = elempart[i];
	      neededvector[elempart[i]] += 1;
	      neededvector[k] -= 1;
	    }
	  }	
	  sharings++;
	}
      }
      
      somethingdone += sharings;
      if(info && sharings) printf("Changed the ownership of %d nodes\n",sharings);
      
    } while (sharings > 0 && m < 3);
  
    if(info) {
      if(sharings) 
	printf("%d problematic sharings may still exist\n",sharings);
      else 
	printf("There shouldn't be any problematic sharings, knock, knock...\n");
    }
  }

  /* This seems to work also iteratively */
  if(!noopt && m+n > 10 && optimize < 50) {
    optimize++;
    printf("Performing ownership optimization round %d\n",optimize);
    goto optimizeownership;
  }

  free_Ivector(neededvector,1,partitions);

  if(!noopt) {
    free_Ivector(probnodes,1,noknots);   
    if(partitions > 2) {
      free_Ivector(elemparts,1,partitions);
      free_Ivector(invelemparts,1,maxshared);
      free_Imatrix(knows,1,maxshared,1,maxshared);
    }
  }


  if(somethingdone) {    
    if( info ) {
      printf("The partitioning was optimized: %d\n",somethingdone); 
      printf("Checking partitioning after optimization\n");
    }
    CheckPartitioning(data,info);
  }
  else {
    if(info) printf("Partitioning was not altered\n");
  }

  return(0);
}


int SaveElmerInputPartitioned(struct FemType *data,struct BoundaryType *bound,
			      char *prefix,int decimals,int *parthalo,int indirect,
			      int parthypre,int subparts,int nooverwrite, int info)
/* Saves the mesh in a form that may be used as input 
   in Elmer calculations in parallel platforms. 
   */
{
  int noknots,noelements,sumsides,partitions,hit,found,parent,parent2,maxnosides;
  int nodesd2,nodesd1,discont,maxelemtype,minelemtype,sidehits,elemsides,side,bctype;
  int part,otherpart,part2,part3,elemtype,sideelemtype,*needednodes,*neededtwice;
  int **bulktypes,*sidetypes,tottypes,splitsides;
  int i,j,k,l,l2,l3,m,n,ind,ind2,sideind[MAXNODESD1],elemhit[MAXNODESD2];
  char filename[MAXFILESIZE],outstyle[MAXFILESIZE];
  char directoryname[MAXFILESIZE],subdirectoryname[MAXFILESIZE];
  int *neededtimes,*elempart,*elementsinpart,*indirectinpart,*sidesinpart;
  int maxneededtimes,indirecttype,bcneeded,trueparent,trueparent2,*ownerpart;
  int *sharednodes,*ownnodes,reorder,*order=NULL,*invorder=NULL;
  int *bcnode,*bcnodedummy,*elementhalo,*neededtimes2;
  int partstart,partfin,filesetsize,nofile,nofile2,nobcnodes,hasbcnodes,halomode;
  int halobulkelems,halobcs,savethis,fail=0,cdstat,immersed,halocopies,anyparthalo;

  FILE *out,*outfiles[MAXPARTITIONS+1];
  int sumelementsinpart,sumownnodes,sumsharednodes,sumsidesinpart,sumindirect;


  
  if(info) {
    printf("Saving Elmer mesh in partitioned format\n");
    if( subparts ) printf("There are %d subpartitions\n",subparts);
  }

  if(!data->created) {
    printf("You tried to save points that were never created.\n");
    bigerror("No Elmer mesh files saved!");
  }

  partitions = data->nopartitions;
  if(!partitions) {
    printf("Tried to save partitioned format without partitions!\n");
    bigerror("No Elmer mesh files saved!");
  }

  anyparthalo = FALSE; 
  for(i=0;i<MAXHALOMODES;i++) {
    if( parthalo[i] ) {
      printf("Saving halo elements mode active: %d\n",i);
      anyparthalo = TRUE;
    }
  }    
  
  if( subparts < 1 && parthalo[3] ) {
    printf("There can be no layer halo since there are no layers!\n");
    bigerror("No Elmer mesh files saved!");
  }


  elempart = data->elempart;
  ownerpart = data->nodepart;
  noelements = data->noelements;
  noknots = data->noknots;

  minelemtype = 101;
  maxelemtype = GetMaxElementType(data);
  indirecttype = 0;
  halobulkelems = 0;

  needednodes = Ivector(1,partitions);
  neededtwice = Ivector(1,partitions);
  sharednodes = Ivector(1,partitions);
  ownnodes = Ivector(1,partitions);
  sidetypes = Ivector(minelemtype,maxelemtype);
  bulktypes =  Imatrix(1,partitions,minelemtype,maxelemtype);
  bcnodedummy = Ivector(1,noknots);

  /* Order the nodes so that the different partitions have a continuous interval of nodes.
     This information is used only just before the saving of node indexes in each instance. 
     This feature was coded for collaboration with Hypre library that assumes this. */
  reorder = parthypre;
  if(reorder) {
    order = Ivector(1,noknots);
    k = 0;
    for(j=1;j<=partitions;j++)
      for(i=1; i <= noknots; i++) 
	if(ownerpart[i] == j) {
	  k++;
	  order[i] = k; 
	}
    invorder = Ivector(1,noknots);
    for(i=1;i<=noknots;i++) 
      invorder[order[i]] = i;
    if(info) printf("Created continuous parallel ordering\n");
  } 

  /* Mark the nodes that are on some boundary in order to create boundary halos */
  bcnode = Ivector(1,noknots);
  for(i=1;i<=noknots;i++)
    bcnode[i] = FALSE;
  for(j=0;j < MAXBOUNDARIES;j++) {
    for(i=1; i <= bound[j].nosides; i++) {		
      GetBoundaryElement(i,&bound[j],data,sideind,&sideelemtype); 
      nodesd1 = sideelemtype%100;
      for(l=0;l<nodesd1;l++) 
	bcnode[sideind[l]] = 1;
    }
  }
  nobcnodes = 0;
  for(i=1;i<=noknots;i++) {
    if( bcnode[i] == 1) nobcnodes += 1;
  }
  if(info) printf("Number of boundary nodes at the boundary: %d\n",nobcnodes);
  
  
  for(i=1;i<=noknots;i++)
    bcnodedummy[i] = 0;
  for(i=1;i<=noelements;i++) {
    k = data->material[i];
    elemtype = data->elementtypes[i];
    nodesd2 = elemtype%100;

    for(j=0;j < nodesd2;j++) {
      ind = data->topology[i][j];
      if( bcnode[ind] == 1) continue;
      if( bcnodedummy[ind] ) {
	if( bcnodedummy[ind] != k ) bcnode[ind] = 2;
      } else {
	bcnodedummy[ind] = k;
      }
    }
  }
  nobcnodes = 0;
  for(i=1;i<=noknots;i++) {
    if( bcnode[i] == 2) nobcnodes += 1;
  }
  if(info) printf("Number of additional interface nodes: %d\n",nobcnodes);

  
  sprintf(directoryname,"%s",prefix);
  sprintf(subdirectoryname,"%s.%d","partitioning",partitions);


#ifdef MINGW32
  mkdir(directoryname);
#else
  fail = mkdir(directoryname,0750);
#endif
  if(info && !fail) printf("Created mesh directory: %s\n",directoryname);
  fail = chdir(directoryname);
#ifdef MINGW32
  fail = mkdir(subdirectoryname);
#else
  fail = mkdir(subdirectoryname,0750);
#endif
  if(fail) {
    if(info) printf("Reusing existing subdirectory: %s\n",subdirectoryname);
    if(nooverwrite) {
      printf("Mesh seems to already exist, writing is cancelled!\n"); 
      return(0);
    }
  }
  else {
    if(info) printf("Created subdirectory: %s\n",subdirectoryname);
  }
  cdstat = chdir(subdirectoryname);



  if(info) printf("Saving mesh in parallel ElmerSolver format to directory %s/%s.\n",
		  directoryname,subdirectoryname);

  filesetsize = MAXPARTITIONS;
  if(partitions > filesetsize) 
    if(info) printf("Saving %d partitions in maximum sets of %d\n",partitions,filesetsize);

  elementsinpart = Ivector(1,partitions);
  indirectinpart = Ivector(1,partitions);
  sidesinpart = Ivector(1,partitions);
  elementhalo = Ivector(1,partitions);
  for(i=1;i<=partitions;i++)
    elementsinpart[i] = indirectinpart[i] = sidesinpart[i] = elementhalo[i] = 0;

  for(j=1;j<=partitions;j++)
    for(i=minelemtype;i<=maxelemtype;i++)
      bulktypes[j][i] = 0;

  /* Compute how many times a node may be needed at maximum */
  maxneededtimes = data->maxpartitiontable;
  neededtimes = Ivector(1,noknots);
  for(i=1;i<=noknots;i++) {
    neededtimes[i] = 0;
    for(j=1;j<=maxneededtimes;j++)
      if(data->partitiontable[j][i]) neededtimes[i] += 1;
  }
  if(info) printf("Nodes belong to %d partitions in maximum\n",maxneededtimes);


  /*********** part.n.elements *********************/
  /* Save elements in all partitions and where they are needed */

  
  partstart = 1;
  partfin = MIN( partitions, filesetsize );

 next_elements_set:

  for(part=partstart;part<=partfin;part++) {
    sprintf(filename,"%s.%d.%s","part",part,"elements");
    nofile = part - partstart + 1;
    outfiles[nofile] = fopen(filename,"w");
  }

  for(i=1;i<=noelements;i++) {
    part = elempart[i];

    elemtype = data->elementtypes[i];
    nodesd2 = elemtype%100;

    if(part >= partstart && part <= partfin) {
      nofile = part - partstart + 1;
      bulktypes[part][elemtype] += 1;
      elementsinpart[part] += 1;
    
      fprintf(outfiles[nofile],"%d %d %d ",i,data->material[i],elemtype);

      for(j=0;j < nodesd2;j++) {
	ind = data->topology[i][j];
	if(reorder) ind = order[ind];
	fprintf(outfiles[nofile],"%d ",ind);
      }
      fprintf(outfiles[nofile],"\n");    
    }

    
    if(!anyparthalo) continue;
    halocopies = 0;



    for( halomode = 1; halomode <= 4; halomode++) {
      
      if(!parthalo[halomode]) continue;
      
      
      /* This creates a simple halo when the elements have been partitioned such
	 that they are in number of "subparts" intervals of z coordinate. */
      if( halomode == 3 ) {
	if( part <= subparts ) {
	  int leftright;
	  for(leftright=-1;leftright <=1;leftright += 2) {
	    part2 = part+leftright;	  
	    if( part2 < 1 || part2 > subparts ) continue;		
	    elementhalo[part2] = i;
	    halocopies += 1;
	  }
	}
      }
      
      /* For the following halos we study the shared and boundary nodes. */
      otherpart = 0;
      hasbcnodes = FALSE;
      for(j=0;j < nodesd2;j++) {
	ind = data->topology[i][j];
	if(neededtimes[ind] > 1) otherpart++;
	if(bcnode[ind]) hasbcnodes = TRUE;
      }
            
      /* Classical halo needed by discontinuous Galerkin.
	 Elements that are attached by joined face element are added to halo */
      if( halomode == 1 && otherpart ) {
	/* If the saving of halo is requested check it for elements which have at least 
	   two nodes in shared partitions. First make this quick test. */
      
	elemsides = elemtype / 100;
	if(elemsides == 8) {
	  immersed = (otherpart >= 4);
	  elemsides = 6;
	}
	else if(elemsides == 7) {
	  immersed = (otherpart >= 3);
	  elemsides = 5;
	}
	else if(elemsides == 6) {
	  immersed = (otherpart >= 3);
	  elemsides = 5;
	}
	else if(elemsides == 5) {
	  immersed = (otherpart >= 3);
	  elemsides = 4;
	}      
	else {
	  immersed = (otherpart >= 2);
	}
	
	if( immersed ) {
      
	  /* In order for the halo to be present the element should have a boundary 
	     fully immersed in the other partition. This test takes more time. */
	
	  for(side=0;side<elemsides;side++) {

	    GetElementSide(i,side,1,data,&sideind[0],&sideelemtype);

	    /* Because every node must be on the boundary use the 1st index as the 
	       first test */
	    for(l=1;l<=neededtimes[sideind[0]];l++) {
	      part2 = data->partitiontable[l][sideind[0]];
	    
	      /* We did already save this in partition part */
	      if(part2 == part) continue;
	      if(elementhalo[part2] == i) continue;	    
	      if(part2 < partstart || part2 > partfin) continue;

	      /* Now check that also the other nodes on the face are at the interface */
	      sidehits = 1;	 	  
	      for(k=1;k<sideelemtype%100;k++) {
		for(l2=1;l2<=neededtimes[sideind[k]];l2++) {
		  part3 = data->partitiontable[l2][sideind[k]];
		  if(part2 == part3) sidehits++;
		}
	      }
	      if( sidehits <  sideelemtype % 100 ) continue;
	      
	      if(0) printf("Adding halo for partition %d and element %d\n",part2,i);
	      
	      /* Remember that this element is saved for this partition */
	      elementhalo[part2] = i;
	      halocopies += 1;
	    }
	  }
	}            
      }
      
      /* This creates a halo that includes all elements attached to boundary or interface nodes */
      if( halomode == 2 && hasbcnodes ) {
	for(j=0;j < nodesd2;j++) {
	  ind = data->topology[i][j];
	    
	  /* If the node is not on the boundary then cycle */
	  if( !bcnode[ind] ) continue;
	    
	  /* Check to which partitions this element is associated with */
	  for(k=1;k<=neededtimes[ind];k++) {
	    part2 = data->partitiontable[k][ind]; 
	      
	    /* This element is already saved to its primary partition */
	    if( part == part2 ) continue;
	    if( elementhalo[part2] == i) continue;
	    if( part2 < partstart || part2 > partfin ) continue;  	  
	      
	    if(0) printf("Adding bc halo for partition %d and element %d\n",part2,i);
	      
	    /* Remember that this element is saved for this partition */
	    elementhalo[part2] = i;
	    halocopies += 1;
	  }
	}
      }
      
      /* This greedy routine makes the halo even if there is just one node in the shared boundary. 
	 Currently not active. */
      if( halomode == 4 && otherpart ) {	
	for(j=0;j < nodesd2;j++) {
	  ind = data->topology[i][j];
	  if(neededtimes[ind] == 1) continue;
	  
	  for(l=1;l<=neededtimes[ind];l++) {
	    part2 = data->partitiontable[l][ind];
	    
	    /* We did already save this in partition part */
	    if(part2 == part) continue;	    
	    if( elementhalo[part2] == i) continue;
	    if(part2 < partstart || part2 > partfin) continue;
	    
	    /* Remember that this element is saved for this partition */
	    elementhalo[part2] = i;
	    halocopies += 1;
	  }
	}
      }
    }

      
    if( halocopies ) {
      halobulkelems += halocopies;
      
      for(part2=partstart; part2 <= partfin; part2++) {
	if( part == part2) continue;
	if( elementhalo[part2] != i ) continue;

	nofile2 = part2 - partstart + 1;
	fprintf(outfiles[nofile2],"%d/%d %d %d ",i,part,data->material[i],elemtype);
	      
	for(j=0;j < nodesd2;j++) {
	  ind = data->topology[i][j];
	  if(reorder) ind = order[ind];
	  fprintf(outfiles[nofile2],"%d ",ind);
	}
	fprintf(outfiles[nofile2],"\n");    	    
	
	bulktypes[part2][elemtype] += 1;
	elementsinpart[part2] += 1;	
	
	/* Add the halo on-the-fly to the partitiontable of the nodes */	    
	for(j=0;j < nodesd2;j++) {
	  ind = data->topology[i][j];
	  hit = FALSE;
	  for(k=1;k<=maxneededtimes;k++) {
	    part3 = data->partitiontable[k][ind];
	    if(part3 == part2) hit = TRUE;
	    if(!part3) break;
	  }
	  if(!hit) {
	    if(k <= maxneededtimes) {
	      data->partitiontable[k][ind] = part2;
	    } 
	    else {
	      maxneededtimes++;
	      if(0) printf("Allocating new column %d in partitiontable\n",k);
	      data->partitiontable[k] = Ivector(1,noknots);
	      for(m=1;m<=noknots;m++)
		data->partitiontable[k][m] = 0;
	      data->partitiontable[k][ind] = part2;
	    }
	  }
	}
      }
    } /* halocopies */
    
  } /* noelements */

  if( anyparthalo ) {
    if(info) printf("There are %d bulk elements in the halo.\n",halobulkelems);
  }
	  

  for(part=partstart;part<=partfin;part++) {
    nofile = part - partstart + 1;
    fclose(outfiles[nofile]);
  }
  if(partfin < partitions) {
    partstart = partfin + 1;
    partfin = MIN( partfin + filesetsize, partitions);
    goto next_elements_set;
  }
  /* part.n.elements saved */


  /* The partitiontable has been changed to include the halo elements. The need for saving the 
     halo nodes may be checked by looking whether the number of how many partitions needs 
     the element has changed. */
  if(anyparthalo) {
    int halonodes;
    neededtimes2 = Ivector(1,noknots);
    halonodes = 0;

    k = 0;
    l = 0;
    for(i=1;i<=noknots;i++) {
      neededtimes2[i] = 0;
      for(j=1;j<=maxneededtimes;j++) 
	if(data->partitiontable[j][i]) {
	  neededtimes2[i] += 1;
	  k++;
	}
      halonodes += neededtimes2[i] - neededtimes[i];
      l += neededtimes[i];
    }
    if(data->maxpartitiontable < maxneededtimes) {
      data->maxpartitiontable = maxneededtimes;
      if(info) printf("With the halos nodes belong to %d partitions in maximum\n",maxneededtimes);
    }
    if(info) printf("There are %d additional shared nodes resulting from the halo.\n",halonodes);
  }
  else {
    neededtimes2 = neededtimes;
  }

  /* Define new BC numbers for indirect connections. These should not be mixed with
     existing BCs as they only serve the purpose of automatically creating the matrix structure. */
  if(indirect) {
    indirecttype = 0;
    for(j=0;j < MAXBOUNDARIES;j++) 
      for(i=1; i <= bound[j].nosides; i++) 
	if(bound[j].types[i] > indirecttype) indirecttype = bound[j].types[i];
    indirecttype++;
    if(info) printf("Indirect connections given index %d and elementtype 102.\n",indirecttype);
  }

  /* The output format is the same for all partitions */
  if(data->dim == 2) 
    sprintf(outstyle,"%%d %%d %%.%dg %%.%dg 0.0\n",decimals,decimals);
  else 
    sprintf(outstyle,"%%d %%d %%.%dg %%.%dg %%.%dg\n",decimals,decimals,decimals);

  if(info) printf("Saving mesh for %d partitions\n",partitions);


  /*********** part.n.nodes *********************/

  partstart = 1;
  partfin = MIN( partitions, filesetsize);

  for(i=1;i<=partitions;i++) {
    needednodes[i] = 0;
    neededtwice[i] = 0;
    sharednodes[i] = 0;
    ownnodes[i] = 0;
  }    

 next_nodes_set:

  for(part=partstart;part<=partfin;part++) {
    sprintf(filename,"%s.%d.%s","part",part,"nodes");
    nofile = part - partstart + 1;
    outfiles[nofile] = fopen(filename,"w");
  }
  
  for(l=1; l <= noknots; l++) {      
    i = l;
    if(reorder) i=invorder[l];

    /*    for(j=1;j<=neededtimes2[i];j++) { */
    for(j=1;j<=maxneededtimes;j++) {

      k = data->partitiontable[j][i];
      if(!k) break;

      if(k < partstart || k > partfin) continue;
      nofile = k - partstart + 1;

      ind = i;
      if(reorder) ind=order[i];

      if(data->dim == 2)
	fprintf(outfiles[nofile],outstyle,ind,-1,data->x[i],data->y[i]);
      else if(data->dim == 3)
	fprintf(outfiles[nofile],outstyle,ind,-1,data->x[i],data->y[i],data->z[i]);	  	    
      
      needednodes[k] += 1;
      if(k == ownerpart[i]) 
	ownnodes[k] += 1;
      else 
	sharednodes[k] += 1;
    }
  }

  for(part=partstart;part<=partfin;part++) {
    nofile = part - partstart + 1;
    fclose(outfiles[nofile]);
  }
  if(partfin < partitions) {
    partstart = partfin + 1;
    partfin = MIN( partfin + filesetsize, partitions);
    goto next_nodes_set;
  }
  /* part.n.nodes saved */
      

  /*********** part.n.shared *********************/

  partstart = 1;
  partfin = MIN( partitions, filesetsize );

 next_shared_set:

  for(part=partstart;part<=partfin;part++) {
    sprintf(filename,"%s.%d.%s","part",part,"shared");
    nofile = part - partstart + 1;
    outfiles[nofile] = fopen(filename,"w");
  }

  for(l=1; l <= noknots; l++) {      
    i = l;
    if(reorder) i = invorder[l];

    if(neededtimes2[i] <= 1) continue;

    for(j=1;j<=neededtimes2[i];j++) {
      k = data->partitiontable[j][i];

      /* if( parthalo[2] && neededtimes[i] == 1 && ownerpart[i] == k ) continue;  */
      
      if(k < partstart || k > partfin) continue;
      nofile = k - partstart + 1;

      ind = i;
      if(reorder) ind = order[i];
      neededtwice[k] += 1; 

      /* if( parthalo[2] ) {
	fprintf(outfiles[nofile],"%d %d %d",ind,MAX(1,neededtimes[i]),ownerpart[i]);      
	for(m=1;m<=neededtimes[i];m++) 
	  if(data->partitiontable[m][i] != ownerpart[i]) fprintf(outfiles[nofile]," %d",data->partitiontable[m][i]);
	fprintf(outfiles[nofile],"\n");
      }
      else { */
      fprintf(outfiles[nofile],"%d %d %d",ind,neededtimes2[i],ownerpart[i]);      
      for(m=1;m<=neededtimes2[i];m++) 
	if(data->partitiontable[m][i] != ownerpart[i]) fprintf(outfiles[nofile]," %d",data->partitiontable[m][i]);
      fprintf(outfiles[nofile],"\n");
      /* } */
    }
  }

  for(part=partstart;part<=partfin;part++) {
    nofile = part - partstart + 1;
    fclose(outfiles[nofile]);
  }
  if(partfin < partitions) {
    partstart = partfin + 1;
    partfin = MIN( partfin + filesetsize, partitions);
    goto next_shared_set;
  }
  /* part.n.shared saved */




   
  /*********** part.n.boundary *********************/
  /* This is still done in partition loop as the subroutines are quite complicated */
  discont = FALSE;
  splitsides = 0;

  maxnosides = 0;
  for(j=0;j < MAXBOUNDARIES;j++) 
    maxnosides = MAX( maxnosides, bound[j].nosides );

  halobcs = 0;
  for(part=1;part<=partitions;part++) { 
    int bcneeded2,closeparent,closeparent2,haloelem,saveelem;

    sprintf(filename,"%s.%d.%s","part",part,"boundary");
    out = fopen(filename,"w");
   
    for(i=minelemtype;i<=maxelemtype;i++)
      sidetypes[i] = 0;
    
    sumsides = 0;

    /* Loop over boundaries and boundary elements therein. */
    for(j=0;j < MAXBOUNDARIES;j++) {

      if(bound[j].nosides == 0 ) continue;

      /* Normal boundary conditions */
      for(i=1; i <= bound[j].nosides; i++) {
	  
	GetBoundaryElement(i,&bound[j],data,sideind,&sideelemtype); 

	bctype = bound[j].types[i];
	nodesd1 = sideelemtype%100;
	  
	parent = bound[j].parent[i];
	parent2 = bound[j].parent2[i];

	
	/* Check whether the side is such that it belongs to the domain.
	   Then it will be always saved - no matter of the halos. */
	trueparent = trueparent2 = FALSE;
	if( parent ) trueparent = (elempart[parent] == part);
	if( parent2 ) trueparent2 = (elempart[parent2] == part);	   
	
	saveelem = FALSE;
	part2 = 0;

	if(trueparent || trueparent2) {
	  /* Either parent must be associated with this partition, otherwise do not save this (except for halo nodes) */
	  if( parent && !trueparent ) {	  
	    splitsides++;
	    /* For halo elements also the 2nd true parent exists by construction. */
	    if(!parthalo[1] && !parthalo[2] && !parthalo[4]) parent = 0;
	  }
	  else if( parent2 && !trueparent2 ) {
	    splitsides++;
	    if(!parthalo[1] && !parthalo[2] && !parthalo[4]) parent2 = 0;
	  }
	  saveelem = TRUE;
	}

	/* If boundary element is not saved by default then study different halos.
	   First look for the extruded halo and offset +/-1 in partition index. */
	if(!saveelem && parthalo[3] && part <= subparts ) {
	  int leftright;
	  for(leftright=-1;leftright <=1;leftright += 2) {	    
	    part2 = part+leftright;
	    if( part2 < 1 || part2 > subparts ) continue;		
	    if( parent ) trueparent = (elempart[parent] == part2);
	    if( parent2 ) trueparent2 = (elempart[parent2] == part2);

	    /* If other parent is not true parent then set the parent to zero because it is not 
	       in this partition. */
	    if(trueparent || trueparent2) {
	      if( parent && !trueparent )
		parent = 0;
	      else if( parent2 && !trueparent2 )
		parent2 = 0;
	      saveelem = TRUE;
	      break;
	    }
	  }
	  if(saveelem) halobcs += 1;
	}
	  
	/* These halo strategies save elements even with just one hit on the original partition 
	   since the aggressive halos ensure that all necessary bulk elements are saved. */ 
	if(!saveelem && (parthalo[2] || parthalo[4]) ) {
	  bcneeded = 0;
	  for(l=0;l<nodesd1;l++) {
	    ind = sideind[l];
	    for(k=1;k<=neededtimes[ind];k++)
	      if(part == data->partitiontable[k][ind]) bcneeded++;	   
	  }
	  saveelem = bcneeded;
	  if(saveelem) halobcs +=1;
	}

	/* If we have nothing to save then just continue to next boundary element. */
	if(!saveelem) continue;
	   
	sumsides++;	
	sidetypes[sideelemtype] += 1;
	    
	if( part2 ) 
	  fprintf(out,"%d/%d %d %d %d %d",
		  sumsides,part2,bctype,parent,parent2,sideelemtype);	    
	else 
	  fprintf(out,"%d %d %d %d %d",
		  sumsides,bctype,parent,parent2,sideelemtype);	  
	    
	if(reorder) {
	  for(l=0;l<nodesd1;l++)
	    fprintf(out," %d",order[sideind[l]]);
	} else {
	  for(l=0;l<nodesd1;l++)
	    fprintf(out," %d",sideind[l]);	  
	}
	fprintf(out,"\n");
      }
    }

    /* The discontinuous stuff is more or less obsolete as these things can now be made also
       within ElmerSolver. */
    /* The second side for discontinuous boundary conditions.
       Note that this has not been treated for orphan control. */
    for(j=0;j < MAXBOUNDARIES;j++) {
      for(i=1; i <= bound[j].nosides; i++) {
	if(bound[j].ediscont) 
	  discont = bound[j].discont[i];
	else 
	  discont = FALSE;
	
	if(!bound[j].parent2[i] || !discont) continue;
	
	GetElementSide(bound[j].parent2[i],bound[j].side2[i],-bound[j].normal[i],
		       data,sideind,&sideelemtype); 
	nodesd1 = sideelemtype%100;	
	
	bcneeded = 0;
	for(l=0;l<nodesd1;l++) {
	  ind = sideind[l];
	  for(k=1;k<=neededtimes[ind];k++)
	    if(part == data->partitiontable[k][ind]) bcneeded++;
	}
	if(bcneeded < nodesd1) continue;
	
	
	trueparent = (elempart[bound[j].parent2[i]] == part);
	if(!trueparent) continue;
	
	sumsides++;
	fprintf(out,"%d %d %d %d ",
		sumsides,bound[j].types[i],bound[j].parent2[i],bound[j].parent[i]);
	
	fprintf(out,"%d ",sideelemtype);
	sidetypes[sideelemtype] += 1;
	if(reorder) {
	  for(l=0;l<nodesd1;l++)
	    fprintf(out,"%d ",order[sideind[l]]);
	} 
	else {
	  for(l=0;l<nodesd1;l++)
	    fprintf(out,"%d ",sideind[l]);	  
	} 
	fprintf(out,"\n");
      }
    }
    sidesinpart[part] = sumsides;
        

    /* Boundary nodes that express indirect couplings between different partitions.
       This makes it possible for ElmerSolver to create a matrix connection that 
       is known to exist. */

    if (indirect) {
      int maxsides,nodesides,maxnodeconnections,connectednodes,m;
      int **nodepairs=NULL,*nodeconnections,**indpairs;      

      nodeconnections = bcnodedummy;
      l = 0;
      maxsides = 0;
      nodesides = 0;

  findindirect:

      /* First calculate the maximum number of additional sides */
      for(i=1;i<=noelements;i++) {

	/* owner partition cannot cause an indirect coupling */
	if(elempart[i] == part) continue;
	
	elemtype = data->elementtypes[i];
	nodesd2 = elemtype%100;
	
	/* Check how many nodes still belong to this partition, 
	   if more than one there may be indirect coupling. */
	for(j=0;j < nodesd2;j++) {
	  elemhit[j] = FALSE;
	  ind = data->topology[i][j];
	  for(k=1;k<=neededtimes[ind];k++) 
	    if(part == data->partitiontable[k][ind]) elemhit[j] = TRUE;
	}
	bcneeded = 0;
	for(j=0;j < nodesd2;j++) 
	  if(elemhit[j]) bcneeded++;
	if(bcneeded <= 1) continue;
	
	if(l == 0) {
	  maxsides += (bcneeded-1)*bcneeded/2;
	} 
	else {
	  for(j=0;j < nodesd2;j++) {	  
	    for(k=j+1;k < nodesd2;k++) {
	      if(elemhit[j] && elemhit[k]) {
		nodesides += 1;

		/* The minimum index always first */
		if(data->topology[i][j] <= data->topology[i][k]) {
		  nodepairs[nodesides][1] = data->topology[i][j];
		  nodepairs[nodesides][2] = data->topology[i][k];
		}
		else {
		  nodepairs[nodesides][1] = data->topology[i][k];
		  nodepairs[nodesides][2] = data->topology[i][j];		  
		}
	      }
	    }
	  }
	}
      }

      /* After first round allocate enough space to memorize all indirect non-element couplings. */      
      if(l == 0) {
	nodepairs = Imatrix(1,maxsides,1,2);
	for(i=1;i<=maxsides;i++)
	  nodepairs[i][1] = nodepairs[i][2] = 0;
	l++;
	goto findindirect;
      }
      if(0) printf("Number of non-element connections is %d\n",nodesides);
      
      
      for(i=1;i<=noknots;i++)
	nodeconnections[i] = 0;
      
      for(i=1;i<=nodesides;i++)
	nodeconnections[nodepairs[i][1]] += 1;
      
      maxnodeconnections = 0;
      for(i=1;i<=noknots;i++)
	maxnodeconnections = MAX(maxnodeconnections, nodeconnections[i]);     
      if(0) printf("Maximum number of node-to-node connections %d\n",maxnodeconnections);

      connectednodes = 0;
      for(i=1;i<=noknots;i++) {
	if(nodeconnections[i] > 0) {
	  connectednodes++;
	  nodeconnections[i] = connectednodes;
	}
      }
      if(0) printf("Number of nodes with non-element connections %d\n",connectednodes);

      indpairs = Imatrix(1,connectednodes,1,maxnodeconnections);
      for(i=1;i<=connectednodes;i++)
	for(j=1;j<=maxnodeconnections;j++)
	  indpairs[i][j] = 0;
      
      for(i=1;i<=nodesides;i++) {
	ind = nodeconnections[nodepairs[i][1]];
	for(j=1;j<=maxnodeconnections;j++) {
	  if(indpairs[ind][j] == 0) {
	    indpairs[ind][j] = i;	    
	    break;
	  }
	}
      }

      /* Remove duplicate connections */
      l = 0;
      for(i=1;i<=connectednodes;i++) {
	for(j=1;j<=maxnodeconnections;j++)
	  for(k=j+1;k<=maxnodeconnections;k++) {
	    ind = indpairs[i][j];
	    ind2 = indpairs[i][k];
	    if(!ind || !ind2) continue;
	    
	    if(!nodepairs[ind][1] || !nodepairs[ind][2]) continue;

	    if(nodepairs[ind][2] == nodepairs[ind2][2]) {
	      nodepairs[ind2][1] = nodepairs[ind2][2] = 0;
	      l++;
	    }
	  }
      }
      if(0) printf("Removed %d duplicate connections\n",l);

      
      /* Remove connections that already exist */
      m = 0;
      for(i=1;i<=noelements;i++) {
	if(elempart[i] != part) continue;
	
	elemtype = data->elementtypes[i];
	nodesd2 = elemtype%100;
	
	for(j=0;j < nodesd2;j++) {
	  ind = nodeconnections[data->topology[i][j]];
	  if(!ind) continue;
	  
	  for(k=0;k < nodesd2;k++) {
	    if(j==k) continue;
	    
	    for(l=1;l<=maxnodeconnections;l++) {
	      ind2 = indpairs[ind][l];
	      if(!ind2) break;

	      if(nodepairs[ind2][1] == data->topology[i][j] && nodepairs[ind2][2] == data->topology[i][k]) {
		nodepairs[ind2][1] = nodepairs[ind2][2] = 0;	    
		m++;
	      }
	    }
	  }
	}
      }
      if(0) printf("Removed %d connections that already exists in other elements\n",m);
      
      for(i=1; i <= nodesides; i++) {
	ind = nodepairs[i][1]; 
	ind2 = nodepairs[i][2];
	if(!ind || !ind2) continue;	
	sumsides++;

	sideelemtype = 102;
	if(reorder) {
	  fprintf(out,"%d %d %d %d %d %d %d\n",
		  sumsides,indirecttype,0,0,sideelemtype,order[ind],order[ind2]);
	} else {
	  fprintf(out,"%d %d %d %d %d %d %d\n",
		  sumsides,indirecttype,0,0,sideelemtype,ind,ind2);	  
	}
	sidetypes[sideelemtype] += 1;
	indirectinpart[part] += 1;	
      }

      /* Finally free some extra space that was allocated */
      free_Imatrix(indpairs,1,connectednodes,1,maxnodeconnections);
      free_Imatrix(nodepairs,1,maxsides,1,2);
    }
    /* End of indirect couplings */


    fclose(out);
    /*********** end of part.n.boundary *********************/



    
    /*********** start of part.n.header *********************/
    tottypes = 0;
    for(i=minelemtype;i<=maxelemtype;i++) {
      if(bulktypes[part][i]) tottypes++;
      if(sidetypes[i]) tottypes++;
    }

    sprintf(filename,"%s.%d.%s","part",part,"header");
    out = fopen(filename,"w");
    fprintf(out,"%-6d %-6d %-6d\n",
	    needednodes[part],elementsinpart[part],sumsides);  

    fprintf(out,"%-6d\n",tottypes);
    for(i=minelemtype;i<=maxelemtype;i++) 
      if(bulktypes[part][i]) 
	fprintf(out,"%-6d %-6d\n",i,bulktypes[part][i]);

    for(i=minelemtype;i<=maxelemtype;i++) 
      if(sidetypes[i]) 
	fprintf(out,"%-6d %-6d\n",i,sidetypes[i]);

    fprintf(out,"%-6d %-6d\n",neededtwice[part],0);
    fclose(out);

    if(info) {
      if( indirect ) {
	if(part == 1) {
	  printf("   %-5s %-10s %-10s %-8s %-8s %-8s\n",
		 "part","elements","nodes","shared","bc elems","indirect");
	}
	printf("   %-5d %-10d %-10d %-8d %-8d %-8d\n",
	       part,elementsinpart[part],ownnodes[part],sharednodes[part],sidesinpart[part],
	       indirectinpart[part]);
      }
      else {
	if( part == 1 ) {
	  printf("   %-5s %-10s %-10s %-8s %-8s\n",
		 "part","elements","nodes","shared","bc elems");
	}
	printf("   %-5d %-10d %-10d %-8d %-8d\n",
	       part,elementsinpart[part],ownnodes[part],sharednodes[part],sidesinpart[part]);	
      }
    }
  }
  /*********** end of part.n.header *********************/

  
  sumelementsinpart = sumownnodes = sumsharednodes = sumsidesinpart = sumindirect = 0;
  for(i=1;i<=partitions;i++) {
    sumelementsinpart += elementsinpart[i];
    sumownnodes += ownnodes[i];
    sumsharednodes += sharednodes[i];
    sumsidesinpart += sidesinpart[i];
    sumindirect += indirectinpart[i];
  }
  n = partitions;
  printf("----------------------------------------------------------------------------------------------\n");
  printf("   ave   %-10.1f %-10.1f %-8.1f %-8.1f %-8.1f\n",
	 1.0*sumelementsinpart/n,1.0*sumownnodes/n,1.0*sumsharednodes/n,
	 1.0*sumsidesinpart/n,1.0*sumindirect/n);

  if( info && halobcs ) {
    printf("Number of boundary elements associated with halo: %d\n",halobcs);
  }
 

#if 0
  {
    int noparents[2][100];

    for(j=0;j<2;j++)
      for(i=0;i<100;i++)
	noparents[j][i] = 0;

    for(j=0;j < MAXBOUNDARIES;j++) {    
      if(bound[j].nosides == 0 ) continue;   
      for(i=1; i <= bound[j].nosides; i++) {      
	GetBoundaryElement(i,&bound[j],data,sideind,&sideelemtype);       	  
	parent = bound[j].parent[i];
	parent2 = bound[j].parent2[i];
	bctype = bound[j].types[i];

	k = 0;
	if(parent) k++;
	if(parent2) k++;

	noparents[k][bctype] += 1;
      }
    }

    printf("Number of BC parents\n");
    for(i=0;i<100;i++) {
      k = 0;
      for(j=0;j<2;j++)
	k = k + noparents[j][i];
      if( k ) {
	printf("BC %d: %d %d %d\n",i,noparents[0][i],noparents[1][i],noparents[2][i]); 
      }
    }
  }
#endif
  
	

  if(splitsides && !anyparthalo) {
    printf("************************* Warning ****************************\n");
    printf("Number or boundary elements split at between parents: %d\n",splitsides);
    printf("This could be a problem for internal jump conditions\n");
    printf("You could try to use '-halobc' flag as remedy with ElmerSolver.\n");
    printf("**************************************************************\n");
  }

  if(anyparthalo) free_Ivector(neededtimes2,1,noknots);
    
  cdstat = chdir("..");
  cdstat = chdir("..");

  if(reorder) free_Ivector(order,1,noknots);
  free_Ivector(needednodes,1,partitions);
  free_Ivector(neededtwice,1,partitions);
  free_Ivector(sharednodes,1,partitions);
  free_Ivector(ownnodes,1,partitions);
  free_Ivector(sidetypes,minelemtype,maxelemtype);
  free_Imatrix(bulktypes,1,partitions,minelemtype,maxelemtype);
  
  if(info) printf("Writing of partitioned mesh finished\n");
  
  return(0);
}


#if USE_METIS 
int ReorderElementsMetis(struct FemType *data,int info)
/* Calls the fill reduction ordering algorithm of Metis library. */
{
  int i,j,k,nn,totcon,maxcon,con,options[8];
  int noelements,noknots,nonodes;
  int *xadj,*adjncy,numflag,*iperm=NULL,*perm=NULL,**newtopology;
  Real *newx,*newy,*newz = NULL;

  noelements = data->noelements;
  noknots = data->noknots;

  if(info) printf("Reordering %d knots and %d elements using Metis reordering routine.\n",
		  noknots,noelements);
  perm = NULL;
  i = CalculateIndexwidth(data,FALSE,perm);
  if(info) printf("Indexwidth of the original node order is %d.\n",i);


  CreateNodalGraph(data,TRUE,info);
  maxcon = data->nodalmaxconnections;

  totcon = 0;
  for(i=1;i<=noknots;i++) {
    for(j=0;j<maxcon;j++) {
      con = data->nodalgraph[j][i];
      if(con) totcon++;
    }
  }
  if(info) printf("There are %d connections altogether\n",totcon);

  xadj = Ivector(0,noknots);
  adjncy = Ivector(0,totcon-1);
  for(i=0;i<totcon;i++) 
    adjncy[i] = 0;

  totcon = 0;
  for(i=1;i<=noknots;i++) {
    xadj[i-1] = totcon;
    for(j=0;j<maxcon;j++) {
      con = data->nodalgraph[j][i];
      if(con) {
	adjncy[totcon] = con-1;
	totcon++;
      }
    }
  }
  xadj[noknots] = totcon;

  nn = noknots;
  numflag = 0;
  for(i=0;i<8;i++) options[i] = 0;
  perm = Ivector(0,noknots-1);
  iperm = Ivector(0,noknots-1);
  
  if(info) printf("Starting Metis reordering routine.\n");

  METIS_NodeND(&nn,xadj,adjncy,&numflag,&options[0],perm,iperm);

  if(info) printf("Finished Metis reordering routine.\n");

  if(info) printf("Moving knots to new positions\n");
  newx = Rvector(1,data->noknots);
  newy = Rvector(1,data->noknots);
  newz = Rvector(1,data->noknots);

  for(i=1;i<=data->noknots;i++) {
    newx[i] = data->x[perm[i-1]+1];
    newy[i] = data->y[perm[i-1]+1];
    newz[i] = data->z[perm[i-1]+1];
  }

  free_Rvector(data->x,1,data->noknots);
  free_Rvector(data->y,1,data->noknots);
  free_Rvector(data->z,1,data->noknots);

  data->x = newx;
  data->y = newy;
  data->z = newz;


  if(info) printf("Changing the element topology\n");

  newtopology = Imatrix(1,noelements,0,data->maxnodes-1);

  for(j=1;j<=noelements;j++) {
    nonodes = data->elementtypes[j]%100;
    for(i=0;i<nonodes;i++) {
      k = data->topology[j][i];
      newtopology[j][i] = iperm[k-1]+1;
    }
  }
  free_Imatrix(data->topology,1,noelements,0,data->maxnodes-1);
  data->topology = newtopology;

  i = CalculateIndexwidth(data,FALSE,perm);
  if(info) printf("Indexwidth of the new node order is %d.\n",i); 

  if(0) printf("Deallocating vectors needed for reordering.\n");
  free_Ivector(iperm,0,noknots-1);
  free_Ivector(perm,0,noknots-1);

  return(0);
}
#endif
