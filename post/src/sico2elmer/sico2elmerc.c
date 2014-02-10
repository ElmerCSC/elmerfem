#include "sico2elmer.h"

#include "../../config.h"

/*
  jv: added fortran function name and char ptr macros to (hopefully) enhance portability
 */

/******************************************************
      Output of SICOPOLIS grid in ELMER POST format
*******************************************************/
void STDCALLBULL FC_FUNC(postgrid,POSTGRID) (float  *xi, /* unscaled x coordinate index i: 0,imax */
	       float  *eta, /* unscaled y coordinate index j from 0,jmax */
	       float  *z_c, /* unscaled z coordinate index i: 0,imax, j: 0,jmax, kc: 0,kcmax */
	       float  *z_t, /* unscaled y coordinate index i: 0,imax, j: 0,jmax, kt: 0,kcmax */
	       float  *deltaX, /* horizontal grod spacing */
	       int    *imax_in, /* grid steps in xi-direction */
	       int    *jmax_in, /* grid steps in eta-direction */
	       int    *kcmax_in, /* grid steps in z-direction in cold ice layer */
	       int    *ktmax_in, /* grid steps in z-direction in temperate ice layer */
	       FC_CHAR_PTR(runname,i1), /*name of run*/
	       FC_CHAR_PTR(ergnum,i2), /*number of file*/
	       int    *maske, /*mask of vertex type */
	       int    *flag){
  register int i, j, k, n, element_nr, boundary_nr,kn;
  int   number_of_layers[2], number_of_elements, elements_in_one_layer, number_of_iced_nodes, number_of_nodes[2], nodes_of_element[8], *nodes_of_side_element, number_of_iced_collums, nodes_in_one_layer, number_of_ice_boundaries, *iced, *boundary, *gridmap, auxiliary;
  int   imax, jmax, kcmax, ktmax;
  float *staggered_grid[2];
  float actual_scaled_coord[3];
  char  groupid[4], filename[80], yes_no, *suffix=".ep";
  FILE  *ptFile;
  
  /* constants */
  imax= *imax_in;
  jmax= *jmax_in;
  kcmax= *kcmax_in;
  ktmax= *ktmax_in;
  nodes_in_one_layer = (imax+2)*(jmax+2);
  elements_in_one_layer = (imax+1)*(jmax+1);
  number_of_layers[0] = kcmax+1;
  number_of_layers[1] = (*flag)*(ktmax);
  number_of_nodes[0] = nodes_in_one_layer*number_of_layers[0];
  number_of_nodes[1] = nodes_in_one_layer*number_of_layers[1];
  /* print out little summary */
  printf("---------------------------------------------------------------\n");
  printf("|        Output of SICOPOLIS Grid vor ELMER POST\n");
  printf("---------------------------------------------------------------\n");
  printf("| imax/jmax/kcmax/ktmax=%d/%d/%d/%d\n",imax, jmax, kcmax, ktmax);
  printf("---------------------------------------------------------------\n");
  printf("| nodes in original grid:\n");
  printf("|       cold layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(kcmax+1), (imax+1)*(jmax+1), (kcmax+1));
  printf("|  temperate layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+1), (imax+1)*(jmax+1), (ktmax+1));
  printf("|                    -------------\n");
  printf("|                    %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+2+kcmax), (imax+1)*(jmax+1), ktmax+2+kcmax);
  printf("---------------------------------------------------------------\n");
  if (*flag)
    printf("| Output of temperate layer enabled\n");
  else  printf("| Output of temperate layer disabled\n");  
  printf("---------------------------------------------------------------\n");
  printf("| nodes in full staggered grid:\n");
  printf("|       cold layer:  %d=%d*%d\n", number_of_nodes[0], nodes_in_one_layer, number_of_layers[0]);
  printf("|  temperate layer:  %d=%d*%d\n", number_of_nodes[1], nodes_in_one_layer, number_of_layers[1]);
  printf("|                    -------------\n");
  printf("|                    %d\n", number_of_nodes[1] +  number_of_nodes[0]);
  printf("---------------------------------------------------------------\n");
  /* allocate memory */
  staggered_grid[0] = (float *) malloc((size_t) (number_of_nodes[0]*3*sizeof(float)));
  if (staggered_grid[0] == NULL){
    printf("ERROR in allocating memory for staggered grid of cold ice layer\n");
    return;
  }
  staggered_grid[1] = (float *) malloc((size_t) (number_of_nodes[1]*3*sizeof(float)));
  if (staggered_grid[1] == NULL){
    printf("ERROR in allocating memory for staggered grid of temperate ice layer\n");
    free(staggered_grid[0]);
    return;
  }
  iced = (int *) malloc((size_t) ((size_t) (imax+1)*(jmax+1)*sizeof(int)));
  if (iced == NULL){
    printf("ERROR in allocating memory for glaciation information\n");    
    free(staggered_grid[0]);free(staggered_grid[1]);
    return;
  }
  /* get staggered grid(s) */
  auxiliary = get_staggered_grid(xi,eta,z_c,imax,jmax,kcmax,deltaX,staggered_grid[0]);
/*   exit(0); */
  if (auxiliary!=number_of_nodes[0]){
    printf("WARNING: number of written %d gridpoints for cold layer does not match calculated %d\n",auxiliary, number_of_nodes[0]);
    printf("ERROR: Interpolation of staggered grid for cold layer failed!\n");
    free(staggered_grid[0]);free(staggered_grid[1]); 
    return;
  }else
    printf("| succeeded in interpolating staggered grid for cold layer\n");
  if (*flag){
    auxiliary = get_staggered_grid(xi,eta,z_t,imax,jmax,(ktmax-1),deltaX,staggered_grid[1]);
    if(auxiliary!=number_of_nodes[1]){
      printf("WARNING: number of written %d gridpoints for cold layer does not match calculated %d\n",auxiliary, number_of_nodes[1]);
      printf("ERROR: Interpolation of staggered grid for temperate layer failed!\n");
      free(staggered_grid[0]);free(staggered_grid[1]); 
      return;
    }else
      printf("| succeeded in interpolating staggered grid for temperate layer\n");
  }else
    printf("| no staggered grid for temperate layer interpolated\n");
  printf("---------------------------------------------------------------\n");
  /* get glaciation info */
  number_of_iced_collums=get_glaciation_info(imax,jmax,iced,maske);
  if (number_of_iced_collums<1){
    printf("| no glaciation\n");
    boundary = (int *) malloc((size_t) 4*sizeof(int)); /* dummy array size */
/*     free(staggered_grid[0]);free(staggered_grid[1]);free(iced);  */
/*     return; */
  }else{
    printf("| number of iced colums:       %d out of %d (%3.2f %%)\n", number_of_iced_collums, (imax+1)*(jmax+1), ((float) number_of_iced_collums)*100.0/((float) (imax+1)*(jmax+1)));
    boundary = (int *) malloc((size_t) number_of_iced_collums*4*sizeof(int));
  }
  if (boundary==NULL){
    printf("ERROR in allocation of memory for ice-boundary information\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced); 
    return;
  }
  gridmap = (int *) malloc((size_t) nodes_in_one_layer*sizeof(int));
  if (gridmap==NULL){
    printf("ERROR in allocation of memory for staggered grid glaciation information\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced); 
    return;
  }
  number_of_iced_nodes = get_grid_map(imax,jmax,iced,gridmap);
  if (number_of_iced_nodes>1){
/*     printf("ERROR in calculation of glaciation information for staggered grid\n"); */
/*     free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(gridmap); */
/*     return;   */
    printf("| number of iced nodes in one layer:  %d out of %d\n",number_of_iced_nodes, nodes_in_one_layer);
  }

  number_of_ice_boundaries = get_glaciation_boundary_info(imax,jmax,iced,boundary);  
  if (number_of_ice_boundaries<1){
    printf("| no glaciation\n");
/*     printf("ERROR in calculation of boundaries\n"); */
/*     free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap); */
/*     return; */
  }else{
    printf("| number of boundary elements lines: %d (%f km)\n", number_of_ice_boundaries, ((float) number_of_ice_boundaries)*(*deltaX));
  }
  /* calculate constants for output mesh*/
  if (*flag){
    number_of_nodes[0] =  number_of_iced_nodes*(kcmax+1);
    number_of_nodes[1] =  number_of_iced_nodes*(ktmax-1)+nodes_in_one_layer;
  }else{
    number_of_nodes[0] =  number_of_iced_nodes*kcmax;
    number_of_nodes[1] =  nodes_in_one_layer;
  }
  number_of_elements = elements_in_one_layer + (1+(*flag))*number_of_iced_collums + (kcmax+(*flag)*ktmax)*number_of_iced_collums + number_of_ice_boundaries*(kcmax+(*flag)*ktmax);
  nodes_of_side_element = (int *) malloc((size_t) number_of_ice_boundaries*4*(kcmax+(*flag)*ktmax)*sizeof(int));
  if (nodes_of_side_element == NULL){
    printf("ERROR in allocation of side element array\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);
    return;
  }
  printf("---------------------------------------------------------------\n");
  printf("| number of nodes:\n");
  printf("|             in cold layer       %d\n",number_of_nodes[0]+(1-(*flag))*nodes_in_one_layer);
  printf("|        in temperate layer       %d\n",number_of_nodes[1]+(*flag)*nodes_in_one_layer);
  printf("|                                 -----------\n");
  printf("|                                 %d\n",number_of_nodes[0]+number_of_nodes[1]);
  printf("---------------------------------------------------------------\n");
  printf("| number of elements:\n");
  printf("|           in cold layer volume  %d\n",kcmax*number_of_iced_collums);
  printf("|      in temperate layer volume  %d\n",(*flag)*ktmax*number_of_iced_collums);
  printf("|                        at base  %d\n",elements_in_one_layer);
  printf("|                         at CTS  %d\n",(*flag)*number_of_iced_collums);
  printf("|                at free surface  %d\n",number_of_iced_collums);
  printf("|      on boundary of cold layer  %d\n",number_of_ice_boundaries*(kcmax));
  printf("| on boundary of temperate layer  %d\n",(*flag)*number_of_ice_boundaries*(ktmax));
  printf("|                                 -----------\n");
  printf("|                                 %d\n",number_of_elements);
  printf("---------------------------------------------------------------\n");
  /* write mesh file header*/
  sprintf(filename,"%5s%2s%s", runname, ergnum, suffix);
  printf("| Writing mesh-file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);
    return;
  }
  fprintf(ptFile, "%d %d 1 1\n", number_of_nodes[0]+number_of_nodes[1],  number_of_elements);
  /* write nodes */
  if (*flag){
    for(j=0,n=0;j<jmax+2;++j){
      for (i=0;i<imax+2;++i){
	fprintf(ptFile, "%6.4e %6.4e %6.4e\n", staggered_grid[1][(j*(imax+2) + i)*3 + 0], staggered_grid[1][(j*(imax+2) + i)*3 + 1], staggered_grid[1][(j*(imax+2) + i)*3 + 2]);++n;
      }
    }    
    for (k=1;k<ktmax;++k){
      for(j=0;j<jmax+2;++j){
	for (i=0;i<imax+2;++i){
	  if (gridmap[j*(imax+2) + i]!=-1){
	    fprintf(ptFile, "%6.4e %6.4e %6.4e\n", staggered_grid[1][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 0], staggered_grid[1][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 1], staggered_grid[1][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]);
	    ++n;
	  }
	}
      }
    }
    for (k=0;k<kcmax+1;++k){
      for(j=0;j<jmax+2;++j){
	for (i=0;i<imax+2;++i){
	  if (gridmap[j*(imax+2) + i]!=-1){
	    fprintf(ptFile, "%6.4e %6.4e %6.4e\n", staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 0], staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 1], staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]);
	    ++n;
	  }
	}
      }
    }
  }else{
    for(j=0,n=0;j<jmax+2;++j){
      for (i=0;i<imax+2;++i){
	fprintf(ptFile, "%6.4e %6.4e %6.4e\n", staggered_grid[0][(j*(imax+2) + i)*3 + 0], staggered_grid[0][(j*(imax+2) + i)*3 + 1], staggered_grid[0][(j*(imax+2) + i)*3 + 2]);++n;
      }      
    }
    for (k=1;k<kcmax+1;++k){
      for(j=0;j<jmax+2;++j){
	for (i=0;i<imax+2;++i){
	  if (gridmap[j*(imax+2) + i]!=-1){
	    fprintf(ptFile, "%6.4e %6.4e %6.4e\n", staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 0], staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 1], staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]);
	    ++n;
	  }
	}
      }
    }
  }
  if (n!=number_of_nodes[0]+number_of_nodes[1]){
    printf("WARNING: Number of written nodes %d does not match calculated %d\n", n, number_of_nodes[0]+number_of_nodes[1]);
    fclose(ptFile);
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);
    return;
  }else printf("| %d nodes written\n",n);
  /* write elements */
  /******** Order of Elements *****
   *     7         6              *
   *     +---------+       E=0    *
   *     |\         \      N=1    *
   *     | \       | \     W=2    *
   *  j  |  \4        \5   S=3    *
   *  ^  |   +-----+---+          *
   *   \ |   |         |          *
   *    \|   |N    |   |          *
   *     + - + - -S+   |          *
   *     3\  |     2\  |          *
   *       \ |        E|          *
   *      W \|        \|          *
   *         +---------+  нннн>i  *
   *         0    S    1          *
   ********************************/
  boundary_nr=0;
  element_nr=0;
  if (*flag){
    sprintf(groupid,"temp");
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	if (iced[j*(imax+1)+i]!=-1){
	  nodes_of_element[0] = j*(imax+2)+i;
	  nodes_of_element[1] = j*(imax+2)+i+1;
	  nodes_of_element[2] = (j+1)*(imax+2)+i+1;
	  nodes_of_element[3] = (j+1)*(imax+2)+i;
	  nodes_of_element[4] = nodes_in_one_layer + gridmap[j*(imax+2)+i];
	  nodes_of_element[5] = nodes_in_one_layer + gridmap[j*(imax+2)+i+1];
	  nodes_of_element[6] = nodes_in_one_layer + gridmap[(j+1)*(imax+2)+i+1];
	  nodes_of_element[7] = nodes_in_one_layer + gridmap[(j+1)*(imax+2)+i];
	  fprintf(ptFile,"%s 808 %d %d %d %d %d %d %d %d\n", groupid,
		  nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3], 
		  nodes_of_element[4], nodes_of_element[5], nodes_of_element[6], nodes_of_element[7]); 
	  ++element_nr;
	  if (boundary[iced[(j*(imax+1)+i)]*4+0]){/* east */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[5];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[6];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+1]){/* north */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[6];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+2]){/* west */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+3]){/* south */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[5];
	    ++boundary_nr;
	  }
	}
      }
    }
    for (k=1;k<ktmax;++k){
      for (j=0;j<jmax+1;++j){
	for (i=0;i<imax+1;++i){
	  if (iced[j*(imax+1)+i]!=-1){
	    for (n=0; n<2; n++){/* lower level: n=0; upper level n=1; each counterclkws */
	      nodes_of_element[n*4] = nodes_in_one_layer +  (k-1+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i];
	      nodes_of_element[n*4+1] = nodes_in_one_layer + (k-1+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i+1];
	      nodes_of_element[n*4+2] = nodes_in_one_layer + (k-1+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i+1];
	      nodes_of_element[n*4+3] = nodes_in_one_layer + (k-1+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i];
	    }
	    fprintf(ptFile,"%s 808 %d %d %d %d %d %d %d %d\n", groupid,
		    nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3], 
		    nodes_of_element[4], nodes_of_element[5], nodes_of_element[6], nodes_of_element[7]); 
	    ++element_nr;
	    if (boundary[iced[(j*(imax+1)+i)]*4+0]){/* east */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[2];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[1];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[5];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[6];
	      ++boundary_nr;
	    }
	    if (boundary[iced[(j*(imax+1)+i)]*4+1]){/* north */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[2];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[6];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	      ++boundary_nr;
	    }
	    if (boundary[iced[(j*(imax+1)+i)]*4+2]){/* west */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	      ++boundary_nr;
	    }
	    if (boundary[iced[(j*(imax+1)+i)]*4+3]){/* south */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[1];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[5];
	      ++boundary_nr;
	    }
	  }
	}
      }
    }
  }
  sprintf(groupid,"cold");
  if (*flag) kn=0;
  else{
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	if (iced[j*(imax+1)+i]!=-1){
	  nodes_of_element[0] = j*(imax+2)+i;
	  nodes_of_element[1] = j*(imax+2)+i+1;
	  nodes_of_element[2] = (j+1)*(imax+2)+i+1;
	  nodes_of_element[3] = (j+1)*(imax+2)+i;
	  nodes_of_element[4] = nodes_in_one_layer + gridmap[j*(imax+2)+i];
	  nodes_of_element[5] = nodes_in_one_layer + gridmap[j*(imax+2)+i+1];
	  nodes_of_element[6] = nodes_in_one_layer + gridmap[(j+1)*(imax+2)+i+1];
	  nodes_of_element[7] = nodes_in_one_layer + gridmap[(j+1)*(imax+2)+i];
	  fprintf(ptFile,"%s 808 %d %d %d %d %d %d %d %d\n", groupid,
		  nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3], 
		  nodes_of_element[4], nodes_of_element[5], nodes_of_element[6], nodes_of_element[7]); 
	  ++element_nr;
	  if (boundary[iced[(j*(imax+1)+i)]*4+0]){/* east */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[5];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[6];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+1]){/* north */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[6];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+2]){/* west */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+3]){/* south */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[5];
	    ++boundary_nr;
	  }
	}
      }
    }
    kn=1;
  }
  for (k=kn;k<kcmax;++k){
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	if (iced[j*(imax+1)+i]!=-1){
	  for (n=0; n<2; n++){/* lower level: n=0; upper level n=1; each counterclkws */
	    nodes_of_element[n*4] = number_of_nodes[1] + (k-kn+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i];
	    nodes_of_element[n*4+1] = number_of_nodes[1] + (k-kn+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i+1];
	    nodes_of_element[n*4+2] = number_of_nodes[1] + (k-kn+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i+1];
	    nodes_of_element[n*4+3] = number_of_nodes[1] + (k-kn+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i];
	  }
	  fprintf(ptFile,"%s 808 %d %d %d %d %d %d %d %d\n", groupid,
		  nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3], 
		  nodes_of_element[4], nodes_of_element[5], nodes_of_element[6], nodes_of_element[7]); 
	  ++element_nr;
	  if (boundary[iced[(j*(imax+1)+i)]*4+0]){/* east */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[5];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[6];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+1]){/* north */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[6];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+2]){/* west */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+3]){/* south */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[5];
	    ++boundary_nr;
	  }
	}
      }
    }
  }
  sprintf(groupid,"base");    
  for (j=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      nodes_of_element[0]=j*(imax+2)+i;
      nodes_of_element[1]=j*(imax+2)+i+1;
      nodes_of_element[2]=(j+1)*(imax+2)+i+1;
      nodes_of_element[3]=(j+1)*(imax+2)+i;
      fprintf(ptFile,"%s 404 %d %d %d %d\n", groupid, nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3]);
      ++element_nr;
    }
  }
  if (*flag){
    sprintf(groupid,"cts"); 
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	if (iced[j*(imax+1)+i]!=-1){
	  nodes_of_element[0] = number_of_nodes[1] + gridmap[j*(imax+2)+i];
	  nodes_of_element[1] = number_of_nodes[1] + gridmap[j*(imax+2)+i+1];
	  nodes_of_element[2] = number_of_nodes[1] + gridmap[(j+1)*(imax+2)+i+1];
	  nodes_of_element[3] = number_of_nodes[1] + gridmap[(j+1)*(imax+2)+i];
	  fprintf(ptFile,"%s 404 %d %d %d %d\n", groupid, nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3]);
	  ++element_nr;
	}
      }      
    }
  }
  sprintf(groupid,"free");
  for (j=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      if (iced[j*(imax+1)+i]!=-1){
	nodes_of_element[0] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[j*(imax+2)+i];
	nodes_of_element[1] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[j*(imax+2)+i+1];
	nodes_of_element[2] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i+1];
	nodes_of_element[3] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i];
	fprintf(ptFile,"%s 404 %d %d %d %d\n", groupid, nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3]);
	++element_nr;
      }
    }
  }  
  if (*flag){  
    sprintf(groupid,"side_t");
    for (n=0;n<number_of_ice_boundaries*ktmax;++n){
      fprintf(ptFile,"%s 404 %d %d %d %d\n", groupid, nodes_of_side_element[n*4+0], nodes_of_side_element[n*4+1], nodes_of_side_element[n*4+2], nodes_of_side_element[n*4+3]);
      ++element_nr;
    }
    sprintf(groupid,"side_c");
    for (n=number_of_ice_boundaries*ktmax;n<number_of_ice_boundaries*(kcmax+ktmax);++n){
      fprintf(ptFile,"%s 404 %d %d %d %d\n", groupid, nodes_of_side_element[n*4+0], nodes_of_side_element[n*4+1], nodes_of_side_element[n*4+2], nodes_of_side_element[n*4+3]);
      ++element_nr;
    }
  }else{
    sprintf(groupid,"side_c");
    for (n=0;n<number_of_ice_boundaries*kcmax;++n){ 
      fprintf(ptFile,"%s 404 %d %d %d %d\n", groupid, nodes_of_side_element[n*4+0], nodes_of_side_element[n*4+1], nodes_of_side_element[n*4+2], nodes_of_side_element[n*4+3]);
      ++element_nr;
    }
  }
  if (element_nr!=number_of_elements)
    printf("WARNING: number of written elements %d does not match calculated %d\n", element_nr, number_of_elements);
  else printf("| %d elements written\n", element_nr);
  printf("---------------------------------------------------------------\n");
  fclose(ptFile);
  free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);
  return;
} 
/*******************************************************
	      write data for ELMER post 
********************************************************/
void STDCALLBULL FC_FUNC(elmerdata,ELMERDATA) (int   *imax_in, /* grid steps in xi-direction */
				    int   *jmax_in, /* grid steps in eta-direction */
				    int   *kcmax_in, /* grid steps in z-direction in cold ice layer */
				    int   *ktmax_in, /* grid steps in z-direction in temperate ice layer */
				    float *z_c, /* z-coordinate in cold layer for given i,j-position in plane  and kc in vertical*/
				    float *z_t, /* z-coordinate in  temperate  layer for given i,j-position in plane  and kc in vertical*/
				    float *vx_c, /* velocity component in x for cold region for given i,j-position in plane and kc in vertical */
				    float *vy_c, /* velocity component in y for cold region for given i,j-position in plane and kc in vertical */
				    float *vz_c, /* velocity component in z for cold region for given i,j-position in plane and kc in vertical */
				    float *age_c, /* age in cold region for given i,j-position in plane and kc in vertical */
				    float *temp_c, /* temperature in cold region for given i,j-position in plane and kt in vertical */
				    float *vx_t, /* velocity component in x for temperated region for given i,j-position in plane and kt in vertical */
				    float *vy_t, /* velocity component in y for temperated region for given i,j-position in plane and kt in vertical */
				    float *vz_t, /* velocity component in z for temperated region for given i,j-position in plane and kt in vertical */
				    float *temp_t_m, /* melting temperature in temperate ice region for given i,j-position in plane and kt in vertical */
				    float *age_t, /* age in temperate region for given i,j-position in plane and kc in vertical */
				    float *omega_t, /* H2O mass fraction in temperated region for given i,j-position in plane and kt in vertical */
				    float *Q_bm, /* production rate of melting-water at bottom for given i,j-position in plane */
				    float *Q_tld, /*  water drainage rate from the temperated region for given i,j-position in plane */
				    float *am_perp, /* ice volume flux through CTS for given i,j-position in plane */
				    float *qx, /* mass-flux in x-direction for given i,j-position in plane */
				    float *qy, /* mass-flux in y-direction for given i,j-position in plane */
				    int   *n_cts, /* polythermal conditions for given i,j-position at base (-1 cold ice; 0 temp. ice base with cold ice above; 1 temp. ice base with temperate ice above; */
				    int   *maske, /* glaciation information for given i,j-position at base (glaciated=0, ice-free=1, 2=sea-point, 3=floating ice) */
				    FC_CHAR_PTR(runname,runname_l), /*name of run*/
				    FC_CHAR_PTR(ergnum,ergnum_l), /*number of file*/
				    int   *flag){
  register int i, j, k, n;
  int   imax, jmax, kcmax, ktmax, kcmin, offset, nodes_in_temp, nodes_in_cold, elements_in_layer, number_of_iced_nodes, number_of_written_nodes, *gridmap, ok;
  int nodes_in_one_layer, number_of_nodes, nodes_in_layer_of_staggered_grid, number_of_iced_collums, *iced, auxiliary, number_of_properties;
  float *array_for_output, *array_for_interpolation, *age, *temperature, *flux[2], *velocity[3], *height, *drainage, *melt, *ice_land_sea, *type_of_base, *float_cts, *float_maske;
  char data_filename[80], *suffix=".dat";
  FILE  *ptFile;
  /* constants */
  imax=*imax_in;
  jmax=*jmax_in;
  kcmax=*kcmax_in;
  ktmax=*ktmax_in;
  elements_in_layer=(imax+1)*(jmax+1);
  nodes_in_layer_of_staggered_grid = (imax+2)*(jmax+2);
  nodes_in_temp = (*flag)*nodes_in_layer_of_staggered_grid*ktmax;
  nodes_in_cold = nodes_in_layer_of_staggered_grid*(kcmax+1);
  number_of_nodes = nodes_in_temp + nodes_in_cold;
  number_of_properties= (NVEC2*2 + NVEC3*3 + NSCAL);
  /* print out little summary */
  printf("---------------------------------------------------------------\n");
  printf("|        Output of SICOPOLIS Data for ELMER POST\n");
  printf("---------------------------------------------------------------\n");
  printf("| imax/jmax/kcmax/ktmax=%d/%d/%d/%d\n",imax, jmax, kcmax, ktmax);
  printf("---------------------------------------------------------------\n");
  printf("| nodes in original grid:\n");
  printf("|       cold layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(kcmax+1), (imax+1)*(jmax+1), (kcmax+1));
  printf("|  temperate layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+1), (imax+1)*(jmax+1), (ktmax+1));
  printf("|                    -------------\n");
  printf("|                    %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+2+kcmax), (imax+1)*(jmax+1), ktmax+2+kcmax);
  printf("---------------------------------------------------------------\n");
  if (*flag)
    printf("| Output of temperate layer enabled\n");
  else  printf("| Output of temperate layer disabled\n");  
  printf("---------------------------------------------------------------\n");
  printf("| nodes in full staggered grid:\n");
  printf("|       cold layer:  %d\n", nodes_in_cold);
  printf("|  temperate layer:  %d\n", nodes_in_temp);
  printf("|                    -------------\n");
  printf("|                    %d\n",number_of_nodes);
  printf("---------------------------------------------------------------\n");
  if ((float_maske = (float *) malloc((size_t) nodes_in_layer_of_staggered_grid*sizeof(float)))==NULL){
    printf("ERROR in allocation of memory!\n");
    return;
  }
  if((float_cts  = (float *) malloc((size_t) nodes_in_layer_of_staggered_grid*sizeof(float)))==NULL){
    printf("ERROR in allocation of memory!\n");
    free(float_maske);
    return;
  }
  if ((array_for_interpolation = (float *) malloc((size_t) number_of_nodes*number_of_properties*sizeof(float)))==NULL){
    printf("ERROR in allocation of memory for interpolating data on staggered grid!\n");
    free(float_cts);
    free(float_maske);
    return;
  }
  if ((iced = (int *) malloc((size_t) nodes_in_layer_of_staggered_grid*sizeof(int)))==NULL){
    printf("ERROR in allocation of memory!\n");
    free(float_cts);    
    free(float_maske);
    free(array_for_interpolation);
  }
  gridmap = (int *) malloc((size_t) nodes_in_layer_of_staggered_grid*sizeof(int));
  if (gridmap==NULL){
    printf("ERROR in allocation of memory!\n");
    free(float_cts);    
    free(float_maske);
    free(array_for_interpolation);
    free(iced);
  }
  number_of_iced_collums=get_glaciation_info(imax,jmax,iced,maske);
  number_of_iced_nodes = get_grid_map(imax,jmax,iced,gridmap);
  if (*flag) number_of_written_nodes=nodes_in_layer_of_staggered_grid+number_of_iced_nodes*(ktmax+kcmax);
  else number_of_written_nodes=nodes_in_layer_of_staggered_grid+number_of_iced_nodes*(kcmax);
  printf("| number of iced colums:       %d out of %d (%3.2f %%)\n", number_of_iced_collums, (imax+1)*(jmax+1), ((float) number_of_iced_collums)*100.0/((float) (imax+1)*(jmax+1)));
  printf("| number of iced nodes in layer of staggered grid: %d out of %d\n", number_of_iced_nodes, nodes_in_layer_of_staggered_grid);
  printf("| number of nodes in written grid: %d of %d\n", number_of_written_nodes, number_of_nodes);
  printf("| number of properties written for each node: %d\n", number_of_properties);
  printf("| size of output array =%d=%d*%d\n", number_of_written_nodes*number_of_properties, number_of_written_nodes,number_of_properties);
  printf("---------------------------------------------------------------\n");
  array_for_output = (float *) malloc((size_t) number_of_written_nodes*number_of_properties*sizeof(float));
  if (array_for_output==NULL){
    printf("ERROR in allocation of memory for output of data!\n");
    free(float_cts);    
    free(float_maske);
    free(array_for_interpolation);
    free(iced);
    free(gridmap);
    return;
  }
  /* assign pointers to array of output */
  height = (float *) &array_for_interpolation[0];
  velocity[0] = (float *) &array_for_interpolation[number_of_nodes];
  velocity[1] = (float *) &array_for_interpolation[2*number_of_nodes];
  velocity[2] = (float *) &array_for_interpolation[3*number_of_nodes];
  drainage = (float *) &array_for_interpolation[4*number_of_nodes];
  melt  = (float *) &array_for_interpolation[5*number_of_nodes];
  ice_land_sea = (float *) &array_for_interpolation[6*number_of_nodes];
  type_of_base = (float *) &array_for_interpolation[7*number_of_nodes];
  temperature = (float *) &array_for_interpolation[8*number_of_nodes];
  age = (float *) &array_for_interpolation[9*number_of_nodes];
  flux[0] = (float *) &array_for_interpolation[10*number_of_nodes];
  flux[1] = (float *) &array_for_interpolation[11*number_of_nodes];
  
  make_float_from_integer_scalar_field(maske, float_maske, nodes_in_layer_of_staggered_grid, 1);
  make_float_from_integer_scalar_field(n_cts, float_cts, nodes_in_layer_of_staggered_grid, 0);

  /* set not used 3d parts of 2d arrays to default value */
  for (n=nodes_in_layer_of_staggered_grid; n<number_of_nodes;++n){
    ice_land_sea[n]=0.0;
    type_of_base[n]=2.0;    
    drainage[n]=0.0;
    melt[n] = 0.0;
    flux[0][n]=0.0;
    flux[1][n]=0.0;
  }
  /* interpolate properties on staggered grid */
  ok=1;
  if (*flag){
    auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, ktmax-1, z_t,  height);
    if (auxiliary!=nodes_in_temp){
      printf("WARNING: number of returned values %d of interpolated height in temperate layer does not match calculated %d\n",auxiliary, nodes_in_temp);
      ok=0;
    }
    auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, ktmax-1, vx_t,  &velocity[0][0]);
    if (auxiliary!=nodes_in_temp){
      printf("WARNING: number of returned values %d of interpolated velocity_x in temperate layer does not match calculated %d\n", auxiliary, nodes_in_temp);
      ok=0;
    }
    auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, ktmax-1, vy_t,  &velocity[1][0]);
    if (auxiliary!=nodes_in_temp){
      printf("WARNING: number of returned values %d of interpolated velocity_y in temperate layer does not match calculated %d\n",auxiliary, nodes_in_temp);
      ok=0;
    }
    auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, ktmax-1, vz_t,  &velocity[2][0]);
    if (auxiliary!=nodes_in_temp){
      printf("WARNING: number of returned values %d of interpolated  velocity_z in temperate layer does not match calculated %d\n",auxiliary, nodes_in_temp);
      ok=0;
    }
    auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, ktmax-1, temp_t_m, temperature);
    if (auxiliary!=nodes_in_temp){
      printf("WARNING: number of returned values %d of interpolated temperature in temperate layer does not match calculated %d\n",auxiliary, nodes_in_temp);
      ok=0;
    }
    auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, ktmax-1, age_t, age);
    if (auxiliary!=nodes_in_temp){
      printf("WARNING: number of returned values %d of interpolated age of ice in temperate layer does not match calculated %d\n",auxiliary, nodes_in_temp);
      ok=0;
    }
  }
  if (!(ok)){
    printf("ERROR in interpolating data on staggered grid in temperate layer!\n");
    free(array_for_interpolation);
    free(float_maske);
    free(float_cts);
    free(iced);
    free(gridmap);
    free(array_for_output);
    return;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, 0, Q_tld, drainage);
  if (auxiliary!=nodes_in_layer_of_staggered_grid){
    printf("WARNING: number of returned values %d of interpolated drainage  at base layer does not match calculated %d\n",auxiliary, nodes_in_layer_of_staggered_grid);
      ok=0;
  }      
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, 0, Q_bm, melt);
  if (auxiliary!=nodes_in_layer_of_staggered_grid){
    printf("WARNING: number of returned values %d of interpolated melting rate  at base layer does not match calculated %d\n",auxiliary, nodes_in_layer_of_staggered_grid);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, 0, float_maske, ice_land_sea);
  if (auxiliary!=nodes_in_layer_of_staggered_grid){
    printf("WARNING: number of returned values %d of interpolated ice-sea-land mask at base does not match calculated %d\n",auxiliary, nodes_in_layer_of_staggered_grid);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, 0, float_cts, type_of_base);
  if (auxiliary!=nodes_in_layer_of_staggered_grid){
    printf("WARNING: number of returned values %d of interpolated ype of base does not match calculated %d\n",auxiliary, nodes_in_layer_of_staggered_grid);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, 0, qx, flux[0]);
  if (auxiliary!=nodes_in_layer_of_staggered_grid){
    printf("WARNING: number of returned values %d of interpolated flux_x at base does not match calculated %d\n",auxiliary, nodes_in_layer_of_staggered_grid);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, 0, qy, flux[1]);
  if (auxiliary!=nodes_in_layer_of_staggered_grid){
    printf("WARNING: number of returned values %d of interpolated flux_y at base does not match calculated %d\n",auxiliary, nodes_in_layer_of_staggered_grid);
      ok=0;
  }
  if (!(ok)){
    printf("ERROR in interpolating data on staggered grid at base!\n");
    free(array_for_interpolation);
    free(float_maske);
    free(float_cts);
    free(iced);
    free(gridmap);
    free(array_for_output);
    return;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, kcmax, z_c,  &height[nodes_in_temp]);
  if (auxiliary!=nodes_in_cold){
    printf("WARNING: number of returned values %d of interpolated height in cold layer does not match calculated %d\n",auxiliary, nodes_in_cold);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, kcmax, vx_c,  &velocity[0][nodes_in_temp]);
  if (auxiliary!=nodes_in_cold){
    printf("WARNING: number of returned values %d of interpolated velocity_x in cold layer does not match calculated %d\n",auxiliary, nodes_in_cold);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, kcmax, vy_c,  &velocity[1][nodes_in_temp]);
  if (auxiliary!=nodes_in_cold){
    printf("WARNING: number of returned values %d of interpolated velocity_y in cold layer does not match calculated %d\n",auxiliary, nodes_in_cold);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, kcmax, vz_c,  &velocity[2][nodes_in_temp]);
  if (auxiliary!=nodes_in_cold){
    printf("WARNING: number of returned values %d of interpolated velocity_z in cold layer does not match calculated %d\n",auxiliary, nodes_in_cold);
      ok=0;
  }  

  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, kcmax, temp_c, &temperature[nodes_in_temp]);
  if (auxiliary!=nodes_in_cold){
    printf("WARNING: number of returned values %d of interpolated temperature in cold layer does not match calculated %d\n",auxiliary, nodes_in_cold);
      ok=0;
  }
  auxiliary=get_interpolated_property_on_staggered_grid(imax, jmax, kcmax, age_c, &age[nodes_in_temp]);
  if (auxiliary!=nodes_in_cold){
    printf("WARNING: number of returned values %d of age of ice interpolated in cold layer does not match calculated %d\n",auxiliary, nodes_in_cold);
      ok=0;
  }
  if (!(ok)){
    printf("ERROR in interpolating data on staggered grid in cold layer!\n");
    free(array_for_interpolation);
    free(float_maske);
    free(float_cts);
    free(iced);
    free(gridmap);
    free(array_for_output);
    return;
  }
  /* write to output array */
  for(j=0,n=0;j<jmax+2;++j){
    for(i=0;i<imax+2;++i){
      array_for_output[(j*(imax+2)+i)*number_of_properties+0]=height[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+1]=velocity[0][j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+2]=velocity[1][j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+3]=velocity[2][j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+4]=drainage[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+5]=melt[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+6]=ice_land_sea[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+7]=type_of_base[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+8]=temperature[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+9]=age[j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+10]=flux[0][j*(imax+2)+i];
      array_for_output[(j*(imax+2)+i)*number_of_properties+11]=flux[1][j*(imax+2)+i];
      ++n;
    }
  }
  if (*flag){
    for (k=1;k<ktmax;++k){
      for(j=0;j<jmax+2;++j){
	for(i=0;i<imax+2;++i){
	  if (gridmap[j*(imax+2) + i]!=-1){
	    array_for_output[n*number_of_properties + 0]=height[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+1]=velocity[0][k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+2]=velocity[1][k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+3]=velocity[2][k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+4]=drainage[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+5]=melt[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+6]=ice_land_sea[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+7]=type_of_base[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+8]=temperature[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+9]=age[k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+10]=flux[0][k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    array_for_output[n*number_of_properties+11]=flux[1][k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	    ++n;
	  }
	}
      }
    }
    kcmin=0;
  }else{    
    kcmin=1;
  }
  for (k=kcmin;k<kcmax+1;++k){
    for(j=0;j<jmax+2;++j){
      for(i=0;i<imax+2;++i){
	if (gridmap[j*(imax+2) + i]!=-1){
	  array_for_output[n*number_of_properties + 0]=height[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+1]=velocity[0][nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+2]=velocity[1][nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+3]=velocity[2][nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+4]=drainage[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+5]=melt[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+6]=ice_land_sea[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+7]=type_of_base[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+8]=temperature[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+9]=age[nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+10]=flux[0][nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  array_for_output[n*number_of_properties+11]=flux[1][nodes_in_temp+k*nodes_in_layer_of_staggered_grid+j*(imax+2)+i];
	  ++n;
	}
      }
    }
  }
  if (n!=number_of_written_nodes){
    printf("ERROR: Number of nodes for data %d does not match calculated %d\n", n, number_of_written_nodes);
  }else{
    /* write binary data on file */
    sprintf(data_filename,"%5s%2s%s",runname, ergnum, suffix);
    printf("| Writing on file %s\n", data_filename);
    if((ptFile=fopen(data_filename, "w"))==NULL){
      printf("ERROR: Could not open Elmer timeslice %s  file for writing!\n", data_filename);  
    }else{  
      fprintf(ptFile, "0 0 %d %d\n", number_of_written_nodes, number_of_properties);
      auxiliary = fwrite((void *) array_for_output, (size_t) number_of_properties*sizeof(float), (size_t) number_of_written_nodes, ptFile);
      if (auxiliary!=number_of_written_nodes){
	printf("ERROR: Only %d of %d items written on file\n",auxiliary,number_of_nodes);
	printf("ERROR: Not able to write binary data on time-slice file %s\n", data_filename);
      }else{
	printf("| Succeeded in writing file\n");
	printf("---------------------------------------------------------------\n");
      }
    }
  }
  fclose(ptFile);
  free(array_for_interpolation);
  free(float_maske);
  free(float_cts);
  free(iced);
  free(gridmap);
  free(array_for_output);
  return;
}
/******************************************************
   Output of SICOPOLIS grid in ELMER Solver format
*******************************************************/
void STDCALLBULL FC_FUNC(pregrid,PREGRID) (float  *xi, /* unscaled x coordinate index i: 0,imax */
			       float  *eta, /* unscaled y coordinate index j from 0,jmax */
			       float  *z_c, /* unscaled z coordinate index i: 0,imax, j: 0,jmax, kc: 0,kcmax */
			       float  *z_t, /* unscaled y coordinate index i: 0,imax, j: 0,jmax, kt: 0,kcmax */
			       int    *imax_in, /* grid steps in xi-direction */
			       int    *jmax_in, /* grid steps in eta-direction */
			       int    *kcmax_in, /* grid steps in z-direction in cold ice layer */
			       int    *ktmax_in, /* grid steps in z-direction in temperate ice layer */
			       FC_CHAR_PTR(runname,runname_l), /*name of run*/
			       FC_CHAR_PTR(ergnum,ergnum_l), /*number of file*/
			       int    *maske, /*mask of vertex type */
			       float  *deltaX, /* stepsize of grid */
			       int    *flag)
{
  register int i, j, k, n, element_nr, boundary_nr,kn;
  int   number_of_layers[2], number_of_bulk_elements, number_of_boundary_elements, elements_in_one_layer, number_of_iced_nodes, number_of_nodes[2], nodes_of_element[8], *nodes_of_side_element, *parent_of_side_element, number_of_iced_collums, nodes_in_one_layer, number_of_ice_boundaries, *iced, *boundary,  *gridmap, auxiliary, idnr, min_max_i_j[2][2];
  int   imax, jmax, kcmax, ktmax;
  float *staggered_grid[2], min_max_xyz[2][3];
  float actual_scaled_coord[3];
  char  groupid[4], filename[80], yes_no, *suffix=".ep";
  FILE  *ptFile;
  
  /* constants */
  imax= *imax_in;
  jmax= *jmax_in;
  kcmax= *kcmax_in;
  ktmax= *ktmax_in;
  nodes_in_one_layer = (imax+2)*(jmax+2);
  elements_in_one_layer = (imax+1)*(jmax+1);
  number_of_layers[0] = kcmax+1;
  number_of_layers[1] = (*flag)*(ktmax);
  number_of_nodes[0] = nodes_in_one_layer*number_of_layers[0];
  number_of_nodes[1] = nodes_in_one_layer*number_of_layers[1];
  /* print out little summary */
  printf("---------------------------------------------------------------\n");
  printf("|        Output of SICOPOLIS Grid for ELMER Solver\n");
  printf("---------------------------------------------------------------\n");
  printf("| imax/jmax/kcmax/ktmax=%d/%d/%d/%d\n",imax, jmax, kcmax, ktmax);
  printf("---------------------------------------------------------------\n");
  printf("| nodes in original grid:\n");
  printf("|       cold layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(kcmax+1), (imax+1)*(jmax+1), (kcmax+1));
  printf("|  temperate layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+1), (imax+1)*(jmax+1), (ktmax+1));
  printf("|                    -------------\n");
  printf("|                    %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+2+kcmax), (imax+1)*(jmax+1), ktmax+2+kcmax);
  printf("| x = %.1f -> %.1f , y = %.1f -> %.1f, dx=%.1f\n", xi[0], xi[imax], eta[0], eta[jmax], *deltaX);
  printf("---------------------------------------------------------------\n");
  if (*flag)
    printf("| Output of temperate layer enabled\n");
  else  printf("| Output of temperate layer disabled\n");  
  printf("---------------------------------------------------------------\n");
  printf("| nodes in full staggered grid:\n");
  printf("|       cold layer:  %d=%d*%d\n", number_of_nodes[0], nodes_in_one_layer, number_of_layers[0]);
  printf("|  temperate layer:  %d=%d*%d\n", number_of_nodes[1], nodes_in_one_layer, number_of_layers[1]);
  printf("|                    -------------\n");
  printf("|                    %d\n", number_of_nodes[1] +  number_of_nodes[0]);
  printf("---------------------------------------------------------------\n");
  /* allocate memory */
  staggered_grid[0] = (float *) malloc((size_t) (number_of_nodes[0]*3*sizeof(float)));
  if (staggered_grid[0] == NULL){
    printf("ERROR in allocating memory for staggered grid of cold ice layer\n");
    return;
  }
  staggered_grid[1] = (float *) malloc((size_t) (number_of_nodes[1]*3*sizeof(float)));
  if (staggered_grid[1] == NULL){
    printf("ERROR in allocating memory for staggered grid of temperate ice layer\n");
    free(staggered_grid[0]);
    return;
  }
  iced = (int *) malloc((size_t) ((size_t) (imax+1)*(jmax+1)*sizeof(int)));
  if (iced == NULL){
    printf("ERROR in allocating memory for glaciation information\n");    
    free(staggered_grid[0]);free(staggered_grid[1]);
    return;
  }
  /* get staggered grid(s) */
  auxiliary = get_staggered_grid(xi,eta,z_c,imax,jmax,kcmax,deltaX,staggered_grid[0]);
  if (auxiliary!=number_of_nodes[0]){
    printf("WARNING: number of written %d gridpoints for cold layer does not match calculated %d\n",auxiliary, number_of_nodes[0]);
    printf("ERROR: Interpolation of staggered grid for cold layer failed!\n");
    free(staggered_grid[0]);free(staggered_grid[1]); 
    return;
  }else
    printf("| succeeded in interpolating staggered grid for cold layer\n");
  if (*flag){
    auxiliary = get_staggered_grid(xi,eta,z_t,imax,jmax,(ktmax-1),deltaX,staggered_grid[1]);
    if(auxiliary!=number_of_nodes[1]){
      printf("WARNING: number of written %d gridpoints for cold layer does not match calculated %d\n",auxiliary, number_of_nodes[1]);
      printf("ERROR: Interpolation of staggered grid for temperate layer failed!\n");
      free(staggered_grid[0]);free(staggered_grid[1]); 
      return;
    }else
      printf("| succeeded in interpolating staggered grid for temperate layer\n");
  }else
    printf("| no staggered grid for temperate layer interpolated\n");
  printf("---------------------------------------------------------------\n");
  printf("| staggered grid:\n");
  printf("| x = %.1f -> %.1f , y = %.1f -> %.1f\n",staggered_grid[0][0], staggered_grid[0][(imax+1)*3], staggered_grid[0][1], staggered_grid[0][(jmax+1)*(imax+2)*3+1]);
  printf("---------------------------------------------------------------\n");
  /* get glaciation info */
  number_of_iced_collums=get_glaciation_info(imax,jmax,iced,maske);
  if (number_of_iced_collums<1){
  printf("| no glaciation\n");
/*     printf("ERROR in calculation of glaciation\n"); */
/*     free(staggered_grid[0]);free(staggered_grid[1]);free(iced);  */
/*     return; */
  }else{
    printf("| number of iced colums:       %d out of %d (%3.2f %%)\n", number_of_iced_collums, (imax+1)*(jmax+1), ((float) number_of_iced_collums)*100.0/((float) (imax+1)*(jmax+1)));
  }
  gridmap = (int *) malloc((size_t) nodes_in_one_layer*sizeof(int));
  if (gridmap==NULL){
    printf("ERROR in allocation of memory for staggered grid glaciation information\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced); 
    return;
  }
  number_of_iced_nodes = get_grid_map(imax,jmax,iced,gridmap);
/*   if (number_of_iced_nodes<1){ */
/*     printf("ERROR in calculation of glaciation information for staggered grid\n"); */
/*     free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(gridmap); */
/*     return; */
/*   } */
  printf("| number of iced nodes in one layer:  %d out of %d\n",number_of_iced_nodes, nodes_in_one_layer);
  boundary = (int *) malloc((size_t) number_of_iced_collums*4*sizeof(int));
  if (boundary==NULL){
    printf("ERROR in allocation of memory for ice-boundary information\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(gridmap); 
    return;
  }
  number_of_ice_boundaries = get_glaciation_boundary_info(imax,jmax,iced,boundary);  
/*   if (number_of_ice_boundaries<1){ */
/*     printf("ERROR in calculation of boundaries\n"); */
/*     free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap); */
/*     return; */
/*   } */
  printf("| number of boundary elements lines: %d (%f km)\n", number_of_ice_boundaries, ((float) number_of_ice_boundaries)*(*deltaX));
   /* calculate constants for output mesh*/
  if (*flag){
    number_of_nodes[0] =  number_of_iced_nodes*(kcmax+1);
    number_of_nodes[1] =  number_of_iced_nodes*(ktmax);
  }else{
    number_of_nodes[0] =  number_of_iced_nodes*(kcmax+1);
    number_of_nodes[1] =  0;
  }
  number_of_bulk_elements =  (kcmax+(*flag)*ktmax)*number_of_iced_collums;
  number_of_boundary_elements = number_of_ice_boundaries*(kcmax+(*flag)*ktmax) + (2+(*flag))*number_of_iced_collums;
  nodes_of_side_element = (int *) malloc((size_t) number_of_ice_boundaries*4*(kcmax+(*flag)*ktmax)*sizeof(int));
  if (nodes_of_side_element == NULL){
    printf("ERROR in allocation of side element array\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);
    return;
  }
  parent_of_side_element = (int *) malloc((size_t) number_of_ice_boundaries*(kcmax+(*flag)*ktmax)*sizeof(int));
  if (parent_of_side_element == NULL){
    printf("ERROR in allocation of side element array\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);
    return;
  }
  printf("---------------------------------------------------------------\n");
  printf("| number of nodes:\n");
  printf("|             in cold layer       %d\n",number_of_nodes[0]);
  printf("|        in temperate layer       %d\n",number_of_nodes[1]);
  printf("|                                 -----------\n");
  printf("|                                 %d\n",number_of_nodes[0]+number_of_nodes[1]);
  printf("---------------------------------------------------------------\n");
  printf("| number of elements:\n");
  printf("|           in cold layer volume  %d\n",kcmax*number_of_iced_collums);
  printf("|      in temperate layer volume  %d\n",(*flag)*ktmax*number_of_iced_collums);
  printf("|                                 -----------\n");
  printf("|                                 %d\n",number_of_bulk_elements);
  printf("|                                 -----------\n");
  printf("|                                 \n");
  printf("|                        at base  %d\n",number_of_iced_collums);
  printf("|                         at cts  %d\n",number_of_iced_collums);
  printf("|                at free surface  %d\n",number_of_iced_collums);
  printf("|      on boundary of cold layer  %d\n",number_of_ice_boundaries*(kcmax));
  printf("| on boundary of temperate layer  %d\n",(*flag)*number_of_ice_boundaries*(ktmax));
  printf("|                                 -----------\n");
  printf("|                                 %d\n",number_of_boundary_elements);
  printf("---------------------------------------------------------------\n");
  /* writing header file */
  sprintf(filename,"mesh.header");
  printf("| Writing mesh header file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
    return;
  }
  fprintf(ptFile,"%d %d %d\n",number_of_nodes[0]+number_of_nodes[1], number_of_bulk_elements,number_of_boundary_elements);
  fprintf(ptFile,"2\n");
  fprintf(ptFile,"808 %d\n", number_of_bulk_elements);
  fprintf(ptFile,"404 %d\n", number_of_boundary_elements);
  printf("| succeeded in writting mesh header file %s.\n",filename);
  printf("---------------------------------------------------------------\n");
  fclose(ptFile);
  /* check min/max values of grid */
  min_max_i_j[0][0]=imax;
  min_max_i_j[1][0]=0;
  min_max_i_j[0][1]=jmax;
  min_max_i_j[1][1]=0;
  for(j=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      if (iced[j*(imax+1)+i]!=-1){
	if (i < min_max_i_j[0][0]) min_max_i_j[0][0]=i;
	if (i > min_max_i_j[1][0]) min_max_i_j[1][0]=i;
	if (j < min_max_i_j[0][1]) min_max_i_j[0][1]=j;
	if (j > min_max_i_j[1][1]) min_max_i_j[1][1]=j;
      }
    }
  }
  printf("| glaciation information:\n");
  printf("| original grid:\n");
  printf("| i = %d -> %d of %d \n", min_max_i_j[0][0], min_max_i_j[1][0], imax);
  printf("| j = %d -> %d of %d \n", min_max_i_j[0][1], min_max_i_j[1][1], jmax);
  min_max_i_j[0][0]=imax+1;
  min_max_i_j[1][0]=0;
  min_max_i_j[0][1]=jmax+1;
  min_max_i_j[1][1]=0;
  for (i=0;i<(imax+1)*(jmax+1);++i){
    if (gridmap[i]!=-1){
      for (n=0;n<3;++n)
	min_max_xyz[0][n]=staggered_grid[0][i*3 + n];
      break;
    }
  }
  for(j=0;j<jmax+2;++j){
    for (i=0;i<imax+2;++i){
      if (gridmap[j*(imax+2) + i]!=-1){	
	if (i < min_max_i_j[0][0]) min_max_i_j[0][0]=i;
	if (i > min_max_i_j[1][0]) min_max_i_j[1][0]=i;
	if (j < min_max_i_j[0][1]) min_max_i_j[0][1]=j;
	if (j > min_max_i_j[1][1]) min_max_i_j[1][1]=j;
	for (n=0;n<2;++n){
	  if (min_max_xyz[0][n] > staggered_grid[0][(j*(imax+2) + i)*3 + n]) min_max_xyz[0][n] = staggered_grid[0][(j*(imax+2) + i)*3 + n];
	  if (min_max_xyz[1][n] < staggered_grid[0][(j*(imax+2) + i)*3 + n]) min_max_xyz[1][n] = staggered_grid[0][(j*(imax+2) + i)*3 + n];
	}
	if (*flag){
	  if (min_max_xyz[0][2] > staggered_grid[0][(j*(imax+2) + i)*3 + 2]) min_max_xyz[0][2] = staggered_grid[1][(j*(imax+2) + i)*3 + 2];
	}
	for (k=0;k<(kcmax+1);++k){
	  if (min_max_xyz[0][2] > staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]) min_max_xyz[0][2] = staggered_grid[0][(j*(imax+2) + i)*3 + 2];
	  if (min_max_xyz[1][2] < staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]) min_max_xyz[1][2] = staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2];
	}
      }
    }
  }
  printf("| staggered grid:\n");
  printf("| i = %d -> %d of %d \n", min_max_i_j[0][0], min_max_i_j[1][0], imax);
  printf("| j = %d -> %d of %d \n", min_max_i_j[0][1], min_max_i_j[1][1], jmax);
  printf("| x = %.8f ->  %.8f\n",min_max_xyz[0][0], min_max_xyz[1][0]);
  printf("| y = %.8f ->  %.8f\n",min_max_xyz[0][1], min_max_xyz[1][1]);
  printf("| z = %.8f ->  %.8f\n",min_max_xyz[0][2], min_max_xyz[1][2]);
  printf("---------------------------------------------------------------\n");
  
  /* writing nodes file */
  sprintf(filename,"mesh.nodes");
  printf("| Writing node-file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
    return;
  }
  n=0;
  if (*flag){
    for (k=0;k<ktmax;++k){
      for(j=0;j<jmax+2;++j){
	for (i=0;i<imax+2;++i){
	  if (gridmap[j*(imax+2) + i]!=-1){
	    fprintf(ptFile, "%d -1 %.8f %.8f %.8f\n", n+1, staggered_grid[1][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 0], staggered_grid[1][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 1], staggered_grid[1][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]);
	    ++n;
	  }
	}
      }
    }
  }
  for (k=0;k<kcmax+1;++k){
    for(j=0;j<jmax+2;++j){
      for (i=0;i<imax+2;++i){
	if (gridmap[j*(imax+2) + i]!=-1){
	  fprintf(ptFile, "%d -1 %.8f %.8f %.8f\n", n+1, staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 0], staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 1], staggered_grid[0][(k*nodes_in_one_layer + j*(imax+2) + i)*3 + 2]);
	  ++n;
	}
      }
    }
  }
  if (n!=number_of_nodes[0]+number_of_nodes[1]){
    printf("ERROR: number of written nodes %d does not match calculated %d\n",n, number_of_nodes[0]+number_of_nodes[1]);
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
    fclose(ptFile);
    return;
  }
  printf("| succeeded in writting node file %s.\n",filename);
  printf("---------------------------------------------------------------\n");
  /* writing element file */
  sprintf(filename,"mesh.elements");
  printf("| Writing bulk element file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
    return;
  }
  n=0;element_nr=0;boundary_nr=0; idnr=1;
  if (*flag){
    for (k=0;k<ktmax;++k){
      for (j=0;j<jmax+1;++j){
	for (i=0;i<imax+1;++i){
	  if (iced[j*(imax+1)+i]!=-1){
	    for (n=0; n<2; n++){/* lower level: n=0; upper level n=1; each counterclkws */
	      nodes_of_element[n*4] = (k+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i] + 1;
	      nodes_of_element[n*4+1] = (k+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i+1] + 1;
	      nodes_of_element[n*4+2] = (k+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i+1] + 1;
	      nodes_of_element[n*4+3] = (k+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i] + 1;
	    }
	    fprintf(ptFile,"%d %d 808 %d %d %d %d %d %d %d %d\n", element_nr+1, idnr,
		    nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3], 
		    nodes_of_element[4], nodes_of_element[5], nodes_of_element[6], nodes_of_element[7]);
	    if (element_nr != k+number_of_iced_collums+ iced[j*(imax+1)+i]) printf("element_nr=%d does not match %d\n", element_nr, k+number_of_iced_collums+ iced[j*(imax+1)+i]);
	    ++element_nr;
	    if (boundary[iced[(j*(imax+1)+i)]*4+0]){/* east */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[2];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[1];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[5];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[6];
	      parent_of_side_element[boundary_nr]=element_nr;
	      ++boundary_nr;
	    }
	    if (boundary[iced[(j*(imax+1)+i)]*4+1]){/* north */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[2];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[6];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	      parent_of_side_element[boundary_nr]=element_nr;
	      ++boundary_nr;
	    }
	    if (boundary[iced[(j*(imax+1)+i)]*4+2]){/* west */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	      parent_of_side_element[boundary_nr]=element_nr;
	      ++boundary_nr;
	    }
	    if (boundary[iced[(j*(imax+1)+i)]*4+3]){/* south */
	      nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[1];
	      nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	      nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	      nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[5];
	      parent_of_side_element[boundary_nr]=element_nr;
	      ++boundary_nr;
	    }
	  }
	}
      }
    }
    idnr=2;
  }
  for (k=0;k<kcmax;++k){
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	if (iced[j*(imax+1)+i]!=-1){
	  for (n=0; n<2; n++){/* lower level: n=0; upper level n=1; each counterclkws */
	    nodes_of_element[n*4] = number_of_nodes[1] + (k+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i] + 1;
	    nodes_of_element[n*4+1] = number_of_nodes[1] + (k+n)*number_of_iced_nodes + gridmap[j*(imax+2)+i+1] + 1;
	    nodes_of_element[n*4+2] = number_of_nodes[1] + (k+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i+1] + 1;
	    nodes_of_element[n*4+3] = number_of_nodes[1] + (k+n)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i] + 1;
	  }
	  fprintf(ptFile,"%d %d 808 %d %d %d %d %d %d %d %d\n", element_nr+1, idnr,
		  nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3], 
		  nodes_of_element[4], nodes_of_element[5], nodes_of_element[6], nodes_of_element[7]); 
	  if (element_nr != k*number_of_iced_collums+ iced[j*(imax+1)+i]) printf("element_nr=%d does not match %d\n", element_nr, k*number_of_iced_collums+ iced[j*(imax+1)+i]);
	  ++element_nr;
	  if (boundary[iced[(j*(imax+1)+i)]*4+0]){/* east */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[5];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[6];
	    parent_of_side_element[boundary_nr]=element_nr;
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+1]){/* north */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[2];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[6];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    parent_of_side_element[boundary_nr]=element_nr;
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+2]){/* west */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[3];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[7];
	    parent_of_side_element[boundary_nr]=element_nr;
	    ++boundary_nr;
	  }
	  if (boundary[iced[(j*(imax+1)+i)]*4+3]){/* south */
	    nodes_of_side_element[boundary_nr*4+0]=nodes_of_element[1];
	    nodes_of_side_element[boundary_nr*4+1]=nodes_of_element[0];
	    nodes_of_side_element[boundary_nr*4+2]=nodes_of_element[4];
	    nodes_of_side_element[boundary_nr*4+3]=nodes_of_element[5];
	    parent_of_side_element[boundary_nr]=element_nr;
	    ++boundary_nr;
	  }
	}
      }
    }
  }
  if (number_of_bulk_elements!=element_nr){
    printf("ERROR: Number of written bulk elements %d does not match calculated %d\n", number_of_bulk_elements, element_nr);
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
    return;
  }
  printf("| succeeded in writting bulk element file %s.\n",filename);
  printf("---------------------------------------------------------------\n");  
  fclose(ptFile);
  /* writing boundary element file */
  sprintf(filename,"mesh.boundary");
  printf("| Writing boundary element file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer boundary element file for writing!\n");
    free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
    return;
  }
  boundary_nr=0;
  /* base */
  idnr=1;
  for (j=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      if (iced[j*(imax+1)+i]!=-1){
	nodes_of_element[0]= gridmap[j*(imax+2)+i]+1;
	nodes_of_element[1]= gridmap[j*(imax+2)+i+1]+1;
	nodes_of_element[2]= gridmap[(j+1)*(imax+2)+i+1]+1;
	nodes_of_element[3]= gridmap[(j+1)*(imax+2)+i]+1;
	fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n", boundary_nr+1, idnr, iced[j*(imax+1)+i]+1, nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3]);
	++boundary_nr;
      }
    }
  }
  /* cts */
  if (*flag){
    ++idnr;
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	if (iced[j*(imax+1)+i]!=-1){
	  nodes_of_element[0] = number_of_nodes[1] + gridmap[j*(imax+2)+i]+1;
	  nodes_of_element[1] = number_of_nodes[1] + gridmap[j*(imax+2)+i+1]+1;
	  nodes_of_element[2] = number_of_nodes[1] + gridmap[(j+1)*(imax+2)+i+1]+1;
	  nodes_of_element[3] = number_of_nodes[1] + gridmap[(j+1)*(imax+2)+i]+1;
	  fprintf(ptFile,"%d %d %d %d 404 %d %d %d %d\n", boundary_nr+1, idnr, (ktmax-2)*number_of_iced_collums+iced[j*(imax+1)+i]+1, (ktmax)*number_of_iced_collums+iced[j*(imax+1)+i]+1, nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3]);
	  ++boundary_nr;
	}
      }      
    }
  }
  /* free surface */
  ++idnr;
  for (j=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      if (iced[j*(imax+1)+i]!=-1){
	nodes_of_element[0] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[j*(imax+2)+i]+1;
	nodes_of_element[1] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[j*(imax+2)+i+1]+1;
	nodes_of_element[2] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i+1]+1;
	nodes_of_element[3] = number_of_nodes[1] + (kcmax)*number_of_iced_nodes + gridmap[(j+1)*(imax+2)+i]+1;
	fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n", boundary_nr+1, idnr, ((ktmax)*(*flag)+kcmax-1)*number_of_iced_collums+iced[j*(imax+1)+i]+1, nodes_of_element[0], nodes_of_element[1], nodes_of_element[2], nodes_of_element[3]);
	++boundary_nr;
      }
    }
  }
  /* sides */
  if (*flag){  
    ++idnr;
    for (n=0;n<number_of_ice_boundaries*ktmax;++n){
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n", boundary_nr+1, parent_of_side_element[n], idnr, nodes_of_side_element[n*4+0], nodes_of_side_element[n*4+1], nodes_of_side_element[n*4+2], nodes_of_side_element[n*4+3]);
      ++boundary_nr;
    }
    ++idnr;
    for (;n<number_of_ice_boundaries*(kcmax+ktmax);++n){
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n", boundary_nr+1, idnr, parent_of_side_element[n], nodes_of_side_element[n*4+0], nodes_of_side_element[n*4+1], nodes_of_side_element[n*4+2], nodes_of_side_element[n*4+3]);
      ++boundary_nr;
    }
  }else{
    ++idnr;
    for (n=0;n<number_of_ice_boundaries*kcmax;++n){ 
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n", boundary_nr+1, idnr, parent_of_side_element[n], nodes_of_side_element[n*4+0], nodes_of_side_element[n*4+1], nodes_of_side_element[n*4+2], nodes_of_side_element[n*4+3]);
      ++boundary_nr;
    }
  }
  if (boundary_nr!=number_of_boundary_elements){
    printf("ERROR: Number of written boundary elements %d does not match calculated %d.\n", boundary_nr, number_of_boundary_elements);
  }else{
    printf("| succeeded in writting boundary element file %s.\n",filename);
    printf("---------------------------------------------------------------\n");
  }
  fclose(ptFile);
  free(staggered_grid[0]);free(staggered_grid[1]);free(iced);free(boundary);free(gridmap);free(nodes_of_side_element);free(parent_of_side_element);
  return;  
}

void STDCALLBULL FC_FUNC(asciidata,ASCIIDATA) (float  *xi, /* unscaled x coordinate index i: 0,imax */
		float  *eta, /* unscaled y coordinate index j from 0,jmax */
		int   *imax_in, /* grid steps in xi-direction */
		int   *jmax_in, /* grid steps in eta-direction */
		int   *kcmax_in, /* grid steps in z-direction in cold ice layer */
		int   *ktmax_in, /* grid steps in z-direction in temperate ice layer */
		float *z_c, /* z-coordinate in cold layer for given i,j-position in plane  and kc in vertical*/
		float *z_t, /* z-coordinate in  temperate  layer for given i,j-position in plane  and kc in vertical*/
		float *vx_c, /* velocity component in x for cold region for given i,j-position in plane and kc in vertical */
		float *vy_c, /* velocity component in y for cold region for given i,j-position in plane and kc in vertical */
		float *vz_c, /* velocity component in z for cold region for given i,j-position in plane and kc in vertical */
		float *age_c, /* age in cold region for given i,j-position in plane and kc in vertical */
		float *temp_c, /* temperature in cold region for given i,j-position in plane and kt in vertical */
		float *vx_t, /* velocity component in x for temperated region for given i,j-position in plane and kt in vertical */
		float *vy_t, /* velocity component in y for temperated region for given i,j-position in plane and kt in vertical */
		float *vz_t, /* velocity component in z for temperated region for given i,j-position in plane and kt in vertical */
		float *temp_t_m, /* melting temperature in temperate ice region for given i,j-position in plane and kt in vertical */
		float *age_t, /* age in temperate region for given i,j-position in plane and kc in vertical */
		float *omega_t, /* H2O mass fraction in temperated region for given i,j-position in plane and kt in vertical */
		float *Q_bm, /* production rate of melting-water at bottom for given i,j-position in plane */
		float *Q_tld, /*  water drainage rate from the temperated region for given i,j-position in plane */
		float *am_perp, /* ice volume flux through CTS for given i,j-position in plane */
		float *qx, /* mass-flux in x-direction for given i,j-position in plane */
		float *qy, /* mass-flux in y-direction for given i,j-position in plane */
		int   *n_cts, /* polythermal conditions for given i,j-position at base (-1 cold ice; 0 temp. ice base with cold ice above; 1 temp. ice base with temperate ice above; */
		int   *maske, /* glaciation information for given i,j-position at base (glaciated=0, ice-free=1, 2=sea-point, 3=floating ice) */
		FC_CHAR_PTR(runname,p1), /*name of run*/
		FC_CHAR_PTR(ergnum,p2),  /*number of file*/
		int   *flag)
{ 
  register int i, j, k, n, current_index;
  int   imax, jmax, kcmax, ktmax, kcmin, offset, nodes_in_temp, nodes_in_cold, elements_in_layer, number_of_iced_nodes, number_of_written_nodes, *gridmap, ok;
  int nodes_in_one_layer, number_of_nodes, nodes_in_layer_of_staggered_grid, number_of_iced_collums, *iced, auxiliary, number_of_properties, outputflags[12], failed=0;
  char filename[80], inputstring[40], dummy_string[39], user_name[10],*suffix=".asc";
  struct stat buf;
  time_t  how_late_is_it;
  FILE  *ptFile;
  /* constants */
  imax=*imax_in;
  jmax=*jmax_in;
  kcmax=*kcmax_in;
  ktmax=*ktmax_in;
  elements_in_layer=(imax+1)*(jmax+1);
  printf("---------------------------------------------------------------\n");
  printf("|        Output of SICOPOLIS Data in ASCII-format - Hi Olli :-)\n");
  printf("---------------------------------------------------------------\n");
  printf("| imax/jmax/kcmax/ktmax=%d/%d/%d/%d\n",imax, jmax, kcmax, ktmax);
  printf("---------------------------------------------------------------\n");
  printf("| nodes in grid:\n");
  printf("|       cold layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(kcmax+1), (imax+1)*(jmax+1), (kcmax+1));
  printf("|  temperate layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+1), (imax+1)*(jmax+1), (ktmax+1));
  printf("|                    -------------\n");
  printf("|                    %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+2+kcmax), (imax+1)*(jmax+1), ktmax+2+kcmax);
  printf("---------------------------------------------------------------\n");
  printf("| Trying to load preferences from .sico2elmer.rc\n");
  /* read preference file */
  if (stat(".sico2elmer.rc" , &buf)!=0){
    if (errno ==  ENOENT){
      printf("WARNING: The preference-file .sico2elmer.rc does not exist. Writing all output-components by default\n");
      failed=1;
    }
    else if (errno ==  EACCES)
      printf("WARNING: No permission to read file .sico2elmer.rc. Writing all output-components by default\n");
    else printf("WARNING: Error occured during reading  preference-file .sico2elmer.rc. Writing all output-components by default\n");
    failed=1;
  }else{
    if((ptFile=fopen(".sico2elmer.rc", "r"))==NULL){
      printf("WARNING: Error while opening file .sico2elmer.rc . Writing all output-components\n");
      failed=1;
    }else{
      /* read in comment lines */
      for (n=0;n<7;++n){
	if (fgets(inputstring, 40, ptFile)==NULL){ 
	  printf("WARNING: Error while reading .sico2elmer.rc . Writing all output-components\n");
	  failed=1;
	}
      }
      /* read: flag for output of temperate layer */
      if (fgets(inputstring, 40, ptFile)==NULL){ 
	printf("WARNING: Error while reading .sico2elmer.rc . Writing all output-components\n");
	failed=1;
      }else{
	if (sscanf(inputstring, "%1d%s",&outputflags[0], dummy_string)<0){
	  printf("WARNING: Error while copying inputstring no 0 from .sico2elmer.rc\n");
	  outputflags[0]=1;
	}
      }
      /* read: flag for output on not iced vertices*/
      if (fgets(inputstring, 40, ptFile)==NULL){ 
	printf("WARNING: Error while reading .sico2elmer.rc . Writing all output-components\n");
	failed=1;
      }else{
	if (sscanf(inputstring, "%1d%s",&outputflags[1], dummy_string)<0){
	  printf("WARNING: Error while copying inputstring no 0 from .sico2elmer.rc\n");
	  outputflags[1]=0;
	}
      }
      if (fgets(inputstring, 40, ptFile)==NULL){ 
	printf("WARNING: Error while reading .sico2elmer.rc. Writing all output-components\n");
	failed=1;
      }
      /* read: 3d variable flags */
      for (n=0;n<4;++n){
	if (fgets(inputstring, 40, ptFile)==NULL){ 
	  printf("WARNING: Error while reading .sico2elmer.rc. Writing all output-components\n");
	  failed=1;
	}else{
	  if (sscanf(inputstring, "%1d%s",&outputflags[2+n], dummy_string)<0){
	    printf("WARNING: Error while copying inputstring no %d from .sico2elmer.rc\n", n+1);
	    outputflags[2+n]=1;
	  }
	}
      }
      if (fgets(inputstring, 40, ptFile)==NULL){ 
	printf("WARNING: Error while reading .sico2elmer.rc. Writing all output-components\n");
	failed=1;
      } 
      /* read: 2d variable flags */
      for (n=0;n<6;++n){
	if (fgets(inputstring, 40, ptFile)==NULL){ 
	  printf("WARNING: Error while reading .sico2elmer.rc. Writing all output-components\n");
	  failed=1;
	}else{
	  if (sscanf(inputstring, "%1d%s",&outputflags[6+n], dummy_string)<0){
	    printf("WARNING: Error while copying inputstring no %d from .sico2elmer.rc\n", n+1);
	    outputflags[6+n]=1;
	  }
	}
      }
    }
    fclose(ptFile);
  }
  if (failed){
    for (n=0;n<12;++n) outputflags[n]=1;
    printf("| Error occured during reading of file\n");
    printf("| Taking default values\n");
  }else printf("| Loading of preferences succeessful\n");
  printf("| Following values for output are set (1=yes,0=no):\n");
  printf("| \n");
  if (outputflags[0] && ((*flag)==0)){
    outputflags[0] = 0;
    printf("|  Output temperate layer       %d (reset, since no temperate layer has been loaded)\n", outputflags[0]);    
  }else
    printf("|  Output temperate layer       %d\n", outputflags[0]);
  printf("|  Output at not iced vertices  %d\n", outputflags[1]);
  printf("| \n");
  printf("|  Output 3d coordinates        %d\n", outputflags[2]);
  printf("|  Output velocities            %d\n", outputflags[3]);
  printf("|  Output temperature           %d\n", outputflags[4]);
  printf("|  Output age                   %d\n", outputflags[5]);
  printf("| \n");
  printf("|  Output bedrock and ice depht %d\n", outputflags[6]);
  printf("|  Output water drainage        %d\n", outputflags[7]);
  printf("|  Output mask for glaciation   %d\n", outputflags[8]);
  printf("|  Output type of base          %d\n", outputflags[9]);
  printf("|  Output flux                  %d\n", outputflags[10]);
  printf("|  Output 2d coordinates        %d\n", outputflags[11]);
  printf("---------------------------------------------------------------\n");
  /* open file for writing of 3d variables */
  sprintf(filename,"%5s%2s_3d%s", runname, ergnum, suffix);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open file for writing of 3d ASCII-data\n");
    return;
  }else{
    printf("| writing on 3d output file %s\n", filename);
    fprintf(ptFile,"# *********************************************************\n");
    fprintf(ptFile,"# ASCCI output file of run %s, file no. %s\n", runname, ergnum);
    if (time(&how_late_is_it)==-1){
      printf("WARNING: Could not evaluate current time\n");
      fprintf(ptFile,"# no specific date was able to be inquiered\n");
    }else{      
      fprintf(ptFile,"# written on  %s", ctime(&how_late_is_it));
    }
    /* replace with getpwuid(geteuid())  */
#if defined(__APPLE__) || defined(MINGW32) || defined(WIN32) || defined(BSD)
    if ( 1 )
#else
    if (cuserid(user_name)==NULL)
#endif
      fprintf(ptFile,"# no user id was able to be inquiered\n");
    else 
      fprintf(ptFile,"# by %10s\n", user_name);
    fprintf(ptFile,"# *********************************************************\n");
    fprintf(ptFile,"# 3d variables:\n");
    fprintf(ptFile,"# ");
    if (outputflags[1])
      fprintf(ptFile,"x           y           z           ");
    if (outputflags[3])
      fprintf(ptFile,"v_x         v_y         v_z         ");
    if (outputflags[4])
      fprintf(ptFile,"T           ");
    if (outputflags[5])
      fprintf(ptFile,"Age         ");
    if ((outputflags[2] + outputflags[3] + outputflags[4] + outputflags[5]) !=0){
      fprintf(ptFile,"\n# --------- temperate layer --------------------------------------------------------------------\n");      
      if (outputflags[0]){
	for (k=0;k<ktmax;++k){
	  for (j=0;j<jmax+1;++j){
	    for (i=0;i<imax+1;++i){
	      if( (outputflags[1]) && !((maske[j*(imax+1) + i]==0) || (maske[j*(imax+1) + i]==3)) ) continue;
	      current_index = k*elements_in_layer + j*(imax+1) + i;
	      if (outputflags[2])
		fprintf(ptFile,"% .4e % .4e % .4e ",xi[i], eta[j], z_t[current_index]);
	      if (outputflags[3])
		fprintf(ptFile,"% .4e % .4e % .4e ",vx_t[current_index], vy_t[current_index], vz_t[current_index]);
	      if (outputflags[4])
		fprintf(ptFile,"% .4e ",temp_t_m[current_index]);
	      if (outputflags[5])
		fprintf(ptFile,"% .4e ",age_t[current_index]); 
	      fprintf(ptFile,"\n");
	    }
	    if (BLANK) fprintf(ptFile,"\n\n");
	  }
	}
      }
      fprintf(ptFile,"# --------- cold layer -----------------------------------------------------------------------\n");
      for (k=0;k<kcmax+1;++k){
	for (j=0;j<jmax+1;++j){
	  for (i=0;i<imax+1;++i){	
	    if( (outputflags[1]) && !((maske[j*(imax+1) + i]==0) || (maske[j*(imax+1) + i]==3)) ) continue;
	    current_index = k*elements_in_layer + j*(imax+1) + i;
	    if (outputflags[2])
	      fprintf(ptFile,"%+.4e %+.4e %+.4e ",xi[i], eta[j], z_c[current_index]);
	    if (outputflags[3])
	      fprintf(ptFile,"% .4e % .4e % .4e ",vx_c[current_index], vy_c[current_index], vz_c[current_index]);
	    if (outputflags[4])
	      fprintf(ptFile,"% .4e ",temp_c[current_index]);
	    if (outputflags[5])
	      fprintf(ptFile,"% .4e ",age_c[current_index]); 
	    fprintf(ptFile,"\n");
	  }
	  if (BLANK) fprintf(ptFile,"\n\n\n");
	  else fprintf(ptFile,"\n");	
	}
      }
    }
    fclose(ptFile);
  }
  /* open file for writing of 2d variables */
  sprintf(filename,"%5s%2s_2d%s", runname, ergnum, suffix);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open file for writing of 3d ASCII-data\n");
    return;
  }else{
    printf("| writing on 2d output file %s\n", filename);
    fprintf(ptFile,"# *********************************************************\n");
    fprintf(ptFile,"# ASCCI output file of run %s, file no. %s\n", runname, ergnum);
    if (time(&how_late_is_it)==-1){
      printf("WARNING: Could not evaluate current time\n");
      fprintf(ptFile,"# no specific date was able to be inquiered\n");
    }else{      
      fprintf(ptFile,"# written on  %s", ctime(&how_late_is_it));
    }

#if defined(__APPLE__) || defined(MINGW32) || defined(WIN32) || defined(BSD)
    if ( 1 )
#else
    if (cuserid(user_name)==NULL)
#endif
      fprintf(ptFile,"# no user id was able to be inquiered\n");
    else 
      fprintf(ptFile,"# by %10s\n", user_name);
    fprintf(ptFile,"# *********************************************************\n");
    fprintf(ptFile,"# 2d variables:\n");
    fprintf(ptFile,"# ");
    if (outputflags[11])
      fprintf(ptFile,"x           y           ");
    if (outputflags[6])
      fprintf(ptFile,"b           s           ");
    if (outputflags[7])
      fprintf(ptFile,"Drainage    ");
    if (outputflags[8])
      fprintf(ptFile,"glaciation  ");
    if (outputflags[9])
      fprintf(ptFile,"type_base   ");
    if (outputflags[10])
      fprintf(ptFile,"flux_x      flux_y      ");
    fprintf(ptFile,"\n# --------------------------------------------------------------------------------------------------------\n");
    for (j=0;j<jmax+1;++j){
      for (i=0;i<imax+1;++i){
	current_index = j*(imax+1) + i;
/* 	printf("%d ++ %d ==  %d,  %d\n", outputflags[1], !((maske[current_index]==0) || (maske[current_index]==3)), ( (outputflags[1]==0) && !((maske[current_index]==0) || (maske[current_index]==3)) ), maske[current_index]); */
	if( (outputflags[1]==0) && !((maske[current_index]==0) || (maske[current_index]==3)) ) continue;
	if (outputflags[11])
	  fprintf(ptFile,"% .4e % .4e", xi[i], eta[j]);
	if (outputflags[6]){
/* 	  printf("%d %d\n", i,j); */
	  if (outputflags[0]){
	    fprintf(ptFile,"% .4e % .4e", z_t[current_index], z_c[kcmax*elements_in_layer + current_index]);
	  }
	  else
	    fprintf(ptFile,"% .4e % .4e", z_c[current_index], z_c[kcmax*elements_in_layer + current_index]);
	}
	if (outputflags[7])
	  fprintf(ptFile,"% .4e ",Q_tld[current_index]);
	if (outputflags[8])
	  fprintf(ptFile,"% 1d          ",maske[current_index]);
	if (outputflags[9])
	  fprintf(ptFile,"% 1d          ",n_cts[current_index]);
	if (outputflags[10])
	  fprintf(ptFile,"% .4e % .4e",qx[current_index],qy[current_index]);
	fprintf(ptFile,"\n");
      }
      if (BLANK) fprintf(ptFile,"\n\n");
    }
    fclose(ptFile);
  }
  printf("---------------------------------------------------------------\n");
}
/*******************************************************
              calculate staggered grid 
********************************************************/
int get_staggered_grid(float  *xi, /* unscaled x coordinate index i: 0,imax */
		       float  *eta, /* unscaled y coordinate index j from 0,jmax */
		       float  *z_in, /* unscaled z coordinate index i: 0,imax, j: 0,jmax, kc: 0,kcmax */
		       int    imax, /* grid steps in xi-direction */
		       int    jmax, /* grid steps in eta-direction */
		       int    kmax, /* grid steps in z-direction in cold ice layer */
		       float  *deltaX, /* stepsize of original grid */
		       float  *staggered_grid){ /* coords of staggered FEM grid */
  
  register int i,j,k,n;
  int nodes_in_layer_of_staggered_grid, nodes_in_staggered_grid, auxiliary;
  float delta_x, *x,*y,*z;

  nodes_in_layer_of_staggered_grid = (imax+2)*(jmax+2);
  nodes_in_staggered_grid =nodes_in_layer_of_staggered_grid*(kmax+1);

  x=(float *) malloc((size_t) (imax+2)*sizeof(float));
  y=(float *) malloc((size_t) (jmax+2)*sizeof(float));
  z=(float *) malloc((size_t) nodes_in_staggered_grid*sizeof(float));
  
  delta_x= 0.5*(*deltaX);

  for (i=0;i<imax+1;++i)
    x[i]=xi[i] - delta_x;
  x[imax+1]=xi[imax] + delta_x;
  for (j=0;j<jmax+1;++j)
    y[j]=eta[j] - delta_x;
  y[jmax+1]=eta[jmax] + delta_x;
  
  auxiliary = get_interpolated_property_on_staggered_grid(imax, jmax, kmax, z_in, z);
  if (auxiliary != nodes_in_staggered_grid){
    printf("WARNING: Returned value of interpolated z-values %d does not match calculated %d\n", auxiliary, nodes_in_staggered_grid); 
    free(x);free(y);free(z);
    return -1;
  }
  for (k=0,n=0;k<kmax+1;++k){
    for(j=0;j<jmax+2;++j){
      for (i=0;i<imax+2;++i,++n){
	staggered_grid[(k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i)*3]=x[i];
	staggered_grid[(k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i)*3+1]=y[j];
	staggered_grid[(k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i)*3+2]=z[k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i];
/* 	printf("%f %f %f\n", staggered_grid[(k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i)*3], staggered_grid[(k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i)*3+1], staggered_grid[(k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i)*3 + 2]); */
      }
    }
/*     printf("\n"); */
  }
  free(x);free(y);free(z);
  return n; 
}

/*******************************************************
   interpolate property on staggered grid
********************************************************/
int get_interpolated_property_on_staggered_grid(int imax,
						int jmax,
						int kmax,
						float *property_in,
						float *property_out){  
  register int i,j,k,n;
  int nodes_in_layer_of_staggered_grid, nodes_in_layer;

  nodes_in_layer = (imax+1)*(jmax+1);
  nodes_in_layer_of_staggered_grid = (imax+2)*(jmax+2);
  
  /* inside array */
  for (k=0,n=0;k<kmax+1;++k){
    for (j=1;j<jmax+1;++j){
      for (i=1;i<imax+1;++i){
	property_out[k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i]=0.25*(property_in[k*nodes_in_layer + j*(imax+1) + i-1] + property_in[k*nodes_in_layer + (j-1)*(imax+1) + i-1] + property_in[k*nodes_in_layer + j*(imax+1) + i] + property_in[k*nodes_in_layer + (j-1)*(imax+1) + i]);
	++n;
      }
    }
  }      
  /* at borders */
  for (k=0;k<kmax+1;++k){
    for (j=1;j<jmax+1;++j){
      property_out[k*nodes_in_layer_of_staggered_grid + j*(imax+2)]=0.5*(property_in[k*nodes_in_layer + j*(imax+1)] + property_in[k*nodes_in_layer + (j-1)*(imax+1)]);
      ++n;
      property_out[k*nodes_in_layer_of_staggered_grid + j*(imax+2) + imax+1]=0.5*(property_in[k*nodes_in_layer + j*(imax+1) + imax] + property_in[k*nodes_in_layer + (j-1)*(imax+1) + imax]);
      ++n;
    }
    for (i=1;i<imax+1;++i){
      property_out[k*nodes_in_layer_of_staggered_grid + i]= 0.5*(property_in[k*nodes_in_layer + i] + property_in[k*nodes_in_layer + i-1]);
      ++n;
      property_out[k*nodes_in_layer_of_staggered_grid + (jmax+1)*(imax+2) + i]=0.5*(property_in[k*nodes_in_layer + jmax*(imax+1) + i] + property_in[k*nodes_in_layer + jmax*(imax+1) + i-1]);
      ++n;
    }
  }
  /* at corners */
  for (k=0;k<kmax+1;++k){
    property_out[k*nodes_in_layer_of_staggered_grid] = property_in[k*nodes_in_layer];
    ++n;
    property_out[k*nodes_in_layer_of_staggered_grid + imax+1] =  property_in[k*nodes_in_layer+ imax];
    ++n;
    property_out[k*nodes_in_layer_of_staggered_grid + (jmax+1)*(imax+2)] = property_in[k*nodes_in_layer+jmax*(imax+1)];
    ++n;
    property_out[k*nodes_in_layer_of_staggered_grid + j*(imax+2) + i] = property_in[k*nodes_in_layer+jmax*(imax+1)+imax];
    ++n;
  }
  return n;
}
/************************************************
   get information on glaciation of ice-sheet(s)
*************************************************/
int get_glaciation_info(int imax,
			int jmax,
			int *iced,
			int *maske){
  register int i,j,n;

  for (j=0, n=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      if ((maske[j*(imax+1) + i]==0) || (maske[j*(imax+1) + i]==3)){
	iced[j*(imax+1) + i]=n;
	++n;
      }else{
	iced[j*(imax+1) + i]=-1;  
      }
    }
  }  
  return n;
}
/************************************************
   get 2d glaciation map for staggered grid
*************************************************/
int get_grid_map(int imax,
		 int jmax,
		 int *iced,
		 int *gridmap){
  register int i,j,n;
  int nodes_in_layer_of_staggered_grid, nodes_in_layer, *staggered_iced;
  nodes_in_layer = (imax+1)*(jmax+1);
  nodes_in_layer_of_staggered_grid = (imax+2)*(jmax+2);
  if ((staggered_iced = (int *) malloc((size_t)  nodes_in_layer_of_staggered_grid*sizeof(int)))==NULL){
    printf("ERROR during allocation of memory for glaciation information of staggered grid\n");
    return -1;
  }  
  for (j=0;j<jmax+2;++j){
    for (i=0;i<imax+2;++i){
      staggered_iced[j*(imax+2) + i]=0;
    }
  }
  for (j=0;j<jmax+1;++j){
    for (i=0;i<imax+1;++i){
      if (iced[j*(imax+1) + i]!=-1){
	staggered_iced[j*(imax+2) + i]=1;
	staggered_iced[(j+1)*(imax+2) + i]=1;
	staggered_iced[j*(imax+2) + i+1]=1;
	staggered_iced[(j+1)*(imax+2) + i+1]=1;
      }
    }
  }
  for (j=0, n=0;j<jmax+2;++j){
    for (i=0;i<imax+2;++i){
      if (staggered_iced[j*(imax+2) + i]){
	gridmap[j*(imax+2) + i]=n;
	++n;
      }else{
	gridmap[j*(imax+2) + i]=-1;
      }
    }
  }
  free(staggered_iced);
  return n;
}
/**************************************************
   get 2d information on boundaries of ice-sheet(s)
***************************************************/
int get_glaciation_boundary_info(int imax,
				 int jmax,
				 int *iced,
				 int *boundary){
  register int i,j,n,nn;
  /* inside array */
  for (j=1,n=0,nn=0;j<jmax;++j){
    for (i=1;i<imax;++i){
      if (iced[j*(imax+1) + i]!=-1){
	if (iced[j*(imax+1) + i+1]==-1){ /* east */
	  boundary[iced[j*(imax+1) + i]*4+0]=1;
	  ++n;
	}else boundary[iced[j*(imax+1) + i]*4+0]=0;
	if (iced[(j+1)*(imax+1) + i]==-1){ /* north */
	  boundary[iced[j*(imax+1) + i]*4+1]=1;
	  ++n;
	}else boundary[iced[j*(imax+1) + i]*4+1]=0;
	if (iced[j*(imax+1) + i-1]==-1){ /* west */
	  boundary[iced[j*(imax+1) + i]*4+2]=1;
	  ++n;
	}else boundary[iced[j*(imax+1) + i]*4+2]=0;
	if (iced[(j-1)*(imax+1) + i]==-1){ /* south */
	  boundary[iced[j*(imax+1) + i]*4+3]=1;
	  ++n;
	}else boundary[iced[j*(imax+1) + i]*4+3]=0;
      }
    }
  }
  /* at borders */
  for (j=1;j<jmax;++j){
    /* most western row (i==0)*/
    if (iced[j*(imax+1)]!=-1){
      boundary[iced[j*(imax+1)]*4+2]=1; /* west is always boundary */
      ++n;
      if (iced[j*(imax+1) +1]==-1){ /* east */
	boundary[iced[j*(imax+1)]*4+0]=1;
	++n;
      }else boundary[iced[j*(imax+1)]*4+0]=0;
      if (iced[(j+1)*(imax+1)]==-1){ /* north */
	boundary[iced[j*(imax+1)]*4+1]=1;
	++n;
      }else boundary[iced[j*(imax+1)]*4+1]=0;
      if (iced[(j-1)*(imax+1)]==-1){ /* south */
	boundary[iced[j*(imax+1)]*4+3]=1;
	++n;
      }else boundary[iced[j*(imax+1)]*4+3]=0;
    }
    /* most eastern row (i==imax)*/
    if (iced[j*(imax+1)+imax]!=-1){
      boundary[iced[j*(imax+1) + imax]*4+0]=1;/* east is always boundary */
      ++n;
      if (iced[(j+1)*(imax+1) + imax]==-1){ /* north */
	boundary[iced[j*(imax+1) + imax]*4+1]=1;
	++n;
      }else boundary[iced[j*(imax+1) + imax]*4+1]=0;
      if (iced[j*(imax+1) + imax-1]==-1){ /* west */
	boundary[iced[j*(imax+1) + imax]*4+2]=1;
	++n;
      }else boundary[iced[j*(imax+1) + imax]*4+2]=0;
      if (iced[(j-1)*(imax+1) + imax]==-1){ /* south */
	boundary[iced[j*(imax+1) + imax]*4+3]=1;
	++n;
      }else boundary[iced[j*(imax+1) + imax]*4+3]=0;
    }
  }
  for (i=1;i<imax;++i){
    /* most northern row (j==jmax)*/
    if (iced[jmax*(imax+1) + i]!=-1){
      boundary[iced[jmax*(imax+1) + i]*4+1]=1;/* north is always boundary */
      ++n;
      if (iced[jmax*(imax+1) + i+1]==-1){ /* east */
	boundary[iced[jmax*(imax+1) + i]*4+0]=1;
	++n;
      }else boundary[iced[jmax*(imax+1) + i]*4+0]=0;
      if (iced[jmax*(imax+1) + i-1]==-1){ /* west */
	boundary[iced[jmax*(imax+1) + i]*4+2]=1;
	++n;
      }else boundary[iced[jmax*(imax+1) + i]*4+2]=0;
      if (iced[(jmax-1)*(imax+1) + i]==-1){ /* south */
	boundary[iced[jmax*(imax+1) + i]*4+3]=1;
	++n;
      }else boundary[iced[jmax*(imax+1) + i]*4+3]=0;
    }
    /* most southern (j==0)*/
    if (iced[i]!=-1){
      boundary[iced[i]*4+3]=1; /* south is always boundary */
      ++n;
      if (iced[i+1]==-1){ /* east */
	boundary[iced[i]*4+0]=1;
	++n;
      }else boundary[iced[i]*4+0]=0;
      if (iced[1*(imax+1) + i]==-1){ /* north */
	boundary[iced[i]*4+1]=1;
	++n;
      }else boundary[iced[i]*4+1]=0;
      if (iced[i-1]==-1){ /* west */
	boundary[iced[i]*4+2]=1;
	++n;
      }else boundary[iced[i]*4+2]=0;
    }
  }
  /* at corners */
  /* northeast */
  if (iced[jmax*(imax+1) + imax]!=-1){
    boundary[iced[jmax*(imax+1) + imax]*4+0]=1; /* east is always boundary */
    ++n;
    boundary[iced[jmax*(imax+1) + i]*4+1]=1;/* north is always boundary */
    ++n;
    if (iced[jmax*(imax+1) + imax-1]==-1){ /* west */
      boundary[iced[jmax*(imax+1) + imax]*4+2]=1;
      ++n;
    }else boundary[iced[jmax*(imax+1) + imax]*4+2]=0;
    if (iced[(jmax-1)*(imax+1) + imax]==-1){ /* south */
      boundary[iced[jmax*(imax+1) + imax]*4+3]=1;
      ++n;
    }else boundary[iced[jmax*(imax+1) + imax]*4+3]=0;
  }
  /* northwest */
  if (iced[jmax*(imax+1)]!=-1){
    if (iced[jmax*(imax+1)+1]==-1){ /* east */
      boundary[iced[jmax*(imax+1)]*4+0]=1;
      ++n;
    }else boundary[iced[jmax*(imax+1)]*4+0]=0;
    boundary[iced[jmax*(imax+1)]*4+1]=1;/* north is always boundary */
    ++n;
    boundary[iced[jmax*(imax+1)]*4+2]=1;/* west is always boundary */
    ++n;
    if (iced[(jmax-1)*(imax+1)]==-1){ /* south */
      boundary[iced[jmax*(imax+1)]*4+3]=1;
      ++n;
    }else boundary[iced[jmax*(imax+1)]*4+3]=0;
  }
  /* southwest */
  if (iced[0]!=-1){
    if (iced[1]==-1){ /* east */
      boundary[iced[0]*4+0]=1;
      ++n;
    }else boundary[iced[0]*4+0]=0;
    if (iced[1*(imax+1)]==-1){ /* north */
      boundary[iced[0]*4+1]=1;
      ++n;
    }else boundary[iced[0]*4+1]=0;
    boundary[iced[0]*4+2]=1;/* west is always boundary */
    ++n;
    boundary[iced[0]*4+3]=1;/* south is always boundary */
    ++n;
  }
  /* southeast */
  if (iced[imax]!=-1){
    boundary[iced[imax]*4+0]=1;/* east is always boundary */
    ++n;
    if (iced[1*(imax+1) + imax]==-1){ /* north */
      boundary[iced[imax]*4+1]=1;
      ++n;
    }else boundary[iced[imax]*4+1]=0;
    if (iced[imax-1]==-1){ /* west */
      boundary[iced[imax]*4+2]=1;
      ++n;
    }else boundary[iced[imax]*4+2]=0;
    boundary[iced[imax]*4+3]=1;/* south is always boundary */
    ++n;
  }
  return n;
}
/*******************************************************
              transform integer into float array
********************************************************/
void  make_float_from_integer_scalar_field(int   *input_property,
					   float *output_property, 
					   int   number_of_nodes,
					   int   reorder_ice_land_sea_mask){
  register int n;
  if (reorder_ice_land_sea_mask){
    for (n=0;n<number_of_nodes;++n){
      if ((input_property[n]==0) || (input_property[n]==3)) output_property[n] = 0; /* ice */
      else if (input_property[n]==1) output_property[n] = 1; /* land */
      else output_property[n] = 2; /* sea */
    }
  }else{
    for (n=0;n<number_of_nodes;++n)
      output_property[n] = (float) input_property[n];    
  }
  return;
}






/****************************************************************
   Output of SICOPOLIS grid (not staggerd) in ELMER Solver format
*****************************************************************/
void STDCALLBULL FC_FUNC(solvergrid,SOLVERGRID)(float  *xi, /* unscaled x coordinate index i: 0,imax */
				    float  *eta, /* unscaled y coordinate index j from 0,jmax */
				    float  *z_c, /* unscaled z coordinate index i: 0,imax, j: 0,jmax, kc: 0,kcmax */
				    float  *z_t, /* unscaled y coordinate index i: 0,imax, j: 0,jmax, kt: 0,kcmax */
				    int    *imax_in, /* grid steps in xi-direction */
				    int    *jmax_in, /* grid steps in eta-direction */
				    int    *kcmax_in, /* grid steps in z-direction in cold ice layer */
				    int    *ktmax_in, /* grid steps in z-direction in temperate ice layer */
				    FC_CHAR_PTR(runname,runname_l), /*name of run*/
				    FC_CHAR_PTR(ergnum,ergnum_l), /*number of file*/
				    int    *maske, /*mask of vertex type */
				    float  *deltaX, /* stepsize of grid */
				    int    *flag)
{
  register int i, j, k, l, m, n;
  int   number_of_bulk_elements, number_of_boundary_elements, number_of_iced_collums, elements_in_one_layer, number_of_nodes, nodes_of_element[8], *nodes_of_side_element, *parent_of_side_element, nodes_in_one_layer, *iced, idnr, parent_element, sidebulk=0;
  int   imax, jmax, kcmax, ktmax, kmax = 0;
  float  min_max_xyz[2][3], *bottom, freesurf, *delta_z;
  float actual_scaled_coord[3];
  char  groupid[4], filename[80], yes_no;
  FILE  *ptFile;
  
  /* constants */
  imax= *imax_in;
  jmax= *jmax_in;
  kcmax= *kcmax_in;
  ktmax= *ktmax_in;

  /* print out little summary */
  printf("---------------------------------------------------------------\n");
  printf("|        Output of SICOPOLIS Grid for ELMER Solver\n");
  printf("---------------------------------------------------------------\n");
  printf("| imax/jmax/kcmax/ktmax=%d/%d/%d/%d\n",imax, jmax, kcmax, ktmax);
  printf("---------------------------------------------------------------\n");
  printf("| nodes in original grid:\n");
  printf("|       cold layer:  %d=%d * %d\n", (imax+1)*(jmax+1)*(kcmax+1), (imax+1)*(jmax+1), (kcmax+1));
  printf("|  temperate layer:  %d=%d * (%d - 1)\n", (imax+1)*(jmax+1)*(ktmax), (imax+1)*(jmax+1), (ktmax+1));
  printf("|                    -------------\n");
  printf("|                    %d=%d * %d\n", (imax+1)*(jmax+1)*(ktmax+1+kcmax), (imax+1)*(jmax+1), ktmax+1+kcmax);
  printf("---------------------------------------------------------------\n");
  printf("| elements in original grid:\n");
  printf("|       cold layer:  %d\n", (imax)*(jmax)*kcmax);
  printf("|  temperate layer: %d\n", (*flag)*(imax)*(jmax)*(ktmax));
  printf("---------------------------------------------------------------\n");
  printf("| x = %.1f -> %.1f , y = %.1f -> %.1f, dx=%.1f\n", xi[0], xi[imax], eta[0], eta[jmax], *deltaX);
  printf("---------------------------------------------------------------\n");

  /* inquire number of grid levels */
  /* ------------------------------*/
  while (kmax < 1){
    printf("| How many vertical grid layers the grid shall contain? \n");
    scanf("%d", &kmax);
    printf("| \n");
    if (kmax < 1)  printf("| No. of vertical (current value %d) levels must exceed 1!\n", kmax);
  }
  printf("Thanks. Taking %d  vertical layers\n", kmax);
  printf("---------------------------------------------------------------\n");

  /* inquire if sidebulk on free surface is desired */
  /* -----------------------------------------------*/
  printf("| sidebulk elements on free surface (y/n)? \n");
  scanf("%1s", &yes_no);
  printf("| \n");
  if ((yes_no == 'y') || (yes_no == 'Y') ){
    printf("| Thanks. Writing sidebulk elements on free surface\n");
    sidebulk = 1;
  }else{
    printf("| Thanks. No output of sidebulk elements on free surface\n");
    sidebulk = 0;
  }
  printf("---------------------------------------------------------------\n");

  /* calculate constants for output mesh*/
  /* -----------------------------------*/
  nodes_in_one_layer = (imax+1)*(jmax+1);
  elements_in_one_layer = (imax)*(jmax);
  number_of_nodes =  nodes_in_one_layer*(kmax+1);
  number_of_bulk_elements =  elements_in_one_layer*kmax;
  number_of_boundary_elements = 2*elements_in_one_layer + 2*(imax+jmax)*kmax;

  /* allocate arrays */
  /* ----------------*/
  iced = (int *) malloc((size_t) (imax+1)*(jmax+1)*sizeof(int));
  if (iced == NULL){
    printf("ERROR in allocating memory for glaciation information\n");    
    free(iced);
    return;
  }  
  delta_z = (float *) malloc((size_t) (imax+1)*(jmax+1)*sizeof(float));
  if (delta_z == NULL){
    printf("ERROR in allocating memory for vertical grid step sizes\n");    
    free(iced);free(delta_z);
    return;
  }
  bottom = (float *) malloc((size_t) (imax+1)*(jmax+1)*sizeof(float));
  if (bottom == NULL){
    printf("ERROR in allocating memory for bottom topography\n");    
    free(iced);free(delta_z);free(bottom);
    return;
  }
  
  /* get glaciation info */
  /* --------------------*/  
  number_of_iced_collums = get_glaciation_info(imax,jmax,iced,maske);
  if (number_of_iced_collums<1){
    printf("WARNING: domain seems to be ice-free\n");
  }
  printf("| number of iced colums:       %d out of %d (%3.2f %%)\n", number_of_iced_collums, (imax+1)*(jmax+1), ((float) number_of_iced_collums)*100.0/((float) (imax+1)*(jmax+1)));
  printf("---------------------------------------------------------------\n");
  printf("| number of nodes:             %d=%d * %d\n",number_of_nodes, nodes_in_one_layer, kmax+1);
  printf("| number of bulk elements:     %d=%d * %d\n",number_of_bulk_elements, elements_in_one_layer, kmax);
  printf("| number of sidebulk elements: %d\n", sidebulk*elements_in_one_layer);
  printf("| number of boundary elements: %d=2*[%d + (%d + %d) *%d] + %d\n", number_of_boundary_elements, elements_in_one_layer, imax, jmax, kmax, 2*(imax+jmax)*sidebulk);
  printf("| sidebulk boundary elements:  %d=%d*2*(%d + %d)\n", sidebulk*2*(imax+jmax), sidebulk, imax, jmax);
  printf("---------------------------------------------------------------\n");

  /* writing header file */
  /* --------------------*/
  sprintf(filename,"mesh.header");
  printf("| Writing mesh header file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(iced);free(delta_z);free(bottom);
    return;
  }
  fprintf(ptFile, "%d %d %d\n", number_of_nodes,  number_of_bulk_elements + sidebulk*elements_in_one_layer, number_of_boundary_elements + sidebulk*2*(imax+jmax));
  fprintf(ptFile,"%d\n",2+sidebulk);
  if (sidebulk)
    fprintf(ptFile,"202 %d\n", 2*(imax+jmax));
  fprintf(ptFile,"404 %d\n", number_of_boundary_elements + sidebulk*elements_in_one_layer);
  fprintf(ptFile,"808 %d\n", number_of_bulk_elements);

  printf("| succeeded in writting mesh header file %s.\n",filename);
  printf("---------------------------------------------------------------\n");
  fclose(ptFile);

  /* init min/max */
  /* -------------*/
  min_max_xyz[0][2]= min_max_xyz[1][2] = z_c[0];
  min_max_xyz[0][0]= xi[0];
  min_max_xyz[1][0] = xi[imax];
  min_max_xyz[0][1]=  eta[0];
  min_max_xyz[1][1] = eta[jmax];

  /* get delta_z for specific i,j value */
  /* -----------------------------------*/
  for (j=0;j<jmax+1;++j){ /* loop over all j-values */
    for (i=0;i<imax+1;++i) { /* loop over all i-values */
      freesurf = z_c[kcmax*nodes_in_one_layer + j*(imax+1) + i];
      if (*flag){
	bottom[j*(imax+1) + i] = z_t[j*(imax+1) + i]; 
      }else{
	bottom[j*(imax+1) + i] = z_c[j*(imax+1) + i];
      }
      if (min_max_xyz[0][2]>bottom[j*(imax+1) + i]) min_max_xyz[0][2]=bottom[j*(imax+1) + i];
      if (min_max_xyz[1][2]<freesurf) min_max_xyz[0][2]=freesurf;
      delta_z[j*(imax+1) + i] = (freesurf - bottom[j*(imax+1) + i])/((float) kmax);
/*       printf("dz=%f = %f - %f \n",   delta_z[j*(imax+1) + i], freesurf, bottom[j*(imax+1) + i]); */
/*       printf("%d %d % f\n", i,j,delta_z[j*(imax+1) + i]); */
      if (delta_z[j*(imax+1) + i] < 0.0){
	printf("\a delta z = %f - %f for (%d, %d) = %f < 0\n", freesurf, bottom[j*(imax+1) + i], i, j, delta_z[j*(imax+1) + i]);
	free(iced);free(delta_z);free(bottom);
	fclose(ptFile);
	return;
      }else if((delta_z[j*(imax+1) + i] < MINHEIGHT/((float) kmax))|| ((maske[j*(imax+1) + i]!=0) && (maske[j*(imax+1) + i]!=3))){
	delta_z[j*(imax+1) + i] = MINHEIGHT/((float) kmax);
      }
    }
  }
  printf("| succeeded in calculating grid step sizes %s.\n",filename);
  printf("---------------------------------------------------------------\n");

  /* writing nodes file */
  /* -------------------*/
  sprintf(filename,"mesh.nodes");
  printf("| Writing node-file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(iced);free(delta_z);free(bottom);
    return;
  }
  for (n=0,k=0;k<kmax+1;++k){ /* loop over all levels */
    for (j=0;j<jmax+1;++j){ /* loop over all j-values */
      for (i=0;i<imax+1;++i) { /* loop over all i-values */
	++n;
	fprintf(ptFile, "%d -1 %.8f %.8f %.8f\n", n, eta[i], xi[j], bottom[j*(imax+1) + i] + delta_z[j*(imax+1) + i] * ((float) k));
      }
    }
  }
  fclose(ptFile);
  if (n != number_of_nodes){
    printf("\a %d written nodes does not match %d calculated nodes in grid\n", n, number_of_nodes);
    return;	   
  }else{
    printf("| succeeded in writing node file %s.\n",filename);
    printf("---------------------------------------------------------------\n");
  }
   
  /* writing element file */
  /* ---------------------*/
  sprintf(filename,"mesh.elements");
  printf("| Writing bulk element file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer mesh-file file for writing!\n");
    free(iced);free(delta_z);free(bottom);
    return;
  }
  /* write bulk elements */
  idnr = 1;
  for (n=0,k=0;k<kmax;++k){ /* loop over all levels */
    for (j=0;j<jmax;++j){ /* loop over all j-values */
      for (i=0;i<imax;++i) { /* loop over all i-values */
	n++;
	for (l=0; l<2; ++l){/* lower level: n=0; upper level n=1; each counterclkws */
	  nodes_of_element[l*4] = (k+l)*nodes_in_one_layer + j*(imax+1) + i;
	  nodes_of_element[l*4+1] = (k+l)*nodes_in_one_layer + j*(imax+1) + i+1;
	  nodes_of_element[l*4+2] = (k+l)*nodes_in_one_layer + (j+1)*(imax+1) + i+1;
	  nodes_of_element[l*4+3] = (k+l)*nodes_in_one_layer + (j+1)*(imax+1) + i;
	}
	fprintf(ptFile,"%d %d 808 %d %d %d %d %d %d %d %d\n", n, idnr,
		nodes_of_element[0]+1, nodes_of_element[1]+1, nodes_of_element[2]+1, nodes_of_element[3]+1, 
		nodes_of_element[4]+1, nodes_of_element[5]+1, nodes_of_element[6]+1, nodes_of_element[7]+1);
      }
    }
  }
  /* write sidebulk elements (if demanded) */
  if (sidebulk){
    idnr = 2;
    for (m=0,j=0;j<jmax;++j){ /* loop over all j-values */
      for (i=0;i<imax;++i) { /* loop over all i-values */
	++m;
	nodes_of_element[0] = kmax*nodes_in_one_layer + j*(imax+1) + i;
	nodes_of_element[1] = kmax*nodes_in_one_layer + j*(imax+1) + i+1;
	nodes_of_element[2] = kmax*nodes_in_one_layer + (j+1)*(imax+1) + i+1;
	nodes_of_element[3] = kmax*nodes_in_one_layer + (j+1)*(imax+1) + i;
	fprintf(ptFile,"%d %d 404 %d %d %d %d\n", n+m, idnr,
		nodes_of_element[0]+1, nodes_of_element[1]+1, nodes_of_element[2]+1, nodes_of_element[3]+1);
      }
    }
  }else
    m=0;
  fclose(ptFile);
  if (n+m != number_of_bulk_elements + sidebulk*elements_in_one_layer){
    printf("\a %d=%d + %d  written bulk elements does not match %d= %d + %d  calculated elements in grid\n",
	   n+m, n,m, number_of_bulk_elements+elements_in_one_layer, number_of_bulk_elements, elements_in_one_layer);
    return;	   
  }else{
    printf("| %d number of bulk elements written.\n", n);
    printf("| %d number of sidebulk elements written.\n", m);
    printf("| succeeded in writting bulk element file %s.\n",filename);
    printf("---------------------------------------------------------------\n");
  }
  
  /* writing boundary element file */
  /* ------------------------------*/
  sprintf(filename,"mesh.boundary");
  printf("| Writing boundary element file %s.\n",filename);
  if((ptFile=fopen(filename, "w"))==NULL){
    printf("\a Could not open Elmer boundary element file for writing!\n");
    free(iced);free(delta_z);free(bottom);
    return;
  }
  /* base */
  idnr = 1;
  printf("| lower (base, z=0) boundary for bulk; idnr=%d\n", idnr);
  for (n=0,j=0;j<jmax;++j){ /* loop over all j-values */
    for (i=0;i<imax;++i) { /* loop over all i-values */
      n++;
      nodes_of_element[0] = j*(imax+1) + i;
      nodes_of_element[1] = j*(imax+1) + i+1;
      nodes_of_element[2] = (j+1)*(imax+1) + i+1;
      nodes_of_element[3] = (j+1)*(imax+1) + i;
      parent_element = j*imax + i;
      if (parent_element+1 > number_of_bulk_elements){
	printf("parent element %d > number of bulk elments %d\n",parent_element+1, number_of_bulk_elements);
	return;
      }
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n",
	      n, idnr, parent_element+1, nodes_of_element[0]+1, nodes_of_element[1]+1, nodes_of_element[2]+1, nodes_of_element[3]+1);
    }
  }
  /* free surface */
  idnr = 2;
  printf("| upper (free surface, z=max) boundary for bulk; idnr=%d\n", idnr);
  for (j=0;j<jmax;++j){ /* loop over all j-values */
    for (i=0;i<imax;++i) { /* loop over all i-values */
      n++;
      nodes_of_element[0] = kmax*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = kmax*nodes_in_one_layer + j*(imax+1) + i+1;
      nodes_of_element[2] = kmax*nodes_in_one_layer + (j+1)*(imax+1) + i+1;
      nodes_of_element[3] = kmax*nodes_in_one_layer + (j+1)*(imax+1) + i;
      parent_element = (kmax-1)*elements_in_one_layer + j*imax + i; /* same numbering as in bulk */
      if (parent_element > number_of_bulk_elements){
	printf("parent element %d > number of bulk elments %d\n",parent_element+1, number_of_bulk_elements);
	return;
      }
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n",
	      n, idnr, parent_element+1, nodes_of_element[0]+1, nodes_of_element[1]+1, nodes_of_element[2]+1, nodes_of_element[3]+1);
    }
  }
  /* side faces: */
  /* south */
  idnr = 3;
  j=0;
  printf("| southern (y=0) boundary for bulk; idnr=%d\n", idnr);
  for (k=0;k<kmax;++k){ /* loop over all levels */
    for (i=0;i<imax;++i) { /* loop over all i-values */  
      ++n;
      nodes_of_element[0] = k*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = k*nodes_in_one_layer + j*(imax+1) + i+1;
      nodes_of_element[2] = (k+1)*nodes_in_one_layer + j*(imax+1) + i+1;
      nodes_of_element[3] = (k+1)*nodes_in_one_layer + j*(imax+1) + i;
      parent_element = k*elements_in_one_layer + j*imax + i; /* same numbering as in bulk */
      if (parent_element+1 > number_of_bulk_elements){
	printf("parent element %d > number of bulk elments %d\n",parent_element+1, number_of_bulk_elements);
	return;
      }
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n",
	      n, idnr, parent_element+1, nodes_of_element[0]+1, nodes_of_element[1]+1, nodes_of_element[2]+1, nodes_of_element[3]+1);
    }
  }
  /* west */
  idnr = 4;
  i=0;
  printf("| western (x=0) boundary for bulk; idnr=%d\n", idnr);
  for (k=0;k<kmax;++k){ /* loop over all levels */
    for (j=0;j<jmax;++j){ /* loop over all j-values */
      ++n;
      nodes_of_element[0] = k*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = k*nodes_in_one_layer + (j+1)*(imax+1) + i;
      nodes_of_element[2] = (k+1)*nodes_in_one_layer + (j+1)*(imax+1) + i;
      nodes_of_element[3] = (k+1)*nodes_in_one_layer + j*(imax+1) + i;
      parent_element = k*elements_in_one_layer + j*imax + i; /* same numbering as in bulk */ 
      if (parent_element+1 > number_of_bulk_elements){
	printf("parent element %d > number of bulk elments %d\n",parent_element+1, number_of_bulk_elements);
	return;
      }
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n",
	      n, idnr, parent_element+1, nodes_of_element[0]+1, nodes_of_element[1]+1, nodes_of_element[2]+1, nodes_of_element[3]+1);
    }
  }
  /* north */ 
  idnr = 5;
  j=jmax;
  printf("| northern (y=max) boundary for bulk; idnr=%d\n", idnr);
  for (k=0;k<kmax;++k){ /* loop over all levels */
    for (i=0;i<imax;++i) { /* loop over all i-values */  
      nodes_of_element[0] = k*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = k*nodes_in_one_layer + j*(imax+1) + i+1;
      nodes_of_element[2] = (k+1)*nodes_in_one_layer + j*(imax+1) + i+1;
      nodes_of_element[3] = (k+1)*nodes_in_one_layer + j*(imax+1) + i;
      ++n;
      parent_element = k*elements_in_one_layer + (jmax-1)*imax + i; /* same numbering as in bulk */
      if (parent_element+1 > number_of_bulk_elements){
	printf("parent element %d > number of bulk elments %d\n",parent_element+1, number_of_bulk_elements);
	return;
      }
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n",
	      n, idnr, parent_element + 1, nodes_of_element[0] + 1, nodes_of_element[1] + 1, nodes_of_element[2] + 1, nodes_of_element[3] + 1);
    }
  }
  /* east */ 
  idnr = 6;
  i=imax;
  printf("| eastern (x=max) boundary for bulk; idnr=%d\n", idnr);
  for (k=0;k<kmax;++k){ /* loop over all levels */
    for (j=0;j<jmax;++j){ /* loop over all j-values */
      nodes_of_element[0] = k*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = k*nodes_in_one_layer + (j+1)*(imax+1) + i;
      nodes_of_element[2] = (k+1)*nodes_in_one_layer + (j+1)*(imax+1) + i;
      nodes_of_element[3] = (k+1)*nodes_in_one_layer + j*(imax+1) + i;
      ++n;
      parent_element = k*elements_in_one_layer + j*imax + i-1; /* same numbering as in bulk */
      if (parent_element >= number_of_bulk_elements){
	printf("parent element %d > number of bulk elments %d\n",parent_element, number_of_bulk_elements);
/* 	return; */
      }
      fprintf(ptFile,"%d %d %d 0 404 %d %d %d %d\n",
	      n, idnr, parent_element + 1, nodes_of_element[0] + 1, nodes_of_element[1] + 1, nodes_of_element[2] + 1, nodes_of_element[3] + 1);
    }
  }
  m=0;
  if (sidebulk){
    /* frame of sidebulk (i.e. free surface): */
    /* south */
    ++idnr;
    j=0;
    printf("| southern (y=0) boundary for sidebulk; idnr=%d\n", idnr);
    for (i=0;i<imax;++i) { /* loop over all i-values */
      ++m;
      nodes_of_element[0] = kmax*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = kmax*nodes_in_one_layer + j*(imax+1) + i+1;
      parent_element = number_of_bulk_elements + j*imax + i;
      fprintf(ptFile,"%d %d %d 0 202 %d %d\n",
	      n+m, idnr, parent_element + 1, nodes_of_element[0] + 1, nodes_of_element[1] + 1);
    }
    /* west */
    ++idnr;
    i=0;
    printf("| western (x=0) boundary for sidebulk; idnr=%d\n", idnr);
    for (j=0;j<jmax;++j){ /* loop over all j-values */
      ++m;
      nodes_of_element[0] = kmax*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = kmax*nodes_in_one_layer + (j+1)*(imax+1) + i;
      parent_element = number_of_bulk_elements + j*imax + i;
      fprintf(ptFile,"%d %d %d 0 202 %d %d\n",
	      n+m, idnr, parent_element + 1, nodes_of_element[0] + 1, nodes_of_element[1] + 1);
    }
    /* north */
    ++idnr;
    j=jmax;
    printf("| northern (y=max) boundary for sidebulk; idnr=%d\n", idnr); 
    for (i=0;i<imax;++i) { /* loop over all i-values */
      ++m;
      nodes_of_element[0] = (kmax-1)*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = (kmax-1)*nodes_in_one_layer + j*(imax+1) + i + 1;
      parent_element = number_of_bulk_elements + (j-1)*imax + i;
      fprintf(ptFile,"%d %d %d 0 202 %d %d\n",
	      n+m, idnr, parent_element + 1, nodes_of_element[0] + 1, nodes_of_element[1] + 1);
    }
    /* east */
    ++idnr;
    i=imax;
    printf("| western (x=max) boundary for sidebulk; idnr=%d\n", idnr);
    for (j=0;j<jmax;++j){ /* loop over all j-values */
      ++m;
      nodes_of_element[0] = (kmax-1)*nodes_in_one_layer + j*(imax+1) + i;
      nodes_of_element[1] = (kmax-1)*nodes_in_one_layer + j*(imax+1) + i+1;
      parent_element = number_of_bulk_elements + j*imax + i-1;
      fprintf(ptFile,"%d %d %d 0 202 %d %d\n",
	      n+m, idnr, parent_element + 1, nodes_of_element[0] + 1, nodes_of_element[1] + 1);
    }
  }
  fclose(ptFile);
  if (n+m !=  number_of_boundary_elements + sidebulk*2*(imax+jmax)){
    printf("\a %d=%d + %d  written boundary elements does not match %d= %d + %d  calculated elements in grid\n",
	   n+m, n,m, number_of_boundary_elements+sidebulk*2*(imax+jmax), number_of_boundary_elements, sidebulk*2*(imax+jmax));
    return;	   
  }else{
    printf("| %d number of boundary elements for bulk written.\n", n);
    printf("| %d number of boundary elements for sidebulk written.\n", m);
    printf("| succeeded in writting boundary element file %s.\n",filename);
    printf("---------------------------------------------------------------\n");
  }

  free(iced);free(delta_z);free(bottom);
  return;
}


/****************************************************************
   Get parameters from log file of SICOPOLIS run
*****************************************************************/
void STDCALLBULL FC_FUNC_(readlog_c,READLOG_C) (FC_CHAR_PTR(runname,l1), /*name of run*/
				    int    *imax, /* grid steps in xi-direction */
				    int    *jmax, /* grid steps in eta-direction */
				    int    *kcmax, /* grid steps in z-direction cold ice layer */
				    int    *ktmax,/* grid steps in z-direction temperate ice layer */
				    int    *krmax,/* grid steps in z-direction in bedrock */
				    float  *deform, /* parameter for vertical scaling */
				    float  *deltaX, /* horizontal grod spacing */
				    int    *gotit){
  /* variable declaration */
  int   gotparameter=0, i=0, inputint;
  char  filename[80], inputstring[100], parameter[6], chardummy[1];
  float inputfloat;
  FILE  *ptFile;

  
  /* compose filename */
  sprintf(filename,"%s%s", runname, ".log");

  /* print out little summary */
  printf("---------------------------------------------------------------\n");
  printf("|        Reading in parameters from log file\n");
  printf("---------------------------------------------------------------\n");

  if((ptFile = fopen(filename, "r")) == NULL){
    printf("Cannot open file %s\n", filename);
    *gotit = 0;
    return;
  }
  
  while((fgets(inputstring, LINESIZE - 1, ptFile) != NULL) && (gotparameter <7)){
    ++i;
    /* printf("%s\n", inputstring); */
    if (sscanf(inputstring, "%s %c %i3", parameter, chardummy, &inputint) == 3){
      /* printf("%d: %s, %c, %f\n", i, parameter, chardummy, inputfloat); */
      if (strcmp(parameter, "imax") == 0){
	*imax = inputint;
	gotparameter++;
      }
      else if (strcmp(parameter, "jmax") == 0){
	*jmax = inputint;
	gotparameter++;
      }
      else if (strcmp(parameter, "kcmax") == 0){
	*kcmax = inputint;
	gotparameter++;
      }
      else if (strcmp(parameter, "ktmax") == 0){
	*ktmax = inputint;
	gotparameter++;
      }
      else if (strcmp(parameter, "krmax") == 0){
	*ktmax = inputint;
	gotparameter++;
      }
    }
    if (sscanf(inputstring, "%s %c %f", parameter, chardummy, &inputfloat) == 3){
      if (strcmp(parameter, "a") == 0){
	*deform = inputfloat;
	gotparameter++;
      }
      else if (strcmp(parameter, "dx") == 0){
	*deltaX = inputfloat;
	gotparameter++;
      }
    }
  }
  printf("| imax/jmax/kcmax/ktmax=%d/%d/%d/%d\n",*imax, *jmax, *kcmax, *ktmax);
  printf("| deform = %f, dx = %f\n", *deform, *deltaX);
  printf("---------------------------------------------------------------\n");
  fclose(ptFile);
  *gotit=1;
  return;
}
