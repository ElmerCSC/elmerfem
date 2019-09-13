//
//  ExtrudeMesh - extrudes 2d mesh on a constant number of levels and
//                interpolates lower and upper surface from values
//                given by a DEM
//
//
//  Authors: thomas Zwinger, Thorsten Malm
//  Email: Thomas.Zwinger@csc.fi
//  Address: CSC - IT Center for Science Ltd.
//  Keilaranta 14
//  02101 Espoo, Finland
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor,
//  Boston, MA  02110-1301, USA.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/file.h>
#include <sys/param.h>

#define PERMS 0666
#define BUF_LEN 100000


#ifdef WIN32
#include <direct.h>
#define MKDIR(x) _mkdir(x)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(x) mkdir(x, S_IRWXU)
#endif


double f(double x,double GL, double dZ, double BL) {
  return dZ - (BL+dZ)*pow(x,GL)+BL*pow(x,GL+1);
}  /* f */

double f_prime(double x,double GL, double dZ, double BL) {
  return -(BL+dZ)*GL*pow(x,GL-1)+BL*(GL+1)*pow(x,GL);
}  /* f_prime */



double newton(double x_0, double tol, int max_iters, 
	      int* iters_p, int* converged_p,double GL, double dZ, double BL) {
   double x = x_0;
   double x_prev;
   int    iter = 0;

   do {
      iter++;
      x_prev = x;
      x = x_prev - f(x_prev,GL,dZ,BL)/f_prime(x_prev,GL,dZ,BL);
   } while ((fabs(x - x_prev) > tol && iter < max_iters));

   if (fabs(x - x_prev) <= tol)
      *converged_p = 1;
   else{
      *converged_p = 0;
      printf("not converged, %d\n",converged_p[0]);
   }
   *iters_p = iter;



   
   return x;
}  /* newton algorithm */





// Compute layer coordinate Z for geometric refinement at BL and uniform:
//-----------------------------------------------------
double refinedHeightBL (double S, double B, int levels, int level, int GL, double p, double r)
{

  
   int    iters;     /* Actual number of iterations  */
   int    converged; /* Whether iteration converged  */
   double z;


  int n = GL -1; 


  if(GL < 1){
    double dZ = (S - B)/((double) levels -1);
    z = B + level*dZ;
  }else{
    double BL = (S - B)*fabs(p);
    double dZ = (S - B - BL)/((double) levels - GL -1);
    if(r==0)  r=newton(10,1e-3,1000,&iters,&converged,GL,dZ,BL);
    //for r==0 r is calculated automatically, so that dZ=r*_last_boundary_layer_thickness_; but it's not always possible to find such an r > 1, pay attention to this!!!!

  
    double a;
    if(r>1) a = BL *(1-r)/(1-pow(r,n+1));
    else{
      a= BL/(GL-1);
      printf("no good solution found, since r is equal or smaller as 1: r= %f \n",r);
      printf("solution r %f dZ %f dZlast %f a %f\n",r,dZ,a*pow(r,GL-1),a);  
      printf("equation of the polynom for testing: poly  konst %f GL %i c1  %f c2 %f \n",dZ,GL,-BL-dZ,BL);  
    }


    int k = 0;
    
    if(p>0){
      z = B;
      for(k = 0; k < level; ++k) {
	if(k<GL){
	  z = z + a * pow(r, k);
	}
	else
	  z = z +  dZ;
      }
    }else{
      z = B;
      for(k = 0; k < level; ++k) {
	if(k<levels-GL-1){
	  z = z +  dZ;
	}
	else
	  z = z + a * pow(r, levels-k-2);
      }
    }
  }
  return z;
}

int file_exists(const char * filename){
  FILE * fileptr;
    if (fileptr = fopen(filename, "r")) 
    {
        fclose(fileptr);
        return 1;
    }
    return 0;
}


int checkFileEntries(char *filename, char *argv[], int *entriesPerLine){


  FILE *fileptr;
  char dummys[BUF_LEN], *charptr;
  int lines,entries;


  fileptr = fopen(filename, "r");
  if (fileptr==NULL){
    printf("(%s) Error opening file %s\n",argv[0],filename);
    return -1;
  }
    
  lines = 0;
  entries = 0;
  fgets(dummys, BUF_LEN,fileptr);
  charptr = strtok(dummys, "\t ");
  while(charptr != NULL) {
    entries++;
    charptr = strtok(NULL, "\t ");
  }
  while( fgets(dummys, BUF_LEN, fileptr) != NULL ){
    lines++;
  }
  fclose(fileptr);
  *entriesPerLine = entries;
  return ++lines;
}

int readDEM(double *inputfield, int *isnoval, int pointsInDEM, char *filename, char *argv[]){
  int i,validPointsInDEM;
  double noval = atof(argv[13]);
  FILE *fileptr;
  
  fileptr = fopen(filename, "r");
  if (fileptr==NULL){
    printf("(%s) Error opening file %s\n",argv[0],filename);
    return -1;
  }
  for (i=0,validPointsInDEM=0;i<pointsInDEM;++i){
    fscanf(fileptr,"%lf %lf %lf",&inputfield[3*validPointsInDEM],&inputfield[3*validPointsInDEM+1],&inputfield[3*validPointsInDEM+2]);
    if (inputfield[3*validPointsInDEM+2] == noval) 
      isnoval[validPointsInDEM] = 1;
      // inputfield[3*validPointsInDEM] = inputfield[3*validPointsInDEM+1] = inputfield[3*validPointsInDEM+2] = 0.0;
    else{
      validPointsInDEM++;
      isnoval[validPointsInDEM] = 0;
    }
  }
  fclose(fileptr);
  return validPointsInDEM;
}

int interpolatePoint(double X, double Y, double *B, double *S, double *bed, double *surf,  double *thick, int *pointsInDEM, double minrad, double wexp, int interpolationScheme, char *argv[], int *isnoval0, int *isnoval1, int *isnoval2){
  int i, j, k, rounds;
  double weightsum, weight, radius, localthickness;
  
  //printf("minrad=%e wexp=%e X=%e Y=%e \n",minrad,wexp, X, Y);

  // do we need an interpolated thickness?
  if (interpolationScheme > 11){
    rounds = 0;    
    do{
      localthickness = 0.0;
      for (i=0,weightsum=0.0;i<pointsInDEM[2];++i){
	radius = sqrt(pow((X-thick[3*i]),2.0) + pow((Y-thick[(3*i)+1]),2.0));
    
    
	//printf("radius=%e, bed=%e,%e\n",radius,bed[3*i], bed[3*i+1]); 

	if (radius <= minrad*(1.0 + (double)(rounds+1)/100.0)){
	  if (isnoval2[i] == 1) {
	    weight = 0.0;
	      }else{
	    weight = 1.0/(pow(radius,wexp) + 1.0E-09);
	  }
	  weight = (weight > 1.0E05) ? 1.0E05 : weight;
	  weightsum += weight;
	  localthickness += weight * thick[(3*i)+2];
	}
      }
      localthickness /= weightsum;
      if (localthickness != localthickness){
	rounds++;    
      }
    }while ((rounds < 100) && (*B != *B));

    if (rounds >= 99){
      fprintf(stderr,"(%s) (thickness) no interpolation point within range of double the cut-off value %e\n", argv[0], minrad);
	return EXIT_FAILURE;
    }else if (rounds > 0){
      fprintf(stderr,"(%s) (thickness) interpolation needed increasing cut-off value by %4e percent\n",argv[0],  rounds*1.0); 
    }
  }
  

  // bedrock interpolation from DEM
  if ((interpolationScheme == 11) || (interpolationScheme == 101)){
    rounds = 0;
    do{
      *B = 0.0;
      for (i=0,weightsum=0.0;i<pointsInDEM[0];++i){
	radius = sqrt(pow((X-bed[3*i]),2.0) + pow((Y-bed[(3*i)+1]),2.0));
    
    
	//printf("radius=%e, bed=%e,%e\n",radius,bed[3*i], bed[3*i+1]); 

	if (radius <= minrad*(1.0 + (double)(rounds+1)/100.0)){
	  if (isnoval0[i] == 1) {
	    weight = 0.0;
	      }else{
	    weight = 1.0/(pow(radius,wexp) + 1.0E-09);
	  }
	  weight = (weight > 1.0E05) ? 1.0E05 : weight;
	  weightsum += weight;
	  *B += weight * bed[(3*i)+2];
	}
      }
      *B /= weightsum;
      if (*B != *B){
	rounds++;    
      }
    }while ((rounds < 100) && (*B != *B));

    if (rounds >= 99){
      fprintf(stderr,"(%s) no interpolation point within range of double the cut-off value %e\n", argv[0], minrad);
      return EXIT_FAILURE;
    }else if (rounds > 0){
      fprintf(stderr,"(%s) interpolation needed increasing cut-off value by %4e percent\n",argv[0],  rounds*1.0); 
    }
  }
  
  //surface
  if ((interpolationScheme == 11) || (interpolationScheme == 110)){
    rounds = 0;
    do{
      *S=0.0;
      for (i=0,weightsum=0.0;i<pointsInDEM[1];++i){
	radius = sqrt(pow((X-surf[3*i]),2.0) + pow((Y-surf[(3*i)+1]),2.0));
	if (radius <= minrad*(1.0 + (double)(rounds+1)/100.0)){
	  if (isnoval1[i] == 1) {
	    weight = 0.0;
	      }else{
	    weight = 1.0/(pow(radius,wexp) + 1.0E-09);
	  }
	  weightsum += weight;
	  *S += weight * surf[(3*i)+2];
	}
      }
      if (weightsum <= 0.0){
	rounds++;    
      }else *S /= weightsum;    
    }while ((rounds < 100) && (weightsum <= 0.0));

    if (rounds >= 99){
      fprintf(stderr,"(%s)  (surface) no interpolation point within range of double the cut-off value %e\n", argv[0], minrad);
      return EXIT_FAILURE;
    }else if (rounds > 0){
      fprintf(stderr,"(%s) (surface) interpolation needed increasing cut-off value by %4e percent\n",argv[0],  rounds*1.0); 
    }
  }

  // get surface from bed and thickness
  if (interpolationScheme == 101) *S = *B + localthickness;
  // get bedrock from surface and thickness
  if (interpolationScheme == 110) *B = *S - localthickness;
    
//-----

  return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------------------------
// extrudes every single partition
//---------------------------------------------------------------------------------------------------
int extrudepartition( char *argv[], char *inputfilename, char *outputfilename, int partition, int partitions, int levels, double depth, int pointsinlevel, int elementsinlevel, int belementsinlevel, int interpolationScheme,  double *bed, double *surf, double *thick, int *pointsInDEM, double cutOffValue, double wexp, int *isnoval0, int *isnoval1, int *isnoval2)
{// extrudepartition
  int i,j,k,l,level, dummyint, nodesinpartition,  elementsinpartition, belementsinpartition, sharednodes, points, quads, lines, triangles, bricks, wedges, parent[2], elementtype, belementtype;
  char instring[BUF_LEN],dummys[BUF_LEN], numberstring[BUF_LEN], dummyc, *charptr;
  double dummydbl, dZ;
  // arrays to be allocated
  double *X,*Y,*S,*B;
  int *nodeinfo, *elementinfo, *belementinfo, *sharednodeinfo;
  char  **inputfile, **outputfile;
  FILE **outfids, **infids, *infofid;

  // allocate and initiate needed arrays for input
  //----------------------------------------------
  infids = (FILE **) malloc(5 * sizeof(FILE *));
  inputfile = (char **) malloc(5 * sizeof(char *)); 
  for (i=0;i<5;++i)
    inputfile[i] = (char *) malloc((MAXPATHLEN+1) * sizeof(char));

  for (i=0;i<5;++i){
    strcpy(inputfile[i],inputfilename);
  }
  strcat(inputfile[0],".header");//header
  strcat(inputfile[1],".nodes");//elements
  strcat(inputfile[2],".elements");//elements
  strcat(inputfile[3],".boundary");//boundary
  strcat(inputfile[4],".shared");//shared nodes
  for (i=0;i<5;++i){
    infids[i] = fopen(inputfile[i], "r");      
    if (infids[i] != NULL) continue;
    printf("(%s) Could not open file %s\n",argv[0],inputfile[i]);
    return EXIT_FAILURE;
  }
  // allocate and initiate needed arrays for output
  //-----------------------------------------------
  outfids = (FILE **) malloc(5 * sizeof(FILE *));
  outputfile = (char **) malloc(5 * sizeof(char *)); 
  for (i=0;i<5;++i)
    outputfile[i] = (char *) malloc((MAXPATHLEN+1) * sizeof(char));

  for (i=0;i<5;++i){
    strcpy(outputfile[i],outputfilename);
  }
  strcat(outputfile[0],".header");//header
  strcat(outputfile[1],".nodes");//elements
  strcat(outputfile[2],".elements");//elements
  strcat(outputfile[3],".boundary");//boundary
  strcat(outputfile[4],".shared");//shared nodes
  for (i=0;i<5;++i){
    outfids[i] = fopen(outputfile[i], "w");      
    if (outfids[i] != NULL) continue;
    printf("(%s) Could not open file %s\n",argv[0],outputfile[i]);
    return EXIT_FAILURE;
  }
  for (i=0;i<5;++i){
    printf("<-%s\n",inputfile[i]);
    printf("->%s\n",outputfile[i]);
  }
  // read in header in order to 
  //  inquire mesh sizes and open input files
  //-----------------------------------------
  fscanf(infids[0],"%i ",&nodesinpartition);
  fscanf(infids[0],"%i ",&elementsinpartition);
  fscanf(infids[0],"%i ",&belementsinpartition);

  fscanf(infids[0],"%i ",&j);
  points = lines = triangles = quads = 0;
  for (i=0;i<j;++i){
    fscanf(infids[0],"%i %i", &elementtype, &dummyint);
    if (elementtype == 101) points = dummyint;
    else if (elementtype == 202) lines = dummyint;
    else if (elementtype == 303) triangles = dummyint;
    else if (elementtype == 404) quads = dummyint;
    else{
      printf("(%s) element type %3i of input file %s not defined\n",argv[0],elementtype,inputfile[0]);
      return EXIT_FAILURE;
    }
  }
  fscanf(infids[0],"%i %i", &sharednodes, &dummyint);
  printf("(%s) Partition %i of %i input: %i nodes, %i elements, %i boundary elements, %i shared nodes\n",argv[0], partition+1, partitions, nodesinpartition,elementsinpartition,belementsinpartition,sharednodes);
  // now we also know the amount of wedges and bricks
  wedges = triangles*(levels - 1); // extruded triangles
  bricks = quads*(levels -1); // extruded quads
  triangles *= 2; // double, as free surface and bottom
  quads *= 2; // double, as free surface and bottom
  quads += (levels -1)*lines; // plus the extruded lines
  lines = (levels - 1)*points; // all extruded points

  // write partition header
  //printf("%i %i + %i  %i\n",points, wedges, bricks, quads+triangles+lines);
  fprintf(outfids[0],"%i %i %i\n",nodesinpartition*levels, wedges+bricks, quads+triangles+lines); // no. points, no- elements, no. boundary elements
  fprintf(outfids[0],"%i\n",(quads > 0) + (triangles > 0) + (bricks > 0) + (wedges > 0) + (lines > 0) +  (points > 0)); // no of different element types
  if (lines > 0) fprintf(outfids[0],"202 %i\n", lines);
  if (triangles > 0) fprintf(outfids[0],"303 %i\n", triangles);
  if (quads > 0)  fprintf(outfids[0],"404 %i\n", quads);
  if (wedges > 0) fprintf(outfids[0],"706 %i\n", wedges);
  if (bricks > 0) fprintf(outfids[0],"808 %i\n", bricks);
  fprintf(outfids[0],"%i 0 \n",sharednodes*levels); 
  printf("(%s) Partition %i of %i output: %i nodes, %i elements, %i boundary elements, %i shared nodes\n\n",argv[0], partition+1, partitions, nodesinpartition*levels, wedges+bricks, quads+triangles+lines, sharednodes*levels);

  
  X = (double *) malloc(nodesinpartition * sizeof(double));
  Y = (double *) malloc(nodesinpartition * sizeof(double));
  S = (double *) malloc(nodesinpartition * sizeof(double));
  B = (double *) malloc(nodesinpartition * sizeof(double));
  nodeinfo = (int *) malloc(nodesinpartition * sizeof(int));
  elementinfo = (int *) malloc(elementsinpartition * 8 * sizeof(int));
  belementinfo = (int *) malloc(belementsinpartition * 8 * sizeof(int));
  sharednodeinfo = (int *) malloc((partitions + 2) * sizeof(int));


  // load original mesh (footprint) nodes
  // and compute lower as well as upper
  //  Z-coordinate
  //-------------------------------------
  for (i=0; i<nodesinpartition;++i){
    fscanf(infids[1],"%i %i %le %le %le",&nodeinfo[i], &dummyint, &X[i], &Y[i], &B[i]);
    if (interpolationScheme > 0){
      if (interpolatePoint(X[i], Y[i], &B[i], &S[i], bed, surf, thick, pointsInDEM, cutOffValue, wexp, interpolationScheme, argv, isnoval0, isnoval1, isnoval2) == EXIT_FAILURE){
	fprintf(stderr,"(%s) Failed interpolating values for point %i %e %e\n",argv[0],i,X[i],Y[i]);
	return EXIT_FAILURE;
      }else{
	//printf("S,B= %e %e\n", B[i], S[i]);

	S[i] = (S[i] <= B[i] - depth) ? (B[i] + depth) :  S[i];
      }
    }else{
      S[i] = B[i] + depth;
      //printf("X(%i)= (%f, %f) Z=%f -> %f) \n ",nodeinfo[i],X[i],Y[i],B[i],S[i]);
    }
  }
/*   for (i=0; i<nodesinpartition;++i){ */
/*     fscanf(infids[1],"%i %i %le %le %le",&nodeinfo[i], &dummyint, &X[i], &Y[i], &B[i]); */
/*     if (interpolateDEM==1){ */
/*       printf("(%s) DEM interpolation yet not implemented\n",argv[0]); */
/*       return EXIT_FAILURE; */
/*     }else{ */
/*       S[i] = B[i] + depth; */
/*       printf("X(%i)= (%f, %f) Z=%f -> %f) \n ",nodeinfo[i],X[i],Y[i],B[i],S[i]); */
/*     } */
/*   } */



  // load original mesh (footprint) elements
  // elementinfo[8*k+0]/elementinfo[8*k+7] ... number/partition (if halo)
  // elementinfo[8*k+1] ... body
  // elementinfo[8*k+2] ... element type
  // elementinfo[8*k+3..6] ... element nodes 
  // elementinfo[8*k+7] ... partition owing the element(see first line)
  //---------------------------------------------------------------------
  for (k=0,quads=0,triangles=0; k<elementsinpartition;++k){
    fgets(instring, BUF_LEN-1,infids[2]);

    //printf("%i: %s\n",k,instring);
    sscanf(instring,"%s %i %i %i %i %i %i",&numberstring[0], &elementinfo[8*k+1],&elementinfo[8*k+2],
	   &elementinfo[8*k+3],&elementinfo[8*k+4],&elementinfo[8*k+5], &elementinfo[8*k+6]);
    sscanf(numberstring,"%i%c%i",&elementinfo[8*k+0], &dummyc, &elementinfo[8*k+7]);
    if (dummyc != '/')  elementinfo[8*k+7] = 0;
    if (elementinfo[8*k+2] == 404) quads++;
    else if (elementinfo[8*k+2] == 303){
      elementinfo[8*k+6] = -1;
      triangles++;
    } else {
      printf("(%s): element type %i not recognised for entry %i in partition %i\n", argv[0], elementinfo[8*k+2],k,partition);
    }
  }



  // load original mesh (footprint) boundary elements
  // belementinfo[8*k+0]/elementinfo[8*k+7] ... num  
  // belementinfo[8*k+1] ... boundary
  // belementinfo[8*k+2] ... parent 1
  // belementinfo[8*k+3] ... parent 2 
  // belementinfo[8*k+4] ... element type
  // belementinfo[8*k+5..6] ... element nodes
  // belementinfo[8*k+7] ... partition owing the element(see first line)
  //--------------------------------------------------------------------
  for (k=0; k<belementsinpartition;++k){
/*     fscanf(infids[3],"%s %i %i %i %i", &numberstring,&belementinfo[8*k+1],&belementinfo[8*k+2],&belementinfo[8*k+3],&belementinfo[8*k+4]);  */ 
    fgets(instring, BUF_LEN-1,infids[3]);
    //printf("%i: %s\n",k,instring);
    sscanf(instring,"%s %i %i %i %i %i %i",&numberstring[0], &belementinfo[8*k+1],&belementinfo[8*k+2],
	   &belementinfo[8*k+3],&belementinfo[8*k+4],&belementinfo[8*k+5], &belementinfo[8*k+6]);
    sscanf(numberstring,"%i%c%i",&belementinfo[8*k+0], &dummyc, &belementinfo[8*k+7]);
    //printf("%s: %i ,%c, %i\n",numberstring, belementinfo[8*k+0], dummyc, belementinfo[8*k+7]);
    if (dummyc != '/')  belementinfo[8*k+7] = 0;
    if (belementinfo[8*k+4] == 202){ 
      j=2;
    }else if (belementinfo[8*k+4] == 101){
      j=1;
      belementinfo[8*k+6] = -1;
    }else {
      printf("(%s): boundary element type %i not recognised for entry %i in partition %i\n",
	     argv[0], belementinfo[8*k+4],k,partition);
    }
/*     printf("R %i/%i %i %i %i %i %i %i\n",belementinfo[8*k+0],belementinfo[8*k+7], */
/* 	     belementinfo[8*k+1],belementinfo[8*k+2],belementinfo[8*k+3], */
/* 	     belementinfo[8*k+4],belementinfo[8*k+5],belementinfo[8*k+6]); */

  }



  // read and write shared nodes file
  //---------------------------------
  k=0;
  //printf("BUF_LEN = %i\n",BUF_LEN);
  while(fgets(instring, BUF_LEN-1, infids[4]) != NULL) { 
    ++k;
    //printf("%i: %s\n", k, instring);
    charptr = strtok(instring," ");
    i=0;
    while (charptr != NULL)
      {
	sharednodeinfo[i++] = atoi(charptr);
	//printf("%s->%i ", charptr,sharednodeinfo[i-1]);
	charptr = strtok(NULL," ");
      }
    //printf("\n");
    for(level=0;level<levels;++level){
      fprintf(outfids[4],"%i %i", pointsinlevel*level + sharednodeinfo[0],sharednodeinfo[1]);
      for (j=2;j<i;++j){
	fprintf(outfids[4],"  %i", sharednodeinfo[j]);
      }
      fprintf(outfids[4],"\n");
    }
  }

  // write extruded mesh
  //--------------------
  for (j=0,level=0;level<levels;++level){
    for (k=0;k<nodesinpartition;++k){
      l=100*k/pointsinlevel;
      if (l>j){
	j=l;
	//printf("done: %3i percent\n", j);
      }
            dZ = (S[k] - B[k])/((double) levels - 1);


      //printf("%i: S=%f, B=%f, dZ=%f\n",k,S[k],B[k],dZ);
 
            fprintf(outfids[1],"%i -1 %e %e %e\n", pointsinlevel*level + nodeinfo[k], X[k], Y[k], B[k]+ ((double) level)*dZ ); //

    }
  }

  // write elements
  // elementinfo[8*k+0]/elementinfo[8*k+7] ... number/partition (if halo)
  // elementinfo[8*k+1] ... body
  // elementinfo[8*k+2] ... element type
  // elementinfo[8*k+3..6] ... element nodes 
  // elementinfo[8*k+7] ... partition owing the element(see first line)
  //--------------------------------------------------------------------
  for (k=0;k<elementsinpartition;++k){
    if (elementinfo[8*k+2] == 303) {elementtype = 706;l=3;}
    else if (elementinfo[8*k+2] == 404) {elementtype = 808;l=4;}
    for (level=0;level<levels-1;++level){
      fprintf(outfids[2],"%i", elementsinlevel*level + elementinfo[8*k+0]); // number
      if (elementinfo[8*k+7] > 0)
	fprintf(outfids[2],"/%i ", elementinfo[8*k+7]); // if halo, then owner partition
      fprintf(outfids[2]," %i %i", elementinfo[8*k+1], elementtype); //  body, elementinfo[8*k+ 2]
      for (j=0;j<l;++j) fprintf(outfids[2]," %i", pointsinlevel*level + elementinfo[8*k+ 3 + j]); //lower level points
      for (j=0;j<l;++j) fprintf(outfids[2]," %i", pointsinlevel*(level + 1) + elementinfo[8*k+ 3 + j]);//upper level points
      fprintf(outfids[2],"\n");
    }
  }

  // write boundary elements
  // belementinfo[8*k+0]/belementinfo[8*k+7] ... num 
  // belementinfo[8*k+1] ... boundary (off-setted by 2 for original ones)
  // belementinfo[8*k+2] ... parent 1
  // belementinfo[8*k+3] ... parent 2 
  // belementinfo[8*k+4] ... element type
  // belementinfo[8*k+5..6] ... element nodes
  // belementinfo[8*k+7] ... partition owing the element(see first line)
  //--------------------------------------------------------------------
  //  bottom and surface (new, bcids=1,2)
  parent[1]=0; 
  for (i=0,lines=0;i<2;++i){ // i==0 ... bottom, i==1 ... top 
    for (k=0;k<elementsinpartition;++k){ 
      parent[0] =  elementsinlevel*(levels-2)*i + elementinfo[8*k+0];
      fprintf(outfids[3],"%i",  i*elementsinpartition + k +1); // element number	      
      if (elementinfo[8*k+7] > 0)
	fprintf(outfids[3],"/%i", elementinfo[8*k+7]); // if halo, then owner partition
      fprintf(outfids[3]," %1i %i %i %3i", i+1, parent[0], parent[1], elementinfo[8*k+2]); // bcid, parents, type
      if (elementinfo[8*k+2] == 303) {l=3;}
      else if (elementinfo[8*k+2] == 404) {l=4;}
      for(j=0;j<l;++j){
	fprintf(outfids[3]," %i", pointsinlevel*(levels-1)*i + elementinfo[8*k+ 3 + j]);
      }
      fprintf(outfids[3],"\n");
    }
  }

  // already existing boundaries (extruded sides, bcids = original + 2)
  for (k=0;k<belementsinpartition;++k){ //belementsinpartition
    for (level=0;level<levels-1;++level){
      for (i=0;i<2;++i)
	parent[i] = (belementinfo[8*k+2+i] > 0) ? (level*elementsinlevel + belementinfo[8*k+2+i]) : 0;
      if (belementinfo[8*k+4] == 202){
	belementtype = 404;
	++quads;
	j=2;
      }else if(belementinfo[8*k+4] == 101){
	belementtype = 202;
	++lines;
	j=1;
      }
      else{
	printf("(%s) boundary element type %3i of input file %s (line %i) not defined\n",argv[0],belementinfo[8*k+4],inputfile[3],k);
	return EXIT_FAILURE;
      }
/*       printf("W %i/%i %i %i %i %i %i %i\n",belementinfo[8*k+0],belementinfo[8*k+7], */
/* 	     belementinfo[8*k+1]+2,belementinfo[8*k+2],belementinfo[8*k+3], */
/* 	     belementinfo[8*k+4],belementinfo[8*k+5],belementinfo[8*k+6]); */
      fprintf(outfids[3],"%i",  2*elementsinpartition + belementsinpartition*level + belementinfo[8*k+0]); // number
      if (belementinfo[8*k+7] > 0)
	fprintf(outfids[3],"/%i",belementinfo[8*k+7]); //partition (if halo)
      fprintf(outfids[3]," %i %i %i %3i",  
	      belementinfo[8*k+1]+2, parent[0], parent[1], belementtype); // bcid, parents, type
      for (l=0;l<2;++l){
	if (l==0) {
	  for (i=0;i<j;++i){
	    fprintf(outfids[3]," %i", pointsinlevel*level + belementinfo[8*k+5+i]);
	    //printf("low: %i + %i\n", pointsinlevel*(level+l) , belementinfo[8*k+5+i]);
	  }
	}else{
	  for (i=j-1;i>-1;--i){
	    fprintf(outfids[3]," %i", pointsinlevel*(level+1) + belementinfo[8*k+5+i]);
	    //printf(" %i + %i\n", pointsinlevel*(level+l) , belementinfo[8*k+5+i]);
	  }
	}
      }
      fprintf(outfids[3],"\n");
    }
  }


  // close files
  //------------
  for (i=0;i<5;++i){
    fclose(infids[i]);
    fclose(outfids[i]);
  }
  // free memory
  //------------
  free(X);
  free(Y);
  free(S);
  free(B);
  free(nodeinfo);
  free(elementinfo);
  free(belementinfo);
  free(sharednodeinfo);

  free(infids); 
  for (i=0;i<5;++i)
    free(inputfile[i]);
  free(inputfile); 
  free(outfids); 
  for (i=0;i<5;++i)
    free(outputfile[i]);
  free(outputfile); 

  return EXIT_SUCCESS;
}
  



//---------------------------------------------------------------------------------------------------
// extrudes serial mesh
//---------------------------------------------------------------------------------------------------
int extrudeserial( char *argv[], char *inputfilename, char *outputfilename, int levels, double depth, int pointsinlevel, int elementsinlevel, int belementsinlevel, int interpolationScheme,  double *bed, double *surf, double *thick, int *pointsInDEM, double cutOffValue, double wexp, int GL, double percentage, double ratio, int baseline, int corrbed, int *isnoval0, int *isnoval1, int *isnoval2 )
{//extrudeserial
  int i,j,k,l,level, dummyint, nodesinpartition,  elementsinpartition, belementsinpartition, points, quads, lines, triangles, bricks, wedges, parent[2], elementtype, belementtype, bcoffset, maxorigBCno;
  char instring[BUF_LEN],dummys[BUF_LEN], numberstring[BUF_LEN],  *charptr;
  double dummydbl, dZ;
  // arrays to be allocated
  double *X,*Y,*S,*B, Z;
  int *nodeinfo, *elementinfo, *belementinfo;
  char  **inputfile, **outputfile;
  FILE **outfids, **infids, *infofid;

  // allocate and initiate needed arrays for input
  //----------------------------------------------
  infids = (FILE **) malloc(5 * sizeof(FILE *));
  inputfile = (char **) malloc(5 * sizeof(char *)); 
  for (i=0;i<5;++i)
    inputfile[i] = (char *) malloc((MAXPATHLEN+1) * sizeof(char));

  for (i=0;i<5;++i){
    strcpy(inputfile[i],inputfilename);
  }
  strcat(inputfile[0],"/mesh.header");//header
  strcat(inputfile[1],"/mesh.nodes");//elements
  strcat(inputfile[2],"/mesh.elements");//elements
  strcat(inputfile[3],"/mesh.boundary");//boundary
  for (i=0;i<5;++i){
    infids[i] = fopen(inputfile[i], "r");      
    if (infids[i] != NULL) continue;
    printf("(%s) Could not open file %s\n",argv[0],inputfile[i]);
    return EXIT_FAILURE;
  }
  // allocate and initiate needed arrays for output
  //-----------------------------------------------
  outfids = (FILE **) malloc(4 * sizeof(FILE *));
  outputfile = (char **) malloc(4 * sizeof(char *)); 
  for (i=0;i<4;++i)
    outputfile[i] = (char *) malloc((MAXPATHLEN+1) * sizeof(char));

  for (i=0;i<4;++i){
    strcpy(outputfile[i],outputfilename);
  }
  strcat(outputfile[0],"/mesh.header");//header
  strcat(outputfile[1],"/mesh.nodes");//elements
  strcat(outputfile[2],"/mesh.elements");//elements
  strcat(outputfile[3],"/mesh.boundary");//boundary
  for (i=0;i<4;++i){
    outfids[i] = fopen(outputfile[i], "w");      
    if (outfids[i] != NULL) continue;
    printf("(%s) Could not open file %s\n",argv[0],outputfile[i]);
    return EXIT_FAILURE;
  }
  for (i=0;i<4;++i){
    printf("<-%s\n",inputfile[i]);
    printf("->%s\n",outputfile[i]);
  }
  // read in header in order to 
  //  inquire mesh sizes and open input files
  //-----------------------------------------
  fscanf(infids[0],"%i ",&nodesinpartition);
  fscanf(infids[0],"%i ",&elementsinpartition);
  fscanf(infids[0],"%i ",&belementsinpartition);

  fscanf(infids[0],"%i ",&j);
  points = lines = triangles = quads = 0;
  for (i=0;i<j;++i){
    fscanf(infids[0],"%i %i", &elementtype, &dummyint);
    if (elementtype == 101) points = dummyint;
    else if (elementtype == 202) lines = dummyint;
    else if (elementtype == 303) triangles = dummyint;
    else if (elementtype == 404) quads = dummyint;
    else{
      printf("(%s) element type %3i of input file %s not defined\n",argv[0],elementtype,inputfile[0]);
      return EXIT_FAILURE;
    }
  }

  printf("(%s) Serial mesh input: %i nodes, %i elements, %i boundary elements\n",argv[0], nodesinpartition,elementsinpartition,belementsinpartition);
  // now we also know the amount of wedges and bricks
  wedges = triangles*(levels - 1); // extruded triangles
  bricks = quads*(levels -1); // extruded quads
  triangles *= 2; // double, as free surface and bottom
  quads *= 2; // double, as free surface and bottom
  quads += (levels -1)*lines; // plus the extruded lines
  lines = (levels - 1)*points + baseline*lines; // all extruded points and included baselines

  // write partition header
  //printf("%i %i + %i  %i\n",points, wedges, bricks, quads+triangles+lines);
  fprintf(outfids[0],"%i %i %i\n",nodesinpartition*levels, wedges+bricks, quads+triangles+lines); // no. points, no- elements, no. boundary elements
  fprintf(outfids[0],"%i\n",(quads > 0) + (triangles > 0) + (bricks > 0) + (wedges > 0) + (lines > 0)); // no of different element types
  if (lines > 0) fprintf(outfids[0],"202 %i\n", lines);
  if (triangles > 0) fprintf(outfids[0],"303 %i\n", triangles);
  if (quads > 0)  fprintf(outfids[0],"404 %i\n", quads);
  if (wedges > 0) fprintf(outfids[0],"706 %i\n", wedges);
  if (bricks > 0) fprintf(outfids[0],"808 %i\n", bricks);
  printf("(%s) Serial mesh output: %i nodes, %i elements, %i boundary elements\n\n",argv[0],  nodesinpartition*levels, wedges+bricks, quads+triangles+lines);

  
  X = (double *) malloc(nodesinpartition * sizeof(double));
  Y = (double *) malloc(nodesinpartition * sizeof(double));
  S = (double *) malloc(nodesinpartition * sizeof(double));
  B = (double *) malloc(nodesinpartition * sizeof(double));
  nodeinfo = (int *) malloc(nodesinpartition * sizeof(int));
  elementinfo = (int *) malloc(elementsinpartition * 7 * sizeof(int));
  belementinfo = (int *) malloc(belementsinpartition * 7 * sizeof(int));

 
  // load original mesh (footprint) nodes
  // and compute lower as well as upper
  //  Z-coordinate
  //-------------------------------------
  for (i=0; i<nodesinpartition;++i){
    fscanf(infids[1],"%i %i %le %le %le",&nodeinfo[i], &dummyint, &X[i], &Y[i], &B[i]);
    if (interpolationScheme > 0){
      if (interpolatePoint(X[i], Y[i], &B[i], &S[i], bed, surf, thick, pointsInDEM, cutOffValue, wexp, interpolationScheme, argv, isnoval0, isnoval1, isnoval2) == EXIT_FAILURE){
	fprintf(stderr,"(%s) Failed interpolating values for point %i %e %e\n",argv[0],i,X[i],Y[i]);
	return EXIT_FAILURE;
      }else{
	//printf("S,B= %e %e\n", B[i], S[i]);
	if (S[i] -  B[i] < depth){
	  printf("corrected  S[%i] %e %e %e\n", i, X[i], Y[i], S[i]);
	  if (corrbed == 1){
	    B[i] = S[i]  - depth;
	  }else{
	    S[i]  = B[i] + depth;
	  }
	  printf(" ->  %e %e %e\n", X[i], Y[i], S[i]);
	}
	//S[i] = (S[i] <= B[i] - depth) ? (B[i] + depth) :  S[i];
      }
    }else{ // constant extrusion offset
      S[i] = B[i] + depth;
      //printf("X(%i)= (%f, %f) Z=%f -> %f) \n ",nodeinfo[i],X[i],Y[i],B[i],S[i]);
    }
  }


  // load original mesh (footprint) elements
  // elementinfo[7*k+0] ... number
  // elementinfo[7*k+1] ... body
  // elementinfo[7*k+2] ... element type
  // elementinfo[7*k+3..6] ... element nodes 
  //---------------------------------------------------------------------
  for (k=0,quads=0,triangles=0; k<elementsinpartition;++k){
    fgets(instring, BUF_LEN-1,infids[2]);

    //printf("%i: %s\n",k,instring);
    sscanf(instring,"%i %i %i %i %i %i %i",&elementinfo[7*k+0], &elementinfo[7*k+1],&elementinfo[7*k+2],
	   &elementinfo[7*k+3],&elementinfo[7*k+4],&elementinfo[7*k+5], &elementinfo[7*k+6]);
    if (elementinfo[7*k+2] == 404) quads++;
    else if (elementinfo[7*k+2] == 303){
      elementinfo[7*k+6] = -1;
      triangles++;
    } else {
      printf("(%s): element type %i not recognised for entry %i\n", argv[0], elementinfo[7*k+2],k);
    }
  }



  // load original mesh (footprint) boundary elements
  // belementinfo[7*k+0] ... number 
  // belementinfo[7*k+1] ... boundary
  // belementinfo[7*k+2] ... parent 1
  // belementinfo[7*k+3] ... parent 2 
  // belementinfo[7*k+4] ... element type
  // belementinfo[7*k+5..6] ... element nodes
  //--------------------------------------------------------------------
  maxorigBCno = 0;
  for (k=0; k<belementsinpartition;++k){
/*     fscanf(infids[3],"%s %i %i %i %i", &numberstring,&belementinfo[7*k+1],&belementinfo[7*k+2],&belementinfo[7*k+3],&belementinfo[7*k+4]);  */ 
    fgets(instring, BUF_LEN-1,infids[3]);
    //printf("%i: %s\n",k,instring);
    sscanf(instring,"%i %i %i %i %i %i %i",&belementinfo[7*k+0], &belementinfo[7*k+1],&belementinfo[7*k+2],
	   &belementinfo[7*k+3],&belementinfo[7*k+4],&belementinfo[7*k+5], &belementinfo[7*k+6]);
    if (belementinfo[7*k+4] == 202){ 
      j=2;
    }else if (belementinfo[7*k+4] == 101){
      j=1;
      belementinfo[7*k+6] = -1;
    }else {
      printf("(%s): boundary element type %i not recognised for entry %i\n",
	     argv[0], belementinfo[7*k+4],k);
    }
/*     printf("R %i/%i %i %i %i %i %i %i\n",belementinfo[7*k+0],belementinfo[7*k+7], */
/* 	     belementinfo[7*k+1],belementinfo[7*k+2],belementinfo[7*k+3], */
/* 	     belementinfo[7*k+4],belementinfo[7*k+5],belementinfo[7*k+6]); */
    maxorigBCno = (belementinfo[7*k+1] > maxorigBCno) ? belementinfo[7*k+1]:maxorigBCno;
  }


  // write extruded mesh
  //--------------------
  for (j=0,level=0;level<levels;++level){
    for (k=0;k<nodesinpartition;++k){
      l=100*k/pointsinlevel;
      //if (l>j){
      //j=l;
	//printf("done: %3i percent\n", j);
      //}
      //      dZ = (S[k] - B[k])/((double) levels - 1);
      Z = refinedHeightBL(S[k], B[k], levels, level, GL, percentage, ratio);
      //if (percentage < 0) Z=S[k] + Z + B[k];

      //printf("%i: S=%f, B=%f, dZ=%f\n",k,S[k],B[k],dZ);
 
      fprintf(outfids[1],"%i -1 %e %e %e\n", pointsinlevel*level + nodeinfo[k], X[k], Y[k], Z ); //

 //    fprintf(outfids[1],"%i -1 %e %e %e\n", pointsinlevel*level + nodeinfo[k], X[k], Y[k], B[k]+ ((double) level)*dZ ); //
    }
  }

  // write elements
  // elementinfo[7*k+0] ... number
  // elementinfo[7*k+1] ... body
  // elementinfo[7*k+2] ... element type
  // elementinfo[7*k+3..6] ... element nodes 
  //--------------------------------------------------------------------
  for (k=0,bcoffset=1;k<elementsinpartition;++k){
    if (elementinfo[7*k+2] == 303) {elementtype = 706;l=3;}
    else if (elementinfo[7*k+2] == 404) {elementtype = 808;l=4;}
    for (level=0;level<levels-1;++level){
      fprintf(outfids[2],"%i", elementsinlevel*level + elementinfo[7*k+0]); // number
      fprintf(outfids[2]," %i %i", elementinfo[7*k+1], elementtype); //  body, elementinfo[7*k+ 2]
      bcoffset = (bcoffset > elementinfo[7*k+1]) ? bcoffset : elementinfo[7*k+1];
      for (j=0;j<l;++j) fprintf(outfids[2]," %i", pointsinlevel*level + elementinfo[7*k+ 3 + j]); //lower level points
      for (j=0;j<l;++j) fprintf(outfids[2]," %i", pointsinlevel*(level + 1) + elementinfo[7*k+ 3 + j]);//upper level points
      fprintf(outfids[2],"\n");
    }
  }

  printf("bodies: %i\n", bcoffset);

  // write boundary elements
  // belementinfo[7*k+0] ... number
  // belementinfo[7*k+1] ... boundary (off-setted by 2 for original ones)
  // belementinfo[7*k+2] ... parent 1
  // belementinfo[7*k+3] ... parent 2 
  // belementinfo[7*k+4] ... element type
  // belementinfo[7*k+5..6] ... element nodes
  //--------------------------------------------------------------------
  //  bottom and surface (new, bcids=1..bcoffset,bcoffset+1..2*bcoffset)
  parent[1]=0; 
  for (i=0,lines=0;i<2;++i){ // i==0 ... bottom, i==1 ... top 
    for (k=0;k<elementsinpartition;++k){ 
      parent[0] =  elementsinlevel*(levels-2)*i + elementinfo[7*k+0];
      fprintf(outfids[3],"%i",  i*elementsinpartition + k +1); // element number
      fprintf(outfids[3]," %i %i %i %3i",  elementinfo[7*k+1]+i*bcoffset, parent[0], parent[1], elementinfo[7*k+2]); // bcid, parents, type
      if (elementinfo[7*k+2] == 303) {l=3;}
      else if (elementinfo[7*k+2] == 404) {l=4;}
      for(j=0;j<l;++j){
	fprintf(outfids[3]," %i", pointsinlevel*(levels-1)*i + elementinfo[7*k+ 3 + j]);
      }
      fprintf(outfids[3],"\n");
    }
  }

  // already existing boundaries (extruded sides, bcids = original + 2)
  for (k=0;k<belementsinpartition;++k){ //belementsinpartition
    for (level=0;level<levels-1;++level){
      for (i=0;i<2;++i)
	parent[i] = (belementinfo[7*k+2+i] > 0) ? (level*elementsinlevel + belementinfo[7*k+2+i]) : 0;
      if (belementinfo[7*k+4] == 202){
	belementtype = 404;
	++quads;
	j=2;
      }else if(belementinfo[7*k+4] == 101){
	belementtype = 202;
	++lines;
	j=1;
      }
      else{
	printf("(%s) boundary element type %3i of input file %s (line %i) not defined\n",argv[0],belementinfo[7*k+4],inputfile[3],k);
	return EXIT_FAILURE;
      }
      fprintf(outfids[3],"%i",  2*elementsinpartition + belementsinpartition*level + belementinfo[7*k+0]); // number
      fprintf(outfids[3]," %i %i %i %3i",  
	      belementinfo[7*k+1]+2*bcoffset, parent[0], parent[1], belementtype); // bcid, parents, type
      for (l=0;l<2;++l){
	if (l==0) {
	  for (i=0;i<j;++i){
	    fprintf(outfids[3]," %i", pointsinlevel*level + belementinfo[7*k+5+i]);
	    //printf("low: %i + %i\n", pointsinlevel*(level+l) , belementinfo[7*k+5+i]);
	  }
	}else{
	  for (i=j-1;i>-1;--i){
	    fprintf(outfids[3]," %i", pointsinlevel*(level+1) + belementinfo[7*k+5+i]);
	    //printf(" %i + %i\n", pointsinlevel*(level+l) , belementinfo[7*k+5+i]);
	  }
	}
      }
      fprintf(outfids[3],"\n");
    }
  }

  // include the frame of the 2d footprint (e.g., in order to place BC for DIM-1 problem)
  if (baseline) {
    for (k=0;k<belementsinpartition;++k){ //belementsinpartition
      parent[0] = parent[1] = 0;
      //	for (i=0;i<2;++i)
      //	  parent[i] = (belementinfo[7*k+2+i] > 0) ? (belementinfo[7*k+2+i]) : 0;
	belementtype = belementinfo[7*k+4];
	if (belementtype == 202)  {++lines; j=2;}
	else if (belementtype == 101)  {++points; j=2;}
	fprintf(outfids[3],"%i",  2*elementsinpartition + belementsinpartition*(levels-1) + belementinfo[7*k+0]); //number
	fprintf(outfids[3]," %i %i %i %3i",  
		belementinfo[7*k+1] + (2*bcoffset) + maxorigBCno, parent[0], parent[1], belementtype); // bcid, parents, type
	for (i=0;i<j;++i){
	    fprintf(outfids[3]," %i", belementinfo[7*k+5+i]);
	}
	fprintf(outfids[3],"\n");
    }
  }


  // close files
  //------------
  for (i=0;i<4;++i){
    fclose(infids[i]);
    fclose(outfids[i]);
  }
  // free memory
  //------------
  free(X);
  free(Y);
  free(S);
  free(B);
  free(nodeinfo);
  free(elementinfo);
  free(belementinfo);

  free(infids); 
  for (i=0;i<4;++i)
    free(inputfile[i]);
  free(inputfile); 
  free(outfids); 
  for (i=0;i<4;++i)
    free(outputfile[i]);
  free(outputfile); 

  return EXIT_SUCCESS;
}



int main(int argc, char *argv[])
{
  int i,j,k,l,levels,partitions,partition,pointsinlevel, elementsinlevel, belementsinlevel, stat, interpolateDEM, fileExists[3], interpolationScheme, noelementtypes,triangles,quads,lines,points, pointsInDEM[3],dummyint,GL, baseline, *isnoval0, *isnoval1, *isnoval2, corrbed=0;
  char inputdirectoryname[MAXPATHLEN+1],outputdirectoryname[MAXPATHLEN+1],directoryname[MAXPATHLEN+1], cpartition[11], inputfilename[MAXPATHLEN+1], outputfilename[MAXPATHLEN+1], filename[MAXPATHLEN+1], *charptr;
  double depth, cutOffValue, *bed, *surf, *thick, wexp, noval, percentage, ratio;
  FILE *infofid, *headerfid;

  // failure/usage message
  //----------------------
  if((argc != 10) && ((argc != 14)) && (argc != 15)){
    printf("%s usage:\a\n\n",argv[0]);
    printf("%s inputdir outputdir levels extrudedepth N baseline GL percentage ratio DEM cutoff wexp noval\n\n",argv[0]);
    printf("     inputdir ... directory containing the 2D footrpint to be extruded\n");
    printf("    outputdir ... directory (will be created if not existing) containing\n");
    printf("                  the extruded mesh\n");
    printf("                  WARNING: Files in outputdir will be overwritten!\n");
    printf("       levels ... levels in direction of extrushion (>2)\n");
    printf(" extrudedepth ... depth in direction of extrushion in unit-length of\n");    
    printf("                  input mesh\n");
    printf("           N  ... partitions\n");
    printf("                  inputdir should contain an already partitioned mesh)\n");
    printf("                  give N=1 for serial mesh\n");
    printf("   baseline  ...  1/0\n");
    printf("                  1 if the existing side boundary for footprint (point/line in 2d/3d)\n");
    printf("                  shall be included, else excluded\n\n");
    printf("         GL  ... number of layers inside boundary layer\n");
    printf("                 has to be smaller than levels\n");
    printf("                 give 0, if no boundary layer wanted\n");
    printf("  percentage ... percentage height of boundary layer\n");
    printf("                 with respect to the local height (as 0..1) \n"); 
    printf("                 if positive BL at bottom, if negative at the surface \n"); 
    printf("       ratio ... ratio between two adjacent layers in boundary layer(>1) \n"); 
    printf("                  has to be >1, if it's 0 automatic calculation, but does not work always \n"); 

    printf("         DEM  ... (optional) directory of digital elevation model\n");
    printf("                  for interpolation of bed/surface (see below)\n");
    printf("      cutoff  ... (optional, but mandatory if previous was given)\n");
    printf("                  cutoff radius for interpolation\n");
    printf("      wexp  ... (optional, but mandatory if previous was given)\n");
    printf("                  exponent of weight for interpolation w=1/r^wexp\n");
    printf("      noval  ... (optional, but mandatory if previous was given)\n");
    printf("                  value indicating void/unvalid entry in DEM (e.g. -999.99)\n");
    printf("     corrbed ... (optional) value: 1 or any other number\n");
    printf("                 corrections induced by minimum flowdpeth are by default applied\n");
    printf("                 correctiong the side of the free surface, only of value 1 is defined\n");
    printf("                 the surface is kept constant and the bedrock is adjusted\n"); 
    printf("outputdir contains:\n");
    printf(" mesh.{nodes,elements,boundary} ... mesh files (if N==1)\n");
    printf(" partitioning.N/mesh.{nodes,elements,boundary} ... mesh files (if N>1)\n");
    printf(" info.dat ... file containing information on the extruded mesh\n"); 
    printf("DEM (if chosen) should provide (two of these three):\n");
    printf(" DEM/bed.xyz  (hard-coded name!)\n");
    printf(" DEM/surf.xyz (hard-coded name!)\n");
    printf(" DEM/thick.xyz (hard-coded name!)\n");
    printf(" mind, that those filenames are hardcoded!\n");
    printf(" the xyz-files can contain a set of irregullary distributed (x,y,z)\n");
    printf(" coordinates, between those the mesh will be interpolated.\n");
    printf(" In case of a DEM being provided, the earlier given\n");
    printf(" parameter \"extrudedepth\" is taken as the minimum flow depth\n\n");
   printf(" Boundary layer implemented only in serial\n\n");
    return EXIT_FAILURE;
  } 

  // inquire parameters
  //-------------------
  levels = atoi(argv[3]);
  depth =  atof(argv[4]);
  partitions = atoi(argv[5]);
  baseline = atoi(argv[6]);
  if (baseline != 1) baseline = 0;
  GL = atoi(argv[7]);
  if (GL >= levels) {
    fprintf(stderr, "GL = %i >= levels = %i !\n", GL, levels);
    return EXIT_FAILURE;
  }
  percentage =  atof(argv[8]);
  if ((percentage < -0.99) || (percentage > 0.99)) {
    fprintf(stderr, "0.01 < percentage = %f10.4 < 0.99 !\n", percentage);
    return EXIT_FAILURE;
  }
  ratio = atof(argv[9]); 
     if (ratio < 1.0 && ratio != 0) {
       fprintf(stderr, "ratio < 1.0 = %f10.4 < 1.0 or equal 0!\n", ratio);
     return EXIT_FAILURE;
   }

  printf("(%s) executing with following inputs ...\n",argv[0]);
  printf("    ... input directory: %s \n",argv[1]);
  printf("    ... output directory: %s\n",argv[2]);
  printf("    ... levels: %i\n",levels);
  printf("    ... depth: %10.4f\n",depth);
  printf("    ... partitions: %i\n",partitions);
  printf("    ... include baseline(1=yes,0=no): %i\n",baseline);
  printf("    ... GL: %i\n",GL);
  printf("    ... percentage: %10.4f\n",percentage);
  printf("    ... ratio: %10.4f\n",ratio);
  if (partitions < 1){
    printf("(%s) Wrong number of paritions %i\n",argv[0],partitions);
    printf("(%s) Assuming serial mesh (N=1)\n",argv[0]);
    partitions=1;
  }
  interpolateDEM = 0;
  if ((argc == 14) || (argc == 15)){// we have a DEM to interpolate
    for (i=0;i<argc;++i) {
      printf("argv[%i] = %s\n",i, argv[i]);
    }

    cutOffValue =  atof(argv[11]);
    wexp = atof(argv[12]);
    noval = atof(argv[13]);
    printf("    ... DEM directory: %s\n", argv[10]);
    printf("    ... cutoff radius: %e\n", cutOffValue);
    printf("    ... exponent m of weight (1/r^m): %e\n", wexp);
    printf("    ... noval: %10.4f\n",noval);
    if (corrbed == 1){
      printf("    ... depth corrections of bedrock\n");
    }else{
      printf("    ... depth corrections of surface (default)\n");
    }

    int entriesPerLine;
 
    // the bedrock DEM
    printf("(%s) Checking for input files of DEM in directory %s\n", argv[0], argv[7]);
    strcpy(filename,argv[10]); strcat(filename,"/bed.xyz");
    if (file_exists(filename)){
      pointsInDEM[0] = checkFileEntries(filename, argv, &entriesPerLine);
      if ((pointsInDEM[0] < 1) || (entriesPerLine < 3)){
	fprintf(stderr,"(%s) Wrong file entries (%i,%i) found in file %s\n", argv[0], pointsInDEM[0], entriesPerLine, filename);
	return EXIT_FAILURE;
      }else{
	printf("(%s) %i file entries found in file %s\n", argv[0], pointsInDEM[0], filename);
	bed = (double *) malloc(pointsInDEM[0] * 3 * sizeof(double));
	//      for (i=0;i<pointsInDEM[0];++i){ bed[3*i]=(double) 3*i, bed[(3*i)+1]=2.0, bed[(3*i)+2]=3.0; printf("-> %f %f %f\n",bed[3*i], bed[(3*i)+1], bed[(3*i)+2]);}
	isnoval0 = (int *) malloc( pointsInDEM[0] * sizeof(int));
	dummyint = readDEM(bed, isnoval0, pointsInDEM[0], filename, argv);
	if (dummyint > 1){
	  printf("(%s) %i out of %i valid points in DEM %s\n",argv[0],dummyint,pointsInDEM[0],filename);
	  pointsInDEM[0] = dummyint;
	  // for (i=0;i<pointsInDEM[0];++i) printf("->%i / %i %f %f %f\n",i, pointsInDEM[0], bed[3*i], bed[(3*i)+1], bed[(3*i)+2]);
	} else {
	  fprintf(stderr,"(%s) Error reading file %s\n",argv[0],filename);
	  return EXIT_FAILURE;
	}
	interpolateDEM = 1;
	interpolationScheme = 1;
      }
    }else{
      interpolationScheme = 0;
      printf("(%s) no file %s\n", argv[0], filename);
    }

    // the surface DEM
    strcpy(filename,argv[10]); strcat(filename,"/surf.xyz");
    if (file_exists(filename)){
      pointsInDEM[1] = checkFileEntries(filename, argv, &entriesPerLine);
      if ((pointsInDEM[1] < 1) || (entriesPerLine < 3)){
	fprintf(stderr,"(%s) Wrong file entries (%i,%i) found in file %s\n", argv[0], pointsInDEM[1], entriesPerLine, filename);
	return EXIT_FAILURE;
      }else{
	printf("(%s) %i file entries found in file %s\n", argv[0], pointsInDEM[1], filename);
	surf = (double *) malloc(pointsInDEM[1] * 3 * sizeof(double));
	for (i=0;i<pointsInDEM[1];++i){ surf[3*i]=0.0, surf[3*i+1]=0.0, surf[3*i+2]=0.0;}
	isnoval1 = (int *) malloc( pointsInDEM[1] * sizeof(int));
	dummyint  = readDEM(surf, isnoval1, pointsInDEM[1], filename, argv);
	if (dummyint > 1){
	  printf("(%s) %i out of %i valid points in DEM %s\n",argv[0],dummyint,pointsInDEM[1],filename);
	  pointsInDEM[1] = dummyint;
	} else {
	  fprintf(stderr,"(%s) Error reading file %s\n",argv[0],filename);
	  return EXIT_FAILURE;
	}
      }
      interpolateDEM = 1;
      interpolationScheme += 10;
    }else{      
      printf("(%s) no file %s\n", argv[0], filename);
      if (interpolationScheme == 0) {
	fprintf(stderr,"(%s) We either need a surface or a bedrock\n",argv[0]);
	return EXIT_FAILURE;
      }
    }

    if (interpolationScheme < 11){
      // the thickness DEM
      strcpy(filename,argv[10]); strcat(filename,"/thick.xyz");
      if (file_exists(filename)){
	pointsInDEM[2] = checkFileEntries(filename, argv, &entriesPerLine);
	if ((pointsInDEM[2] < 1) || (entriesPerLine < 3)){
	  fprintf(stderr,"(%s) Wrong file entries (%i,%i) found in file %s\n", argv[0], pointsInDEM[2], entriesPerLine, filename);
	  return EXIT_FAILURE;
	}else{
	  printf("(%s) %i file entries found in file %s\n", argv[0], pointsInDEM[2], filename);
	  thick = (double *) malloc(pointsInDEM[2] * 3 * sizeof(double));
	  for (i=0;i<pointsInDEM[2];++i){ thick[3*i]=0.0, thick[3*i+1]=0.0, thick[3*i+2]=0.0;}
	  isnoval2 = (int *) malloc( pointsInDEM[0] * sizeof(int));
	  dummyint  = readDEM(thick, isnoval2, pointsInDEM[2], filename, argv);
	  if (dummyint > 1){
	    printf("(%s) %i out of %i valid points in DEM %s\n",argv[0],dummyint,pointsInDEM[2],filename);
	    pointsInDEM[2] = dummyint;
	  }else{
	    fprintf(stderr,"(%s) Error reading file %s\n",argv[0],filename);
	    return EXIT_FAILURE;
	  }
	}
	interpolateDEM = 1;
	interpolationScheme += 100;
      }else{
	fprintf(stderr,"(%s) no file %s\n", argv[0], filename);
	fprintf(stderr,"(%s) We either need a thickness file\n",argv[0]);
	return EXIT_FAILURE;
      }
    }
  }else{
    interpolationScheme = 0;
  }

 

  strcpy(inputdirectoryname,argv[1]);
  // read in global header file
  //---------------------------
  strcpy(filename,inputdirectoryname); strcat(filename,"/mesh.header");
  headerfid = fopen(filename, "r");
  if (headerfid==NULL){
    printf("(%s) Error opening file %s\n",argv[0],filename);
    return EXIT_FAILURE;
  }
 
  printf("(%s) reading global header file %s\n",argv[0],filename);
  fscanf(headerfid,"%i %i %i",&pointsinlevel, &elementsinlevel,&belementsinlevel);
  fscanf(headerfid,"%i",&noelementtypes);
  points = lines = triangles = quads = 0;
  for (i=0;i<noelementtypes;++i){
    fscanf(headerfid,"%i %i",&l, &k);
    if (l == 101) points = k;
    else if (l == 202) lines = k;
    else if ( l == 303) triangles = k;
    else if (l == 404) quads = k;
    else {
      printf("(%s) element type %i in header file %s not recognised\n",argv[0],l,filename);
      return EXIT_FAILURE;
    } 
  }
  //compose input and output directory name
  //---------------------------------------
  if (partitions > 1){// we have a parallel mesh
    if (partitions < 10)
      sprintf(cpartition,"%1i", partitions);
    else if (partitions  < 100)
      sprintf(cpartition, "%2i", partitions);
    else if (partitions  < 1000)
      sprintf(cpartition, "%3i", partitions);
    else if (partitions  < 10000)
      sprintf(cpartition, "%4i", partitions);      
    else{
      printf("(%s) More than 9999 partitions? Get real or (if you insist) change the source code of this executable!\n",argv[0]);
      return EXIT_FAILURE;
    }
    strcat(inputdirectoryname,"/partitioning.");strcat(inputdirectoryname,cpartition);
    printf("(%s) Reading from directory: %s\n",argv[0], inputdirectoryname);
    
    strcpy(outputdirectoryname, argv[2]);
    printf("(%s) creating output file directory %s\n", argv[0], outputdirectoryname);
    MKDIR(outputdirectoryname);

    strcat(outputdirectoryname,"/partitioning.");
    strcat(outputdirectoryname,cpartition);  
    printf("(%s) creating output file partition directory %s\n", argv[0], outputdirectoryname);
    MKDIR(outputdirectoryname);

 
  }else{     
    printf("(%s) Reading from directory: %s\n",argv[0], inputdirectoryname);
    
    strcpy(outputdirectoryname, argv[2]);
    printf("(%s) creating output file directory %s\n", argv[0], outputdirectoryname);
    MKDIR(outputdirectoryname);
  }

  strcpy(filename,outputdirectoryname); strcat(filename,"/info.dat");
  infofid = fopen(filename, "w");
  if (infofid==NULL){
    printf("(%s) Error opening file %s\n",argv[0],filename);
    return EXIT_FAILURE;
  }
  // write input to info file
  //-------------------------
  fprintf(infofid,"command line was: %s %s %s %i %8.2f %i %i %i  %5.2f %4.2f  %s %8.2f %3.2f %8.2f [%i]\n",argv[0], argv[1],argv[2],levels,depth, partitions,baseline,GL,percentage,ratio,argv[10],cutOffValue,wexp,noval,corrbed);
  fprintf(infofid,"(%s) was invoked with following inputs ...\n",argv[0]);
  fprintf(infofid,"    ... input directory: %s \n",argv[1]);
  fprintf(infofid,"    ... output directory: %s\n",argv[2]);
  fprintf(infofid,"    ... levels: %i\n",levels);
  fprintf(infofid,"    ... depth: %10.4f\n\n",depth);
  fprintf(infofid,"    ... partitions: %i\n",partitions);
  fprintf(infofid,"    ... baseline: %i\n",baseline);
  fprintf(infofid,"    ... GL: %i\n\n",GL);
  fprintf(infofid,"    ... percentage: %10.4f\n\n",percentage);
  fprintf(infofid,"    ... ratio: %10.4f\n\n",ratio);
  if (interpolateDEM){   
    fprintf(infofid,"    ... DEM directory: %s\n",argv[10]);
    fprintf(infofid,"    ... cutoff radius: %10.4f\n", cutOffValue);
    fprintf(infofid,"    ... exponent m of weight (1/r^m): %10.4f\n",wexp);
    fprintf(infofid,"    ... noval: %10.4f\n\n",noval); 
    fprintf(infofid,"    ... corrbed: %1i\n\n",corrbed);
  }

  fprintf(infofid,"    ... total poitns in 2d mesh: %i \n",pointsinlevel);
  fprintf(infofid,"    ... total elements in 2d mesh: %i \n",elementsinlevel);
  fprintf(infofid,"    ... total boundary elements in 2d mesh: %i \n",belementsinlevel);
  fprintf(infofid,"    ... no of element type 101: %i \n",points);
  fprintf(infofid,"    ... no of element type 202: %i \n",lines);
  fprintf(infofid,"    ... no of element type 303: %i \n",triangles);
  fprintf(infofid,"    ... no of element type 404: %i \n",quads);
  printf("    ... total points in 2d mesh: %i \n",pointsinlevel);
  printf("    ... total elements in 2d mesh: %i \n",elementsinlevel);
  printf("    ... total boundary elements in 2d mesh: %i \n",belementsinlevel);
  printf("    ... no of element type 101: %9i \n",points);
  printf("    ... no of element type 202: %9i \n",lines);
  printf("    ... no of element type 303: %9i \n",triangles);
  printf("    ... no of element type 404: %9i \n",quads);
  fclose(headerfid);

  // write global header file
  //-------------------------
  strcpy(filename,outputdirectoryname); strcat(filename,"/mesh.header.global");
  headerfid = fopen(filename, "w");
  if (headerfid==NULL){
    printf("(%s) Error opening file %s\n",argv[0],filename);
    return EXIT_FAILURE;
  }
  fprintf(headerfid,"%9i %9i %9i\n",pointsinlevel*levels,elementsinlevel*(levels-1),elementsinlevel*2 + belementsinlevel*(levels-1));
  //  fprintf(headerfid,"%i\n",noelementtypes);
  fprintf(headerfid,"%9i\n", (points > 0) + 2*(triangles > 0) + (quads > 0) +1);
  if (points > 0)  fprintf(headerfid,"202 %i\n",points * (levels-1)) + lines * baseline;
  if ((lines > 0) || (quads > 0))  fprintf(headerfid,"404 %i\n",lines * (levels-1)  +  2 * quads);
  if (triangles > 0)  fprintf(headerfid,"706 %i\n", triangles * (levels-1));
  if (quads > 0)  fprintf(headerfid,"808 %i\n", quads * (levels-1));


  if (partitions > 1){
    // extrude meshes for each partition
    //----------------------------------
    for (partition=0;partition<partitions;++partition){
      if (partition < 10)
	sprintf(cpartition,"/part.%1i", partition+1);
      else if (partition  < 100)
	sprintf(cpartition, "/part.%2i", partition+1);
      else if (partition  < 1000)
	sprintf(cpartition, "/part.%3i", partition+1);
      else
	sprintf(cpartition, "/part.%4i", partition+1);  
      strcpy(inputfilename,inputdirectoryname); strcat(inputfilename,cpartition);
      printf("(%s) reading partition %i: base file name %s\n",argv[0],partition+1,inputfilename);
      strcpy(outputfilename,outputdirectoryname); strcat(outputfilename,cpartition);
      stat = extrudepartition(argv, inputfilename, outputfilename, partition, partitions, levels, depth, pointsinlevel, elementsinlevel, belementsinlevel, interpolationScheme, bed, surf, thick, pointsInDEM, cutOffValue, wexp, isnoval0, isnoval1, isnoval2);
      //stat = extrudepartition(argv, inputfilename, outputfilename, partition, partitions, levels, depth, pointsinlevel, elementsinlevel, belementsinlevel, interpolateDEM=0, charptr);
      
    } // end of loop over partitions
  }else{
    printf("(%s) Doing serial extrusion\n",argv[0]);    
    strcpy(inputfilename,inputdirectoryname);
    strcpy(outputfilename,outputdirectoryname);
    stat = extrudeserial(argv, inputfilename, outputfilename, levels, depth, pointsinlevel, elementsinlevel, belementsinlevel, interpolationScheme, bed, surf, thick, pointsInDEM, cutOffValue, wexp, GL, percentage, ratio, baseline, corrbed, isnoval0, isnoval1, isnoval2);
    
  }
    
  

  fclose(headerfid);
  fclose(infofid);
  if (interpolationScheme > 0){
    free(bed);
    free(surf);
  }
  return stat;

}
