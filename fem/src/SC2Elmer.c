/*
 *  A program to convert Silicomp mesh file to elmer mesh.
 *
 *  Usage: SC2Elmer silicomp_file
 *
 *  the result is files: mesh.header, mesh.nodes, mesh.elements, mesh.boundary
 *  containing the mesh in elmer format.
 *
 *  boundaries may be grouped based on boundary normals if file groups.dat
 *  is found in current directory. The groups.dat file consist of following:
 * 
 *  n tolerance_angle
 *  nx_1 ny_1
 *  ...
 *  ...
 * 
 *  nx_n ny_n
 *
 * tolerance angle is given in degrees.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_LEN 1000
char buf[MAX_LEN+1];

typedef struct edge
{
   struct edge *next;

   int n1,n2,n3,p1,p2;
} edge_t;

double *nx, *ny, *fx;

/*
 *  return next line from file.
 */
char *scan( FILE *fp )
{
   char *ptr;

   *buf = ' ';
   while( 1 )
   {
      if ( feof(fp) || ferror(fp) ) break;
      fgets( buf, MAX_LEN, fp );

      ptr = buf;
      while( *ptr != '\0' && *ptr != ' ' ) ptr++;

      if ( *ptr != '\0' ) break;
   }

   return buf;
}

/*
 *  Find mesh edges.
 */
int findedges( int n, int *el[8], edge_t **edge )
{
   int i, j, k, e1,e2, e3, edges, n1,n2,n3;
   edge_t *ptr, *ptr0;

   int map[4][3] = { { 0,1,4 }, { 1,2,5 }, { 2,3,6 }, { 3,0,7 } };

   edges = 0;
   for( i=0; i<n; i++ )
   {
      for( j=0; j<4; j++ ) {
         n1 = el[map[j][0]][i];
         n2 = el[map[j][1]][i];
         n3 = el[map[j][2]][i];

         if ( n2 < n1 ) { k=n1; n1=n2; n2=k; }

         ptr0 = edge[n1];
         ptr  = edge[n1];
         while( ptr )
         {
            if ( ptr->n2 == n2 ) {
              ptr->p2 = i;
              break;
            }
            ptr0 = ptr;
            ptr = ptr->next;
         }

         if ( !ptr  ) {
            if ( !ptr0 ) {
               ptr = edge[n1] = (edge_t *)calloc( sizeof( edge_t ), 1 );
            } else {
               ptr = ptr0->next = (edge_t *)calloc( sizeof( edge_t ), 1 );
            }

            ptr->n1 = n1;
            ptr->n2 = n2;
            ptr->n3 = n3;
            ptr->p1 =  i;
            ptr->p2 = -1;

            edges++;
         }
      }
   }
   fprintf( stderr, "found edges: %d\n", edges );

   return edges;
}

/*
 *  Compute boundary element outer normal
 */
void normal( double *n, int p1, int p2, double cx, double cy )
{
    double dxdu, dydu, detA;

    dxdu = nx[p2] - nx[p1];
    dydu = ny[p2] - ny[p1];

    detA = dxdu*dxdu + dydu*dydu;
    detA = 1.0 / sqrt(detA);

    n[0] = -dydu * detA;
    n[1] =  dxdu * detA;

    cx = cx - ( nx[p1] + nx[p2] ) / 2;
    cy = cy - ( ny[p1] + ny[p2] ) / 2;

    if ( n[0]*cx + n[1]*cy > 0 ) { n[0] = -n[0]; n[1] = -n[1]; }
}

int main( int argc, char **argv )
{
   FILE *fp = fopen( argv[1], "r" ), *fp_out, *fp_grp;

   char *line;

   int i,j,k,i1,i2,i3,i4,i5,i6,i7,i8,nodes,elements,*ref, *renumber, *el[8], grp, groups;
   edge_t **edge, *ptr;
   double x,y,n[2],cx,cy, g, *groups_x, *groups_y, s, twopi, ang;

   line = scan( fp );
   sscanf( line, "%d %d", &elements, &nodes );

   for( i=0; i<8; i++ )
      el[i] = (int *)malloc( elements*sizeof(int) );

   edge = (edge_t **)calloc( nodes,sizeof(edge_t *) );

   ref = (int *)calloc( nodes, sizeof(int) );
   renumber = (int *)calloc( nodes, sizeof(int) );

   nx = (double *)calloc( nodes, sizeof(double) );
   ny = (double *)calloc( nodes, sizeof(double) );
   fx = (double *)calloc( nodes, sizeof(double) );

   /* 
    * Read the elements from silicomp format file, store in memory
    */
   for( i=0; i<elements; i++ )
   {
      line = scan( fp );
      j = sscanf( line, "%d %d %d %d %d %d %d %d %d", &k, &i1,&i2,&i3,&i4,&i5,&i6,&i7,&i8 );
      if ( j != 9 ) break;

      ref[i1] = 1;
      ref[i2] = 1;
      ref[i3] = 1;
      ref[i4] = 1;
      ref[i5] = 1;
      ref[i6] = 1;
      ref[i7] = 1;
      ref[i8] = 1;

      el[0][i] = i1;
      el[1][i] = i2;
      el[2][i] = i3;
      el[3][i] = i4;
      el[4][i] = i5;
      el[5][i] = i6;
      el[6][i] = i7;
      el[7][i] = i8;
   }
   if ( i>=elements) scan( fp );
   elements = i;

   /* 
    * Do a continuous renumbering of nodal points.
    */
   k = 0;
   for( i=0; i<nodes; i++ )
   {
      if ( ref[i] ) { renumber[i] = k++; }
   }

   /* 
    * Output elements in ELMER format.
    */
   fp_out = fopen( "mesh.elements", "w" );
   for( i=0; i<elements; i++ )
   {
      for ( j=0; j<8; j++ ) el[j][i] = renumber[ el[j][i] ];

      fprintf( fp_out, "%d 1 408 %d %d %d %d %d %d %d %d\n", i+1,
               el[0][i]+1, 
               el[1][i]+1, 
               el[2][i]+1, 
               el[3][i]+1, 
               el[4][i]+1, 
               el[5][i]+1, 
               el[6][i]+1, 
               el[7][i]+1 );
   }
   fclose( fp_out );

   /* 
    * Read nodal points from Silicomp format file.
    */
   for( i=0; i<nodes; i++ )
   {
      if ( 3 != sscanf( line, "%d %lf %lf", &i1, &x, &y ) ) break;
      line = scan( fp );

      nx[renumber[i1]] = x;
      ny[renumber[i1]] = y;
   }
   if ( i>=nodes ) scan( fp );
   nodes = i;

   /* 
    * Write nodal points in ELMER format.
    */
   fp_out = fopen( "mesh.nodes", "w" );
   for( i=0; i<nodes; i++ )
   {
      fprintf( fp_out, "%d -1 %g %g 0.0\n", i+1, nx[i], ny[i] );
   }
   fclose( fp_out );

   /* 
    * Read temperature field from Silicomp format file.
    */
   for( i=0; i<nodes; i++ )
   {
      sscanf( line, "%d %lf", &i1, &y );
      fx[renumber[i1]] = y;
      line = scan( fp );
   }

   /* 
    * Write temperature field in ELMER format.
    */
   fp_out = fopen( "result.dat", "w" );
   fprintf( fp_out,
          "Degrees of freedom:\n"
          "temperature  1  :heat equation\n"
          "Total DOFs:   1\n"
          "Time:       1      1 0.100000000000E+001\n"
          "temperature\n" );

   for( i=0; i<nodes; i++ )
   {
      fprintf( fp_out,  "%d %d %g\n", i+1,i+1,fx[i] );
   }
   fclose( fp_out );

   /* 
    * Find all edges of the mesh. Boundaries are identified as edges
    * who have only one parent element.
    */
   k = findedges( elements, el,edge );

   /* 
    * Read in groups.dat file which gives normals to organize
    * boundaries to separate groups:
    * 
    * * if given group normal and outward boundary normal are 
    *    within given angle add boundary element to boundary group.
    */
   fp_grp = fopen( "groups.dat", "r" );
   groups = 0;
   if ( fp_grp ) {
      twopi = 4 * acos(0.0);
      ang = 0.0;
      fscanf( fp_grp, "%d %lf", &groups, &ang );
      ang = 1.0 - cos( ang * twopi / 360.0 );
      groups_x = (double *)malloc( groups*sizeof(double) );
      groups_y = (double *)malloc( groups*sizeof(double) );
      for( i=0; i<groups; i++ )
      {
        fscanf( fp_grp, "%lf %lf", &groups_x[i], &groups_y[i] );
        s = sqrt( groups_x[i]*groups_x[i] + groups_y[i]*groups_y[i] );
        groups_x[i] /= s;
        groups_y[i] /= s;
      }
      fclose(fp_grp);
   }

   /*
    * Output boundary elements in ELMER format.
    */
   fp_out = fopen( "mesh.boundary", "w" );
   j = 0;
   for( i=0; i<nodes; i++ )
   {
      ptr = edge[i];
      while( ptr ) {
        if ( ptr->p2 < 0 ) {
           j++;
           k = ptr->p1;
           cx = (nx[el[0][k]] + nx[el[1][k]] + nx[el[2][k]] + nx[el[3][k]]) / 4.0;
           cy = (ny[el[0][k]] + ny[el[1][k]] + ny[el[2][k]] + ny[el[3][k]]) / 4.0;

           normal( n, ptr->n1, ptr->n2, cx,cy  );
           for( k=0; k<groups; k++ )
           {
              if ( fabs( 1 - n[0]*groups_x[k] + n[1]*groups_y[k] ) < ang+1.0e-10 ) break;
           }
           fprintf( fp_out, "%d %d %d 0 203 %d %d %d\n", j,k+1, ptr->p1+1,ptr->n1+1,ptr->n2+1,ptr->n3+1 );
        }
        ptr = ptr->next;
      }
   }
   fclose( fp_out );

   /*
    * Output mesh header in ELMER format.
    */
   fp_out = fopen( "mesh.header", "w" );
   fprintf( fp_out, "%d %d %d\n", nodes, elements, j );

   fprintf( fp_out, "2\n408 %d\n203 %d\n", elements, j );
   fclose( fp_out );
}
