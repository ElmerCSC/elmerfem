#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define MAX_LENGTH 255

int main(int argc, char **argv)
{
   FILE *fp = fopen( argv[1], "r" );
   char line[MAX_LENGTH+1], *ptr;

   static double norm, f, target_nrm, target_eps = 1e-5, REALT, CPUT;
   int n, success;

   success = 0;
   while( fgets( line, MAX_LENGTH, fp ) )
   {
      if ( strstr( line, "END TEST CASE" ) ) {
         ptr = strstr( line, "NRM=" );         
         if ( ptr ) sscanf( ptr,"NRM=%lf", &target_nrm );
         ptr = strstr( line, "EPS=" );         
         if ( ptr ) sscanf( ptr,"EPS=%lf", &target_eps );         
         success = compare( norm, target_nrm, target_eps );
         if ( !success ) {
           n = strlen(line)-1;
           while( line[n]==10 || line[n]==13 ) line[n--] = '\0';
           fprintf( stderr, "[FAILED]: %s, Computed NRM=%g: ", line, norm );
           break;
         }
      }
      else if ( ptr=strstr( line, "(NRM,RELC)" ) ) {
        while( *ptr != '\0' && *ptr!='+' && *ptr != '-' && *ptr != '.' && !isdigit(*ptr) ) ptr++;
        n = sscanf( ptr, "%lf", &f );
        if ( n==1 && f != 0.0 ) norm = f;
      } else if (ptr = strstr( line, "Check NRM=") ) {
        while( *ptr != '\0' && *ptr!='+' && *ptr != '-' && *ptr != '.' && !isdigit(*ptr) ) ptr++;
        n = sscanf( ptr, "%lf", &f );
        if ( n==1 && f != 0.0 ) norm = f;
      } 
      else if ( ptr = strstr( line, "(CPU,REAL):" )  ) {
        sscanf(  ptr+11, "%lf %lf", &CPUT, &REALT); 
      }
   }

   if ( argc<=2 ) {
     fprintf( stdout, "%d\n", success );
   } else {
     fprintf( stdout, "%g\n", atof(argv[2])+CPUT );
   }
}



int compare( double norm1, double norm2, double eps )
{
   int n;

   if ( norm1 != -1 ) {
      if ( eps < 0 ) { 
         if ( norm2 < norm1 )
            return 1;
         else
            return 0;
      } else  if ( 2 * fabs(norm1-norm2) / (norm1+norm2) < eps )
         return 1;
      else
         return 0;
   } else {
      return 1;
   }
}
