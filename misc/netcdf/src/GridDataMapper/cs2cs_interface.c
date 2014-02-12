/***** Vili Forsell, 28.6.2011
* Function for interfacing with cs2cs (i.e. PROJ.4)
* See "/user/include/proj_api.h" for cs2cs application interface
*/
#include <stdlib.h>
#include <proj_api.h>
#include <stdio.h>

char const MSG_HEADER[] = "GridDataMapper: CS2CS coordinate system transformation:";

// Add boolean values
enum bool {
	true = 1,
	false = 0
	};
typedef enum bool bool;

/* Function: cs2cs_transform()
*** Transforms the given input from a coordinate system to another with cs2cs
@param coord The input (Elmer) coordinates to be transformed
@param hasZ If 0, the coord parameter contains 3rd dimension, otherwise 2D
@param isRad If 0, then uses radians during transformation, otherwise degrees
@param elmer_proj The cs2cs parameter string for Elmer grid coordinate system (from)
@param netcdf_proj The cs2cs parameter string for NetCDF grid coordinate system (to)
@param res The resulting (NetCDF) data values for further indexing
*/
void cs2cs_transform( double coord[3], int hasZ, int isRad, char elmer_proj[], char netcdf_proj[], double res[3] ) {
	projPJ pj_elmer = NULL, pj_netcdf = NULL; // Coordinate system definitions
	bool hasErrors = false;
	int transf_val = 0; // For return value from the transformation predefined by cs2cs as 0, if no errors, and as error code otherwise
	bool isZset = true; // True, if z coordinate is used
	bool isInputRad = false; // True, if the given input is in radians

	// Set constants
	long point_count = 1; // Number of processed points (x,y,z can be arrays)
	int point_offset = 1; // Step size, if multiple points (f.ex. every second) processed

	// Input information
	if ( hasZ == 0 ) isZset = false;
	if ( isRad != 0 ) isRad = true;

/* DEBUG PRINTOUTS
	printf("(%.2f, %.2f, %.2f) ; hasZ = %d ; isRad = %d\n",coord[0],coord[1],coord[2],hasZ,isRad);
	fprintf(stdout,"\n\"%s\"\n",elmer_proj);
	fprintf(stdout,"\n\"%s\"\n",netcdf_proj);
*/

	// Initializes the coordinate systems for the read Elmer data point and the NetCDF data
	pj_elmer = pj_init_plus(elmer_proj);
	if ( pj_elmer == NULL ) {
		printf("%s Failed to initialize Elmer coordinate system with parameters: %s!\n", MSG_HEADER, elmer_proj);
		hasErrors = true;
	}
	pj_netcdf = pj_init_plus(netcdf_proj);
	if ( pj_netcdf == NULL ) {
		printf("%s Failed to initialize NetCDF coordinate system with parameters: %s!\n", MSG_HEADER, netcdf_proj);
		hasErrors = true;
	}

	// Transformation from degrees to radians, if necessary
	if ( !isInputRad ) {
		printf("%s Converting values from degrees to radians.\n", MSG_HEADER);
		coord[0] *= DEG_TO_RAD;
		coord[1] *= DEG_TO_RAD;
		coord[2] *= DEG_TO_RAD;
	}

	// Performs the transformation (according to input dimensions)
	if ( hasErrors != true ) {
		if (isZset) transf_val = pj_transform(elmer_proj,netcdf_proj,point_count,point_offset,&coord[0],&coord[1],&coord[2]);
		else transf_val = pj_transform(elmer_proj,netcdf_proj,point_count,point_offset,&coord[0],&coord[1],NULL);
		if ( transf_val != 0 ) {
			printf("%s Failed to transform coordinates; Proj.4 error code %d.\n",MSG_HEADER,transf_val);
			hasErrors = true;
		}
	}

	// Saves output, if no errors
	if ( hasErrors == true ) {
		res[0] = 0;
		res[1] = 0;
		res[2] = 0;
		abort(); // Immediate abort, since this is likely not intended and affects all results
        } else {
		res[0] = coord[0];
		res[1] = coord[1];
		res[2] = coord[2];
        }
  
/* DEBUG PRINTOUTS
	printf("Done with (%.2f, %.2f)\n\n", coord[0], coord[1]);
*/
}
