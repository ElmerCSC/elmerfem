#include <stdio.h>
#include "mpi.h"

void x_(long int *c)
{
  	fprintf( stdout, "eh                 %ld %ld\n", c, MPI_IN_PLACE );
}
