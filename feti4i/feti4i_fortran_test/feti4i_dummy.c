
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../feti4i.h"

struct FETI4IStructMatrix {
  FETI4IInt     indexBase;
  FETI4IInt     size;
  FETI4IInt*    indices;
  FETI4IReal*   values;
	FETI4IInt     nelems;
};

struct FETI4IStructInstance {
  FETI4IMatrix  matrix;
  FETI4IReal* 	rhs;
  FETI4IInt* 		l2g;                     /* length of both rhs and l2g is size */ 
  FETI4IMPIInt 	neighbours_size;
  FETI4IMPIInt*	neighbours;
  FETI4IInt 		dirichlet_size;
  FETI4IInt* 		dirichlet_indices;  //TODO which numbering? we prefer global numbering
  FETI4IReal* 	dirichlet_values;
};


void FETI4ICreateStiffnessMatrix(FETI4IMatrix *matrix, FETI4IInt indexBase)
{
  *matrix = (FETI4IMatrix) malloc(sizeof(struct FETI4IStructMatrix));
	(*matrix)->indexBase = indexBase;
	(*matrix)->size = 1;
	(*matrix)->indices = NULL;
	(*matrix)->values = NULL;
	(*matrix)->nelems = 0;
	printf("%s: created matrix %p, indexBase %d\n", __func__, (void*) *matrix, indexBase);
}


void FETI4IAddElement(FETI4IMatrix matrix, FETI4IInt size, FETI4IInt* indices, FETI4IReal* values)
{
  matrix->size += (size-1);
	matrix->indices = indices;
	matrix->values = values;
	printf("%s: added element %d to matrix %p\n", __func__, matrix->nelems++, (void*) matrix);
}

void FETI4ICreateInstance(FETI4IInstance *instance, FETI4IMatrix matrix, FETI4IInt size,
  FETI4IReal* rhs, FETI4IInt* l2g, FETI4IMPIInt neighbours_size, FETI4IMPIInt* neighbours,
  FETI4IInt dirichlet_size, FETI4IInt* dirichlet_indices, FETI4IReal* dirichlet_values)
{
	printf("%s: setting rhs %p (%d), l2g %p (%d), " 
    "neighbours %p (%d), dirichlet_indices %p (%d), "
    "dirichlet_values %p (%d)\n",
    __func__, (void*) rhs, size, (void*) l2g, size,
    (void*) neighbours, neighbours_size, (void*) dirichlet_indices, dirichlet_size,
    (void*) dirichlet_indices, dirichlet_size);

  *instance = (FETI4IInstance) malloc(sizeof(struct FETI4IStructInstance));
  (*instance)->matrix = matrix;
  assert(size == matrix->size);
  (*instance)->rhs = rhs;
  (*instance)->l2g = l2g;
  (*instance)->neighbours_size = neighbours_size;
  (*instance)->neighbours = neighbours;
  (*instance)->dirichlet_size = dirichlet_size;
  (*instance)->dirichlet_indices = dirichlet_indices;
  (*instance)->dirichlet_values = dirichlet_values;
}

void FETI4ISolve(FETI4IInstance instance, FETI4IInt solution_size, FETI4IReal* solution)
{
  FETI4IInt i;
  FETI4IMatrix matrix = instance->matrix;
  FETI4IReal *rhs = instance->rhs;

	printf("%s: solving with matrix %p, rhs %p\n", __func__, (void*) instance->matrix, (void*) instance->rhs);
  assert(solution_size == matrix->size);
  for (i=0; i<solution_size; i++) {
    printf("#%d %g\n",i,rhs[i]);
  }
  printf("\n");
  for (i=0; i<solution_size; i++) {
    printf("#%d %g ",i,solution[i]);
    solution[i]=-i;
    printf(" %g\n",solution[i]);
  }
}

void FETI4IDestroy(void *ptr)
{
  FETI4IInstance instance = (FETI4IInstance) ptr;
  free(instance->matrix);
  free(instance);
}

