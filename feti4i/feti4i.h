
#ifndef FETI4I_H_
#define FETI4I_H_

#include "mpi.h"

/*-----------------------------------------------------------------------------
 Set data-types used in FETI4I

 Possible values for FETI4I_INT_WIDTH:
   32: use 32 bit signed integer
   64: use 64 bit signed integer

 Possible values for FETI4I_REAL_WIDTH:
   64: ESPRESO supports only 64 bit real values
------------------------------------------------------------------------------*/
#ifndef FETI4I_INT_WIDTH
#define FETI4I_INT_WIDTH 32
#endif

#ifndef FETI4I_REAL_WIDTH
#define FETI4I_REAL_WIDTH 64
#endif

#if FETI4I_INT_WIDTH == 32
	typedef int FETI4IInt;
#elif FETI4I_INT_WIDTH == 64
	typedef long FETI4IInt;
#else
#error "Incorrect user-supplied value of FETI4I_INT_WIDTH"
#endif

/* MPI integer (e.g. rank) are always 32-bit */
	typedef int FETI4IMPIInt;

#if FETI4I_REAL_WIDTH == 64
	typedef double FETI4IReal;
#else
#error "Incorrect user-supplied value of FETI4I_REAL_WIDTH"
#endif


#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 Definitions of internal structures used in FETI4I
------------------------------------------------------------------------------*/
typedef struct FETI4IStructMatrix* FETI4IMatrix;
typedef struct FETI4IStructInstance* FETI4IInstance;

/*-----------------------------------------------------------------------------
 Functions for manipulating with FETI4I internal structures
------------------------------------------------------------------------------*/

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix, //TODO size?
		FETI4IInt		indexBase
);

void FETI4IAddElement(   //TODO CRS or CCS ?
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IInt* 		indices,
		FETI4IReal* 	values
);

/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

void FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IReal* 	rhs,
		FETI4IInt* 		l2g,                     /* length of both rhs and l2g is size */ 
		FETI4IMPIInt 	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,  //TODO which numbering? we prefer global numbering
		FETI4IReal* 	dirichlet_values
);

void FETI4ISolve(
		FETI4IInstance 	instance,
		FETI4IInt 		solution_size,
		FETI4IReal*		solution
);

/*-----------------------------------------------------------------------------
 Functions for updating a created instance
------------------------------------------------------------------------------*/

void FETI4IUpdateStiffnessMatrix(
		FETI4IInstance 	instance,
		FETI4IMatrix 	stiffnessMatrix
);

void FETI4IUpdateRhs(
		FETI4IInstance 	instance,
		FETI4IInt 		rhs_size,
		FETI4IReal* 	rhs_values
);

//TODO VH: neighbours perhaps should not be passed here and they are FETI4IMPIInt
void FETI4IUpdateDirichlet(
		FETI4IInstance 	instance,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IReal* 	l2g,
		FETI4IInt 		neighbour_size,
		FETI4IReal* 	neighbour
);


/*-----------------------------------------------------------------------------
 Destroy an arbitrary internal structure
------------------------------------------------------------------------------*/

void FETI4IDestroy(
		void* 			ptr
);


#ifdef __cplusplus
}
#endif


#endif /* FETI4I_H_ */

