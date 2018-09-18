SUBROUTINE OptimizeMesh_MMG3D(Model,Solver,dt,TransientSimulation)

#ifdef HAVE_MMG
  USE MMG3DTools
#endif

  IMPLICIT NONE
#ifdef HAVE_MMG
#include "mmg/mmg3d/libmmgtypesf.h"
#endif

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

#ifdef HAVE_MMG

  TYPE(Mesh_t),POINTER :: Mesh,NewMesh
  TYPE(ValueList_t), POINTER :: SolverParams
  REAL(KIND=dp) :: hsiz(3)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SolverName = "OptimizeMesh_MMG3D"
  SolverParams => GetSolverParams()

  hsiz = (/2000.0,2000.0,150.0/)
  Mesh => Solver % Mesh
  mmgMesh = 0
  mmgSol  = 0

  CALL MMG3D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  CALL SET_MMG3D_MESH(Mesh)

  IF (DEBUG) PRINT *,'--**-- SET MMG3D PARAMETERS '
  CALL SET_MMG3D_PARAMETERS(SolverParams)

  CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ier)
  IF ( ier == MMG5_STRONGFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
    STOP MMG5_STRONGFAILURE
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG3DLIB"
  ENDIF
  IF (DEBUG) PRINT *,'--**-- MMG3D_mmg3dlib DONE'

  !! GET THE NEW MESH
  CALL GET_MMG3D_MESH(NewMesh)

#else
  CALL Fatal(SolverName, "Remeshing utility MMG3D has not been installed")
#endif
END SUBROUTINE OptimizeMesh_MMG3D
