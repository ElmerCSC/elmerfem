SUBROUTINE TimeIntTest( Model,Solver,dt,TransientSimulation )
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element

  LOGICAL :: Found

  INTEGER :: n
  REAL(KIND=dp) :: Norm

  TYPE(ValueList_t), POINTER :: BodyForce
  REAL(KIND=dp) :: STIFF(1,1), MASS(1,1), LOAD(1), FORCE(1)
!------------------------------------------------------------------------------

   !Initialize the system and do the assembly:
   !------------------------------------------
   CALL DefaultInitialize()

   Element => GetActiveElement(1)
   n = GetElementNOFNodes()
   LOAD = 0.0d0

   BodyForce => GetBodyForce()
   IF ( ASSOCIATED(BodyForce) ) &
      Load(1:n) = GetReal( BodyForce, 'Source', Found )

   !Get element local matrix and rhs vector:
   !----------------------------------------
   STIFF = 0.0d0
   MASS  = 1.0d0
   FORCE = LOAD

   !Update global matrix and rhs vector from local matrix & vector:
   !---------------------------------------------------------------
   CALL Default1stOrderTime( MASS, STIFF, FORCE )
   CALL DefaultUpdateEquations( STIFF, FORCE )
   CALL DefaultFinishBulkAssembly()
   CALL DefaultFinishBoundaryAssembly()
   CALL DefaultFinishAssembly()
!
!  Solve the system and we are done:
!  ---------------------------------
   Norm = DefaultSolve()
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE TimeIntTest
!------------------------------------------------------------------------------
