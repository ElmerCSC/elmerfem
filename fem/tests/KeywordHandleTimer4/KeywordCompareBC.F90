!------------------------------------------------------------------------------
!> Compare keyword strategies with and without a Handle. 
!> This test for logical, integer, and string.
!------------------------------------------------------------------------------
SUBROUTINE KeywordCompare( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  INTEGER :: t, NoActive
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Material, BodyForce
  TYPE(ValueHandle_t) :: A_h, B_h, C_h, D_h
  LOGICAL :: ValL, StrL
  INTEGER :: ValI
  REAL(KIND=dp) :: ValR 
  CHARACTER(LEN=MAX_NAME_LEN) :: ValC

  INTEGER :: i, TimeSweeps
  INTEGER :: OldCounter(4),NewCounter(4)
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp) :: OldTime, NewTime, PseudoNorm
  TYPE(ValueList_t), POINTER :: BC

  
  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  CALL Info( 'KeywordCompare','')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')
  CALL Info( 'KeywordCompare','Comparing keyword strategies for BCs')

  Mesh => GetMesh()

  
  NoActive = GetNOFBoundaryElements()
  TimeSweeps = ListGetInteger( Solver % Values,'Timer Sweeps')

  OldTime = 0.0_dp
  NewTime = 0.0_dp
  OldCounter = 0
  NewCounter = 0

  
  CALL Info('KeywordCompare','Starting bc keyword comparision')
  OldTime = RealTime()
  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetBoundaryElement(t)
      IF( .NOT. ActiveBoundaryElement()) CYCLE
      
      CurrentModel % CurrentElement => Element
      BC => GetBC(Element)
      
      IF(.NOT. ASSOCIATED( BC ) ) CYCLE
      
      ValL = ListGetLogical( BC,'Logical Value A',Found )
      ValI = ListGetInteger( BC,'Integer Value B',Found )
      ValR = ListGetCReal( BC,'Float Value C',Found )
      ValC = ListGetString( BC,'String Value D',Found )
      StrL = ( ValC == 'string target d') 
      
      IF( ValL ) OldCounter(1) = OldCounter(1) + 1
      OldCounter(2) = OldCounter(2) + ValI
      OldCounter(3) = OldCounter(3) + ValR
      IF( StrL ) OldCounter(4) = OldCounter(4) + 1 
    END DO
  END DO
  OldTime = RealTime() - OldTime
  CALL Info('KeywordCompare','finished testing old keyword style')


  NewTime = RealTime()  
  CALL ListInitElementKeyword( A_h,'Boundary Condition','Logical Value A' )
  CALL ListInitElementKeyword( B_h,'Boundary Condition','Integer Value B' )
  CALL ListInitElementKeyword( C_h,'Boundary Condition','Float Value C' )
  CALL ListInitElementKeyword( D_h,'Boundary Condition','String Value D' )

  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetBoundaryElement(t)
      IF( .NOT. ActiveBoundaryElement()) CYCLE
     
      ValL = ListGetElementLogical( A_h, Element, Found )
      ValI = ListGetElementInteger( B_h, Element, Found )
      ValR = ListGetElementReal( C_h, Found = Found )
      StrL = ListCompareElementString( D_h, 'string target d',Element, Found )
      
      IF( ValL ) NewCounter(1) = NewCounter(1) + 1
      NewCounter(2) = NewCounter(2) + ValI
      NewCounter(3) = NewCounter(3) + ValR
      IF( StrL ) NewCounter(4) = NewCounter(4) + 1
    END DO
  END DO
  NewTime = RealTime() - NewTime


  PRINT *,'Old Counter:',OldCounter
  PRINT *,'New Counter:',NewCounter

  PRINT *,'Speed-up for BC operations:           ',OldTime / NewTime

  PseudoNorm = 1.0_dp * SUM(NewCounter) / SUM(OldCounter)
  PRINT *,'PseudoNorm:',PseudoNorm

  Solver % Variable % Values = PseudoNorm
  Solver % Variable % Norm = PseudoNorm
        
  CALL Info( 'KeywordCompare','Add done')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')


END SUBROUTINE KeywordCompare
