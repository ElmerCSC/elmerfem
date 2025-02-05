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
  TYPE(ValueHandle_t) :: LogA_h, LogB_h, LogC_h, LogD_h
  TYPE(ValueHandle_t) :: IntA_h, IntB_h, IntC_h, IntD_h
  TYPE(ValueHandle_t) :: StrA_h, StrB_h
  LOGICAL :: LogA, LogB, LogC, LogD
  INTEGER :: IntA, IntB, IntC, IntD
  CHARACTER(LEN=MAX_NAME_LEN) :: StrA, StrB     

  INTEGER :: i, TimeSweeps
  INTEGER :: OldCounter(4),NewCounter(4)
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp) :: OldTimes(4), NewTimes(4), PseudoNorm
  
  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  CALL Info( 'KeywordCompare','')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')
  CALL Info( 'KeywordCompare','Comparing keyword strategies')

  Mesh => GetMesh()

  
  NoActive = GetNOFActive()
  TimeSweeps = ListGetInteger( Solver % Values,'Timer Sweeps')

  OldTimes = 0.0_dp
  NewTimes = 0.0_dp
  OldCounter = 0
  NewCounter = 0


  CALL Info('KeywordCompare','Starting logical comparision')
  OldTimes(1) = RealTime()
  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)
      Material => GetMaterial()
      BodyForce => GetBodyForce()

      LogA = ListGetLogical( Material,'Logical Value A',Found )
      LogB = ListGetLogical( Material,'Logical Value B',Found )
      LogC = ListGetLogical( Material,'Logical Value C',Found )
      LogD = ListGetLogical( BodyForce,'Logical Value D',Found )

      IF( LogA ) OldCounter(1) = OldCounter(1) + 1
      IF( LogB ) OldCounter(1) = OldCounter(1) + 1
      IF( LogC ) OldCounter(1) = OldCounter(1) + 1
      IF( LogD ) OldCounter(1) = OldCounter(1) + 1
    END DO
  END DO
  OldTimes(1) = RealTime() - OldTimes(1)


  NewTimes(1) = RealTime()  
  CALL ListInitElementKeyword( LogA_h,'Material','Logical Value A' )
  CALL ListInitElementKeyword( LogB_h,'Material','Logical Value B' )
  CALL ListInitElementKeyword( LogC_h,'Material','Logical Value C')
  CALL ListInitElementKeyword( LogD_h,'Body Force','Logical Value D' )

  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)

      LogA = ListGetElementLogical( LogA_h, Element, Found )
      LogB = ListGetElementLogical( LogB_h, Element, Found )
      LogC = ListGetElementLogical( LogC_h, Found = Found )
      LogD = ListGetElementLogical( LogD_h )

      IF( LogA ) NewCounter(1) = NewCounter(1) + 1
      IF( LogB ) NewCounter(1) = NewCounter(1) + 1
      IF( LogC ) NewCounter(1) = NewCounter(1) + 1
      IF( LogD ) NewCounter(1) = NewCounter(1) + 1
    END DO
  END DO
  NewTimes(1) = RealTime() - NewTimes(1)



  CALL Info('KeywordCompare','Starting integer comparision')
  OldTimes(2) = RealTime()
  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)
      Material => GetMaterial()
      BodyForce => GetBodyForce()

      IntA = ListGetInteger( Material,'Integer Value A',Found )
      IntB = ListGetInteger( Material,'Integer Value B',Found )
      IntC = ListGetInteger( Material,'Integer Value C',Found )
      IntD = ListGetInteger( BodyForce,'Integer Value D',Found )

      OldCounter(2) = OldCounter(2) + IntA + IntB + IntC + IntD
    END DO
  END DO
  OldTimes(2) = RealTime() - OldTimes(2)


  NewTimes(2) = RealTime() 
  CALL ListInitElementKeyword( LogA_h,'Material','Integer Value A' )
  CALL ListInitElementKeyword( LogB_h,'Material','Integer Value B' )
  CALL ListInitElementKeyword( LogC_h,'Material','Integer Value C')
  CALL ListInitElementKeyword( LogD_h,'Body Force','Integer Value D' )

  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)

      IntA = ListGetElementInteger( LogA_h, Element, Found )
      IntB = ListGetElementInteger( LogB_h, Element, Found )
      IntC = ListGetElementInteger( LogC_h, Found = Found )
      IntD = ListGetElementInteger( LogD_h )

      NewCounter(2) = NewCounter(2) + IntA + IntB + IntC + IntD
    END DO
  END DO
  NewTimes(2) = RealTime() - NewTimes(2)



  CALL Info('KeywordCompare','Starting string comparision')
  OldTimes(3) = RealTime()
  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)
      Material => GetMaterial()
      BodyForce => GetBodyForce()

      StrA = ListGetString( Material,'String Value A',Found )
      StrB = ListGetString( Material,'String Value B',Found )

      IF( StrA == 'str target a' ) OldCounter(3) = OldCounter(3) + 1 
      IF( StrB == 'str target b' ) OldCounter(3) = OldCounter(3) + 1 
    END DO
  END DO
  OldTimes(3) = RealTime() - OldTimes(3)



  NewTimes(3) = RealTime()
  CALL ListInitElementKeyword( StrA_h,'Material','String Value A' )
  CALL ListInitElementKeyword( StrB_h,'Material','String Value B' )

  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)

      StrA = ListGetElementString( StrA_h, Element, Found )
      StrB = ListGetElementString( StrB_h, Element, Found )

      IF( StrA == 'str target a' ) NewCounter(3) = NewCounter(3) + 1 
      IF( StrB == 'str target b' ) NewCounter(3) = NewCounter(3) + 1 
    END DO
  END DO
  NewTimes(3) = RealTime() - NewTimes(3)


  NewTimes(4) = RealTime()   
  DO i=1,TimeSweeps
    DO t=1,NoActive
      Element => GetActiveElement(t)

      LogA = ListCompareElementString( StrA_h, 'str target a',Element, Found )
      LogB = ListCompareElementString( StrB_h, 'str target b',Element, Found )

      IF( LogA ) NewCounter(4) = NewCounter(4) + 1 
      IF( LogB ) NewCounter(4) = NewCounter(4) + 1 
    END DO
  END DO
  NewTimes(4) = RealTime() - NewTimes(4) 

  ! There is no old comparison for ListCompareelementString
  OldTimes(4) = OldTimes(3)
  OldCounter(4) = OldCounter(3)

  !PRINT *,'Old Counter:',OldCounter
  !PRINT *,'New Counter:',NewCounter

  PRINT *,'Speed-up for Logical:           ',OldTimes(1)/NewTimes(1)
  PRINT *,'Speed-up for Integer:           ',OldTimes(2)/NewTimes(2)
  PRINT *,'Speed-up for String:            ',OldTimes(3)/NewTimes(3)
  PRINT *,'Speed-up for String comparison: ',OldTimes(4)/NewTimes(4)

  PseudoNorm = 1.0_dp * SUM(NewCounter) / SUM(OldCounter)
  PRINT *,'PseudoNorm:',PseudoNorm

  Solver % Variable % Values = PseudoNorm
  Solver % Variable % Norm = PseudoNorm
   
      
  CALL Info( 'KeywordCompare','Add done')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')


END SUBROUTINE KeywordCompare
