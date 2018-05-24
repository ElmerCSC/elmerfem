!------------------------------------------------------------------------------
!> Compare keyword strategies with and without a Handle. 
!> This test for real valued keywords.
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
  TYPE(ValueHandle_t) :: RealVal_h
  INTEGER :: i, TimeSweeps, key
  REAL(KIND=dp) :: OldCounter(10),NewCounter(10)
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp) :: RealVal_AtIp, OldTimes(10), NewTimes(10), PseudoNorm, Val


  REAL(KIND=dp), ALLOCATABLE :: Basis(:), RealVal(:)
  REAL(KIND=dp) :: DetJ
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  INTEGER :: elem, n
  INTEGER, POINTER :: Indexes(:)
  LOGICAL :: Stat
  CHARACTER(LEN=20) :: KeywordName
  
  SAVE Nodes


  
  CALL Info( 'KeywordCompare','')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')
  CALL Info( 'KeywordCompare','Comparing keyword strategies for real values')

  Mesh => GetMesh()

  n = 20
  ALLOCATE( Basis(n), RealVal(n) )
  
  
  NoActive = GetNOFActive()
  TimeSweeps = ListGetInteger( Solver % Values,'Timer Sweeps')

  OldTimes = 0.0_dp
  NewTimes = 0.0_dp
  OldCounter = 0.0_dp
  NewCounter = 0.0_dp


  CALL Info('KeywordCompare','Starting real value keyword comparision')
  OldTimes(1) = RealTime()

  DO key = 1, 10
    PRINT *,'Testing keyword number:',key
    KeywordName = 'Float Value '//TRIM(I2S(key))
    
    OldTimes(key) = RealTime()  
    DO i=1,TimeSweeps
      DO elem=1,NoActive
        Element => GetActiveElement(elem)
        Material => GetMaterial()
        
        CALL GetElementNodes( Nodes )
        n  = GetElementNOFNodes()
        indexes => Element % NodeIndexes
        
        RealVal(1:n) = ListGetReal( Material,KeywordName,n,Indexes,Found )

        IP = GaussPoints( Element )
        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )
          RealVal_AtIp = SUM( Basis(1:n) * RealVal(1:n) )
          OldCounter(key) = OldCounter(key) +  RealVal_AtIp * IP % s(t)
        END DO
      END DO
    END DO
    OldTimes(key) = RealTime() - OldTimes(key)
    

    NewTimes(key) = RealTime()  
    
    CALL ListInitElementKeyword( RealVal_h,'Material',KeywordName )
      
    DO i=1,TimeSweeps
      DO elem=1,NoActive
        Element => GetActiveElement(elem)

        CALL GetElementNodes( Nodes )
        n  = GetElementNOFNodes()      

        IP = GaussPoints( Element )
        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )
          RealVal_AtIp = ListGetElementReal( RealVal_h, Basis, Element, Found )
          NewCounter(key) = NewCounter(key) + RealVal_AtIp * IP % s(t)
        END DO
      END DO
    END DO
    NewTimes(key) = RealTime() - NewTimes(key)

    PRINT *,'Compare integrals:', OldCounter(key), NewCounter(key), &
        NewCounter(key)/OldCounter(key)
    PRINT *,'Compare time consumptions:',OldTimes(key),NewTimes(key),&
        OldTimes(key)/NewTimes(key)
  END DO


  PseudoNorm = 1.0_dp * SUM(NewCounter) / SUM(OldCounter)
  PRINT *,'PseudoNorm:',PseudoNorm

  Solver % Variable % Values = PseudoNorm
  Solver % Variable % Norm = PseudoNorm
   
      
  CALL Info( 'KeywordCompare','Add done')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')


END SUBROUTINE KeywordCompare
