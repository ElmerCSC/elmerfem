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
  REAL(KIND=dp) :: OldTimes(10), NewTimes(10), PseudoNorm, Val
  
  
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), OldValAtIp(:,:)
  REAL(KIND=dp), POINTER :: RealVal(:,:,:) => NULL(), NewValAtIp(:,:)
  REAL(KIND=dp) :: DetJ
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  INTEGER :: elem, n, fdim, n1, n2, i1, i2, rdim
  INTEGER, POINTER :: Indexes(:)
  LOGICAL :: Stat
  CHARACTER(LEN=20) :: KeywordName
  
  SAVE Nodes, RealVal

  ALLOCATE( RealVal(1,1,1) )
  
  CALL Info( 'KeywordCompare','')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')
  CALL Info( 'KeywordCompare','Comparing keyword strategies for real values')

  Mesh => GetMesh()

  n = 20
  ALLOCATE( Basis(n), OldValAtIp(3,3) )
  
  
  NoActive = GetNOFActive()
  TimeSweeps = ListGetInteger( Solver % Values,'Timer Sweeps')

  OldTimes = 0.0_dp
  NewTimes = 0.0_dp
  OldCounter = 0.0_dp
  NewCounter = 0.0_dp


  CALL Info('KeywordCompare','Starting real value keyword comparision')
  OldTimes(1) = RealTime()

  DO key = 1, 10
    KeywordName = 'Float Value '//TRIM(I2S(key))
    PRINT *,'Testing keyword: '//TRIM(KeywordName)
    
    OldTimes(key) = RealTime()  
    DO i=1,TimeSweeps
      DO elem=1,NoActive
        Element => GetActiveElement(elem)
        Material => GetMaterial()
        
        CALL GetElementNodes( Nodes )
        n  = GetElementNOFNodes()
        indexes => Element % NodeIndexes

        
        CALL ListGetRealArray( Material,KeywordName,RealVal,N,Indexes,Found )

        IF(.NOT. Found ) CYCLE
        
        n1 = SIZE( RealVal, 1 )
        n2 = SIZE( RealVal, 2 ) 

        
        IP = GaussPoints( Element )
        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

          DO i1=1,n1
            DO i2=1,n2
              OldValAtIp(i1,i2) = SUM( Basis(1:n) * RealVal(i1,i2,1:n) )
            END DO
          END DO

          OldCounter(key) = OldCounter(key) +  SUM(OldValAtIp(1:n1,1:n2)) * IP % s(t) 
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
          Val = ListGetElementReal( RealVal_h, Basis, Element, Found, Rdim = Rdim, Rtensor = NewValAtIp )          
          IF(.NOT. Found ) CYCLE
          
          IF( Rdim == 0 ) CALL Fatal('KeywordCompare','All test keywords here should be tensors!')
          
          NewCounter(key) = NewCounter(key) + SUM( NewValAtIp ) * IP % s(t)
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
