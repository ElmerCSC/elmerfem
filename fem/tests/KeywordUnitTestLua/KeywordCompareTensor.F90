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
  TYPE(ValueHandle_t) :: RealValLua_h, RealValMatc_h
  INTEGER :: i, key
  REAL(KIND=dp) :: OldCounter(10),NewCounter(10), DiffCounter(10)
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp) :: PseudoNorm, Val
  
  
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), OldValAtIp_lua(:,:), OldValAtIp_matc(:,:)
  REAL(KIND=dp), POINTER :: RealVal_lua(:,:,:) => NULL(), realval_matc(:,:,:), &
      ValAtIp_lua(:,:) => NULL(), ValAtIp_matc(:,:) => NULL()
  REAL(KIND=dp) :: DetJ
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  INTEGER :: elem, n, fdim, n1, n2, i1, i2, rdim
  INTEGER, POINTER :: Indexes(:)
  LOGICAL :: Stat
  CHARACTER(LEN=20) :: KeywordName_lua, keywordname_matc
  
  SAVE Nodes

  
  CALL Info( 'KeywordCompare','')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')
  CALL Info( 'KeywordCompare','Comparing keyword strategies for real values')

  Mesh => GetMesh()

  n = 20
  ALLOCATE( Basis(n), OldValAtIp_lua(3,3), OldValAtIp_matc(3,3), &
      realval_lua(1,1,1), realval_matc(1,1,1) ) 
  
  
  NoActive = GetNOFActive()

  OldCounter = 0.0_dp
  NewCounter = 0.0_dp
  DiffCounter = 0.0_dp


  CALL Info('KeywordCompare','Starting real value keyword comparision')

  DO key = 1, 4
    KeywordName_lua = 'Float Value '//TRIM(I2S(key))
    KeywordName_matc = 'Float Value '//TRIM(I2S(key+4))
    PRINT *,'Testing keywords: '//TRIM(KeywordName_lua) // ' and ' // trim(KeywordName_matc)
    ! CALL ListInitElementKeyword( RealVal_h,'Material',KeywordName_lua )

    DO elem=1,NoActive
      Element => GetActiveElement(elem)
      Material => GetMaterial()

      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()
      indexes => Element % NodeIndexes


      CALL ListGetRealArray( Material,KeywordName_lua,RealVal_lua,N,Indexes,Found )
      CALL ListGetRealArray( Material,KeywordName_matc,RealVal_matc,N,Indexes,Found )

      IF(.NOT. Found ) CYCLE

      n1 = SIZE( RealVal_lua, 1 )
      n2 = SIZE( RealVal_lua, 2 ) 


      IP = GaussPoints( Element )
      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )

        DO i1=1,n1
          DO i2=1,n2
            OldValAtIp_lua(i1,i2) = SUM( Basis(1:n) * RealVal_lua(i1,i2,1:n) )
            OldValAtIp_matc(i1,i2) = SUM( Basis(1:n) * RealVal_matc(i1,i2,1:n) )
          END DO
        END DO

        newcounter(key) = Newcounter(key)+  SUM(OldValAtIp_lua(1:n1,1:n2)) * IP % s(t) 
        oldcounter(key) = oldcounter(key)+  SUM(OldValAtIp_matc(1:n1,1:n2)) * IP % s(t) 
        DiffCounter(key) = DiffCounter(key) +  SUM(OldValAtIp_lua(1:n1,1:n2)-OldValAtIp_matc(1:n1,1:n2)) * IP % s(t) 
      END DO
    END DO

    PRINT *,'Compare integrals:', 2.0_dp*DiffCounter(key)  / (NewCounter(key) + OldCounter(key))
  END DO
  PseudoNorm = sum(2.0_dp*DiffCounter(1:4) / (NewCounter(1:4) + OldCounter(1:4)))

  DO key = 1, 4
    KeywordName_lua = 'Float Value '//TRIM(I2S(key))
    KeywordName_matc = 'Float Value '//TRIM(I2S(key+4))
    PRINT *,'Testing keywords using handles: '//TRIM(KeywordName_lua) // ' and ' // trim(KeywordName_matc)

    CALL ListInitElementKeyword( RealValLua_h,'Material',KeywordName_lua )
    CALL ListInitElementKeyword( RealValMatc_h,'Material',KeywordName_matc )

    DO elem=1,NoActive
      Element => GetActiveElement(elem)
      Material => GetMaterial()

      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()
      indexes => Element % NodeIndexes



      IF(.NOT. Found ) CYCLE

      n1 = SIZE( RealVal_lua, 1 )
      n2 = SIZE( RealVal_lua, 2 ) 

      IP = GaussPoints( Element )
      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )
        Val = ListGetElementReal( RealValLua_h, Basis, Element, Found, Rdim = Rdim, Rtensor = ValAtIp_lua )
        Val = ListGetElementReal( RealValMatc_h, Basis, Element, Found, Rdim = Rdim, Rtensor = ValAtIp_matc )

        DO i1=1,n1
          DO i2=1,n2
            OldValAtIp_lua(i1,i2) = SUM( Basis(1:n) * RealVal_lua(i1,i2,1:n) )
            OldValAtIp_matc(i1,i2) = SUM( Basis(1:n) * RealVal_matc(i1,i2,1:n) )
          END DO
        END DO

        newcounter(key) = Newcounter(key)+  SUM(OldValAtIp_lua(1:n1,1:n2)) * IP % s(t) 
        oldcounter(key) = oldcounter(key)+  SUM(OldValAtIp_matc(1:n1,1:n2)) * IP % s(t) 
        DiffCounter(key) = DiffCounter(key) +  SUM(OldValAtIp_lua(1:n1,1:n2)-OldValAtIp_matc(1:n1,1:n2)) * IP % s(t) 
      END DO
    END DO

    PRINT *,'Compare integrals:', 2.0_dp*DiffCounter(key)  / (NewCounter(key) + OldCounter(key))
  END DO


  PseudoNorm = PseudoNorm + sum(2.0_dp*DiffCounter(1:4) / (NewCounter(1:4) + OldCounter(1:4)))
  PRINT *,'PseudoNorm:',PseudoNorm

  Solver % Variable % Values = PseudoNorm
  Solver % Variable % Norm = PseudoNorm
   
      
  CALL Info( 'KeywordCompare','Add done')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')


END SUBROUTINE KeywordCompare
