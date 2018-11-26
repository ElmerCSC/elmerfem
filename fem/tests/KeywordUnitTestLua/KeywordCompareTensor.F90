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
      ValAtIp_lua(:,:) => NULL(), ValAtIp_matc(:,:) => NULL(), realval_lua_vec(:,:) => NULL(), &
      RealVal_matc_vec(:,:) => NULL()
  REAL(KIND=dp) :: DetJ
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  INTEGER :: elem, n, fdim, n1, n2, i1, i2, rdim
  INTEGER, POINTER :: Indexes(:)
  LOGICAL :: Stat
  CHARACTER(LEN=13) :: KeywordName_lua, keywordname_matc
  CHARACTER(len=8) :: suffixes = 'abcdefgh'
  
  SAVE Nodes

  
  CALL Info( 'KeywordCompare','')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')
  CALL Info( 'KeywordCompare','Comparing keyword strategies for real values')

  Mesh => GetMesh()

  n = 20
  ALLOCATE( Basis(n), OldValAtIp_lua(3,3), OldValAtIp_matc(3,3), &
      realval_lua(1,1,1), realval_matc(1,1,1), realval_lua_vec(1,1), realval_matc_vec(1,1) ) 
  
  
  NoActive = GetNOFActive()

  PseudoNorm = 0.0_dp

  CALL Info('KeywordCompare','Starting real value keyword comparision')

  OldCounter = 0.0_dp
  NewCounter = 0.0_dp
  DiffCounter = 0.0_dp
  DO key = 1, 4, 2 ! Test odd keys with ListGetRealVec TODO: should fail {{{
    call set_keyword_names()
    PRINT *,'Testing keywords: '//TRIM(KeywordName_lua) // ' and ' // trim(KeywordName_matc)
    ! CALL ListInitElementKeyword( RealVal_h,'Material',KeywordName_lua )

    ! Skipping for now since it fails
    
    DO elem=1,NoActive
      Element => GetActiveElement(elem)
      Material => GetMaterial()

      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()
      indexes => Element % NodeIndexes
      ALLOCATE(realval_matc_vec(3,N), realval_lua_vec(3,N))

      Found = .false.
      CALL ListGetRealVector( Material, KeywordName_lua,  RealVal_lua_vec,  N, Indexes, Found )
      IF(.NOT. Found ) then
        print *, '>'// KeywordName_Lua //'< not found'
        CYCLE
      end if
      CALL ListGetRealVector( Material, KeywordName_matc, RealVal_matc_vec, N, Indexes, Found )
      IF(.NOT. Found ) then
        print *, '>'// KeywordName_matc //'< not found'
        CYCLE
      end if


      n1 = SIZE( RealVal_lua_vec, 1 )
      n2 = 1

      IP = GaussPoints( Element )
      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )

        DO i1=1,n1
          DO i2=1,n2
            OldValAtIp_lua(i1,i2) = SUM( Basis(1:n) * RealVal_lua_vec(i1,1:n) )
            OldValAtIp_matc(i1,i2) = SUM( Basis(1:n) * RealVal_matc_vec(i1,1:n) )
          END DO
        END DO

        newcounter(key) = Newcounter(key)+  SUM(OldValAtIp_lua(1:n1,1:n2)) * IP % s(t) 
        oldcounter(key) = oldcounter(key)+  SUM(OldValAtIp_matc(1:n1,1:n2)) * IP % s(t) 
        DiffCounter(key) = DiffCounter(key) +  SUM(OldValAtIp_lua(1:n1,1:n2)-OldValAtIp_matc(1:n1,1:n2)) * IP % s(t) 
      END DO
    END DO

    call printit()
    PseudoNorm = PseudoNorm + (2.0_dp*DiffCounter(key) / (NewCounter(key) + OldCounter(key)))
  END DO ! }}}

  DiffCounter = 0.0_dp
  NewCounter = 0.0_dp
  OldCOunter = 0.0_dp
  DO key = 1, 4 ! test all values with ListGetRealArray {{{
    call set_keyword_names()
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

    call printit()
    PseudoNorm = PseudoNorm + (2.0_dp*DiffCounter(key) / (NewCounter(key) + OldCounter(key)))
  END DO ! }}}
  ! PseudoNorm = PseudoNorm + sum(2.0_dp*DiffCounter(1:4) / (NewCounter(1:4) + OldCounter(1:4)))

  DiffCounter = 0.0_dp
  NewCounter = 0.0_dp
  OldCOunter = 0.0_dp
  DO key = 1, 4 ! use handles {{{ 
    call set_keyword_names()
    ! KeywordName_lua = 'Float Value '//TRIM(I2S(key))
    ! KeywordName_matc = 'Float Value '//TRIM(I2S(key+4))
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

    call printit()
    PseudoNorm = PseudoNorm + (2.0_dp*DiffCounter(key) / (NewCounter(key) + OldCounter(key)))
  END DO ! }}}


  ! PseudoNorm = PseudoNorm + sum(2.0_dp*DiffCounter(1:4) / (NewCounter(1:4) + OldCounter(1:4)))
  PRINT *,'PseudoNorm:',PseudoNorm

  Solver % Variable % Values = PseudoNorm
  Solver % Variable % Norm = PseudoNorm
   
      
  CALL Info( 'KeywordCompare','Add done')
  CALL Info( 'KeywordCompare','----------------------------------------------------------')

  contains
    subroutine set_keyword_names()
      KeywordName_lua = 'Float Value '//suffixes(key:key)
      KeywordName_matc = 'Float Value '//suffixes(key+4:key+4)
    end subroutine

    subroutine printit()
      print *,'Compare integrals (old, new, rel diff):', oldcounter(key), newcounter(key), &
          2.0_dp*DiffCounter(key)  / (NewCounter(key) + OldCounter(key))
    end subroutine

END SUBROUTINE KeywordCompare
